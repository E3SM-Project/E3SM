module ActiveLayerMod
  
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines for calculation of active layer dynamics
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_const_mod   , only : SHR_CONST_TKFRZ
  use elm_varctl      , only : iulog, spinup_state, use_polygonal_tundra
  use TemperatureType , only : temperature_type
  use CanopyStateType , only : canopystate_type
  use GridcellType    , only : grc_pp       
  use ColumnType      , only : col_pp
  use ColumnDataType  , only : col_es, col_ws
  use LandunitType    , only : lun_pp
  use landunit_varcon , only : ilowcenpoly, iflatcenpoly, ihighcenpoly
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: alt_calc
  !-----------------------------------------------------------------------
  
contains

  !-----------------------------------------------------------------------
  subroutine alt_calc(num_soilc, filter_soilc, &
       temperature_vars, canopystate_vars) 
    !
    ! !DESCRIPTION:
    !  define active layer thickness similarly to frost_table, except set as deepest thawed layer and define on nlevgrnd
    !  also update annual maxima, and keep track of prior year for rooting memory
    !
    !  Note (WJS, 6-12-13): This routine just operates over active soil columns. However,
    !  because of its placement in the driver sequence (it is called very early in each
    !  timestep, before weights are adjusted and filters are updated), it effectively operates
    !  over the active columns from the previous time step. (This isn't really an issue now,
    !  but could become one when we have dynamic landunits.) As long as the output variables
    !  computed here are initialized to valid values over non-active points, this shouldn't be
    !  a problem - it will just mean that values are not quite what they should be in the
    !  first timestep a point becomes active. Currently, it seems that these variables ARE
    !  initialized to valid values, so I think this is okay.
    !
    ! !USES:
    use shr_const_mod    , only : SHR_CONST_TKFRZ
    use elm_varpar       , only : nlevgrnd
    use elm_time_manager , only : get_curr_date, get_step_size
    use elm_varctl       , only : iulog
    use elm_varcon       , only : zsoi, dzsoi, zisoi
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(temperature_type) , intent(in)    :: temperature_vars
    type(canopystate_type) , intent(inout) :: canopystate_vars
    !
    ! !LOCAL VARIABLES:
    integer  :: c, j, fc, g                             ! counters
    integer  :: alt_ind                                 ! index of base of active layer
    integer  :: year                                    ! year (0, ...) for nstep+1
    integer  :: mon                                     ! month (1, ..., 12) for nstep+1
    integer  :: day                                     ! day of month (1, ..., 31) for nstep+1
    integer  :: sec                                     ! seconds into current date for nstep+1
    integer  :: dtime                                   ! time step length in seconds
    integer  :: k_frz                                   ! index of first nonfrozen soil layer
    logical  :: found_thawlayer                         ! used to break loop when first unfrozen layer reached
    real(r8) :: t1, t2, z1, z2, orig_excess, old_mfrac  ! temporary variables
    real(r8), dimension(nlevgrnd) :: melt_profile       ! profile of melted excess ice
    !-----------------------------------------------------------------------

    associate(                                                                &
         t_soisno             =>    col_es%t_soisno        ,    & ! Input:   [real(r8) (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

         alt                  =>    canopystate_vars%alt_col             ,      & ! Output:  [real(r8) (:)   ]  current depth of thaw
         altmax               =>    canopystate_vars%altmax_col          ,      & ! Output:  [real(r8) (:)   ]  maximum annual depth of thaw
         altmax_lastyear      =>    canopystate_vars%altmax_lastyear_col ,      & ! Output:  [real(r8) (:)   ]  prior year maximum annual depth of thaw
         altmax_1989          =>    canopystate_vars%altmax_1989_col     ,      & ! Output:  [real(r8) (:)   ]  maximum ALT in 1989
         altmax_ever          =>    canopystate_vars%altmax_ever_col     ,      & ! Output:  [real(r8) (:)   ]  maximum thaw depth since initialization
         alt_indx             =>    canopystate_vars%alt_indx_col        ,      & ! Output:  [integer  (:)   ]  current depth of thaw
         altmax_indx          =>    canopystate_vars%altmax_indx_col     ,      & ! Output:  [integer  (:)   ]  maximum annual depth of thaw
         altmax_lastyear_indx =>    canopystate_vars%altmax_lastyear_indx_col , & ! Output:  [integer  (:)   ]  prior year maximum annual depth of thaw
         altmax_1989_indx     =>    canopystate_vars%altmax_1989_indx_col,      & ! Output:  [integer  (:)   ]  index of maximum ALT in 1989
         altmax_ever_indx     =>    canopystate_vars%altmax_ever_indx_col,      & ! Output:  [integer  (:)   ]  maximum thaw depth since initialization
         excess_ice           =>    col_ws%excess_ice                    ,      & ! Input/output:[real(r8) (:,:)]  depth variable excess ice content in soil column (-)
         rmax                 =>    col_ws%iwp_microrel                  ,      & ! Output:  [real(r8) (:)   ]  ice wedge polygon microtopographic relief (m)
         vexc                 =>    col_ws%iwp_exclvol                   ,      & ! Output:  [real(r8) (:)   ]  ice wedge polygon excluded volume (m)
         ddep                 =>    col_ws%iwp_ddep                      ,      & ! Output:  [real(r8) (:)   ]  ice wedge polygon depression depth (m)
         subsidence           =>    col_ws%iwp_subsidence                ,      & ! Input/output:[real(r8)(:)]  ice wedge polygon subsidence (m)
         frac_melted          =>    col_ws%frac_melted                          & ! Input/output:[real(r8)(:)]  fraction of layer that has ever melted (-)
         )

      ! on a set annual timestep, update annual maxima
      ! make this 1 January for NH columns, 1 July for SH columns
      call get_curr_date(year, mon, day, sec)
      dtime =  get_step_size()
      if ( (mon .eq. 1) .and. (day .eq. 1) .and. ( sec / dtime .eq. 1) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            g = col_pp%gridcell(c)
            if ( grc_pp%lat(g) > 0. ) then
               
               altmax_lastyear(c) = altmax(c)
               altmax_lastyear_indx(c) = altmax_indx(c)
               altmax(c) = 0.
               altmax_indx(c) = 0
            endif
         end do
      endif
      if ( (mon .eq. 7) .and. (day .eq. 1) .and. ( sec / dtime .eq. 1) ) then
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            g = col_pp%gridcell(c)
            if ( grc_pp%lat(g) <= 0. ) then 
               altmax_lastyear(c) = altmax(c)
               altmax_lastyear_indx(c) = altmax_indx(c)
               altmax(c) = 0.
               altmax_indx(c) = 0
            endif
         end do
      endif

      do fc = 1,num_soilc
         c = filter_soilc(fc)

         ! calculate alt for a given timestep
         ! start from base of soil and search upwards for first thawed layer.
         ! note that this will put talik in with active layer
         ! a different way of doing this could be to keep track of how long a given layer has ben frozen for,
         ! and define ALT as the first layer that has been frozen for less than 2 years.
         if (t_soisno(c,nlevgrnd) > SHR_CONST_TKFRZ ) then
            alt(c) = zsoi(nlevgrnd)
            alt_indx(c) = nlevgrnd
         else
            k_frz=0
            found_thawlayer = .false.
            do j=nlevgrnd-1,1,-1
               if ( ( t_soisno(c,j) > SHR_CONST_TKFRZ ) .and. .not. found_thawlayer ) then
                  k_frz=j
                  found_thawlayer = .true.
               endif
            end do

            if ( k_frz > 0 ) then
               ! define active layer as the depth at which the linearly interpolated temperature line intersects with zero
               z1 = zsoi(k_frz)
               z2 = zsoi(k_frz+1)
               t1 = t_soisno(c,k_frz)
               t2 = t_soisno(c,k_frz+1)
               alt(c) = z1 + (t1-SHR_CONST_TKFRZ)*(z2-z1)/(t1-t2)
               alt_indx(c) = k_frz
            else
               alt(c)=0._r8
               alt_indx(c) = 0
            endif
         endif

         ! if appropriate, update maximum annual active layer thickness
         if (alt(c) > altmax(c)) then
            altmax(c) = alt(c)
            altmax_indx(c) = alt_indx(c)
         endif
         if (alt(c) > altmax_ever(c)) then
            if (spinup_state .eq. 0) then
                altmax_ever(c) = alt(c)
                altmax_ever_indx(c) = alt_indx(c)
            else !overwrite if in spinup
                altmax_ever(c) = 0._r8
                altmax_ever_indx(c) = 0
            endif
         endif

         if (use_polygonal_tundra) then
            ! special case for year = 1989. For now it is the assumed baseline year for 
            ! changes in polygonal ground
            if (year .eq. 1989) then
               altmax_1989(c) = altmax(c)
               altmax_1989_indx(c) = altmax_indx(c)
            endif

           ! update subsidence based on change in ALT
           ! melt_profile stores the amount of excess_ice
           ! melted in this timestep.
           ! note that this may cause some unexpected results
           ! for taliks

           ! initialize melt_profile as zero
           melt_profile(:) = 0._r8

           do j = nlevgrnd,1,-1 ! note, this will go from bottom to top
              if (j .gt. k_frz + 1) then ! all layers below k_frz + 1 remain frozen
                melt_profile(j) = 0.0_r8
              else if (j .eq. k_frz + 1) then ! first layer below the 'thawed' layer
                ! need to check to see if the active layer thickness is is actually
                ! in this layer (and not between the midpoint of j_frz and bottom interface
                ! or else inferred melt will be negative
                ! also note: only have ice to melt if alt has never been this deep, otherwise
                ! ice will continue to be removed each time step the alt remains in this layer
                if ((alt(c)-zisoi(j-1)) .ge. 0._r8 .and. (alt(c) .eq. altmax_ever(c)) .and. (frac_melted(c,j) .lt. 1._r8)) then
                  orig_excess = (1._r8/(1._r8-frac_melted(c,j))) * excess_ice(c,j)
                  old_mfrac = frac_melted(c,j)
                  ! update frac melted
                  frac_melted(c,j) = min(max(frac_melted(c,j), (alt(c)-zisoi(j-1))/dzsoi(j)),1._r8)
                  melt_profile(j) = orig_excess*(frac_melted(c,j) - old_mfrac)
                  excess_ice(c,j) = excess_ice(c,j) - melt_profile(j)
                else
                  melt_profile(j) = 0._r8 ! no melt
                end if
              else if (j .eq. k_frz) then
                if (alt(c) .eq. altmax_ever(c) .and. (frac_melted(c,j) .lt. 1._r8)) then
                  orig_excess = (1._r8/(1._r8 - frac_melted(c,j))) * excess_ice(c,j)
                  old_mfrac = frac_melted(c,j)
                  ! update frac_melted:
                  frac_melted(c,j) = min(max(frac_melted(c,j), (alt(c)-zsoi(j-1))/dzsoi(j)),1._r8)
                  ! remove ice, only if alt has never been this deep before:
                  melt_profile(j) = orig_excess*(frac_melted(c,j) - old_mfrac)
                  excess_ice(c,j) = excess_ice(c,j) - melt_profile(j)
                else
                  melt_profile(j) = 0._r8
                end if
              else !
                 melt_profile(j) = excess_ice(c,j)
                 ! remove melted excess ice
                 excess_ice(c,j) = 0._r8
              end if
              ! calculate subsidence at this layer:
              melt_profile(j) = melt_profile(j) * dzsoi(j)
           end do

           ! subsidence is integral of melt profile:
           if ((year .ge. 1989) .and. (altmax_ever(c) .ge. altmax_1989(c))) then
              subsidence(c) = subsidence(c) + sum(melt_profile)
           end if

           ! limit subsidence to 0.4 m
           subsidence(c) = min(0.4_r8, subsidence(c))

           ! update ice wedge polygon microtopographic parameters if in polygonal ground
           if (lun_pp%ispolygon(col_pp%landunit(c))) then
             if (lun_pp%polygontype(col_pp%landunit(c)) .eq. ilowcenpoly) then
               rmax(c) = 0.4_r8
               vexc(c) = 0.2_r8
               ddep(c) = max(0.05_r8, 0.15_r8 - 0.25_r8*subsidence(c))
             elseif (lun_pp%polygontype(col_pp%landunit(c)) .eq. iflatcenpoly) then
               rmax(c) = min(0.4_r8, 0.1_r8 + 0.75_r8*subsidence(c))
               vexc(c) = min(0.2_r8, 0.05_r8 + 0.375_r8*subsidence(c))
               ddep(c) = min(0.05_r8, 0.01_r8 + 0.1_r8*subsidence(c))
             elseif (lun_pp%polygontype(col_pp%landunit(c)) .eq. ihighcenpoly) then
               rmax(c) = 0.4_r8
               vexc(c) = 0.2_r8
               ddep(c) = 0.05_r8
             else
               !call endrun !<- TODO: needed? Potential way to prevent unintended updating of microtopography
               ! if polygonal ground is misspecified on surface file.
             endif
           endif
         endif
       end do

    end associate 

  end subroutine alt_calc
  
end module ActiveLayerMod
