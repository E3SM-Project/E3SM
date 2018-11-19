module RootDynamicsMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module holding routines used for determining fine root distribution for all pfts.
  ! Includes dynamic root depth for crops
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use clm_time_manager    , only : get_step_size
  use clm_varpar          , only : nlevsoi, nlevgrnd
  use clm_varctl          , only : use_vertsoilc
  use decompMod           , only : bounds_type
  use pftvarcon           , only : noveg, npcropmin, roota_par, rootb_par, root_dmx, evergreen
  use ColumnType          , only : col_pp 
  use VegetationType      , only : veg_pp
  use CanopyStateType     , only: canopystate_type
  use CNStateType         , only : cnstate_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use EnergyFluxType      , only: energyflux_type
  use SoilStateType       , only : soilstate_type
  use CropType            , only : crop_type
  use SimpleMathMod       , only : array_normalization
  use RootBiophysMod      , only : init_vegrootfr

  ! !PUBLIC TYPES:
  implicit none
  save
  private
  public :: RootDynamics
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  !
  subroutine RootDynamics(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       canopystate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars,  &
       cnstate_vars, crop_vars,  energyflux_vars, soilstate_vars)
    !
    ! !DESCRIPTION:
    ! This routine determine the fine root distribution
    ! Needs to be called after the photosynthesis calculation
    ! May need to update other subroutines that use the fixed root profile for calculations
    ! i.e. VerticalProfileMod
    !
    ! !USES:


    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds             ! bounds
    integer                  , intent(in)    :: num_soilc
    integer                  , intent(in)    :: filter_soilc(:)
    integer                  , intent(in)    :: num_soilp          ! number of soil pfts in filter
    integer                  , intent(in)    :: filter_soilp(:)    ! filter for soil pfts
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(carbonflux_type)    , intent(in)    :: carbonflux_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars !
    type(crop_type)          , intent(in)    :: crop_vars
    type(energyflux_type)    , intent(in)    :: energyflux_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars

    !
    ! !LOCAL VARIABLES:

    integer  :: f,c,p,lev,j                                    ! indices
    real(r8) :: dt                                             ! radiation time step delta t (seconds)
    real(r8) :: w_limit(bounds%begp:bounds%endp)               ! soil water weighting factor
    real(r8) :: rswa(bounds%begp:bounds%endp,1:nlevgrnd)       ! soil water availability in each soil layer
    real(r8) :: rsmn(bounds%begp:bounds%endp,1:nlevgrnd)       ! soil nitrogen availability in each soil layer
    real(r8) :: sumrswa(bounds%begp:bounds%endp)               ! scaling soil water availability in each soil layer
    real(r8) :: sumrsmn(bounds%begp:bounds%endp)               ! scaling  soil mineral N availability in each soil layer
    real(r8) :: frootc_dz(bounds%begp:bounds%endp, 1:nlevgrnd) ! root carbon in each soil layer (gC)
    real(r8) :: sumfrootc(bounds%begp:bounds%endp)             ! fine root carbon total before turnover in each step
    real(r8) :: rootfr_coarse(bounds%begp:bounds%endp, 1:nlevgrnd) ! coarse root distribution
    real(r8) :: psi                                            ! soil moisture potential
    real(r8) :: maxpsi                                         ! maximum soil moisture potential
    real(r8) :: new_growth, new_croot_growth                   ! new carbon allocated to roots this timestep
    real(r8), parameter :: minpsi = -1.5_r8                    ! minimum soil moisture potential - permanent wilting point (MPa)
    real(r8), parameter :: soil_water_factor_min = 0.9_r8
    real(r8), parameter :: exp_decay_factor = 3._r8

    !-----------------------------------------------------------------------
    ! Assign local pointers to derived type arrays (in)
    associate(&
         ivt                    => veg_pp%itype                                , & ! Input  :  [integer (:)]  pft vegetation type
         pcolumn                => veg_pp%column                               , & ! Input  :  [integer (:)]  pft's column index
         croplive               => crop_vars%croplive_patch                    , & ! Input  :  [logical (:)]  flag, true if planted, not harvested
         cpool_to_frootc        => carbonflux_vars%cpool_to_frootc_patch       , & ! Input  :  [real(r8) (:)] allocation to fine root C (gC/m2/s)
         cpool_to_frootc_storage=> carbonflux_vars%cpool_to_frootc_storage_patch, & ! Input:  [real(r8) (:)] allocation to fine root C storage (gC/m2/s)
         frootc_xfer_to_frootc  => carbonflux_vars%frootc_xfer_to_frootc_patch , & ! Input  :  [real(r8) (:)] fine root C growth from storage (gC/m2/s)
         onset_flag             => cnstate_vars%onset_flag_patch               , & ! Input  :  [real(r8) (:)] onset flag
         dormant_flag           => cnstate_vars%dormant_flag_patch             , & ! Input  :  [real(r8) (:)]  dormancy flag
         root_depth             => soilstate_vars%root_depth_patch             , & ! InOut  :  [real(r8) (:)] current root depth
         dz                     => col_pp%dz                                   , & ! Input  :  layer thickness (m)  (-nlevsno+1:nlevgrnd)
         zi                     => col_pp%zi                                   , & ! Input  :  interface level below a "z" level (m) (-nlevsno+0:nlevgrnd)
         nlevbed                => col_pp%nlevbed                              , & ! Input  :  # levels to bedrock
         rootfr                 => soilstate_vars%rootfr_patch                 , & ! Output :  [real(r8) (:,:)]  fraction of roots in each soil layer
         sucsat                 => soilstate_vars%sucsat_col                   , & ! Input  :  minimum soil suction (mm)
         soilpsi                => soilstate_vars%soilpsi_col                  , & ! Input  :  soil water potential in each soil layer (MPa)
         rresis                 => energyflux_vars%rresis_patch                , & ! Input  :  [real(r8) (:,:) ]  root soil water stress (resistance) by layer (0-1)
         sminn_vr               => nitrogenstate_vars%sminn_vr_col             , & ! Iniput :  [real(r8) (:,:)]  (gN/m3) soil mineral N
         frootc                 => carbonstate_vars%frootc_patch               , & ! Input  :  [real(r8) (:)]  (gC/m2) fine root C
         hui                    => crop_vars%gddplant_patch                    , & ! Input  :  [real(r8) (:)]  =gdd since planting (gddplant)
         huigrain               => cnstate_vars%huigrain_patch                 , & ! Input  :  [real(r8) (:)]  same to reach vegetative maturity
         livecrootc             => carbonstate_vars%livecrootc_patch           , & !
         deadcrootc             => carbonstate_vars%deadcrootc_patch           , &
         cpool_to_livecrootc    => carbonflux_vars%cpool_to_livecrootc_patch   , & ! Input  :  [real(r8) (:)] allocation to coarse root C (gC/m2/s)
         cpool_to_livecrootc_storage => carbonflux_vars%cpool_to_livecrootc_storage_patch, & ! Input:  [real(r8) (:)] allocation to coarse root C storage (gC/m2/s)
         cpool_to_deadcrootc    => carbonflux_vars%cpool_to_deadcrootc_patch   , & ! Input  :  [real(r8) (:)] allocation to dead coarse root C (gC/m2/s)
         cpool_to_deadcrootc_storage => carbonflux_vars%cpool_to_deadcrootc_storage_patch, & ! Input:  [real(r8) (:)] allocation to dead coarse root C storage (gC/m2/s)
         livecrootc_xfer_to_livecrootc => carbonflux_vars%livecrootc_xfer_to_livecrootc_patch , & ! Input  :  [real(r8) (:)] coarse root C growth from storage (gC/m2/s)
         deadcrootc_xfer_to_deadcrootc => carbonflux_vars%deadcrootc_xfer_to_deadcrootc_patch ,  & ! Input  :  [real(r8) (:)] dead coarse root C growth from storage (gC/m2/s)
         altmax_lastyear         => canopystate_vars%altmax_lastyear_col         & ! Input: [real(r8) (:)   ]  maximum annual depth of thaw
         )

      ! set time steps
      dt = get_step_size()

      !initialize to 0
      w_limit(bounds%begp:bounds%endp)              = 0._r8
      sumrswa(bounds%begp:bounds%endp)              = 0._r8
      sumrsmn(bounds%begp:bounds%endp)              = 0._r8
      sumfrootc(bounds%begp:bounds%endp)            = 0._r8
      rswa(bounds%begp:bounds%endp,1:nlevgrnd)      = 0._r8
      rsmn(bounds%begp:bounds%endp,1:nlevgrnd)      = 0._r8
      frootc_dz(bounds%begp:bounds%endp,1:nlevgrnd) = 0._r8
      root_depth(bounds%begp:bounds%endp)           = 0._r8
      rootfr_coarse(bounds%begp:bounds%endp,1:nlevgrnd) = 0._r8

      !---------------------------------------------------------------
      ! Set root depth, dynamic for crops, fixed for other vegetation
      !---------------------------------------------------------------

      do f = 1, num_soilp
         p = filter_soilp(f)
         c = pcolumn(p)
         if (ivt(p) /= noveg) then
            if ((ivt(p)) >= npcropmin) then !skip generic crop types
               if (huigrain(p) > 0._r8) then
                  root_depth(p) = max(zi(c,2), min(hui(p)/huigrain(p)* root_dmx(ivt(p)), root_dmx(ivt(p))))
               end if
            else
               ! this can be changed to any depth (i.e. the maximum soil depth)
               root_depth(p) = min(altmax_lastyear(c), min(zi(c,nlevsoi), root_dmx(ivt(p))))
            end if
         else
            root_depth(p) = 0._r8
         end if
      end do

      !----------------------------------------------------------------
      !  ! calculate a weighting function by soil depth that depends on the
      ! fine root distribution per pft and depth and the pft weight on the column.
      ! This will be used to weight the temperature and water potential scalars
      ! for decomposition control.

      do f = 1,num_soilp
         p = filter_soilp(f)
         c = pcolumn(p)
         do j = 1,nlevsoi
            maxpsi = sucsat(c,j) * (-9.8e-6_r8)
            psi = min(soilpsi(c,j),maxpsi)
            if (psi > minpsi) then

               ! First calculate water in the root zone
               if (root_depth(p) >  zi(c,3) .and. (zi(c,j) <= root_depth(p) .or. &
                    (zi(c,j-1) < root_depth(p) .and. zi(c,j) > root_depth(p)))) then
                  w_limit(p) = w_limit(p) + max(rresis(p,j)*rootfr(p,j),0._r8)
                  w_limit(p) = min(soil_water_factor_min,w_limit(p))
               end if
               ! Calculate the water in each soil layer
               if (root_depth(p) >= zi(c,j) .or. &
                    (zi(c,j-1) < root_depth(p) .and. zi(c,j) > root_depth(p))) then
                  rswa(p,j) = max(0._r8, rresis(p,j))*dz(c,j)
               end if
            end if
            sumrswa(p) = sumrswa(p) + rswa(p,j)

            ! Calculate the nitrogen profile in each layer
            ! For now, the profile for each PFT is equivilent to the
            ! column profile, in the future, this could be changed to a weighted profile
            if (use_vertsoilc) then !for vertical soil profile
               rsmn(p,j) = sminn_vr(c,j)*dz(c,j)
            else ! need to calculate a profile, an exponential decay
               rsmn(p,j) = dz(c,j)*exp(-exp_decay_factor*zi(c,j))
            end if
            if (root_depth(p) >= zi(c,j).or. &
                 (zi(c,j-1) < root_depth(p) .and. zi(c,j) > root_depth(p))) then
               sumrsmn(p) = sumrsmn(p) + rsmn(p,j)
            end if
         end do
      end do

      !------------------------------------------------------------------
      ! Calculate coarse root fraction from standard root distribution
      !------------------------------------------------------------------

      call init_vegrootfr(bounds, nlevsoi, nlevgrnd, &
                   nlevbed(bounds%begc:bounds%endc), &
                   rootfr_coarse(bounds%begp:bounds%endp,1:nlevgrnd))

      !--------------------------------------------------------------------
      ! Now calculate the density of roots in each soil layer for each pft
      ! based on this timesteps growth
      !--------------------------------------------------------------------

      do f = 1, num_soilp
         p = filter_soilp(f) 
         c = pcolumn(p)
         new_growth = 0._r8
         new_croot_growth = 0._r8
         if (onset_flag(p) == 0._r8) then
            new_growth = (cpool_to_frootc(p) + frootc_xfer_to_frootc(p))*dt
            new_croot_growth = (cpool_to_livecrootc(p) + cpool_to_deadcrootc(p) + &
                                deadcrootc_xfer_to_deadcrootc(p) + livecrootc_xfer_to_livecrootc(p))*dt
         end if
         if (evergreen(ivt(p)) == 0._r8) then
            new_growth = new_growth + cpool_to_frootc_storage(p)*dt
            new_croot_growth = new_croot_growth + (cpool_to_livecrootc_storage(p) + &
                                                   cpool_to_deadcrootc_storage(p))*dt
         end if
         do lev=1,nlevsoi
            if (sumrswa(p) <= 0._r8 .or. sumrsmn(p) <= 0._r8) then
               ! when sumrswa or sumrsmn are less than or equal to 0 rootfr will not be updated
            else
               frootc_dz(p,lev) = (livecrootc(p) + deadcrootc(p) + frootc(p))*rootfr(p,lev) &
                    + new_croot_growth * rootfr_coarse(p,lev) &
                    + new_growth * ((1._r8 - w_limit(p)) * rswa(p,lev) / sumrswa(p) &
                    + w_limit(p) * rsmn(p,lev) / sumrsmn(p))
            end if

            sumfrootc(p) = sumfrootc(p) + frootc_dz(p,lev)

         end do
      end do

      !----------------------------------
      !Calculate root fraction
      !----------------------------------

      ! normalize the root fraction for each pft

      do f = 1, num_soilp
         p = filter_soilp(f)
         c = pcolumn(p)
         do lev = 1, nlevgrnd
            if (sumfrootc(p) > 0._r8) then
               rootfr(p,lev) = frootc_dz(p,lev)/sumfrootc(p)
            end if
         end do
      end do


    end associate

  end subroutine RootDynamics

end module RootDynamicsMod
