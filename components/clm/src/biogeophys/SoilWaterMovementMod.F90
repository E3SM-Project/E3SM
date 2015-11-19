module SoilWaterMovementMod

  !-----------------------------------------------------------------------
  ! DESCRIPTION
  ! module contains different subroutines to couple soil and root water interactions
  !
  ! created by Jinyun Tang, Mar 12, 2014
  implicit none
  save 
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: SoilWater            ! Calculate soil hydrology   
  public :: init_soilwater_movement
  public :: init_vsfm_condition_ids
  !
  ! !PRIVATE DATA MEMBERS:
  integer, parameter :: zengdecker_2009 = 0
  integer, parameter :: vsfm = 1
  integer :: soilroot_water_method     !0: use the Zeng and deck method, this will be readin from namelist in the future

  ! IDs to indentify the conditions for VSFM
  integer :: vsfm_cond_id_for_infil
  integer :: vsfm_cond_id_for_et
  integer :: vsfm_cond_id_for_dew
  integer :: vsfm_cond_id_for_drainage
  integer :: vsfm_cond_id_for_snow
  integer :: vsfm_cond_id_for_sublimation
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine init_soilwater_movement()
    !
    !DESCRIPTION
    !specify method for doing soil&root water interactions
    !
    use clm_varctl, only : use_vsfm
    use spmdMod,    only : mpicom, MPI_LOGICAL
    ! !ARGUMENTS:
    implicit none
    integer :: ier ! error status
    !------------------------------------------------------------------------------

    soilroot_water_method = zengdecker_2009

    ! Initialize condition id for VSFM
    vsfm_cond_id_for_infil         = -1
    vsfm_cond_id_for_et            = -1
    vsfm_cond_id_for_dew           = -1
    vsfm_cond_id_for_drainage      = -1
    vsfm_cond_id_for_snow          = -1
    vsfm_cond_id_for_sublimation   = -1

    ! GB-FIX-ME: The call to control_spmd() [in subroutine control_init()] before
    !            call to init_hydrology() would avoid the mpi broadcast
    call mpi_bcast (use_vsfm, 1, MPI_LOGICAL, 0, mpicom, ier)
    if (use_vsfm) soilroot_water_method = vsfm

  end subroutine init_soilwater_movement

  !-----------------------------------------------------------------------
  subroutine SoilWater(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, &
       waterflux_vars, waterstate_vars, temperature_vars, soil_water_retention_curve)
    !
    ! DESCRIPTION
    ! select one subroutine to do the soil and root water coupling
    !
    !USES
    use clm_varctl        , only : use_betr
    use shr_kind_mod      , only : r8 => shr_kind_r8
    use clm_varpar        , only : nlevsoi    
    use decompMod         , only : bounds_type   
    use abortutils        , only : endrun   
    use SoilHydrologyType , only : soilhydrology_type
    use SoilStateType     , only : soilstate_type
    use TemperatureType   , only : temperature_type
    use WaterFluxType     , only : waterflux_type
    use WaterStateType    , only : waterstate_type
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use clm_varcon        , only : denh2o, denice, watmin
    use ColumnType        , only : col        
    !
    ! !ARGUMENTS:
    implicit none     
    type(bounds_type)        , intent(in)    :: bounds                ! bounds
    integer                  , intent(in)    :: num_hydrologyc        ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_hydrologyc(:)  ! column filter for soil points
    integer                  , intent(in)    :: num_urbanc            ! number of column urban points in column filter
    integer                  , intent(in)    :: filter_urbanc(:)      ! column filter for urban points
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(soilstate_type)     , intent(inout) :: soilstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    character(len=32)                        :: subname = 'SoilWater'       ! subroutine name
    real(r8)                                 :: xs(bounds%begc:bounds%endc) !excess soil water above urban ponding limit

    integer  :: fc, c, j
    
    
    !------------------------------------------------------------------------------
    associate(                                                         &
      wa                 =>    soilhydrology_vars%wa_col             , & ! Input:  [real(r8) (:)   ] water in the unconfined aquifer (mm)
      dz                 =>    col%dz                                , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)    
      h2osoi_ice         =>    waterstate_vars%h2osoi_ice_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
      h2osoi_vol         =>    waterstate_vars%h2osoi_vol_col        , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
      h2osoi_liq         =>    waterstate_vars%h2osoi_liq_col          & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
    )

    select case(soilroot_water_method)

    case (zengdecker_2009)

       call soilwater_zengdecker2009(bounds, num_hydrologyc, filter_hydrologyc, &
            num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, &
            waterflux_vars, waterstate_vars, temperature_vars, soil_water_retention_curve)

    case (vsfm)
       call soilwater_vsfm(bounds, num_hydrologyc, filter_hydrologyc, &
            num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, &
            waterflux_vars, waterstate_vars, temperature_vars)
    case default

       call endrun(subname // ':: a SoilWater implementation must be specified!')          

    end select

    if(use_betr)then
    !a work around of the negative liquid water embarrassment, which is
    !critical for a meaningufl tracer transport in betr. Jinyun Tang, Jan 14, 2015

    do j = 1, nlevsoi-1
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          if (h2osoi_liq(c,j) < 0._r8) then
             xs(c) = watmin - h2osoi_liq(c,j)
          else
             xs(c) = 0._r8
          end if
          h2osoi_liq(c,j  ) = h2osoi_liq(c,j  ) + xs(c)
          h2osoi_liq(c,j+1) = h2osoi_liq(c,j+1) - xs(c)
       end do
    end do

    j = nlevsoi
    do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if (h2osoi_liq(c,j) < watmin) then
          xs(c) = watmin-h2osoi_liq(c,j)
        else
          xs(c) = 0._r8
       end if
       h2osoi_liq(c,j) = h2osoi_liq(c,j) + xs(c)
       wa(c) = wa(c) - xs(c)
    end do
    
    !update volumetric soil moisture for bgc calculation
    do j = 1, nlevsoi
       do fc = 1, num_hydrologyc
          c = filter_hydrologyc(fc)
          h2osoi_vol(c,j) = h2osoi_liq(c,j)/(dz(c,j)*denh2o) &
                            + h2osoi_ice(c,j)/(dz(c,j)*denice)
       enddo
    enddo
    endif
  end associate

  end subroutine SoilWater

  !-----------------------------------------------------------------------   
  subroutine soilwater_zengdecker2009(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, &
       waterflux_vars, waterstate_vars, temperature_vars, soil_water_retention_curve)
    !
    ! !DESCRIPTION:
    ! Soil hydrology
    ! Soil moisture is predicted from a 10-layer model (as with soil
    ! temperature), in which the vertical soil moisture transport is governed
    ! by infiltration, runoff, gradient diffusion, gravity, and root
    ! extraction through canopy transpiration.  The net water applied to the
    ! surface layer is the snowmelt plus precipitation plus the throughfall
    ! of canopy dew minus surface runoff and evaporation.
    ! CLM3.5 uses a zero-flow bottom boundary condition.
    !
    ! The vertical water flow in an unsaturated porous media is described by
    ! Darcy's law, and the hydraulic conductivity and the soil negative
    ! potential vary with soil water content and soil texture based on the work
    ! of Clapp and Hornberger (1978) and Cosby et al. (1984). The equation is
    ! integrated over the layer thickness, in which the time rate of change in
    ! water mass must equal the net flow across the bounding interface, plus the
    ! rate of internal source or sink. The terms of water flow across the layer
    ! interfaces are linearly expanded by using first-order Taylor expansion.
    ! The equations result in a tridiagonal system equation.
    !
    ! Note: length units here are all millimeter
    ! (in temperature subroutine uses same soil layer
    ! structure required but lengths are m)
    !
    ! Richards equation:
    !
    ! d wat      d     d wat d psi
    ! ----- = - -- [ k(----- ----- - 1) ] + S
    !   dt      dz       dz  d wat
    !
    ! where: wat = volume of water per volume of soil (mm**3/mm**3)
    ! psi = soil matrix potential (mm)
    ! dt  = time step (s)
    ! z   = depth (mm)
    ! dz  = thickness (mm)
    ! qin = inflow at top (mm h2o /s)
    ! qout= outflow at bottom (mm h2o /s)
    ! s   = source/sink flux (mm h2o /s)
    ! k   = hydraulic conductivity (mm h2o /s)
    !
    !                       d qin                  d qin
    ! qin[n+1] = qin[n] +  --------  d wat(j-1) + --------- d wat(j)
    !                       d wat(j-1)             d wat(j)
    !                ==================|=================
    !                                  < qin
    !
    !                 d wat(j)/dt * dz = qin[n+1] - qout[n+1] + S(j)
    !
    !                                  > qout
    !                ==================|=================
    !                        d qout               d qout
    ! qout[n+1] = qout[n] + --------- d wat(j) + --------- d wat(j+1)
    !                        d wat(j)             d wat(j+1)
    !
    !
    ! Solution: linearize k and psi about d wat and use tridiagonal
    ! system of equations to solve for d wat,
    ! where for layer j
    !
    !
    ! r_j = a_j [d wat_j-1] + b_j [d wat_j] + c_j [d wat_j+1]
    !
    ! !USES:
    use shr_kind_mod         , only : r8 => shr_kind_r8     
    use shr_const_mod        , only : SHR_CONST_TKFRZ, SHR_CONST_LATICE, SHR_CONST_G
    use decompMod            , only : bounds_type        
    use clm_varcon           , only : wimp,grav,hfus,tfrz
    use clm_varcon           , only : e_ice,denh2o, denice
    use clm_varpar           , only : nlevsoi, max_patch_per_col, nlevgrnd
    use clm_time_manager     , only : get_step_size
    use column_varcon        , only : icol_roof, icol_road_imperv
    use TridiagonalMod       , only : Tridiagonal
    use abortutils           , only : endrun     
    use SoilStateType        , only : soilstate_type
    use SoilHydrologyType    , only : soilhydrology_type
    use TemperatureType      , only : temperature_type
    use WaterFluxType        , only : waterflux_type
    use WaterStateType       , only : waterstate_type
    use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
    use PatchType            , only : pft
    use ColumnType           , only : col
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type)       , intent(in)    :: bounds               ! bounds
    integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                 , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                 , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    type(soilhydrology_type), intent(inout) :: soilhydrology_vars
    type(soilstate_type)    , intent(inout) :: soilstate_vars
    type(waterflux_type)    , intent(inout) :: waterflux_vars
    type(waterstate_type)   , intent(inout) :: waterstate_vars
    type(temperature_type)  , intent(in)    :: temperature_vars
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    integer  :: p,c,fc,j                                     ! do loop indices
    integer  :: jtop(bounds%begc:bounds%endc)                ! top level at each column
    real(r8) :: dtime                                        ! land model time step (sec)
    real(r8) :: hk(bounds%begc:bounds%endc,1:nlevsoi)        ! hydraulic conductivity [mm h2o/s]
    real(r8) :: dhkdw(bounds%begc:bounds%endc,1:nlevsoi)     ! d(hk)/d(vol_liq)
    real(r8) :: amx(bounds%begc:bounds%endc,1:nlevsoi+1)     ! "a" left off diagonal of tridiagonal matrix
    real(r8) :: bmx(bounds%begc:bounds%endc,1:nlevsoi+1)     ! "b" diagonal column for tridiagonal matrix
    real(r8) :: cmx(bounds%begc:bounds%endc,1:nlevsoi+1)     ! "c" right off diagonal tridiagonal matrix
    real(r8) :: rmx(bounds%begc:bounds%endc,1:nlevsoi+1)     ! "r" forcing term of tridiagonal matrix
    real(r8) :: zmm(bounds%begc:bounds%endc,1:nlevsoi+1)     ! layer depth [mm]
    real(r8) :: dzmm(bounds%begc:bounds%endc,1:nlevsoi+1)    ! layer thickness [mm]
    real(r8) :: den                                          ! used in calculating qin, qout
    real(r8) :: dqidw0(bounds%begc:bounds%endc,1:nlevsoi+1)  ! d(qin)/d(vol_liq(i-1))
    real(r8) :: dqidw1(bounds%begc:bounds%endc,1:nlevsoi+1)  ! d(qin)/d(vol_liq(i))
    real(r8) :: dqodw1(bounds%begc:bounds%endc,1:nlevsoi+1)  ! d(qout)/d(vol_liq(i))
    real(r8) :: dqodw2(bounds%begc:bounds%endc,1:nlevsoi+1)  ! d(qout)/d(vol_liq(i+1))
    real(r8) :: dsmpdw(bounds%begc:bounds%endc,1:nlevsoi+1)  ! d(smp)/d(vol_liq)
    real(r8) :: num                                          ! used in calculating qin, qout
    real(r8) :: qin(bounds%begc:bounds%endc,1:nlevsoi+1)     ! flux of water into soil layer [mm h2o/s]
    real(r8) :: qout(bounds%begc:bounds%endc,1:nlevsoi+1)    ! flux of water out of soil layer [mm h2o/s]
    real(r8) :: s_node                                       ! soil wetness
    real(r8) :: s1                                           ! "s" at interface of layer
    real(r8) :: s2                                           ! k*s**(2b+2)
    real(r8) :: smp(bounds%begc:bounds%endc,1:nlevsoi)       ! soil matrix potential [mm]
    real(r8) :: sdamp                                        ! extrapolates soiwat dependence of evaporation
    integer  :: pi                                           ! pft index
    real(r8) :: temp(bounds%begc:bounds%endc)                ! accumulator for rootr weighting
    integer  :: jwt(bounds%begc:bounds%endc)                 ! index of the soil layer right above the water table (-)
    real(r8) :: smp1,dsmpdw1,wh,wh_zwt,ka
    real(r8) :: dwat2(bounds%begc:bounds%endc,1:nlevsoi+1)
    real(r8) :: dzq                                          ! used in calculating qin, qout (difference in equilbirium matric potential)
    real(r8) :: zimm(bounds%begc:bounds%endc,0:nlevsoi)      ! layer interface depth [mm]
    real(r8) :: zq(bounds%begc:bounds%endc,1:nlevsoi+1)      ! equilibrium matric potential for each layer [mm]
    real(r8) :: vol_eq(bounds%begc:bounds%endc,1:nlevsoi+1)  ! equilibrium volumetric water content
    real(r8) :: tempi                                        ! temp variable for calculating vol_eq
    real(r8) :: temp0                                        ! temp variable for calculating vol_eq
    real(r8) :: voleq1                                       ! temp variable for calculating vol_eq
    real(r8) :: zwtmm(bounds%begc:bounds%endc)               ! water table depth [mm]
    real(r8) :: imped(bounds%begc:bounds%endc,1:nlevsoi)             
    real(r8) :: vol_ice(bounds%begc:bounds%endc,1:nlevsoi)
    real(r8) :: z_mid
    real(r8) :: vwc_zwt(bounds%begc:bounds%endc)
    real(r8) :: vwc_liq(bounds%begc:bounds%endc,1:nlevsoi+1) ! liquid volumetric water content
    real(r8) :: smp_grad(bounds%begc:bounds%endc,1:nlevsoi+1)
    real(r8) :: dsmpds                                       !temporary variable
    real(r8) :: dhkds                                        !temporary variable
    real(r8) :: hktmp                                        !temporary variable
    !-----------------------------------------------------------------------

    associate(& 
         z                 =>    col%z                              , & ! Input:  [real(r8) (:,:) ]  layer depth (m)                                 
         zi                =>    col%zi                             , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)           
         dz                =>    col%dz                             , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)                             

         origflag          =>    soilhydrology_vars%origflag        , & ! Input:  constant
         qcharge           =>    soilhydrology_vars%qcharge_col     , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)                      
         zwt               =>    soilhydrology_vars%zwt_col         , & ! Input:  [real(r8) (:)   ]  water table depth (m)                             
         fracice           =>    soilhydrology_vars%fracice_col     , & ! Input:  [real(r8) (:,:) ]  fractional impermeability (-)                   
         icefrac           =>    soilhydrology_vars%icefrac_col     , & ! Input:  [real(r8) (:,:) ]  fraction of ice                                 
         hkdepth           =>    soilhydrology_vars%hkdepth_col     , & ! Input:  [real(r8) (:)   ]  decay factor (m)                                  

         smpmin            =>    soilstate_vars%smpmin_col          , & ! Input:  [real(r8) (:)   ]  restriction for min of soil potential (mm)        
         watsat            =>    soilstate_vars%watsat_col          , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)  
         hksat             =>    soilstate_vars%hksat_col           , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity at saturation (mm H2O /s)
         bsw               =>    soilstate_vars%bsw_col             , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"                        
         sucsat            =>    soilstate_vars%sucsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)                       
         eff_porosity      =>    soilstate_vars%eff_porosity_col    , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice         
         rootr_col         =>    soilstate_vars%rootr_col           , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer  
         smp_l             =>    soilstate_vars%smp_l_col           , & ! Input:  [real(r8) (:,:) ]  soil matrix potential [mm]                      
         hk_l              =>    soilstate_vars%hk_l_col            , & ! Input:  [real(r8) (:,:) ]  hydraulic conductivity (mm/s)                   
         rootr_pft         =>    soilstate_vars%rootr_patch         , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer  

         h2osoi_ice        =>    waterstate_vars%h2osoi_ice_col     , & ! Input:  [real(r8) (:,:) ]  ice water (kg/m2)
         h2osoi_liq        =>    waterstate_vars%h2osoi_liq_col     , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_vol        =>    waterstate_vars%h2osoi_vol_col     , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]

         qflx_deficit      =>    waterflux_vars%qflx_deficit_col    , & ! Input:  [real(r8) (:)   ]  water deficit to keep non-negative liquid water content
         qflx_infl         =>    waterflux_vars%qflx_infl_col       , & ! Input:  [real(r8) (:)   ]  infiltration (mm H2O /s)                          
         qflx_tran_veg_col =>    waterflux_vars%qflx_tran_veg_col   , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)  
         qflx_tran_veg_pft =>    waterflux_vars%qflx_tran_veg_patch , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)  
         qflx_rootsoi      =>    waterflux_vars%qflx_rootsoi_col    , & ! Output: [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/s) (+ = to atm)

         t_soisno          =>    temperature_vars%t_soisno_col        & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)                       
         )

      ! Get time step
      
      dtime = get_step_size()

      ! Because the depths in this routine are in mm, use local
      ! variable arrays instead of pointers

      do j = 1, nlevsoi
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            zmm(c,j) = z(c,j)*1.e3_r8
            dzmm(c,j) = dz(c,j)*1.e3_r8
            zimm(c,j) = zi(c,j)*1.e3_r8

            ! calculate icefrac up here
            vol_ice(c,j) = min(watsat(c,j), h2osoi_ice(c,j)/(dz(c,j)*denice))
            icefrac(c,j) = min(1._r8,vol_ice(c,j)/watsat(c,j))
            vwc_liq(c,j) = max(h2osoi_liq(c,j),1.0e-6_r8)/(dz(c,j)*denh2o)
         end do
      end do

      do fc = 1, num_hydrologyc 
         c = filter_hydrologyc(fc)
         zimm(c,0) = 0.0_r8
         zwtmm(c)  = zwt(c)*1.e3_r8
      end do

      ! First step is to calculate the column-level effective rooting
      ! fraction in each soil layer. This is done outside the usual
      ! PFT-to-column averaging routines because it is not a simple
      ! weighted average of the PFT level rootr arrays. Instead, the
      ! weighting depends on both the per-unit-area transpiration
      ! of the PFT and the PFTs area relative to all PFTs.

      temp(bounds%begc : bounds%endc) = 0._r8

      do j = 1, nlevsoi
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            rootr_col(c,j) = 0._r8
         end do
      end do

      do pi = 1,max_patch_per_col
         do j = 1,nlevsoi
            do fc = 1, num_hydrologyc
               c = filter_hydrologyc(fc)
               if (pi <= col%npfts(c)) then
                  p = col%pfti(c) + pi - 1
                  if (pft%active(p)) then
                     rootr_col(c,j) = rootr_col(c,j) + rootr_pft(p,j) * qflx_tran_veg_pft(p) * pft%wtcol(p)
                  end if
               end if
            end do
         end do
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if (pi <= col%npfts(c)) then
               p = col%pfti(c) + pi - 1
               if (pft%active(p)) then
                  temp(c) = temp(c) + qflx_tran_veg_pft(p) * pft%wtcol(p)
               end if
            end if
         end do
      end do

      do j = 1, nlevsoi
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if (temp(c) /= 0._r8) then
               rootr_col(c,j) = rootr_col(c,j)/temp(c)
            end if
         end do
      end do

      !compute jwt index
      ! The layer index of the first unsaturated layer, i.e., the layer right above
      ! the water table

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         jwt(c) = nlevsoi
         ! allow jwt to equal zero when zwt is in top layer
         do j = 1,nlevsoi
            if(zwt(c) <= zi(c,j)) then
               jwt(c) = j-1
               exit
            end if
         enddo

         ! compute vwc at water table depth (mainly for case when t < tfrz)
         !     this will only be used when zwt is below the soil column
         vwc_zwt(c) = watsat(c,nlevsoi)
         if(t_soisno(c,jwt(c)+1) < tfrz) then
            vwc_zwt(c) = vwc_liq(c,nlevsoi)
            do j = nlevsoi,nlevgrnd
               if(zwt(c) <= zi(c,j)) then
                  smp1 = hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
                  !smp1 = max(0._r8,smp1)
                  smp1 = max(sucsat(c,nlevsoi),smp1)
                  vwc_zwt(c) = watsat(c,nlevsoi)*(smp1/sucsat(c,nlevsoi))**(-1._r8/bsw(c,nlevsoi))
                  ! for temperatures close to tfrz, limit vwc to total water content 
                  vwc_zwt(c) = min(vwc_zwt(c), 0.5*(watsat(c,nlevsoi) + h2osoi_vol(c,nlevsoi)) )
                  exit
               endif
            enddo
         endif
      end do

      ! calculate the equilibrium water content based on the water table depth

      do j=1,nlevsoi 
         do fc=1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if ((zwtmm(c) <= zimm(c,j-1))) then 
               vol_eq(c,j) = watsat(c,j)

               ! use the weighted average from the saturated part (depth > wtd) and the equilibrium solution for the
               ! rest of the layer, the equilibrium solution is based on Clapp-Hornberg parameterization
               ! and no extension to full range swrc is needed

            else if ((zwtmm(c) .lt. zimm(c,j)) .and. (zwtmm(c) .gt. zimm(c,j-1))) then
               tempi = 1.0_r8
               temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
               voleq1 = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zwtmm(c)-zimm(c,j-1))*(tempi-temp0)
               vol_eq(c,j) = (voleq1*(zwtmm(c)-zimm(c,j-1)) + watsat(c,j)*(zimm(c,j)-zwtmm(c)))/(zimm(c,j)-zimm(c,j-1))
               vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
               vol_eq(c,j) = max(vol_eq(c,j),0.0_r8)
            else
               tempi = (((sucsat(c,j)+zwtmm(c)-zimm(c,j))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
               temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j-1))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
               vol_eq(c,j) = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zimm(c,j)-zimm(c,j-1))*(tempi-temp0)
               vol_eq(c,j) = max(vol_eq(c,j),0.0_r8)
               vol_eq(c,j) = min(watsat(c,j),vol_eq(c,j))
            endif
            zq(c,j) = -sucsat(c,j)*(max(vol_eq(c,j)/watsat(c,j),0.01_r8))**(-bsw(c,j))
            zq(c,j) = max(smpmin(c), zq(c,j))
         end do
      end do

      ! If water table is below soil column calculate zq for the 11th layer
      j = nlevsoi
      do fc=1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if(jwt(c) == nlevsoi) then 
            tempi = 1._r8
            temp0 = (((sucsat(c,j)+zwtmm(c)-zimm(c,j))/sucsat(c,j)))**(1._r8-1._r8/bsw(c,j))
            vol_eq(c,j+1) = -sucsat(c,j)*watsat(c,j)/(1._r8-1._r8/bsw(c,j))/(zwtmm(c)-zimm(c,j))*(tempi-temp0)
            vol_eq(c,j+1) = max(vol_eq(c,j+1),0.0_r8)
            vol_eq(c,j+1) = min(watsat(c,j),vol_eq(c,j+1))
            zq(c,j+1) = -sucsat(c,j)*(max(vol_eq(c,j+1)/watsat(c,j),0.01_r8))**(-bsw(c,j))
            zq(c,j+1) = max(smpmin(c), zq(c,j+1))
         end if
      end do

      ! Hydraulic conductivity and soil matric potential and their derivatives

      sdamp = 0._r8
      do j = 1, nlevsoi
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            ! compute hydraulic conductivity based on liquid water content only

            if (origflag == 1) then
               s1 = 0.5_r8*(h2osoi_vol(c,j) + h2osoi_vol(c,min(nlevsoi, j+1))) / &
                    (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
            else
               s1 = 0.5_r8*(vwc_liq(c,j) + vwc_liq(c,min(nlevsoi, j+1))) / &
                    (0.5_r8*(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))
            endif
            s1 = min(1._r8, s1)
            s2 = hksat(c,j)*s1**(2._r8*bsw(c,j)+2._r8)

            ! replace fracice with impedance factor, as in zhao 97,99
            if (origflag == 1) then
               imped(c,j)=(1._r8-0.5_r8*(fracice(c,j)+fracice(c,min(nlevsoi, j+1))))
            else
               imped(c,j)=10._r8**(-e_ice*(0.5_r8*(icefrac(c,j)+icefrac(c,min(nlevsoi, j+1)))))
            endif
            hk(c,j) = imped(c,j)*s1*s2
            dhkdw(c,j) = imped(c,j)*(2._r8*bsw(c,j)+3._r8)*s2* &
                 (1._r8/(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))

            !compute un-restricted hydraulic conductivity
            !call soil_water_retention_curve%soil_hk(hksat(c,j), imped(c,j), s1, bsw(c,j), hktmp, dhkds)
            !if(hktmp/=hk(c,j))write(10,*)'diff',hktmp,hk(c,j)
            !    call endrun('bad in hk')
            !endif    
            !apply ice impedance
            !hk(c,j) = imped(c,j)*hk(c,j)          
            !dhkdw(c,j) = imped(c,j) * dhkds * (1._r8/(watsat(c,j)+watsat(c,min(nlevsoi, j+1))))


            ! compute matric potential and derivative based on liquid water content only
            if (origflag == 1) then
               s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
            else
               s_node = max(vwc_liq(c,j)/watsat(c,j), 0.01_r8)
            endif
            s_node = min(1.0_r8, s_node)

            !call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp(c,j), dsmpds)

            smp(c,j) = -sucsat(c,j)*s_node**(-bsw(c,j))
            smp(c,j) = max(smpmin(c), smp(c,j))
            !do not turn on the line below, which will cause bit to bit error, jyt, 2014 Mar 6
            !dsmpdw(c,j) = dsmpds/watsat(c,j)

            if (origflag == 1) then             
               dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/(s_node*watsat(c,j))
            else
               dsmpdw(c,j) = -bsw(c,j)*smp(c,j)/vwc_liq(c,j)
            endif

            smp_l(c,j) = smp(c,j)
            hk_l(c,j) = hk(c,j)

         end do
      end do

      ! aquifer (11th) layer
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         zmm(c,nlevsoi+1) = 0.5*(1.e3_r8*zwt(c) + zmm(c,nlevsoi))
         if(jwt(c) < nlevsoi) then
            dzmm(c,nlevsoi+1) = dzmm(c,nlevsoi)
         else
            dzmm(c,nlevsoi+1) = (1.e3_r8*zwt(c) - zmm(c,nlevsoi))
         end if
      end do

      ! Set up r, a, b, and c vectors for tridiagonal solution

      ! Node j=1 (top)

      j = 1
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         qin(c,j)    = qflx_infl(c)
         den    = (zmm(c,j+1)-zmm(c,j))
         dzq    = (zq(c,j+1)-zq(c,j))
         num    = (smp(c,j+1)-smp(c,j)) - dzq
         qout(c,j)   = -hk(c,j)*num/den
         dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
         dqodw2(c,j) = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
         rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c) * rootr_col(c,j)
         amx(c,j) =  0._r8
         bmx(c,j) =  dzmm(c,j)*(sdamp+1._r8/dtime) + dqodw1(c,j)
         cmx(c,j) =  dqodw2(c,j)
      end do

      ! Nodes j=2 to j=nlevsoi-1

      do j = 2, nlevsoi - 1
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            den    = (zmm(c,j) - zmm(c,j-1))
            dzq    = (zq(c,j)-zq(c,j-1))
            num    = (smp(c,j)-smp(c,j-1)) - dzq
            qin(c,j)    = -hk(c,j-1)*num/den
            dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
            dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
            den    = (zmm(c,j+1)-zmm(c,j))
            dzq    = (zq(c,j+1)-zq(c,j))
            num    = (smp(c,j+1)-smp(c,j)) - dzq
            qout(c,j)   = -hk(c,j)*num/den
            dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
            dqodw2(c,j) = -( hk(c,j)*dsmpdw(c,j+1) + num*dhkdw(c,j))/den
            rmx(c,j)    =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c)*rootr_col(c,j)
            amx(c,j)    = -dqidw0(c,j)
            bmx(c,j)    =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
            cmx(c,j)    =  dqodw2(c,j)
         end do
      end do

      ! Node j=nlevsoi (bottom)

      j = nlevsoi
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         if(j > jwt(c)) then !water table is in soil column
            den    = (zmm(c,j) - zmm(c,j-1))
            dzq    = (zq(c,j)-zq(c,j-1))
            num    = (smp(c,j)-smp(c,j-1)) - dzq
            qin(c,j)    = -hk(c,j-1)*num/den
            dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
            dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
            qout(c,j)   =  0._r8
            dqodw1(c,j) =  0._r8
            rmx(c,j)    =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c)*rootr_col(c,j)
            amx(c,j)    = -dqidw0(c,j)
            bmx(c,j)    =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
            cmx(c,j)    =  0._r8

            ! next set up aquifer layer; hydrologically inactive
            rmx(c,j+1) = 0._r8
            amx(c,j+1) = 0._r8
            bmx(c,j+1) = dzmm(c,j+1)/dtime
            cmx(c,j+1) = 0._r8
         else ! water table is below soil column

            ! compute aquifer soil moisture as average of layer 10 and saturation
            if(origflag == 1) then
               s_node = max(0.5*(1.0_r8+h2osoi_vol(c,j)/watsat(c,j)), 0.01_r8)
            else
               s_node = max(0.5*((vwc_zwt(c)+vwc_liq(c,j))/watsat(c,j)), 0.01_r8)
            endif
            s_node = min(1.0_r8, s_node)

            ! compute smp for aquifer layer
            !call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp1, dsmpds)
            smp1 = -sucsat(c,j)*s_node**(-bsw(c,j))
            smp1 = max(smpmin(c), smp1)

            ! compute dsmpdw for aquifer layer
            !dsmpdw1 = dsmpds/watsat(c,j)
            dsmpdw1 = -bsw(c,j)*smp1/(s_node*watsat(c,j))

            ! first set up bottom layer of soil column
            den    = (zmm(c,j) - zmm(c,j-1))
            dzq    = (zq(c,j)-zq(c,j-1))
            num    = (smp(c,j)-smp(c,j-1)) - dzq
            qin(c,j)    = -hk(c,j-1)*num/den
            dqidw0(c,j) = -(-hk(c,j-1)*dsmpdw(c,j-1) + num*dhkdw(c,j-1))/den
            dqidw1(c,j) = -( hk(c,j-1)*dsmpdw(c,j)   + num*dhkdw(c,j-1))/den
            den    = (zmm(c,j+1)-zmm(c,j))
            dzq    = (zq(c,j+1)-zq(c,j))
            num    = (smp1-smp(c,j)) - dzq
            qout(c,j)   = -hk(c,j)*num/den
            dqodw1(c,j) = -(-hk(c,j)*dsmpdw(c,j)   + num*dhkdw(c,j))/den
            dqodw2(c,j) = -( hk(c,j)*dsmpdw1 + num*dhkdw(c,j))/den

            rmx(c,j) =  qin(c,j) - qout(c,j) - qflx_tran_veg_col(c)*rootr_col(c,j)
            amx(c,j) = -dqidw0(c,j)
            bmx(c,j) =  dzmm(c,j)/dtime - dqidw1(c,j) + dqodw1(c,j)
            cmx(c,j) =  dqodw2(c,j)

            ! next set up aquifer layer; den/num unchanged, qin=qout
            qin(c,j+1)    = qout(c,j)
            dqidw0(c,j+1) = -(-hk(c,j)*dsmpdw(c,j) + num*dhkdw(c,j))/den
            dqidw1(c,j+1) = -( hk(c,j)*dsmpdw1   + num*dhkdw(c,j))/den
            qout(c,j+1)   =  0._r8  ! zero-flow bottom boundary condition
            dqodw1(c,j+1) =  0._r8  ! zero-flow bottom boundary condition
            rmx(c,j+1) =  qin(c,j+1) - qout(c,j+1)
            amx(c,j+1) = -dqidw0(c,j+1)
            bmx(c,j+1) =  dzmm(c,j+1)/dtime - dqidw1(c,j+1) + dqodw1(c,j+1)
            cmx(c,j+1) =  0._r8
         endif
      end do

      ! Solve for dwat

      jtop(bounds%begc : bounds%endc) = 1
      call Tridiagonal(bounds, 1, nlevsoi+1, &
           jtop(bounds%begc:bounds%endc), &
           num_hydrologyc, filter_hydrologyc, &
           amx(bounds%begc:bounds%endc, :), &
           bmx(bounds%begc:bounds%endc, :), &
           cmx(bounds%begc:bounds%endc, :), &
           rmx(bounds%begc:bounds%endc, :), &
           dwat2(bounds%begc:bounds%endc, :) )

      ! Renew the mass of liquid water
      ! also compute qcharge from dwat in aquifer layer
      ! update in drainage for case jwt < nlevsoi

      do fc = 1,num_hydrologyc
         c = filter_hydrologyc(fc)
         do j = 1, nlevsoi
            h2osoi_liq(c,j) = h2osoi_liq(c,j) + dwat2(c,j)*dzmm(c,j)
         end do

         ! calculate qcharge for case jwt < nlevsoi
         if(jwt(c) < nlevsoi) then
            wh_zwt = 0._r8   !since wh_zwt = -sucsat - zq_zwt, where zq_zwt = -sucsat

            ! Recharge rate qcharge to groundwater (positive to aquifer)
            s_node = max(h2osoi_vol(c,jwt(c)+1)/watsat(c,jwt(c)+1), 0.01_r8)
            s1 = min(1._r8, s_node)

            !scs: this is the expression for unsaturated hk
            ka = imped(c,jwt(c)+1)*hksat(c,jwt(c)+1) &
                 *s1**(2._r8*bsw(c,jwt(c)+1)+3._r8)

            !compute unsaturated hk, this shall be tested later, because it
            !is not bit for bit
            !call soil_water_retention_curve%soil_hk(hksat(c,jwt(c)+1), s1, bsw(c,jwt(c)+1), ka)
            !apply ice impedance
            !ka = imped(c,jwt(c)+1) * ka 
            ! Recharge rate qcharge to groundwater (positive to aquifer)
            smp1 = max(smpmin(c), smp(c,max(1,jwt(c))))
            wh      = smp1 - zq(c,max(1,jwt(c)))

            !scs: original formulation
            if(jwt(c) == 0) then
               qcharge(c) = -ka * (wh_zwt-wh)  /((zwt(c)+1.e-3)*1000._r8)
            else
               !             qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._r8)
               !scs: 1/2, assuming flux is at zwt interface, saturation deeper than zwt
               qcharge(c) = -ka * (wh_zwt-wh)/((zwt(c)-z(c,jwt(c)))*1000._r8*2.0)
            endif

            ! To limit qcharge  (for the first several timesteps)
            qcharge(c) = max(-10.0_r8/dtime,qcharge(c))
            qcharge(c) = min( 10.0_r8/dtime,qcharge(c))
         else
            ! if water table is below soil column, compute qcharge from dwat2(11)
            qcharge(c) = dwat2(c,nlevsoi+1)*dzmm(c,nlevsoi+1)/dtime
         endif
      end do

      ! compute the water deficit and reset negative liquid water content
      !  Jinyun Tang
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         qflx_deficit(c) = 0._r8
         do j = 1, nlevsoi
            if(h2osoi_liq(c,j)<0._r8)then
               qflx_deficit(c) = qflx_deficit(c) - h2osoi_liq(c,j)
            endif
         enddo
      enddo

      do j = 1, nlevsoi
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            qflx_rootsoi(c,j) = qflx_tran_veg_col(c) * rootr_col(c,j) * 1.e-3_r8       ![m H2O/s]
         enddo
      enddo

    end associate 

  end subroutine soilwater_zengdecker2009

  !-----------------------------------------------------------------------
  subroutine init_vsfm_condition_ids()
    !
    !DESCRIPTION
    !  Determines the IDs of various source-sink conditions in VSFM
    !
    use clm_varctl, only : use_vsfm
#ifdef USE_PETSC_LIB
    use MultiPhysicsProbVSFM             , only : vsfm_mpp
    use MultiPhysicsProbConstants        , only : COND_SS
    use MultiPhysicsProbConstants        , only : COND_NULL
    use clm_varctl                       , only : iulog
    use abortutils                       , only : endrun
    use shr_log_mod                      , only : errMsg => shr_log_errMsg
#endif
    ! !ARGUMENTS:
    implicit none
    integer :: ier ! error status
    character (len=256), pointer :: cond_names(:)
    integer                      :: num_conds
    integer                      :: nn
    integer                      :: kk
    character (len=256)          :: cond_name
    !------------------------------------------------------------------------------

#ifdef USE_PETSC_LIB

    ! Get the number of
    call vsfm_mpp%sysofeqns%GetConditionNames(COND_SS, COND_NULL, num_conds, cond_names)

    if (num_conds /= 6) then
      write(iulog,*)'In init_vsfm_condition_ids: Source-sink conditions /= 6'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    do nn = 1, num_conds
       select case(trim(cond_names(nn)))
       case ("Infiltration_Flux")
         vsfm_cond_id_for_infil = nn
       case ("Evapotranspiration_Flux")
         vsfm_cond_id_for_et = nn
       case ("Dew_Flux")
         vsfm_cond_id_for_dew = nn
       case ("Drainage_Flux")
         vsfm_cond_id_for_drainage = nn
       case ("Snow_Disappearance_Flux")
         vsfm_cond_id_for_snow = nn
       case ("Sublimation_Flux")
         vsfm_cond_id_for_sublimation = nn
       case default
         write(iulog,*)'In init_vsfm_condition_ids: Unknown flux.'
         call endrun(msg=errMsg(__FILE__, __LINE__))
       end select
    enddo

    if (vsfm_cond_id_for_infil == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_infil not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_cond_id_for_et == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_et not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_cond_id_for_dew == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_dew not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_cond_id_for_drainage == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_drainage not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_cond_id_for_snow == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_snow not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    if (vsfm_cond_id_for_sublimation == -1) then
      write(iulog,*)'In init_vsfm_condition_ids: vsfm_cond_id_for_sublimation not defined.'
      call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    deallocate(cond_names)

#endif

  end subroutine init_vsfm_condition_ids

  !-----------------------------------------------------------------------
  subroutine soilwater_vsfm(bounds, num_hydrologyc, filter_hydrologyc, &
       num_urbanc, filter_urbanc, soilhydrology_vars, soilstate_vars, &
       waterflux_vars, waterstate_vars, temperature_vars)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use decompMod                 , only : bounds_type
    use clm_varcon                , only : denh2o
    use clm_varpar                , only : nlevsoi, max_patch_per_col, nlevgrnd
    use clm_time_manager          , only : get_step_size, get_nstep
    use abortutils                , only : endrun
    use SoilStateType             , only : soilstate_type
    use SoilHydrologyType         , only : soilhydrology_type
    use TemperatureType           , only : temperature_type
    use WaterFluxType             , only : waterflux_type
    use WaterStateType            , only : waterstate_type
    use PatchType                 , only : pft
    use ColumnType                , only : col
    use clm_varcon                , only : watmin
    use LandunitType              , only : lun
    use landunit_varcon           , only : istsoil, istcrop
    use clm_varctl                , only : iulog
    use shr_log_mod               , only : errMsg => shr_log_errMsg
#ifdef USE_PETSC_LIB
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_FRAC_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
#endif
    !
    ! !ARGUMENTS:
    implicit none
#ifdef USE_PETSC_LIB
#include "finclude/petscsys.h"
#endif
    !
    type(bounds_type)       , intent(in)    :: bounds               ! bounds
    integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    integer                 , intent(in)    :: num_urbanc           ! number of column urban points in column filter
    integer                 , intent(in)    :: filter_urbanc(:)     ! column filter for urban points
    type(soilhydrology_type), intent(inout) :: soilhydrology_vars
    type(soilstate_type)    , intent(inout) :: soilstate_vars
    type(waterflux_type)    , intent(inout) :: waterflux_vars
    type(waterstate_type)   , intent(inout) :: waterstate_vars
    type(temperature_type)  , intent(in)    :: temperature_vars
    !
    ! !LOCAL VARIABLES:
    integer              :: p,c,fc,j                                                         ! do loop indices
    real(r8)             :: dtime                                                            ! land model time step (sec)
    real(r8)             :: temp(bounds%begc:bounds%endc)                                    ! accumulator for rootr weighting
    integer              :: pi                                                               ! pft index
    real(r8)             :: dzsum                                                            ! summation of dzmm of layers below water table (mm)
    real(r8)             :: frac_ice                    (bounds%begc:bounds%endc,1:nlevgrnd) ! fraction of ice
    real(r8)             :: total_mass_flux_col         (bounds%begc:bounds%endc)            ! Sum of all source-sinks conditions for VSFM solver at column level
    real(r8)             :: total_mass_flux_et_col      (bounds%begc:bounds%endc)            ! ET sink for VSFM solver at column level
    real(r8)             :: total_mass_flux_infl_col    (bounds%begc:bounds%endc)            ! Infiltration source for VSFM solver at column level
    real(r8)             :: total_mass_flux_dew_col     (bounds%begc:bounds%endc)            ! Dew source for VSFM solver at column level
    real(r8)             :: total_mass_flux_drain_col   (bounds%begc:bounds%endc)            ! Drainage sink for VSFM solver at column level
    real(r8)             :: total_mass_flux_snowlyr_col (bounds%begc:bounds%endc)            ! Flux due to disappearance of snow for VSFM solver at column level
    real(r8)             :: total_mass_flux_sub_col     (bounds%begc:bounds%endc)            ! Sublimation sink for VSFM solver at column level
    real(r8)             :: vsfm_mass_prev_col          (bounds%begc:bounds%endc,1:nlevgrnd) ! Mass of water before a VSFM solve
    real(r8)             :: vsfm_dmass_col              (bounds%begc:bounds%endc)            ! Change in mass of water after a VSFM solve
    real(r8)             :: mass_beg_col                (bounds%begc:bounds%endc)            ! Total mass before a VSFM solve
    real(r8)             :: mass_end_col                (bounds%begc:bounds%endc)            ! Total mass after a VSFM solve
    integer              :: ier                                                              ! error status

#ifdef USE_PETSC_LIB
    PetscInt             :: jwt                                                              ! index of first unsaturated soil layer
    PetscInt             :: idx                                                              ! 1D index for (c,j)
    PetscInt             :: soe_auxvar_id                                                    ! Index of system-of-equation's (SoE's) auxvar
    PetscReal            :: flux_unit_conversion                                             ! [mm/s] ---> [kg/s]
    PetscReal            :: area                                                             ! [m^2]
    PetscReal            :: z_up, z_dn                                                       ! [m]
    PetscReal            :: qflx_drain_layer                                                 ! Drainage flux from a soil layer (mm H2O/s)
    PetscReal            :: qflx_drain_tot                                                   ! Cummulative drainage flux from soil layers within a column (mm H2O/s)
    PetscErrorCode       :: ierr                                                             ! PETSc return error code

    PetscBool            :: converged                                                        ! Did VSFM solver converge to a solution with given PETSc SNES tolerances
    PetscReal            :: atol_default                                                     ! Default SNES absolute convergance tolerance
    PetscReal            :: rtol_default                                                     ! Default SNES relative convergance tolerance
    PetscReal            :: stol_default                                                     ! Default SNES solution convergance tolerance
    PetscInt             :: max_it_default                                                   ! Default SNES maximum number of iteration
    PetscInt             :: max_f_default                                                    ! Default SNES maximum number of function evaluation
    PetscReal            :: stol                                                             ! solution convergance tolerance
    PetscReal,parameter  :: stol_alternate = 1.d-10                                          ! Alternate solution convergance tolerance

    PetscReal            :: mass_beg                                                         ! Sum of mass of water for all active soil columns before VSFM is called
    PetscReal            :: mass_end                                                         ! Sum of mass of water for all active soil columns after VSFM is called
    PetscReal            :: total_mass_flux_et                                               ! Sum of mass ET mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_infl                                             ! Sum of mass infiltration mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_dew                                              ! Sum of mass dew mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_drain                                            ! Sum of mass drainage mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_snowlyr                                          ! Sum of mass snow layer disappearance mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_sub                                              ! Sum of mass sublimation mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux                                                  ! Sum of mass ALL mass flux of water for all active soil columns
    PetscInt             :: iter_count                                                       ! How many times VSFM solver is called

    PetscInt, parameter  :: max_iter_count = 10                                              ! Maximum number of times VSFM can be called
    PetscInt             :: diverged_count                                                   ! Number of time VSFM solver diverged
    PetscInt             :: mass_bal_err_count                                               ! Number of time VSFM solver returns a solution that isn't within acceptable mass balance error threshold
    PetscReal            :: abs_mass_error_col                                               ! Maximum absolute error for any active soil column
    PetscReal, parameter :: max_abs_mass_error_col  = 1.e-5                                  ! Acceptable mass balance error
    PetscBool            :: successful_step                                                  ! Is the solution return by VSFM acceptable
#endif
    !-----------------------------------------------------------------------

    associate(& 
         z                         =>    col%z                                      , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         zi                        =>    col%zi                                     , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         dz                        =>    col%dz                                     , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         snl                       =>    col%snl                                    , & ! Input:  [integer  (:)   ]  minus number of snow layers

         qcharge                   =>    soilhydrology_vars%qcharge_col             , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)
         zwt                       =>    soilhydrology_vars%zwt_col                 , & ! Input:  [real(r8) (:)   ]  water table depth (m)

         watsat                    =>    soilstate_vars%watsat_col                  , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         rootr_col                 =>    soilstate_vars%rootr_col                   , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer
         rootr_pft                 =>    soilstate_vars%rootr_patch                 , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer
         smp_l                     =>    soilstate_vars%smp_l_col                   , & ! Output: [real(r8) (:,:) ]  soil matrix potential [mm]
         hk_l                      =>    soilstate_vars%hk_l_col                    , & ! Output: [real(r8) (:,:) ]  hydraulic conductivity (mm/s)

         h2osoi_ice                =>    waterstate_vars%h2osoi_ice_col             , & ! Input:  [real(r8) (:,:) ]  ice water (kg/m2)
         h2osoi_liq                =>    waterstate_vars%h2osoi_liq_col             , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_vol                =>    waterstate_vars%h2osoi_vol_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         frac_h2osfc               =>    waterstate_vars%frac_h2osfc_col            , & ! Input:  [real(r8) (:)   ]
         frac_sno                  =>    waterstate_vars%frac_sno_eff_col           , & ! Input:  [real(r8) (:)   ]  fraction of ground covered by snow (0 to 1)
         vsfm_fliq_col_1d          =>    waterstate_vars%vsfm_fliq_col_1d           , & ! Output: [real(r8) (:)   ]  1D fraction of liquid saturation for VSFM [-]
         vsfm_sat_col_1d           =>    waterstate_vars%vsfm_sat_col_1d            , & ! Output: [real(r8) (:)   ]  1D liquid saturation from VSFM [-]
         vsfm_mass_col_1d          =>    waterstate_vars%vsfm_mass_col_1d           , & ! Output: [real(r8) (:)   ]  1D liquid mass per unit area from VSFM [kg H2O/m^2]
         vsfm_smpl_col_1d          =>    waterstate_vars%vsfm_smpl_col_1d           , & ! Output: [real(r8) (:)   ]  1D soil matrix potential liquid from VSFM [m]
         vsfm_soilp_col_1d         =>    waterstate_vars%vsfm_soilp_col_1d          , & ! Output: [real(r8) (:)   ]  1D soil water pressure from VSFM [Pa]
         soilp_col                 =>    waterstate_vars%soilp_col                  , & ! Output: [real(r8) (:,:) ]  soil water pressure (Pa)

         qflx_deficit              =>    waterflux_vars%qflx_deficit_col            , & ! Input:  [real(r8) (:)   ]  water deficit to keep non-negative liquid water content
         qflx_infl                 =>    waterflux_vars%qflx_infl_col               , & ! Input:  [real(r8) (:)   ]  infiltration (mm H2O /s)
         qflx_tran_veg_col         =>    waterflux_vars%qflx_tran_veg_col           , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         qflx_tran_veg_pft         =>    waterflux_vars%qflx_tran_veg_patch         , & ! Input:  [real(r8) (:)   ]  vegetation transpiration (mm H2O/s) (+ = to atm)
         qflx_dew_snow             =>    waterflux_vars%qflx_dew_snow_col           , & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s) [+]
         qflx_dew_grnd             =>    waterflux_vars%qflx_dew_grnd_col           , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
         qflx_sub_snow             =>    waterflux_vars%qflx_sub_snow_col           , & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
         qflx_drain                =>    waterflux_vars%qflx_drain_col              , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)
         qflx_drain_perched        =>    waterflux_vars%qflx_drain_perched_col      , & ! Input:  [real(r8) (:)   ]  perched wt sub-surface runoff (mm H2O /s)
         mflx_infl_col_1d          =>    waterflux_vars%mflx_infl_col_1d            , & ! Input:  [real(r8) (:)   ]  infiltration source in top soil control volume (kg H2O /s)
         mflx_dew_col_1d           =>    waterflux_vars%mflx_dew_col_1d             , & ! Input:  [real(r8) (:)   ]  (liquid+snow) dew source in top soil control volume (kg H2O /s)
         mflx_et_col_1d            =>    waterflux_vars%mflx_et_col_1d              , & ! Input:  [real(r8) (:)   ]  evapotranspiration sink from all soil coontrol volumes (kg H2O /s) (+ = to atm)
         mflx_snowlyr_col_1d       =>    waterflux_vars%mflx_snowlyr_col_1d         , & ! Input:  [real(r8) (:)   ]  mass flux to top soil layer due to disappearance of snow (kg H2O /s)
         mflx_sub_snow_col_1d      =>    waterflux_vars%mflx_sub_snow_col_1d        , & ! Output: [real(r8) (:)   ]  mass flux from top soil layer due to sublimation of snow (kg H2O /s)
         mflx_snowlyr_col          =>    waterflux_vars%mflx_snowlyr_col            , & ! Input:  [real(r8) (:)   ]  mass flux to top soil layer due to disappearance of snow (kg H2O /s)
         mflx_drain_col_1d         =>    waterflux_vars%mflx_drain_col_1d           , & ! Input:  [real(r8) (:)   ]  drainage from groundwater and perched water table (kg H2O /s)
         mflx_drain_perched_col_1d =>    waterflux_vars%mflx_drain_perched_col_1d   , & ! Input:  [real(r8) (:)   ]  drainage from perched water table (kg H2O /s)
         mflx_neg_snow_col_1d      =>    waterflux_vars%mflx_neg_snow_col_1d        , & ! Input:  [real(r8) (:)   ]  mass flux from top soil layer due to negative water content in snow layers (kg H2O /s)

         t_soil_col_1d             =>    temperature_vars%t_soil_col_1d             , & ! Input:  [real(r8) (:)   ]  1D soil temperature (Kelvin)
         t_soisno                  =>    temperature_vars%t_soisno_col                & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)
         )

      ! Get time step

      dtime = get_step_size()


      ! First step is to calculate the column-level effective rooting
      ! fraction in each soil layer. This is done outside the usual
      ! PFT-to-column averaging routines because it is not a simple
      ! weighted average of the PFT level rootr arrays. Instead, the
      ! weighting depends on both the per-unit-area transpiration
      ! of the PFT and the PFTs area relative to all PFTs.

      temp(bounds%begc : bounds%endc) = 0._r8

      do j = 1, nlevsoi
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            rootr_col(c,j) = 0._r8
         end do
      end do

      do pi = 1,max_patch_per_col
         do j = 1,nlevsoi
            do fc = 1, num_hydrologyc
               c = filter_hydrologyc(fc)
               if (pi <= col%npfts(c)) then
                  p = col%pfti(c) + pi - 1
                  if (pft%active(p)) then
                     rootr_col(c,j) = rootr_col(c,j) + rootr_pft(p,j) * qflx_tran_veg_pft(p) * pft%wtcol(p)
                  end if
               end if
            end do
         end do
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if (pi <= col%npfts(c)) then
               p = col%pfti(c) + pi - 1
               if (pft%active(p)) then
                  temp(c) = temp(c) + qflx_tran_veg_pft(p) * pft%wtcol(p)
               end if
            end if
         end do
      end do

      do j = 1, nlevsoi
         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)
            if (temp(c) /= 0._r8) then
               rootr_col(c,j) = rootr_col(c,j)/temp(c)
            end if
         end do
      end do

#ifdef USE_PETSC_LIB

      area = 1.d0 ! [m^2]

      ! initialize
      mflx_et_col_1d(:)                = 0.d0
      mflx_infl_col_1d(:)              = 0.d0
      mflx_dew_col_1d(:)               = 0.d0
      mflx_drain_col_1d(:)             = 0.d0
      mflx_sub_snow_col_1d(:)          = 0.d0
      mflx_snowlyr_col_1d(:)           = 0.d0
      t_soil_col_1d(:)                 = 298.15d0

      mass_beg                         = 0.d0
      mass_end                         = 0.d0
      total_mass_flux                  = 0.d0
      total_mass_flux_et               = 0.d0
      total_mass_flux_infl             = 0.d0
      total_mass_flux_dew              = 0.d0
      total_mass_flux_drain            = 0.d0
      total_mass_flux_snowlyr          = 0.d0
      total_mass_flux_sub              = 0.d0

      mass_beg_col(:)                  = 0.d0
      mass_end_col(:)                  = 0.d0
      total_mass_flux_col(:)           = 0.d0
      total_mass_flux_et_col(:)        = 0.d0
      total_mass_flux_infl_col(:)      = 0.d0
      total_mass_flux_dew_col(:)       = 0.d0
      total_mass_flux_drain_col(:)     = 0.d0
      total_mass_flux_snowlyr_col(:)   = 0.d0
      total_mass_flux_sub_col(:)       = 0.d0

      vsfm_mass_prev_col(:,:)          = 0.d0
      vsfm_dmass_col(:)                = 0.d0

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         ! [mm/s] --> [kg/s]   [m^2] [kg/m^3]  [m/mm]
         flux_unit_conversion     = area * denh2o * 1.0d-3

         do j = 1, nlevsoi
            ! ET sink
            idx = (c-bounds%begc)*nlevgrnd + j
            mflx_et_col_1d(idx) = -qflx_tran_veg_col(c)*rootr_col(c,j)*flux_unit_conversion
         end do

         !do j = 1, nlevgrnd
         !   ! Temperature
         !   idx = (c-1)*nlevgrnd + j
         !   t_soil_col_1d(idx) = t_soisno(c,j)
         !end do

         ! Infiltration source term
         idx = c-bounds%begc+1
         mflx_infl_col_1d(idx) = qflx_infl(c)*flux_unit_conversion

         ! Dew and snow sublimation source/sink term
         if (snl(c) >= 0) then
            mflx_dew_col_1d(idx)       = (qflx_dew_snow(c) + qflx_dew_grnd(c))* &
                                         (1._r8 - frac_h2osfc(c))*              &
                                         flux_unit_conversion

            mflx_sub_snow_col_1d(idx)  = -qflx_sub_snow(c)*          &
                                          (1._r8 - frac_h2osfc(c))*  &
                                          flux_unit_conversion
         end if


         if (qflx_drain(c) > 0.d0) then

            ! Find soil layer just above water table
            jwt = nlevgrnd

            ! allow jwt to equal zero when zwt is in top layer
            do j = 1,nlevgrnd
               if (zwt(c) <= zi(c,j)) then
                  jwt = j-1
                  exit
               end if
            enddo

            ! Now ensure the soil layer index corresponding to the water table
            ! is greater than or equal to the first soil layer.
            jwt = max(jwt,1)

            dzsum = 0.d0
            do j = jwt, nlevgrnd
               dzsum = dzsum + dz(c,j)
            end do

            qflx_drain_tot = 0.d0
            do j = jwt, nlevgrnd
               qflx_drain_layer = qflx_drain(c) * dz(c,j)/dzsum

               ! if the amount of water being drained from a given layer
               ! exceeds the allowable water, limit the drainage
               if (qflx_drain_layer*dtime > (h2osoi_liq(c,j)-watmin)) then
                  qflx_drain_layer = (h2osoi_liq(c,j)-watmin)/dtime
               endif
               qflx_drain_tot = qflx_drain_tot + qflx_drain_layer

               idx = (c-bounds%begc)*nlevgrnd + j
               mflx_drain_col_1d(idx) = -qflx_drain_layer*flux_unit_conversion

           end do
           qflx_drain(c) = qflx_drain_tot

         endif

         ! The mass flux associated with disapperance of snow layer over the
         ! last time step.
         idx = c-bounds%begc+1
         mflx_snowlyr_col_1d(c-bounds%begc+1) = mflx_snowlyr_col(c) + &
                                                mflx_neg_snow_col_1d(c-bounds%begc+1)
         mflx_snowlyr_col(c) = 0._r8

      end do

      ! Set temperature
      soe_auxvar_id = 1;
      call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_INTERNAL ,      &
                                             VAR_TEMPERATURE ,      &
                                             soe_auxvar_id   ,      &
                                             t_soil_col_1d          &
                                            )

      ! Set Infiltration
      soe_auxvar_id = vsfm_cond_id_for_infil;
      call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS           ,  &
                                             VAR_BC_SS_CONDITION ,  &
                                             soe_auxvar_id       ,  &
                                             mflx_infl_col_1d       &
                                            )

      ! Set ET
      soe_auxvar_id = vsfm_cond_id_for_et;
      call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS           ,  &
                                             VAR_BC_SS_CONDITION ,  &
                                             soe_auxvar_id       ,  &
                                             mflx_et_col_1d         &
                                            )

      ! Set Dew
      soe_auxvar_id = vsfm_cond_id_for_dew;
      call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS           ,  &
                                             VAR_BC_SS_CONDITION ,  &
                                             soe_auxvar_id       ,  &
                                             mflx_dew_col_1d        &
                                            )

      ! Set Drainage sink
      soe_auxvar_id = vsfm_cond_id_for_drainage;
      mflx_drain_col_1d(:) = mflx_drain_col_1d(:) + mflx_drain_perched_col_1d(:)
      call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS           ,  &
                                             VAR_BC_SS_CONDITION ,  &
                                             soe_auxvar_id       ,  &
                                             mflx_drain_col_1d      &
                                            )

      ! Set mass flux associated with disappearance of snow layer
      ! from last time step
      soe_auxvar_id = vsfm_cond_id_for_snow;
      call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS           ,  &
                                             VAR_BC_SS_CONDITION ,  &
                                             soe_auxvar_id       ,  &
                                             mflx_snowlyr_col_1d    &
                                            )

      ! Set mass flux associated with sublimation of snow
      soe_auxvar_id = vsfm_cond_id_for_sublimation;
      call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_SS            , &
                                             VAR_BC_SS_CONDITION  , &
                                             soe_auxvar_id        , &
                                             mflx_sub_snow_col_1d   &
                                            )

      ! Get total mass
      soe_auxvar_id = 1;
      call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL ,       &
                                            VAR_MASS        ,       &
                                            soe_auxvar_id   ,       &
                                            vsfm_mass_col_1d        &
                                           )

      frac_ice(:,:) = 0.d0
      vsfm_fliq_col_1d(:) = 1.d0
      do fc = 1,num_hydrologyc
         c = filter_hydrologyc(fc)
         do j = 1, nlevgrnd

            frac_ice(c,j) = h2osoi_ice(c,j)/(h2osoi_liq(c,j) + h2osoi_ice(c,j))

            idx = (c-bounds%begc)*nlevgrnd + j
            vsfm_fliq_col_1d(idx) = 1._r8 - frac_ice(c,j)
         end do
      end do

      ! Set frac_liq
      soe_auxvar_id = 1;
      call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_INTERNAL  , &
                                             VAR_FRAC_LIQ_SAT , &
                                             soe_auxvar_id    , &
                                             vsfm_fliq_col_1d   &
                                            )

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         do j = 1, nlevgrnd

            idx = (c-bounds%begc)*nlevgrnd + j
            total_mass_flux_et           = total_mass_flux_et           + mflx_et_col_1d(idx)
            total_mass_flux_et_col(c)    = total_mass_flux_et_col(c)    + mflx_et_col_1d(idx)

            total_mass_flux_drain        = total_mass_flux_drain        + mflx_drain_col_1d(idx)
            total_mass_flux_drain_col(c) = total_mass_flux_drain_col(c) + mflx_drain_col_1d(idx)

            mass_beg                     = mass_beg                     + vsfm_mass_col_1d(idx)
            mass_beg_col(c)              = mass_beg_col(c)              + vsfm_mass_col_1d(idx)
            vsfm_mass_prev_col(c,j)      = vsfm_mass_col_1d(idx)
         end do

         idx = c-bounds%begc+1
         total_mass_flux_dew            = total_mass_flux_dew            + mflx_dew_col_1d(idx)
         total_mass_flux_dew_col(c)     = total_mass_flux_dew_col(c)     + mflx_dew_col_1d(idx)

         total_mass_flux_infl           = total_mass_flux_infl           + mflx_infl_col_1d(idx)
         total_mass_flux_infl_col(c)    = total_mass_flux_infl_col(c)    + mflx_infl_col_1d(idx)

         total_mass_flux_snowlyr        = total_mass_flux_snowlyr        + mflx_snowlyr_col_1d(idx)
         total_mass_flux_snowlyr_col(c) = total_mass_flux_snowlyr_col(c) + mflx_snowlyr_col_1d(idx)

         total_mass_flux_sub            = total_mass_flux_sub            + mflx_sub_snow_col_1d(idx)
         total_mass_flux_sub_col(c)     = total_mass_flux_sub_col(c)     + mflx_sub_snow_col_1d(idx)

         total_mass_flux_col(c) = total_mass_flux_et_col(c)      + &
                                  total_mass_flux_infl_col(c)    + &
                                  total_mass_flux_dew_col(c)     + &
                                  total_mass_flux_drain_col(c)   + &
                                  total_mass_flux_snowlyr_col(c) + &
                                  total_mass_flux_sub_col(c)
      end do
      total_mass_flux        = total_mass_flux_et        + &
                               total_mass_flux_infl      + &
                               total_mass_flux_dew       + &
                               total_mass_flux_drain     + &
                               total_mass_flux_snowlyr   + &
                               total_mass_flux_sub

      ! Preform Pre-StepDT operations
      call vsfm_mpp%sysofeqns%PreStepDT()

      ! Get default SNES settings
      call SNESGetTolerances(vsfm_mpp%sysofeqns%snes , &
                             atol_default            , &
                             rtol_default            , &
                             stol_default            , &
                             max_it_default          , &
                             max_f_default           , &
                             ierr                      &
                            )
      CHKERRQ(ierr)

      stol = stol_default

      !
      ! Solve the VSFM.
      !
      iter_count           = 0
      diverged_count       = 0
      mass_bal_err_count   = 0
      abs_mass_error_col   = 0.d0
      successful_step      = PETSC_FALSE

      do

         iter_count = iter_count + 1

         call SNESSetTolerances(vsfm_mpp%sysofeqns%snes , &
                                atol_default            , &
                                rtol_default            , &
                                stol                    , &
                                max_it_default          , &
                                max_f_default           , &
                                ierr                      &
                               );
         CHKERRQ(ierr)

         call vsfm_mpp%sysofeqns%StepDT(dtime, converged, ierr); CHKERRQ(ierr)

         if (.not. converged) then

            ! VSFM solver did not converge, so let's try again with different
            ! solver settings.

            stol             = stol_alternate
            diverged_count   = diverged_count + 1
            successful_step  = PETSC_FALSE

            ! Reduce total run length time by the amount VSFM ran successfully
            ! with previous solver settings
            dtime = dtime - vsfm_mpp%sysofeqns%time

            if (diverged_count > 1) then
               ! Set frac_liq
               vsfm_fliq_col_1d(:) = 1.d0
               soe_auxvar_id = 1;

               call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_INTERNAL  , &
                                                      VAR_FRAC_LIQ_SAT , &
                                                      soe_auxvar_id    , &
                                                      vsfm_fliq_col_1d   &
                                                     )
            end if
         else

            ! Solver converged, so let's copy data from VSFM model to
            ! CLM's data structure.

            ! Get Liquid saturation
            soe_auxvar_id = 1;
            call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL , &
                                                  VAR_LIQ_SAT     , &
                                                  soe_auxvar_id   , &
                                                  vsfm_sat_col_1d   &
                                                 )

            ! Get total mass
            soe_auxvar_id = 1;
            call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL  , &
                                                  VAR_MASS         , &
                                                  soe_auxvar_id    , &
                                                  vsfm_mass_col_1d   &
                                                 )

            ! Get liquid soil matrix potential
            soe_auxvar_id = 1;
            call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL     , &
                                                  VAR_SOIL_MATRIX_POT , &
                                                  soe_auxvar_id       , &
                                                  vsfm_smpl_col_1d      &
                                                 )

            ! Get soil liquid pressure. This is the prognostic state of VSFM
            ! and needs to be saved in the restart file.
            soe_auxvar_id = 1;
            call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL   , &
                                                  VAR_PRESSURE      , &
                                                  soe_auxvar_id     , &
                                                  vsfm_soilp_col_1d   &
                                                 )

            ! Put the data in CLM's data structure
            mass_end        = 0.d0
            do fc = 1,num_hydrologyc
               c = filter_hydrologyc(fc)

               ! initialization
               jwt = -1

               ! Loops in decreasing j so WTD can be computed in the same loop
               do j = nlevgrnd, 1, -1
                  idx = (c-bounds%begc)*nlevgrnd + j

                  h2osoi_liq(c,j) = (1.d0 - frac_ice(c,j))*vsfm_mass_col_1d(idx)
                  h2osoi_ice(c,j) = frac_ice(c,j)         *vsfm_mass_col_1d(idx)

                  mass_end        = mass_end        + vsfm_mass_col_1d(idx)
                  mass_end_col(c) = mass_end_col(c) + vsfm_mass_col_1d(idx)

                  vsfm_dmass_col(c) = vsfm_dmass_col(c) + &
                                      (vsfm_mass_col_1d(idx)-vsfm_mass_prev_col(c,j))

                  smp_l(c,j)        = vsfm_smpl_col_1d(idx)*1000.0_r8      ! [m] --> [mm]

                  if (jwt == -1) then
                     ! Find the first soil that is unsaturated
                     if (smp_l(c,j) < 0._r8) jwt = j
                  end if

               end do

               ! Find maximum water balance error over the column
               abs_mass_error_col = max(abs_mass_error_col,                     &
                                        abs(mass_beg_col(c) - mass_end_col(c) + &
                                            total_mass_flux_col(c)*get_step_size()))

               qcharge(c) = 0._r8

               if (jwt == -1 .or. jwt == nlevgrnd) then
                  ! Water table below or in the last layer
                  zwt(c) = zi(c,nlevgrnd)
               else
                  z_dn = (zi(c,jwt-1) + zi(c,jwt  ))/2._r8
                  z_up = (zi(c,jwt ) + zi(c,jwt+1))/2._r8
                  zwt(c) = (0._r8 - smp_l(c,jwt))/(smp_l(c,jwt) - smp_l(c,jwt+1))*(z_dn - z_up) + z_dn
               endif
            end do

            ! Save soil liquid pressure from VSFM for all (active+nonactive) cells.
            ! soilp_col is used for restarting VSFM.
            do c = bounds%begc, bounds%endc
               do j = 1, nlevgrnd
                  idx = (c-bounds%begc)*nlevgrnd + j
                  soilp_col(c,j) = vsfm_soilp_col_1d(idx)
               end do
            end do

            ! For the solution that did converge, is the mass error acceptable?
            if (abs_mass_error_col >= max_abs_mass_error_col) then

               ! For the solution that converged, the mass error
               ! is unacceptable. So let's try again with tighter
               ! solution tolerance (stol) for SNES.

               mass_bal_err_count  = mass_bal_err_count + 1
               stol                = stol/10._r8
               dtime               = get_step_size()
               successful_step     = PETSC_FALSE
               abs_mass_error_col  = 0._r8
               mass_end_col(:)     = 0._r8

               ! Perform Pre-StepDT operations
               call vsfm_mpp%sysofeqns%PreStepDT()

            else

               successful_step  = PETSC_TRUE

            endif

         endif

         if (successful_step) exit

         if (iter_count >= max_iter_count) then
            write(iulog,*)'In soilwater_vsfm: VSFM failed to converge after multiple attempts.'
            call endrun(msg=errMsg(__FILE__, __LINE__))
         end if

      end do

      call SNESSetTolerances(vsfm_mpp%sysofeqns%snes, atol_default, rtol_default, stol_default, &
                             max_it_default, max_f_default, ierr); CHKERRQ(ierr)

      call vsfm_mpp%sysofeqns%PostStepDT()

#if VSFM_DEBUG
      write(iulog,*)'VSFM-DEBUG: nstep                      = ',get_nstep()
      write(iulog,*)'VSFM-DEBUG: dtime                      = ',get_step_size()
      write(iulog,*)'VSFM-DEBUG: change in mass between dt  = ',-(mass_beg - mass_end)
      write(iulog,*)'VSFM-DEBUG: change in mass due to flux = ',total_mass_flux*get_step_size()
      write(iulog,*)'VSFM-DEBUG: Error in mass conservation = ',mass_beg - mass_end + total_mass_flux*get_step_size()
      write(iulog,*)'VSFM-DEBUG: et_flux    * dtime         = ',total_mass_flux_et*get_step_size()
      write(iulog,*)'VSFM-DEBUG: infil_flux * dtime         = ',total_mass_flux_infl*get_step_size()
      write(iulog,*)'VSFM-DEBUG: dew_flux   * dtime         = ',total_mass_flux_dew*get_step_size()
      write(iulog,*)'VSFM-DEBUG: drain_flux * dtime         = ',total_mass_flux_drain*get_step_size()
      write(iulog,*)'VSFM-DEBUG: snow_flux  * dtime         = ',total_mass_flux_snowlyr*get_step_size()
      write(iulog,*)'VSFM-DEBUG: sub_flux   * dtime         = ',total_mass_flux_sub*get_step_size()
      write(iulog,*)'VSFM-DEBUG: total_mass_flux            = ',total_mass_flux/flux_unit_conversion
      write(iulog,*)'VSFM-DEBUG: et_flux                    = ',total_mass_flux_et
      write(iulog,*)'VSFM-DEBUG: infil_flux                 = ',total_mass_flux_infl
      write(iulog,*)'VSFM-DEBUG: dew_flux                   = ',total_mass_flux_dew
      write(iulog,*)'VSFM-DEBUG: drain_flux                 = ',total_mass_flux_drain
      write(iulog,*)'VSFM-DEBUG: snow_flux                  = ',total_mass_flux_snowlyr
      write(iulog,*)'VSFM-DEBUG: sub_flux                   = ',total_mass_flux_sub
      write(iulog,*)''
#endif

#endif

      ! compute the water deficit and reset negative liquid water content
      !  Jinyun Tang
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         qflx_deficit(c) = 0._r8
      enddo

    end associate

   end subroutine soilwater_vsfm

 end module SoilWaterMovementMod
