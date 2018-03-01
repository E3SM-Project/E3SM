module MPPVSFMALM_Driver

#ifdef USE_PETSC_LIB

  implicit none
  
  public :: MPPVSFMALM_Solve

  !------------------------------------------------------------------------
contains

  !-----------------------------------------------------------------------
  subroutine MPPVSFMALM_Solve(bounds, num_hydrologyc, filter_hydrologyc, &
       soilhydrology_vars, soilstate_vars, &
       waterflux_vars, waterstate_vars, temperature_vars)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
#include <petsc/finclude/petsc.h>
    use petscsnes

    use decompMod                 , only : bounds_type
    use shr_kind_mod              , only : r8 => shr_kind_r8
    use mpp_varcon                , only : denh2o
    use mpp_varpar                , only : nlevsoi, max_patch_per_col, nlevgrnd
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
    use mpp_varcon                , only : istsoil, istcrop
    use mpp_varctl                , only : iulog
    use shr_log_mod               , only : errMsg => shr_log_errMsg
    use mpp_varctl                , only : lateral_connectivity
    use mpp_varctl                , only : vsfm_lateral_model_type
    use mpp_varctl                , only : vsfm_include_seepage_bc
    use MultiPhysicsProbVSFM      , only : vsfm_mpp
    use MultiPhysicsProbConstants , only : VAR_BC_SS_CONDITION
    use MultiPhysicsProbConstants , only : VAR_TEMPERATURE
    use MultiPhysicsProbConstants , only : VAR_PRESSURE
    use MultiPhysicsProbConstants , only : VAR_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_FRAC_LIQ_SAT
    use MultiPhysicsProbConstants , only : VAR_MASS
    use MultiPhysicsProbConstants , only : VAR_SOIL_MATRIX_POT
    use MultiPhysicsProbConstants , only : VAR_LATERAL_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : VAR_BC_MASS_EXCHANGED
    use MultiPhysicsProbConstants , only : AUXVAR_INTERNAL
    use MultiPhysicsProbConstants , only : AUXVAR_BC
    use MultiPhysicsProbConstants , only : AUXVAR_SS
    use domainLateralMod          , only : ExchangeColumnLevelGhostData, ldomain_lateral
    use MPPVSFMALM_Initialize     , only : vsfm_cond_id_for_infil
    use MPPVSFMALM_Initialize     , only : vsfm_cond_id_for_et
    use MPPVSFMALM_Initialize     , only : vsfm_cond_id_for_dew
    use MPPVSFMALM_Initialize     , only : vsfm_cond_id_for_drainage
    use MPPVSFMALM_Initialize     , only : vsfm_cond_id_for_snow
    use MPPVSFMALM_Initialize     , only : vsfm_cond_id_for_sublimation
    use MPPVSFMALM_Initialize     , only : vsfm_cond_id_for_lateral_flux
    use mpp_bounds                , only : bounds_proc_begc_all, bounds_proc_endc_all
    use mpp_bounds                , only : bounds_proc_begc, bounds_proc_endc
    !
    ! !ARGUMENTS:
    implicit none
    !
    type(bounds_type)       , intent(in)    :: bounds               ! bounds
    integer                 , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
    integer                 , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
    type(soilhydrology_type), intent(inout) :: soilhydrology_vars
    type(soilstate_type)    , intent(inout) :: soilstate_vars
    type(waterflux_type)    , intent(inout) :: waterflux_vars
    type(waterstate_type)   , intent(inout) :: waterstate_vars
    type(temperature_type)  , intent(in)    :: temperature_vars
    !
    ! !LOCAL VARIABLES:
    integer              :: p,c,fc,j,g                                                       ! do loop indices
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
    real(r8)             :: total_mass_flux_lateral_col (bounds%begc:bounds%endc)            ! Lateral flux computed by VSFM solver at column level
    real(r8)             :: total_mass_flux_seepage_col (bounds%begc:bounds%endc)            ! Seepage flux computed by VSFM solver at column level
    real(r8)             :: qflx_seepage                (bounds%begc:bounds%endc)            ! Seepage flux computed by VSFM solver at column level
    real(r8)             :: vsfm_mass_prev_col          (bounds%begc:bounds%endc,1:nlevgrnd) ! Mass of water before a VSFM solve
    real(r8)             :: vsfm_dmass_col              (bounds%begc:bounds%endc)            ! Change in mass of water after a VSFM solve
    real(r8)             :: mass_beg_col                (bounds%begc:bounds%endc)            ! Total mass before a VSFM solve
    real(r8)             :: mass_end_col                (bounds%begc:bounds%endc)            ! Total mass after a VSFM solve
    integer              :: ier                                                              ! error status

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
    PetscInt             :: converged_reason                                                 ! SNES converged due to which criteria
    PetscReal            :: atol_default                                                     ! Default SNES absolute convergance tolerance
    PetscReal            :: rtol_default                                                     ! Default SNES relative convergance tolerance
    PetscReal            :: stol_default                                                     ! Default SNES solution convergance tolerance
    PetscInt             :: max_it_default                                                   ! Default SNES maximum number of iteration
    PetscInt             :: max_f_default                                                    ! Default SNES maximum number of function evaluation
    PetscReal            :: stol                                                             ! solution convergance tolerance
    PetscReal            :: rtol                                                             ! relative convergance tolerance
    PetscReal,parameter  :: stol_alternate = 1.d-10                                          ! Alternate solution convergance tolerance

    PetscReal            :: mass_beg                                                         ! Sum of mass of water for all active soil columns before VSFM is called
    PetscReal            :: mass_end                                                         ! Sum of mass of water for all active soil columns after VSFM is called
    PetscReal            :: total_mass_flux_et                                               ! Sum of mass ET mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_infl                                             ! Sum of mass infiltration mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_dew                                              ! Sum of mass dew mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_drain                                            ! Sum of mass drainage mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_snowlyr                                          ! Sum of mass snow layer disappearance mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_sub                                              ! Sum of mass sublimation mass flux of water for all active soil columns
    PetscReal            :: total_mass_flux_lateral                                          ! Sum of lateral mass flux for all active soil columns
    PetscReal            :: total_mass_flux                                                  ! Sum of mass ALL mass flux of water for all active soil columns
    PetscInt             :: iter_count                                                       ! How many times VSFM solver is called

    PetscInt, parameter  :: max_iter_count = 10                                              ! Maximum number of times VSFM can be called
    PetscInt             :: diverged_count                                                   ! Number of time VSFM solver diverged
    PetscInt             :: mass_bal_err_count                                               ! Number of time VSFM solver returns a solution that isn't within acceptable mass balance error threshold
    PetscReal            :: abs_mass_error_col                                               ! Maximum absolute error for any active soil column
    PetscReal, parameter :: max_abs_mass_error_col  = 1.e-5                                  ! Acceptable mass balance error
    PetscBool            :: successful_step                                                  ! Is the solution return by VSFM acceptable
    PetscReal, pointer   :: vsfm_soilp_col_ghosted_1d(:)
    PetscReal, pointer   :: vsfm_fliq_col_ghosted_1d(:)
    PetscReal, pointer   :: mflx_lateral_col_1d(:)
    PetscReal, pointer   :: lat_mass_exc_col_1d(:)
    PetscReal, pointer   :: seepage_mass_exc_col_1d(:)
    PetscReal, pointer   :: seepage_press_1d(:)
    !-----------------------------------------------------------------------

    associate(& 
         z                         =>    col%z                                      , & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         zi                        =>    col%zi                                     , & ! Input:  [real(r8) (:,:) ]  interface level below a "z" level (m)
         dz                        =>    col%dz                                     , & ! Input:  [real(r8) (:,:) ]  layer thickness (m)
         snl                       =>    col%snl                                    , & ! Input:  [integer  (:)   ]  minus number of snow layers

         qcharge                   =>    soilhydrology_vars%qcharge_col             , & ! Input:  [real(r8) (:)   ]  aquifer recharge rate (mm/s)
         zwt                       =>    soilhydrology_vars%zwt_col                 , & ! Input:  [real(r8) (:)   ]  water table depth (m)

         rootr_col                 =>    soilstate_vars%rootr_col                   , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer
         rootr_pft                 =>    soilstate_vars%rootr_patch                 , & ! Input:  [real(r8) (:,:) ]  effective fraction of roots in each soil layer
         smp_l                     =>    soilstate_vars%smp_l_col                   , & ! Output: [real(r8) (:,:) ]  soil matrix potential [mm]

         h2osoi_ice                =>    waterstate_vars%h2osoi_ice_col             , & ! Input:  [real(r8) (:,:) ]  ice water (kg/m2)
         h2osoi_liq                =>    waterstate_vars%h2osoi_liq_col             , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2)
         h2osoi_vol                =>    waterstate_vars%h2osoi_vol_col             , & ! Input:  [real(r8) (:,:) ]  volumetric soil water (0<=h2osoi_vol<=watsat) [m3/m3]
         frac_h2osfc               =>    waterstate_vars%frac_h2osfc_col            , & ! Input:  [real(r8) (:)   ]

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
         qflx_lateral              =>    waterflux_vars%qflx_lateral_col            , & ! Input:  [real(r8) (:)   ]  lateral flux of water to neighboring column (mm H2O /s)
         qflx_surf                 =>    waterflux_vars%qflx_surf_col               , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)
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
      total_mass_flux_lateral          = 0.d0

      mass_beg_col(:)                  = 0.d0
      mass_end_col(:)                  = 0.d0
      total_mass_flux_col(:)           = 0.d0
      total_mass_flux_et_col(:)        = 0.d0
      total_mass_flux_infl_col(:)      = 0.d0
      total_mass_flux_dew_col(:)       = 0.d0
      total_mass_flux_drain_col(:)     = 0.d0
      total_mass_flux_snowlyr_col(:)   = 0.d0
      total_mass_flux_sub_col(:)       = 0.d0
      total_mass_flux_lateral_col(:)   = 0.d0

      vsfm_mass_prev_col(:,:)          = 0.d0
      vsfm_dmass_col(:)                = 0.d0

      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         if (lateral_connectivity) then
            g    = col%gridCell(c)
            area = ldomain_lateral%ugrid%areaGrid_ghosted(g)
         endif

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
         mflx_snowlyr_col_1d(c-bounds%begc+1) = mflx_snowlyr_col(c)*area + &
                                                mflx_neg_snow_col_1d(c-bounds%begc+1)*area
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

      if (vsfm_lateral_model_type == 'source_sink') then


         allocate(vsfm_soilp_col_ghosted_1d((bounds_proc_endc_all - bounds_proc_begc_all+1)*nlevgrnd))
         allocate(vsfm_fliq_col_ghosted_1d( (bounds_proc_endc_all - bounds_proc_begc_all+1)*nlevgrnd))
         allocate(mflx_lateral_col_1d( (bounds_proc_endc - bounds_proc_begc+1)*nlevgrnd))

         soe_auxvar_id = 1;
         call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL   , &
                                               VAR_PRESSURE      , &
                                               soe_auxvar_id     , &
                                               vsfm_soilp_col_1d   &
                                               )

         call ExchangeColumnLevelGhostData(bounds, nlevgrnd, vsfm_soilp_col_1d, vsfm_soilp_col_ghosted_1d)
         call ExchangeColumnLevelGhostData(bounds, nlevgrnd, vsfm_fliq_col_1d,  vsfm_fliq_col_ghosted_1d )

         soe_auxvar_id = 1;
         call vsfm_mpp%sysofeqns%SetDataFromCLMForGhost(AUXVAR_INTERNAL           , &
                                                        VAR_PRESSURE              , &
                                                        soe_auxvar_id             , &
                                                        vsfm_soilp_col_ghosted_1d   &
                                                       )

         soe_auxvar_id = 1;
         call vsfm_mpp%sysofeqns%SetDataFromCLMForGhost(AUXVAR_INTERNAL          , &
                                                        VAR_FRAC_LIQ_SAT         , &
                                                        soe_auxvar_id            , &
                                                        vsfm_fliq_col_ghosted_1d   &
                                                       )

         call vsfm_mpp%sysofeqns%ComputeLateralFlux(dtime)

         soe_auxvar_id = vsfm_cond_id_for_lateral_flux;
         call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_SS   , &
                                               VAR_BC_SS_CONDITION      , &
                                               soe_auxvar_id     , &
                                               mflx_lateral_col_1d   &
                                               )

         do fc = 1, num_hydrologyc
            c = filter_hydrologyc(fc)

            g    = col%gridCell(c)
            area = ldomain_lateral%ugrid%areaGrid_ghosted(g)

            ! [mm/s] --> [kg/s]   [m^2] [kg/m^3]  [m/mm]
            flux_unit_conversion     = area * denh2o * 1.0d-3

            qflx_lateral(c) = 0._r8
            do j = 1, nlevgrnd
               idx = (c-bounds%begc)*nlevgrnd + j

               total_mass_flux_lateral_col(c) =      &
                    total_mass_flux_lateral_col(c) + &
                    mflx_lateral_col_1d(idx)

               qflx_lateral(c) = qflx_lateral(c) - &
                    mflx_lateral_col_1d(idx)/flux_unit_conversion
            enddo

            total_mass_flux_lateral = total_mass_flux_lateral + &
                 total_mass_flux_lateral_col(c)
         enddo

         deallocate(vsfm_soilp_col_ghosted_1d )
         deallocate(vsfm_fliq_col_ghosted_1d  )
         deallocate(mflx_lateral_col_1d       )

         if (vsfm_include_seepage_bc) then
            allocate(seepage_press_1d( (bounds_proc_endc - bounds_proc_begc+1)))
            seepage_press_1d(:) = 101325.d0
            soe_auxvar_id = 1
            call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_BC,  &
                 VAR_BC_SS_CONDITION, soe_auxvar_id, seepage_press_1d)
            deallocate(seepage_press_1d)
         endif

      else if (vsfm_lateral_model_type == 'three_dimensional') then

         if (vsfm_include_seepage_bc) then
            allocate(seepage_press_1d( (bounds_proc_endc - bounds_proc_begc+1)))
            seepage_press_1d(:) = 101325.d0
            soe_auxvar_id = 1
            call vsfm_mpp%sysofeqns%SetDataFromCLM(AUXVAR_BC,  &
                 VAR_BC_SS_CONDITION, soe_auxvar_id, seepage_press_1d)
            deallocate(seepage_press_1d)
         endif

      endif

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
                                  total_mass_flux_sub_col(c)     + &
                                  total_mass_flux_lateral_col(c)
      end do
      total_mass_flux        = total_mass_flux_et        + &
                               total_mass_flux_infl      + &
                               total_mass_flux_dew       + &
                               total_mass_flux_drain     + &
                               total_mass_flux_snowlyr   + &
                               total_mass_flux_sub       + &
                               total_mass_flux_lateral

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
      rtol = rtol_default

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
                                rtol                    , &
                                stol                    , &
                                max_it_default          , &
                                max_f_default           , &
                                ierr                      &
                               );
         CHKERRQ(ierr)

         call vsfm_mpp%sysofeqns%StepDT(dtime, get_nstep(), &
              converged, converged_reason, ierr); CHKERRQ(ierr)

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

            qflx_seepage(:) = 0._r8

            if (vsfm_lateral_model_type == 'source_sink') then

               ! Get following fluxes from VSFM:
               ! (i) seepage mass exchanged.

               allocate(seepage_mass_exc_col_1d( (bounds_proc_endc - bounds_proc_begc+1)         ))
               seepage_mass_exc_col_1d = 0.d0

               if (vsfm_include_seepage_bc) then
                  soe_auxvar_id = 1;
                  call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_BC              ,  &
                                                        VAR_BC_MASS_EXCHANGED  ,  &
                                                        soe_auxvar_id          ,  &
                                                        seepage_mass_exc_col_1d   &
                                                        )
               endif

               do fc = 1, num_hydrologyc
                  c = filter_hydrologyc(fc)

                  g    = col%gridCell(c)
                  area = ldomain_lateral%ugrid%areaGrid_ghosted(g)

                  ! [mm/s] --> [kg/s]   [m^2] [kg/m^3]  [m/mm]
                  flux_unit_conversion     = area * denh2o * 1.0d-3

                  idx = (c-bounds%begc) + 1

                  qflx_seepage(c)                = seepage_mass_exc_col_1d(idx)/flux_unit_conversion/dtime
                  total_mass_flux_seepage_col(c) = -seepage_mass_exc_col_1d(idx)/dtime

                  total_mass_flux_col(c) = total_mass_flux_et_col(c)      + &
                                           total_mass_flux_infl_col(c)    + &
                                           total_mass_flux_dew_col(c)     + &
                                           total_mass_flux_drain_col(c)   + &
                                           total_mass_flux_snowlyr_col(c) + &
                                           total_mass_flux_sub_col(c)     + &
                                           total_mass_flux_lateral_col(c) + &
                                           total_mass_flux_seepage_col(c)
               enddo

            else if (vsfm_lateral_model_type == 'three_dimensional') then

               ! Get following fluxes from VSFM:
               ! (i) lateral mass exchanged, and
               ! (ii) seepage mass exchanged.

               allocate(lat_mass_exc_col_1d(     (bounds_proc_endc - bounds_proc_begc+1)*nlevgrnd))
               allocate(seepage_mass_exc_col_1d( (bounds_proc_endc - bounds_proc_begc+1)         ))

               lat_mass_exc_col_1d(:) = 0.d0
               seepage_mass_exc_col_1d(:) = 0.d0

               soe_auxvar_id = 1
               call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_INTERNAL            , &
                                                     VAR_LATERAL_MASS_EXCHANGED , &
                                                     soe_auxvar_id              , &
                                                     lat_mass_exc_col_1d          &
                                                     )

               if (vsfm_include_seepage_bc) then
                  soe_auxvar_id = 1;
                  call vsfm_mpp%sysofeqns%GetDataForCLM(AUXVAR_BC              ,  &
                                                        VAR_BC_MASS_EXCHANGED  ,  &
                                                        soe_auxvar_id          ,  &
                                                        seepage_mass_exc_col_1d   &
                                                        )
               endif

               total_mass_flux_lateral_col(:)   = 0.d0
               total_mass_flux_lateral          = 0.d0

               do fc = 1, num_hydrologyc
                  c = filter_hydrologyc(fc)

                  g    = col%gridCell(c)
                  area = ldomain_lateral%ugrid%areaGrid_ghosted(g)

                  ! [mm/s] --> [kg/s]   [m^2] [kg/m^3]  [m/mm]
                  flux_unit_conversion     = area * denh2o * 1.0d-3

                  qflx_lateral(c) = 0._r8
                  do j = 1, nlevgrnd
                     idx = (c-bounds%begc)*nlevgrnd + j

                     total_mass_flux_lateral_col(c) = &
                          total_mass_flux_lateral_col(c) + &
                          lat_mass_exc_col_1d(idx)/dtime

                     qflx_lateral(c)  = qflx_lateral(c) - &
                          lat_mass_exc_col_1d(idx)/flux_unit_conversion/dtime

                  enddo

                  idx = (c-bounds%begc) + 1
                  qflx_seepage(c)                = seepage_mass_exc_col_1d(idx)/flux_unit_conversion/dtime
                  total_mass_flux_seepage_col(c) = -seepage_mass_exc_col_1d(idx)/dtime

                  total_mass_flux_col(c) = total_mass_flux_et_col(c)      + &
                                           total_mass_flux_infl_col(c)    + &
                                           total_mass_flux_dew_col(c)     + &
                                           total_mass_flux_drain_col(c)   + &
                                           total_mass_flux_snowlyr_col(c) + &
                                           total_mass_flux_sub_col(c)     + &
                                           total_mass_flux_lateral_col(c) + &
                                           total_mass_flux_seepage_col(c)

                  total_mass_flux_lateral = total_mass_flux_lateral + &
                       total_mass_flux_lateral_col(c)
               enddo

               deallocate(lat_mass_exc_col_1d)
            endif

            ! Put the data in CLM's data structure
            mass_end        = 0.d0
            area            = 1.d0 ! [m^2]

            do fc = 1,num_hydrologyc
               c = filter_hydrologyc(fc)

               if (lateral_connectivity) then
                  g    = col%gridCell(c)
                  area = ldomain_lateral%ugrid%areaGrid_ghosted(g)
               endif

               ! initialization
               jwt = -1

               ! Loops in decreasing j so WTD can be computed in the same loop
               do j = nlevgrnd, 1, -1
                  idx = (c-bounds%begc)*nlevgrnd + j

                  h2osoi_liq(c,j) = (1.d0 - frac_ice(c,j))*vsfm_mass_col_1d(idx)/area
                  h2osoi_ice(c,j) = frac_ice(c,j)         *vsfm_mass_col_1d(idx)/area

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

               if (converged_reason == SNES_CONVERGED_FNORM_RELATIVE) then
                  rtol = rtol/10._r8
               else if (converged_reason == SNES_CONVERGED_SNORM_RELATIVE) then
                  stol = stol/10._r8
               endif

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

      ! Add seepage flux from VSFM to surface runoff
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)

         qflx_surf(c) = qflx_surf(c) + qflx_seepage(c)
      enddo

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
      write(iulog,*)'VSFM-DEBUG: lat_flux   * dtime         = ',total_mass_flux_lateral*get_step_size()
      write(iulog,*)'VSFM-DEBUG: total_mass_flux            = ',total_mass_flux/flux_unit_conversion
      write(iulog,*)'VSFM-DEBUG: et_flux                    = ',total_mass_flux_et
      write(iulog,*)'VSFM-DEBUG: infil_flux                 = ',total_mass_flux_infl
      write(iulog,*)'VSFM-DEBUG: dew_flux                   = ',total_mass_flux_dew
      write(iulog,*)'VSFM-DEBUG: drain_flux                 = ',total_mass_flux_drain
      write(iulog,*)'VSFM-DEBUG: snow_flux                  = ',total_mass_flux_snowlyr
      write(iulog,*)'VSFM-DEBUG: sub_flux                   = ',total_mass_flux_sub
      write(iulog,*)''
#endif

      ! compute the water deficit and reset negative liquid water content
      !  Jinyun Tang
      do fc = 1, num_hydrologyc
         c = filter_hydrologyc(fc)
         qflx_deficit(c) = 0._r8
      enddo

    end associate

  end subroutine MPPVSFMALM_Solve

#endif

end module MPPVSFMALM_Driver
