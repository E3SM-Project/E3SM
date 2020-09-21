module dynSubgridDriverMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! High-level routines for dynamic subgrid areas (prescribed transient Patches and
  ! dynamic landunits).
  !
  ! !USES:
  use dynSubgridControlMod, only : get_flanduse_timeseries
  use dynSubgridControlMod, only : get_do_transient_pfts, get_do_transient_crops
  use dynSubgridControlMod, only : get_do_harvest
  use dynPriorWeightsMod  , only : prior_weights_type
  use dynPatchStateUpdaterMod      , only : patch_state_updater_type
  use dynColumnStateUpdaterMod     , only : column_state_updater_type
  use UrbanParamsType     , only : urbanparams_type
  use CanopyStateType     , only : canopystate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNStateType         , only : cnstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use EnergyFluxType      , only : energyflux_type
  use LakeStateType       , only : lakestate_type
  use PhotosynthesisType  , only : photosyns_type
  use SoilHydrologyType   , only : soilhydrology_type  
  use SoilStateType       , only : soilstate_type
  use WaterfluxType       , only : waterflux_type
  use WaterstateType      , only : waterstate_type
  use TemperatureType     , only : temperature_type
  use glc2lndMod          , only : glc2lnd_type
  use dynLandunitAreaMod  , only : update_landunit_weights
  use CropType            , only : crop_type
  use PhosphorusStateType , only : phosphorusstate_type
  use PhosphorusFluxType  , only : phosphorusflux_type
  use dyncropFileMod      , only : dyncrop_init, dyncrop_interp
  use filterMod           , only : filter, filter_inactive_and_active

  use GridcellDataType    , only : gridcell_carbon_state, gridcell_carbon_flux
  use ColumnDataType      , only : column_carbon_state, column_carbon_flux
  use ColumnDataType      , only : col_ns, col_ps
  use VegetationDataType  , only : vegetation_carbon_state, veg_ns, veg_ps

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save

  public :: dynSubgrid_init             ! initialize transient land cover
  public :: dynSubgrid_driver           ! top-level driver for transient land cover
  public :: dynSubgrid_wrapup_weight_changes ! reconcile various variables after subgrid weights change
  !
  ! !PRIVATE TYPES:
  ! saved weights from before the subgrid weight updates
  type(prior_weights_type) :: prior_weights

  ! object used to update patch-level states after subgrid weight updates
  type(patch_state_updater_type), target :: patch_state_updater

  ! object used to update column-level states after subgrid weight updates
  type(column_state_updater_type), target :: column_state_updater
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_init(bounds, glc2lnd_vars, crop_vars)
    !
    ! !DESCRIPTION:
    ! Determine initial subgrid weights for prescribed transient Patches and/or
    ! dynamic landunits. Note that these weights will be overwritten in a restart run.
    !
    ! This should be called from initialization. 
    !
    ! Note that dynpft_init / dynpft_interp need to be called from outside any loops over
    ! clumps - so this routine needs to be called from outside any loops over clumps.
    !
    ! !USES:
    use decompMod         , only : bounds_type, BOUNDS_LEVEL_PROC
    use decompMod         , only : get_proc_clumps, get_clump_bounds
    use dynpftFileMod     , only : dynpft_init
    use dynHarvestMod     , only : dynHarvest_init
    use dynpftFileMod     , only : dynpft_interp
    !
    ! !ARGUMENTS:
    type(bounds_type) , intent(in)    :: bounds  ! processor-level bounds
    type(glc2lnd_type), intent(inout) :: glc2lnd_vars
    type(crop_type)   , intent(inout) :: crop_vars
    !
    ! !LOCAL VARIABLES:
    integer           :: nclumps      ! number of clumps on this processor
    integer           :: nc           ! clump index
    type(bounds_type) :: bounds_clump ! clump-level bounds
    character(len=*), parameter :: subname = 'dynSubgrid_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    nclumps = get_proc_clumps()

    prior_weights        = prior_weights_type(bounds)
    patch_state_updater  = patch_state_updater_type(bounds)
    column_state_updater = column_state_updater_type(bounds, nclumps)

    ! Initialize stuff for prescribed transient Patches
    if (get_do_transient_pfts()) then
       call dynpft_init(bounds, dynpft_filename=get_flanduse_timeseries())
    end if

    ! Initialize stuff for harvest (currently shares the flanduse_timeseries file)
    if (get_do_harvest()) then
       call dynHarvest_init(bounds, harvest_filename=get_flanduse_timeseries())
    end if

    ! Initialize stuff for prescribed transient crops
    if (get_do_transient_crops()) then
       call dyncrop_init(bounds, dyncrop_filename=get_flanduse_timeseries())
    end if

    ! ------------------------------------------------------------------------
    ! Set initial subgrid weights for aspects that are read from file. This is relevant
    ! for cold start and use_init_interp-based initialization.
    ! ------------------------------------------------------------------------

    if (get_do_transient_pfts()) then
       call dynpft_interp(bounds)
    end if

    if (get_do_transient_crops()) then
       call dyncrop_interp(bounds, crop_vars)
    end if

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)
       call dynSubgrid_wrapup_weight_changes(bounds_clump, glc2lnd_vars)
    end do
    !$OMP END PARALLEL DO
    
  end subroutine dynSubgrid_init

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_driver(bounds_proc, &
       urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
       waterstate_vars, waterflux_vars, temperature_vars, energyflux_vars, &
       canopystate_vars, photosyns_vars, cnstate_vars, &
       veg_cs, c13_veg_cs, c14_veg_cs, &
       col_cs, c13_col_cs, c14_col_cs, col_cf, &
       grc_cs, grc_cf, &
       carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
       nitrogenstate_vars, nitrogenflux_vars, glc2lnd_vars,&
       phosphorusstate_vars,phosphorusflux_vars, crop_vars)
    !
    ! !DESCRIPTION:
    ! Update subgrid weights for prescribed transient Patches and/or dynamic
    ! landunits. Also do related adjustments (water & energy, carbon & nitrogen).
    !
    ! This should be called every time step in CLM's run loop.
    !
    ! Note that this routine operates partly at the proc-level (outside an OMP region),
    ! and partly at the clump level (inside OMP regions). Thus, this must be called from
    ! OUTSIDE any loops over clumps in the driver.
    !
    ! !USES:
    use elm_varctl           , only : use_cn, create_glacier_mec_landunit, use_fates
    use decompMod            , only : bounds_type, get_proc_clumps, get_clump_bounds
    use decompMod            , only : BOUNDS_LEVEL_PROC
    use dynInitColumnsMod    , only : initialize_new_columns
    use dynConsBiogeophysMod , only : dyn_hwcontent_init, dyn_hwcontent_final
    use dynConsBiogeochemMod , only : dyn_cnbal_patch, dyn_cnbal_column
    use dynpftFileMod        , only : dynpft_interp
    use dynHarvestMod        , only : dynHarvest_interp
    use dynEDMod             , only : dyn_ED
    use reweightMod          , only : reweight_wrapup
    use subgridWeightsMod    , only : compute_higher_order_weights, set_subgrid_diagnostic_fields
    use CarbonStateUpdate1Mod   , only : CarbonStateUpdateDynPatch
    use NitrogenStateUpdate1Mod   , only : NitrogenStateUpdateDynPatch
    use PhosphorusStateUpdate1Mod     , only : PhosphorusStateUpdateDynPatch
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds_proc  ! processor-level bounds
    type(urbanparams_type)   , intent(in)    :: urbanparams_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(soilhydrology_type) , intent(inout) :: soilhydrology_vars
    type(lakestate_type)     , intent(in)    :: lakestate_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(energyflux_type)    , intent(inout) :: energyflux_vars
    type(canopystate_type)   , intent(inout) :: canopystate_vars
    type(photosyns_type)     , intent(inout) :: photosyns_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(vegetation_carbon_state), intent(inout) :: veg_cs
    type(vegetation_carbon_state), intent(inout) :: c13_veg_cs
    type(vegetation_carbon_state), intent(inout) :: c14_veg_cs
    type(column_carbon_state)    , intent(inout) :: col_cs
    type(column_carbon_state)    , intent(inout) :: c13_col_cs
    type(column_carbon_state)    , intent(inout) :: c14_col_cs
    type(column_carbon_flux)     , intent(inout) :: col_cf
    type(gridcell_carbon_state)  , intent(inout) :: grc_cs
    type(gridcell_carbon_flux)   , intent(inout) :: grc_cf
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(glc2lnd_type)       , intent(inout) :: glc2lnd_vars

    type(phosphorusstate_type) , intent(inout)    :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    type(crop_type)          , intent(inout) :: crop_vars

    !
    ! !LOCAL VARIABLES:
    integer           :: nclumps      ! number of clumps on this processor
    integer           :: nc           ! clump index
    type(bounds_type) :: bounds_clump ! clump-level bounds

    character(len=*), parameter :: subname = 'dynSubgrid_driver'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds_proc%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    nclumps = get_proc_clumps()
    
    ! ==========================================================================
    ! Do initialization, prior to land cover change
    ! ==========================================================================

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)

       call dyn_hwcontent_init(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_lakec, filter(nc)%lakec, &
            urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
            waterstate_vars, waterflux_vars, temperature_vars, energyflux_vars)

       call prior_weights%set_prior_weights(bounds_clump)
       call patch_state_updater%set_old_weights(bounds_clump)
       call column_state_updater%set_old_weights(bounds_clump)
    end do
    !$OMP END PARALLEL DO

    ! ==========================================================================
    ! Do land cover change that requires I/O, and thus must be outside a threaded region
    ! ==========================================================================

    if (get_do_transient_pfts()) then
       call dynpft_interp(bounds_proc)
    end if

    if (get_do_transient_crops()) then
       call dyncrop_interp(bounds_proc,crop_vars)
    end if

    if (get_do_harvest()) then
       call dynHarvest_interp(bounds_proc)
    end if

    ! ==========================================================================
    ! Do everything else related to land cover change
    ! ==========================================================================

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)

       if (use_fates) then
          call dyn_ED(bounds_clump)
       end if

       if (create_glacier_mec_landunit) then
          call glc2lnd_vars%update_glc2lnd(bounds_clump)
       end if

       ! Everything following this point in this loop only needs to be called if we have
       ! actually changed some weights in this time step. This is also required in the
       ! first time step of the run to update filters to reflect state of CISM
       ! (particularly mask that is past through coupler).

       call dynSubgrid_wrapup_weight_changes(bounds_clump, glc2lnd_vars)
       call patch_state_updater%set_new_weights(bounds_clump)
       call column_state_updater%set_new_weights(bounds_clump, nc)

       call set_subgrid_diagnostic_fields(bounds_clump)

       call initialize_new_columns(bounds_clump, &
            prior_weights%cactive(bounds_clump%begc:bounds_clump%endc), &
            temperature_vars, waterstate_vars, soilhydrology_vars)

       call dyn_hwcontent_final(bounds_clump, &
            filter(nc)%num_nolakec, filter(nc)%nolakec, &
            filter(nc)%num_lakec, filter(nc)%lakec, &
            urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
            waterstate_vars, waterflux_vars, temperature_vars, energyflux_vars)

       if (use_cn) then
          call dyn_cnbal_patch(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilp, filter_inactive_and_active(nc)%soilp, &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc, &
               prior_weights, &
               patch_state_updater, &
               canopystate_vars, photosyns_vars, cnstate_vars, &
               veg_cs, c13_veg_cs, c14_veg_cs, &
               carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
               nitrogenstate_vars, nitrogenflux_vars, &
               veg_ns, &
               phosphorusstate_vars,phosphorusflux_vars, veg_ps)

          call CarbonStateUpdateDynPatch(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc)

          call NitrogenStateUpdateDynPatch(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc)

          call PhosphorusStateUpdateDynPatch(bounds_clump, &
               filter_inactive_and_active(nc)%num_soilc, filter_inactive_and_active(nc)%soilc)

          call dyn_cnbal_column(bounds_clump, nc, column_state_updater, &
               col_cs, c13_col_cs, c14_col_cs, &
               phosphorusstate_vars, col_ns, col_ps )
       end if

    end do
    !$OMP END PARALLEL DO

  end subroutine dynSubgrid_driver

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_wrapup_weight_changes(bounds_clump, glc2lnd_vars)
    !
    ! !DESCRIPTION:
    ! Reconcile various variables after subgrid weights change
    !
    ! !USES:
    use decompMod         , only : bounds_type
    use subgridWeightsMod , only : compute_higher_order_weights
    use reweightMod       , only : reweight_wrapup
    use decompMod         , only : BOUNDS_LEVEL_CLUMP
    !
    ! !ARGUMENTS:
    type(bounds_type)  , intent(in)    :: bounds_clump ! clump-level bounds
    type(glc2lnd_type) , intent(inout) :: glc2lnd_vars
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'dynSubgrid_wrapup_weight_changes'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds_clump%level == BOUNDS_LEVEL_CLUMP, subname // ': argument must be CLUMP-level bounds')

    call update_landunit_weights(bounds_clump)

    call compute_higher_order_weights(bounds_clump)

    ! Here: filters are re-made
    !
    ! This call requires clump-level bounds, which is why we need to ensure that the
    ! argument to this routine is clump-level bounds
    call reweight_wrapup(bounds_clump, glc2lnd_vars%icemask_grc(bounds_clump%begg:bounds_clump%endg))

  end subroutine dynSubgrid_wrapup_weight_changes

end module dynSubgridDriverMod
