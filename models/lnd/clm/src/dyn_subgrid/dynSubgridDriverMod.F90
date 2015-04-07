module dynSubgridDriverMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! High-level routines for dynamic subgrid areas (prescribed transient Patches, CNDV, and
  ! dynamic landunits).
  !
  ! !USES:
  use dynPriorWeightsMod           , only : prior_weights_type
  use UrbanParamsType              , only : urbanparams_type
  use CNDVType                     , only : dgvs_type
  use CanopyStateType              , only : canopystate_type
  use CNVegStateType               , only : cnveg_state_type
  use CNVegCarbonStateType         , only : cnveg_carbonstate_type
  use CNVegCarbonFluxType          , only : cnveg_carbonflux_type
  use CNVegNitrogenStateType       , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType        , only : cnveg_nitrogenflux_type
  use SoilBiogeochemStateType      , only : soilBiogeochem_state_type
  use SoilBiogeochemCarbonFluxType , only : soilBiogeochem_carbonflux_type
  use EnergyFluxType               , only : energyflux_type
  use LakeStateType                , only : lakestate_type
  use PhotosynthesisMod            , only : photosyns_type
  use SoilHydrologyType            , only : soilhydrology_type  
  use SoilStateType                , only : soilstate_type
  use WaterfluxType                , only : waterflux_type
  use WaterstateType               , only : waterstate_type
  use TemperatureType              , only : temperature_type
  use glc2lndMod                   , only : glc2lnd_type
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  !
  public :: dynSubgrid_init             ! initialize transient land cover
  public :: dynSubgrid_driver           ! top-level driver for transient land cover
  !
  ! !PRIVATE TYPES:
  type(prior_weights_type) :: prior_weights ! saved weights from before the subgrid weight updates
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_init(bounds, dgvs_inst)
    !
    ! !DESCRIPTION:
    ! Determine initial subgrid weights for prescribed transient Patches, CNDV, and/or
    ! dynamic landunits. Note that these weights will be overwritten in a restart run.
    !
    ! This should be called from initialization. 
    !
    ! Note that dynpft_init / dynpft_interp need to be called from outside any loops over
    ! clumps - so this routine needs to be called from outside any loops over clumps.
    !
    ! !USES:
    use clm_varctl        , only : flanduse_timeseries, use_cndv
    use decompMod         , only : bounds_type, BOUNDS_LEVEL_PROC
    use dynpftFileMod     , only : dynpft_init
    use dynHarvestMod     , only : dynHarvest_init
    use dynCNDVMod        , only : dynCNDV_init
    use subgridWeightsMod , only : compute_higher_order_weights
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in)    :: bounds  ! processor-level bounds
    type(dgvs_type)  , intent(inout) :: dgvs_inst
    !
    ! !LOCAL VARIABLES:
    character(len=*), parameter :: subname = 'dynSubgrid_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    prior_weights = prior_weights_type(bounds)

    ! Initialize stuff for prescribed transient Patches
    if (flanduse_timeseries /= ' ') then
       call dynpft_init(bounds)
    end if

    ! Initialize stuff for harvest (currently shares the flanduse_timeseries file)
    if (flanduse_timeseries /= ' ') then
       call dynHarvest_init(bounds)
    end if

    if (use_cndv) then
       call dynCNDV_init(bounds, dgvs_inst)
    end if

    call compute_higher_order_weights(bounds)
    
  end subroutine dynSubgrid_init

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_driver(bounds_proc,                                            &
       urbanparams_inst, soilstate_inst, soilhydrology_inst, lakestate_inst,           &
       waterstate_inst, waterflux_inst, temperature_inst, energyflux_inst,             &
       canopystate_inst, photosyns_inst, dgvs_inst, glc2lnd_inst, cnveg_state_inst,    &
       cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst, &
       cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,    &
       cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst,                              &
       soilbiogeochem_state_inst, soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! Update subgrid weights for prescribed transient Patches, CNDV, and/or dynamic
    ! landunits. Also do related adjustments (water & energy, carbon & nitrogen).
    !
    ! This should be called every time step in CLM's run loop.
    !
    ! Note that this routine operates partly at the proc-level (outside an OMP region),
    ! and partly at the clump level (inside OMP regions). Thus, this must be called from
    ! OUTSIDE any loops over clumps in the driver.
    !
    ! !USES:
    use clm_varctl           , only : flanduse_timeseries, use_cndv, use_cn, create_glacier_mec_landunit, use_ed
    use decompMod            , only : bounds_type, get_proc_clumps, get_clump_bounds
    use decompMod            , only : BOUNDS_LEVEL_PROC
    use dynLandunitAreaMod   , only : update_landunit_weights
    use dynInitColumnsMod    , only : initialize_new_columns
    use dynConsBiogeophysMod , only : dyn_hwcontent_init, dyn_hwcontent_final
    use dynConsBiogeochemMod , only : dyn_cnbal_patch
    use dynpftFileMod        , only : dynpft_interp
    use dynHarvestMod        , only : dynHarvest_interp
    use dynCNDVMod           , only : dynCNDV_interp
    use dynEDMod             , only : dyn_ED
    use reweightMod          , only : reweight_wrapup
    use subgridWeightsMod    , only : compute_higher_order_weights, set_subgrid_diagnostic_fields
    !
    ! !ARGUMENTS:
    type(bounds_type)                    , intent(in)    :: bounds_proc  ! processor-level bounds
    type(urbanparams_type)               , intent(in)    :: urbanparams_inst
    type(soilstate_type)                 , intent(in)    :: soilstate_inst
    type(soilhydrology_type)             , intent(in)    :: soilhydrology_inst
    type(lakestate_type)                 , intent(in)    :: lakestate_inst
    type(waterstate_type)                , intent(inout) :: waterstate_inst
    type(waterflux_type)                 , intent(inout) :: waterflux_inst
    type(temperature_type)               , intent(inout) :: temperature_inst
    type(energyflux_type)                , intent(inout) :: energyflux_inst
    type(canopystate_type)               , intent(inout) :: canopystate_inst
    type(photosyns_type)                 , intent(inout) :: photosyns_inst
    type(dgvs_type)                      , intent(inout) :: dgvs_inst
    type(glc2lnd_type)                   , intent(inout) :: glc2lnd_inst
    type(cnveg_state_type)               , intent(inout) :: cnveg_state_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: c13_cnveg_carbonstate_inst
    type(cnveg_carbonstate_type)         , intent(inout) :: c14_cnveg_carbonstate_inst
    type(cnveg_carbonflux_type)          , intent(inout) :: cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)          , intent(inout) :: c13_cnveg_carbonflux_inst
    type(cnveg_carbonflux_type)          , intent(inout) :: c14_cnveg_carbonflux_inst
    type(cnveg_nitrogenstate_type)       , intent(inout) :: cnveg_nitrogenstate_inst
    type(cnveg_nitrogenflux_type)        , intent(inout) :: cnveg_nitrogenflux_inst
    type(soilbiogeochem_state_type)      , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
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
            urbanparams_inst, soilstate_inst, soilhydrology_inst, lakestate_inst, &
            waterstate_inst, waterflux_inst, temperature_inst, energyflux_inst)

       call prior_weights%set_prior_weights(bounds_clump)
    end do
    !$OMP END PARALLEL DO

    ! ==========================================================================
    ! Do land cover change that requires I/O, and thus must be outside a threaded region
    ! ==========================================================================

    if (flanduse_timeseries /= ' ') then
       call dynpft_interp(bounds_proc)
    end if

    if (flanduse_timeseries /= ' ') then
       call dynHarvest_interp(bounds_proc)
    end if

    ! ==========================================================================
    ! Do everything else related to land cover change
    ! ==========================================================================

    !$OMP PARALLEL DO PRIVATE (nc, bounds_clump)
    do nc = 1, nclumps
       call get_clump_bounds(nc, bounds_clump)

       if (use_cndv) then
          call dynCNDV_interp(bounds_clump, dgvs_inst)
       end if
       
       if (use_ed) then
          call dyn_ED(bounds_clump)
       end if

       if (create_glacier_mec_landunit) then
          call glc2lnd_inst%update_glc2lnd(bounds_clump)
       end if

       ! Everything following this point in this loop only needs to be called if we have
       ! actually changed some weights in this time step. This is also required in the
       ! first time step of the run to update filters to reflect state of CISM
       ! (particularly mask that is past through coupler).

       call update_landunit_weights(bounds_clump)

       call compute_higher_order_weights(bounds_clump)

       ! Here: filters are re-made
       call reweight_wrapup(bounds_clump, &
            glc2lnd_inst%icemask_grc(bounds_clump%begg:bounds_clump%endg))

       call set_subgrid_diagnostic_fields(bounds_clump)

       call initialize_new_columns(bounds_clump, &
            prior_weights%cactive(bounds_clump%begc:bounds_clump%endc), &
            temperature_inst)

       call dyn_hwcontent_final(bounds_clump, &
            urbanparams_inst, soilstate_inst, soilhydrology_inst, lakestate_inst, &
            waterstate_inst, waterflux_inst, temperature_inst, energyflux_inst)

       if (use_cn) then
          call dyn_cnbal_patch(bounds_clump, prior_weights,                                       &
               canopystate_inst, photosyns_inst, cnveg_state_inst,                                &
               cnveg_carbonstate_inst, c13_cnveg_carbonstate_inst, c14_cnveg_carbonstate_inst,    &
               cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,       &
               cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, soilbiogeochem_carbonflux_inst, &
               soilbiogeochem_state_inst)
       end if

    end do
    !$OMP END PARALLEL DO

  end subroutine dynSubgrid_driver

end module dynSubgridDriverMod
