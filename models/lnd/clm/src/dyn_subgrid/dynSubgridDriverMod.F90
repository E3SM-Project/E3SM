module dynSubgridDriverMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  !
  ! !DESCRIPTION:
  ! High-level routines for dynamic subgrid areas (prescribed transient Patches, CNDV, and
  ! dynamic landunits).
  !
  ! !USES:
  use dynPriorWeightsMod  , only : prior_weights_type
  use UrbanParamsType     , only : urbanparams_type
  use CNDVType            , only : dgvs_type
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
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  private
  save

  public :: dynSubgrid_init             ! initialize transient land cover
  public :: dynSubgrid_driver           ! top-level driver for transient land cover
  !
  ! !PRIVATE TYPES:
  type(prior_weights_type) :: prior_weights ! saved weights from before the subgrid weight updates
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_init(bounds, dgvs_vars)
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
    type(dgvs_type)  , intent(inout) :: dgvs_vars
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
       call dynCNDV_init(bounds, dgvs_vars)
    end if

    call compute_higher_order_weights(bounds)
    
  end subroutine dynSubgrid_init

  !-----------------------------------------------------------------------
  subroutine dynSubgrid_driver(bounds_proc, &
       urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
       waterstate_vars, waterflux_vars, temperature_vars, energyflux_vars, &
       canopystate_vars, photosyns_vars, cnstate_vars, dgvs_vars, &
       carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, &
       carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
       nitrogenstate_vars, nitrogenflux_vars, glc2lnd_vars)
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
    type(bounds_type)        , intent(in)    :: bounds_proc  ! processor-level bounds
    type(urbanparams_type)   , intent(in)    :: urbanparams_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(lakestate_type)     , intent(in)    :: lakestate_vars
    type(waterstate_type)    , intent(inout) :: waterstate_vars
    type(waterflux_type)     , intent(inout) :: waterflux_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(energyflux_type)    , intent(inout) :: energyflux_vars
    type(canopystate_type)   , intent(inout) :: canopystate_vars
    type(photosyns_type)     , intent(inout) :: photosyns_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(dgvs_type)          , intent(inout) :: dgvs_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars
    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(glc2lnd_type)       , intent(inout) :: glc2lnd_vars
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
            urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
            waterstate_vars, waterflux_vars, temperature_vars, energyflux_vars)

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
          call dynCNDV_interp(bounds_clump, dgvs_vars)
       end if
       
       if (use_ed) then
          call dyn_ED(bounds_clump)
       end if

       if (create_glacier_mec_landunit) then
          call glc2lnd_vars%update_glc2lnd(bounds_clump)
       end if

       ! Everything following this point in this loop only needs to be called if we have
       ! actually changed some weights in this time step. This is also required in the
       ! first time step of the run to update filters to reflect state of CISM
       ! (particularly mask that is past through coupler).

       call update_landunit_weights(bounds_clump)

       call compute_higher_order_weights(bounds_clump)

       ! Here: filters are re-made
       call reweight_wrapup(bounds_clump, &
            glc2lnd_vars%icemask_grc(bounds_clump%begg:bounds_clump%endg))

       call set_subgrid_diagnostic_fields(bounds_clump)

       call initialize_new_columns(bounds_clump, &
            prior_weights%cactive(bounds_clump%begc:bounds_clump%endc), &
            temperature_vars)

       call dyn_hwcontent_final(bounds_clump, &
            urbanparams_vars, soilstate_vars, soilhydrology_vars, lakestate_vars, &
            waterstate_vars, waterflux_vars, temperature_vars, energyflux_vars)

       if (use_cn) then
          call dyn_cnbal_patch(bounds_clump, prior_weights, &
               canopystate_vars, photosyns_vars, cnstate_vars, &
               carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, &
               carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
               nitrogenstate_vars, nitrogenflux_vars)
       end if

    end do
    !$OMP END PARALLEL DO

  end subroutine dynSubgrid_driver

end module dynSubgridDriverMod
