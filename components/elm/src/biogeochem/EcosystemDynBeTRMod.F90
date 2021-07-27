module EcosystemDynBeTRMod

  !
  ! DESCRIPTION
  ! betr based aboveground belowground coupling
  !
  ! Created by Jinyun Tang
  ! Now it is only for generic carbon coupling no isotope is attempted below, but will
  ! be enabled gradually.
  use shr_kind_mod              , only : r8 => shr_kind_r8
  use shr_sys_mod               , only : shr_sys_flush
  use elm_varctl                , only : use_c13, use_c14, use_fates, use_dynroot
  use elm_varctl                , only : use_pheno_flux_limiter, iulog, use_erosion
  use elm_varpar                , only : nlevsoi
  use elm_varctl                , only : use_erosion
  use decompMod                 , only : bounds_type
  use perf_mod                  , only : t_startf, t_stopf
  use spmdMod                   , only : masterproc
  use elm_varctl                , only : use_century_decomp
  use CNStateType               , only : cnstate_type
  use CNCarbonFluxType          , only : carbonflux_type
  use CNCarbonStateType         , only : carbonstate_type
  use CNNitrogenFluxType        , only : nitrogenflux_type
  use CNNitrogenStateType       , only : nitrogenstate_type
  use CanopyStateType           , only : canopystate_type
  use SoilStateType             , only : soilstate_type
  use TemperatureType           , only : temperature_type
  use WaterstateType            , only : waterstate_type
  use WaterfluxType             , only : waterflux_type
  use atm2lndType               , only : atm2lnd_type
  use CanopyStateType           , only : canopystate_type
  use PhotosynthesisType        , only : photosyns_type
  use CH4Mod                    , only : ch4_type
  use EnergyFluxType            , only : energyflux_type
  use SoilHydrologyType         , only : soilhydrology_type
  use FrictionVelocityType      , only : frictionvel_type
  use tracerfluxType            , only : tracerflux_type
  use tracerstatetype           , only : tracerstate_type
  use BetrTracerType            , only : betrtracer_type
  use PhosphorusFluxType        , only : phosphorusflux_type
  use PhosphorusStateType       , only : phosphorusstate_type
  use dynSubgridControlMod      , only : get_do_harvest
  use clm_time_manager          , only : get_nstep
  use ColumnType                , only : column_physical_properties
  use VegetationType            , only : vegetation_physical_properties
  use ColumnDataType            , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType            , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType            , only : col_ns, col_nf, col_ws
  use ColumnDataType            , only : col_ps, col_pf, col_es
  use VegetationDataType        , only : veg_cs, c13_veg_cs, c14_veg_cs
  use VegetationDataType        , only : veg_cf, c13_veg_cf, c14_veg_cf
  use VegetationDataType        , only : veg_ns, veg_nf
  use VegetationDataType        , only : veg_ps, veg_pf
  use AllocationMod             , only : nu_com_nfix, nu_com_phosphatase
  use SoilLittDecompMod         , only : SoilLittDecompAlloc
  use SoilLittDecompMod         , only : SoilLittDecompAlloc2 !after SoilLittDecompAlloc
  use timeinfoMod
  use VegetationDataType        , only : veg_cf_summary, veg_cf_summary_for_ch4, veg_cf_summary_rr
  use VegetationDataType        , only : veg_nf_summary, veg_ns_summary, veg_cs_Summary
  use VegetationDataType        , only : veg_pf_summary, veg_ps_summary
  use VegetationDataType        , only : veg_cf_setvalues, veg_nf_setvalues, veg_pf_setvalues

  use ColumnDataType            , only : col_cf_summary, col_nf_summary, col_pf_Summary
  use ColumnDataType            , only : col_cs_summary, col_ns_summary, col_ps_summary
  use ColumnDataType            , only : col_cf_summary_for_ch4
  use ColumnDataType            , only : col_cf_setvalues, col_nf_setvalues, col_pf_setvalues
  use timeinfoMod
  use perfMod_GPU
  use ErosionMod                , only : ErosionFluxes
  implicit none

  private
  public :: CNEcosystemDynBeTR0
  public :: CNEcosystemDynBeTR1
  public :: CNFluxStateBeTR0Summary
  public :: CNFluxStateBeTR1Summary
  contains

  !-----------------------------------------------------------------------
  subroutine CNEcosystemDynBeTR0(bounds,                             &
         num_soilc, filter_soilc,                                        &
         num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb,      &
         cnstate_vars, atm2lnd_vars, canopystate_vars, soilstate_vars,   &
         crop_vars, photosyns_vars, soilhydrology_vars, energyflux_vars, &
         PlantMicKinetics_vars, ch4_vars, sedflux_vars)

    ! Description:
    ! Update vegetation related state variables and
    ! setup fluxes and parameters for plant-microbe coupling in soibgc
    !
    ! !USES:
    use PhenologyMod                     , only : Phenology, CNLitterToColumn
    use GrowthRespMod                    , only : GrowthResp
    use GapMortalityMod                  , only : GapMortality
    use CarbonStateUpdate1BeTRMod        , only : CarbonStateUpdate1,CarbonStateUpdate0
    use CarbonStateUpdate2BeTRMod        , only : CarbonStateUpdate2, CarbonStateUpdate2h
    use CarbonStateUpdate3BeTRMod        , only : CarbonStateUpdate3
    use NitrogenStateUpdate1BeTRMod      , only : NitrogenStateUpdate1
    use NitrogenStateUpdate2BeTRMod      , only : NitrogenStateUpdate2, NitrogenStateUpdate2h
    use NitrogenStateUpdate3BeTRMod      , only : NitrogenStateUpdate3
    use PhosphorusStateUpdate1BeTRMod    , only : PhosphorusStateUpdate1
    use PhosphorusStateUpdate2BeTRMod    , only : PhosphorusStateUpdate2, PhosphorusStateUpdate2h
    use PhosphorusStateUpdate3BeTRMod    , only : PhosphorusStateUpdate3
    use FireMod                          , only: FireArea, FireFluxes
    use CarbonIsoFluxMod                 , only : CarbonIsoFlux1, CarbonIsoFlux2, CarbonIsoFlux2h, CarbonIsoFlux3
    use C14DecayMod                      , only : C14Decay, C14BombSpike
    use ErosionMod                       , only : ErosionFluxes
    use WoodProductsMod                  , only : WoodProducts
    use CropType                         , only : crop_type
    use dynHarvestMod                    , only : CNHarvest
    use elm_varpar                       , only : crop_prog
    use NitrogenDynamicsMod              , only : NitrogenLeaching
    use CropHarvestPoolsMod              , only : CropHarvestPools
    use PlantMicKineticsMod              , only : PlantMicKinetics_type
    use PhosphorusDynamicsMod            , only : PhosphorusDeposition,PhosphorusWeathering,PhosphorusAdsportion
    use PhosphorusDynamicsMod            , only : PhosphorusLeaching, PhosphorusOcclusion,PhosphorusDesoprtion
    use VerticalProfileMod               , only : decomp_vertprofiles
    use SoilLittVertTranspMod            , only : SoilLittVertTransp
    use RootDynamicsMod                  , only : RootDynamics
    use PhenologyFLuxLimitMod            , only : phenology_flux_limiter
    use abortutils                       , only : endrun
    use shr_log_mod                      , only : errMsg => shr_log_errMsg
    use EcosystemDynMod                  , only : EcosystemDynNoLeaching1
    use SedFluxType                      , only : sedflux_type
    implicit none


    !
    ! !ARGUMENTS:
    type(bounds_type)                , intent(in)    :: bounds
    integer                          , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                          , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                          , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                          , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    integer                          , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                          , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    logical                          , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(cnstate_type)               , intent(inout) :: cnstate_vars
    type(atm2lnd_type)               , intent(in)    :: atm2lnd_vars
    type(canopystate_type)           , intent(in)    :: canopystate_vars
    type(soilstate_type)             , intent(inout) :: soilstate_vars
    type(crop_type)                  , intent(inout) :: crop_vars
    type(ch4_type)                    , intent(in)    :: ch4_vars
    type(photosyns_type)             , intent(in)    :: photosyns_vars
    type(soilhydrology_type)         , intent(in)    :: soilhydrology_vars
    type(energyflux_type)            , intent(in)    :: energyflux_vars
    type(PlantMicKinetics_type)      , intent(inout) :: PlantMicKinetics_vars
    type(sedflux_type)               , intent(in)    :: sedflux_vars

    real(r8) :: dt
    integer :: c13, c14
    c13 = 0
    c14 = 1

    dt = dtime_mod

    call EcosystemDynNoLeaching1(bounds,                                  &
         num_soilc, filter_soilc,                                            &
         num_soilp, filter_soilp,                            &
         cnstate_vars, atm2lnd_vars,               &
         canopystate_vars, soilstate_vars, crop_vars,   &
         ch4_vars, photosyns_vars)


    ! Call the main CN routines
    ! only do if ed is off
    if( .not. use_fates ) then

       call t_startf('SoilLittDecompAlloc')
       !----------------------------------------------------------------
            ! directly run clm-bgc
            ! if (use_clm_interface & use_clm_bgc), then CNDecomAlloc is called in clm_driver
       call SoilLittDecompAlloc (bounds, num_soilc, filter_soilc,    &
                       num_soilp, filter_soilp,                     &
                       canopystate_vars, soilstate_vars,            &
                       cnstate_vars, ch4_vars,                      &
                       dt)
       !----------------------------------------------------------------
       ! SoilLittDecompAlloc2 is called by both clm-bgc & pflotran
       ! pflotran: call 'SoilLittDecompAlloc2' to calculate some diagnostic variables and 'fpg' for plant N uptake
       ! pflotran & clm-bgc : 'Allocation3_AG' and vertically integrate net and gross mineralization fluxes
       call SoilLittDecompAlloc2 (bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,           &
                photosyns_vars, canopystate_vars, soilstate_vars,   &
                cnstate_vars, ch4_vars,                                        &
                crop_vars, atm2lnd_vars,                 &
                dt)

       !----------------------------------------------------------------
       call t_stopf('SoilLittDecompAlloc')
       !----------------------------------------------------------------

       !--------------------------------------------
       ! Phenology
       !--------------------------------------------

       ! Phenology needs to be called after SoilLittDecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call t_startf('Phenology')
       call Phenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            num_pcropp, filter_pcropp, doalb, atm2lnd_vars, &
            crop_vars, canopystate_vars, soilstate_vars, &
            cnstate_vars )
       call t_stopf('Phenology')

       !--------------------------------------------
       ! Growth respiration
       !--------------------------------------------

       call t_startf('GrowthResp')
       call GrowthResp(num_soilp, filter_soilp)
       call t_stopf('GrowthResp')

       call veg_cf_summary_rr(veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, col_cf)
       if(use_c13) then
         call veg_cf_summary_rr(c13_veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, c13_col_cf)
       endif
       if(use_c14) then
         call veg_cf_summary_rr(c14_veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, c14_col_cf)
       endif

       !--------------------------------------------
       ! Dynamic Roots
       !--------------------------------------------

       if( use_dynroot ) then
          call t_startf('RootDynamics')

          call RootDynamics(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               canopystate_vars,  &
               cnstate_vars, crop_vars, energyflux_vars, soilstate_vars)
          call t_stopf('RootDynamics')
       end if

       !--------------------------------------------
       ! CNUpdate0
       !--------------------------------------------

       call t_startf('CNUpdate0')
       call CarbonStateUpdate0(num_soilp, filter_soilp, veg_cs, veg_cf)
       if ( use_c13 ) then
          call CarbonStateUpdate0(num_soilp, filter_soilp, c13_veg_cs, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate0(num_soilp, filter_soilp, c14_veg_cs, c14_veg_cf)
       end if
       call t_stopf('CNUpdate0')

       !--------------------------------------------
       if(use_pheno_flux_limiter)then
         call t_startf('phenology_flux_limiter')
         call phenology_flux_limiter(bounds, num_soilc, filter_soilc,&
           num_soilp, filter_soilp, crop_vars, cnstate_vars,  &
           veg_cf, veg_cs, &
           c13_veg_cf, c13_veg_cs, &
           c14_veg_cf, c14_veg_cs, &
           veg_nf, veg_ns, veg_pf, veg_ps)
         call t_stopf('phenology_flux_limiter')
       endif
       call t_startf('CNLitterToColumn')
       call CNLitterToColumn(num_soilc, filter_soilc,  cnstate_vars)

       call t_stopf('CNLitterToColumn')
       !--------------------------------------------
       ! Update1
       !--------------------------------------------

       call t_startf('CNUpdate1')

       if ( use_c13 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, &
               isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if

       if ( use_c14 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars , &
               isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            crop_vars, col_cs, veg_cs, col_cf, veg_cf, dt)

       if ( use_c13 ) then
          call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf, dt)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, dt)
       end if

       call NitrogenStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, dt)

       call PhosphorusStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, dt)

       call t_stopf('CNUpdate1')

       call t_startf('SoilLittVertTransp')
       call SoilLittVertTransp(bounds, &
            num_soilc, filter_soilc, &
            canopystate_vars, cnstate_vars)
       call t_stopf('SoilLittVertTransp')

       call t_startf('CNGapMortality')
       call GapMortality( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars)
       call t_stopf('CNGapMortality')

       !--------------------------------------------
       ! Update2
       !--------------------------------------------

       call t_startf('CNUpdate2')

       if ( use_c13 ) then
           call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars,  &
                isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if

       if ( use_c14 ) then
           call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                cnstate_vars, &
                isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate2( num_soilc, filter_soilc, num_soilp, filter_soilp, &
             col_cs, veg_cs, col_cf, veg_cf, dt)

       if ( use_c13 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf, dt)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, dt)
       end if
       call NitrogenStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)

       call PhosphorusStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp,dt)

       if (get_do_harvest()) then
          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars)
       end if

       if ( use_c13 ) then
          call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, &
               isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars , &
               isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate2h( num_soilc, filter_soilc,  num_soilp, filter_soilp, &
             col_cs, veg_cs, col_cf, veg_cf)
       if ( use_c13 ) then
          call CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf)
       end if

       call NitrogenStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ns, veg_ns, col_nf, veg_nf)

       call PhosphorusStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ps, veg_ps, col_pf, veg_pf)

       call WoodProducts(num_soilc, filter_soilc)

       call CropHarvestPools(num_soilc, filter_soilc, dt)

       call FireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            atm2lnd_vars, energyflux_vars, soilhydrology_vars,&
            cnstate_vars)

       call FireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars)

       if ( use_erosion ) then
            call ErosionFluxes(bounds, num_soilc, filter_soilc, soilstate_vars, sedflux_vars)
       end if

       call t_stopf('CNUpdate2')

       !--------------------------------------------
       ! Update3
       !--------------------------------------------

       if ( use_c13 ) then
          call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnstate_vars, &
             isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars , &
               isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_cs, veg_cs, col_cf, veg_cf, dt)

       if ( use_c13 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf, dt)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, dt)
       end if


       if ( use_c14 ) then
          call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars)

          call C14BombSpike(num_soilp, filter_soilp, &
               cnstate_vars)
       end if

       call veg_cf_summary_for_ch4(veg_cf,bounds, num_soilp, filter_soilp)
       if( use_c13 ) then
          call col_cf_summary_for_ch4(c13_col_cf,bounds, num_soilc, filter_soilc)
          call veg_cf_summary_for_ch4(c13_veg_cf,bounds, num_soilp, filter_soilp)
       endif
       if( use_c14 ) then
          call col_cf_summary_for_ch4(c14_col_cf,bounds, num_soilc, filter_soilc)
          call veg_cf_summary_for_ch4(c14_veg_cf,bounds, num_soilp, filter_soilp)
       endif
    end if !end of if not use_fates block

    call col_cf_summary_for_ch4(col_cf,bounds, num_soilc, filter_soilc)
  end subroutine CNEcosystemDynBeTR0

  !-----------------------------------------------------------------------
  subroutine CNEcosystemDynBeTR1(bounds, col, pft,                       &
         num_soilc, filter_soilc,                                        &
         num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb,      &
         cnstate_vars, atm2lnd_vars, canopystate_vars, soilstate_vars,   &
         crop_vars, photosyns_vars, soilhydrology_vars, energyflux_vars, &
         PlantMicKinetics_vars, ch4_vars,sedflux_vars,                   &
         ep_betr, soil_water_retention_curve)

    ! Description:
    ! Update vegetation related state variables and
    ! setup fluxes and parameters for plant-microbe coupling in soibgc
    !
    ! !USES:
    use PhenologyMod                      , only : Phenology, CNLitterToColumn
    use GrowthRespMod                     , only : GrowthResp
    use CarbonStateUpdate1BeTRMod         , only : CarbonStateUpdate1,CarbonStateUpdate0
    use CarbonStateUpdate2BeTRMod         , only : CarbonStateUpdate2, CarbonStateUpdate2h
    use CarbonStateUpdate3BeTRMod         , only : CarbonStateUpdate3
    use NitrogenStateUpdate1BeTRMod       , only : NitrogenStateUpdate1
    use NitrogenStateUpdate2BeTRMod       , only : NitrogenStateUpdate2, NitrogenStateUpdate2h
    use NitrogenStateUpdate3BeTRMod       , only : NitrogenStateUpdate3Veg
    use PhosphorusStateUpdate1BeTRMod     , only : PhosphorusStateUpdate1
    use PhosphorusStateUpdate2BeTRMod     , only : PhosphorusStateUpdate2, PhosphorusStateUpdate2h
    use PhosphorusStateUpdate3BeTRMod     , only : PhosphorusStateUpdate3Veg
    use GapMortalityMod                   , only : GapMortality
    use FireMod                           , only : FireArea, FireFluxes
    use CarbonIsoFluxMod                  , only : CarbonIsoFlux1, CarbonIsoFlux2, CarbonIsoFlux2h, CarbonIsoFlux3
    use C14DecayMod                       , only : C14Decay, C14BombSpike
    use WoodProductsMod                   , only : WoodProducts
    use CropType                          , only : crop_type
    use dynHarvestMod                     , only : CNHarvest
    use elm_varpar                        , only : crop_prog
    use CropHarvestPoolsMod               , only : CropHarvestPools
    use PlantMicKineticsMod               , only : PlantMicKinetics_type
    use AllocationMod                     , only : update_PlantMicKinetics_pars,Allocation3_PlantCNPAlloc
    use PhosphorusDynamicsMod             , only : PhosphorusDeposition,PhosphorusWeathering,PhosphorusLeaching
    use PhosphorusDynamicsMod             , only : PhosphorusAdsportion,PhosphorusDesoprtion,PhosphorusOcclusion
    use NitrogenDynamicsMod               , only : NitrogenLeaching
    use VerticalProfileMod                , only : decomp_vertprofiles
    use RootDynamicsMod                   , only : RootDynamics
    use PhenologyFLuxLimitMod             , only : phenology_flux_limiter
    use abortutils                        , only : endrun
    use shr_log_mod                       , only : errMsg => shr_log_errMsg
    use EcosystemDynMod                   , only : EcosystemDynNoLeaching1
    use SoilLittVertTranspMod             , only : SoilLittVertTransp
    use SoilLittDecompMod                 , only : SoilLittDecompAlloc2
    use BeTRSimulationALM                 , only : betr_simulation_alm_type
    use SoilWaterRetentionCurveMod        , only : soil_water_retention_curve_type
    use subgridAveMod                     , only : p2c
    use AllocationMod                     , only : calc_nfix_stress
    use SedFluxType         , only : sedflux_type
    implicit none


    !
    ! !ARGUMENTS:
    type(bounds_type)                , intent(in)    :: bounds
    type(column_physical_properties)        , intent(in)    :: col
    type(vegetation_physical_properties)         , intent(in)    :: pft
    integer                          , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                          , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                          , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                          , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    integer                          , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                          , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    logical                          , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(cnstate_type)               , intent(inout) :: cnstate_vars
    type(atm2lnd_type)               , intent(in)    :: atm2lnd_vars
    type(canopystate_type)           , intent(in)    :: canopystate_vars
    type(soilstate_type)             , intent(inout) :: soilstate_vars
    type(crop_type)                  , intent(inout) :: crop_vars
    type(ch4_type)                   , intent(in)    :: ch4_vars
    type(photosyns_type)             , intent(in)    :: photosyns_vars
    type(soilhydrology_type)         , intent(in)    :: soilhydrology_vars
    type(energyflux_type)            , intent(in)    :: energyflux_vars
    type(PlantMicKinetics_type)      , intent(inout) :: PlantMicKinetics_vars
    type(sedflux_type)       , intent(in)    :: sedflux_vars
    class(betr_simulation_alm_type)  , intent(inout) :: ep_betr
    class(soil_water_retention_curve_type) , intent(in) :: soil_water_retention_curve
    character(len=64) :: event
    real(r8) :: dt
    integer :: c13, c14
    c13 = 0
    c14 = 1

    dt = dtime_mod

    call t_startf('EcosystemDynNoLeaching1')
    call EcosystemDynNoLeaching1(bounds,                                  &
         num_soilc, filter_soilc,                                         &
         num_soilp, filter_soilp,                                         &
         cnstate_vars,                  &
         atm2lnd_vars,                   &
         canopystate_vars, soilstate_vars,  crop_vars,   &
         ch4_vars, photosyns_vars)
    call t_stopf('EcosystemDynNoLeaching1')


    !----------------------------------------------------------------
    !call decomposition method from betr
    !----------------------------------------------------------------
    call t_startf('betr type1 soil bgc')

    call update_PlantMicKinetics_pars(bounds, num_soilc, filter_soilc  , &
        cnstate_vars, PlantMicKinetics_vars)

    call ep_betr%CalcSmpL(bounds, 1, nlevsoi, num_soilc, filter_soilc, &
            col_es%t_soisno(bounds%begc:bounds%endc,1:nlevsoi), &
            soilstate_vars, col_ws, soil_water_retention_curve)

    col_ws%finundated(filter_soilc(1:num_soilc))=0._r8


    call veg_cf_summary_rr(veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, col_cf)
    if(use_c13) then
       call veg_cf_summary_rr(c13_veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, c13_col_cf)
    endif
    if(use_c14) then
       call veg_cf_summary_rr(c14_veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, c14_col_cf)
    endif


    call ep_betr%SetBiophysForcing(bounds, col, pft,  &
         waterstate_vars=col_ws, temperature_vars=col_es,&
         atm2lnd_vars=atm2lnd_vars, soilstate_vars=soilstate_vars, carbonflux_vars=col_cf)

    !pass in parameters into betr
    call ep_betr%EnterOutLoopBGC(bounds, col, pft, &
       num_soilc, filter_soilc, &
       col_cs, col_cf, &
       c13_col_cs, c14_col_cs,  &
       col_ns,  col_ps, PlantMicKinetics_vars)

    !do bgc reaction inside betr
    call ep_betr%OutLoopSoilBGC(bounds, col, pft)

    !gathering outputs from betr
    call ep_betr%ExitOutLoopBGC(bounds, col, pft, &
        col_cs, col_cf, &
        c13_col_cs, c13_col_cf, &
        c14_col_cs, c14_col_cf, &
        col_ns, col_nf, veg_nf, &
        col_ps, col_pf, veg_pf)

   !To avoid erroneous zero out, the following two calls should not be bundled
   !with the summary subroutine
    call p2c(bounds, num_soilc, filter_soilc, &
         veg_nf%sminn_to_plant(bounds%begp:bounds%endp), &
         col_nf%sminn_to_plant(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         veg_pf%sminp_to_plant(bounds%begp:bounds%endp), &
         col_pf%sminp_to_plant(bounds%begc:bounds%endc))

    call calc_nfix_stress(num_soilc, filter_soilc, cnstate_vars, col_cf, veg_cs, veg_ns, veg_nf)

    !resolve nutrient allocation
    call Allocation3_PlantCNPAlloc (bounds            , &
        num_soilc, filter_soilc, num_soilp, filter_soilp    , &
        canopystate_vars                                    , &
        cnstate_vars, crop_vars, dt)
       !--------------------------------------------
       ! Phenology
       !--------------------------------------------

       ! Phenology needs to be called after SoilLittDecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

     call t_startf('Phenology')

        call Phenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             num_pcropp, filter_pcropp, doalb, atm2lnd_vars, &
             crop_vars, canopystate_vars, soilstate_vars, &
             cnstate_vars )
     call t_stopf('Phenology')

     !--------------------------------------------
     ! Growth respiration
     !--------------------------------------------

     call t_startf('GrowthResp')
     call GrowthResp(num_soilp, filter_soilp)
     call t_stopf('GrowthResp')

     call veg_cf%SummaryRR(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, col_cf, cnstate_vars)

     if(use_c13) then
       call c13_veg_cf%SummaryRR(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, c13_col_cf)
     endif

     if(use_c14) then
       call c14_veg_cf%SummaryRR(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, c14_col_cf)
     endif
     !--------------------------------------------
     ! Dynamic Roots
     !--------------------------------------------

     if( use_dynroot ) then
        call t_startf('RootDynamics')

        call RootDynamics(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               canopystate_vars, &
               cnstate_vars, crop_vars, energyflux_vars, soilstate_vars)
        call t_stopf('RootDynamics')
     end if

     !--------------------------------------------
     ! CNUpdate0
     !--------------------------------------------

     call t_startf('CNUpdate0')
     call CarbonStateUpdate0(&
            num_soilp, filter_soilp, &
            veg_cs, veg_cf)

     if ( use_c13 ) then
        call CarbonStateUpdate0(&
               num_soilp, filter_soilp, &
               c13_veg_cs, c13_veg_cf)
     end if

     if ( use_c14 ) then
        call CarbonStateUpdate0(&
               num_soilp, filter_soilp, &
               c14_veg_cs, c14_veg_cf)
     end if
     call t_stopf('CNUpdate0')
     !--------------------------------------------
     if(use_pheno_flux_limiter)then
       call t_startf('phenology_flux_limiter')
       call phenology_flux_limiter(bounds, num_soilc, filter_soilc,&
           num_soilp, filter_soilp, crop_vars, cnstate_vars,  &
           veg_cf, veg_cs , c13_veg_cf, c13_veg_cs , c14_veg_cf, c14_veg_cs , &
           veg_nf, veg_ns, veg_pf, veg_ps)
     endif
     call CNLitterToColumn(num_soilc, filter_soilc,  cnstate_vars)

     call t_stopf('phenology_flux_limiter')
     !--------------------------------------------
     ! Update1
     !--------------------------------------------

     call t_startf('CNUpdate1')


     if ( use_c13 ) then
       call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, &
            isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
     end if

     if ( use_c14 ) then
       call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars , &
            isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
     end if


     call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            crop_vars, col_cs, veg_cs, col_cf, veg_cf, dt, ldecomp_on=.false.)

     if ( use_c13 ) then
        call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf, dt, ldecomp_on=.false.)
     end if
     if ( use_c14 ) then
        call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, dt, ldecomp_on=.false.)
     end if

     call NitrogenStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, col_ns, veg_ns, col_nf, veg_nf, dt, ldecomp_on=.false.)

     call PhosphorusStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, col_ps, veg_ps, col_pf, veg_pf, dt, ldecomp_on=.false.)
     call t_stopf('CNUpdate1')

     call t_startf('SoilLittVertTransp')
     call SoilLittVertTransp(bounds, &
            num_soilc, filter_soilc, &
            canopystate_vars, cnstate_vars)
     call t_stopf('SoilLittVertTransp')

     call t_startf('CNGapMortality')
     call GapMortality( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars)
     call t_stopf('CNGapMortality')

     !--------------------------------------------
     ! Update2
     !--------------------------------------------

     call t_startf('CNUpdate2')

    if ( use_c13 ) then
        call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnstate_vars,  &
             isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
    end if

    if ( use_c14 ) then
        call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnstate_vars, &
             isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
     end if

     call CarbonStateUpdate2( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_cs, veg_cs, col_cf, veg_cf, dt)

     if ( use_c13 ) then
        call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf, dt)
     end if
     if ( use_c14 ) then
        call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, dt)
     end if
     call NitrogenStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)


     call PhosphorusStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp,dt)

     if (get_do_harvest()) then
        call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars)
     end if

     if ( use_c13 ) then
        call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnstate_vars, &
             isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
     end if
     if ( use_c14 ) then
        call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
             cnstate_vars , &
             isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
     end if

     call CarbonStateUpdate2h( num_soilc, filter_soilc,  num_soilp, filter_soilp, &
            col_cs, veg_cs, col_cf, veg_cf)
     if ( use_c13 ) then
        call CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf)
     end if
     if ( use_c14 ) then
        call CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf)
     end if

     call NitrogenStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ns, veg_ns, col_nf, veg_nf)

     call PhosphorusStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ps, veg_ps, col_pf, veg_pf)

     call WoodProducts(num_soilc, filter_soilc)

     call CropHarvestPools(num_soilc, filter_soilc, dt)

     call FireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            atm2lnd_vars, energyflux_vars, soilhydrology_vars,  &
            cnstate_vars)

     call FireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars)

     call t_stopf('CNUpdate2')


     if ( use_erosion ) then
         event = 'ErosionFluxes'
         call t_start_lnd(event)
         call ErosionFluxes(bounds, num_soilc, filter_soilc, soilstate_vars, sedflux_vars )
         call t_stop_lnd(event)
     end if

     !--------------------------------------------
     ! Update3
     !--------------------------------------------

     if ( use_c13 ) then
     call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, cnstate_vars, &
          isotope=c13, isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
     end if
     if ( use_c14 ) then
        call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, cnstate_vars , &
             isotope=c14, isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
     end if

     call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_cs, veg_cs, col_cf, veg_cf, dt)

     if ( use_c13 ) then
        call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf, dt)
     end if
     if ( use_c14 ) then
        call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, dt)
     end if

     if ( use_c14 ) then
        call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars)

        call C14BombSpike(num_soilp, filter_soilp, &
               cnstate_vars)
     end if

     call t_startf('PhosphorusWeathering')
     call PhosphorusWeathering(num_soilc, filter_soilc, &
     cnstate_vars, dt)
     call t_stopf('PhosphorusWeathering')

     call t_startf('PhosphorusAdsportion')
     call PhosphorusAdsportion(num_soilc, filter_soilc, &
            cnstate_vars, dt)
     call t_stopf('PhosphorusAdsportion')

     call t_startf('PhosphorusDesoprtion')
     call PhosphorusDesoprtion(num_soilc, filter_soilc, &
            cnstate_vars, dt)
     call t_stopf('PhosphorusDesoprtion')

     call t_startf('PhosphorusOcclusion')
     call PhosphorusOcclusion(num_soilc, filter_soilc, &
             cnstate_vars, dt)
     call t_stopf('PhosphorusOcclusion')
!     write(iulog,*)'NitrogenLeaching'
!     call t_startf('NitrogenLeaching')
!     call NitrogenLeaching(bounds, num_soilc, filter_soilc, &
!            waterstate_vars, waterflux_vars, nitrogenstate_vars, nitrogenflux_vars)
!     call t_startf('NitrogenLeaching')

!     call t_startf('PhosphorusLeaching')
!     call PhosphorusLeaching(bounds, num_soilc, filter_soilc, &
!            waterstate_vars, waterflux_vars, phosphorusstate_vars, phosphorusflux_vars)
!     call t_stopf('PhosphorusLeaching')

     call t_startf('CNUpdate3Veg')
     call NitrogenStateUpdate3Veg(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)
     call t_stopf('CNUpdate3Veg')

     call t_startf('PUpdate3Veg')
     call PhosphorusStateUpdate3Veg(bounds,num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars,dt)
     call t_stopf('PUpdate3Veg')

     call veg_cf%SummaryCH4(bounds, num_soilp, filter_soilp)
  end subroutine CNEcosystemDynBeTR1
  !-----------------------------------------------------------------------
  subroutine update_plant_nutrient_buffer(bounds,col, pft, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       nitrogenflux_vars, nitrogenstate_vars, phosphorusflux_vars, phosphorusstate_vars)
    !
    ! DESCRIPTION
    ! calculate gpp downregulation factor
    use clm_time_manager         , only : get_step_size
    use ColumnType               , only : column_physical_properties
    use VegetationType           , only : vegetation_physical_properties
    use pftvarcon                , only : noveg
    implicit none
    type(bounds_type)        , intent(in)    :: bounds
    type(column_physical_properties)        , intent(in)    :: col
    type(vegetation_physical_properties)         , intent(in)    :: pft
    integer                  , intent(in)    :: num_soilc                                      ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)                                ! filter for soil columns
    integer                  , intent(in)    :: num_soilp
    integer                  , intent(in)    :: filter_soilp(:)
    type(nitrogenflux_type)  , intent(in)    :: nitrogenflux_vars    !
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars   !
    type(phosphorusflux_type), intent(in)    :: phosphorusflux_vars
    type(phosphorusstate_type),intent(inout) :: phosphorusstate_vars

    integer :: fc, c, p
    real(r8) :: dtime
    associate(&
         plant_n_buffer_patch        => nitrogenstate_vars%plant_n_buffer_patch            , & ! Inout:  [real(r8) (:)   ] gN/m2
         plant_p_buffer_patch        => phosphorusstate_vars%plant_p_buffer_patch          , & ! Inout:  [real(r8) (:)   ] gN/m2
         smin_nh4_to_plant_patch     => nitrogenflux_vars%smin_nh4_to_plant_patch          , &
         smin_no3_to_plant_patch     => nitrogenflux_vars%smin_no3_to_plant_patch          , &
         sminp_to_plant_patch        => phosphorusflux_vars%sminp_to_plant_patch             &
    )
    dtime =  get_step_size()

    do fc=1,num_soilc
       c = filter_soilc(fc)
       do p = col%pfti(c), col%pftf(c)
         if (pft%active(p).and. (pft%itype(p) .ne. noveg)) then
           plant_n_buffer_patch(p) = plant_n_buffer_patch(p) + dtime * &
               (smin_nh4_to_plant_patch(p) + smin_no3_to_plant_patch(p))

           plant_p_buffer_patch(p) = plant_p_buffer_patch(p) + dtime * &
               sminp_to_plant_patch(p)
         endif
       enddo
    enddo
    end associate
  end subroutine update_plant_nutrient_buffer

  !-----------------------------------------------------------------------

  subroutine CNFluxStateBeTR1Summary(bounds, col, pft, num_soilc, filter_soilc, &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, &
       cnstate_vars, frictionvel_vars, canopystate_vars)
    !
    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use spmdMod              , only: masterproc
    use PhosphorusDynamicsMod         , only: PhosphorusWeathering,PhosphorusAdsportion
    use PhosphorusDynamicsMod         , only: PhosphorusDesoprtion,PhosphorusOcclusion
    use NitrogenDynamicsMod       , only: NitrogenLeaching
    use NitrogenStateUpdate3BeTRMod   , only: NitrogenStateUpdate3Soil
    use PhosphorusStateUpdate3BeTRMod     , only: PhosphorusStateUpdate3Soil
    use PrecisionControlMod  , only: PrecisionControl
    use perf_mod             , only: t_startf, t_stopf
    use shr_sys_mod          , only: shr_sys_flush
    use PhosphorusDynamicsMod         , only: PhosphorusLeaching

    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    type(column_physical_properties)        , intent(in)    :: col
    type(vegetation_physical_properties)         , intent(in)    :: pft
    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                  , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    integer                  , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                  , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    logical                  , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(frictionvel_type)   , intent(in)    :: frictionvel_vars
    type(canopystate_type)   , intent(inout) :: canopystate_vars

    character(len=64) :: event
    real(r8) :: dt
    !-----------------------------------------------------------------------
    dt = dtime_mod;


    !-----------------------------------------------------------------------
    ! only do if ed is off
    if( .not. use_fates) then

       call t_startf('NitrogenLeaching')
       call NitrogenLeaching(bounds, num_soilc, filter_soilc, dt)
       call t_startf('NitrogenLeaching')

       call t_startf('PhosphorusLeaching')
       call PhosphorusLeaching(bounds, num_soilc, filter_soilc, dt)
       call t_stopf('PhosphorusLeaching')

       call t_startf('CNUpdate3')
       call NitrogenStateUpdate3Soil(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ns,  col_nf)
       call t_stopf('CNUpdate3')

       call t_startf('PUpdate3')
       call PhosphorusStateUpdate3Soil(bounds,num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars,col_ps, col_pf, is_decomp_on= .false.)
       call t_stopf('PUpdate3')

       call t_startf('CNPsum')
       call PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call col_cf_summary_for_ch4(col_cf,bounds, num_soilc, filter_soilc)

       call veg_cf_Summary(veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'bulk', col_cf)
       if ( use_c13 ) then
          call veg_cf_Summary(veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'c13', c13_col_cf)
          call col_cf_Summary(col_cf,bounds, num_soilc, filter_soilc, 'c13')
       end if
       if ( use_c14 ) then
          call veg_cf_Summary(veg_cf,bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'c14', c14_col_cf)
          call col_cf_Summary(col_cf,bounds, num_soilc, filter_soilc, 'c14')
       end if
       call veg_cs_Summary(veg_cs,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_cs)
       if ( use_c13 ) then
          call veg_cs_Summary(c13_veg_cs,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, c13_col_cs)
          call col_cs_Summary(c13_col_cs,bounds, num_soilc, filter_soilc)
       end if
       if ( use_c14 ) then
          call veg_cs_Summary(c14_veg_cs,bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, c14_col_cs)
          call col_cs_Summary(c14_col_cs,bounds, num_soilc, filter_soilc)
       end if
       call veg_nf_Summary(veg_nf, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_nf)
       call veg_ns_Summary(veg_ns, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ns)
       call veg_pf_Summary(veg_pf, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_pf)
       call veg_ps_Summary(veg_ps, bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ps)


       call col_cf_Summary(col_cf,bounds, num_soilc, filter_soilc, 'bulk')
       call col_cs_Summary(col_cs,bounds, num_soilc, filter_soilc)

       call col_nf_Summary(col_nf,bounds, num_soilc, filter_soilc)
       call col_ns_Summary(col_ns,bounds, num_soilc, filter_soilc)

       call col_pf_Summary(col_pf,bounds, num_soilc, filter_soilc)
       call col_ps_Summary(col_ps,bounds, num_soilc, filter_soilc)

       call t_stopf('CNPsum')

    end if !end of if not use_fates block

  end subroutine CNFluxStateBeTR1Summary

  !-----------------------------------------------------------------------

  subroutine CNFluxStateBeTR0Summary(bounds, col, pft, num_soilc, filter_soilc, &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, &
       cnstate_vars, frictionvel_vars, canopystate_vars)
    !
    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use spmdMod              , only: masterproc
    use PhosphorusDynamicsMod         , only: PhosphorusWeathering,PhosphorusAdsportion
    use PhosphorusDynamicsMod         , only: PhosphorusDesoprtion,PhosphorusOcclusion
    use NitrogenDynamicsMod       , only: NitrogenLeaching
    use NitrogenStateUpdate3BeTRMod   , only: NitrogenStateUpdate3
    use PhosphorusStateUpdate3BeTRMod     , only: PhosphorusStateUpdate3
    use PrecisionControlMod  , only: PrecisionControl
    use perf_mod             , only: t_startf, t_stopf
    use shr_sys_mod          , only: shr_sys_flush
    use PhosphorusDynamicsMod         , only: PhosphorusLeaching
    use PhosphorusDynamicsMod         , only: PhosphorusBiochemMin
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    type(column_physical_properties)        , intent(in)    :: col
    type(vegetation_physical_properties)         , intent(in)    :: pft
    integer                  , intent(in)    :: num_soilc         ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)   ! filter for soil columns
    integer                  , intent(in)    :: num_soilp         ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)   ! filter for soil patches
    integer                  , intent(in)    :: num_pcropp        ! number of prog. crop patches in filter
    integer                  , intent(in)    :: filter_pcropp(:)  ! filter for prognostic crop patches
    logical                  , intent(in)    :: doalb             ! true = surface albedo calculation time step
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(frictionvel_type)   , intent(in)    :: frictionvel_vars
    type(canopystate_type)   , intent(inout) :: canopystate_vars


    character(len=64) :: event
    real(r8) :: dt
    !-----------------------------------------------------------------------
    dt = dtime_mod;

    ! only do if ed is off

    call t_startf('PhosphorusWeathering')
    call PhosphorusWeathering(num_soilc, filter_soilc, &
         cnstate_vars, dt)
    call t_stopf('PhosphorusWeathering')

    call t_startf('PhosphorusAdsportion')
    call PhosphorusAdsportion(num_soilc, filter_soilc, &
         cnstate_vars, dt)
    call t_stopf('PhosphorusAdsportion')

    call t_startf('PhosphorusDesoprtion')
    call PhosphorusDesoprtion(num_soilc, filter_soilc, &
           cnstate_vars,dt)
    call t_stopf('PhosphorusDesoprtion')

    call t_startf('PhosphorusOcclusion')
    call PhosphorusOcclusion(num_soilc, filter_soilc, &
         cnstate_vars,dt)
    call t_stopf('PhosphorusOcclusion')

    if (.not. nu_com_phosphatase) then
       call t_startf('PhosphorusBiochemMin')
       call PhosphorusBiochemMin(bounds,num_soilc, filter_soilc, &
                cnstate_vars,dt)
       call t_stopf('PhosphorusBiochemMin')

    end if

    !-----------------------------------------------------------------------
    ! pflotran: when both 'pf-bgc' and 'pf-h' on, no need to call CLM-CN's N leaching module

    call NitrogenLeaching(bounds, num_soilc, filter_soilc, dt)

    call PhosphorusLeaching(bounds, num_soilc, filter_soilc, dt)
       !-----------------------------------------------------------------------

    call t_startf('CNUpdate3')

    call NitrogenStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, dt)
    call t_stopf('CNUpdate3')


    call t_startf('PUpdate3')
    call PhosphorusStateUpdate3(bounds,num_soilc, filter_soilc, num_soilp, filter_soilp, dt)
    call t_stopf('PUpdate3')


    !-----------------------------------------------------------------------
    ! only do if ed is off
    if( .not. use_fates) then
       call t_startf('CNPsum')
       call PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call col_cf%SummaryCH4(bounds, num_soilc, filter_soilc)
       call veg_cf%SummaryCH4(bounds, num_soilp, filter_soilp)

       call veg_cf%Summary(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'bulk', col_cf)
       call col_cf%Summary(bounds, num_soilc, filter_soilc, 'bulk')
       if ( use_c13 ) then
          call c13_veg_cf%Summary(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'c13', c13_col_cf)
          call c13_col_cf%Summary(bounds, num_soilc, filter_soilc, 'c13')
       end if
       if ( use_c14 ) then
          call c14_veg_cf%Summary(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'c14', c14_col_cf)
          call c14_col_cf%Summary(bounds, num_soilc, filter_soilc, 'c14')
       end if

       call veg_cs%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_cs)
       call col_cs%Summary(bounds, num_soilc, filter_soilc)
       if ( use_c13 ) then
          call c13_veg_cs%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, c13_col_cs)
          call c13_col_cs%Summary(bounds, num_soilc, filter_soilc)
       end if
       if ( use_c14 ) then
          call c14_veg_cs%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, c14_col_cs)
          call c14_col_cs%Summary(bounds, num_soilc, filter_soilc)

       end if

       call veg_nf%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_nf)
       call col_nf%Summary(bounds, num_soilc, filter_soilc)

       call veg_ns%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ns)
       call col_ns%Summary(bounds, num_soilc, filter_soilc)

       call veg_pf%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_pf)
       call col_pf%Summary(bounds, num_soilc, filter_soilc)

       call veg_ps%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ps)
       call col_ps%Summary(bounds, num_soilc, filter_soilc)
       call t_stopf('CNPsum')

    end if !end of if not use_fates block
  end subroutine CNFluxStateBeTR0Summary

end module EcosystemDynBeTRMod
