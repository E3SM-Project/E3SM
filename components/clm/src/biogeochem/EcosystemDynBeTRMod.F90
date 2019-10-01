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
  use clm_varctl                , only : use_c13, use_c14, use_fates, use_dynroot
  use clm_varpar                , only : nlevsoi
  use decompMod                 , only : bounds_type
  use perf_mod                  , only : t_startf, t_stopf
  use spmdMod                   , only : masterproc
  use clm_varctl                , only : use_century_decomp
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
  use ColumnDataType      , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType      , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType      , only : col_ns, col_nf
  use ColumnDataType      , only : col_ps, col_pf
  use VegetationDataType  , only : veg_cs, c13_veg_cs, c14_veg_cs
  use VegetationDataType  , only : veg_cf, c13_veg_cf, c14_veg_cf
  use VegetationDataType  , only : veg_ns, veg_nf
  use VegetationDataType  , only : veg_ps, veg_pf

  implicit none

  private
  public :: CNEcosystemDynBeTR0
  public :: CNEcosystemDynBeTR1
  public :: CNEcosystemDynBeTR2
  public :: CNFluxStateBeTR1Summary
  public :: CNFluxStateBeTR2Summary
  contains

  !-----------------------------------------------------------------------
  subroutine CNEcosystemDynBeTR0(bounds,                             &
         num_soilc, filter_soilc,                                        &
         num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb,      &
         cnstate_vars, carbonflux_vars, carbonstate_vars,                &
         c13_carbonflux_vars, c13_carbonstate_vars,                      &
         c14_carbonflux_vars, c14_carbonstate_vars,                      &
         nitrogenflux_vars, nitrogenstate_vars,                          &
         atm2lnd_vars, waterstate_vars, waterflux_vars,                  &
         canopystate_vars, soilstate_vars, temperature_vars, crop_vars,  &
         photosyns_vars, soilhydrology_vars, energyflux_vars, &
         PlantMicKinetics_vars, ch4_vars,                                &
         phosphorusflux_vars, phosphorusstate_vars, sedflux_vars)

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
    use WoodProductsMod                  , only : WoodProducts
    use CropType                         , only : crop_type
    use dynHarvestMod                    , only : CNHarvest
    use clm_varpar                       , only : crop_prog
    use NitrogenDynamicsMod              , only : NitrogenLeaching
    use CropHarvestPoolsMod              , only : CropHarvestPools
    use PlantMicKineticsMod              , only : PlantMicKinetics_type
    use PhosphorusDynamicsMod            , only : PhosphorusDeposition,PhosphorusWeathering,PhosphorusAdsportion
    use PhosphorusDynamicsMod            , only : PhosphorusLeaching, PhosphorusOcclusion,PhosphorusDesoprtion
    use VerticalProfileMod               , only : decomp_vertprofiles
     use SoilLittVertTranspMod           , only : SoilLittVertTransp
    use RootDynamicsMod                  , only : RootDynamics
    use PhenologyFLuxLimitMod            , only : phenology_flux_limiter
    use abortutils                       , only : endrun
    use shr_log_mod                      , only : errMsg => shr_log_errMsg
    use EcosystemDynMod                  , only : EcosystemDynNoLeaching1, EcosystemDynNoLeaching2
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
    type(carbonflux_type)            , intent(inout) :: carbonflux_vars
    type(carbonstate_type)           , intent(inout) :: carbonstate_vars
    type(carbonflux_type)            , intent(inout) :: c13_carbonflux_vars
    type(carbonstate_type)           , intent(inout) :: c13_carbonstate_vars
    type(carbonflux_type)            , intent(inout) :: c14_carbonflux_vars
    type(carbonstate_type)           , intent(inout) :: c14_carbonstate_vars
    type(nitrogenflux_type)          , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)         , intent(inout) :: nitrogenstate_vars
    type(atm2lnd_type)               , intent(in)    :: atm2lnd_vars
    type(waterstate_type)            , intent(inout) :: waterstate_vars
    type(waterflux_type)             , intent(in)    :: waterflux_vars
    type(canopystate_type)           , intent(in)    :: canopystate_vars
    type(soilstate_type)             , intent(inout) :: soilstate_vars
    type(temperature_type)           , intent(inout) :: temperature_vars
    type(crop_type)                  , intent(inout) :: crop_vars
    type(ch4_type)                    , intent(in)    :: ch4_vars
    type(photosyns_type)             , intent(in)    :: photosyns_vars
    type(soilhydrology_type)         , intent(in)    :: soilhydrology_vars
    type(energyflux_type)            , intent(in)    :: energyflux_vars
    type(PlantMicKinetics_type)      , intent(inout) :: PlantMicKinetics_vars
    type(phosphorusflux_type)        , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type)       , intent(inout) :: phosphorusstate_vars
    type(sedflux_type)               , intent(in)    :: sedflux_vars

    call EcosystemDynNoLeaching1(bounds,                                  &
         num_soilc, filter_soilc,                                            &
         num_soilp, filter_soilp,                                            &
         cnstate_vars, carbonflux_vars, carbonstate_vars,                 &
         c13_carbonflux_vars,                                             &
         c14_carbonflux_vars,                                             &
         nitrogenflux_vars, nitrogenstate_vars,                           &
         atm2lnd_vars, waterstate_vars, waterflux_vars,                   &
         canopystate_vars, soilstate_vars, temperature_vars, crop_vars,   &
         ch4_vars, photosyns_vars,                                        &
         phosphorusflux_vars,phosphorusstate_vars)

    call EcosystemDynNoLeaching2(bounds,                                   &
       num_soilc, filter_soilc,                                                 &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb,               &
       cnstate_vars, carbonflux_vars, carbonstate_vars,                         &
       c13_carbonflux_vars, c13_carbonstate_vars,                               &
       c14_carbonflux_vars, c14_carbonstate_vars,                               &
       nitrogenflux_vars, nitrogenstate_vars,                                   &
       atm2lnd_vars, waterstate_vars, waterflux_vars,                           &
       canopystate_vars, soilstate_vars, temperature_vars, crop_vars, ch4_vars, &
       photosyns_vars, soilhydrology_vars, energyflux_vars,          &
       phosphorusflux_vars, phosphorusstate_vars, sedflux_vars)

  end subroutine CNEcosystemDynBeTR0

  !-----------------------------------------------------------------------
  subroutine CNEcosystemDynBeTR1(bounds, col, pft,                       &
         num_soilc, filter_soilc,                                        &
         num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb,      &
         cnstate_vars, carbonflux_vars, carbonstate_vars,                &
         c13_carbonflux_vars, c13_carbonstate_vars,                      &
         c14_carbonflux_vars, c14_carbonstate_vars,                      &
         nitrogenflux_vars, nitrogenstate_vars,                          &
         atm2lnd_vars, waterstate_vars, waterflux_vars,                  &
         canopystate_vars, soilstate_vars, temperature_vars, crop_vars,  &
         photosyns_vars, soilhydrology_vars, energyflux_vars, &
         PlantMicKinetics_vars, ch4_vars,                                &
         phosphorusflux_vars, phosphorusstate_vars, ep_betr, soil_water_retention_curve)

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
    use NitrogenStateUpdate3BeTRMod       , only : NitrogenStateUpdate3
    use PhosphorusStateUpdate1BeTRMod     , only : PhosphorusStateUpdate1
    use PhosphorusStateUpdate2BeTRMod     , only : PhosphorusStateUpdate2, PhosphorusStateUpdate2h
    use PhosphorusStateUpdate3BeTRMod     , only : PhosphorusStateUpdate3
    use GapMortalityMod                   , only : GapMortality
    use FireMod                           , only : FireArea, FireFluxes
    use CarbonIsoFluxMod                  , only : CarbonIsoFlux1, CarbonIsoFlux2, CarbonIsoFlux2h, CarbonIsoFlux3
    use C14DecayMod                       , only : C14Decay, C14BombSpike
    use WoodProductsMod                   , only : WoodProducts
    use CropType                          , only : crop_type
    use dynHarvestMod                     , only : CNHarvest
    use clm_varpar                        , only : crop_prog
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
    use subgridAveMod                     , only: p2c
    use AllocationMod                     , only : calc_nfix_stress
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
    type(carbonflux_type)            , intent(inout) :: carbonflux_vars
    type(carbonstate_type)           , intent(inout) :: carbonstate_vars
    type(carbonflux_type)            , intent(inout) :: c13_carbonflux_vars
    type(carbonstate_type)           , intent(inout) :: c13_carbonstate_vars
    type(carbonflux_type)            , intent(inout) :: c14_carbonflux_vars
    type(carbonstate_type)           , intent(inout) :: c14_carbonstate_vars
    type(nitrogenflux_type)          , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)         , intent(inout) :: nitrogenstate_vars
    type(atm2lnd_type)               , intent(in)    :: atm2lnd_vars
    type(waterstate_type)            , intent(inout) :: waterstate_vars
    type(waterflux_type)             , intent(in)    :: waterflux_vars
    type(canopystate_type)           , intent(in)    :: canopystate_vars
    type(soilstate_type)             , intent(inout) :: soilstate_vars
    type(temperature_type)           , intent(inout) :: temperature_vars
    type(crop_type)                  , intent(inout) :: crop_vars
    type(ch4_type)                   , intent(in)    :: ch4_vars
    type(photosyns_type)             , intent(in)    :: photosyns_vars
    type(soilhydrology_type)         , intent(in)    :: soilhydrology_vars
    type(energyflux_type)            , intent(in)    :: energyflux_vars
    type(PlantMicKinetics_type)      , intent(inout) :: PlantMicKinetics_vars
    type(phosphorusflux_type)        , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type)       , intent(inout) :: phosphorusstate_vars
    class(betr_simulation_alm_type)  , intent(inout) :: ep_betr
    class(soil_water_retention_curve_type) , intent(in) :: soil_water_retention_curve

    call t_startf('EcosystemDynNoLeaching1')
    call EcosystemDynNoLeaching1(bounds,                                  &
         num_soilc, filter_soilc,                                         &
         num_soilp, filter_soilp,                                         &
         cnstate_vars, carbonflux_vars, carbonstate_vars,                 &
         c13_carbonflux_vars,                                             &
         c14_carbonflux_vars,                                             &
         nitrogenflux_vars, nitrogenstate_vars,                           &
         atm2lnd_vars, waterstate_vars, waterflux_vars,                   &
         canopystate_vars, soilstate_vars, temperature_vars, crop_vars,   &
         ch4_vars, photosyns_vars,                                        &
         phosphorusflux_vars,phosphorusstate_vars)
     call t_stopf('EcosystemDynNoLeaching1')
       !----------------------------------------------------------------
       !call decomposition method from betr
       !----------------------------------------------------------------
      call t_startf('betr type1 soil bgc')

      call update_PlantMicKinetics_pars(bounds, num_soilc, filter_soilc  , &
        cnstate_vars, carbonstate_vars, carbonflux_vars,  nitrogenstate_vars, &
        nitrogenflux_vars, phosphorusstate_vars, phosphorusflux_vars, &
        PlantMicKinetics_vars)

      call ep_betr%CalcSmpL(bounds, 1, nlevsoi, num_soilc, filter_soilc, &
            temperature_vars%t_soisno_col(bounds%begc:bounds%endc,1:nlevsoi), &
            soilstate_vars, waterstate_vars, soil_water_retention_curve)


      call ep_betr%SetBiophysForcing(bounds, col, pft,  &
         waterstate_vars=waterstate_vars, temperature_vars=temperature_vars,&
         atm2lnd_vars=atm2lnd_vars, soilstate_vars=soilstate_vars, carbonflux_vars=col_cf)

      !pass in parameters into betr
      call ep_betr%EnterOutLoopBGC(bounds, col, pft, &
       num_soilc, filter_soilc, &
       col_cs, col_cf, &
       c13_col_cs, c14_col_cs,  &
       col_ns,  col_ps, PlantMicKinetics_vars)

      !do bgc inside betr
      call ep_betr%OutLoopSoilBGC(bounds, col, pft)

      !gathering outputs from betr
      call ep_betr%ExitOutLoopBGC(bounds, col, pft, &
        col_cs, col_cf, &
        c13_col_cs, c13_col_cf, &
        c14_col_cs, c14_col_cf, &
        col_ns, col_nf, veg_nf, &
        col_ps, col_pf, veg_pf)

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
        cnstate_vars, carbonstate_vars, carbonflux_vars     , &
        c13_carbonflux_vars, c14_carbonflux_vars            , &
        nitrogenstate_vars, nitrogenflux_vars               , &
        phosphorusstate_vars, phosphorusflux_vars, crop_vars)
       !--------------------------------------------
       ! Phenology
       !--------------------------------------------

       ! Phenology needs to be called after SoilLittDecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call t_startf('Phenology')
       call Phenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            num_pcropp, filter_pcropp, doalb, atm2lnd_vars, &
            waterstate_vars, temperature_vars, crop_vars, canopystate_vars, soilstate_vars, &
            cnstate_vars, carbonstate_vars, carbonflux_vars, &
            nitrogenstate_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('Phenology')

       !--------------------------------------------
       ! Growth respiration
       !--------------------------------------------

       call t_startf('GrowthResp')
       call GrowthResp(num_soilp, filter_soilp, &
            carbonflux_vars)
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
               canopystate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars,  &
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

       call t_startf('phenology_flux_limiter')
       call phenology_flux_limiter(bounds, num_soilc, filter_soilc,&
           num_soilp, filter_soilp, crop_vars, cnstate_vars,  &
           veg_cf, veg_cs , c13_veg_cf, c13_veg_cs , c14_veg_cf, c14_veg_cs , &
           veg_nf, veg_ns, veg_pf, veg_ps)

       call CNLitterToColumn(num_soilc, filter_soilc, &
         cnstate_vars, carbonflux_vars, nitrogenflux_vars,phosphorusflux_vars)

       call t_stopf('phenology_flux_limiter')
       !--------------------------------------------
       ! Update1
       !--------------------------------------------

       call t_startf('CNUpdate1')

       if ( use_c13 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13',isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if

       if ( use_c14 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14',isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            crop_vars, col_cs, veg_cs, col_cf, veg_cf, ldecomp_on=.false.)

       if ( use_c13 ) then
          call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf, ldecomp_on=.false.)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf, ldecomp_on=.false.)
       end if

       call NitrogenStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, col_ns, veg_ns, col_nf, veg_nf, ldecomp_on=.false.)

       call PhosphorusStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, col_ps, veg_ps, col_pf, veg_pf, ldecomp_on=.false.)
       call t_stopf('CNUpdate1')

       call t_startf('SoilLittVertTransp')
       call SoilLittVertTransp(bounds, &
            num_soilc, filter_soilc, &
            canopystate_vars, cnstate_vars,                               &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, &
            carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,    &
            nitrogenstate_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('SoilLittVertTransp')

       call t_startf('CNGapMortality')
       call GapMortality( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars )
       call t_stopf('CNGapMortality')

       !--------------------------------------------
       ! Update2
       !--------------------------------------------

       call t_startf('CNUpdate2')

       if ( use_c13 ) then
          call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13',isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if

       if ( use_c14 ) then
          call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14',isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate2( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_cs, veg_cs, col_cf, veg_cf)

       if ( use_c13 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf)
       end if
       call NitrogenStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ns, veg_ns, col_nf, veg_nf)

       call PhosphorusStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ps, veg_ps, col_pf, veg_pf)

       if (get_do_harvest()) then
          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
               phosphorusstate_vars, phosphorusflux_vars)
       end if

       if ( use_c13 ) then
          call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13', isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14',  isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
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

       call WoodProducts(num_soilc, filter_soilc, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, &
            carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)

       call CropHarvestPools(num_soilc, filter_soilc, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, &
            phosphorusstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
            nitrogenflux_vars, phosphorusflux_vars)

       call FireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            atm2lnd_vars, temperature_vars, energyflux_vars, soilhydrology_vars, waterstate_vars, &
            cnstate_vars, carbonstate_vars)

       call FireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, carbonstate_vars, nitrogenstate_vars, &
            carbonflux_vars,nitrogenflux_vars,phosphorusstate_vars,phosphorusflux_vars)

       call t_stopf('CNUpdate2')

       !--------------------------------------------
       ! Update3
       !--------------------------------------------

       if ( use_c13 ) then
          call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13',  isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14', isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_cs, veg_cs, col_cf, veg_cf)

       if ( use_c13 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf)
       end if


       if ( use_c14 ) then
          call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, c14_carbonstate_vars)

          call C14BombSpike(num_soilp, filter_soilp, &
               cnstate_vars)
       end if

       call t_startf('PhosphorusWeathering')
       call PhosphorusWeathering(num_soilc, filter_soilc, &
       cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('PhosphorusWeathering')

       call t_startf('PhosphorusAdsportion')
       call PhosphorusAdsportion(num_soilc, filter_soilc, &
            cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('PhosphorusAdsportion')

       call t_startf('PhosphorusDesoprtion')
       call PhosphorusDesoprtion(num_soilc, filter_soilc, &
            cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('PhosphorusDesoprtion')

       call t_startf('PhosphorusOcclusion')
       call PhosphorusOcclusion(num_soilc, filter_soilc, &
             cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('PhosphorusOcclusion')

       !-----------------------------------------------------------------------
       ! in type 1 bgc, leaching will be done in betr, evenutally.
       !-----------------------------------------------------------------------
       ! in type 1 bgc, leaching will be done in betr, evenutally.
!       call t_startf('NitrogenLeaching')
!       call NitrogenLeaching(bounds, num_soilc, filter_soilc, &
!            waterstate_vars, waterflux_vars, nitrogenstate_vars, nitrogenflux_vars)
!       call t_stopf('NitrogenLeaching')

       call t_startf('PhosphorusLeaching')
       call PhosphorusLeaching(bounds, num_soilc, filter_soilc, &
            waterstate_vars, waterflux_vars, phosphorusstate_vars, phosphorusflux_vars)
       call t_stopf('PhosphorusLeaching')

       call t_startf('CNUpdate3')
       call NitrogenStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ns, veg_ns, col_nf, veg_nf)
       call t_stopf('CNUpdate3')

       call t_startf('PUpdate3')
       call PhosphorusStateUpdate3(bounds,num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars,col_ps, veg_ps, col_pf, veg_pf,ldecomp_on=.false.)
       call t_stopf('PUpdate3')

  end subroutine CNEcosystemDynBeTR1
  !-----------------------------------------------------------------------
  subroutine CNEcosystemDynBeTR2(bounds,                             &
         num_soilc, filter_soilc,                                        &
         num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb,      &
         cnstate_vars, carbonflux_vars, carbonstate_vars,                &
         c13_carbonflux_vars, c13_carbonstate_vars,                      &
         c14_carbonflux_vars, c14_carbonstate_vars,                      &
         nitrogenflux_vars, nitrogenstate_vars,                          &
         atm2lnd_vars, waterstate_vars, waterflux_vars,                  &
         canopystate_vars, soilstate_vars, temperature_vars, crop_vars,  &
         photosyns_vars, soilhydrology_vars, energyflux_vars, &
         PlantMicKinetics_vars,                                          &
         phosphorusflux_vars, phosphorusstate_vars)

    ! Description:
    ! Update vegetation related state variables and
    ! setup fluxes and parameters for plant-microbe coupling in soibgc
    !
    ! !USES:
    use NitrogenDynamicsMod         , only: NitrogenDeposition, NitrogenFert, CNSoyfix
    use NitrogenDynamicsMod         , only: NitrogenFixation_balance
    use MaintenanceRespMod             , only: MaintenanceResp
    use PhenologyMod            , only : Phenology, CNLitterToColumn
    use GrowthRespMod             , only: GrowthResp
    use CarbonStateUpdate1BeTRMod    , only : CarbonStateUpdate1,CarbonStateUpdate0
    use NitrogenStateUpdate1BeTRMod    , only : NitrogenStateUpdate1
    use GapMortalityMod         , only : GapMortality
    use CarbonStateUpdate2BeTRMod        , only : CarbonStateUpdate2, CarbonStateUpdate2h
    use NitrogenStateUpdate2BeTRMod    , only : NitrogenStateUpdate2, NitrogenStateUpdate2h
    use FireMod              , only: FireArea, FireFluxes
    use CarbonStateUpdate3BeTRMod        , only : CarbonStateUpdate3
    use CarbonIsoFluxMod             , only : CarbonIsoFlux1, CarbonIsoFlux2, CarbonIsoFlux2h, CarbonIsoFlux3
    use C14DecayMod             , only : C14Decay, C14BombSpike
    use WoodProductsMod         , only : WoodProducts
    use CropType                  , only : crop_type
    use dynHarvestMod             , only : CNHarvest
    use clm_varpar                , only : crop_prog
    use CropHarvestPoolsMod     , only : CropHarvestPools
    use PlantMicKineticsMod       , only : PlantMicKinetics_type
    use CNAllocationBetrMod       , only : SetPlantMicNPDemand, Allocation3_PlantCNPAlloc
    use PhosphorusStateUpdate3BeTRMod          , only : PhosphorusStateUpdate3
    use NitrogenStateUpdate3BeTRMod    , only : NitrogenStateUpdate3
    use PhosphorusStateUpdate1BeTRMod      , only : PhosphorusStateUpdate1
    use PhosphorusStateUpdate2BeTRMod     , only : PhosphorusStateUpdate2, PhosphorusStateUpdate2h
    use PhosphorusDynamicsMod              , only : PhosphorusDeposition,PhosphorusWeathering,PhosphorusBiochemMin_Ptaseact
    use VerticalProfileMod      , only : decomp_vertprofiles
    use RootDynamicsMod           , only: RootDynamics
    use PhenologyFLuxLimitMod     , only : phenology_flux_limiter
    use abortutils          , only : endrun
    use shr_log_mod         , only : errMsg => shr_log_errMsg
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
    type(carbonflux_type)            , intent(inout) :: carbonflux_vars
    type(carbonstate_type)           , intent(inout) :: carbonstate_vars
    type(carbonflux_type)            , intent(inout) :: c13_carbonflux_vars
    type(carbonstate_type)           , intent(inout) :: c13_carbonstate_vars
    type(carbonflux_type)            , intent(inout) :: c14_carbonflux_vars
    type(carbonstate_type)           , intent(inout) :: c14_carbonstate_vars
    type(nitrogenflux_type)          , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)         , intent(inout) :: nitrogenstate_vars
    type(atm2lnd_type)               , intent(in)    :: atm2lnd_vars
    type(waterstate_type)            , intent(in)    :: waterstate_vars
    type(waterflux_type)             , intent(in)    :: waterflux_vars
    type(canopystate_type)           , intent(in)    :: canopystate_vars
    type(soilstate_type)             , intent(inout) :: soilstate_vars
    type(temperature_type)           , intent(inout) :: temperature_vars
    type(crop_type)                  , intent(inout) :: crop_vars
    type(photosyns_type)             , intent(in)    :: photosyns_vars
    type(soilhydrology_type)         , intent(in)    :: soilhydrology_vars
    type(energyflux_type)            , intent(in)    :: energyflux_vars
    type(PlantMicKinetics_type)      , intent(inout) :: PlantMicKinetics_vars
    type(phosphorusflux_type)        , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type)       , intent(inout) :: phosphorusstate_vars

    if(.not. use_fates)then
       ! --------------------------------------------------
       ! zero the column-level C and N fluxes
       ! --------------------------------------------------

       call t_startf('CNZero')

       call carbonflux_vars%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)

       if ( use_c13 ) then
          call c13_carbonflux_vars%SetValues( &
               num_soilp, filter_soilp, 0._r8, &
               num_soilc, filter_soilc, 0._r8)
       end if

       if ( use_c14 ) then
          call c14_carbonflux_vars%SetValues( &
               num_soilp, filter_soilp, 0._r8, &
               num_soilc, filter_soilc, 0._r8)
       end if
       call nitrogenflux_vars%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)

       call phosphorusflux_vars%SetValues( &
            num_soilp, filter_soilp, 0._r8, &
            num_soilc, filter_soilc, 0._r8)

       call t_stopf('CNZero')
!       call phosphorusstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,'loc-1')
       ! --------------------------------------------------
       ! Nitrogen Deposition, Fixation and Respiration, phosphorus dynamics
       ! --------------------------------------------------

       call t_startf('CNDeposition')
       call NitrogenDeposition(bounds, &
            atm2lnd_vars, nitrogenflux_vars)
       call t_stopf('CNDeposition')

       call t_startf('MaintenanceResp')
       if (crop_prog) then
          call  NitrogenFert(bounds, num_soilc,filter_soilc, &
               nitrogenflux_vars)

       end if
       call MaintenanceResp(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            canopystate_vars, soilstate_vars, temperature_vars, photosyns_vars, &
            carbonflux_vars, carbonstate_vars, nitrogenstate_vars)
       call t_stopf('MaintenanceResp')
!       call phosphorusstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,'afresp')
       ! for P competition purpose, calculate P fluxes that will potentially increase solution P pool
       ! then competitors take up solution P
       call t_startf('PhosphorusWeathering')
       call PhosphorusWeathering(num_soilc, filter_soilc, &
               cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('PhosphorusWeathering')


       ! --------------------------------------------------
       ! Phosphorus Deposition ! X.SHI
       ! --------------------------------------------------

       call t_startf('PhosphorusDeposition')
       call PhosphorusDeposition(bounds, &
            atm2lnd_vars, phosphorusflux_vars)
       call t_stopf('PhosphorusDeposition')

       !This specifies the vertical distribution of deposition fluxes and
       !root exudates
       call decomp_vertprofiles(bounds,                      &
           num_soilc, filter_soilc, num_soilp, filter_soilp, &
           soilstate_vars, canopystate_vars, cnstate_vars)
!!--------------------------------------------------------------
!       call phosphorusstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,'bfdemand')
       call t_startf('CNAllocation - phase-1')
       call SetPlantMicNPDemand (bounds                                     , &
                num_soilc, filter_soilc, num_soilp, filter_soilp,temperature_vars, &
                photosyns_vars, crop_vars, canopystate_vars,soilstate_vars,cnstate_vars   , &
                carbonstate_vars, carbonflux_vars, c13_carbonflux_vars      , &
                c14_carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars  , &
                phosphorusstate_vars, phosphorusflux_vars, PlantMicKinetics_vars)

       call t_stopf('CNAllocation - phase-1')
!       call phosphorusstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,'afdemand')
       call t_startf('CNFixation')
       !nfixation comes after SetPlantMicNPDemand because it needs cnp ratio
       !computed first
       call NitrogenFixation_balance( num_soilc, filter_soilc, &
               cnstate_vars, carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, &
               temperature_vars, waterstate_vars, carbonstate_vars, phosphorusstate_vars)
       call t_stopf('CNFixation')

       ! nu_com_phosphatase is true
       call t_startf('PhosphorusBiochemMin')
       call PhosphorusBiochemMin_Ptaseact(bounds,num_soilc, filter_soilc, &
                  cnstate_vars,nitrogenstate_vars,phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('PhosphorusBiochemMin')
!       call phosphorusstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,'afptase')
       if (crop_prog) then
          !be careful about CNSoyfix, it is coded by using CTC-RD formulation
          !of CN interactions
          call CNSoyfix(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               waterstate_vars, crop_vars, cnstate_vars, &
               nitrogenstate_vars, nitrogenflux_vars)
      endif
      call t_startf('CNAllocation - phase-3')
      call Allocation3_PlantCNPAlloc (bounds                      , &
                num_soilc, filter_soilc, num_soilp, filter_soilp    , &
                canopystate_vars                                    , &
                cnstate_vars, carbonstate_vars, carbonflux_vars     , &
                c13_carbonflux_vars, c14_carbonflux_vars            , &
                nitrogenstate_vars, nitrogenflux_vars               , &
                phosphorusstate_vars, phosphorusflux_vars, crop_vars)
      call t_stopf('CNAllocation - phase-3')

       !--------------------------------------------
       ! Phenology
       !--------------------------------------------
!       call phosphorusstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,'loc0')
       ! CNphenology needs to be called after CNdecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call t_startf('CNPhenology')
       call Phenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            num_pcropp, filter_pcropp, doalb, atm2lnd_vars, &
            waterstate_vars, temperature_vars, crop_vars, canopystate_vars, soilstate_vars, &
            cnstate_vars, carbonstate_vars, carbonflux_vars, &
            nitrogenstate_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)!, litfall_on=(get_nstep()/=23))
       call t_stopf('CNPhenology')

      !--------------------------------------------
       ! Growth respiration
       !--------------------------------------------

       call t_startf('CNGResp')
       call GrowthResp(num_soilp, filter_soilp, &
            carbonflux_vars)
       call t_stopf('CNGResp')

       call carbonflux_vars%summary_rr(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)

       if(use_c13) then
         call c13_carbonflux_vars%summary_rr(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)
       endif

       if(use_c14) then
         call c14_carbonflux_vars%summary_rr(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)
       endif
       !--------------------------------------------
       ! Dynamic Roots
       !--------------------------------------------

       if( use_dynroot ) then
          call t_startf('CNRootDyn')

          call RootDynamics(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               canopystate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars,  &
               cnstate_vars, crop_vars, energyflux_vars, soilstate_vars)
          call t_stopf('CNRootDyn')
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
       ! Update1
       !--------------------------------------------

       call phenology_flux_limiter(bounds, num_soilc, filter_soilc,&
           num_soilp, filter_soilp, crop_vars, cnstate_vars,  &
           veg_cf, veg_cs , c13_veg_cf, c13_veg_cs , c14_veg_cf, c14_veg_cs , &
           veg_nf, veg_ns, veg_pf, veg_ps)

      ! gather all patch-level litterfall fluxes to the column for litter C and N inputs

       call CNLitterToColumn(num_soilc, filter_soilc, &
         cnstate_vars, carbonflux_vars, nitrogenflux_vars,phosphorusflux_vars)

       call t_startf('CNUpdate1')

       if ( use_c13 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13',isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if

       if ( use_c14 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14',isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            crop_vars, col_cs, veg_cs, col_cf, veg_cf)

       if ( use_c13 ) then
          call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               crop_vars, c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf)
       end if

       call NitrogenStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, col_ns, veg_ns, col_nf, veg_nf)

       call PhosphorusStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, col_ps, veg_ps, col_pf, veg_pf)

       call t_stopf('CNUpdate1')

       call t_startf('CNGapMortality')

       call GapMortality( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars )
       call t_stopf('CNGapMortality')

       !--------------------------------------------
       ! Update2
       !--------------------------------------------

       call t_startf('CNUpdate2')

       if ( use_c13 ) then
          call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13',isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if

       if ( use_c14 ) then
          call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14',isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate2( num_soilc, filter_soilc, num_soilp, filter_soilp, &
             col_cs, veg_cs, col_cf, veg_cf)

       if ( use_c13 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf)
       end if
       call NitrogenStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ns, veg_ns, col_nf, veg_nf)

       call PhosphorusStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ps, veg_ps, col_pf, veg_pf)

       if (get_do_harvest()) then
          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
               phosphorusstate_vars, phosphorusflux_vars)
       endif
       if ( use_c13 ) then
          call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13',isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14',isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
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


       call WoodProducts(num_soilc, filter_soilc, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, &
            carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)

       call CropHarvestPools(num_soilc, filter_soilc, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, &
            phosphorusstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
            nitrogenflux_vars, phosphorusflux_vars)


       call FireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            atm2lnd_vars, temperature_vars, energyflux_vars, soilhydrology_vars, waterstate_vars, &
            cnstate_vars, carbonstate_vars)

       call FireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, carbonstate_vars, nitrogenstate_vars, &
            carbonflux_vars,nitrogenflux_vars,phosphorusstate_vars,phosphorusflux_vars)

       call t_stopf('CNUpdate2')

       !--------------------------------------------
       ! Update3
       !--------------------------------------------

       if ( use_c13 ) then
          call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13', isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14', isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_cs, veg_cs, col_cf, veg_cf)

       if ( use_c13 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf)
       end if


       if ( use_c14 ) then
          call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, c14_carbonstate_vars)

          call C14BombSpike(num_soilp, filter_soilp, &
               cnstate_vars)
       end if

       call t_startf('CNUpdate3')
       call NitrogenStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            col_ns, veg_ns, col_nf, veg_nf)
       call t_stopf('CNUpdate3')

       call t_startf('PUpdate3')
       call PhosphorusStateUpdate3(bounds,num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars,col_ps, veg_ps, col_pf, veg_pf)
       call t_stopf('PUpdate3')

    endif

  end subroutine CNEcosystemDynBeTR2


  !-----------------------------------------------------------------------
  subroutine CNFluxStateBetr2Summary(bounds, col, pft, num_soilc, filter_soilc, &
       num_soilp, filter_soilp,                                      &
       carbonflux_vars, carbonstate_vars,                            &
       c13_carbonflux_vars, c13_carbonstate_vars,                    &
       c14_carbonflux_vars, c14_carbonstate_vars,                    &
       nitrogenflux_vars, nitrogenstate_vars,                        &
       phosphorusflux_vars, phosphorusstate_vars)
    !
    ! DESCRIPTION
    ! summarize all fluxes and state varaibles, prepare for mass balance analysis
    !
    use PrecisionControlMod, only: PrecisionControl
    implicit none
    type(bounds_type)        , intent(in)    :: bounds
    type(column_physical_properties)        , intent(in)    :: col
    type(vegetation_physical_properties)         , intent(in)    :: pft
    integer                  , intent(in)    :: num_soilc            ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:)      ! filter for soil columns
    integer                  , intent(in)    :: num_soilp            ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:)      ! filter for soil patches
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars      !
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars     !
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars  !
    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars !
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars  !
    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars !
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars    !
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars   !
    type(phosphorusflux_type), intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type),intent(inout) :: phosphorusstate_vars

    call t_startf('CNsumBetr')

    call PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars,phosphorusstate_vars)

    call carbonflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp,'bulk')
    call carbonstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)

    if ( use_c13 ) then
       call c13_carbonflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       call c13_carbonstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    end if

    if ( use_c14 ) then
       call c14_carbonflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')
       call c14_carbonstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    end if

    call update_plant_nutrient_buffer(bounds, col, pft, num_soilc, filter_soilc, num_soilp, filter_soilp, &
      nitrogenflux_vars, nitrogenstate_vars, phosphorusflux_vars, phosphorusstate_vars)

    call nitrogenflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    call nitrogenstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)

    call phosphorusflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    call phosphorusstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)

    call t_stopf('CNsumBetr')

  end subroutine CNFluxStateBetr2Summary

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
       cnstate_vars, carbonflux_vars, carbonstate_vars, &
       c13_carbonflux_vars, c13_carbonstate_vars, &
       c14_carbonflux_vars, c14_carbonstate_vars, &
       nitrogenflux_vars, nitrogenstate_vars, &
       waterstate_vars, waterflux_vars, frictionvel_vars, canopystate_vars,&
       phosphorusflux_vars,phosphorusstate_vars)
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
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: c13_carbonflux_vars
    type(carbonstate_type)   , intent(inout) :: c13_carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: c14_carbonflux_vars
    type(carbonstate_type)   , intent(inout) :: c14_carbonstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(frictionvel_type)   , intent(in)    :: frictionvel_vars
    type(canopystate_type)   , intent(inout) :: canopystate_vars

    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars

    !-----------------------------------------------------------------------
    ! only do if ed is off
    if( .not. use_fates) then

       call t_startf('CNPsum')
       call PrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars,phosphorusstate_vars)

       call carbonflux_vars%summary_cflux_for_ch4(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc)

       call carbonflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')
       call carbonstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       if ( use_c13 ) then
          call c13_carbonflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')
          call c13_carbonstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       end if
       if ( use_c14 ) then
          call c14_carbonflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')
          call c14_carbonstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       end if
       call nitrogenflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       call nitrogenstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)

       call phosphorusflux_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       call phosphorusstate_vars%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)

       call t_stopf('CNPsum')

    end if !end of if not use_fates block

  end subroutine CNFluxStateBeTR1Summary

end module EcosystemDynBeTRMod
