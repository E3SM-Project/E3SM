module CNEcosystemDynBetrMod

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
  use TemperatureType           , only : temperature_type
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
  use ColumnDataType            , only : col_cs, c13_col_cs, c14_col_cs
  use ColumnDataType            , only : col_cf, c13_col_cf, c14_col_cf
  use ColumnDataType            , only : col_ns, col_nf
  use ColumnDataType            , only : col_ps, col_pf
  use VegetationDataType        , only : veg_cs, c13_veg_cs, c14_veg_cs
  use VegetationDataType        , only : veg_cf, c13_veg_cf, c14_veg_cf
  use VegetationDataType        , only : veg_ns, veg_nf
  use VegetationDataType        , only : veg_ps, veg_pf

  implicit none

  private
  public :: CNEcosystemDynBeTR
  public :: CNFluxStateBetrSummary
  contains


  !-----------------------------------------------------------------------
  subroutine CNEcosystemDynBetr(bounds,                             &
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
    use NitrogenDynamicsMod            , only : NitrogenDeposition,NitrogenFixation, NitrogenFert, CNSoyfix
    use MaintenanceRespMod                , only : MaintenanceResp
    use SoilLittDecompMod               , only : SoilLittDecompAlloc
    use CNPhenologyBeTRMod            , only : CNPhenology
    use GrowthRespMod                , only : GrowthResp
    use CarbonStateUpdate1Mod        , only : CarbonStateUpdate1,CarbonStateUpdate0
    use CNNStateUpdate1BeTRMod        , only : NStateUpdate1
    use CNGapMortalityBeTRMod         , only : CNGapMortality
    use CarbonStateUpdate2Mod        , only : CarbonStateUpdate2, CarbonStateUpdate2h
    use CNNStateUpdate2BeTRMod        , only : NStateUpdate2, NStateUpdate2h
    use FireMod                 , only : FireArea, FireFluxes
    use CarbonStateUpdate3Mod        , only : CarbonStateUpdate3
    use CarbonIsoFluxMod             , only : CarbonIsoFlux1, CarbonIsoFlux2, CarbonIsoFlux2h, CarbonIsoFlux3
    use C14DecayMod             , only : C14Decay, C14BombSpike
    use WoodProductsMod         , only : WoodProducts
    use DecompCascadeBGCMod     , only : decomp_rate_constants_bgc
    use DecompCascadeCNMod      , only : decomp_rate_constants_cn
    use CropType                  , only : crop_type
    use dynHarvestMod             , only : CNHarvest
    use elm_varpar                , only : crop_prog
    use CropHarvestPoolsMod       , only : CropHarvestPools
    use PlantMicKineticsMod       , only : PlantMicKinetics_type
    use CNAllocationBetrMod       , only : SetPlantMicNPDemand, Allocation3_PlantCNPAlloc
    use CNNStateUpdate3BeTRMod        , only : NStateUpdate3
    use NitrogenDynamicsMod            , only : NitrogenFixation_balance
    use PhosphorusStateUpdate1Mod          , only : PhosphorusStateUpdate1
    use PhosphorusStateUpdate2Mod          , only : PhosphorusStateUpdate2, PhosphorusStateUpdate2h
    use PhosphorusDynamicsMod              , only : PhosphorusBiochemMin_balance,PhosphorusDeposition,PhosphorusWeathering
    use VerticalProfileMod      , only : decomp_vertprofiles
    use RootDynamicsMod              , only : RootDynamics
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

       call col_cf%SetValues(num_soilc, filter_soilc, 0._r8)
       call veg_cf%SetValues(num_soilp, filter_soilp, 0._r8)

       if ( use_c13 ) then
          call c13_col_cf%SetValues(num_soilc, filter_soilc, 0._r8)
          call c13_veg_cf%SetValues(num_soilp, filter_soilp, 0._r8)
       end if

       if ( use_c14 ) then
          call c14_col_cf%SetValues(num_soilc, filter_soilc, 0._r8)
          call c14_veg_cf%SetValues(num_soilp, filter_soilp, 0._r8)
       end if
       
       call col_nf%SetValues (num_soilc, filter_soilc, 0._r8)
       call veg_nf%SetValues (num_soilp, filter_soilp, 0._r8)

       call col_pf%SetValues (num_soilc, filter_soilc, 0._r8)
       call veg_pf%SetValues (num_soilp, filter_soilp, 0._r8)

       call t_stopf('CNZero')

       ! --------------------------------------------------
       ! Nitrogen Deposition, Fixation and Respiration, phosphorus dynamics
       ! --------------------------------------------------

       call t_startf('CNDeposition')
       call NitrogenDeposition(bounds, &
            atm2lnd_vars, nitrogenflux_vars)
       call t_stopf('CNDeposition')

       call t_startf('MaintenanceResp')
       if (crop_prog) then
          call NitrogenFert(bounds, num_soilc,filter_soilc, &
               nitrogenflux_vars)

       end if
       call MaintenanceResp(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            canopystate_vars, soilstate_vars, temperature_vars, photosyns_vars, &
            carbonflux_vars, carbonstate_vars, nitrogenstate_vars)
       call t_stopf('MaintenanceResp')

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

       call t_startf('CNAllocation - phase-1')
       call SetPlantMicNPDemand (bounds                                     , &
                num_soilc, filter_soilc, num_soilp, filter_soilp            , &
                photosyns_vars, crop_vars, canopystate_vars, cnstate_vars   , &
                carbonstate_vars, carbonflux_vars, c13_carbonflux_vars      , &
                c14_carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars  , &
                phosphorusstate_vars, phosphorusflux_vars, PlantMicKinetics_vars)

       call t_stopf('CNAllocation - phase-1')

       call t_startf('CNFixation')
       !nfixation comes after SetPlantMicNPDemand because it needs cnp ratio
       !computed first
       call NitrogenFixation_balance( num_soilc, filter_soilc, &
               cnstate_vars, carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, &
               temperature_vars, waterstate_vars, carbonstate_vars, phosphorusstate_vars)
       call t_stopf('CNFixation')

       ! nu_com_phosphatase is true
       call t_startf('PhosphorusBiochemMin')
       call PhosphorusBiochemMin_balance(bounds,num_soilc, filter_soilc, &
                  cnstate_vars,nitrogenstate_vars,phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('PhosphorusBiochemMin')

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

       ! CNphenology needs to be called after SoilLittDecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call t_startf('CNPhenology')
       call CNPhenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            num_pcropp, filter_pcropp, doalb, &
            waterstate_vars, temperature_vars, crop_vars, canopystate_vars, soilstate_vars, &
            cnstate_vars, carbonstate_vars, carbonflux_vars, &
            nitrogenstate_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('CNPhenology')


      !--------------------------------------------
       ! Growth respiration
       !--------------------------------------------

       call t_startf('GrowthResp')
       call GrowthResp(num_soilp, filter_soilp, &
            carbonflux_vars)
       call t_stopf('CNGResp')
       
       call veg_cf%SummaryRR(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, col_cf)
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
       ! C State Update 0
       !--------------------------------------------

       call t_startf('CarbonStateUpdate0')
       call CarbonStateUpdate0(num_soilp, filter_soilp, veg_cs, veg_cf)
       if ( use_c13 ) then
          call CarbonStateUpdate0(num_soilp, filter_soilp, c13_veg_cs, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate0(num_soilp, filter_soilp, c14_veg_cs, c14_veg_cf)
       end if
       call t_stopf('CarbonStateUpdate0')

       !--------------------------------------------
       ! Update1
       !--------------------------------------------

       call t_startf('CNUpdate1')

       if ( use_c13 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13', isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if

       if ( use_c14 ) then
          call CarbonIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14', isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
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

       call NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, nitrogenflux_vars, nitrogenstate_vars)

       call PhosphorusStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, phosphorusflux_vars, phosphorusstate_vars)

       call t_stopf('CNUpdate1')

       call t_startf('CNGapMortality')
       call CNGapMortality( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, &
            carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
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
               isotope='c13', isocol_cs=c13_col_cs, isoveg_cs=c13_veg_cs, isocol_cf=c13_col_cf, isoveg_cf=c13_veg_cf)
       end if

       if ( use_c14 ) then
          call CarbonIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14', isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate2( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            carbonflux_vars, carbonstate_vars, col_cs, veg_cs, col_cf, veg_cf)

       if ( use_c13 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars, c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars, c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf)
       end if
       call NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            nitrogenflux_vars, nitrogenstate_vars)

       call PhosphorusStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            phosphorusflux_vars, phosphorusstate_vars)

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
               isotope='c14', isocol_cs=c14_col_cs, isoveg_cs=c14_veg_cs, isocol_cf=c14_col_cf, isoveg_cf=c14_veg_cf)
       end if

       call CarbonStateUpdate2h( num_soilc, filter_soilc,  num_soilp, filter_soilp, &
            carbonflux_vars, carbonstate_vars, col_cs, veg_cs, col_cf, veg_cf)
       if ( use_c13 ) then
          call CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars, c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars, c14_col_cs, c13_veg_cs, c14_col_cf, c14_veg_cf)
       end if

       call NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            nitrogenflux_vars, nitrogenstate_vars)

       call PhosphorusStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            phosphorusflux_vars, phosphorusstate_vars)

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
            carbonflux_vars, carbonstate_vars, col_cs, veg_cs, col_cf, veg_cf)

       if ( use_c13 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars, c13_col_cs, c13_veg_cs, c13_col_cf, c13_veg_cf)
       end if
       if ( use_c14 ) then
          call CarbonStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars, c14_col_cs, c14_veg_cs, c14_col_cf, c14_veg_cf)
       end if


       if ( use_c14 ) then
          call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, c14_carbonstate_vars)

          call C14BombSpike(num_soilp, filter_soilp, &
               cnstate_vars)
       end if

    endif

  end subroutine CNEcosystemDynBetr


  !-----------------------------------------------------------------------
  subroutine CNFluxStateBetrSummary(bounds, col, pft, num_soilc, filter_soilc, &
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
    use ColumnType           , only : column_physical_properties
    use VegetationType       , only : vegetation_physical_properties
    use PrecisionControlMod  , only: PrecisionControl
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

    call veg_cf%Summary(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'bulk', col_cf)
    call col_cf%Summary(bounds, num_soilc, filter_soilc, 'bulk')
    if ( use_c13 ) then
       call c13_veg_cf%Summary(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'c13', c13_col_cf)
       call c13_col_cf%Summary(bounds, num_soilc, filter_soilc, 'c13')

    end if
    if ( use_c14 ) then
       call C14_veg_cf%Summary(bounds, num_soilp, filter_soilp, num_soilc, filter_soilc, 'c14', c14_col_cf)
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

    call update_plant_nutrient_buffer(bounds, col, pft, num_soilc, filter_soilc, num_soilp, filter_soilp, &
      nitrogenflux_vars, nitrogenstate_vars, phosphorusflux_vars, phosphorusstate_vars)

    call veg_nf%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_nf)
    call col_nf%Summary(bounds, num_soilc, filter_soilc)

    call veg_ns%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ns)
    call col_ns%Summary(bounds, num_soilc, filter_soilc)

    call veg_pf%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_pf)
    call col_pf%Summary(bounds, num_soilc, filter_soilc)

    call veg_ps%Summary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, col_ps)
    call col_ps%Summary(bounds, num_soilc, filter_soilc)

    call t_stopf('CNsumBetr')

  end subroutine CNFluxStateBetrSummary

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
         plant_n_buffer_patch        => veg_ns%plant_n_buffer            , & ! Inout:  [real(r8) (:)   ] gN/m2
         plant_p_buffer_patch        => veg_ps%plant_p_buffer          , & ! Inout:  [real(r8) (:)   ] gN/m2
         smin_nh4_to_plant_patch     => veg_nf%smin_nh4_to_plant          , &
         smin_no3_to_plant_patch     => veg_nf%smin_no3_to_plant          , &
         sminp_to_plant_patch        => veg_pf%sminp_to_plant             &
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

end module CNEcosystemDynBetrMod
