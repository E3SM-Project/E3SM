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
  use clm_varctl                , only : flanduse_timeseries, use_c13, use_c14, use_ed, use_dynroot
  use decompMod                 , only : bounds_type
  use perf_mod                  , only : t_startf, t_stopf
  use spmdMod                   , only : masterproc
  use clm_varctl                , only : use_century_decomp
  use CNStateType               , only : cnstate_type
  use CNCarbonFluxType          , only : carbonflux_type
  use CNCarbonStateType         , only : carbonstate_type
  use CNNitrogenFluxType        , only : nitrogenflux_type
  use CNNitrogenStateType       , only : nitrogenstate_type
  use CNDVType                  , only : dgvs_type
  use CanopyStateType           , only : canopystate_type
  use SoilStateType             , only : soilstate_type
  use TemperatureType           , only : temperature_type
  use WaterstateType            , only : waterstate_type
  use WaterfluxType             , only : waterflux_type
  use atm2lndType               , only : atm2lnd_type
  use CanopyStateType           , only : canopystate_type
  use TemperatureType           , only : temperature_type
  use PhotosynthesisType        , only : photosyns_type
  use ch4Mod                    , only : ch4_type
  use EnergyFluxType            , only : energyflux_type
  use SoilHydrologyType         , only : soilhydrology_type
  use FrictionVelocityType      , only : frictionvel_type
  use tracerfluxType            , only : tracerflux_type
  use tracerstatetype           , only : tracerstate_type
  use BetrTracerType            , only : betrtracer_type
  use PhosphorusFluxType        , only : phosphorusflux_type
  use PhosphorusStateType       , only : phosphorusstate_type

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
         dgvs_vars, photosyns_vars, soilhydrology_vars, energyflux_vars, &
         PlantMicKinetics_vars,                                          &
         phosphorusflux_vars, phosphorusstate_vars)

    ! Description:
    ! Update vegetation related state variables and
    ! setup fluxes and parameters for plant-microbe coupling in soibgc
    !
    ! !USES:
    use CNNDynamicsMod            , only : CNNDeposition,CNNFixation, CNNFert, CNSoyfix
    use CNMRespMod                , only : CNMResp
    use CNDecompMod               , only : CNDecompAlloc
    use CNPhenologyMod            , only : CNPhenology
    use CNGRespMod                , only : CNGResp
    use CNCStateUpdate1Mod        , only : CStateUpdate1,CStateUpdate0
    use CNNStateUpdate1Mod        , only : NStateUpdate1
    use CNGapMortalityMod         , only : CNGapMortality
    use CNCStateUpdate2Mod        , only : CStateUpdate2, CStateUpdate2h
    use CNNStateUpdate2Mod        , only : NStateUpdate2, NStateUpdate2h
    use CNFireMod                 , only : CNFireArea, CNFireFluxes
    use CNCStateUpdate3Mod        , only : CStateUpdate3
    use CNCIsoFluxMod             , only : CIsoFlux1, CIsoFlux2, CIsoFlux2h, CIsoFlux3
    use CNC14DecayMod             , only : C14Decay, C14BombSpike
    use CNWoodProductsMod         , only : CNWoodProducts
    use CNDecompCascadeBGCMod     , only : decomp_rate_constants_bgc
    use CNDecompCascadeCNMod      , only : decomp_rate_constants_cn
    use CropType                  , only : crop_type
    use dynHarvestMod             , only : CNHarvest
    use clm_varpar                , only : crop_prog
    use CNCropHarvestPoolsMod     , only : CNCropHarvestPools
    use CNAllocationMod           , only : update_plant_stoichiometry
    use PlantMicKineticsMod       , only : PlantMicKinetics_type
    use CNAllocationBetrMod       , only : SetPlantMicNPDemand, CNAllocation3_PlantCNPAlloc
    use CNNStateUpdate3Mod        , only : NStateUpdate3
    use CNNDynamicsMod            , only : CNNFixation_balance
    use PStateUpdate1Mod          , only : PStateUpdate1
    use PStateUpdate2Mod          , only : PStateUpdate2, PStateUpdate2h
    use PDynamicsMod              , only : PBiochemMin_balance,PDeposition,PWeathering
    use CNVerticalProfileMod      , only : decomp_vertprofiles
    use CNRootDynMod              , only : CNRootDyn
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
    type(dgvs_type)                  , intent(inout) :: dgvs_vars
    type(photosyns_type)             , intent(in)    :: photosyns_vars
    type(soilhydrology_type)         , intent(in)    :: soilhydrology_vars
    type(energyflux_type)            , intent(in)    :: energyflux_vars
    type(PlantMicKinetics_type)      , intent(inout) :: PlantMicKinetics_vars
    type(phosphorusflux_type)        , intent(inout) :: phosphorusflux_vars
    type(phosphorusstate_type)       , intent(inout) :: phosphorusstate_vars

    if(.not. use_ed)then
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

       ! --------------------------------------------------
       ! Nitrogen Deposition, Fixation and Respiration, phosphorus dynamics
       ! --------------------------------------------------

       call t_startf('CNDeposition')
       call CNNDeposition(bounds, &
            atm2lnd_vars, nitrogenflux_vars)
       call t_stopf('CNDeposition')

       call t_startf('CNFixation')
       call CNNFixation_balance( num_soilc, filter_soilc, &
               cnstate_vars, carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, &
               temperature_vars, waterstate_vars, carbonstate_vars, phosphorusstate_vars)
       call t_stopf('CNFixation')

       call t_startf('CNMResp')
       if (crop_prog) then
          call CNNFert(bounds, num_soilc,filter_soilc, &
               nitrogenflux_vars)

          call CNSoyfix(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               waterstate_vars, crop_vars, cnstate_vars, &
               nitrogenstate_vars, nitrogenflux_vars)
       end if
       call CNMResp(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            canopystate_vars, soilstate_vars, temperature_vars, photosyns_vars, &
            carbonflux_vars, nitrogenstate_vars)
       call t_stopf('CNMResp')

       ! for P competition purpose, calculate P fluxes that will potentially increase solution P pool
       ! then competitors take up solution P
       call t_startf('PWeathering')
       call PWeathering(num_soilc, filter_soilc, &
               cnstate_vars,phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('PWeathering')

       ! nu_com_phosphatase is true
       call t_startf('PBiochemMin')
       call PBiochemMin_balance(bounds,num_soilc, filter_soilc, &
                  cnstate_vars,nitrogenstate_vars,phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('PBiochemMin')

       ! --------------------------------------------------
       ! Phosphorus Deposition ! X.SHI
       ! --------------------------------------------------

       call t_startf('PDeposition')
       call PDeposition(bounds, &
            atm2lnd_vars, phosphorusflux_vars)
       call t_stopf('PDeposition')

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


      call t_startf('CNAllocation - phase-3')
      call CNAllocation3_PlantCNPAlloc (bounds                      , &
                num_soilc, filter_soilc, num_soilp, filter_soilp    , &
                canopystate_vars                                    , &
                cnstate_vars, carbonstate_vars, carbonflux_vars     , &
                c13_carbonflux_vars, c14_carbonflux_vars            , &
                nitrogenstate_vars, nitrogenflux_vars               , &
                phosphorusstate_vars, phosphorusflux_vars)
      call t_stopf('CNAllocation - phase-3')

       !--------------------------------------------
       ! Phenology
       !--------------------------------------------

       ! CNphenology needs to be called after CNdecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call t_startf('CNPhenology')
       call CNPhenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            num_pcropp, filter_pcropp, doalb, &
            waterstate_vars, temperature_vars, crop_vars, canopystate_vars, soilstate_vars, &
            dgvs_vars, cnstate_vars, carbonstate_vars, carbonflux_vars, &
            nitrogenstate_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)
       call t_stopf('CNPhenology')


      !--------------------------------------------
       ! Growth respiration
       !--------------------------------------------

       call t_startf('CNGResp')
       call CNGResp(num_soilp, filter_soilp, &
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

          call CNRootDyn(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               carbonstate_vars, nitrogenstate_vars, carbonflux_vars,  &
               cnstate_vars, crop_vars,  soilstate_vars)
          call t_stopf('CNRootDyn')
       end if

       !--------------------------------------------
       ! CNUpdate0
       !--------------------------------------------

       call t_startf('CNUpdate0')
       call CStateUpdate0(&
            num_soilp, filter_soilp, &
            carbonflux_vars, carbonstate_vars)

       if ( use_c13 ) then
          call CStateUpdate0(&
               num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars)
       end if

       if ( use_c14 ) then
          call CStateUpdate0(&
               num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars)
       end if
       call t_stopf('CNUpdate0')

       !--------------------------------------------
       ! Update1
       !--------------------------------------------

       call t_startf('CNUpdate1')

       if ( use_c13 ) then
          call CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13')
       end if

       if ( use_c14 ) then
          call CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14')
       end if

       call CStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, carbonflux_vars, carbonstate_vars)

       if ( use_c13 ) then
          call CStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, c13_carbonflux_vars, c13_carbonstate_vars)
       end if
       if ( use_c14 ) then
          call CStateUpdate1(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, c14_carbonflux_vars, c14_carbonstate_vars)
       end if

       call NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, nitrogenflux_vars, nitrogenstate_vars)

       call PStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, phosphorusflux_vars, phosphorusstate_vars)

       call t_stopf('CNUpdate1')

       call t_startf('CNGapMortality')
       call CNGapMortality( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            dgvs_vars, cnstate_vars, &
            carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars )
       call t_stopf('CNGapMortality')

       !--------------------------------------------
       ! Update2
       !--------------------------------------------

       call t_startf('CNUpdate2')

       if ( use_c13 ) then
          call CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13')
       end if

       if ( use_c14 ) then
          call CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14')
       end if

       call CStateUpdate2( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            carbonflux_vars, carbonstate_vars)

       if ( use_c13 ) then
          call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars)
       end if
       if ( use_c14 ) then
          call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars)
       end if
       call NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            nitrogenflux_vars, nitrogenstate_vars)

       call PStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            phosphorusflux_vars, phosphorusstate_vars)

       if (flanduse_timeseries /= ' ') then
          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
               phosphorusstate_vars, phosphorusflux_vars)
       end if

       if ( use_c13 ) then
          call CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13')
       end if
       if ( use_c14 ) then
          call CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14')
       end if

       call CStateUpdate2h( num_soilc, filter_soilc,  num_soilp, filter_soilp, &
            carbonflux_vars, carbonstate_vars)
       if ( use_c13 ) then
          call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars)
       end if
       if ( use_c14 ) then
          call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars)
       end if

       call NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            nitrogenflux_vars, nitrogenstate_vars)

       call PStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            phosphorusflux_vars, phosphorusstate_vars)

       call CNWoodProducts(num_soilc, filter_soilc, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, &
            carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)

       call CNCropHarvestPools(num_soilc, filter_soilc, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, &
            phosphorusstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
            nitrogenflux_vars, phosphorusflux_vars)

       call CNFireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            atm2lnd_vars, temperature_vars, energyflux_vars, soilhydrology_vars, waterstate_vars, &
            cnstate_vars, carbonstate_vars)

       call CNFireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            dgvs_vars, cnstate_vars, carbonstate_vars, nitrogenstate_vars, &
            carbonflux_vars,nitrogenflux_vars,phosphorusstate_vars,phosphorusflux_vars)

       call t_stopf('CNUpdate2')

       !--------------------------------------------
       ! Update3
       !--------------------------------------------

       if ( use_c13 ) then
          call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c13_carbonflux_vars, isotopestate_vars=c13_carbonstate_vars, &
               isotope='c13')
       end if
       if ( use_c14 ) then
          call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonflux_vars, carbonstate_vars, &
               isotopeflux_vars=c14_carbonflux_vars, isotopestate_vars=c14_carbonstate_vars, &
               isotope='c14')
       end if

       call CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            carbonflux_vars, carbonstate_vars)

       if ( use_c13 ) then
          call CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c13_carbonflux_vars, c13_carbonstate_vars)
       end if
       if ( use_c14 ) then
          call CStateUpdate3( num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_carbonflux_vars, c14_carbonstate_vars)
       end if


       if ( use_c14 ) then
          call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               c14_carbonstate_vars)

          call C14BombSpike(num_soilp, filter_soilp, &
               cnstate_vars)
       end if

       call update_plant_stoichiometry(num_soilp, filter_soilp, &
               carbonstate_vars, nitrogenstate_vars, phosphorusstate_vars)

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
    use ColumnType               , only : column_type
    use PatchType                , only : patch_type
    implicit none
    type(bounds_type)        , intent(in)    :: bounds
    type(column_type)        , intent(in)    :: col
    type(patch_type)         , intent(in)    :: pft
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

    call update_plant_nutrient_buffer(bounds, col, pft, num_soilc, filter_soilc, num_soilp, filter_soilp, &
      nitrogenflux_vars, nitrogenstate_vars, phosphorusflux_vars, phosphorusstate_vars)

    call t_stopf('CNsum')

  end subroutine CNFluxStateBetrSummary

  !-----------------------------------------------------------------------
  subroutine update_plant_nutrient_buffer(bounds,col, pft, num_soilc, filter_soilc, num_soilp, filter_soilp, &
       nitrogenflux_vars, nitrogenstate_vars, phosphorusflux_vars, phosphorusstate_vars)
    !
    ! DESCRIPTION
    ! calculate gpp downregulation factor
    use clm_time_manager         , only : get_step_size
    use ColumnType               , only : column_type
    use PatchType                , only : patch_type
    use pftvarcon                , only : noveg
    implicit none
    type(bounds_type)        , intent(in)    :: bounds
    type(column_type)        , intent(in)    :: col
    type(patch_type)         , intent(in)    :: pft
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

end module CNEcosystemDynBetrMod
