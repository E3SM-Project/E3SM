module CNEcosystemDynMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Ecosystem dynamics: phenology, vegetation
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_sys_mod         , only : shr_sys_flush
  use clm_varctl          , only : flanduse_timeseries, use_c13, use_c14, use_ed
  use decompMod           , only : bounds_type
  use perf_mod            , only : t_startf, t_stopf
  use spmdMod             , only : masterproc
  use clm_varctl          , only : use_century_decomp
  use CNStateType         , only : cnstate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNDVType            , only : dgvs_type
  use CanopyStateType     , only : canopystate_type
  use SoilStateType       , only : soilstate_type
  use TemperatureType     , only : temperature_type
  use WaterstateType      , only : waterstate_type
  use WaterfluxType       , only : waterflux_type
  use atm2lndType         , only : atm2lnd_type
  use SoilStateType       , only : soilstate_type
  use CanopyStateType     , only : canopystate_type
  use TemperatureType     , only : temperature_type 
  use PhotosynthesisType  , only : photosyns_type
  use ch4Mod              , only : ch4_type
  use EnergyFluxType      , only : energyflux_type
  use SoilHydrologyType   , only : soilhydrology_type
  use FrictionVelocityType, only : frictionvel_type
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNEcosystemDynInit   ! Ecosystem dynamics initialization
  public :: CNEcosystemDynNoLeaching       ! Ecosystem dynamics: phenology, vegetation, before doing N leaching
  public :: CNEcosystemDynLeaching       ! Ecosystem dynamics: phenology, vegetation, doing N leaching
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNEcosystemDynInit(bounds)
    !
    ! !DESCRIPTION:
    ! Initialzation of the CN Ecosystem dynamics.
    !
    ! !USES:
    use CNAllocationMod, only : CNAllocationInit
    use CNPhenologyMod , only : CNPhenologyInit
    use CNFireMod      , only : CNFireInit
    use CNC14DecayMod  , only : C14_init_BombSpike
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds      
    !-----------------------------------------------------------------------

    call CNAllocationInit (bounds)
    call CNPhenologyInit  (bounds)
    call CNFireInit       (bounds)
    
    if ( use_c14 ) then
       call C14_init_BombSpike()
    end if

  end subroutine CNEcosystemDynInit

  !-----------------------------------------------------------------------
  subroutine CNEcosystemDynNoLeaching(bounds, &
       num_soilc, filter_soilc, &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, &
       cnstate_vars, carbonflux_vars, carbonstate_vars, &
       c13_carbonflux_vars, c13_carbonstate_vars, &
       c14_carbonflux_vars, c14_carbonstate_vars, &
       nitrogenflux_vars, nitrogenstate_vars, &
       atm2lnd_vars, waterstate_vars, waterflux_vars, &
       canopystate_vars, soilstate_vars, temperature_vars, crop_vars, ch4_vars, &
       dgvs_vars, photosyns_vars, soilhydrology_vars, energyflux_vars)
    !
    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use CNNDynamicsMod         , only: CNNDeposition,CNNFixation, CNNFert, CNSoyfix
    use CNMRespMod             , only: CNMResp
    use CNDecompMod            , only: CNDecompAlloc
    use CNPhenologyMod         , only: CNPhenology
    use CNGRespMod             , only: CNGResp
    use CNCStateUpdate1Mod     , only: CStateUpdate1,CStateUpdate0
    use CNNStateUpdate1Mod     , only: NStateUpdate1
    use CNGapMortalityMod      , only: CNGapMortality
    use CNCStateUpdate2Mod     , only: CStateUpdate2, CStateUpdate2h
    use CNNStateUpdate2Mod     , only: NStateUpdate2, NStateUpdate2h
    use CNFireMod              , only: CNFireArea, CNFireFluxes
    use CNCStateUpdate3Mod     , only: CStateUpdate3
    use CNCIsoFluxMod          , only: CIsoFlux1, CIsoFlux2, CIsoFlux2h, CIsoFlux3
    use CNC14DecayMod          , only: C14Decay, C14BombSpike
    use CNWoodProductsMod      , only: CNWoodProducts
    use CNSoilLittVertTranspMod, only: CNSoilLittVertTransp
    use CNDecompCascadeBGCMod  , only: decomp_rate_constants_bgc
    use CNDecompCascadeCNMod   , only: decomp_rate_constants_cn
    use CropType               , only: crop_type
    use dynHarvestMod          , only: CNHarvest
    use clm_varpar             , only: crop_prog
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
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
    type(atm2lnd_type)       , intent(in)    :: atm2lnd_vars 
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(crop_type)          , intent(in)    :: crop_vars
    type(ch4_type)           , intent(in)    :: ch4_vars
    type(dgvs_type)          , intent(inout) :: dgvs_vars
    type(photosyns_type)     , intent(in)    :: photosyns_vars
    type(soilhydrology_type) , intent(in)    :: soilhydrology_vars
    type(energyflux_type)    , intent(in)    :: energyflux_vars
    !-----------------------------------------------------------------------

    ! Call the main CN routines

    ! only do if ed is off
    if( .not. use_ed ) then

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

       call t_stopf('CNZero')

       ! --------------------------------------------------
       ! Nitrogen Deposition, Fixation and Respiration
       ! --------------------------------------------------

       call t_startf('CNDeposition')
       call CNNDeposition(bounds, &
            atm2lnd_vars, nitrogenflux_vars)
       call t_stopf('CNDeposition')

       call t_startf('CNFixation')
       call CNNFixation( num_soilc, filter_soilc, &
            carbonflux_vars, nitrogenflux_vars)
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


       call t_startf('CNDecompAlloc')

       if (use_century_decomp) then
          call decomp_rate_constants_bgc(bounds, num_soilc, filter_soilc, &
               canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars)
       else
          call decomp_rate_constants_cn(bounds, num_soilc, filter_soilc, &
               canopystate_vars, soilstate_vars, temperature_vars, ch4_vars, carbonflux_vars)
       end if

       call CNDecompAlloc(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            photosyns_vars, canopystate_vars, soilstate_vars, temperature_vars, waterstate_vars, &
            cnstate_vars, ch4_vars, &
            carbonstate_vars, carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, &
            nitrogenstate_vars, nitrogenflux_vars, crop_vars)

       call t_stopf('CNDecompAlloc')

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
            nitrogenstate_vars, nitrogenflux_vars)
       call t_stopf('CNPhenology')

       !--------------------------------------------
       ! Growth respiration
       !--------------------------------------------

       call t_startf('CNGResp')
       call CNGResp(num_soilp, filter_soilp, &
            carbonflux_vars)
       call t_stopf('CNGResp')

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

       call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, carbonflux_vars, carbonstate_vars)

       if ( use_c13 ) then
          call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, c13_carbonflux_vars, c13_carbonstate_vars)
       end if
       if ( use_c14 ) then
          call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, c14_carbonflux_vars, c14_carbonstate_vars)
       end if

       call NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            cnstate_vars, nitrogenflux_vars, nitrogenstate_vars)
       call t_stopf('CNUpdate1')

       call t_startf('CNSoilLittVertTransp')
       call CNSoilLittVertTransp(bounds, &
            num_soilc, filter_soilc, &
            canopystate_vars, cnstate_vars,                               &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, &
            carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars,    &
            nitrogenstate_vars, nitrogenflux_vars)
       call t_stopf('CNSoilLittVertTransp')

       call t_startf('CNGapMortality')
       call CNGapMortality( num_soilc, filter_soilc, num_soilp, filter_soilp, &
            dgvs_vars, cnstate_vars, &
            carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars)
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


       if (flanduse_timeseries /= ' ') then
          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp, &
               cnstate_vars, carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars)
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

       call CNWoodProducts(num_soilc, filter_soilc, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars, &
            carbonflux_vars, c13_carbonflux_vars, c14_carbonflux_vars, nitrogenflux_vars)

       call CNFireArea(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, &
            atm2lnd_vars, temperature_vars, energyflux_vars, soilhydrology_vars, waterstate_vars, &
            cnstate_vars, carbonstate_vars)

       call CNFireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            dgvs_vars, cnstate_vars, carbonstate_vars, nitrogenstate_vars, &
            carbonflux_vars, nitrogenflux_vars)

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

    end if !end of if not use_ed block

  end subroutine CNEcosystemDynNoLeaching
  
  !-----------------------------------------------------------------------
  subroutine CNEcosystemDynLeaching(bounds, num_soilc, filter_soilc, &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb, &
       cnstate_vars, carbonflux_vars, carbonstate_vars, &
       c13_carbonflux_vars, c13_carbonstate_vars, &
       c14_carbonflux_vars, c14_carbonstate_vars, dgvs_vars, &
       nitrogenflux_vars, nitrogenstate_vars, &
       waterstate_vars, waterflux_vars, frictionvel_vars, canopystate_vars)
    !
    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use spmdMod              , only: masterproc
    use CNNDynamicsMod       , only: CNNLeaching
    use CNNStateUpdate3Mod   , only: NStateUpdate3
    use CNPrecisionControlMod, only: CNPrecisionControl
    use perf_mod             , only: t_startf, t_stopf
    use shr_sys_mod          , only: shr_sys_flush
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds  
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
    type(dgvs_type)          , intent(in)    :: dgvs_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(waterflux_type)     , intent(in)    :: waterflux_vars
    type(frictionvel_type)   , intent(in)    :: frictionvel_vars
    type(canopystate_type)   , intent(inout) :: canopystate_vars
    !-----------------------------------------------------------------------
  
    ! only do if ed is off
    if( .not. use_ed ) then

       call CNNLeaching(bounds, num_soilc, filter_soilc, &
            waterstate_vars, waterflux_vars, nitrogenstate_vars, nitrogenflux_vars)

       call t_startf('CNUpdate3')

       call NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            nitrogenflux_vars, nitrogenstate_vars)
       call t_stopf('CNUpdate3')

       call t_startf('CNsum')
       call CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp, &
            carbonstate_vars, c13_carbonstate_vars, c14_carbonstate_vars, nitrogenstate_vars)

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

       call t_stopf('CNsum')

    end if !end of if not use_ed block

  end subroutine CNEcosystemDynLeaching
  
end  module CNEcosystemDynMod
