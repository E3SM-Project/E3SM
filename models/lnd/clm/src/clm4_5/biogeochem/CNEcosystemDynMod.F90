module CNEcosystemDynMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Ecosystem dynamics: phenology, vegetation
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varctl  , only: fpftdyn, use_c13, use_c14
  use decompMod   , only: bounds_type
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
    type(bounds_type), intent(in) :: bounds      ! bounds
    !-----------------------------------------------------------------------

    call CNAllocationInit (bounds)
    call CNPhenologyInit  (bounds)
    call CNFireInit       (bounds)
    
    if ( use_c14 ) then
       call C14_init_BombSpike()
    end if

  end subroutine CNEcosystemDynInit

  !-----------------------------------------------------------------------

  subroutine CNEcosystemDynNoLeaching(bounds, num_soilc, filter_soilc, &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb)
    !
    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use clmtype
    use spmdMod                , only: masterproc
    use CNSetValueMod          , only: CNZeroFluxes
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
    use CNAnnualUpdateMod      , only: CNAnnualUpdate
    use CNCIsoFluxMod          , only: CIsoFlux1, CIsoFlux2, CIsoFlux2h, CIsoFlux3
    use CNC14DecayMod          , only: C14Decay, C14BombSpike
    use dynHarvestMod          , only: CNHarvest
    use CNWoodProductsMod      , only: CNWoodProducts
    use CNSoilLittVertTranspMod, only: CNSoilLittVertTransp
    use perf_mod               , only: t_startf, t_stopf
    use clm_varpar             , only: crop_prog
    use shr_sys_mod            , only: shr_sys_flush
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_soilc         ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:)   ! filter for soil columns
    integer, intent(in) :: num_soilp         ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:)   ! filter for soil pfts
    integer, intent(in) :: num_pcropp        ! number of prog. crop pfts in filter
    integer, intent(in) :: filter_pcropp(:)  ! filter for prognostic crop pfts
    logical, intent(in) :: doalb             ! true = surface albedo calculation time step
    !-----------------------------------------------------------------------

       ! Call the main CN routines
       call t_startf('CNZero')
       call CNZeroFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)
       call t_stopf('CNZero')

       call t_startf('CNDeposition')
       call CNNDeposition(bounds)
       call t_stopf('CNDeposition')

       call t_startf('CNFixation')
       call CNNFixation(num_soilc,filter_soilc)
       call t_stopf('CNFixation')

       call t_startf('CNMResp')
       if (crop_prog) call CNNFert(bounds, num_soilc,filter_soilc)

       if (crop_prog) call CNSoyfix(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNMResp(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
       call t_stopf('CNMResp')

       call t_startf('CNDecompAlloc')
       call CNDecompAlloc(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp )
       call t_stopf('CNDecompAlloc')

       ! CNphenology needs to be called after CNdecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call t_startf('CNPhenology')
       call CNPhenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                        num_pcropp, filter_pcropp, doalb)
       call t_stopf('CNPhenology')

       call t_startf('CNUpdate0')
       call CNGResp(num_soilp, filter_soilp)
       
       call CStateUpdate0(num_soilp, filter_soilp, 'bulk')

       if ( use_c13 ) call CStateUpdate0(num_soilp, filter_soilp, 'c13')
       call t_stopf('CNUpdate0')

       call t_startf('CNUpdate1')
       if ( use_c13 ) call CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       if ( use_c14 ) call CStateUpdate0(num_soilp, filter_soilp, 'c14')

       if ( use_c14 ) call CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')

       call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')

       if ( use_c13 ) call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')
       
       if ( use_c14 ) call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')
       
       call NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call t_startf('CNSoilLittVertTransp')
       call CNSoilLittVertTransp(bounds, num_soilc, filter_soilc)
       call t_stopf('CNSoilLittVertTransp')

       call t_startf('CNGapMortality')
       call CNGapMortality(num_soilc, filter_soilc, num_soilp, filter_soilp)
       call t_stopf('CNGapMortality')

       call t_stopf('CNUpdate1')

       call t_startf('CNUpdate2')
       if ( use_c13 ) call CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       if ( use_c14 ) call CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')

       call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')

       if ( use_c13 ) call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       if ( use_c14 ) call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')

       call NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       if (fpftdyn /= ' ') then
          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp)
       end if 

       if ( use_c13 ) call CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       if ( use_c14 ) call CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')

       call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')

       if ( use_c13 ) call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       if ( use_c14 ) call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')

       call NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       call CNWoodProducts(num_soilc, filter_soilc)
       
       call CNFireArea(bounds, num_soilc, filter_soilc,num_soilp, filter_soilp)

       call CNFireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call t_stopf('CNUpdate2')


       if ( use_c13 ) call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       if ( use_c14 ) call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')

       call CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')

       if ( use_c13 ) call CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       if ( use_c14 ) call CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')

       if ( use_c14 ) call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp)

       if ( use_c14 ) call C14BombSpike(num_soilp, filter_soilp)

  end subroutine CNEcosystemDynNoLeaching
  
  !-----------------------------------------------------------------------
  subroutine CNEcosystemDynLeaching(bounds, num_soilc, filter_soilc, &
       num_soilp, filter_soilp, num_pcropp, filter_pcropp, doalb)
    !
    ! !DESCRIPTION:
    ! The core CN code is executed here. Calculates fluxes for maintenance
    ! respiration, decomposition, allocation, phenology, and growth respiration.
    ! These routines happen on the radiation time step so that canopy structure
    ! stays synchronized with albedo calculations.
    !
    ! !USES:
    use clmtype
    use spmdMod              , only: masterproc
    use CNNDynamicsMod       , only: CNNLeaching
    use CNNStateUpdate3Mod   , only: NStateUpdate3
    use CNPrecisionControlMod, only: CNPrecisionControl
    use CNVegStructUpdateMod , only: CNVegStructUpdate
    use CNSummaryMod         , only: CSummary, NSummary
    use perf_mod             , only: t_startf, t_stopf
    use shr_sys_mod          , only: shr_sys_flush
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    integer, intent(in) :: num_soilc         ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(:)   ! filter for soil columns
    integer, intent(in) :: num_soilp         ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(:)   ! filter for soil pfts
    integer, intent(in) :: num_pcropp        ! number of prog. crop pfts in filter
    integer, intent(in) :: filter_pcropp(:)  ! filter for prognostic crop pfts
    logical, intent(in) :: doalb             ! true = surface albedo calculation time step
    !-----------------------------------------------------------------------
  
    call CNNLeaching(bounds, num_soilc, filter_soilc)

    call t_startf('CNUpdate3')
       
    call NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)
    call t_stopf('CNUpdate3')
    
    call t_startf('CNsum')
    call CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)
    
    if (doalb) then   
       call CNVegStructUpdate(num_soilp, filter_soilp)
    end if
   
    call CSummary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')
    
    if ( use_c13 ) call CSummary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')
    
    if ( use_c14 ) call CSummary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')
    
    call NSummary(bounds, num_soilc, filter_soilc, num_soilp, filter_soilp)
    call t_stopf('CNsum')
    
  end subroutine CNEcosystemDynLeaching
  
end  module CNEcosystemDynMod
