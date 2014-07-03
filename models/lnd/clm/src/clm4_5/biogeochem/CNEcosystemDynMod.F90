module CNEcosystemDynMod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNEcosystemDynMod
!
! !DESCRIPTION:
! Ecosystem dynamics: phenology, vegetation
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varctl  , only: fpftdyn, use_c13, use_c14
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNEcosystemDynInit   ! Ecosystem dynamics initialization
  public :: CNEcosystemDyn       ! Ecosystem dynamics: phenology, vegetation
!
!
! !REVISION HISTORY:
! Created by Peter Thornton
! 19 May 2009: PET - modified to include call to harvest routine
!F. Li and S. Levis (11/06/12)
!
! !PRIVATE MEMBER FUNCTIONS:
!
! !PRIVATE TYPES:
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNEcosystemDynInit
!
! !INTERFACE:
  subroutine CNEcosystemDynInit(lbg, ubg, lbc, ubc, lbp, ubp )
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
    integer, intent(in) :: lbg, ubg        ! gridcell bounds
    integer, intent(in) :: lbc, ubc        ! column bounds
    integer, intent(in) :: lbp, ubp        ! pft bounds
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 04/05/11, Erik Kluzek creation
!
! !LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------
     call CNAllocationInit ( lbc, ubc, lbp, ubp )
     call CNPhenologyInit  ( lbp, ubp )
     call CNFireInit       ( lbg, ubg )

     if ( use_c14 ) call C14_init_BombSpike()

  end subroutine CNEcosystemDynInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNEcosystemDyn
!
! !INTERFACE:
  subroutine CNEcosystemDyn(lbc, ubc, lbp, ubp, num_soilc, filter_soilc, &
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
    use CNSetValueMod        , only: CNZeroFluxes
    use CNNDynamicsMod       , only: CNNDeposition,CNNFixation, CNNLeaching, CNNFert, CNSoyfix
    use CNMRespMod           , only: CNMResp
    use CNDecompMod          , only: CNDecompAlloc
    use CNPhenologyMod       , only: CNPhenology
    use CNGRespMod           , only: CNGResp
    use CNCStateUpdate1Mod   , only: CStateUpdate1,CStateUpdate0
    use CNNStateUpdate1Mod   , only: NStateUpdate1
    use CNGapMortalityMod    , only: CNGapMortality
    use CNCStateUpdate2Mod   , only: CStateUpdate2, CStateUpdate2h
    use CNNStateUpdate2Mod   , only: NStateUpdate2, NStateUpdate2h
    use CNFireMod            , only: CNFireArea, CNFireFluxes
    use CNCStateUpdate3Mod   , only: CStateUpdate3
    use CNNStateUpdate3Mod   , only: NStateUpdate3
    use CNPrecisionControlMod, only: CNPrecisionControl
    use CNVegStructUpdateMod , only: CNVegStructUpdate
    use CNAnnualUpdateMod    , only: CNAnnualUpdate
    use CNSummaryMod         , only: CSummary, NSummary
    use CNCIsoFluxMod        , only: CIsoFlux1, CIsoFlux2, CIsoFlux2h, CIsoFlux3
    use CNC14DecayMod        , only: C14Decay, C14BombSpike
    use pftdynMod              , only: CNHarvest
    use CNWoodProductsMod      , only: CNWoodProducts
    use CNSoilLittVertTranspMod, only: CNSoilLittVertTransp
    use perf_mod               , only: t_startf, t_stopf
    use surfrdMod              , only: crop_prog
    use shr_sys_mod            , only: shr_sys_flush
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: lbc, ubc        ! column bounds
    integer, intent(in) :: lbp, ubp        ! pft bounds
    integer, intent(in) :: num_soilc       ! number of soil columns in filter
    integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
    integer, intent(in) :: num_soilp       ! number of soil pfts in filter
    integer, intent(in) :: filter_soilp(ubp-lbp+1) ! filter for soil pfts
    integer, intent(in) :: num_pcropp      ! number of prog. crop pfts in filter
    integer, intent(in) :: filter_pcropp(ubp-lbp+1)! filter for prognostic crop pfts
    logical, intent(in) :: doalb           ! true = surface albedo calculation time step
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 10/22/03, Peter Thornton: created from EcosystemDyn during migration to
!                           new vector code.
! 11/3/03, Peter Thornton: removed update of elai, esai, frac_veg_nosno_alb.
!     These are now done in CNVegStructUpdate(), which is called
!     prior to SurfaceAlbedo().
! 11/13/03, Peter Thornton: switched from nolake to soil filtering.
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arguments
!
! local pointers to implicit out arguments
!
! !OTHER LOCAL VARIABLES:
!
!EOP
!-----------------------------------------------------------------------

 !   if (doalb) then

       ! Call the main CN routines
       call t_startf('CNZero')
       call CNZeroFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)
       call t_stopf('CNZero')

       call t_startf('CNDeposition')
       call CNNDeposition(lbc, ubc)
       call t_stopf('CNDeposition')

       call t_startf('CNFixation')
       call CNNFixation(num_soilc,filter_soilc)
       call t_stopf('CNFixation')

       call t_startf('CNMResp')
       if (crop_prog) call CNNFert(num_soilc,filter_soilc)

       if (crop_prog) call CNSoyfix(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNMResp(lbc, ubc, num_soilc, filter_soilc, num_soilp, filter_soilp)
       call t_stopf('CNMResp')

       call t_startf('CNDecompAlloc')
       call CNDecompAlloc(lbp, ubp, lbc, ubc, num_soilc, filter_soilc, &
                          num_soilp, filter_soilp )
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
       call CNSoilLittVertTransp(lbc, ubc, num_soilc, filter_soilc)
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
       
       call CNFireArea(num_soilc, filter_soilc,num_soilp, filter_soilp)

       call CNFireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNNLeaching(lbc, ubc, num_soilc, filter_soilc)
       call t_stopf('CNUpdate2')

       call t_startf('CNUpdate3')
       if ( use_c13 ) call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       if ( use_c14 ) call CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')

       call CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')

       if ( use_c13 ) call CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       if ( use_c14 ) call CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')

       if ( use_c14 ) call C14Decay(num_soilc, filter_soilc, num_soilp, filter_soilp)

       if ( use_c14 ) call C14BombSpike(num_soilp, filter_soilp)

       call NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)
       call t_stopf('CNUpdate3')

       call t_startf('CNsum')
       call CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)

    if (doalb) then   
       call CNVegStructUpdate(num_soilp, filter_soilp)
    end if

!       call CNAnnualUpdate(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       call CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp, 'bulk')
       
       if ( use_c13 ) call CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c13')

       if ( use_c14 ) call CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp, 'c14')
       
       call NSummary(num_soilc, filter_soilc, num_soilp, filter_soilp)
       call t_stopf('CNsum')

!    end if  !end of if-doalb block

  end subroutine CNEcosystemDyn
#endif
!-----------------------------------------------------------------------
end  module CNEcosystemDynMod
