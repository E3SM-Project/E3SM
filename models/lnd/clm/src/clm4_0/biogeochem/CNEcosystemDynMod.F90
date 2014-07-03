module CNEcosystemDynMod

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
  use clm_varctl  , only: fpftdyn, use_c13
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNEcosystemDynInit   ! Ecosystem dynamics initialization
  public :: CNEcosystemDyn       ! Ecosystem dynamics: phenology, vegetation
!
! !REVISION HISTORY:
! Created by Peter Thornton
! 19 May 2009: PET - modified to include call to harvest routine
!
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
  subroutine CNEcosystemDynInit(lbc, ubc, lbp, ubp )
!
! !DESCRIPTION:
! Initialzation of the CN Ecosystem dynamics.
!
! !USES:
    use CNAllocationMod, only : CNAllocationInit
    use CNPhenologyMod , only : CNPhenologyInit
!
! !ARGUMENTS:
    implicit none
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
    use CNNDynamicsMod       , only: CNNDeposition,CNNFixation, CNNLeaching
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
    use CNBalanceCheckMod    , only: CBalanceCheck, NBalanceCheck
    use CNPrecisionControlMod, only: CNPrecisionControl
    use CNVegStructUpdateMod , only: CNVegStructUpdate
    use CNAnnualUpdateMod    , only: CNAnnualUpdate
    use CNSummaryMod         , only: CSummary, NSummary

    use CNC13StateUpdate1Mod , only: C13StateUpdate1,C13StateUpdate0
    use CNC13StateUpdate2Mod , only: C13StateUpdate2, C13StateUpdate2h
    use CNC13StateUpdate3Mod , only: C13StateUpdate3
    use CNC13FluxMod         , only: C13Flux1, C13Flux2, C13Flux2h, C13Flux3
    use C13SummaryMod        , only: C13Summary

    use pftdynMod               , only: CNHarvest
    use CNWoodProductsMod    , only: CNWoodProducts
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
    integer, intent(in) :: filter_pcropp(:)! filter for prognostic crop pfts
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
       call CNZeroFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNNDeposition(lbc, ubc)

       call CNNFixation(num_soilc,filter_soilc)

       call CNMResp(lbc, ubc, num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNDecompAlloc(lbp, ubp, lbc, ubc, num_soilc, filter_soilc, &
                          num_soilp, filter_soilp, num_pcropp)

       ! CNphenology needs to be called after CNdecompAlloc, because it
       ! depends on current time-step fluxes to new growth on the last
       ! litterfall timestep in deciduous systems

       call CNPhenology(num_soilc, filter_soilc, num_soilp, filter_soilp, &
                        num_pcropp, filter_pcropp, doalb)

       call CNGResp(num_soilp, filter_soilp)
       
       call CStateUpdate0(num_soilp, filter_soilp)

       if (use_c13) then
          call C13StateUpdate0(num_soilp, filter_soilp)

          call C13Flux1(num_soilc, filter_soilc, num_soilp, filter_soilp)
       endif

       call CStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)

       if (use_c13) then
          call C13StateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)
       endif
       
       call NStateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNGapMortality(num_soilc, filter_soilc, num_soilp, filter_soilp)

       if (use_c13) then
          call C13Flux2(num_soilc, filter_soilc, num_soilp, filter_soilp)
       endif

       call CStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)

       if (use_c13) then
          call C13StateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)
       endif

       call NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       if (fpftdyn /= ' ') then
          call CNHarvest(num_soilc, filter_soilc, num_soilp, filter_soilp)
       end if 

       if (use_c13) then
          call C13Flux2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
       endif

       call CStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)

       if (use_c13) then
          call C13StateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
       endif

       call NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       call CNWoodProducts(num_soilc, filter_soilc)
       
       call CNFireArea(num_soilc, filter_soilc)

       call CNFireFluxes(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNNLeaching(lbc, ubc, num_soilc, filter_soilc)

       if (use_c13) then
          call C13Flux3(num_soilc, filter_soilc, num_soilp, filter_soilp)
       endif

       call CStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)

       if (use_c13) then
          call C13StateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)
       endif

       call NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)

       call CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)

    if (doalb) then   
       call CNVegStructUpdate(num_soilp, filter_soilp)
    end if

!       call CNAnnualUpdate(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       call CSummary(num_soilc, filter_soilc, num_soilp, filter_soilp)
       
       if (use_c13) then
          call C13Summary(num_soilc, filter_soilc, num_soilp, filter_soilp)
       endif
       
       call NSummary(num_soilc, filter_soilc, num_soilp, filter_soilp)

!    end if  !end of if-doalb block

  end subroutine CNEcosystemDyn

!-----------------------------------------------------------------------
end  module CNEcosystemDynMod
