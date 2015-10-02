module CNC13FluxMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: C13FluxMod
!
! !DESCRIPTION:
! Module for 13-carbon flux variable update, non-mortality fluxes.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
!
! !PUBLIC MEMBER FUNCTIONS:
    public:: C13Flux1
    public:: C13Flux2
    public:: C13Flux2h
    public:: C13Flux3
    private:: CNC13LitterToColumn
    private:: CNC13GapPftToColumn
    private:: CNC13HarvestPftToColumn
    private:: C13FluxCalc
!
! !REVISION HISTORY:
! 4/21/2005: Created by Peter Thornton and Neil Suits
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13Flux1
!
! !INTERFACE:
subroutine C13Flux1(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, set the 13-carbon flux
! variables (except for gap-phase mortality and fire fluxes)
!
! !USES:
   use clmtype
   use clm_varctl, only : use_c13
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   type(column_type), pointer :: c
   integer :: fp,pi
!
!EOP
!-----------------------------------------------------------------------

   if (.not. use_c13) then
      RETURN
   end if

   ! set local pointers
   p => pft
   c => col
	
   ! pft-level non-mortality fluxes
   
   call C13FluxCalc(pc13f%leafc_xfer_to_leafc, pcf%leafc_xfer_to_leafc, &
                    pc13s%leafc_xfer, pcs%leafc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%frootc_xfer_to_frootc, pcf%frootc_xfer_to_frootc, &
                    pc13s%frootc_xfer, pcs%frootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%livestemc_xfer_to_livestemc, pcf%livestemc_xfer_to_livestemc, &
                    pc13s%livestemc_xfer, pcs%livestemc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%deadstemc_xfer_to_deadstemc, pcf%deadstemc_xfer_to_deadstemc, &
                    pc13s%deadstemc_xfer, pcs%deadstemc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%livecrootc_xfer_to_livecrootc, pcf%livecrootc_xfer_to_livecrootc, &
                    pc13s%livecrootc_xfer, pcs%livecrootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%deadcrootc_xfer_to_deadcrootc, pcf%deadcrootc_xfer_to_deadcrootc, &
                    pc13s%deadcrootc_xfer, pcs%deadcrootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%leafc_to_litter, pcf%leafc_to_litter, &
                    pc13s%leafc, pcs%leafc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%frootc_to_litter, pcf%frootc_to_litter, &
                    pc13s%frootc, pcs%frootc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%livestemc_to_deadstemc, pcf%livestemc_to_deadstemc, &
                    pc13s%livestemc, pcs%livestemc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%livecrootc_to_deadcrootc, pcf%livecrootc_to_deadcrootc, &
                    pc13s%livecrootc, pcs%livecrootc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%leaf_curmr, pcf%leaf_curmr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%froot_curmr, pcf%froot_curmr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%livestem_curmr, pcf%livestem_curmr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%livecroot_curmr, pcf%livecroot_curmr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%leaf_xsmr, pcf%leaf_xsmr, &
                    pc13s%totvegc, pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%froot_xsmr, pcf%froot_xsmr, &
                    pc13s%totvegc, pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%livestem_xsmr, pcf%livestem_xsmr, &
                    pc13s%totvegc, pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%livecroot_xsmr, pcf%livecroot_xsmr, &
                    pc13s%totvegc, pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_xsmrpool, pcf%cpool_to_xsmrpool, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_leafc, pcf%cpool_to_leafc, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_leafc_storage, pcf%cpool_to_leafc_storage, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_frootc, pcf%cpool_to_frootc, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_frootc_storage, pcf%cpool_to_frootc_storage, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_livestemc, pcf%cpool_to_livestemc, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_livestemc_storage, pcf%cpool_to_livestemc_storage, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_deadstemc, pcf%cpool_to_deadstemc, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_deadstemc_storage, pcf%cpool_to_deadstemc_storage, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_livecrootc, pcf%cpool_to_livecrootc, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_livecrootc_storage, pcf%cpool_to_livecrootc_storage, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_deadcrootc, pcf%cpool_to_deadcrootc, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_deadcrootc_storage, pcf%cpool_to_deadcrootc_storage, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_leaf_gr, pcf%cpool_leaf_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_froot_gr, pcf%cpool_froot_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_livestem_gr, pcf%cpool_livestem_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_deadstem_gr, pcf%cpool_deadstem_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_livecroot_gr, pcf%cpool_livecroot_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_deadcroot_gr, pcf%cpool_deadcroot_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_leaf_storage_gr, pcf%cpool_leaf_storage_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_froot_storage_gr, pcf%cpool_froot_storage_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_livestem_storage_gr, pcf%cpool_livestem_storage_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_deadstem_storage_gr, pcf%cpool_deadstem_storage_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_livecroot_storage_gr, pcf%cpool_livecroot_storage_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_deadcroot_storage_gr, pcf%cpool_deadcroot_storage_gr, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%cpool_to_gresp_storage, pcf%cpool_to_gresp_storage, &
                    pc13s%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%transfer_leaf_gr, pcf%transfer_leaf_gr, &
                    pc13s%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%transfer_froot_gr, pcf%transfer_froot_gr, &
                    pc13s%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%transfer_livestem_gr, pcf%transfer_livestem_gr, &
                    pc13s%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%transfer_deadstem_gr, pcf%transfer_deadstem_gr, &
                    pc13s%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%transfer_livecroot_gr, pcf%transfer_livecroot_gr, &
                    pc13s%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%transfer_deadcroot_gr, pcf%transfer_deadcroot_gr, &
                    pc13s%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%leafc_storage_to_xfer, pcf%leafc_storage_to_xfer, &
                    pc13s%leafc_storage, pcs%leafc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%frootc_storage_to_xfer, pcf%frootc_storage_to_xfer, &
                    pc13s%frootc_storage, pcs%frootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%livestemc_storage_to_xfer, pcf%livestemc_storage_to_xfer, &
                    pc13s%livestemc_storage, pcs%livestemc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%deadstemc_storage_to_xfer, pcf%deadstemc_storage_to_xfer, &
                    pc13s%deadstemc_storage, pcs%deadstemc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%livecrootc_storage_to_xfer, pcf%livecrootc_storage_to_xfer, &
                    pc13s%livecrootc_storage, pcs%livecrootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%deadcrootc_storage_to_xfer, pcf%deadcrootc_storage_to_xfer, &
                    pc13s%deadcrootc_storage, pcs%deadcrootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	call C13FluxCalc(pc13f%gresp_storage_to_xfer, pcf%gresp_storage_to_xfer, &
                    pc13s%gresp_storage, pcs%gresp_storage, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	! call routine to shift pft-level litterfall fluxes to column, for isotopes
   ! the non-isotope version of this routine is called in CNPhenologyMod.F90
   ! For later clean-up, it would be possible to generalize this function to operate on a single 
   ! pft-to-column flux.
   
   call CNC13LitterToColumn(num_soilc, filter_soilc)
   
   ! column-level non-mortality fluxes
   
   call C13FluxCalc(cc13f%cwdc_to_litr2c, ccf%cwdc_to_litr2c, &
                     cc13s%cwdc, ccs%cwdc, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%cwdc_to_litr3c, ccf%cwdc_to_litr3c, &
                     cc13s%cwdc, ccs%cwdc, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%litr1_hr, ccf%litr1_hr, &
                     cc13s%litr1c, ccs%litr1c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%litr1c_to_soil1c, ccf%litr1c_to_soil1c, &
                     cc13s%litr1c, ccs%litr1c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%litr2_hr, ccf%litr2_hr, &
                     cc13s%litr2c, ccs%litr2c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%litr2c_to_soil2c, ccf%litr2c_to_soil2c, &
                     cc13s%litr2c, ccs%litr2c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%litr3_hr, ccf%litr3_hr, &
                     cc13s%litr3c, ccs%litr3c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%litr3c_to_soil3c, ccf%litr3c_to_soil3c, &
                     cc13s%litr3c, ccs%litr3c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%soil1_hr, ccf%soil1_hr, &
                     cc13s%soil1c, ccs%soil1c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%soil1c_to_soil2c, ccf%soil1c_to_soil2c, &
                     cc13s%soil1c, ccs%soil1c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%soil2_hr, ccf%soil2_hr, &
                     cc13s%soil2c, ccs%soil2c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%soil2c_to_soil3c, ccf%soil2c_to_soil3c, &
                     cc13s%soil2c, ccs%soil2c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%soil3_hr, ccf%soil3_hr, &
                     cc13s%soil3c, ccs%soil3c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%soil3c_to_soil4c, ccf%soil3c_to_soil4c, &
                     cc13s%soil3c, ccs%soil3c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   call C13FluxCalc(cc13f%soil4_hr, ccf%soil4_hr, &
                     cc13s%soil4c, ccs%soil4c, &
                     num_soilc, filter_soilc, 1._r8, 0)
   
   
!	call C13FluxCalc(pc13f%fx, pcf%fx, &
!                    pc13s%sx, pcs%sx, &
!                    num_soilp, filter_soilp, 1._r8, 0)
                    
end subroutine C13Flux1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13Flux2
!
! !INTERFACE:
subroutine C13Flux2(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, set the 13-carbon fluxes for gap mortality
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   type(column_type), pointer :: c
   integer :: fp,pi
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => pft
   c => col
	
   ! pft-level gap mortality fluxes
   
   call C13FluxCalc(pc13f%m_leafc_to_litter, pcf%m_leafc_to_litter, &
                     pc13s%leafc, pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_leafc_storage_to_litter, pcf%m_leafc_storage_to_litter, &
                     pc13s%leafc_storage, pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_leafc_xfer_to_litter, pcf%m_leafc_xfer_to_litter, &
                     pc13s%leafc_xfer, pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_frootc_to_litter, pcf%m_frootc_to_litter, &
                     pc13s%frootc, pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_frootc_storage_to_litter, pcf%m_frootc_storage_to_litter, &
                     pc13s%frootc_storage, pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_frootc_xfer_to_litter, pcf%m_frootc_xfer_to_litter, &
                     pc13s%frootc_xfer, pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livestemc_to_litter, pcf%m_livestemc_to_litter, &
                     pc13s%livestemc, pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livestemc_storage_to_litter, pcf%m_livestemc_storage_to_litter, &
                     pc13s%livestemc_storage, pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livestemc_xfer_to_litter, pcf%m_livestemc_xfer_to_litter, &
                     pc13s%livestemc_xfer, pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadstemc_to_litter, pcf%m_deadstemc_to_litter, &
                     pc13s%deadstemc, pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadstemc_storage_to_litter, pcf%m_deadstemc_storage_to_litter, &
                     pc13s%deadstemc_storage, pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadstemc_xfer_to_litter, pcf%m_deadstemc_xfer_to_litter, &
                     pc13s%deadstemc_xfer, pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livecrootc_to_litter, pcf%m_livecrootc_to_litter, &
                     pc13s%livecrootc, pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livecrootc_storage_to_litter, pcf%m_livecrootc_storage_to_litter, &
                     pc13s%livecrootc_storage, pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livecrootc_xfer_to_litter, pcf%m_livecrootc_xfer_to_litter, &
                     pc13s%livecrootc_xfer, pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadcrootc_to_litter, pcf%m_deadcrootc_to_litter, &
                     pc13s%deadcrootc, pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadcrootc_storage_to_litter, pcf%m_deadcrootc_storage_to_litter, &
                     pc13s%deadcrootc_storage, pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadcrootc_xfer_to_litter, pcf%m_deadcrootc_xfer_to_litter, &
                     pc13s%deadcrootc_xfer, pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_gresp_storage_to_litter, pcf%m_gresp_storage_to_litter, &
                     pc13s%gresp_storage, pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_gresp_xfer_to_litter, pcf%m_gresp_xfer_to_litter, &
                     pc13s%gresp_xfer, pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
	! call routine to shift pft-level gap mortality fluxes to column, for isotopes
   ! the non-isotope version of this routine is in CNGapMortalityMod.F90.

	call CNC13GapPftToColumn(num_soilc, filter_soilc)
   
end subroutine C13Flux2
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13Flux2h
!
! !INTERFACE:
subroutine C13Flux2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! set the 13-carbon fluxes for harvest mortality
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   type(column_type), pointer :: c
   integer :: fp,pi
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => pft
   c => col
	
   ! pft-level gap mortality fluxes
   
   call C13FluxCalc(pc13f%hrv_leafc_to_litter, pcf%hrv_leafc_to_litter, &
                     pc13s%leafc, pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_leafc_storage_to_litter, pcf%hrv_leafc_storage_to_litter, &
                     pc13s%leafc_storage, pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_leafc_xfer_to_litter, pcf%hrv_leafc_xfer_to_litter, &
                     pc13s%leafc_xfer, pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_frootc_to_litter, pcf%hrv_frootc_to_litter, &
                     pc13s%frootc, pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_frootc_storage_to_litter, pcf%hrv_frootc_storage_to_litter, &
                     pc13s%frootc_storage, pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_frootc_xfer_to_litter, pcf%hrv_frootc_xfer_to_litter, &
                     pc13s%frootc_xfer, pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_livestemc_to_litter, pcf%hrv_livestemc_to_litter, &
                     pc13s%livestemc, pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_livestemc_storage_to_litter, pcf%hrv_livestemc_storage_to_litter, &
                     pc13s%livestemc_storage, pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_livestemc_xfer_to_litter, pcf%hrv_livestemc_xfer_to_litter, &
                     pc13s%livestemc_xfer, pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_deadstemc_to_prod10c, pcf%hrv_deadstemc_to_prod10c, &
                     pc13s%deadstemc, pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_deadstemc_to_prod100c, pcf%hrv_deadstemc_to_prod100c, &
                     pc13s%deadstemc, pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_deadstemc_storage_to_litter, pcf%hrv_deadstemc_storage_to_litter, &
                     pc13s%deadstemc_storage, pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_deadstemc_xfer_to_litter, pcf%hrv_deadstemc_xfer_to_litter, &
                     pc13s%deadstemc_xfer, pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_livecrootc_to_litter, pcf%hrv_livecrootc_to_litter, &
                     pc13s%livecrootc, pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_livecrootc_storage_to_litter, pcf%hrv_livecrootc_storage_to_litter, &
                     pc13s%livecrootc_storage, pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_livecrootc_xfer_to_litter, pcf%hrv_livecrootc_xfer_to_litter, &
                     pc13s%livecrootc_xfer, pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_deadcrootc_to_litter, pcf%hrv_deadcrootc_to_litter, &
                     pc13s%deadcrootc, pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_deadcrootc_storage_to_litter, pcf%hrv_deadcrootc_storage_to_litter, &
                     pc13s%deadcrootc_storage, pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_deadcrootc_xfer_to_litter, pcf%hrv_deadcrootc_xfer_to_litter, &
                     pc13s%deadcrootc_xfer, pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_gresp_storage_to_litter, pcf%hrv_gresp_storage_to_litter, &
                     pc13s%gresp_storage, pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%hrv_gresp_xfer_to_litter, pcf%hrv_gresp_xfer_to_litter, &
                     pc13s%gresp_xfer, pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
	call C13FluxCalc(pc13f%hrv_xsmrpool_to_atm, pcf%hrv_xsmrpool_to_atm, &
                    pc13s%totvegc, pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0)
                    
	! call routine to shift pft-level gap mortality fluxes to column, for isotopes
   ! the non-isotope version of this routine is in CNGapMortalityMod.F90.

	call CNC13HarvestPftToColumn(num_soilc, filter_soilc)
   
end subroutine C13Flux2h
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13Flux3
!
! !INTERFACE:
subroutine C13Flux3(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, set the 13-carbon fluxes for fire mortality
!
! !USES:
   use clmtype
   use pft2colMod, only: p2c
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   type(column_type), pointer :: c
   integer :: fp,pi
   real(r8), pointer :: ptrp(:)         ! pointer to input pft array
   real(r8), pointer :: ptrc(:)         ! pointer to output column array
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => pft
   c => col
	
   ! pft-level fire mortality fluxes
   
   call C13FluxCalc(pc13f%m_leafc_to_fire, pcf%m_leafc_to_fire, &
                     pc13s%leafc, pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_leafc_storage_to_fire, pcf%m_leafc_storage_to_fire, &
                     pc13s%leafc_storage, pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_leafc_xfer_to_fire, pcf%m_leafc_xfer_to_fire, &
                     pc13s%leafc_xfer, pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_frootc_to_fire, pcf%m_frootc_to_fire, &
                     pc13s%frootc, pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_frootc_storage_to_fire, pcf%m_frootc_storage_to_fire, &
                     pc13s%frootc_storage, pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_frootc_xfer_to_fire, pcf%m_frootc_xfer_to_fire, &
                     pc13s%frootc_xfer, pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livestemc_to_fire, pcf%m_livestemc_to_fire, &
                     pc13s%livestemc, pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livestemc_storage_to_fire, pcf%m_livestemc_storage_to_fire, &
                     pc13s%livestemc_storage, pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livestemc_xfer_to_fire, pcf%m_livestemc_xfer_to_fire, &
                     pc13s%livestemc_xfer, pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadstemc_to_fire, pcf%m_deadstemc_to_fire, &
                     pc13s%deadstemc, pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadstemc_to_litter_fire, pcf%m_deadstemc_to_litter_fire, &
                     pc13s%deadstemc, pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadstemc_storage_to_fire, pcf%m_deadstemc_storage_to_fire, &
                     pc13s%deadstemc_storage, pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadstemc_xfer_to_fire, pcf%m_deadstemc_xfer_to_fire, &
                     pc13s%deadstemc_xfer, pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livecrootc_to_fire, pcf%m_livecrootc_to_fire, &
                     pc13s%livecrootc, pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livecrootc_storage_to_fire, pcf%m_livecrootc_storage_to_fire, &
                     pc13s%livecrootc_storage, pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_livecrootc_xfer_to_fire, pcf%m_livecrootc_xfer_to_fire, &
                     pc13s%livecrootc_xfer, pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadcrootc_to_fire, pcf%m_deadcrootc_to_fire, &
                     pc13s%deadcrootc, pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadcrootc_to_litter_fire, pcf%m_deadcrootc_to_litter_fire, &
                     pc13s%deadcrootc, pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadcrootc_storage_to_fire, pcf%m_deadcrootc_storage_to_fire, &
                     pc13s%deadcrootc_storage, pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_deadcrootc_xfer_to_fire, pcf%m_deadcrootc_xfer_to_fire, &
                     pc13s%deadcrootc_xfer, pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_gresp_storage_to_fire, pcf%m_gresp_storage_to_fire, &
                     pc13s%gresp_storage, pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
   call C13FluxCalc(pc13f%m_gresp_xfer_to_fire, pcf%m_gresp_xfer_to_fire, &
                     pc13s%gresp_xfer, pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0)
   
	! use routine p2c to calculate the column-level flux of deadstem and deadcrootc to
   ! cwdc as the result of fire mortality.
   call p2c(num_soilc, filter_soilc, pc13f%m_deadstemc_to_litter_fire, cc13f%m_deadstemc_to_cwdc_fire)
   call p2c(num_soilc, filter_soilc, pc13f%m_deadcrootc_to_litter_fire, cc13f%m_deadcrootc_to_cwdc_fire)

   call C13FluxCalc(cc13f%m_litr1c_to_fire, ccf%m_litr1c_to_fire, &
                     cc13s%litr1c, ccs%litr1c, &
                     num_soilc, filter_soilc, 1._r8, 0)
                    
   call C13FluxCalc(cc13f%m_litr2c_to_fire, ccf%m_litr2c_to_fire, &
                     cc13s%litr2c, ccs%litr2c, &
                     num_soilc, filter_soilc, 1._r8, 0)
                    
   call C13FluxCalc(cc13f%m_litr3c_to_fire, ccf%m_litr3c_to_fire, &
                     cc13s%litr3c, ccs%litr3c, &
                     num_soilc, filter_soilc, 1._r8, 0)
                    
   call C13FluxCalc(cc13f%m_cwdc_to_fire, ccf%m_cwdc_to_fire, &
                     cc13s%cwdc, ccs%cwdc, &
                     num_soilc, filter_soilc, 1._r8, 0)
                    
!	call C13FluxCalc(pc13f%fx, pcf%fx, &
!                    pc13s%sx, pcs%sx, &
!                    num_soilc, filter_soilc, 1._r8, 0)
                    
end subroutine C13Flux3
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNC13LitterToColumn
!
! !INTERFACE:
subroutine CNC13LitterToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! called at the end of cn_phenology to gather all pft-level litterfall fluxes
! to the column level and assign them to the three litter pools
!
! !USES:
  use clmtype
  use clm_varpar, only : max_pft_per_col
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)          ! pft vegetation type
   real(r8), pointer :: wtcol(:)        ! weight (relative to column) for this pft (0-1)
   real(r8), pointer :: pwtgcell(:)     ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: leafc_to_litter(:)
   real(r8), pointer :: frootc_to_litter(:)
   real(r8), pointer :: lf_flab(:)      ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)      ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)      ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)      ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)      ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)      ! fine root litter lignin fraction
   integer , pointer :: npfts(:)        ! number of pfts for each column
   integer , pointer :: pfti(:)         ! beginning pft index for each column
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: leafc_to_litr1c(:)
   real(r8), pointer :: leafc_to_litr2c(:)
   real(r8), pointer :: leafc_to_litr3c(:)
   real(r8), pointer :: frootc_to_litr1c(:)
   real(r8), pointer :: frootc_to_litr2c(:)
   real(r8), pointer :: frootc_to_litr3c(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
    integer :: fc,c,pi,p
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => pft%itype
    wtcol                          => pft%wtcol
    pwtgcell                       => pft%wtgcell  
    leafc_to_litter                => pc13f%leafc_to_litter
    frootc_to_litter               => pc13f%frootc_to_litter
    npfts                          => col%npfts
    pfti                           => col%pfti
    lf_flab                        => pftcon%lf_flab
    lf_fcel                        => pftcon%lf_fcel
    lf_flig                        => pftcon%lf_flig
    fr_flab                        => pftcon%fr_flab
    fr_fcel                        => pftcon%fr_fcel
    fr_flig                        => pftcon%fr_flig

   ! assign local pointers to derived type arrays (out)
    leafc_to_litr1c                => cc13f%leafc_to_litr1c
    leafc_to_litr2c                => cc13f%leafc_to_litr2c
    leafc_to_litr3c                => cc13f%leafc_to_litr3c
    frootc_to_litr1c               => cc13f%frootc_to_litr1c
    frootc_to_litr2c               => cc13f%frootc_to_litr2c
    frootc_to_litr3c               => cc13f%frootc_to_litr3c

   do pi = 1,max_pft_per_col
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if ( pi <=  npfts(c) ) then
            p = pfti(c) + pi - 1
            if (pwtgcell(p)>0._r8) then

               ! leaf litter carbon fluxes
               leafc_to_litr1c(c) = leafc_to_litr1c(c) + leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               leafc_to_litr2c(c) = leafc_to_litr2c(c) + leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               leafc_to_litr3c(c) = leafc_to_litr3c(c) + leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root litter carbon fluxes
               frootc_to_litr1c(c) = frootc_to_litr1c(c) + frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               frootc_to_litr2c(c) = frootc_to_litr2c(c) + frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               frootc_to_litr3c(c) = frootc_to_litr3c(c) + frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNC13LitterToColumn
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNC13GapPftToColumn
!
! !INTERFACE:
subroutine CNC13GapPftToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! gather all pft-level gap mortality fluxes
! to the column level and assign them to the three litter pools (+ cwd pool)
!
! !USES:
  use clmtype
  use clm_varpar, only : max_pft_per_col, maxpatch_pft
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! soil column filter
!
! !CALLED FROM:
! subroutine CNphenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)      ! pft vegetation type
   real(r8), pointer :: wtcol(:)    ! pft weight relative to column (0-1)
   real(r8), pointer :: pwtgcell(:) ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: lf_flab(:)  ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)  ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)  ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)  ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)  ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)  ! fine root litter lignin fraction
   integer , pointer :: npfts(:)    ! number of pfts for each column
   integer , pointer :: pfti(:)     ! beginning pft index for each column
   real(r8), pointer :: m_leafc_to_litter(:)
   real(r8), pointer :: m_frootc_to_litter(:)
   real(r8), pointer :: m_livestemc_to_litter(:)
   real(r8), pointer :: m_deadstemc_to_litter(:)
   real(r8), pointer :: m_livecrootc_to_litter(:)
   real(r8), pointer :: m_deadcrootc_to_litter(:)
   real(r8), pointer :: m_leafc_storage_to_litter(:)
   real(r8), pointer :: m_frootc_storage_to_litter(:)
   real(r8), pointer :: m_livestemc_storage_to_litter(:)
   real(r8), pointer :: m_deadstemc_storage_to_litter(:)
   real(r8), pointer :: m_livecrootc_storage_to_litter(:)
   real(r8), pointer :: m_deadcrootc_storage_to_litter(:)
   real(r8), pointer :: m_gresp_storage_to_litter(:)
   real(r8), pointer :: m_leafc_xfer_to_litter(:)
   real(r8), pointer :: m_frootc_xfer_to_litter(:)
   real(r8), pointer :: m_livestemc_xfer_to_litter(:)
   real(r8), pointer :: m_deadstemc_xfer_to_litter(:)
   real(r8), pointer :: m_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: m_gresp_xfer_to_litter(:)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: m_leafc_to_litr1c(:)
   real(r8), pointer :: m_leafc_to_litr2c(:)
   real(r8), pointer :: m_leafc_to_litr3c(:)
   real(r8), pointer :: m_frootc_to_litr1c(:)
   real(r8), pointer :: m_frootc_to_litr2c(:)
   real(r8), pointer :: m_frootc_to_litr3c(:)
   real(r8), pointer :: m_livestemc_to_cwdc(:)
   real(r8), pointer :: m_deadstemc_to_cwdc(:)
   real(r8), pointer :: m_livecrootc_to_cwdc(:)
   real(r8), pointer :: m_deadcrootc_to_cwdc(:)
   real(r8), pointer :: m_leafc_storage_to_litr1c(:)
   real(r8), pointer :: m_frootc_storage_to_litr1c(:)
   real(r8), pointer :: m_livestemc_storage_to_litr1c(:)
   real(r8), pointer :: m_deadstemc_storage_to_litr1c(:)
   real(r8), pointer :: m_livecrootc_storage_to_litr1c(:)
   real(r8), pointer :: m_deadcrootc_storage_to_litr1c(:)
   real(r8), pointer :: m_gresp_storage_to_litr1c(:)
   real(r8), pointer :: m_leafc_xfer_to_litr1c(:)
   real(r8), pointer :: m_frootc_xfer_to_litr1c(:)
   real(r8), pointer :: m_livestemc_xfer_to_litr1c(:)
   real(r8), pointer :: m_deadstemc_xfer_to_litr1c(:)
   real(r8), pointer :: m_livecrootc_xfer_to_litr1c(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_litr1c(:)
   real(r8), pointer :: m_gresp_xfer_to_litr1c(:)
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: fc,c,pi,p               ! indices
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers
   lf_flab                        => pftcon%lf_flab
   lf_fcel                        => pftcon%lf_fcel
   lf_flig                        => pftcon%lf_flig
   fr_flab                        => pftcon%fr_flab
   fr_fcel                        => pftcon%fr_fcel
   fr_flig                        => pftcon%fr_flig

   ! assign local pointers to column-level arrays
   npfts                          => col%npfts
   pfti                           => col%pfti
   m_leafc_to_litr1c              => cc13f%m_leafc_to_litr1c
   m_leafc_to_litr2c              => cc13f%m_leafc_to_litr2c
   m_leafc_to_litr3c              => cc13f%m_leafc_to_litr3c
   m_frootc_to_litr1c             => cc13f%m_frootc_to_litr1c
   m_frootc_to_litr2c             => cc13f%m_frootc_to_litr2c
   m_frootc_to_litr3c             => cc13f%m_frootc_to_litr3c
   m_livestemc_to_cwdc            => cc13f%m_livestemc_to_cwdc
   m_deadstemc_to_cwdc            => cc13f%m_deadstemc_to_cwdc
   m_livecrootc_to_cwdc           => cc13f%m_livecrootc_to_cwdc
   m_deadcrootc_to_cwdc           => cc13f%m_deadcrootc_to_cwdc
   m_leafc_storage_to_litr1c      => cc13f%m_leafc_storage_to_litr1c
   m_frootc_storage_to_litr1c     => cc13f%m_frootc_storage_to_litr1c
   m_livestemc_storage_to_litr1c  => cc13f%m_livestemc_storage_to_litr1c
   m_deadstemc_storage_to_litr1c  => cc13f%m_deadstemc_storage_to_litr1c
   m_livecrootc_storage_to_litr1c => cc13f%m_livecrootc_storage_to_litr1c
   m_deadcrootc_storage_to_litr1c => cc13f%m_deadcrootc_storage_to_litr1c
   m_gresp_storage_to_litr1c      => cc13f%m_gresp_storage_to_litr1c
   m_leafc_xfer_to_litr1c         => cc13f%m_leafc_xfer_to_litr1c
   m_frootc_xfer_to_litr1c        => cc13f%m_frootc_xfer_to_litr1c
   m_livestemc_xfer_to_litr1c     => cc13f%m_livestemc_xfer_to_litr1c
   m_deadstemc_xfer_to_litr1c     => cc13f%m_deadstemc_xfer_to_litr1c
   m_livecrootc_xfer_to_litr1c    => cc13f%m_livecrootc_xfer_to_litr1c
   m_deadcrootc_xfer_to_litr1c    => cc13f%m_deadcrootc_xfer_to_litr1c
   m_gresp_xfer_to_litr1c         => cc13f%m_gresp_xfer_to_litr1c

   ! assign local pointers to pft-level arrays
   ivt                            => pft%itype
   wtcol                          => pft%wtcol
   pwtgcell                       => pft%wtgcell  
   m_leafc_to_litter              => pc13f%m_leafc_to_litter
   m_frootc_to_litter             => pc13f%m_frootc_to_litter
   m_livestemc_to_litter          => pc13f%m_livestemc_to_litter
   m_deadstemc_to_litter          => pc13f%m_deadstemc_to_litter
   m_livecrootc_to_litter         => pc13f%m_livecrootc_to_litter
   m_deadcrootc_to_litter         => pc13f%m_deadcrootc_to_litter
   m_leafc_storage_to_litter      => pc13f%m_leafc_storage_to_litter
   m_frootc_storage_to_litter     => pc13f%m_frootc_storage_to_litter
   m_livestemc_storage_to_litter  => pc13f%m_livestemc_storage_to_litter
   m_deadstemc_storage_to_litter  => pc13f%m_deadstemc_storage_to_litter
   m_livecrootc_storage_to_litter => pc13f%m_livecrootc_storage_to_litter
   m_deadcrootc_storage_to_litter => pc13f%m_deadcrootc_storage_to_litter
   m_gresp_storage_to_litter      => pc13f%m_gresp_storage_to_litter
   m_leafc_xfer_to_litter         => pc13f%m_leafc_xfer_to_litter
   m_frootc_xfer_to_litter        => pc13f%m_frootc_xfer_to_litter
   m_livestemc_xfer_to_litter     => pc13f%m_livestemc_xfer_to_litter
   m_deadstemc_xfer_to_litter     => pc13f%m_deadstemc_xfer_to_litter
   m_livecrootc_xfer_to_litter    => pc13f%m_livecrootc_xfer_to_litter
   m_deadcrootc_xfer_to_litter    => pc13f%m_deadcrootc_xfer_to_litter
   m_gresp_xfer_to_litter         => pc13f%m_gresp_xfer_to_litter

   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1

            if (pwtgcell(p)>0._r8) then

               ! leaf gap mortality carbon fluxes
               m_leafc_to_litr1c(c) = m_leafc_to_litr1c(c) + &
                  m_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               m_leafc_to_litr2c(c) = m_leafc_to_litr2c(c) + &
                  m_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               m_leafc_to_litr3c(c) = m_leafc_to_litr3c(c) + &
                  m_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root gap mortality carbon fluxes
               m_frootc_to_litr1c(c) = m_frootc_to_litr1c(c) + &
                  m_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               m_frootc_to_litr2c(c) = m_frootc_to_litr2c(c) + &
                  m_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               m_frootc_to_litr3c(c) = m_frootc_to_litr3c(c) + &
                  m_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood gap mortality carbon fluxes
               m_livestemc_to_cwdc(c)  = m_livestemc_to_cwdc(c)  + &
                  m_livestemc_to_litter(p)  * wtcol(p)
               m_deadstemc_to_cwdc(c)  = m_deadstemc_to_cwdc(c)  + &
                  m_deadstemc_to_litter(p)  * wtcol(p)
               m_livecrootc_to_cwdc(c) = m_livecrootc_to_cwdc(c) + &
                  m_livecrootc_to_litter(p) * wtcol(p)
               m_deadcrootc_to_cwdc(c) = m_deadcrootc_to_cwdc(c) + &
                  m_deadcrootc_to_litter(p) * wtcol(p)

               ! storage gap mortality carbon fluxes
               m_leafc_storage_to_litr1c(c)      = m_leafc_storage_to_litr1c(c)      + &
                  m_leafc_storage_to_litter(p)      * wtcol(p)
               m_frootc_storage_to_litr1c(c)     = m_frootc_storage_to_litr1c(c)     + &
                  m_frootc_storage_to_litter(p)     * wtcol(p)
               m_livestemc_storage_to_litr1c(c)  = m_livestemc_storage_to_litr1c(c)  + &
                  m_livestemc_storage_to_litter(p)  * wtcol(p)
               m_deadstemc_storage_to_litr1c(c)  = m_deadstemc_storage_to_litr1c(c)  + &
                  m_deadstemc_storage_to_litter(p)  * wtcol(p)
               m_livecrootc_storage_to_litr1c(c) = m_livecrootc_storage_to_litr1c(c) + &
                  m_livecrootc_storage_to_litter(p) * wtcol(p)
               m_deadcrootc_storage_to_litr1c(c) = m_deadcrootc_storage_to_litr1c(c) + &
                  m_deadcrootc_storage_to_litter(p) * wtcol(p)
               m_gresp_storage_to_litr1c(c)      = m_gresp_storage_to_litr1c(c)      + &
                  m_gresp_storage_to_litter(p)      * wtcol(p)

               ! transfer gap mortality carbon fluxes
               m_leafc_xfer_to_litr1c(c)      = m_leafc_xfer_to_litr1c(c)      + &
                  m_leafc_xfer_to_litter(p)      * wtcol(p)
               m_frootc_xfer_to_litr1c(c)     = m_frootc_xfer_to_litr1c(c)     + &
                  m_frootc_xfer_to_litter(p)     * wtcol(p)
               m_livestemc_xfer_to_litr1c(c)  = m_livestemc_xfer_to_litr1c(c)  + &
                  m_livestemc_xfer_to_litter(p)  * wtcol(p)
               m_deadstemc_xfer_to_litr1c(c)  = m_deadstemc_xfer_to_litr1c(c)  + &
                  m_deadstemc_xfer_to_litter(p)  * wtcol(p)
               m_livecrootc_xfer_to_litr1c(c) = m_livecrootc_xfer_to_litr1c(c) + &
                  m_livecrootc_xfer_to_litter(p) * wtcol(p)
               m_deadcrootc_xfer_to_litr1c(c) = m_deadcrootc_xfer_to_litr1c(c) + &
                  m_deadcrootc_xfer_to_litter(p) * wtcol(p)
               m_gresp_xfer_to_litr1c(c)      = m_gresp_xfer_to_litr1c(c)      + &
                  m_gresp_xfer_to_litter(p)      * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNC13GapPftToColumn
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNC13HarvestPftToColumn
!
! !INTERFACE:
subroutine CNC13HarvestPftToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! gather all pft-level harvest mortality fluxes
! to the column level and assign them to the litter, cwd, and wood product pools
!
! !USES:
  use clmtype
  use clm_varpar, only : max_pft_per_col, maxpatch_pft
!
! !ARGUMENTS:
  implicit none
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! soil column filter
!
! !CALLED FROM:
! subroutine CNphenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)      ! pft vegetation type
   real(r8), pointer :: wtcol(:)    ! pft weight relative to column (0-1)
   real(r8), pointer :: pwtgcell(:) ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: lf_flab(:)  ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)  ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)  ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)  ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)  ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)  ! fine root litter lignin fraction
   integer , pointer :: npfts(:)    ! number of pfts for each column
   integer , pointer :: pfti(:)     ! beginning pft index for each column
   real(r8), pointer :: hrv_leafc_to_litter(:)
   real(r8), pointer :: hrv_frootc_to_litter(:)
   real(r8), pointer :: hrv_livestemc_to_litter(:)
   real(r8), pointer :: phrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: phrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: hrv_livecrootc_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_to_litter(:)
   real(r8), pointer :: hrv_leafc_storage_to_litter(:)
   real(r8), pointer :: hrv_frootc_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litter(:)
   real(r8), pointer :: hrv_gresp_storage_to_litter(:)
   real(r8), pointer :: hrv_leafc_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litter(:)
   real(r8), pointer :: hrv_gresp_xfer_to_litter(:)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: hrv_leafc_to_litr1c(:)
   real(r8), pointer :: hrv_leafc_to_litr2c(:)
   real(r8), pointer :: hrv_leafc_to_litr3c(:)
   real(r8), pointer :: hrv_frootc_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_to_litr2c(:)
   real(r8), pointer :: hrv_frootc_to_litr3c(:)
   real(r8), pointer :: hrv_livestemc_to_cwdc(:)
   real(r8), pointer :: chrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: chrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: hrv_livecrootc_to_cwdc(:)
   real(r8), pointer :: hrv_deadcrootc_to_cwdc(:)
   real(r8), pointer :: hrv_leafc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_livestemc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_deadstemc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_livecrootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_deadcrootc_storage_to_litr1c(:)
   real(r8), pointer :: hrv_gresp_storage_to_litr1c(:)
   real(r8), pointer :: hrv_leafc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_frootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_livestemc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_deadstemc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_livecrootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_deadcrootc_xfer_to_litr1c(:)
   real(r8), pointer :: hrv_gresp_xfer_to_litr1c(:)
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: fc,c,pi,p               ! indices
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers
   lf_flab                        => pftcon%lf_flab
   lf_fcel                        => pftcon%lf_fcel
   lf_flig                        => pftcon%lf_flig
   fr_flab                        => pftcon%fr_flab
   fr_fcel                        => pftcon%fr_fcel
   fr_flig                        => pftcon%fr_flig

   ! assign local pointers to column-level arrays
   npfts                          => col%npfts
   pfti                           => col%pfti
   hrv_leafc_to_litr1c              => cc13f%hrv_leafc_to_litr1c
   hrv_leafc_to_litr2c              => cc13f%hrv_leafc_to_litr2c
   hrv_leafc_to_litr3c              => cc13f%hrv_leafc_to_litr3c
   hrv_frootc_to_litr1c             => cc13f%hrv_frootc_to_litr1c
   hrv_frootc_to_litr2c             => cc13f%hrv_frootc_to_litr2c
   hrv_frootc_to_litr3c             => cc13f%hrv_frootc_to_litr3c
   hrv_livestemc_to_cwdc            => cc13f%hrv_livestemc_to_cwdc
   chrv_deadstemc_to_prod10c        => cc13f%hrv_deadstemc_to_prod10c
   chrv_deadstemc_to_prod100c       => cc13f%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_cwdc           => cc13f%hrv_livecrootc_to_cwdc
   hrv_deadcrootc_to_cwdc           => cc13f%hrv_deadcrootc_to_cwdc
   hrv_leafc_storage_to_litr1c      => cc13f%hrv_leafc_storage_to_litr1c
   hrv_frootc_storage_to_litr1c     => cc13f%hrv_frootc_storage_to_litr1c
   hrv_livestemc_storage_to_litr1c  => cc13f%hrv_livestemc_storage_to_litr1c
   hrv_deadstemc_storage_to_litr1c  => cc13f%hrv_deadstemc_storage_to_litr1c
   hrv_livecrootc_storage_to_litr1c => cc13f%hrv_livecrootc_storage_to_litr1c
   hrv_deadcrootc_storage_to_litr1c => cc13f%hrv_deadcrootc_storage_to_litr1c
   hrv_gresp_storage_to_litr1c      => cc13f%hrv_gresp_storage_to_litr1c
   hrv_leafc_xfer_to_litr1c         => cc13f%hrv_leafc_xfer_to_litr1c
   hrv_frootc_xfer_to_litr1c        => cc13f%hrv_frootc_xfer_to_litr1c
   hrv_livestemc_xfer_to_litr1c     => cc13f%hrv_livestemc_xfer_to_litr1c
   hrv_deadstemc_xfer_to_litr1c     => cc13f%hrv_deadstemc_xfer_to_litr1c
   hrv_livecrootc_xfer_to_litr1c    => cc13f%hrv_livecrootc_xfer_to_litr1c
   hrv_deadcrootc_xfer_to_litr1c    => cc13f%hrv_deadcrootc_xfer_to_litr1c
   hrv_gresp_xfer_to_litr1c         => cc13f%hrv_gresp_xfer_to_litr1c

   ! assign local pointers to pft-level arrays
   ivt                            => pft%itype
   wtcol                          => pft%wtcol
   pwtgcell                       => pft%wtgcell  
   hrv_leafc_to_litter              => pc13f%hrv_leafc_to_litter
   hrv_frootc_to_litter             => pc13f%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => pc13f%hrv_livestemc_to_litter
   phrv_deadstemc_to_prod10c        => pc13f%hrv_deadstemc_to_prod10c
   phrv_deadstemc_to_prod100c       => pc13f%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => pc13f%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => pc13f%hrv_deadcrootc_to_litter
   hrv_leafc_storage_to_litter      => pc13f%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => pc13f%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => pc13f%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => pc13f%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => pc13f%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => pc13f%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => pc13f%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => pc13f%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => pc13f%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => pc13f%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => pc13f%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => pc13f%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => pc13f%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => pc13f%hrv_gresp_xfer_to_litter

   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1

            if (pwtgcell(p)>0._r8) then

               ! leaf harvest mortality carbon fluxes
               hrv_leafc_to_litr1c(c) = hrv_leafc_to_litr1c(c) + &
                  hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               hrv_leafc_to_litr2c(c) = hrv_leafc_to_litr2c(c) + &
                  hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               hrv_leafc_to_litr3c(c) = hrv_leafc_to_litr3c(c) + &
                  hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root harvest mortality carbon fluxes
               hrv_frootc_to_litr1c(c) = hrv_frootc_to_litr1c(c) + &
                  hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               hrv_frootc_to_litr2c(c) = hrv_frootc_to_litr2c(c) + &
                  hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               hrv_frootc_to_litr3c(c) = hrv_frootc_to_litr3c(c) + &
                  hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood harvest mortality carbon fluxes
               hrv_livestemc_to_cwdc(c)  = hrv_livestemc_to_cwdc(c)  + &
                  hrv_livestemc_to_litter(p)  * wtcol(p)
               chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
                  phrv_deadstemc_to_prod10c(p)  * wtcol(p)
               chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
                  phrv_deadstemc_to_prod100c(p)  * wtcol(p)
               hrv_livecrootc_to_cwdc(c) = hrv_livecrootc_to_cwdc(c) + &
                  hrv_livecrootc_to_litter(p) * wtcol(p)
               hrv_deadcrootc_to_cwdc(c) = hrv_deadcrootc_to_cwdc(c) + &
                  hrv_deadcrootc_to_litter(p) * wtcol(p)

               ! storage harvest mortality carbon fluxes
               hrv_leafc_storage_to_litr1c(c)      = hrv_leafc_storage_to_litr1c(c)      + &
                  hrv_leafc_storage_to_litter(p)      * wtcol(p)
               hrv_frootc_storage_to_litr1c(c)     = hrv_frootc_storage_to_litr1c(c)     + &
                  hrv_frootc_storage_to_litter(p)     * wtcol(p)
               hrv_livestemc_storage_to_litr1c(c)  = hrv_livestemc_storage_to_litr1c(c)  + &
                  hrv_livestemc_storage_to_litter(p)  * wtcol(p)
               hrv_deadstemc_storage_to_litr1c(c)  = hrv_deadstemc_storage_to_litr1c(c)  + &
                  hrv_deadstemc_storage_to_litter(p)  * wtcol(p)
               hrv_livecrootc_storage_to_litr1c(c) = hrv_livecrootc_storage_to_litr1c(c) + &
                  hrv_livecrootc_storage_to_litter(p) * wtcol(p)
               hrv_deadcrootc_storage_to_litr1c(c) = hrv_deadcrootc_storage_to_litr1c(c) + &
                  hrv_deadcrootc_storage_to_litter(p) * wtcol(p)
               hrv_gresp_storage_to_litr1c(c)      = hrv_gresp_storage_to_litr1c(c)      + &
                  hrv_gresp_storage_to_litter(p)      * wtcol(p)

               ! transfer harvest mortality carbon fluxes
               hrv_leafc_xfer_to_litr1c(c)      = hrv_leafc_xfer_to_litr1c(c)      + &
                  hrv_leafc_xfer_to_litter(p)      * wtcol(p)
               hrv_frootc_xfer_to_litr1c(c)     = hrv_frootc_xfer_to_litr1c(c)     + &
                  hrv_frootc_xfer_to_litter(p)     * wtcol(p)
               hrv_livestemc_xfer_to_litr1c(c)  = hrv_livestemc_xfer_to_litr1c(c)  + &
                  hrv_livestemc_xfer_to_litter(p)  * wtcol(p)
               hrv_deadstemc_xfer_to_litr1c(c)  = hrv_deadstemc_xfer_to_litr1c(c)  + &
                  hrv_deadstemc_xfer_to_litter(p)  * wtcol(p)
               hrv_livecrootc_xfer_to_litr1c(c) = hrv_livecrootc_xfer_to_litr1c(c) + &
                  hrv_livecrootc_xfer_to_litter(p) * wtcol(p)
               hrv_deadcrootc_xfer_to_litr1c(c) = hrv_deadcrootc_xfer_to_litr1c(c) + &
                  hrv_deadcrootc_xfer_to_litter(p) * wtcol(p)
               hrv_gresp_xfer_to_litr1c(c)      = hrv_gresp_xfer_to_litr1c(c)      + &
                  hrv_gresp_xfer_to_litter(p)      * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNC13HarvestPftToColumn
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13FluxCalc
!
! !INTERFACE:
subroutine C13FluxCalc(c13_flux, ctot_flux, c13_state, ctot_state, &
	                    num, filter, frax, diag)
!
! !DESCRIPTION:
! On the radiation time step, set the 13-carbon flux
! variables (except for gap-phase mortality and fire fluxes)
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   real(r8), pointer   :: c13_flux(:)   !OUTPUT 13C flux
   real(r8), pointer   :: ctot_flux(:)  !INPUT  totC flux
   real(r8), pointer   :: c13_state(:)  !INPUT  13C state, upstream pool
   real(r8), pointer   :: ctot_state(:) !INPUT  totC state, upstream pool
   real(r8), intent(in):: frax          !fractionation factor (1 = no fractionation)
   integer, intent(in) :: num           ! number of filter members
   integer, intent(in) :: filter(:)     ! filter indices
   integer, intent(in) :: diag          !0=no diagnostics, 1=print diagnostics
!
! !CALLED FROM:
! subroutine C13Flux1
!
! !REVISION HISTORY:
!
! !OTHER LOCAL VARIABLES:
   integer :: i,f     ! indices
   real(r8) :: temp
!
   ! loop over the supplied filter
   do f = 1,num
      i = filter(f)
      if (ctot_state(i) /= 0._r8) then
      	c13_flux(i) = ctot_flux(i) * (c13_state(i)/ctot_state(i)) * frax
      else
      	c13_flux(i) = 0._r8
      end if
      
      if (diag == 1) then
      ! put diagnostic print statements here for 13C flux calculations
      end if
   end do
end subroutine C13FluxCalc
!-----------------------------------------------------------------------

end module CNC13FluxMod
 
