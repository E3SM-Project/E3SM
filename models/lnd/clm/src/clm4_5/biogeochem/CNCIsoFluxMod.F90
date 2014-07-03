module CNCIsoFluxMod
#if (defined CN)

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CIsoFluxMod
!
! !DESCRIPTION:
! Module for carbon isotopic flux variable update, non-mortality fluxes.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varpar   , only: ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools
    use abortutils  , only: endrun
    implicit none
    save
    private
!
! !PUBLIC MEMBER FUNCTIONS:
    public:: CIsoFlux1
    public:: CIsoFlux2
    public:: CIsoFlux2h
    public:: CIsoFlux3
    private:: CNCIsoLitterToColumn
    private:: CNCIsoGapPftToColumn
    private:: CNCIsoHarvestPftToColumn
    private:: CIsoFluxCalc
!
! !REVISION HISTORY:
! 4/21/2005: Created by Peter Thornton and Neil Suits
! 8/15/2011: Modified by C. Koven from C13Flux to CisoFlux to allow treatment of both C13 and C14 
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CIsoFlux1
!
! !INTERFACE:
subroutine CIsoFlux1(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
!
! !DESCRIPTION:
! On the radiation time step, set the carbon isotopic flux
! variables (except for gap-phase mortality and fire fluxes)
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc         ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:)   ! filter for soil columns
   integer, intent(in) :: num_soilp         ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:)   ! filter for soil pfts
   character(len=*), intent(in) :: isotope  ! 'c13' or 'c14'
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
   type(pft_cflux_type), pointer :: pcisof
   type(pft_cstate_type), pointer :: pcisos
   type(column_cflux_type), pointer :: ccisof
   type(column_cstate_type), pointer :: ccisos
   integer :: fp,pi,l
   integer :: fc,cc,j

   integer,  pointer :: cascade_donor_pool(:)             ! which pool is C taken from for a given decomposition step
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => clm3%g%l%c%p
   c => clm3%g%l%c
   select case (isotope)
   case ('c14')
      pcisof => p%pc14f
      pcisos => p%pc14s
      ccisof => c%cc14f
      ccisos => c%cc14s
   case ('c13')
      pcisof => p%pc13f
      pcisos => p%pc13s
      ccisof => c%cc13f
      ccisos => c%cc13s
   case default
      call endrun('CNCIsoFluxMod: iso must be either c13 or c14')
   end select
   cascade_donor_pool                => decomp_cascade_con%cascade_donor_pool
	
   ! pft-level non-mortality fluxes
   
   call CIsoFluxCalc(pcisof%leafc_xfer_to_leafc, p%pcf%leafc_xfer_to_leafc, &
                    pcisos%leafc_xfer, p%pcs%leafc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%frootc_xfer_to_frootc, p%pcf%frootc_xfer_to_frootc, &
                    pcisos%frootc_xfer, p%pcs%frootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livestemc_xfer_to_livestemc, p%pcf%livestemc_xfer_to_livestemc, &
                    pcisos%livestemc_xfer, p%pcs%livestemc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%deadstemc_xfer_to_deadstemc, p%pcf%deadstemc_xfer_to_deadstemc, &
                    pcisos%deadstemc_xfer, p%pcs%deadstemc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livecrootc_xfer_to_livecrootc, p%pcf%livecrootc_xfer_to_livecrootc, &
                    pcisos%livecrootc_xfer, p%pcs%livecrootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%deadcrootc_xfer_to_deadcrootc, p%pcf%deadcrootc_xfer_to_deadcrootc, &
                    pcisos%deadcrootc_xfer, p%pcs%deadcrootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%leafc_to_litter, p%pcf%leafc_to_litter, &
                    pcisos%leafc, p%pcs%leafc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%frootc_to_litter, p%pcf%frootc_to_litter, &
                    pcisos%frootc, p%pcs%frootc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livestemc_to_deadstemc, p%pcf%livestemc_to_deadstemc, &
                    pcisos%livestemc, p%pcs%livestemc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livecrootc_to_deadcrootc, p%pcf%livecrootc_to_deadcrootc, &
                    pcisos%livecrootc, p%pcs%livecrootc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%leaf_curmr, p%pcf%leaf_curmr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%froot_curmr, p%pcf%froot_curmr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livestem_curmr, p%pcf%livestem_curmr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livecroot_curmr, p%pcf%livecroot_curmr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%leaf_xsmr, p%pcf%leaf_xsmr, &
                    pcisos%totvegc, p%pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%froot_xsmr, p%pcf%froot_xsmr, &
                    pcisos%totvegc, p%pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livestem_xsmr, p%pcf%livestem_xsmr, &
                    pcisos%totvegc, p%pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livecroot_xsmr, p%pcf%livecroot_xsmr, &
                    pcisos%totvegc, p%pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_xsmrpool, p%pcf%cpool_to_xsmrpool, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_leafc, p%pcf%cpool_to_leafc, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_leafc_storage, p%pcf%cpool_to_leafc_storage, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_frootc, p%pcf%cpool_to_frootc, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_frootc_storage, p%pcf%cpool_to_frootc_storage, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_livestemc, p%pcf%cpool_to_livestemc, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_livestemc_storage, p%pcf%cpool_to_livestemc_storage, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_deadstemc, p%pcf%cpool_to_deadstemc, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_deadstemc_storage, p%pcf%cpool_to_deadstemc_storage, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_livecrootc, p%pcf%cpool_to_livecrootc, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_livecrootc_storage, p%pcf%cpool_to_livecrootc_storage, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_deadcrootc, p%pcf%cpool_to_deadcrootc, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_deadcrootc_storage, p%pcf%cpool_to_deadcrootc_storage, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_leaf_gr, p%pcf%cpool_leaf_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_froot_gr, p%pcf%cpool_froot_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_livestem_gr, p%pcf%cpool_livestem_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_deadstem_gr, p%pcf%cpool_deadstem_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_livecroot_gr, p%pcf%cpool_livecroot_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_deadcroot_gr, p%pcf%cpool_deadcroot_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_leaf_storage_gr, p%pcf%cpool_leaf_storage_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_froot_storage_gr, p%pcf%cpool_froot_storage_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_livestem_storage_gr, p%pcf%cpool_livestem_storage_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_deadstem_storage_gr, p%pcf%cpool_deadstem_storage_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_livecroot_storage_gr, p%pcf%cpool_livecroot_storage_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_deadcroot_storage_gr, p%pcf%cpool_deadcroot_storage_gr, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_gresp_storage, p%pcf%cpool_to_gresp_storage, &
                    pcisos%cpool, p%pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_leaf_gr, p%pcf%transfer_leaf_gr, &
                    pcisos%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_froot_gr, p%pcf%transfer_froot_gr, &
                    pcisos%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_livestem_gr, p%pcf%transfer_livestem_gr, &
                    pcisos%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_deadstem_gr, p%pcf%transfer_deadstem_gr, &
                    pcisos%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_livecroot_gr, p%pcf%transfer_livecroot_gr, &
                    pcisos%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_deadcroot_gr, p%pcf%transfer_deadcroot_gr, &
                    pcisos%gresp_xfer, p%pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%leafc_storage_to_xfer, p%pcf%leafc_storage_to_xfer, &
                    pcisos%leafc_storage, p%pcs%leafc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%frootc_storage_to_xfer, p%pcf%frootc_storage_to_xfer, &
                    pcisos%frootc_storage, p%pcs%frootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livestemc_storage_to_xfer, p%pcf%livestemc_storage_to_xfer, &
                    pcisos%livestemc_storage, p%pcs%livestemc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%deadstemc_storage_to_xfer, p%pcf%deadstemc_storage_to_xfer, &
                    pcisos%deadstemc_storage, p%pcs%deadstemc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livecrootc_storage_to_xfer, p%pcf%livecrootc_storage_to_xfer, &
                    pcisos%livecrootc_storage, p%pcs%livecrootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%deadcrootc_storage_to_xfer, p%pcf%deadcrootc_storage_to_xfer, &
                    pcisos%deadcrootc_storage, p%pcs%deadcrootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%gresp_storage_to_xfer, p%pcf%gresp_storage_to_xfer, &
                    pcisos%gresp_storage, p%pcs%gresp_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	! call routine to shift pft-level litterfall fluxes to column, for isotopes
   ! the non-isotope version of this routine is called in CNPhenologyMod.F90
   ! For later clean-up, it would be possible to generalize this function to operate on a single 
   ! pft-to-column flux.
   
   call CNCIsoLitterToColumn(num_soilc, filter_soilc, isotope)
   
   ! column-level non-mortality fluxes

   do fc = 1,num_soilc
      cc = filter_soilc(fc)
      do j = 1, nlevdecomp
         do l = 1, ndecomp_cascade_transitions
            if ( c%ccs%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) /= 0._r8) then
               ccisof%decomp_cascade_hr_vr(cc,j,l)  =  c%ccf%decomp_cascade_hr_vr(cc,j,l) * &
                    (ccisos%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) / c%ccs%decomp_cpools_vr(cc,j,cascade_donor_pool(l))) * 1._r8
            else
               ccisof%decomp_cascade_hr_vr(cc,j,l) = 0._r8
            end if
         end do
      end do
   end do

   do fc = 1,num_soilc
      cc = filter_soilc(fc)
      do j = 1, nlevdecomp
         do l = 1, ndecomp_cascade_transitions
            if ( c%ccs%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) /= 0._r8) then
               ccisof%decomp_cascade_ctransfer_vr(cc,j,l)  =  c%ccf%decomp_cascade_ctransfer_vr(cc,j,l) * &
                    (ccisos%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) / c%ccs%decomp_cpools_vr(cc,j,cascade_donor_pool(l))) * 1._r8
            else
               ccisof%decomp_cascade_ctransfer_vr(cc,j,l) = 0._r8
            end if
         end do
      end do
   end do


end subroutine CIsoFlux1
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CIsoFlux2
!
! !INTERFACE:
subroutine CIsoFlux2(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
!
! !DESCRIPTION:
! On the radiation time step, set the carbon isotopic fluxes for gap mortality
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc         ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:)   ! filter for soil columns
   integer, intent(in) :: num_soilp         ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:)   ! filter for soil pfts
   character(len=*), intent(in) :: isotope  ! 'c13' or 'c14'
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   type(pft_cflux_type), pointer :: pcisof
   type(pft_cstate_type), pointer :: pcisos
   integer :: fp,pi
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => clm3%g%l%c%p
   select case (isotope)
   case ('c14')
      pcisof => p%pc14f
      pcisos => p%pc14s
   case ('c13')
      pcisof => p%pc13f
      pcisos => p%pc13s
   case default
      call endrun('CNCIsoFluxMod: iso must be either c13 or c14')
   end select

   ! pft-level gap mortality fluxes
   
   call CIsoFluxCalc(pcisof%m_leafc_to_litter, p%pcf%m_leafc_to_litter, &
                     pcisos%leafc, p%pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_leafc_storage_to_litter, p%pcf%m_leafc_storage_to_litter, &
                     pcisos%leafc_storage, p%pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_leafc_xfer_to_litter, p%pcf%m_leafc_xfer_to_litter, &
                     pcisos%leafc_xfer, p%pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_to_litter, p%pcf%m_frootc_to_litter, &
                     pcisos%frootc, p%pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_storage_to_litter, p%pcf%m_frootc_storage_to_litter, &
                     pcisos%frootc_storage, p%pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_xfer_to_litter, p%pcf%m_frootc_xfer_to_litter, &
                     pcisos%frootc_xfer, p%pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_to_litter, p%pcf%m_livestemc_to_litter, &
                     pcisos%livestemc, p%pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_storage_to_litter, p%pcf%m_livestemc_storage_to_litter, &
                     pcisos%livestemc_storage, p%pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_xfer_to_litter, p%pcf%m_livestemc_xfer_to_litter, &
                     pcisos%livestemc_xfer, p%pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_to_litter, p%pcf%m_deadstemc_to_litter, &
                     pcisos%deadstemc, p%pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_storage_to_litter, p%pcf%m_deadstemc_storage_to_litter, &
                     pcisos%deadstemc_storage, p%pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_xfer_to_litter, p%pcf%m_deadstemc_xfer_to_litter, &
                     pcisos%deadstemc_xfer, p%pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_to_litter, p%pcf%m_livecrootc_to_litter, &
                     pcisos%livecrootc, p%pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_storage_to_litter, p%pcf%m_livecrootc_storage_to_litter, &
                     pcisos%livecrootc_storage, p%pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_xfer_to_litter, p%pcf%m_livecrootc_xfer_to_litter, &
                     pcisos%livecrootc_xfer, p%pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_to_litter, p%pcf%m_deadcrootc_to_litter, &
                     pcisos%deadcrootc, p%pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_storage_to_litter, p%pcf%m_deadcrootc_storage_to_litter, &
                     pcisos%deadcrootc_storage, p%pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_xfer_to_litter, p%pcf%m_deadcrootc_xfer_to_litter, &
                     pcisos%deadcrootc_xfer, p%pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_gresp_storage_to_litter, p%pcf%m_gresp_storage_to_litter, &
                     pcisos%gresp_storage, p%pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_gresp_xfer_to_litter, p%pcf%m_gresp_xfer_to_litter, &
                     pcisos%gresp_xfer, p%pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
	! call routine to shift pft-level gap mortality fluxes to column, for isotopes
   ! the non-isotope version of this routine is in CNGapMortalityMod.F90.

	call CNCIsoGapPftToColumn(num_soilc, filter_soilc, isotope)
   
end subroutine CIsoFlux2
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CIsoFlux2h
!
! !INTERFACE:
subroutine CIsoFlux2h(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
!
! !DESCRIPTION:
! set the carbon isotopic fluxes for harvest mortality
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc        ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:)  ! filter for soil columns
   integer, intent(in) :: num_soilp        ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:)  ! filter for soil pfts
   character(len=*), intent(in) :: isotope ! 'c13' or 'c14'
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
!
! !LOCAL VARIABLES:
! !OTHER LOCAL VARIABLES:
   type(pft_type), pointer :: p
   type(pft_cflux_type), pointer :: pcisof
   type(pft_cstate_type), pointer :: pcisos
   integer :: fp,pi
!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => clm3%g%l%c%p
   select case (isotope)
   case ('c14')
      pcisof => p%pc14f
      pcisos => p%pc14s
   case ('c13')
      pcisof => p%pc13f
      pcisos => p%pc13s
   case default
      call endrun('CNCIsoFluxMod: iso must be either c13 or c14')
   end select
	
   ! pft-level gap mortality fluxes
   
   call CIsoFluxCalc(pcisof%hrv_leafc_to_litter, p%pcf%hrv_leafc_to_litter, &
                     pcisos%leafc, p%pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_leafc_storage_to_litter, p%pcf%hrv_leafc_storage_to_litter, &
                     pcisos%leafc_storage, p%pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_leafc_xfer_to_litter, p%pcf%hrv_leafc_xfer_to_litter, &
                     pcisos%leafc_xfer, p%pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_frootc_to_litter, p%pcf%hrv_frootc_to_litter, &
                     pcisos%frootc, p%pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_frootc_storage_to_litter, p%pcf%hrv_frootc_storage_to_litter, &
                     pcisos%frootc_storage, p%pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_frootc_xfer_to_litter, p%pcf%hrv_frootc_xfer_to_litter, &
                     pcisos%frootc_xfer, p%pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livestemc_to_litter, p%pcf%hrv_livestemc_to_litter, &
                     pcisos%livestemc, p%pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livestemc_storage_to_litter, p%pcf%hrv_livestemc_storage_to_litter, &
                     pcisos%livestemc_storage, p%pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livestemc_xfer_to_litter, p%pcf%hrv_livestemc_xfer_to_litter, &
                     pcisos%livestemc_xfer, p%pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadstemc_to_prod10c, p%pcf%hrv_deadstemc_to_prod10c, &
                     pcisos%deadstemc, p%pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadstemc_to_prod100c, p%pcf%hrv_deadstemc_to_prod100c, &
                     pcisos%deadstemc, p%pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadstemc_storage_to_litter, p%pcf%hrv_deadstemc_storage_to_litter, &
                     pcisos%deadstemc_storage, p%pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadstemc_xfer_to_litter, p%pcf%hrv_deadstemc_xfer_to_litter, &
                     pcisos%deadstemc_xfer, p%pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livecrootc_to_litter, p%pcf%hrv_livecrootc_to_litter, &
                     pcisos%livecrootc, p%pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livecrootc_storage_to_litter, p%pcf%hrv_livecrootc_storage_to_litter, &
                     pcisos%livecrootc_storage, p%pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livecrootc_xfer_to_litter, p%pcf%hrv_livecrootc_xfer_to_litter, &
                     pcisos%livecrootc_xfer, p%pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadcrootc_to_litter, p%pcf%hrv_deadcrootc_to_litter, &
                     pcisos%deadcrootc, p%pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadcrootc_storage_to_litter, p%pcf%hrv_deadcrootc_storage_to_litter, &
                     pcisos%deadcrootc_storage, p%pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadcrootc_xfer_to_litter, p%pcf%hrv_deadcrootc_xfer_to_litter, &
                     pcisos%deadcrootc_xfer, p%pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_gresp_storage_to_litter, p%pcf%hrv_gresp_storage_to_litter, &
                     pcisos%gresp_storage, p%pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_gresp_xfer_to_litter, p%pcf%hrv_gresp_xfer_to_litter, &
                     pcisos%gresp_xfer, p%pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_xsmrpool_to_atm, p%pcf%hrv_xsmrpool_to_atm, &
                    pcisos%totvegc, p%pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
   ! call routine to shift pft-level gap mortality fluxes to column, for isotopes
   ! the non-isotope version of this routine is in CNGapMortalityMod.F90.

   call CNCIsoHarvestPftToColumn(num_soilc, filter_soilc, isotope)
   
end subroutine CIsoFlux2h
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CIsoFlux3
!
! !INTERFACE:
subroutine CIsoFlux3(num_soilc, filter_soilc, num_soilp, filter_soilp, isotope)
!
! !DESCRIPTION:
! On the radiation time step, set the carbon isotopic fluxes for fire mortality
!
! !USES:
   use clmtype
   use clm_varpar, only : max_pft_per_col

!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc         ! number of soil columns filter
   integer, intent(in) :: filter_soilc(:)   ! filter for soil columns
   integer, intent(in) :: num_soilp         ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:)   ! filter for soil pfts
   character(len=*), intent(in) :: isotope  ! 'c13' or 'c14'
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
   type(pft_cflux_type), pointer :: pcisof
   type(pft_cstate_type), pointer :: pcisos
   type(column_cflux_type), pointer :: ccisof
   type(column_cstate_type), pointer :: ccisos
   integer :: fp,pi,l,pp
   integer :: fc,cc,j
   real(r8), pointer :: ptrp(:)         ! pointer to input pft array
   real(r8), pointer :: ptrc(:)         ! pointer to output column array

   real(r8), pointer :: croot_prof(:,:) ! (1/m) profile of coarse roots
   real(r8), pointer :: stem_prof(:,:)  ! (1/m) profile of stems
   integer , pointer :: npfts(:)        ! number of pfts for each column
   integer , pointer :: pfti(:)         ! beginning pft index for each column
   real(r8), pointer :: wtcol(:)        ! weight (relative to column) for this pft (0-1)
   logical , pointer :: pactive(:)      ! true=>do computations on this pft (see reweightMod for details)


!
!EOP
!-----------------------------------------------------------------------
	! set local pointers
   p => clm3%g%l%c%p
   c => clm3%g%l%c
   select case (isotope)
   case ('c14')
      pcisof => p%pc14f
      pcisos => p%pc14s
      ccisof => c%cc14f
      ccisos => c%cc14s
   case ('c13')
      pcisof => p%pc13f
      pcisos => p%pc13s
      ccisof => c%cc13f
      ccisos => c%cc13s
   case default
      call endrun('CNCIsoFluxMod: iso must be either c13 or c14')
   end select
   croot_prof                     => clm3%g%l%c%p%pps%croot_prof
   stem_prof                      => clm3%g%l%c%p%pps%stem_prof
   npfts                          => clm3%g%l%c%npfts
   pfti                           => clm3%g%l%c%pfti
   wtcol                          => clm3%g%l%c%p%wtcol
   pactive                        => clm3%g%l%c%p%active

	
   ! pft-level fire mortality fluxes
   
   call CIsoFluxCalc(pcisof%m_leafc_to_fire, p%pcf%m_leafc_to_fire, &
                     pcisos%leafc, p%pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_leafc_storage_to_fire, p%pcf%m_leafc_storage_to_fire, &
                     pcisos%leafc_storage, p%pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_leafc_xfer_to_fire, p%pcf%m_leafc_xfer_to_fire, &
                     pcisos%leafc_xfer, p%pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_to_fire, p%pcf%m_frootc_to_fire, &
                     pcisos%frootc, p%pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_storage_to_fire, p%pcf%m_frootc_storage_to_fire, &
                     pcisos%frootc_storage, p%pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_xfer_to_fire, p%pcf%m_frootc_xfer_to_fire, &
                     pcisos%frootc_xfer, p%pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_to_fire, p%pcf%m_livestemc_to_fire, &
                     pcisos%livestemc, p%pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_storage_to_fire, p%pcf%m_livestemc_storage_to_fire, &
                     pcisos%livestemc_storage, p%pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_xfer_to_fire, p%pcf%m_livestemc_xfer_to_fire, &
                     pcisos%livestemc_xfer, p%pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_to_fire, p%pcf%m_deadstemc_to_fire, &
                     pcisos%deadstemc, p%pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_to_litter_fire, p%pcf%m_deadstemc_to_litter_fire, &
                     pcisos%deadstemc, p%pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_storage_to_fire, p%pcf%m_deadstemc_storage_to_fire, &
                     pcisos%deadstemc_storage, p%pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_xfer_to_fire, p%pcf%m_deadstemc_xfer_to_fire, &
                     pcisos%deadstemc_xfer, p%pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_to_fire, p%pcf%m_livecrootc_to_fire, &
                     pcisos%livecrootc, p%pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_storage_to_fire, p%pcf%m_livecrootc_storage_to_fire, &
                     pcisos%livecrootc_storage, p%pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_xfer_to_fire, p%pcf%m_livecrootc_xfer_to_fire, &
                     pcisos%livecrootc_xfer, p%pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_to_fire, p%pcf%m_deadcrootc_to_fire, &
                     pcisos%deadcrootc, p%pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_to_litter_fire, p%pcf%m_deadcrootc_to_litter_fire, &
                     pcisos%deadcrootc, p%pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_storage_to_fire, p%pcf%m_deadcrootc_storage_to_fire, &
                     pcisos%deadcrootc_storage, p%pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_xfer_to_fire, p%pcf%m_deadcrootc_xfer_to_fire, &
                     pcisos%deadcrootc_xfer, p%pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_gresp_storage_to_fire, p%pcf%m_gresp_storage_to_fire, &
                     pcisos%gresp_storage, p%pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_gresp_xfer_to_fire, p%pcf%m_gresp_xfer_to_fire, &
                     pcisos%gresp_xfer, p%pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   ! calculate the column-level flux of deadstem and deadcrootc to cwdc as the result of fire mortality.
   do pi = 1,max_pft_per_col
      do fc = 1,num_soilc
         cc = filter_soilc(fc)
         if ( pi <=  npfts(cc) ) then
            pp = pfti(cc) + pi - 1
            if (pactive(pp)) then
               do j = 1, nlevdecomp
                  ccisof%fire_mortality_c_to_cwdc(cc,j) = ccisof%fire_mortality_c_to_cwdc(cc,j) + &
                       pcisof%m_deadstemc_to_litter_fire(pp) * wtcol(pp) * stem_prof(pp,j)
                  ccisof%fire_mortality_c_to_cwdc(cc,j) = ccisof%fire_mortality_c_to_cwdc(cc,j) + &
                       pcisof%m_deadcrootc_to_litter_fire(pp) * wtcol(pp) * croot_prof(pp,j)
               end do
            end if
         end if
      end do
   end do


   do fc = 1,num_soilc
      cc = filter_soilc(fc)
      do j = 1, nlevdecomp
         do l = 1, ndecomp_pools
            if ( c%ccs%decomp_cpools_vr(cc,j,l) /= 0._r8) then
               ccisof%m_decomp_cpools_to_fire_vr(cc,j,l)  =  c%ccf%m_decomp_cpools_to_fire_vr(cc,j,l) * &
                    (ccisos%decomp_cpools_vr(cc,j,l) / c%ccs%decomp_cpools_vr(cc,j,l)) * 1._r8
            else
               ccisof%m_decomp_cpools_to_fire_vr(cc,j,l) = 0._r8
            end if
         end do
      end do
   end do


                    
end subroutine CIsoFlux3
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNCIsoLitterToColumn
!
! !INTERFACE:
subroutine CNCIsoLitterToColumn (num_soilc, filter_soilc, isotope)
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
  integer, intent(in) :: num_soilc        ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:)  ! filter for soil columns
  character(len=*), intent(in) :: isotope ! 'c13' or 'c14'
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
   logical , pointer :: pactive(:)      ! true=>do computations on this pft (see reweightMod for details)
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
   real(r8), pointer :: phenology_c_to_litr_met_c(:,:)             ! C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
   real(r8), pointer :: phenology_c_to_litr_cel_c(:,:)             ! C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
   real(r8), pointer :: phenology_c_to_litr_lig_c(:,:)             ! C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)

   real(r8), pointer :: leaf_prof(:,:)          ! (1/m) profile of leaves
   real(r8), pointer :: froot_prof(:,:)         ! (1/m) profile of fine roots
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   type(pft_cflux_type), pointer :: pcisof
   type(column_cflux_type), pointer :: ccisof
   integer :: fc,c,pi,p,j
!EOP
!-----------------------------------------------------------------------
    ! select which isotope
    select case (isotope)
    case ('c14')
       pcisof => clm3%g%l%c%p%pc14f
       ccisof => clm3%g%l%c%cc14f
    case ('c13')
       pcisof => clm3%g%l%c%p%pc13f
       ccisof => clm3%g%l%c%cc13f
    case default
       call endrun('CNCIsoFluxMod: iso must be either c13 or c14')
    end select

   ! assign local pointers to derived type arrays (in)
    ivt                            => clm3%g%l%c%p%itype
    wtcol                          => clm3%g%l%c%p%wtcol
    pactive                        => clm3%g%l%c%p%active
    leafc_to_litter                => pcisof%leafc_to_litter
    frootc_to_litter               => pcisof%frootc_to_litter
    npfts                          => clm3%g%l%c%npfts
    pfti                           => clm3%g%l%c%pfti
    lf_flab                        => pftcon%lf_flab
    lf_fcel                        => pftcon%lf_fcel
    lf_flig                        => pftcon%lf_flig
    fr_flab                        => pftcon%fr_flab
    fr_fcel                        => pftcon%fr_fcel
    fr_flig                        => pftcon%fr_flig

   ! assign local pointers to derived type arrays (out)
    phenology_c_to_litr_met_c         => ccisof%phenology_c_to_litr_met_c
    phenology_c_to_litr_cel_c         => ccisof%phenology_c_to_litr_cel_c
    phenology_c_to_litr_lig_c         => ccisof%phenology_c_to_litr_lig_c

    leaf_prof                      => clm3%g%l%c%p%pps%leaf_prof
    froot_prof                     => clm3%g%l%c%p%pps%froot_prof

    do j = 1, nlevdecomp
       do pi = 1,max_pft_per_col
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             
             if ( pi <=  npfts(c) ) then
                p = pfti(c) + pi - 1
                if (pactive(p)) then
                   ! leaf litter carbon fluxes
                   phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) + leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) + leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) + leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   
                   ! fine root litter carbon fluxes
                   phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) + frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                   phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) + frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                   phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) + frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)
                end if
             end if
             
          end do
       end do
       
    end do
    
  end subroutine CNCIsoLitterToColumn
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNCIsoGapPftToColumn
!
! !INTERFACE:
subroutine CNCIsoGapPftToColumn (num_soilc, filter_soilc, isotope)
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
  integer, intent(in) :: num_soilc         ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:)   ! soil column filter
   character(len=*), intent(in) :: isotope ! 'c13' or 'c14'
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
   logical , pointer :: pactive(:)  ! true=>do computations on this pft (see reweightMod for details)
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
   real(r8), pointer :: gap_mortality_c_to_litr_met_c(:,:)         ! C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
   real(r8), pointer :: gap_mortality_c_to_litr_cel_c(:,:)         ! C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
   real(r8), pointer :: gap_mortality_c_to_litr_lig_c(:,:)         ! C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
   real(r8), pointer :: gap_mortality_c_to_cwdc(:,:)               ! C fluxes associated with gap mortality to CWD pool (gC/m3/s)
   real(r8), pointer :: leaf_prof(:,:)          ! (1/m) profile of leaves
   real(r8), pointer :: froot_prof(:,:)         ! (1/m) profile of fine roots
   real(r8), pointer :: croot_prof(:,:)         ! (1/m) profile of coarse roots
   real(r8), pointer :: stem_prof(:,:)          ! (1/m) profile of stems
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   type(pft_cflux_type), pointer :: pcisof
   type(column_cflux_type), pointer :: ccisof
   integer :: fc,c,pi,p,j               ! indices
!EOP
!-----------------------------------------------------------------------

    ! select which isotope
    select case (isotope)
    case ('c14')
       pcisof => clm3%g%l%c%p%pc14f
       ccisof => clm3%g%l%c%cc14f
    case ('c13')
       pcisof => clm3%g%l%c%p%pc13f
       ccisof => clm3%g%l%c%cc13f
    case default
       call endrun('CNCIsoFluxMod: iso must be either c13 or c14')
    end select

   ! assign local pointers
   lf_flab                        => pftcon%lf_flab
   lf_fcel                        => pftcon%lf_fcel
   lf_flig                        => pftcon%lf_flig
   fr_flab                        => pftcon%fr_flab
   fr_fcel                        => pftcon%fr_fcel
   fr_flig                        => pftcon%fr_flig

   ! assign local pointers to column-level arrays
   npfts                          => clm3%g%l%c%npfts
   pfti                           => clm3%g%l%c%pfti
   gap_mortality_c_to_litr_met_c  => ccisof%gap_mortality_c_to_litr_met_c
   gap_mortality_c_to_litr_cel_c  => ccisof%gap_mortality_c_to_litr_cel_c
   gap_mortality_c_to_litr_lig_c  => ccisof%gap_mortality_c_to_litr_lig_c
   gap_mortality_c_to_cwdc        => ccisof%gap_mortality_c_to_cwdc

   ! assign local pointers to pft-level arrays
   ivt                            => clm3%g%l%c%p%itype
   wtcol                          => clm3%g%l%c%p%wtcol
   pactive                        => clm3%g%l%c%p%active
   m_leafc_to_litter              => pcisof%m_leafc_to_litter
   m_frootc_to_litter             => pcisof%m_frootc_to_litter
   m_livestemc_to_litter          => pcisof%m_livestemc_to_litter
   m_deadstemc_to_litter          => pcisof%m_deadstemc_to_litter
   m_livecrootc_to_litter         => pcisof%m_livecrootc_to_litter
   m_deadcrootc_to_litter         => pcisof%m_deadcrootc_to_litter
   m_leafc_storage_to_litter      => pcisof%m_leafc_storage_to_litter
   m_frootc_storage_to_litter     => pcisof%m_frootc_storage_to_litter
   m_livestemc_storage_to_litter  => pcisof%m_livestemc_storage_to_litter
   m_deadstemc_storage_to_litter  => pcisof%m_deadstemc_storage_to_litter
   m_livecrootc_storage_to_litter => pcisof%m_livecrootc_storage_to_litter
   m_deadcrootc_storage_to_litter => pcisof%m_deadcrootc_storage_to_litter
   m_gresp_storage_to_litter      => pcisof%m_gresp_storage_to_litter
   m_leafc_xfer_to_litter         => pcisof%m_leafc_xfer_to_litter
   m_frootc_xfer_to_litter        => pcisof%m_frootc_xfer_to_litter
   m_livestemc_xfer_to_litter     => pcisof%m_livestemc_xfer_to_litter
   m_deadstemc_xfer_to_litter     => pcisof%m_deadstemc_xfer_to_litter
   m_livecrootc_xfer_to_litter    => pcisof%m_livecrootc_xfer_to_litter
   m_deadcrootc_xfer_to_litter    => pcisof%m_deadcrootc_xfer_to_litter
   m_gresp_xfer_to_litter         => pcisof%m_gresp_xfer_to_litter
   leaf_prof                      => clm3%g%l%c%p%pps%leaf_prof
   froot_prof                     => clm3%g%l%c%p%pps%froot_prof
   croot_prof                     => clm3%g%l%c%p%pps%croot_prof
   stem_prof                      => clm3%g%l%c%p%pps%stem_prof

   do j = 1, nlevdecomp
      do pi = 1,maxpatch_pft
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (pi <=  npfts(c)) then
               p = pfti(c) + pi - 1
               
               if (pactive(p)) then
                  
                  ! leaf gap mortality carbon fluxes
                  gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                       m_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  gap_mortality_c_to_litr_cel_c(c,j) = gap_mortality_c_to_litr_cel_c(c,j) + &
                       m_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  gap_mortality_c_to_litr_lig_c(c,j) = gap_mortality_c_to_litr_lig_c(c,j) + &
                       m_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  
                  ! fine root gap mortality carbon fluxes
                  gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                       m_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  gap_mortality_c_to_litr_cel_c(c,j) = gap_mortality_c_to_litr_cel_c(c,j) + &
                       m_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  gap_mortality_c_to_litr_lig_c(c,j) = gap_mortality_c_to_litr_lig_c(c,j) + &
                       m_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  
                  ! wood gap mortality carbon fluxes
                  gap_mortality_c_to_cwdc(c,j)  = gap_mortality_c_to_cwdc(c,j)  + &
                       m_livestemc_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_c_to_cwdc(c,j)  = gap_mortality_c_to_cwdc(c,j)  + &
                       m_deadstemc_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                       m_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                       m_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  
                  ! storage gap mortality carbon fluxes
                  gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                       m_leafc_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)     = gap_mortality_c_to_litr_met_c(c,j)     + &
                       m_frootc_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                       m_livestemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                       m_deadstemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                       m_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                       m_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                       m_gresp_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  
                  ! transfer gap mortality carbon fluxes
                  gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                       m_leafc_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)     = gap_mortality_c_to_litr_met_c(c,j)     + &
                       m_frootc_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                       m_livestemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                       m_deadstemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                       m_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                       m_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                       m_gresp_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  
               end if
            end if
            
         end do
         
      end do
   end do

end subroutine CNCIsoGapPftToColumn
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNCIsoHarvestPftToColumn
!
! !INTERFACE:
subroutine CNCIsoHarvestPftToColumn (num_soilc, filter_soilc, isotope)
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
  integer, intent(in) :: num_soilc         ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:)   ! soil column filter
   character(len=*), intent(in) :: isotope ! 'c13' or 'c14'
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
   logical , pointer :: pactive(:)  ! true=>do computations on this pft (see reweightMod for details)
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
   real(r8), pointer :: chrv_deadstemc_to_prod10c(:)
   real(r8), pointer :: chrv_deadstemc_to_prod100c(:)
   real(r8), pointer :: harvest_c_to_litr_met_c(:,:)               ! C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
   real(r8), pointer :: harvest_c_to_litr_cel_c(:,:)               ! C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
   real(r8), pointer :: harvest_c_to_litr_lig_c(:,:)               ! C fluxes associated with harvest to litter lignin pool (gC/m3/s)
   real(r8), pointer :: harvest_c_to_cwdc(:,:)                     ! C fluxes associated with harvest to CWD pool (gC/m3/s)
   real(r8), pointer :: leaf_prof(:,:)          ! (1/m) profile of leaves
   real(r8), pointer :: froot_prof(:,:)         ! (1/m) profile of fine roots
   real(r8), pointer :: croot_prof(:,:)         ! (1/m) profile of coarse roots
   real(r8), pointer :: stem_prof(:,:)          ! (1/m) profile of stems
!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   type(pft_cflux_type), pointer :: pcisof
   type(column_cflux_type), pointer :: ccisof
   integer :: fc,c,pi,p,j               ! indices
!EOP
!-----------------------------------------------------------------------
    ! select which isotope
    select case (isotope)
    case ('c14')
       pcisof => clm3%g%l%c%p%pc14f
       ccisof => clm3%g%l%c%cc14f
    case ('c13')
       pcisof => clm3%g%l%c%p%pc13f
       ccisof => clm3%g%l%c%cc13f
    case default
       call endrun('CNCIsoFluxMod: iso must be either c13 or c14')
    end select

   ! assign local pointers
   lf_flab                        => pftcon%lf_flab
   lf_fcel                        => pftcon%lf_fcel
   lf_flig                        => pftcon%lf_flig
   fr_flab                        => pftcon%fr_flab
   fr_fcel                        => pftcon%fr_fcel
   fr_flig                        => pftcon%fr_flig

   ! assign local pointers to column-level arrays
   npfts                          => clm3%g%l%c%npfts
   pfti                           => clm3%g%l%c%pfti
   chrv_deadstemc_to_prod10c        => ccisof%hrv_deadstemc_to_prod10c
   chrv_deadstemc_to_prod100c       => ccisof%hrv_deadstemc_to_prod100c
   harvest_c_to_litr_met_c          => ccisof%harvest_c_to_litr_met_c
   harvest_c_to_litr_cel_c          => ccisof%harvest_c_to_litr_cel_c
   harvest_c_to_litr_lig_c          => ccisof%harvest_c_to_litr_lig_c
   harvest_c_to_cwdc                => ccisof%harvest_c_to_cwdc

   ! assign local pointers to pft-level arrays
   ivt                            => clm3%g%l%c%p%itype
   wtcol                          => clm3%g%l%c%p%wtcol
   pactive                        => clm3%g%l%c%p%active
   hrv_leafc_to_litter              => pcisof%hrv_leafc_to_litter
   hrv_frootc_to_litter             => pcisof%hrv_frootc_to_litter
   hrv_livestemc_to_litter          => pcisof%hrv_livestemc_to_litter
   phrv_deadstemc_to_prod10c        => pcisof%hrv_deadstemc_to_prod10c
   phrv_deadstemc_to_prod100c       => pcisof%hrv_deadstemc_to_prod100c
   hrv_livecrootc_to_litter         => pcisof%hrv_livecrootc_to_litter
   hrv_deadcrootc_to_litter         => pcisof%hrv_deadcrootc_to_litter
   hrv_leafc_storage_to_litter      => pcisof%hrv_leafc_storage_to_litter
   hrv_frootc_storage_to_litter     => pcisof%hrv_frootc_storage_to_litter
   hrv_livestemc_storage_to_litter  => pcisof%hrv_livestemc_storage_to_litter
   hrv_deadstemc_storage_to_litter  => pcisof%hrv_deadstemc_storage_to_litter
   hrv_livecrootc_storage_to_litter => pcisof%hrv_livecrootc_storage_to_litter
   hrv_deadcrootc_storage_to_litter => pcisof%hrv_deadcrootc_storage_to_litter
   hrv_gresp_storage_to_litter      => pcisof%hrv_gresp_storage_to_litter
   hrv_leafc_xfer_to_litter         => pcisof%hrv_leafc_xfer_to_litter
   hrv_frootc_xfer_to_litter        => pcisof%hrv_frootc_xfer_to_litter
   hrv_livestemc_xfer_to_litter     => pcisof%hrv_livestemc_xfer_to_litter
   hrv_deadstemc_xfer_to_litter     => pcisof%hrv_deadstemc_xfer_to_litter
   hrv_livecrootc_xfer_to_litter    => pcisof%hrv_livecrootc_xfer_to_litter
   hrv_deadcrootc_xfer_to_litter    => pcisof%hrv_deadcrootc_xfer_to_litter
   hrv_gresp_xfer_to_litter         => pcisof%hrv_gresp_xfer_to_litter
   leaf_prof                      => clm3%g%l%c%p%pps%leaf_prof
   froot_prof                     => clm3%g%l%c%p%pps%froot_prof
   croot_prof                     => clm3%g%l%c%p%pps%croot_prof
   stem_prof                      => clm3%g%l%c%p%pps%stem_prof

   do j = 1, nlevdecomp
      do pi = 1,maxpatch_pft
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            
            if (pi <=  npfts(c)) then
               p = pfti(c) + pi - 1
               
               if (pactive(p)) then
                  
                  ! leaf harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_cel_c(c,j) = harvest_c_to_litr_cel_c(c,j) + &
                       hrv_leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_lig_c(c,j) = harvest_c_to_litr_lig_c(c,j) + &
                       hrv_leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  
                  ! fine root harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_cel_c(c,j) = harvest_c_to_litr_cel_c(c,j) + &
                       hrv_frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_lig_c(c,j) = harvest_c_to_litr_lig_c(c,j) + &
                       hrv_frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  
                  ! wood harvest mortality carbon fluxes
                  harvest_c_to_cwdc(c,j)  = harvest_c_to_cwdc(c,j)  + &
                       hrv_livestemc_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                       hrv_livecrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_cwdc(c,j) = harvest_c_to_cwdc(c,j) + &
                       hrv_deadcrootc_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  
                  ! storage harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_leafc_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)     = harvest_c_to_litr_met_c(c,j)     + &
                       hrv_frootc_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_livestemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_deadstemc_storage_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_livecrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_deadcrootc_storage_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_gresp_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  
                  ! transfer harvest mortality carbon fluxes
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_leafc_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)     = harvest_c_to_litr_met_c(c,j)     + &
                       hrv_frootc_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_livestemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)  = harvest_c_to_litr_met_c(c,j)  + &
                       hrv_deadstemc_xfer_to_litter(p)  * wtcol(p) * stem_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_livecrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j) = harvest_c_to_litr_met_c(c,j) + &
                       hrv_deadcrootc_xfer_to_litter(p) * wtcol(p) * croot_prof(p,j)
                  harvest_c_to_litr_met_c(c,j)      = harvest_c_to_litr_met_c(c,j)      + &
                       hrv_gresp_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
               end if
            end if
            
         end do
         
      end do
   end do
   
   
   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         
         if (pi <=  npfts(c)) then
            p = pfti(c) + pi - 1
            
            if (pactive(p)) then
               
               
               chrv_deadstemc_to_prod10c(c)  = chrv_deadstemc_to_prod10c(c)  + &
                    phrv_deadstemc_to_prod10c(p)  * wtcol(p)
               chrv_deadstemc_to_prod100c(c)  = chrv_deadstemc_to_prod100c(c)  + &
                    phrv_deadstemc_to_prod100c(p)  * wtcol(p)
               
            end if
         end if
         
      end do
      
   end do
   
 end subroutine CNCIsoHarvestPftToColumn
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CIsoFluxCalc
!
! !INTERFACE, isotope:
subroutine CIsoFluxCalc(ciso_flux, ctot_flux, ciso_state, ctot_state, &
	                    num, filter, frax_c13, diag, isotope)
!
! !DESCRIPTION:
! On the radiation time step, set the carbon isotopic flux
! variables (except for gap-phase mortality and fire fluxes)
!
! !USES:
   use clmtype
!
! !ARGUMENTS:
   implicit none
   real(r8), pointer   :: ciso_flux(:)      !OUTPUT isoC flux
   real(r8), pointer   :: ctot_flux(:)      !INPUT  totC flux
   real(r8), pointer   :: ciso_state(:)     !INPUT  isoC state, upstream pool
   real(r8), pointer   :: ctot_state(:)     !INPUT  totC state, upstream pool
   real(r8), intent(in):: frax_c13          ! fractionation factor (1 = no fractionation) for C13
   integer, intent(in) :: num               ! number of filter members
   integer, intent(in) :: filter(:)         ! filter indices
   integer, intent(in) :: diag              ! 0=no diagnostics, 1=print diagnostics
   character(len=*), intent(in) :: isotope  ! 'c13' or 'c14'

!
! !CALLED FROM:
! subroutine CIsoFlux1
!
! !REVISION HISTORY:
!
! !OTHER LOCAL VARIABLES:
   integer :: i,f     ! indices
   real(r8) :: temp
   real(r8) :: frax
!

   ! if C14, double the fractionation
   select case (isotope)
   case ('c14')
      frax = 1._r8 + (1._r8 - frax_c13) * 2._r8
   case ('c13')
      frax = frax_c13
   case default
      call endrun('CNCIsoFluxMod: iso must be either c13 or c14')
   end select

   ! loop over the supplied filter
   do f = 1,num
      i = filter(f)
      if (ctot_state(i) /= 0._r8) then
      	ciso_flux(i) = ctot_flux(i) * (ciso_state(i)/ctot_state(i)) * frax
      else
      	ciso_flux(i) = 0._r8
      end if
      
      if (diag == 1) then
      ! put diagnostic print statements here for isoC flux calculations
      end if
   end do
end subroutine CIsoFluxCalc
!-----------------------------------------------------------------------





#endif

end module CNCIsoFluxMod
 
