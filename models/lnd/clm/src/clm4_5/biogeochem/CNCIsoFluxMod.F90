module CNCIsoFluxMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for carbon isotopic flux variable update, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod , only: r8 => shr_kind_r8
  use clm_varpar   , only: ndecomp_cascade_transitions, nlevdecomp, ndecomp_pools
  use abortutils   , only: endrun
  use shr_log_mod  , only: errMsg => shr_log_errMsg
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public  :: CIsoFlux1
  public  :: CIsoFlux2
  public  :: CIsoFlux2h
  public  :: CIsoFlux3
  private :: CNCIsoLitterToColumn
  private :: CNCIsoGapPftToColumn
  private :: CNCIsoHarvestPftToColumn
  private :: CIsoFluxCalc
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
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
    ! !LOCAL VARIABLES:
    type(pft_type), pointer :: p
    type(column_type), pointer :: c
    type(pft_cflux_type), pointer :: pcisof
    type(pft_cstate_type), pointer :: pcisos
    type(column_cflux_type), pointer :: ccisof
    type(column_cstate_type), pointer :: ccisos
    integer :: fp,pi,l
    integer :: fc,cc,j
    !-----------------------------------------------------------------------

   p => pft
   c => col
   select case (isotope)
   case ('c14')
      pcisof =>  pc14f
      pcisos =>  pc14s
      ccisof =>  cc14f
      ccisos =>  cc14s
   case ('c13')
      pcisof =>  pc13f
      pcisos =>  pc13s
      ccisof =>  cc13f
      ccisos =>  cc13s
   case default
      call endrun(msg='CNCIsoFluxMod: iso must be either c13 or c14'//errMsg(__FILE__, __LINE__))
   end select

   associate(&
   cascade_donor_pool  =>  decomp_cascade_con%cascade_donor_pool  & !  [integer (:)]  which pool is C taken from for a given decomposition step 
   )

   ! pft-level non-mortality fluxes
   
   ! Note: if the variables which are arguments to CIsoFluxCalc are ever changed to NOT be
   ! pointers, then the CIsoFluxCalc routine will need to be changed to declare the bounds
   ! of each argument, these bounds will need to be passed in, and - importantly for
   ! threading to work properly - the subroutine calls will need to be changed so that
   ! instead of 'call CIsoFluxCalc(foo, ...)' we have 'call CIsoFluxCalc(foo(begp:endp), ...)'.

   call CIsoFluxCalc(pcisof%leafc_xfer_to_leafc, pcf%leafc_xfer_to_leafc, &
                    pcisos%leafc_xfer, pcs%leafc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%frootc_xfer_to_frootc, pcf%frootc_xfer_to_frootc, &
                    pcisos%frootc_xfer, pcs%frootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livestemc_xfer_to_livestemc, pcf%livestemc_xfer_to_livestemc, &
                    pcisos%livestemc_xfer, pcs%livestemc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%deadstemc_xfer_to_deadstemc, pcf%deadstemc_xfer_to_deadstemc, &
                    pcisos%deadstemc_xfer, pcs%deadstemc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livecrootc_xfer_to_livecrootc, pcf%livecrootc_xfer_to_livecrootc, &
                    pcisos%livecrootc_xfer, pcs%livecrootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%deadcrootc_xfer_to_deadcrootc, pcf%deadcrootc_xfer_to_deadcrootc, &
                    pcisos%deadcrootc_xfer, pcs%deadcrootc_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%leafc_to_litter, pcf%leafc_to_litter, &
                    pcisos%leafc, pcs%leafc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%frootc_to_litter, pcf%frootc_to_litter, &
                    pcisos%frootc, pcs%frootc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livestemc_to_deadstemc, pcf%livestemc_to_deadstemc, &
                    pcisos%livestemc, pcs%livestemc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livecrootc_to_deadcrootc, pcf%livecrootc_to_deadcrootc, &
                    pcisos%livecrootc, pcs%livecrootc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%leaf_curmr, pcf%leaf_curmr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%froot_curmr, pcf%froot_curmr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livestem_curmr, pcf%livestem_curmr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livecroot_curmr, pcf%livecroot_curmr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%leaf_xsmr, pcf%leaf_xsmr, &
                    pcisos%totvegc, pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%froot_xsmr, pcf%froot_xsmr, &
                    pcisos%totvegc, pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livestem_xsmr, pcf%livestem_xsmr, &
                    pcisos%totvegc, pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livecroot_xsmr, pcf%livecroot_xsmr, &
                    pcisos%totvegc, pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_xsmrpool, pcf%cpool_to_xsmrpool, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_leafc, pcf%cpool_to_leafc, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_leafc_storage, pcf%cpool_to_leafc_storage, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_frootc, pcf%cpool_to_frootc, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_frootc_storage, pcf%cpool_to_frootc_storage, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_livestemc, pcf%cpool_to_livestemc, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_livestemc_storage, pcf%cpool_to_livestemc_storage, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_deadstemc, pcf%cpool_to_deadstemc, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_deadstemc_storage, pcf%cpool_to_deadstemc_storage, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_livecrootc, pcf%cpool_to_livecrootc, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_livecrootc_storage, pcf%cpool_to_livecrootc_storage, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_deadcrootc, pcf%cpool_to_deadcrootc, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_deadcrootc_storage, pcf%cpool_to_deadcrootc_storage, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_leaf_gr, pcf%cpool_leaf_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_froot_gr, pcf%cpool_froot_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_livestem_gr, pcf%cpool_livestem_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_deadstem_gr, pcf%cpool_deadstem_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_livecroot_gr, pcf%cpool_livecroot_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_deadcroot_gr, pcf%cpool_deadcroot_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_leaf_storage_gr, pcf%cpool_leaf_storage_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_froot_storage_gr, pcf%cpool_froot_storage_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_livestem_storage_gr, pcf%cpool_livestem_storage_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_deadstem_storage_gr, pcf%cpool_deadstem_storage_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_livecroot_storage_gr, pcf%cpool_livecroot_storage_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_deadcroot_storage_gr, pcf%cpool_deadcroot_storage_gr, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%cpool_to_gresp_storage, pcf%cpool_to_gresp_storage, &
                    pcisos%cpool, pcs%cpool, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_leaf_gr, pcf%transfer_leaf_gr, &
                    pcisos%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_froot_gr, pcf%transfer_froot_gr, &
                    pcisos%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_livestem_gr, pcf%transfer_livestem_gr, &
                    pcisos%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_deadstem_gr, pcf%transfer_deadstem_gr, &
                    pcisos%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_livecroot_gr, pcf%transfer_livecroot_gr, &
                    pcisos%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%transfer_deadcroot_gr, pcf%transfer_deadcroot_gr, &
                    pcisos%gresp_xfer, pcs%gresp_xfer, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%leafc_storage_to_xfer, pcf%leafc_storage_to_xfer, &
                    pcisos%leafc_storage, pcs%leafc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%frootc_storage_to_xfer, pcf%frootc_storage_to_xfer, &
                    pcisos%frootc_storage, pcs%frootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livestemc_storage_to_xfer, pcf%livestemc_storage_to_xfer, &
                    pcisos%livestemc_storage, pcs%livestemc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%deadstemc_storage_to_xfer, pcf%deadstemc_storage_to_xfer, &
                    pcisos%deadstemc_storage, pcs%deadstemc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%livecrootc_storage_to_xfer, pcf%livecrootc_storage_to_xfer, &
                    pcisos%livecrootc_storage, pcs%livecrootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%deadcrootc_storage_to_xfer, pcf%deadcrootc_storage_to_xfer, &
                    pcisos%deadcrootc_storage, pcs%deadcrootc_storage, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
	call CIsoFluxCalc(pcisof%gresp_storage_to_xfer, pcf%gresp_storage_to_xfer, &
                    pcisos%gresp_storage, pcs%gresp_storage, &
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
            if ( ccs%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) /= 0._r8) then
               ccisof%decomp_cascade_hr_vr(cc,j,l)  =  ccf%decomp_cascade_hr_vr(cc,j,l) * &
                    (ccisos%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) &
                       / ccs%decomp_cpools_vr(cc,j,cascade_donor_pool(l))) * 1._r8
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
            if ( ccs%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) /= 0._r8) then
               ccisof%decomp_cascade_ctransfer_vr(cc,j,l)  =  ccf%decomp_cascade_ctransfer_vr(cc,j,l) * &
                    (ccisos%decomp_cpools_vr(cc,j,cascade_donor_pool(l)) &
                       / ccs%decomp_cpools_vr(cc,j,cascade_donor_pool(l))) * 1._r8
            else
               ccisof%decomp_cascade_ctransfer_vr(cc,j,l) = 0._r8
            end if
         end do
      end do
   end do


 end associate
end subroutine CIsoFlux1

!-----------------------------------------------------------------------
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
  ! !LOCAL VARIABLES:
  type(pft_type), pointer :: p
  type(pft_cflux_type), pointer :: pcisof
  type(pft_cstate_type), pointer :: pcisos
  integer :: fp,pi
  !-----------------------------------------------------------------------

   p => pft
   select case (isotope)
   case ('c14')
      pcisof =>  pc14f
      pcisos =>  pc14s
   case ('c13')
      pcisof =>  pc13f
      pcisos =>  pc13s
   case default
      call endrun(msg='CNCIsoFluxMod: iso must be either c13 or c14'//errMsg(__FILE__, __LINE__))
   end select

   ! pft-level gap mortality fluxes
   
   call CIsoFluxCalc(pcisof%m_leafc_to_litter, pcf%m_leafc_to_litter, &
                     pcisos%leafc, pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_leafc_storage_to_litter, pcf%m_leafc_storage_to_litter, &
                     pcisos%leafc_storage, pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_leafc_xfer_to_litter, pcf%m_leafc_xfer_to_litter, &
                     pcisos%leafc_xfer, pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_to_litter, pcf%m_frootc_to_litter, &
                     pcisos%frootc, pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_storage_to_litter, pcf%m_frootc_storage_to_litter, &
                     pcisos%frootc_storage, pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_xfer_to_litter, pcf%m_frootc_xfer_to_litter, &
                     pcisos%frootc_xfer, pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_to_litter, pcf%m_livestemc_to_litter, &
                     pcisos%livestemc, pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_storage_to_litter, pcf%m_livestemc_storage_to_litter, &
                     pcisos%livestemc_storage, pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_xfer_to_litter, pcf%m_livestemc_xfer_to_litter, &
                     pcisos%livestemc_xfer, pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_to_litter, pcf%m_deadstemc_to_litter, &
                     pcisos%deadstemc, pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_storage_to_litter, pcf%m_deadstemc_storage_to_litter, &
                     pcisos%deadstemc_storage, pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_xfer_to_litter, pcf%m_deadstemc_xfer_to_litter, &
                     pcisos%deadstemc_xfer, pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_to_litter, pcf%m_livecrootc_to_litter, &
                     pcisos%livecrootc, pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_storage_to_litter, pcf%m_livecrootc_storage_to_litter, &
                     pcisos%livecrootc_storage, pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_xfer_to_litter, pcf%m_livecrootc_xfer_to_litter, &
                     pcisos%livecrootc_xfer, pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_to_litter, pcf%m_deadcrootc_to_litter, &
                     pcisos%deadcrootc, pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_storage_to_litter, pcf%m_deadcrootc_storage_to_litter, &
                     pcisos%deadcrootc_storage, pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_xfer_to_litter, pcf%m_deadcrootc_xfer_to_litter, &
                     pcisos%deadcrootc_xfer, pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_gresp_storage_to_litter, pcf%m_gresp_storage_to_litter, &
                     pcisos%gresp_storage, pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_gresp_xfer_to_litter, pcf%m_gresp_xfer_to_litter, &
                     pcisos%gresp_xfer, pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
	! call routine to shift pft-level gap mortality fluxes to column, for isotopes
   ! the non-isotope version of this routine is in CNGapMortalityMod.F90.

	call CNCIsoGapPftToColumn(num_soilc, filter_soilc, isotope)
   
end subroutine CIsoFlux2

!-----------------------------------------------------------------------
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
  ! !LOCAL VARIABLES:
  type(pft_type), pointer :: p
  type(pft_cflux_type), pointer :: pcisof
  type(pft_cstate_type), pointer :: pcisos
  integer :: fp,pi
  !-----------------------------------------------------------------------
  ! set local pointers
   p => pft
   select case (isotope)
   case ('c14')
      pcisof =>  pc14f
      pcisos =>  pc14s
   case ('c13')
      pcisof =>  pc13f
      pcisos =>  pc13s
   case default
      call endrun(msg='CNCIsoFluxMod: iso must be either c13 or c14'//errMsg(__FILE__, __LINE__))
   end select
	
   ! pft-level gap mortality fluxes
   
   call CIsoFluxCalc(pcisof%hrv_leafc_to_litter, pcf%hrv_leafc_to_litter, &
                     pcisos%leafc, pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_leafc_storage_to_litter, pcf%hrv_leafc_storage_to_litter, &
                     pcisos%leafc_storage, pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_leafc_xfer_to_litter, pcf%hrv_leafc_xfer_to_litter, &
                     pcisos%leafc_xfer, pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_frootc_to_litter, pcf%hrv_frootc_to_litter, &
                     pcisos%frootc, pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_frootc_storage_to_litter, pcf%hrv_frootc_storage_to_litter, &
                     pcisos%frootc_storage, pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_frootc_xfer_to_litter, pcf%hrv_frootc_xfer_to_litter, &
                     pcisos%frootc_xfer, pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livestemc_to_litter, pcf%hrv_livestemc_to_litter, &
                     pcisos%livestemc, pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livestemc_storage_to_litter, pcf%hrv_livestemc_storage_to_litter, &
                     pcisos%livestemc_storage, pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livestemc_xfer_to_litter, pcf%hrv_livestemc_xfer_to_litter, &
                     pcisos%livestemc_xfer, pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadstemc_to_prod10c, pcf%hrv_deadstemc_to_prod10c, &
                     pcisos%deadstemc, pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadstemc_to_prod100c, pcf%hrv_deadstemc_to_prod100c, &
                     pcisos%deadstemc, pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadstemc_storage_to_litter, pcf%hrv_deadstemc_storage_to_litter, &
                     pcisos%deadstemc_storage, pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadstemc_xfer_to_litter, pcf%hrv_deadstemc_xfer_to_litter, &
                     pcisos%deadstemc_xfer, pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livecrootc_to_litter, pcf%hrv_livecrootc_to_litter, &
                     pcisos%livecrootc, pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livecrootc_storage_to_litter, pcf%hrv_livecrootc_storage_to_litter, &
                     pcisos%livecrootc_storage, pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_livecrootc_xfer_to_litter, pcf%hrv_livecrootc_xfer_to_litter, &
                     pcisos%livecrootc_xfer, pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadcrootc_to_litter, pcf%hrv_deadcrootc_to_litter, &
                     pcisos%deadcrootc, pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadcrootc_storage_to_litter, pcf%hrv_deadcrootc_storage_to_litter, &
                     pcisos%deadcrootc_storage, pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_deadcrootc_xfer_to_litter, pcf%hrv_deadcrootc_xfer_to_litter, &
                     pcisos%deadcrootc_xfer, pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_gresp_storage_to_litter, pcf%hrv_gresp_storage_to_litter, &
                     pcisos%gresp_storage, pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_gresp_xfer_to_litter, pcf%hrv_gresp_xfer_to_litter, &
                     pcisos%gresp_xfer, pcs%gresp_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%hrv_xsmrpool_to_atm, pcf%hrv_xsmrpool_to_atm, &
                    pcisos%totvegc, pcs%totvegc, &
                    num_soilp, filter_soilp, 1._r8, 0, isotope)
                    
   ! call routine to shift pft-level gap mortality fluxes to column, for isotopes
   ! the non-isotope version of this routine is in CNGapMortalityMod.F90.

   call CNCIsoHarvestPftToColumn(num_soilc, filter_soilc, isotope)
   
end subroutine CIsoFlux2h

!-----------------------------------------------------------------------
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
  ! !LOCAL VARIABLES:
  type(pft_type), pointer :: p
  type(column_type), pointer :: c
  type(pft_cflux_type), pointer :: pcisof
  type(pft_cstate_type), pointer :: pcisos
  type(column_cflux_type), pointer :: ccisof
  type(column_cstate_type), pointer :: ccisos
  integer :: fp,pi,l,pp
  integer :: fc,cc,j
  !-----------------------------------------------------------------------

   p => pft
   c => col
   select case (isotope)
   case ('c14')
      pcisof =>  pc14f
      pcisos =>  pc14s
      ccisof =>  cc14f
      ccisos =>  cc14s
   case ('c13')
      pcisof =>  pc13f
      pcisos =>  pc13s
      ccisof =>  cc13f
      ccisos =>  cc13s
   case default
      call endrun(msg='CNCIsoFluxMod: iso must be either c13 or c14'//errMsg(__FILE__, __LINE__))
   end select

   associate(&
   croot_prof =>   pps%croot_prof , & !  [real(r8) (:,:)]  (1/m) profile of coarse roots                          
   stem_prof  =>   pps%stem_prof  , & !  [real(r8) (:,:)]  (1/m) profile of stems                                 
   npfts      =>   col%npfts      , & !  [integer (:)]  number of pfts for each column                            
   pfti       =>   col%pfti       , & !  [integer (:)]  beginning pft index for each column                       
   wtcol      =>   pft%wtcol      , & !  [real(r8) (:)]  weight (relative to column) for this pft (0-1)           
   pactive    =>   pft%active       & !  [logical (:)]  true=>do computations on this pft 
   )
	
   ! pft-level fire mortality fluxes
   
   call CIsoFluxCalc(pcisof%m_leafc_to_fire, pcf%m_leafc_to_fire, &
                     pcisos%leafc, pcs%leafc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_leafc_storage_to_fire, pcf%m_leafc_storage_to_fire, &
                     pcisos%leafc_storage, pcs%leafc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_leafc_xfer_to_fire, pcf%m_leafc_xfer_to_fire, &
                     pcisos%leafc_xfer, pcs%leafc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_to_fire, pcf%m_frootc_to_fire, &
                     pcisos%frootc, pcs%frootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_storage_to_fire, pcf%m_frootc_storage_to_fire, &
                     pcisos%frootc_storage, pcs%frootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_frootc_xfer_to_fire, pcf%m_frootc_xfer_to_fire, &
                     pcisos%frootc_xfer, pcs%frootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_to_fire, pcf%m_livestemc_to_fire, &
                     pcisos%livestemc, pcs%livestemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_storage_to_fire, pcf%m_livestemc_storage_to_fire, &
                     pcisos%livestemc_storage, pcs%livestemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livestemc_xfer_to_fire, pcf%m_livestemc_xfer_to_fire, &
                     pcisos%livestemc_xfer, pcs%livestemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_to_fire, pcf%m_deadstemc_to_fire, &
                     pcisos%deadstemc, pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_to_litter_fire, pcf%m_deadstemc_to_litter_fire, &
                     pcisos%deadstemc, pcs%deadstemc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_storage_to_fire, pcf%m_deadstemc_storage_to_fire, &
                     pcisos%deadstemc_storage, pcs%deadstemc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadstemc_xfer_to_fire, pcf%m_deadstemc_xfer_to_fire, &
                     pcisos%deadstemc_xfer, pcs%deadstemc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_to_fire, pcf%m_livecrootc_to_fire, &
                     pcisos%livecrootc, pcs%livecrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_storage_to_fire, pcf%m_livecrootc_storage_to_fire, &
                     pcisos%livecrootc_storage, pcs%livecrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_livecrootc_xfer_to_fire, pcf%m_livecrootc_xfer_to_fire, &
                     pcisos%livecrootc_xfer, pcs%livecrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_to_fire, pcf%m_deadcrootc_to_fire, &
                     pcisos%deadcrootc, pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_to_litter_fire, pcf%m_deadcrootc_to_litter_fire, &
                     pcisos%deadcrootc, pcs%deadcrootc, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_storage_to_fire, pcf%m_deadcrootc_storage_to_fire, &
                     pcisos%deadcrootc_storage, pcs%deadcrootc_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_deadcrootc_xfer_to_fire, pcf%m_deadcrootc_xfer_to_fire, &
                     pcisos%deadcrootc_xfer, pcs%deadcrootc_xfer, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_gresp_storage_to_fire, pcf%m_gresp_storage_to_fire, &
                     pcisos%gresp_storage, pcs%gresp_storage, &
                     num_soilp, filter_soilp, 1._r8, 0, isotope)
   
   call CIsoFluxCalc(pcisof%m_gresp_xfer_to_fire, pcf%m_gresp_xfer_to_fire, &
                     pcisos%gresp_xfer, pcs%gresp_xfer, &
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
            if ( ccs%decomp_cpools_vr(cc,j,l) /= 0._r8) then
               ccisof%m_decomp_cpools_to_fire_vr(cc,j,l)  =  ccf%m_decomp_cpools_to_fire_vr(cc,j,l) * &
                    (ccisos%decomp_cpools_vr(cc,j,l) / ccs%decomp_cpools_vr(cc,j,l)) * 1._r8
            else
               ccisof%m_decomp_cpools_to_fire_vr(cc,j,l) = 0._r8
            end if
         end do
      end do
   end do

 end associate
end subroutine CIsoFlux3

!-----------------------------------------------------------------------
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
  ! !LOCAL VARIABLES:
  type(pft_cflux_type), pointer :: pcisof
  type(column_cflux_type), pointer :: ccisof
  integer :: fc,c,pi,p,j
  !-----------------------------------------------------------------------
    ! select which isotope
    select case (isotope)
    case ('c14')
       pcisof =>  pc14f
       ccisof =>  cc14f
    case ('c13')
       pcisof =>  pc13f
       ccisof =>  cc13f
    case default
       call endrun(msg='CNCIsoFluxMod: iso must be either c13 or c14'//errMsg(__FILE__, __LINE__))
    end select

   associate(& 
   ivt                       =>    pft%itype                        , & ! Input:  [integer (:)]  pft vegetation type                                
   wtcol                     =>    pft%wtcol                        , & ! Input:  [real(r8) (:)]  weight (relative to column) for this pft (0-1)    
   pactive                   =>    pft%active                       , & ! Input:  [logical (:)]  true=>do computations on this pft
   leafc_to_litter           =>    pcisof%leafc_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   frootc_to_litter          =>    pcisof%frootc_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   npfts                     =>    col%npfts                        , & ! Input:  [integer (:)]  number of pfts for each column                     
   pfti                      =>    col%pfti                         , & ! Input:  [integer (:)]  beginning pft index for each column                
   lf_flab                   =>    pftcon%lf_flab                   , & ! Input:  [real(r8) (:)]  leaf litter labile fraction                       
   lf_fcel                   =>    pftcon%lf_fcel                   , & ! Input:  [real(r8) (:)]  leaf litter cellulose fraction                    
   lf_flig                   =>    pftcon%lf_flig                   , & ! Input:  [real(r8) (:)]  leaf litter lignin fraction                       
   fr_flab                   =>    pftcon%fr_flab                   , & ! Input:  [real(r8) (:)]  fine root litter labile fraction                  
   fr_fcel                   =>    pftcon%fr_fcel                   , & ! Input:  [real(r8) (:)]  fine root litter cellulose fraction               
   fr_flig                   =>    pftcon%fr_flig                   , & ! Input:  [real(r8) (:)]  fine root litter lignin fraction                  
   phenology_c_to_litr_met_c =>    ccisof%phenology_c_to_litr_met_c , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
   phenology_c_to_litr_cel_c =>    ccisof%phenology_c_to_litr_cel_c , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
   phenology_c_to_litr_lig_c =>    ccisof%phenology_c_to_litr_lig_c , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)
   leaf_prof                 =>    pps%leaf_prof                    , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of leaves                         
   froot_prof                =>    pps%froot_prof                     & ! InOut:  [real(r8) (:,:)]  (1/m) profile of fine roots                     
   )

    do j = 1, nlevdecomp
       do pi = 1,max_pft_per_col
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             
             if ( pi <=  npfts(c) ) then
                p = pfti(c) + pi - 1
                if (pactive(p)) then
                   ! leaf litter carbon fluxes
                   phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                      + leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                      + leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                      + leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                   
                   ! fine root litter carbon fluxes
                   phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                      + frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                   phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                      + frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                   phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                      + frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)
                end if
             end if
             
          end do
       end do
       
    end do
    
    end associate 
   end subroutine CNCIsoLitterToColumn

   !-----------------------------------------------------------------------
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
     ! !LOCAL VARIABLES:
     type(pft_cflux_type), pointer :: pcisof
     type(column_cflux_type), pointer :: ccisof
     integer :: fc,c,pi,p,j               ! indices
     !-----------------------------------------------------------------------

    ! select which isotope
    select case (isotope)
    case ('c14')
       pcisof =>  pc14f
       ccisof =>  cc14f
    case ('c13')
       pcisof =>  pc13f
       ccisof =>  cc13f
    case default
       call endrun(msg='CNCIsoFluxMod: iso must be either c13 or c14'//errMsg(__FILE__, __LINE__))
    end select

   associate(& 
   lf_flab                             =>    pftcon%lf_flab                              , & ! Input:  [real(r8) (:)]  leaf litter labile fraction                       
   lf_fcel                             =>    pftcon%lf_fcel                              , & ! Input:  [real(r8) (:)]  leaf litter cellulose fraction                    
   lf_flig                             =>    pftcon%lf_flig                              , & ! Input:  [real(r8) (:)]  leaf litter lignin fraction                       
   fr_flab                             =>    pftcon%fr_flab                              , & ! Input:  [real(r8) (:)]  fine root litter labile fraction                  
   fr_fcel                             =>    pftcon%fr_fcel                              , & ! Input:  [real(r8) (:)]  fine root litter cellulose fraction               
   fr_flig                             =>    pftcon%fr_flig                              , & ! Input:  [real(r8) (:)]  fine root litter lignin fraction                  
   npfts                               =>   col%npfts                                    , & ! Input:  [integer (:)]  number of pfts for each column                     
   pfti                                =>   col%pfti                                     , & ! Input:  [integer (:)]  beginning pft index for each column                
   gap_mortality_c_to_litr_met_c       =>    ccisof%gap_mortality_c_to_litr_met_c        , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
   gap_mortality_c_to_litr_cel_c       =>    ccisof%gap_mortality_c_to_litr_cel_c        , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
   gap_mortality_c_to_litr_lig_c       =>    ccisof%gap_mortality_c_to_litr_lig_c        , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
   gap_mortality_c_to_cwdc             =>    ccisof%gap_mortality_c_to_cwdc              , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with gap mortality to CWD pool (gC/m3/s)
   ivt                                 =>   pft%itype                                    , & ! Input:  [integer (:)]  pft vegetation type                                
   wtcol                               =>   pft%wtcol                                    , & ! Input:  [real(r8) (:)]  pft weight relative to column (0-1)               
   pactive                             =>    pft%active                                  , & ! Input:  [logical (:)]  true=>do computations on this pft 
   m_leafc_to_litter                   =>    pcisof%m_leafc_to_litter                    , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_to_litter                  =>    pcisof%m_frootc_to_litter                   , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_to_litter               =>    pcisof%m_livestemc_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_to_litter               =>    pcisof%m_deadstemc_to_litter                , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_to_litter              =>    pcisof%m_livecrootc_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_to_litter              =>    pcisof%m_deadcrootc_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_storage_to_litter           =>    pcisof%m_leafc_storage_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_storage_to_litter          =>    pcisof%m_frootc_storage_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_storage_to_litter       =>    pcisof%m_livestemc_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_storage_to_litter       =>    pcisof%m_deadstemc_storage_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_storage_to_litter      =>    pcisof%m_livecrootc_storage_to_litter       , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_storage_to_litter      =>    pcisof%m_deadcrootc_storage_to_litter       , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_storage_to_litter           =>    pcisof%m_gresp_storage_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   m_leafc_xfer_to_litter              =>    pcisof%m_leafc_xfer_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   m_frootc_xfer_to_litter             =>    pcisof%m_frootc_xfer_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   m_livestemc_xfer_to_litter          =>    pcisof%m_livestemc_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_deadstemc_xfer_to_litter          =>    pcisof%m_deadstemc_xfer_to_litter           , & ! Input:  [real(r8) (:)]                                                    
   m_livecrootc_xfer_to_litter         =>    pcisof%m_livecrootc_xfer_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_deadcrootc_xfer_to_litter         =>    pcisof%m_deadcrootc_xfer_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   m_gresp_xfer_to_litter              =>    pcisof%m_gresp_xfer_to_litter               , & ! Input:  [real(r8) (:)]                                                    
   leaf_prof                           =>    pps%leaf_prof                               , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of leaves                         
   froot_prof                          =>    pps%froot_prof                              , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of fine roots                     
   croot_prof                          =>    pps%croot_prof                              , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of coarse roots                   
   stem_prof                           =>    pps%stem_prof                                 & ! InOut:  [real(r8) (:,:)]  (1/m) profile of stems                          
   )

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

    end associate 
 end subroutine CNCIsoGapPftToColumn

 !-----------------------------------------------------------------------
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
   ! !LOCAL VARIABLES:
   type(pft_cflux_type), pointer :: pcisof
   type(column_cflux_type), pointer :: ccisof
   integer :: fc,c,pi,p,j               ! indices
   !-----------------------------------------------------------------------
    ! select which isotope
    select case (isotope)
    case ('c14')
       pcisof =>  pc14f
       ccisof =>  cc14f
    case ('c13')
       pcisof =>  pc13f
       ccisof =>  cc13f
    case default
       call endrun(msg='CNCIsoFluxMod: iso must be either c13 or c14'//errMsg(__FILE__, __LINE__))
    end select

   associate(& 
   lf_flab                             =>    pftcon%lf_flab                              , & ! Input:  [real(r8) (:)]  leaf litter labile fraction                       
   lf_fcel                             =>    pftcon%lf_fcel                              , & ! Input:  [real(r8) (:)]  leaf litter cellulose fraction                    
   lf_flig                             =>    pftcon%lf_flig                              , & ! Input:  [real(r8) (:)]  leaf litter lignin fraction                       
   fr_flab                             =>    pftcon%fr_flab                              , & ! Input:  [real(r8) (:)]  fine root litter labile fraction                  
   fr_fcel                             =>    pftcon%fr_fcel                              , & ! Input:  [real(r8) (:)]  fine root litter cellulose fraction               
   fr_flig                             =>    pftcon%fr_flig                              , & ! Input:  [real(r8) (:)]  fine root litter lignin fraction                  
   npfts                               =>   col%npfts                                    , & ! Input:  [integer (:)]  number of pfts for each column                     
   pfti                                =>   col%pfti                                     , & ! Input:  [integer (:)]  beginning pft index for each column                
   chrv_deadstemc_to_prod10c           =>    ccisof%hrv_deadstemc_to_prod10c             , & ! InOut:  [real(r8) (:)]                                                    
   chrv_deadstemc_to_prod100c          =>    ccisof%hrv_deadstemc_to_prod100c            , & ! InOut:  [real(r8) (:)]                                                    
   harvest_c_to_litr_met_c             =>    ccisof%harvest_c_to_litr_met_c              , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with harvest to litter metabolic pool (gC/m3/s)
   harvest_c_to_litr_cel_c             =>    ccisof%harvest_c_to_litr_cel_c              , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with harvest to litter cellulose pool (gC/m3/s)
   harvest_c_to_litr_lig_c             =>    ccisof%harvest_c_to_litr_lig_c              , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with harvest to litter lignin pool (gC/m3/s)
   harvest_c_to_cwdc                   =>    ccisof%harvest_c_to_cwdc                    , & ! InOut:  [real(r8) (:,:)]  C fluxes associated with harvest to CWD pool (gC/m3/s)
   ivt                                 =>   pft%itype                                    , & ! Input:  [integer (:)]  pft vegetation type                                
   wtcol                               =>   pft%wtcol                                    , & ! Input:  [real(r8) (:)]  pft weight relative to column (0-1)               
   pactive                             =>    pft%active                                  , & ! Input:  [logical (:)]  true=>do computations on this pft 
   hrv_leafc_to_litter                 =>    pcisof%hrv_leafc_to_litter                  , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_to_litter                =>    pcisof%hrv_frootc_to_litter                 , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_to_litter             =>    pcisof%hrv_livestemc_to_litter              , & ! Input:  [real(r8) (:)]                                                    
   phrv_deadstemc_to_prod10c           =>    pcisof%hrv_deadstemc_to_prod10c             , & ! Input:  [real(r8) (:)]                                                    
   phrv_deadstemc_to_prod100c          =>    pcisof%hrv_deadstemc_to_prod100c            , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_to_litter            =>    pcisof%hrv_livecrootc_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_to_litter            =>    pcisof%hrv_deadcrootc_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafc_storage_to_litter         =>    pcisof%hrv_leafc_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_storage_to_litter        =>    pcisof%hrv_frootc_storage_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_storage_to_litter     =>    pcisof%hrv_livestemc_storage_to_litter      , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_storage_to_litter     =>    pcisof%hrv_deadstemc_storage_to_litter      , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_storage_to_litter    =>    pcisof%hrv_livecrootc_storage_to_litter     , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_storage_to_litter    =>    pcisof%hrv_deadcrootc_storage_to_litter     , & ! Input:  [real(r8) (:)]                                                    
   hrv_gresp_storage_to_litter         =>    pcisof%hrv_gresp_storage_to_litter          , & ! Input:  [real(r8) (:)]                                                    
   hrv_leafc_xfer_to_litter            =>    pcisof%hrv_leafc_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   hrv_frootc_xfer_to_litter           =>    pcisof%hrv_frootc_xfer_to_litter            , & ! Input:  [real(r8) (:)]                                                    
   hrv_livestemc_xfer_to_litter        =>    pcisof%hrv_livestemc_xfer_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadstemc_xfer_to_litter        =>    pcisof%hrv_deadstemc_xfer_to_litter         , & ! Input:  [real(r8) (:)]                                                    
   hrv_livecrootc_xfer_to_litter       =>    pcisof%hrv_livecrootc_xfer_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_deadcrootc_xfer_to_litter       =>    pcisof%hrv_deadcrootc_xfer_to_litter        , & ! Input:  [real(r8) (:)]                                                    
   hrv_gresp_xfer_to_litter            =>    pcisof%hrv_gresp_xfer_to_litter             , & ! Input:  [real(r8) (:)]                                                    
   leaf_prof                           =>    pps%leaf_prof                               , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of leaves                         
   froot_prof                          =>    pps%froot_prof                              , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of fine roots                     
   croot_prof                          =>    pps%croot_prof                              , & ! InOut:  [real(r8) (:,:)]  (1/m) profile of coarse roots                   
   stem_prof                           =>    pps%stem_prof                                 & ! InOut:  [real(r8) (:,:)]  (1/m) profile of stems                          
   )

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
   
    end associate 
  end subroutine CNCIsoHarvestPftToColumn

  !-----------------------------------------------------------------------
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
    real(r8), intent(inout), pointer :: ciso_flux(:)      ! isoC flux
    real(r8), intent(in)   , pointer :: ctot_flux(:)      ! totC flux
    real(r8), intent(in)   , pointer :: ciso_state(:)     ! isoC state, upstream pool
    real(r8), intent(in)   , pointer :: ctot_state(:)     ! totC state, upstream pool
    real(r8), intent(in) :: frax_c13          ! fractionation factor (1 = no fractionation) for C13
    integer , intent(in) :: num               ! number of filter members
    integer , intent(in) :: filter(:)         ! filter indices
    integer , intent(in) :: diag              ! 0=no diagnostics, 1=print diagnostics
    character(len=*), intent(in) :: isotope  ! 'c13' or 'c14'
    !
    ! ! LOCAL VARIABLES:
    integer :: i,f     ! indices
    real(r8) :: temp
    real(r8) :: frax
    !-----------------------------------------------------------------------

   ! if C14, double the fractionation
   select case (isotope)
   case ('c14')
      frax = 1._r8 + (1._r8 - frax_c13) * 2._r8
   case ('c13')
      frax = frax_c13
   case default
      call endrun(msg='CNCIsoFluxMod: iso must be either c13 or c14'//errMsg(__FILE__, __LINE__))
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

end module CNCIsoFluxMod
