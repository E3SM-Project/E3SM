module CNSetValueMod
  !-----------------------------------------------------------------------
  ! !MODULE: CNSetValueMod
  !
  ! !DESCRIPTION:
  ! contains code to set all CN variables to specified value
  ! Used for both initialization of special landunit values, and
  ! setting fluxes to 0.0 at the beginning of each time step
  ! 3/23/09, Peter Thornton: Added new subroutine, CNZeroFluxes_dwt(), 
  !     which initialize flux variables used in the pftdyn
  !     routines. This is called from clm_driver1, as
  !     these variables need to be initialized outside of the clumps loop.
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varpar  , only: nlevgrnd, nlevdecomp_full, ndecomp_pools, &
                          ndecomp_cascade_transitions, nlevdecomp, &
                          crop_prog
  use clm_varctl  , only: iulog, use_c13, use_c14, use_cn, use_cndv, &
                          use_nitrif_denitrif   
  use clmtype
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNZeroFluxes
  public :: CNZeroFluxes_dwt
  public :: CNSetPps
  public :: CNSetPepv
  public :: CNSetPcs
  public :: CNSetPns
  public :: CNSetPcf
  public :: CNSetPnf
  public :: CNSetCps
  public :: CNSetCcs
  public :: CNSetCns
  public :: CNSetCcf
  public :: CNSetCnf
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine CNZeroFluxes(num_filterc, filterc, num_filterp, filterp)
    !
    ! !ARGUMENTS:
    implicit none
    integer, intent(in) :: num_filterc ! number of good values in filterc
    integer, intent(in) :: filterc(:)  ! column filter
    integer, intent(in) :: num_filterp ! number of good values in filterp
    integer, intent(in) :: filterp(:)  ! pft filter
    !-----------------------------------------------------------------------

    ! zero the column-level C and N fluxes
    call CNSetCcf(num_filterc, filterc, 0._r8, ccf)

    if ( use_c13 ) call CNSetCcf(num_filterc, filterc, 0._r8, cc13f)

    if ( use_c14 ) call CNSetCcf(num_filterc, filterc, 0._r8, cc14f)

    call CNSetCnf(num_filterc, filterc, 0._r8, cnf)

    ! zero the column-average pft-level C and N fluxes
    call CNSetPcf(num_filterc, filterc, 0._r8, pcf_a)
    call CNSetPnf(num_filterc, filterc, 0._r8, pnf_a)

    ! zero the pft-level C and N fluxes
    call CNSetPcf(num_filterp, filterp, 0._r8, pcf)

    if ( use_c13 ) call CNSetPcf(num_filterp, filterp, 0._r8, pc13f)

    if ( use_c14 ) call CNSetPcf(num_filterp, filterp, 0._r8, pc14f)

    call CNSetPnf(num_filterp, filterp, 0._r8, pnf)

  end subroutine CNZeroFluxes

  !-----------------------------------------------------------------------
  subroutine CNZeroFluxes_dwt( bounds)
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use decompMod   , only : bounds_type
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds  ! bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: c, p, j          ! indices
    !-----------------------------------------------------------------------

    ! set column-level conversion and product pool fluxes
    ! to 0 at the beginning of every timestep

    do c = bounds%begc,bounds%endc
       ! C fluxes
       ccf%dwt_seedc_to_leaf(c) = 0._r8
       ccf%dwt_seedc_to_deadstem(c) = 0._r8
       ccf%dwt_conv_cflux(c) = 0._r8
       ccf%lf_conv_cflux(c) = 0._r8
       ccf%dwt_prod10c_gain(c) = 0._r8
       ccf%dwt_prod100c_gain(c) = 0._r8

       ! N fluxes
       cnf%dwt_seedn_to_leaf(c) = 0._r8
       cnf%dwt_seedn_to_deadstem(c) = 0._r8
       cnf%dwt_conv_nflux(c) = 0._r8
       cnf%dwt_prod10n_gain(c) = 0._r8
       cnf%dwt_prod100n_gain(c) = 0._r8
    end do
    if ( use_c13 ) then
       do c = bounds%begc,bounds%endc
          cc13f%dwt_seedc_to_leaf(c) = 0._r8
          cc13f%dwt_seedc_to_deadstem(c) = 0._r8
          cc13f%dwt_conv_cflux(c) = 0._r8
          cc13f%dwt_prod10c_gain(c) = 0._r8
          cc13f%dwt_prod100c_gain(c) = 0._r8
       end do
    endif

    if ( use_c14 ) then
       do c = bounds%begc,bounds%endc
          cc14f%dwt_seedc_to_leaf(c) = 0._r8
          cc14f%dwt_seedc_to_deadstem(c) = 0._r8
          cc14f%dwt_conv_cflux(c) = 0._r8
          cc14f%dwt_prod10c_gain(c) = 0._r8
          cc14f%dwt_prod100c_gain(c) = 0._r8
       end do
    endif

    do j = 1, nlevdecomp_full
       do c = bounds%begc,bounds%endc
          ! C fluxes
          ccf%dwt_frootc_to_litr_met_c(c,j) = 0._r8
          ccf%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
          ccf%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
          ccf%dwt_livecrootc_to_cwdc(c,j) = 0._r8
          ccf%dwt_deadcrootc_to_cwdc(c,j) = 0._r8

          ! N fluxes
          cnf%dwt_frootn_to_litr_met_n(c,j) = 0._r8
          cnf%dwt_frootn_to_litr_cel_n(c,j) = 0._r8
          cnf%dwt_frootn_to_litr_lig_n(c,j) = 0._r8
          cnf%dwt_livecrootn_to_cwdn(c,j) = 0._r8
          cnf%dwt_deadcrootn_to_cwdn(c,j) = 0._r8
       end do
    end do
    if ( use_c13 ) then
       do j = 1, nlevdecomp_full
          do c = bounds%begc,bounds%endc
             cc13f%dwt_frootc_to_litr_met_c(c,j) = 0._r8
             cc13f%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
             cc13f%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
             cc13f%dwt_livecrootc_to_cwdc(c,j) = 0._r8
             cc13f%dwt_deadcrootc_to_cwdc(c,j) = 0._r8
          end do
       end do
    endif
    if ( use_c14 ) then
       do j = 1, nlevdecomp_full
          do c = bounds%begc,bounds%endc
             cc14f%dwt_frootc_to_litr_met_c(c,j) = 0._r8
             cc14f%dwt_frootc_to_litr_cel_c(c,j) = 0._r8
             cc14f%dwt_frootc_to_litr_lig_c(c,j) = 0._r8
             cc14f%dwt_livecrootc_to_cwdc(c,j) = 0._r8
             cc14f%dwt_deadcrootc_to_cwdc(c,j) = 0._r8
          end do
       end do
    endif

    if (use_cn) then
       do p = bounds%begp,bounds%endp
          pcs%dispvegc(p)   = 0._r8
          pcs%storvegc(p)   = 0._r8
          pcs%totpftc(p)    = 0._r8

          pns%dispvegn(p)   = 0._r8
          pns%storvegn(p)   = 0._r8
          pns%totvegn(p)    = 0._r8
          pns%totpftn(p)    = 0._r8
       end do
       if ( use_c14 ) then
          do p = bounds%begp,bounds%endp
             pc14s%dispvegc(p) = 0._r8
             pc14s%storvegc(p) = 0._r8
             pc14s%totpftc(p)  = 0._r8
          end do
       endif
       if ( use_c13 ) then
          do p = bounds%begp,bounds%endp
             pc13s%dispvegc(p) = 0._r8
             pc13s%storvegc(p) = 0._r8
             pc13s%totpftc(p)  = 0._r8
          end do
       endif
    end if
    
end subroutine CNZeroFluxes_dwt

!-----------------------------------------------------------------------
subroutine CNSetPps(num, filter, val, pps)
  !
  ! !DESCRIPTION:
  ! Set pft physical state variables
  ! !USES:
  use clm_varpar  , only : numrad
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (pft_pstate_type), intent(inout) :: pps
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i,j     ! loop index
  !------------------------------------------------------------------------

  ! currently NOT used

end subroutine CNSetPps

!-----------------------------------------------------------------------
subroutine CNSetPepv (num, filter, val, pepv)
  !
  ! !DESCRIPTION:
  ! Set pft ecophysiological variables
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (pft_epv_type), intent(inout) :: pepv
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i     ! loop index
  !------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      pepv%dormant_flag(i) = val
      pepv%days_active(i) = val
      pepv%onset_flag(i) = val
      pepv%onset_counter(i) = val
      pepv%onset_gddflag(i) = val
      pepv%onset_fdd(i) = val
      pepv%onset_gdd(i) = val
      pepv%onset_swi(i) = val
      pepv%offset_flag(i) = val
      pepv%offset_counter(i) = val
      pepv%offset_fdd(i) = val
      pepv%offset_swi(i) = val
      pepv%fert_counter(i) = val
      pepv%grain_flag(i) = val
      pepv%lgsf(i) = val
      pepv%bglfr(i) = val
      pepv%bgtr(i) = val
      pepv%annavg_t2m(i) = val
      pepv%tempavg_t2m(i) = val
      pepv%gpp(i) = val
      pepv%availc(i) = val
      pepv%xsmrpool_recover(i) = val
      pepv%alloc_pnow(i) = val
      pepv%c_allometry(i) = val
      pepv%n_allometry(i) = val
      pepv%plant_ndemand(i) = val
      pepv%tempsum_potential_gpp(i) = val
      pepv%annsum_potential_gpp(i) = val
      pepv%tempmax_retransn(i) = val
      pepv%annmax_retransn(i) = val
      pepv%avail_retransn(i) = val
      pepv%plant_nalloc(i) = val
      pepv%plant_calloc(i) = val
      pepv%excess_cflux(i) = val
      pepv%downreg(i) = val
      pepv%prev_leafc_to_litter(i) = val
      pepv%prev_frootc_to_litter(i) = val
      pepv%tempsum_npp(i) = val
      pepv%annsum_npp(i) = val
   end do
   if (use_cndv) then
      do fi = 1,num
         i = filter(fi)
         pepv%tempsum_litfall(i) = val
         pepv%annsum_litfall(i) = val
      end do
   end if
   if ( use_c13 ) then
      do fi = 1,num
         i = filter(fi)
         pepv%xsmrpool_c13ratio(i) = val
      end do
   endif

end subroutine CNSetPepv

!-----------------------------------------------------------------------
subroutine CNSetPcs (num, filter, val, pcs)
  !
  ! !DESCRIPTION:
  ! Set pft carbon state variables
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (pft_cstate_type), intent(inout) :: pcs
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i     ! loop index
  !------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      pcs%leafc(i) = val
      pcs%leafc_storage(i) = val
      pcs%leafc_xfer(i) = val
      pcs%frootc(i) = val
      pcs%frootc_storage(i) = val
      pcs%frootc_xfer(i) = val
      pcs%livestemc(i) = val
      pcs%livestemc_storage(i) = val
      pcs%livestemc_xfer(i) = val
      pcs%deadstemc(i) = val
      pcs%deadstemc_storage(i) = val
      pcs%deadstemc_xfer(i) = val
      pcs%livecrootc(i) = val
      pcs%livecrootc_storage(i) = val
      pcs%livecrootc_xfer(i) = val
      pcs%deadcrootc(i) = val
      pcs%deadcrootc_storage(i) = val
      pcs%deadcrootc_xfer(i) = val
      pcs%gresp_storage(i) = val
      pcs%gresp_xfer(i) = val
      pcs%cpool(i) = val
      pcs%xsmrpool(i) = val
      pcs%pft_ctrunc(i) = val
      pcs%dispvegc(i) = val
      pcs%storvegc(i) = val
      pcs%totvegc(i) = val
      pcs%totpftc(i) = val
      pcs%woodc(i) = val

   end do
   if ( crop_prog ) then
      do fi = 1,num
         i = filter(fi)
         pcs%grainc(i)         = val
         pcs%grainc_storage(i) = val
         pcs%grainc_xfer(i)    = val
      end do
   endif


end subroutine CNSetPcs

!-----------------------------------------------------------------------
subroutine CNSetPns(num, filter, val, pns)
  !
  ! !DESCRIPTION:
  ! Set pft nitrogen state variables
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (pft_nstate_type), intent(inout) :: pns
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i     ! loop index
  !------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      pns%leafn(i) = val
      pns%leafn_storage(i) = val
      pns%leafn_xfer(i) = val
      pns%frootn(i) = val
      pns%frootn_storage(i) = val
      pns%frootn_xfer(i) = val
      pns%livestemn(i) = val
      pns%livestemn_storage(i) = val
      pns%livestemn_xfer(i) = val
      pns%deadstemn(i) = val
      pns%deadstemn_storage(i) = val
      pns%deadstemn_xfer(i) = val
      pns%livecrootn(i) = val
      pns%livecrootn_storage(i) = val
      pns%livecrootn_xfer(i) = val
      pns%deadcrootn(i) = val
      pns%deadcrootn_storage(i) = val
      pns%deadcrootn_xfer(i) = val
      pns%retransn(i) = val
      pns%npool(i) = val
      pns%pft_ntrunc(i) = val
      pns%dispvegn(i) = val
      pns%storvegn(i) = val
      pns%totvegn(i) = val
      pns%totpftn(i) = val
   end do
   if ( crop_prog )then
      do fi = 1,num
         i = filter(fi)
         pns%grainn(i)         = val
         pns%grainn_storage(i) = val
         pns%grainn_xfer(i)    = val
      end do
   end if

end subroutine CNSetPns

!-----------------------------------------------------------------------
subroutine CNSetPcf(num, filter, val, pcf)
  !
  ! !DESCRIPTION:
  ! Set pft carbon flux variables
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (pft_cflux_type), intent(inout) :: pcf
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i     ! loop index
  !------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      pcf%m_leafc_to_litter(i) = val
      pcf%m_frootc_to_litter(i) = val
      pcf%m_leafc_storage_to_litter(i) = val
      pcf%m_frootc_storage_to_litter(i) = val
      pcf%m_livestemc_storage_to_litter(i) = val
      pcf%m_deadstemc_storage_to_litter(i) = val
      pcf%m_livecrootc_storage_to_litter(i) = val
      pcf%m_deadcrootc_storage_to_litter(i) = val
      pcf%m_leafc_xfer_to_litter(i) = val
      pcf%m_frootc_xfer_to_litter(i) = val
      pcf%m_livestemc_xfer_to_litter(i) = val
      pcf%m_deadstemc_xfer_to_litter(i) = val
      pcf%m_livecrootc_xfer_to_litter(i) = val
      pcf%m_deadcrootc_xfer_to_litter(i) = val
      pcf%m_livestemc_to_litter(i) = val
      pcf%m_deadstemc_to_litter(i) = val
      pcf%m_livecrootc_to_litter(i) = val
      pcf%m_deadcrootc_to_litter(i) = val
      pcf%m_gresp_storage_to_litter(i) = val
      pcf%m_gresp_xfer_to_litter(i) = val
      pcf%hrv_leafc_to_litter(i) = val             
      pcf%hrv_leafc_storage_to_litter(i) = val     
      pcf%hrv_leafc_xfer_to_litter(i) = val        
      pcf%hrv_frootc_to_litter(i) = val            
      pcf%hrv_frootc_storage_to_litter(i) = val    
      pcf%hrv_frootc_xfer_to_litter(i) = val       
      pcf%hrv_livestemc_to_litter(i) = val         
      pcf%hrv_livestemc_storage_to_litter(i) = val 
      pcf%hrv_livestemc_xfer_to_litter(i) = val    
      pcf%hrv_deadstemc_to_prod10c(i) = val        
      pcf%hrv_deadstemc_to_prod100c(i) = val       
      pcf%hrv_deadstemc_storage_to_litter(i) = val 
      pcf%hrv_deadstemc_xfer_to_litter(i) = val    
      pcf%hrv_livecrootc_to_litter(i) = val        
      pcf%hrv_livecrootc_storage_to_litter(i) = val
      pcf%hrv_livecrootc_xfer_to_litter(i) = val   
      pcf%hrv_deadcrootc_to_litter(i) = val        
      pcf%hrv_deadcrootc_storage_to_litter(i) = val
      pcf%hrv_deadcrootc_xfer_to_litter(i) = val   
      pcf%hrv_gresp_storage_to_litter(i) = val     
      pcf%hrv_gresp_xfer_to_litter(i) = val        
      pcf%hrv_xsmrpool_to_atm(i) = val
   
      ! fire-related variables changed by F. Li and S. Levis           
      pcf%m_leafc_to_fire(i) = val
      pcf%m_leafc_storage_to_fire(i) = val
      pcf%m_leafc_xfer_to_fire(i) = val
      pcf%m_livestemc_to_fire(i) = val
      pcf%m_livestemc_storage_to_fire(i) = val
      pcf%m_livestemc_xfer_to_fire(i) = val
      pcf%m_deadstemc_to_fire(i) = val
      pcf%m_deadstemc_storage_to_fire(i) = val
      pcf%m_deadstemc_xfer_to_fire(i) = val
      pcf%m_frootc_to_fire(i) = val
      pcf%m_frootc_storage_to_fire(i) = val
      pcf%m_frootc_xfer_to_fire(i) = val
      pcf%m_livecrootc_to_fire(i) = val
      pcf%m_livecrootc_storage_to_fire(i) = val
      pcf%m_livecrootc_xfer_to_fire(i) = val
      pcf%m_deadcrootc_to_fire(i) = val
      pcf%m_deadcrootc_storage_to_fire(i) = val
      pcf%m_deadcrootc_xfer_to_fire(i) = val
      pcf%m_gresp_storage_to_fire(i) = val
      pcf%m_gresp_xfer_to_fire(i) = val

      pcf%m_leafc_to_litter_fire(i) = val
      pcf%m_leafc_storage_to_litter_fire(i) = val
      pcf%m_leafc_xfer_to_litter_fire(i) = val
      pcf%m_livestemc_to_litter_fire(i) = val
      pcf%m_livestemc_storage_to_litter_fire(i) = val
      pcf%m_livestemc_xfer_to_litter_fire(i) = val
      pcf%m_livestemc_to_deadstemc_fire(i) = val
      pcf%m_deadstemc_to_litter_fire(i) = val
      pcf%m_deadstemc_storage_to_litter_fire(i) = val
      pcf%m_deadstemc_xfer_to_litter_fire(i) = val
      pcf%m_frootc_to_litter_fire(i) = val
      pcf%m_frootc_storage_to_litter_fire(i) = val
      pcf%m_frootc_xfer_to_litter_fire(i) = val
      pcf%m_livecrootc_to_litter_fire(i) = val
      pcf%m_livecrootc_storage_to_litter_fire(i) = val
      pcf%m_livecrootc_xfer_to_litter_fire(i) = val
      pcf%m_livecrootc_to_deadcrootc_fire(i) = val
      pcf%m_deadcrootc_to_litter_fire(i) = val
      pcf%m_deadcrootc_storage_to_litter_fire(i) = val
      pcf%m_deadcrootc_xfer_to_litter_fire(i) = val
      pcf%m_gresp_storage_to_litter_fire(i) = val
      pcf%m_gresp_xfer_to_litter_fire(i) = val


      pcf%leafc_xfer_to_leafc(i) = val
      pcf%frootc_xfer_to_frootc(i) = val
      pcf%livestemc_xfer_to_livestemc(i) = val
      pcf%deadstemc_xfer_to_deadstemc(i) = val
      pcf%livecrootc_xfer_to_livecrootc(i) = val
      pcf%deadcrootc_xfer_to_deadcrootc(i) = val
      pcf%leafc_to_litter(i) = val
      pcf%frootc_to_litter(i) = val
      pcf%leaf_mr(i) = val
      pcf%froot_mr(i) = val
      pcf%livestem_mr(i) = val
      pcf%livecroot_mr(i) = val
      pcf%grain_mr(i) = val
      pcf%leaf_curmr(i) = val
      pcf%froot_curmr(i) = val
      pcf%livestem_curmr(i) = val
      pcf%livecroot_curmr(i) = val
      pcf%grain_curmr(i) = val
      pcf%leaf_xsmr(i) = val
      pcf%froot_xsmr(i) = val
      pcf%livestem_xsmr(i) = val
      pcf%livecroot_xsmr(i) = val
      pcf%grain_xsmr(i) = val
      pcf%psnsun_to_cpool(i) = val
      pcf%psnshade_to_cpool(i) = val
      pcf%cpool_to_xsmrpool(i) = val
      pcf%cpool_to_leafc(i) = val
      pcf%cpool_to_leafc_storage(i) = val
      pcf%cpool_to_frootc(i) = val
      pcf%cpool_to_frootc_storage(i) = val
      pcf%cpool_to_livestemc(i) = val
      pcf%cpool_to_livestemc_storage(i) = val
      pcf%cpool_to_deadstemc(i) = val
      pcf%cpool_to_deadstemc_storage(i) = val
      pcf%cpool_to_livecrootc(i) = val
      pcf%cpool_to_livecrootc_storage(i) = val
      pcf%cpool_to_deadcrootc(i) = val
      pcf%cpool_to_deadcrootc_storage(i) = val
      pcf%cpool_to_gresp_storage(i) = val
      pcf%cpool_leaf_gr(i) = val
      pcf%cpool_leaf_storage_gr(i) = val
      pcf%transfer_leaf_gr(i) = val
      pcf%cpool_froot_gr(i) = val
      pcf%cpool_froot_storage_gr(i) = val
      pcf%transfer_froot_gr(i) = val
      pcf%cpool_livestem_gr(i) = val
      pcf%cpool_livestem_storage_gr(i) = val
      pcf%transfer_livestem_gr(i) = val
      pcf%cpool_deadstem_gr(i) = val
      pcf%cpool_deadstem_storage_gr(i) = val
      pcf%transfer_deadstem_gr(i) = val
      pcf%cpool_livecroot_gr(i) = val
      pcf%cpool_livecroot_storage_gr(i) = val
      pcf%transfer_livecroot_gr(i) = val
      pcf%cpool_deadcroot_gr(i) = val
      pcf%cpool_deadcroot_storage_gr(i) = val
      pcf%transfer_deadcroot_gr(i) = val
      pcf%leafc_storage_to_xfer(i) = val
      pcf%frootc_storage_to_xfer(i) = val
      pcf%livestemc_storage_to_xfer(i) = val
      pcf%deadstemc_storage_to_xfer(i) = val
      pcf%livecrootc_storage_to_xfer(i) = val
      pcf%deadcrootc_storage_to_xfer(i) = val
      pcf%gresp_storage_to_xfer(i) = val
      pcf%livestemc_to_deadstemc(i) = val
      pcf%livecrootc_to_deadcrootc(i) = val
      pcf%gpp(i) = val
      pcf%mr(i) = val
      pcf%current_gr(i) = val
      pcf%transfer_gr(i) = val
      pcf%storage_gr(i) = val
      pcf%gr(i) = val
      pcf%ar(i) = val
      pcf%rr(i) = val
      pcf%npp(i) = val
      pcf%agnpp(i) = val
      pcf%bgnpp(i) = val
      pcf%litfall(i) = val
      pcf%vegfire(i) = val
      pcf%wood_harvestc(i) = val
      pcf%pft_cinputs(i) = val
      pcf%pft_coutputs(i) = val
      pcf%pft_fire_closs(i) = val
      pcf%frootc_alloc(i) = val
      pcf%frootc_loss(i) = val
      pcf%leafc_alloc(i) = val
      pcf%leafc_loss(i) = val
      pcf%woodc_alloc(i) = val
      pcf%woodc_loss(i) = val
   end do
   if ( crop_prog )then
      do fi = 1,num
         i = filter(fi)
         pcf%xsmrpool_to_atm(i)         = val
         pcf%livestemc_to_litter(i)     = val
         pcf%grainc_to_food(i)          = val
         pcf%grainc_xfer_to_grainc(i)   = val
         pcf%cpool_to_grainc(i)         = val
         pcf%cpool_to_grainc_storage(i) = val
         pcf%cpool_grain_gr(i)          = val
         pcf%cpool_grain_storage_gr(i)  = val
         pcf%transfer_grain_gr(i)       = val
         pcf%grainc_storage_to_xfer(i)  = val
      end do
   end if

end subroutine CNSetPcf

!-----------------------------------------------------------------------
subroutine CNSetPnf(num, filter, val, pnf)
  !
  ! !DESCRIPTION:
  ! Set pft nitrogen flux variables
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (pft_nflux_type), intent(inout) :: pnf
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i     ! loop index
  !------------------------------------------------------------------------

   do fi = 1,num
      i=filter(fi)
      pnf%m_leafn_to_litter(i) = val
      pnf%m_frootn_to_litter(i) = val
      pnf%m_leafn_storage_to_litter(i) = val
      pnf%m_frootn_storage_to_litter(i) = val
      pnf%m_livestemn_storage_to_litter(i) = val
      pnf%m_deadstemn_storage_to_litter(i) = val
      pnf%m_livecrootn_storage_to_litter(i) = val
      pnf%m_deadcrootn_storage_to_litter(i) = val
      pnf%m_leafn_xfer_to_litter(i) = val
      pnf%m_frootn_xfer_to_litter(i) = val
      pnf%m_livestemn_xfer_to_litter(i) = val
      pnf%m_deadstemn_xfer_to_litter(i) = val
      pnf%m_livecrootn_xfer_to_litter(i) = val
      pnf%m_deadcrootn_xfer_to_litter(i) = val
      pnf%m_livestemn_to_litter(i) = val
      pnf%m_deadstemn_to_litter(i) = val
      pnf%m_livecrootn_to_litter(i) = val
      pnf%m_deadcrootn_to_litter(i) = val
      pnf%m_retransn_to_litter(i) = val
      pnf%hrv_leafn_to_litter(i) = val             
      pnf%hrv_frootn_to_litter(i) = val            
      pnf%hrv_leafn_storage_to_litter(i) = val     
      pnf%hrv_frootn_storage_to_litter(i) = val    
      pnf%hrv_livestemn_storage_to_litter(i) = val 
      pnf%hrv_deadstemn_storage_to_litter(i) = val 
      pnf%hrv_livecrootn_storage_to_litter(i) = val
      pnf%hrv_deadcrootn_storage_to_litter(i) = val
      pnf%hrv_leafn_xfer_to_litter(i) = val        
      pnf%hrv_frootn_xfer_to_litter(i) = val       
      pnf%hrv_livestemn_xfer_to_litter(i) = val    
      pnf%hrv_deadstemn_xfer_to_litter(i) = val    
      pnf%hrv_livecrootn_xfer_to_litter(i) = val   
      pnf%hrv_deadcrootn_xfer_to_litter(i) = val   
      pnf%hrv_livestemn_to_litter(i) = val         
      pnf%hrv_deadstemn_to_prod10n(i) = val        
      pnf%hrv_deadstemn_to_prod100n(i) = val       
      pnf%hrv_livecrootn_to_litter(i) = val        
      pnf%hrv_deadcrootn_to_litter(i) = val        
      pnf%hrv_retransn_to_litter(i) = val    

! fire-related variables changed by F. Li and S. Levis                  
      pnf%m_leafn_to_fire(i) = val
      pnf%m_leafn_storage_to_fire(i) = val
      pnf%m_leafn_xfer_to_fire(i) = val
      pnf%m_livestemn_to_fire(i) = val
      pnf%m_livestemn_storage_to_fire(i) = val
      pnf%m_livestemn_xfer_to_fire(i) = val
      pnf%m_deadstemn_to_fire(i) = val
      pnf%m_deadstemn_storage_to_fire(i) = val
      pnf%m_deadstemn_xfer_to_fire(i) = val
      pnf%m_frootn_to_fire(i) = val
      pnf%m_frootn_storage_to_fire(i) = val
      pnf%m_frootn_xfer_to_fire(i) = val
      pnf%m_livecrootn_to_fire(i) = val
      pnf%m_livecrootn_storage_to_fire(i) = val
      pnf%m_livecrootn_xfer_to_fire(i) = val
      pnf%m_deadcrootn_to_fire(i) = val
      pnf%m_deadcrootn_storage_to_fire(i) = val
      pnf%m_deadcrootn_xfer_to_fire(i) = val
      pnf%m_retransn_to_fire(i) = val
      
      
      pnf%m_leafn_to_litter_fire(i) = val
      pnf%m_leafn_storage_to_litter_fire(i) = val
      pnf%m_leafn_xfer_to_litter_fire(i) = val
      pnf%m_livestemn_to_litter_fire(i) = val
      pnf%m_livestemn_storage_to_litter_fire(i) = val
      pnf%m_livestemn_xfer_to_litter_fire(i) = val
      pnf%m_livestemn_to_deadstemn_fire(i) = val
      pnf%m_deadstemn_to_litter_fire(i) = val
      pnf%m_deadstemn_storage_to_litter_fire(i) = val
      pnf%m_deadstemn_xfer_to_litter_fire(i) = val
      pnf%m_frootn_to_litter_fire(i) = val
      pnf%m_frootn_storage_to_litter_fire(i) = val
      pnf%m_frootn_xfer_to_litter_fire(i) = val
      pnf%m_livecrootn_to_litter_fire(i) = val
      pnf%m_livecrootn_storage_to_litter_fire(i) = val
      pnf%m_livecrootn_xfer_to_litter_fire(i) = val
      pnf%m_livecrootn_to_deadcrootn_fire(i) = val
      pnf%m_deadcrootn_to_litter_fire(i) = val
      pnf%m_deadcrootn_storage_to_litter_fire(i) = val
      pnf%m_deadcrootn_xfer_to_litter_fire(i) = val
      pnf%m_retransn_to_litter_fire(i) = val

      pnf%leafn_xfer_to_leafn(i) = val
      pnf%frootn_xfer_to_frootn(i) = val
      pnf%livestemn_xfer_to_livestemn(i) = val
      pnf%deadstemn_xfer_to_deadstemn(i) = val
      pnf%livecrootn_xfer_to_livecrootn(i) = val
      pnf%deadcrootn_xfer_to_deadcrootn(i) = val
      pnf%leafn_to_litter(i) = val
      pnf%leafn_to_retransn(i) = val
      pnf%frootn_to_litter(i) = val
      pnf%retransn_to_npool(i) = val
      pnf%sminn_to_npool(i) = val
      pnf%npool_to_leafn(i) = val
      pnf%npool_to_leafn_storage(i) = val
      pnf%npool_to_frootn(i) = val
      pnf%npool_to_frootn_storage(i) = val
      pnf%npool_to_livestemn(i) = val
      pnf%npool_to_livestemn_storage(i) = val
      pnf%npool_to_deadstemn(i) = val
      pnf%npool_to_deadstemn_storage(i) = val
      pnf%npool_to_livecrootn(i) = val
      pnf%npool_to_livecrootn_storage(i) = val
      pnf%npool_to_deadcrootn(i) = val
      pnf%npool_to_deadcrootn_storage(i) = val
      pnf%leafn_storage_to_xfer(i) = val
      pnf%frootn_storage_to_xfer(i) = val
      pnf%livestemn_storage_to_xfer(i) = val
      pnf%deadstemn_storage_to_xfer(i) = val
      pnf%livecrootn_storage_to_xfer(i) = val
      pnf%deadcrootn_storage_to_xfer(i) = val
      pnf%livestemn_to_deadstemn(i) = val
      pnf%livestemn_to_retransn(i) = val
      pnf%livecrootn_to_deadcrootn(i) = val
      pnf%livecrootn_to_retransn(i) = val
      pnf%ndeploy(i) = val
      pnf%pft_ninputs(i) = val
      pnf%pft_noutputs(i) = val
      pnf%wood_harvestn(i) = val
      pnf%pft_fire_nloss(i) = val
   end do
   if ( crop_prog )then
      do fi = 1,num
         i=filter(fi)
         pnf%livestemn_to_litter(i)     = val
         pnf%grainn_to_food(i)          = val
         pnf%grainn_xfer_to_grainn(i)   = val
         pnf%npool_to_grainn(i)         = val
         pnf%npool_to_grainn_storage(i) = val
         pnf%grainn_storage_to_xfer(i)  = val
         pnf%soyfixn(i)                 = val
         pnf%frootn_to_retransn(i)      = val
      end do
   end if

end subroutine CNSetPnf

!-----------------------------------------------------------------------
subroutine CNSetCps(num, filter, val, cps)
  !
  ! !DESCRIPTION:
  ! Set column physical state variables
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (column_pstate_type), intent(inout) :: cps
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i,j     ! loop index
  !------------------------------------------------------------------------

   do fi = 1,num
      i = filter(fi)
      cps%coszen(i) = val
      cps%fpi(i) = val
      cps%fpg(i) = val
      cps%annsum_counter(i) = val
      cps%cannsum_npp(i) = val
      cps%cannavg_t2m(i) = val
      
  ! fire related variables changed by F. Li and S. Levis
      cps%wf(i) = val
      cps%wf2(i) = val
      cps%nfire(i) = val
      cps%baf_crop(i) = val
      cps%baf_peatf(i) = val
      cps%fbac(i) = val
      cps%fbac1(i) = val
      cps%farea_burned(i) = val
   end do

   do j = 1,nlevdecomp_full
      do fi = 1,num
         i = filter(fi)
         cps%fpi_vr(i,j) = val
      end do
   end do
   
   do j = 1,nlevgrnd
      do fi = 1,num
         i = filter(fi)
         cps%soilpsi(i,j) = val
      end do
   end do



end subroutine CNSetCps

!-----------------------------------------------------------------------
subroutine CNSetCcs(num, filter, val, ccs)
  !
  ! !DESCRIPTION:
  ! Set column carbon state variables
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (column_cstate_type), intent(inout) :: ccs
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i,j,k     ! loop index
  !------------------------------------------------------------------------

   ! column only
   do fi = 1,num
      i = filter(fi)
      ccs%cwdc(i)  = val
      ccs%col_ctrunc(i) = val
      ccs%totlitc(i) = val
      ccs%totsomc(i) = val
      ccs%totecosysc(i) = val
      ccs%totcolc(i) = val
      ccs%rootc_col(i) = val
      ccs%totvegc_col(i) = val
      ccs%leafc_col(i) = val
      ccs%fuelc(i) = val
      ccs%fuelc_crop(i) = val
      ccs%totlitc_1m(i) = val
      ccs%totsomc_1m(i) = val
   end do

   ! column and levdecomp
   do j = 1,nlevdecomp_full
      do fi = 1,num
         i = filter(fi)
         ccs%col_ctrunc_vr(i,j) = val
      end do
   end do

   ! column and decomp_pools
   do k = 1, ndecomp_pools
      do fi = 1,num
         i = filter(fi)
         ccs%decomp_cpools(i,k) = val
         ccs%decomp_cpools_1m(i,k) = val
      end do
   end do

   ! column, levdecomp, and decomp_pools
   do j = 1,nlevdecomp_full
      do k = 1, ndecomp_pools
         do fi = 1,num
            i = filter(fi)
            ccs%decomp_cpools_vr(i,j,k) = val
         end do
      end do
   end do

end subroutine CNSetCcs

!-----------------------------------------------------------------------
subroutine CNSetCns(num, filter, val, cns)
  !
  ! !DESCRIPTION:
  ! Set column nitrogen state variables
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (column_nstate_type), intent(inout) :: cns
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i,j,k     ! loop index
  !------------------------------------------------------------------------

   ! column only
   do fi = 1,num
      i = filter(fi)
      cns%sminn(i) = val
      cns%col_ntrunc(i) = val
      cns%cwdn(i) = val
      if (use_nitrif_denitrif) then
         cns%smin_no3(i) = val
         cns%smin_nh4(i) = val
      end if
      cns%totlitn(i) = val
      cns%totsomn(i) = val
      cns%totecosysn(i) = val
      cns%totcoln(i) = val
      cns%totsomn_1m(i) = val
      cns%totlitn_1m(i) = val
   end do
   
   ! column and levdecomp
   do j = 1,nlevdecomp_full
      do fi = 1,num
         i = filter(fi)
         cns%sminn_vr(i,j) = val
         cns%col_ntrunc_vr(i,j) = val
         if (use_nitrif_denitrif) then
            cns%smin_no3_vr(i,j) = val
            cns%smin_nh4_vr(i,j) = val
         end if
      end do
   end do
   
   ! column and decomp_pools
   do k = 1, ndecomp_pools
      do fi = 1,num
         i = filter(fi)
         cns%decomp_npools(i,k) = val
         cns%decomp_npools_1m(i,k) = val
      end do
   end do
   
   ! column levdecomp, and decomp_pools
   do j = 1,nlevdecomp_full
      do k = 1, ndecomp_pools
         do fi = 1,num
            i = filter(fi)
            cns%decomp_npools_vr(i,j,k) = val
         end do
      end do
   end do
   
end subroutine CNSetCns

!-----------------------------------------------------------------------
subroutine CNSetCcf(num, filter, val, ccf)
  !
  ! !DESCRIPTION:
  ! Set column carbon flux variables
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (column_cflux_type), intent(inout) :: ccf
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i,j,k,l     ! loop index
  !------------------------------------------------------------------------

   do j = 1, nlevdecomp_full
      do fi = 1,num
         i = filter(fi)
         ! phenology: litterfall and crop fluxes associated wit
         ccf%phenology_c_to_litr_met_c(i,j)                = val
         ccf%phenology_c_to_litr_cel_c(i,j)                = val
         ccf%phenology_c_to_litr_lig_c(i,j)                = val
         ! gap mortality
         ccf%gap_mortality_c_to_litr_met_c(i,j)                = val
         ccf%gap_mortality_c_to_litr_cel_c(i,j)                = val
         ccf%gap_mortality_c_to_litr_lig_c(i,j)                = val
         ccf%gap_mortality_c_to_cwdc(i,j)                      = val
         ! fire
         ccf%fire_mortality_c_to_cwdc(i,j)                 = val
         ccf%m_c_to_litr_met_fire(i,j)             = val
         ccf%m_c_to_litr_cel_fire(i,j)             = val  
         ccf%m_c_to_litr_lig_fire(i,j)             = val
         ! harvest
         ccf%harvest_c_to_litr_met_c(i,j)                  = val             
         ccf%harvest_c_to_litr_cel_c(i,j)                  = val             
         ccf%harvest_c_to_litr_lig_c(i,j)                  = val             
         ccf%harvest_c_to_cwdc(i,j)                        = val          
         ! hr
         ccf%hr_vr(i,j) = val
      end do
   end do

   do k = 1, ndecomp_pools
      do j = 1, nlevdecomp_full
         do fi = 1,num
            i = filter(fi)
            ccf%m_decomp_cpools_to_fire_vr(i,j,k) = val
            ccf%decomp_cpools_sourcesink(i,j,k) = val
            ccf%decomp_cpools_transport_tendency(i,j,k) = val
         end do
      end do
   end do

   do l = 1, ndecomp_cascade_transitions
         do fi = 1,num
            i = filter(fi)
            ccf%decomp_cascade_hr(i,l) = val
            ccf%decomp_cascade_ctransfer(i,l) = val
         end do
   end do

   do l = 1, ndecomp_cascade_transitions
      do j = 1, nlevdecomp_full
         do fi = 1,num
            i = filter(fi)
            ccf%decomp_cascade_hr_vr(i,j,l) = val
            ccf%decomp_cascade_ctransfer_vr(i,j,l) = val
            ccf%decomp_k(i,j,l) = val
         end do
      end do
   end do

   do k = 1, ndecomp_pools
      do fi = 1,num
         i = filter(fi)
         ccf%decomp_cpools_leached(i,k) = val
         ccf%m_decomp_cpools_to_fire(i,k) = val
      end do
   end do
   
   do fi = 1,num
      i = filter(fi)
      ccf%hrv_deadstemc_to_prod10c(i)         = val        
      ccf%hrv_deadstemc_to_prod100c(i)        = val  
      ccf%somc_fire(i)                        = val     ! F. Li and S.Levis
      ccf%prod10c_loss(i)                     = val
      ccf%prod100c_loss(i)                    = val
      ccf%product_closs(i)                    = val
      ccf%somhr(i)                            = val
      ccf%lithr(i)                            = val
      ccf%hr(i)                               = val
      ccf%sr(i)                               = val
      ccf%er(i)                               = val
      ccf%litfire(i)                          = val
      ccf%somfire(i)                          = val
      ccf%totfire(i)                          = val
      ccf%nep(i)                              = val
      ccf%nbp(i)                              = val
      ccf%nee(i)                              = val
      ccf%col_cinputs(i)                      = val
      ccf%col_coutputs(i)                     = val
      ccf%col_fire_closs(i)                   = val
      ccf%cwdc_hr(i)                          = val
      ccf%cwdc_loss(i)                        = val
      ccf%litterc_loss(i)                     = val
      ccf%som_c_leached(i) = val
  end do

end subroutine CNSetCcf

!-----------------------------------------------------------------------
subroutine CNSetCnf(num, filter, val, cnf)
  !
  ! !DESCRIPTION:
  ! Set column nitrogen flux variables
  !
  ! !ARGUMENTS:
  implicit none
  integer , intent(in) :: num
  integer , intent(in) :: filter(:)
  real(r8), intent(in) :: val
  type (column_nflux_type), intent(inout) :: cnf
  !
  ! !LOCAL VARIABLES:
  integer :: fi,i,j,k,l     ! loop index
  !------------------------------------------------------------------------


   do j = 1, nlevdecomp_full
      do fi = 1,num
         i = filter(fi)
         ! phenology: litterfall and crop fluxes associated wit
         cnf%phenology_n_to_litr_met_n(i,j)                = val
         cnf%phenology_n_to_litr_cel_n(i,j)                = val
         cnf%phenology_n_to_litr_lig_n(i,j)                = val
         ! gap mortality
         cnf%gap_mortality_n_to_litr_met_n(i,j)                = val
         cnf%gap_mortality_n_to_litr_cel_n(i,j)                = val
         cnf%gap_mortality_n_to_litr_lig_n(i,j)                = val
         cnf%gap_mortality_n_to_cwdn(i,j)                      = val
         ! fire
         cnf%fire_mortality_n_to_cwdn(i,j)                 = val
         cnf%m_n_to_litr_met_fire(i,j)             = val
         cnf%m_n_to_litr_cel_fire(i,j)             = val  
         cnf%m_n_to_litr_lig_fire(i,j)             = val
         ! harvest
         cnf%harvest_n_to_litr_met_n(i,j)                  = val             
         cnf%harvest_n_to_litr_cel_n(i,j)                  = val             
         cnf%harvest_n_to_litr_lig_n(i,j)                  = val             
         cnf%harvest_n_to_cwdn(i,j)                        = val  
         if (.not. use_nitrif_denitrif) then
            cnf%sminn_to_denit_excess_vr(i,j) = val
            cnf%sminn_leached_vr(i,j) = val
         else
            cnf%f_nit_vr(i,j) = val
            cnf%f_denit_vr(i,j) = val
            cnf%smin_no3_leached_vr(i,j) = val
            cnf%smin_no3_runoff_vr(i,j) = val 
            cnf%n2_n2o_ratio_denit_vr(i,j) = val
            cnf%pot_f_nit_vr(i,j) = val
            cnf%pot_f_denit_vr(i,j) = val
            cnf%actual_immob_no3_vr(i,j) = val
            cnf%actual_immob_nh4_vr(i,j) = val
            cnf%smin_no3_to_plant_vr(i,j) = val
            cnf%smin_nh4_to_plant_vr(i,j) = val
            cnf%f_n2o_denit_vr(i,j)  = val
            cnf%f_n2o_nit_vr(i,j)  = val

            cnf%smin_no3_massdens_vr(i,j) = val
            cnf%k_nitr_t_vr(i,j) = val
            cnf%k_nitr_ph_vr(i,j) = val
            cnf%k_nitr_h2o_vr(i,j) = val
            cnf%k_nitr_vr(i,j) = val 
            cnf%wfps_vr(i,j) = val 
            cnf%fmax_denit_carbonsubstrate_vr(i,j) = val 
            cnf%fmax_denit_nitrate_vr(i,j) = val 
            cnf%f_denit_base_vr(i,j) = val

            cnf%diffus(i,j) = val
            cnf%ratio_k1(i,j) = val
            cnf%ratio_no3_co2(i,j) = val
            cnf%soil_co2_prod(i,j) = val
            cnf%fr_WFPS(i,j) = val
            cnf%soil_bulkdensity(i,j) = val

            cnf%r_psi(i,j) = val
            cnf%anaerobic_frac(i,j) = val
         end if
         cnf%potential_immob_vr(i,j) = val
         cnf%actual_immob_vr(i,j) = val
         cnf%sminn_to_plant_vr(i,j) = val
         cnf%supplement_to_sminn_vr(i,j) = val
         cnf%gross_nmin_vr(i,j) = val
         cnf%net_nmin_vr(i,j) = val
      end do
   end do

   do fi = 1,num
      i = filter(fi)
      cnf%ndep_to_sminn(i) = val
      cnf%nfix_to_sminn(i) = val
      cnf%fert_to_sminn(i) = val
      cnf%soyfixn_to_sminn(i) = val
      cnf%hrv_deadstemn_to_prod10n(i) = val        
      cnf%hrv_deadstemn_to_prod100n(i) = val      
      cnf%prod10n_loss(i) = val
      cnf%prod100n_loss(i) = val
      cnf%product_nloss(i) = val
      cnf%potential_immob(i) = val
      cnf%actual_immob(i) = val
      cnf%sminn_to_plant(i) = val
      cnf%supplement_to_sminn(i) = val
      cnf%gross_nmin(i) = val
      cnf%net_nmin(i) = val
      cnf%denit(i) = val
      if (use_nitrif_denitrif) then
         cnf%f_nit(i) = val
         cnf%pot_f_nit(i) = val
         cnf%f_denit(i) = val
         cnf%pot_f_denit(i) = val
         cnf%f_n2o_denit(i) = val
         cnf%f_n2o_nit(i) = val
         cnf%smin_no3_leached(i) = val
         cnf%smin_no3_runoff(i) = val
      else
         cnf%sminn_to_denit_excess(i) = val
         cnf%sminn_leached(i) = val
      end if
      cnf%col_ninputs(i) = val
      cnf%col_noutputs(i) = val
      cnf%col_fire_nloss(i) = val
      cnf%som_n_leached(i) = val
   end do

   do k = 1, ndecomp_pools
      do fi = 1,num
         i = filter(fi)
         cnf%decomp_npools_leached(i,k) = val
         cnf%m_decomp_npools_to_fire(i,k) = val
      end do
   end do
   
   do k = 1, ndecomp_pools
      do j = 1, nlevdecomp_full
         do fi = 1,num
            i = filter(fi)
            cnf%m_decomp_npools_to_fire_vr(i,j,k) = val
            cnf%decomp_npools_sourcesink(i,j,k) = val
            cnf%decomp_npools_transport_tendency(i,j,k) = val
         end do
      end do
   end do
   
   do l = 1, ndecomp_cascade_transitions
         do fi = 1,num
            i = filter(fi)
            cnf%decomp_cascade_ntransfer(i,l) = val
            cnf%decomp_cascade_sminn_flux(i,l) = val
            if (.not. use_nitrif_denitrif) then
               cnf%sminn_to_denit_decomp_cascade(i,l) = val
            end if
         end do
   end do

   do l = 1, ndecomp_cascade_transitions
      do j = 1, nlevdecomp_full
         do fi = 1,num
            i = filter(fi)
            cnf%decomp_cascade_ntransfer_vr(i,j,l) = val
            cnf%decomp_cascade_sminn_flux_vr(i,j,l) = val
            if (.not. use_nitrif_denitrif) then
               cnf%sminn_to_denit_decomp_cascade_vr(i,j,l) = val
            end if
         end do
      end do
   end do


end subroutine CNSetCnf

end module CNSetValueMod
