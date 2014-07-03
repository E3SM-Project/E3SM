
module CNNStateUpdate2Mod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: NStateUpdate2Mod
!
! !DESCRIPTION:
! Module for nitrogen state variable update, mortality fluxes.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    use clm_varpar   , only: nlevsoi, nlevdecomp
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: NStateUpdate2
    public:: NStateUpdate2h
!
! !REVISION HISTORY:
! 4/23/2004: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NStateUpdate2
!
! !INTERFACE:
subroutine NStateUpdate2(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update all the prognostic nitrogen state
! variables affected by gap-phase mortality fluxes
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varctl  , only: iulog
   use clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 8/1/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   real(r8), pointer :: gap_mortality_n_to_litr_met_n(:,:)         ! N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
   real(r8), pointer :: gap_mortality_n_to_litr_cel_n(:,:)         ! N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
   real(r8), pointer :: gap_mortality_n_to_litr_lig_n(:,:)         ! N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
   real(r8), pointer :: gap_mortality_n_to_cwdn(:,:)               ! N fluxes associated with gap mortality to CWD pool (gN/m3/s)
   real(r8), pointer :: m_deadcrootn_storage_to_litter(:)
   real(r8), pointer :: m_deadcrootn_to_litter(:)
   real(r8), pointer :: m_deadcrootn_xfer_to_litter(:)
   real(r8), pointer :: m_deadstemn_storage_to_litter(:)
   real(r8), pointer :: m_deadstemn_to_litter(:)
   real(r8), pointer :: m_deadstemn_xfer_to_litter(:)
   real(r8), pointer :: m_frootn_storage_to_litter(:)
   real(r8), pointer :: m_frootn_to_litter(:)
   real(r8), pointer :: m_frootn_xfer_to_litter(:)
   real(r8), pointer :: m_leafn_storage_to_litter(:)
   real(r8), pointer :: m_leafn_to_litter(:)
   real(r8), pointer :: m_leafn_xfer_to_litter(:)
   real(r8), pointer :: m_livecrootn_storage_to_litter(:)
   real(r8), pointer :: m_livecrootn_to_litter(:)
   real(r8), pointer :: m_livecrootn_xfer_to_litter(:)
   real(r8), pointer :: m_livestemn_storage_to_litter(:)
   real(r8), pointer :: m_livestemn_to_litter(:)
   real(r8), pointer :: m_livestemn_xfer_to_litter(:)
   real(r8), pointer :: m_retransn_to_litter(:)
!
! local pointers to implicit in/out scalars
   real(r8), pointer :: decomp_npools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p,j,l         ! indices
   integer :: fp,fc       ! lake filter indices
   real(r8):: dt          ! radiation time step (seconds)

!EOP
!-----------------------------------------------------------------------
    ! assign local pointers at the column level
    gap_mortality_n_to_litr_met_n  => clm3%g%l%c%cnf%gap_mortality_n_to_litr_met_n
    gap_mortality_n_to_litr_cel_n  => clm3%g%l%c%cnf%gap_mortality_n_to_litr_cel_n
    gap_mortality_n_to_litr_lig_n  => clm3%g%l%c%cnf%gap_mortality_n_to_litr_lig_n
    gap_mortality_n_to_cwdn        => clm3%g%l%c%cnf%gap_mortality_n_to_cwdn
    decomp_npools_vr                   => clm3%g%l%c%cns%decomp_npools_vr
    ! assign local pointers at the pft level
    m_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter
    m_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%m_deadcrootn_to_litter
    m_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter
    m_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter
    m_deadstemn_to_litter          => clm3%g%l%c%p%pnf%m_deadstemn_to_litter
    m_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter
    m_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%m_frootn_storage_to_litter
    m_frootn_to_litter             => clm3%g%l%c%p%pnf%m_frootn_to_litter
    m_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter
    m_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%m_leafn_storage_to_litter
    m_leafn_to_litter              => clm3%g%l%c%p%pnf%m_leafn_to_litter
    m_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter
    m_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter
    m_livecrootn_to_litter         => clm3%g%l%c%p%pnf%m_livecrootn_to_litter
    m_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter
    m_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter
    m_livestemn_to_litter          => clm3%g%l%c%p%pnf%m_livestemn_to_litter
    m_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter
    m_retransn_to_litter           => clm3%g%l%c%p%pnf%m_retransn_to_litter
    deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    deadstemn                      => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    frootn                         => clm3%g%l%c%p%pns%frootn
    frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    leafn                          => clm3%g%l%c%p%pns%leafn
    leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    livecrootn                     => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    livestemn                      => clm3%g%l%c%p%pns%livestemn
    livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    retransn                       => clm3%g%l%c%p%pns%retransn

   ! set time steps
   dt = real( get_step_size(), r8 )

   do j = 1, nlevdecomp
      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         
         ! column-level nitrogen fluxes from gap-phase mortality
         decomp_npools_vr(c,j,i_met_lit) = decomp_npools_vr(c,j,i_met_lit) + gap_mortality_n_to_litr_met_n(c,j) * dt
         decomp_npools_vr(c,j,i_cel_lit) = decomp_npools_vr(c,j,i_cel_lit) + gap_mortality_n_to_litr_cel_n(c,j) * dt
         decomp_npools_vr(c,j,i_lig_lit) = decomp_npools_vr(c,j,i_lig_lit) + gap_mortality_n_to_litr_lig_n(c,j) * dt
         decomp_npools_vr(c,j,i_cwd) = decomp_npools_vr(c,j,i_cwd) + gap_mortality_n_to_cwdn(c,j)  * dt

      end do ! end of column loop
   end do

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! pft-level nitrogen fluxes from gap-phase mortality
      ! displayed pools
      leafn(p)      = leafn(p)      - m_leafn_to_litter(p)      * dt
      frootn(p)     = frootn(p)     - m_frootn_to_litter(p)     * dt
      livestemn(p)  = livestemn(p)  - m_livestemn_to_litter(p)  * dt
      deadstemn(p)  = deadstemn(p)  - m_deadstemn_to_litter(p)  * dt
      livecrootn(p) = livecrootn(p) - m_livecrootn_to_litter(p) * dt
      deadcrootn(p) = deadcrootn(p) - m_deadcrootn_to_litter(p) * dt
      retransn(p)   = retransn(p)   - m_retransn_to_litter(p)   * dt

      ! storage pools
      leafn_storage(p)      = leafn_storage(p)      - m_leafn_storage_to_litter(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - m_frootn_storage_to_litter(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - m_livestemn_storage_to_litter(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - m_deadstemn_storage_to_litter(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - m_livecrootn_storage_to_litter(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - m_deadcrootn_storage_to_litter(p) * dt

      ! transfer pools
      leafn_xfer(p)      = leafn_xfer(p)      - m_leafn_xfer_to_litter(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - m_frootn_xfer_to_litter(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - m_livestemn_xfer_to_litter(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - m_deadstemn_xfer_to_litter(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - m_livecrootn_xfer_to_litter(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - m_deadcrootn_xfer_to_litter(p) * dt

   end do

end subroutine NStateUpdate2
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NStateUpdate2h
!
! !INTERFACE:
subroutine NStateUpdate2h(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Update all the prognostic nitrogen state
! variables affected by harvest mortality fluxes
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 8/1/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   real(r8), pointer :: harvest_n_to_litr_met_n(:,:)               ! N fluxes associated with harvest to litter metabolic pool (gN/m3/s)
   real(r8), pointer :: harvest_n_to_litr_cel_n(:,:)               ! N fluxes associated with harvest to litter cellulose pool (gN/m3/s)
   real(r8), pointer :: harvest_n_to_litr_lig_n(:,:)               ! N fluxes associated with harvest to litter lignin pool (gN/m3/s)
   real(r8), pointer :: harvest_n_to_cwdn(:,:)                     ! N fluxes associated with harvest to CWD pool (gN/m3/s)
   real(r8), pointer :: hrv_deadcrootn_storage_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_to_litter(:)
   real(r8), pointer :: hrv_deadcrootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_storage_to_litter(:)
   real(r8), pointer :: hrv_deadstemn_to_prod10n(:)
   real(r8), pointer :: hrv_deadstemn_to_prod100n(:)
   real(r8), pointer :: hrv_deadstemn_xfer_to_litter(:)
   real(r8), pointer :: hrv_frootn_storage_to_litter(:)
   real(r8), pointer :: hrv_frootn_to_litter(:)
   real(r8), pointer :: hrv_frootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_leafn_storage_to_litter(:)
   real(r8), pointer :: hrv_leafn_to_litter(:)
   real(r8), pointer :: hrv_leafn_xfer_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_storage_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_to_litter(:)
   real(r8), pointer :: hrv_livecrootn_xfer_to_litter(:)
   real(r8), pointer :: hrv_livestemn_storage_to_litter(:)
   real(r8), pointer :: hrv_livestemn_to_litter(:)
   real(r8), pointer :: hrv_livestemn_xfer_to_litter(:)
   real(r8), pointer :: hrv_retransn_to_litter(:)
!
! local pointers to implicit in/out scalars
   real(r8), pointer :: decomp_npools_vr(:,:,:)    ! (gC/m3)  vertically-resolved decomposing (litter, cwd, soil) N pools
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p,j,l         ! indices
   integer :: fp,fc       ! lake filter indices
   real(r8):: dt          ! radiation time step (seconds)

!EOP
!-----------------------------------------------------------------------
    ! assign local pointers at the column level
    harvest_n_to_litr_met_n          => clm3%g%l%c%cnf%harvest_n_to_litr_met_n
    harvest_n_to_litr_cel_n          => clm3%g%l%c%cnf%harvest_n_to_litr_cel_n
    harvest_n_to_litr_lig_n          => clm3%g%l%c%cnf%harvest_n_to_litr_lig_n
    harvest_n_to_cwdn                => clm3%g%l%c%cnf%harvest_n_to_cwdn
    decomp_npools_vr                           => clm3%g%l%c%cns%decomp_npools_vr

    ! assign local pointers at the pft level
    hrv_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_deadcrootn_storage_to_litter
    hrv_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_deadcrootn_to_litter
    hrv_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_deadcrootn_xfer_to_litter
    hrv_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_deadstemn_storage_to_litter
    hrv_deadstemn_to_prod10n         => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod10n
    hrv_deadstemn_to_prod100n        => clm3%g%l%c%p%pnf%hrv_deadstemn_to_prod100n
    hrv_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_deadstemn_xfer_to_litter
    hrv_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%hrv_frootn_storage_to_litter
    hrv_frootn_to_litter             => clm3%g%l%c%p%pnf%hrv_frootn_to_litter
    hrv_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%hrv_frootn_xfer_to_litter
    hrv_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%hrv_leafn_storage_to_litter
    hrv_leafn_to_litter              => clm3%g%l%c%p%pnf%hrv_leafn_to_litter
    hrv_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%hrv_leafn_xfer_to_litter
    hrv_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%hrv_livecrootn_storage_to_litter
    hrv_livecrootn_to_litter         => clm3%g%l%c%p%pnf%hrv_livecrootn_to_litter
    hrv_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%hrv_livecrootn_xfer_to_litter
    hrv_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%hrv_livestemn_storage_to_litter
    hrv_livestemn_to_litter          => clm3%g%l%c%p%pnf%hrv_livestemn_to_litter
    hrv_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%hrv_livestemn_xfer_to_litter
    hrv_retransn_to_litter           => clm3%g%l%c%p%pnf%hrv_retransn_to_litter
    deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
    deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
    deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
    deadstemn                      => clm3%g%l%c%p%pns%deadstemn
    deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
    deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
    frootn                         => clm3%g%l%c%p%pns%frootn
    frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
    frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
    leafn                          => clm3%g%l%c%p%pns%leafn
    leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
    leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
    livecrootn                     => clm3%g%l%c%p%pns%livecrootn
    livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
    livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
    livestemn                      => clm3%g%l%c%p%pns%livestemn
    livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
    livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
    retransn                       => clm3%g%l%c%p%pns%retransn

   ! set time steps
   dt = real( get_step_size(), r8 )

   do j = 1,nlevdecomp
      ! column loop
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         
         ! column-level nitrogen fluxes from harvest mortality
         decomp_npools_vr(c,j,i_met_lit) = decomp_npools_vr(c,j,i_met_lit) + harvest_n_to_litr_met_n(c,j) * dt
         decomp_npools_vr(c,j,i_cel_lit) = decomp_npools_vr(c,j,i_cel_lit) + harvest_n_to_litr_cel_n(c,j) * dt
         decomp_npools_vr(c,j,i_lig_lit) = decomp_npools_vr(c,j,i_lig_lit) + harvest_n_to_litr_lig_n(c,j) * dt
         decomp_npools_vr(c,j,i_cwd) = decomp_npools_vr(c,j,i_cwd) + harvest_n_to_cwdn(c,j)  * dt
         
      end do ! end of column loop
   end do

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! pft-level nitrogen fluxes from harvest mortality
      ! displayed pools
      leafn(p)      = leafn(p)      - hrv_leafn_to_litter(p)      * dt
      frootn(p)     = frootn(p)     - hrv_frootn_to_litter(p)     * dt
      livestemn(p)  = livestemn(p)  - hrv_livestemn_to_litter(p)  * dt
      deadstemn(p)  = deadstemn(p)  - hrv_deadstemn_to_prod10n(p) * dt
      deadstemn(p)  = deadstemn(p)  - hrv_deadstemn_to_prod100n(p)* dt
      livecrootn(p) = livecrootn(p) - hrv_livecrootn_to_litter(p) * dt
      deadcrootn(p) = deadcrootn(p) - hrv_deadcrootn_to_litter(p) * dt
      retransn(p)   = retransn(p)   - hrv_retransn_to_litter(p)   * dt

      ! storage pools
      leafn_storage(p)      = leafn_storage(p)      - hrv_leafn_storage_to_litter(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - hrv_frootn_storage_to_litter(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - hrv_livestemn_storage_to_litter(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - hrv_deadstemn_storage_to_litter(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - hrv_livecrootn_storage_to_litter(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - hrv_deadcrootn_storage_to_litter(p) * dt

      ! transfer pools
      leafn_xfer(p)      = leafn_xfer(p)      - hrv_leafn_xfer_to_litter(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - hrv_frootn_xfer_to_litter(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - hrv_livestemn_xfer_to_litter(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - hrv_deadstemn_xfer_to_litter(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - hrv_livecrootn_xfer_to_litter(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - hrv_deadcrootn_xfer_to_litter(p) * dt

   end do

end subroutine NStateUpdate2h
!-----------------------------------------------------------------------

#endif

end module CNNStateUpdate2Mod
