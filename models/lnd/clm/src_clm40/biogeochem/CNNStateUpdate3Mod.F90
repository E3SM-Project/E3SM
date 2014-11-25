module CNNStateUpdate3Mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: NStateUpdate3Mod
!
! !DESCRIPTION:
! Module for nitrogen state variable update, mortality fluxes.
! Also, sminn leaching flux.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: NStateUpdate3
!
! !REVISION HISTORY:
! 7/27/2004: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: NStateUpdate3
!
! !INTERFACE:
subroutine NStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update all the prognostic nitrogen state
! variables affected by gap-phase mortality fluxes. Also the Sminn leaching flux.
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
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
   real(r8), pointer :: sminn_leached(:) 
   real(r8), pointer :: m_cwdn_to_fire(:)
   real(r8), pointer :: m_deadcrootn_to_cwdn_fire(:)
   real(r8), pointer :: m_deadstemn_to_cwdn_fire(:)
   real(r8), pointer :: m_litr1n_to_fire(:)
   real(r8), pointer :: m_litr2n_to_fire(:)
   real(r8), pointer :: m_litr3n_to_fire(:)
   real(r8), pointer :: m_deadcrootn_storage_to_fire(:)
   real(r8), pointer :: m_deadcrootn_to_fire(:)
   real(r8), pointer :: m_deadcrootn_to_litter_fire(:)
   real(r8), pointer :: m_deadcrootn_xfer_to_fire(:)
   real(r8), pointer :: m_deadstemn_storage_to_fire(:)
   real(r8), pointer :: m_deadstemn_to_fire(:)
   real(r8), pointer :: m_deadstemn_to_litter_fire(:)
   real(r8), pointer :: m_deadstemn_xfer_to_fire(:)
   real(r8), pointer :: m_frootn_storage_to_fire(:)
   real(r8), pointer :: m_frootn_to_fire(:)
   real(r8), pointer :: m_frootn_xfer_to_fire(:)
   real(r8), pointer :: m_leafn_storage_to_fire(:)
   real(r8), pointer :: m_leafn_to_fire(:)
   real(r8), pointer :: m_leafn_xfer_to_fire(:)
   real(r8), pointer :: m_livecrootn_storage_to_fire(:)
   real(r8), pointer :: m_livecrootn_to_fire(:)
   real(r8), pointer :: m_livecrootn_xfer_to_fire(:)
   real(r8), pointer :: m_livestemn_storage_to_fire(:)
   real(r8), pointer :: m_livestemn_to_fire(:)
   real(r8), pointer :: m_livestemn_xfer_to_fire(:)
   real(r8), pointer :: m_retransn_to_fire(:)
!
! local pointers to implicit in/out scalars
   real(r8), pointer :: sminn(:)              ! (gN/m2) soil mineral N
   real(r8), pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real(r8), pointer :: litr1n(:)             ! (gN/m2) litter labile N
   real(r8), pointer :: litr2n(:)             ! (gN/m2) litter cellulose N
   real(r8), pointer :: litr3n(:)             ! (gN/m2) litter lignin N
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
! !OTHER LOCAL VARIABLES:
   integer :: c,p        ! indices
   integer :: fp,fc      ! lake filter indices
   real(r8):: dt         ! radiation time step (seconds)

!EOP
!-----------------------------------------------------------------------

    ! assign local pointers at the column level
    sminn_leached                  => cnf%sminn_leached
    m_cwdn_to_fire                 => cnf%m_cwdn_to_fire
    m_deadcrootn_to_cwdn_fire      => cnf%m_deadcrootn_to_cwdn_fire
    m_deadstemn_to_cwdn_fire       => cnf%m_deadstemn_to_cwdn_fire
    m_litr1n_to_fire               => cnf%m_litr1n_to_fire
    m_litr2n_to_fire               => cnf%m_litr2n_to_fire
    m_litr3n_to_fire               => cnf%m_litr3n_to_fire
    sminn                          => cns%sminn
    cwdn                           => cns%cwdn
    litr1n                         => cns%litr1n
    litr2n                         => cns%litr2n
    litr3n                         => cns%litr3n

    ! assign local pointers at the pft level
    m_deadcrootn_storage_to_fire   => pnf%m_deadcrootn_storage_to_fire
    m_deadcrootn_to_fire           => pnf%m_deadcrootn_to_fire
    m_deadcrootn_to_litter_fire    => pnf%m_deadcrootn_to_litter_fire
    m_deadcrootn_xfer_to_fire      => pnf%m_deadcrootn_xfer_to_fire
    m_deadstemn_storage_to_fire    => pnf%m_deadstemn_storage_to_fire
    m_deadstemn_to_fire            => pnf%m_deadstemn_to_fire
    m_deadstemn_to_litter_fire     => pnf%m_deadstemn_to_litter_fire
    m_deadstemn_xfer_to_fire       => pnf%m_deadstemn_xfer_to_fire
    m_frootn_storage_to_fire       => pnf%m_frootn_storage_to_fire
    m_frootn_to_fire               => pnf%m_frootn_to_fire
    m_frootn_xfer_to_fire          => pnf%m_frootn_xfer_to_fire
    m_leafn_storage_to_fire        => pnf%m_leafn_storage_to_fire
    m_leafn_to_fire                => pnf%m_leafn_to_fire
    m_leafn_xfer_to_fire           => pnf%m_leafn_xfer_to_fire
    m_livecrootn_storage_to_fire   => pnf%m_livecrootn_storage_to_fire
    m_livecrootn_to_fire           => pnf%m_livecrootn_to_fire
    m_livecrootn_xfer_to_fire      => pnf%m_livecrootn_xfer_to_fire
    m_livestemn_storage_to_fire    => pnf%m_livestemn_storage_to_fire
    m_livestemn_to_fire            => pnf%m_livestemn_to_fire
    m_livestemn_xfer_to_fire       => pnf%m_livestemn_xfer_to_fire
    m_retransn_to_fire             => pnf%m_retransn_to_fire
    deadcrootn                     => pns%deadcrootn
    deadcrootn_storage             => pns%deadcrootn_storage
    deadcrootn_xfer                => pns%deadcrootn_xfer
    deadstemn                      => pns%deadstemn
    deadstemn_storage              => pns%deadstemn_storage
    deadstemn_xfer                 => pns%deadstemn_xfer
    frootn                         => pns%frootn
    frootn_storage                 => pns%frootn_storage
    frootn_xfer                    => pns%frootn_xfer
    leafn                          => pns%leafn
    leafn_storage                  => pns%leafn_storage
    leafn_xfer                     => pns%leafn_xfer
    livecrootn                     => pns%livecrootn
    livecrootn_storage             => pns%livecrootn_storage
    livecrootn_xfer                => pns%livecrootn_xfer
    livestemn                      => pns%livestemn
    livestemn_storage              => pns%livestemn_storage
    livestemn_xfer                 => pns%livestemn_xfer
    retransn                       => pns%retransn

   ! set time steps
   dt = real( get_step_size(), r8 )

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)

      ! mineral N loss due to leaching
      sminn(c) = sminn(c) - sminn_leached(c) * dt

      ! column level nitrogen fluxes from fire
      
      ! pft-level wood to column-level CWD (uncombusted wood)
      cwdn(c) = cwdn(c) + m_deadstemn_to_cwdn_fire(c) * dt
      cwdn(c) = cwdn(c) + m_deadcrootn_to_cwdn_fire(c) * dt

      ! litter and CWD losses to fire
      litr1n(c) = litr1n(c) - m_litr1n_to_fire(c) * dt
      litr2n(c) = litr2n(c) - m_litr2n_to_fire(c) * dt
      litr3n(c) = litr3n(c) - m_litr3n_to_fire(c) * dt
      cwdn(c)   = cwdn(c)   - m_cwdn_to_fire(c)   * dt

   end do ! end of column loop

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! pft-level nitrogen fluxes from fire
      ! displayed pools
      leafn(p)      = leafn(p)      - m_leafn_to_fire(p)             * dt
      frootn(p)     = frootn(p)     - m_frootn_to_fire(p)            * dt
      livestemn(p)  = livestemn(p)  - m_livestemn_to_fire(p)         * dt
      deadstemn(p)  = deadstemn(p)  - m_deadstemn_to_fire(p)         * dt
      deadstemn(p)  = deadstemn(p)  - m_deadstemn_to_litter_fire(p)  * dt
      livecrootn(p) = livecrootn(p) - m_livecrootn_to_fire(p)        * dt
      deadcrootn(p) = deadcrootn(p) - m_deadcrootn_to_fire(p)        * dt
      deadcrootn(p) = deadcrootn(p) - m_deadcrootn_to_litter_fire(p) * dt

      ! storage pools
      leafn_storage(p)      = leafn_storage(p)      - m_leafn_storage_to_fire(p)      * dt
      frootn_storage(p)     = frootn_storage(p)     - m_frootn_storage_to_fire(p)     * dt
      livestemn_storage(p)  = livestemn_storage(p)  - m_livestemn_storage_to_fire(p)  * dt
      deadstemn_storage(p)  = deadstemn_storage(p)  - m_deadstemn_storage_to_fire(p)  * dt
      livecrootn_storage(p) = livecrootn_storage(p) - m_livecrootn_storage_to_fire(p) * dt
      deadcrootn_storage(p) = deadcrootn_storage(p) - m_deadcrootn_storage_to_fire(p) * dt

      ! transfer pools
      leafn_xfer(p)      = leafn_xfer(p)      - m_leafn_xfer_to_fire(p)      * dt
      frootn_xfer(p)     = frootn_xfer(p)     - m_frootn_xfer_to_fire(p)     * dt
      livestemn_xfer(p)  = livestemn_xfer(p)  - m_livestemn_xfer_to_fire(p)  * dt
      deadstemn_xfer(p)  = deadstemn_xfer(p)  - m_deadstemn_xfer_to_fire(p)  * dt
      livecrootn_xfer(p) = livecrootn_xfer(p) - m_livecrootn_xfer_to_fire(p) * dt
      deadcrootn_xfer(p) = deadcrootn_xfer(p) - m_deadcrootn_xfer_to_fire(p) * dt

      ! retranslocated N pool
      retransn(p) = retransn(p) - m_retransn_to_fire(p) * dt

   end do

end subroutine NStateUpdate3
!-----------------------------------------------------------------------

end module CNNStateUpdate3Mod
