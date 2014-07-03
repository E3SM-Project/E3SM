module CNC13StateUpdate3Mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: C13StateUpdate3Mod
!
! !DESCRIPTION:
! Module for carbon state variable update, mortality fluxes.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: C13StateUpdate3
!
! !REVISION HISTORY:
! 4/21/2005: Created by Peter Thornton and Neil Suits - copied from 
!  CNCStateUpdate3 for C13 state variables.
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13StateUpdate3
!
! !INTERFACE:
subroutine C13StateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update all the prognostic carbon state
! variables affected by fire fluxes
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varctl, only : use_c13
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
! 3/29/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real(r8), pointer :: m_cwdc_to_fire(:)
   real(r8), pointer :: m_deadcrootc_to_cwdc_fire(:)
   real(r8), pointer :: m_deadstemc_to_cwdc_fire(:)
   real(r8), pointer :: m_litr1c_to_fire(:)
   real(r8), pointer :: m_litr2c_to_fire(:)
   real(r8), pointer :: m_litr3c_to_fire(:)
   real(r8), pointer :: m_deadcrootc_storage_to_fire(:)
   real(r8), pointer :: m_deadcrootc_to_fire(:)
   real(r8), pointer :: m_deadcrootc_to_litter_fire(:)
   real(r8), pointer :: m_deadcrootc_xfer_to_fire(:)
   real(r8), pointer :: m_deadstemc_storage_to_fire(:)
   real(r8), pointer :: m_deadstemc_to_fire(:)
   real(r8), pointer :: m_deadstemc_to_litter_fire(:)
   real(r8), pointer :: m_deadstemc_xfer_to_fire(:)
   real(r8), pointer :: m_frootc_storage_to_fire(:)
   real(r8), pointer :: m_frootc_to_fire(:)
   real(r8), pointer :: m_frootc_xfer_to_fire(:)
   real(r8), pointer :: m_gresp_storage_to_fire(:)
   real(r8), pointer :: m_gresp_xfer_to_fire(:)
   real(r8), pointer :: m_leafc_storage_to_fire(:)
   real(r8), pointer :: m_leafc_to_fire(:)
   real(r8), pointer :: m_leafc_xfer_to_fire(:)
   real(r8), pointer :: m_livecrootc_storage_to_fire(:)
   real(r8), pointer :: m_livecrootc_to_fire(:)
   real(r8), pointer :: m_livecrootc_xfer_to_fire(:)
   real(r8), pointer :: m_livestemc_storage_to_fire(:)
   real(r8), pointer :: m_livestemc_to_fire(:)
   real(r8), pointer :: m_livestemc_xfer_to_fire(:)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real(r8), pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
!
! local pointers to implicit out arrays
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p      ! indices
   integer :: fp,fc    ! lake filter indices
   real(r8):: dt       ! radiation time step (seconds)

!EOP
!-----------------------------------------------------------------------

   if (.not. use_c13) then
      RETURN
   end if

    ! assign local pointers at the column level
    m_cwdc_to_fire                 => cc13f%m_cwdc_to_fire
    m_deadcrootc_to_cwdc_fire      => cc13f%m_deadcrootc_to_cwdc_fire
    m_deadstemc_to_cwdc_fire       => cc13f%m_deadstemc_to_cwdc_fire
    m_litr1c_to_fire               => cc13f%m_litr1c_to_fire
    m_litr2c_to_fire               => cc13f%m_litr2c_to_fire
    m_litr3c_to_fire               => cc13f%m_litr3c_to_fire
    cwdc                           => cc13s%cwdc
    litr1c                         => cc13s%litr1c
    litr2c                         => cc13s%litr2c
    litr3c                         => cc13s%litr3c

    ! assign local pointers at the column level
    m_deadcrootc_storage_to_fire   => pc13f%m_deadcrootc_storage_to_fire
    m_deadcrootc_to_fire           => pc13f%m_deadcrootc_to_fire
    m_deadcrootc_to_litter_fire    => pc13f%m_deadcrootc_to_litter_fire
    m_deadcrootc_xfer_to_fire      => pc13f%m_deadcrootc_xfer_to_fire
    m_deadstemc_storage_to_fire    => pc13f%m_deadstemc_storage_to_fire
    m_deadstemc_to_fire            => pc13f%m_deadstemc_to_fire
    m_deadstemc_to_litter_fire     => pc13f%m_deadstemc_to_litter_fire
    m_deadstemc_xfer_to_fire       => pc13f%m_deadstemc_xfer_to_fire
    m_frootc_storage_to_fire       => pc13f%m_frootc_storage_to_fire
    m_frootc_to_fire               => pc13f%m_frootc_to_fire
    m_frootc_xfer_to_fire          => pc13f%m_frootc_xfer_to_fire
    m_gresp_storage_to_fire        => pc13f%m_gresp_storage_to_fire
    m_gresp_xfer_to_fire           => pc13f%m_gresp_xfer_to_fire
    m_leafc_storage_to_fire        => pc13f%m_leafc_storage_to_fire
    m_leafc_to_fire                => pc13f%m_leafc_to_fire
    m_leafc_xfer_to_fire           => pc13f%m_leafc_xfer_to_fire
    m_livecrootc_storage_to_fire   => pc13f%m_livecrootc_storage_to_fire
    m_livecrootc_to_fire           => pc13f%m_livecrootc_to_fire
    m_livecrootc_xfer_to_fire      => pc13f%m_livecrootc_xfer_to_fire
    m_livestemc_storage_to_fire    => pc13f%m_livestemc_storage_to_fire
    m_livestemc_to_fire            => pc13f%m_livestemc_to_fire
    m_livestemc_xfer_to_fire       => pc13f%m_livestemc_xfer_to_fire
    deadcrootc                     => pc13s%deadcrootc
    deadcrootc_storage             => pc13s%deadcrootc_storage
    deadcrootc_xfer                => pc13s%deadcrootc_xfer
    deadstemc                      => pc13s%deadstemc
    deadstemc_storage              => pc13s%deadstemc_storage
    deadstemc_xfer                 => pc13s%deadstemc_xfer
    frootc                         => pc13s%frootc
    frootc_storage                 => pc13s%frootc_storage
    frootc_xfer                    => pc13s%frootc_xfer
    gresp_storage                  => pc13s%gresp_storage
    gresp_xfer                     => pc13s%gresp_xfer
    leafc                          => pc13s%leafc
    leafc_storage                  => pc13s%leafc_storage
    leafc_xfer                     => pc13s%leafc_xfer
    livecrootc                     => pc13s%livecrootc
    livecrootc_storage             => pc13s%livecrootc_storage
    livecrootc_xfer                => pc13s%livecrootc_xfer
    livestemc                      => pc13s%livestemc
    livestemc_storage              => pc13s%livestemc_storage
    livestemc_xfer                 => pc13s%livestemc_xfer

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! column loop
    do fc = 1,num_soilc
       c = filter_soilc(fc)

       ! column level carbon fluxes from fire

       ! pft-level wood to column-level CWD (uncombusted wood)
       cwdc(c) = cwdc(c) + m_deadstemc_to_cwdc_fire(c) * dt
       cwdc(c) = cwdc(c) + m_deadcrootc_to_cwdc_fire(c) * dt

       ! litter and CWD losses to fire
       litr1c(c) = litr1c(c) - m_litr1c_to_fire(c) * dt
       litr2c(c) = litr2c(c) - m_litr2c_to_fire(c) * dt
       litr3c(c) = litr3c(c) - m_litr3c_to_fire(c) * dt
       cwdc(c)   = cwdc(c)   - m_cwdc_to_fire(c)   * dt

    end do ! end of columns loop

    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)

       ! pft-level carbon fluxes from fire
       ! displayed pools
       leafc(p)               = leafc(p)               - m_leafc_to_fire(p)               * dt
       frootc(p)              = frootc(p)              - m_frootc_to_fire(p)              * dt
       livestemc(p)           = livestemc(p)           - m_livestemc_to_fire(p)           * dt
       deadstemc(p)           = deadstemc(p)           - m_deadstemc_to_fire(p)           * dt
       deadstemc(p)           = deadstemc(p)           - m_deadstemc_to_litter_fire(p)    * dt
       livecrootc(p)          = livecrootc(p)          - m_livecrootc_to_fire(p)          * dt
       deadcrootc(p)          = deadcrootc(p)          - m_deadcrootc_to_fire(p)          * dt
       deadcrootc(p)          = deadcrootc(p)          - m_deadcrootc_to_litter_fire(p)   * dt

       ! storage pools
       leafc_storage(p)       = leafc_storage(p)       - m_leafc_storage_to_fire(p)       * dt
       frootc_storage(p)      = frootc_storage(p)      - m_frootc_storage_to_fire(p)      * dt
       livestemc_storage(p)   = livestemc_storage(p)   - m_livestemc_storage_to_fire(p)   * dt
       deadstemc_storage(p)   = deadstemc_storage(p)   - m_deadstemc_storage_to_fire(p)   * dt
       livecrootc_storage(p)  = livecrootc_storage(p)  - m_livecrootc_storage_to_fire(p)  * dt
       deadcrootc_storage(p)  = deadcrootc_storage(p)  - m_deadcrootc_storage_to_fire(p)  * dt
       gresp_storage(p)       = gresp_storage(p)       - m_gresp_storage_to_fire(p)       * dt

       ! transfer pools
       leafc_xfer(p)      = leafc_xfer(p)      - m_leafc_xfer_to_fire(p)      * dt
       frootc_xfer(p)     = frootc_xfer(p)     - m_frootc_xfer_to_fire(p)     * dt
       livestemc_xfer(p)  = livestemc_xfer(p)  - m_livestemc_xfer_to_fire(p)  * dt
       deadstemc_xfer(p)  = deadstemc_xfer(p)  - m_deadstemc_xfer_to_fire(p)  * dt
       livecrootc_xfer(p) = livecrootc_xfer(p) - m_livecrootc_xfer_to_fire(p) * dt
       deadcrootc_xfer(p) = deadcrootc_xfer(p) - m_deadcrootc_xfer_to_fire(p) * dt
       gresp_xfer(p)      = gresp_xfer(p)      - m_gresp_xfer_to_fire(p)      * dt

    end do ! end of pft loop

end subroutine C13StateUpdate3
!-----------------------------------------------------------------------

end module CNC13StateUpdate3Mod
