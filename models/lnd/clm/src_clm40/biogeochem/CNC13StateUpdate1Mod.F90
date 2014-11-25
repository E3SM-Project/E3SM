module CNC13StateUpdate1Mod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: C13StateUpdate1Mod
!
! !DESCRIPTION:
! Module for carbon state variable update, non-mortality fluxes.
!
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
!
! !PUBLIC MEMBER FUNCTIONS:
    public:: C13StateUpdate1,C13StateUpdate0
!
! !REVISION HISTORY:
! 4/21/2005: Created by Peter Thornton and Neil Suits - copied from 
!  CNCStateUpdate1 for C13 state variables.
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13StateUpdate0
!
! !INTERFACE:
subroutine C13StateUpdate0(num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update cpool carbon state
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size
   use clm_varctl, only : use_c13
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 7/1/05: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
   real(r8), pointer :: psnshade_to_cpool(:)
   real(r8), pointer :: psnsun_to_cpool(:)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
! !OTHER LOCAL VARIABLES:
   integer :: p     ! indices
   integer :: fp   ! lake filter indices
   real(r8):: dt      ! radiation time step (seconds)
!
!EOP
!-----------------------------------------------------------------------

   if (.not. use_c13) then
      RETURN
   end if

    ! assign local pointers at the pft level
    cpool                          => pc13s%cpool
    psnshade_to_cpool              => pc13f%psnshade_to_cpool
    psnsun_to_cpool                => pc13f%psnsun_to_cpool

    ! set time steps
    dt = real( get_step_size(), r8 )

    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)
       ! gross photosynthesis fluxes
       cpool(p) = cpool(p) + psnsun_to_cpool(p)*dt
       cpool(p) = cpool(p) + psnshade_to_cpool(p)*dt
    end do

end subroutine C13StateUpdate0
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: C13StateUpdate1
!
! !INTERFACE:
subroutine C13StateUpdate1(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! On the radiation time step, update all the prognostic carbon state
! variables (except for gap-phase mortality and fire fluxes)
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_step_size

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
! 8/1/03: Created by Peter Thornton
! 12/5/03, Peter Thornton: Added livewood turnover fluxes
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   real(r8), pointer :: woody(:)       ! binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: cwdc_to_litr2c(:)
   real(r8), pointer :: cwdc_to_litr3c(:)
   real(r8), pointer :: frootc_to_litr1c(:)
   real(r8), pointer :: frootc_to_litr2c(:)
   real(r8), pointer :: frootc_to_litr3c(:)
   real(r8), pointer :: leafc_to_litr1c(:)
   real(r8), pointer :: leafc_to_litr2c(:)
   real(r8), pointer :: leafc_to_litr3c(:)
   real(r8), pointer :: litr1_hr(:)
   real(r8), pointer :: litr1c_to_soil1c(:)
   real(r8), pointer :: litr2_hr(:)
   real(r8), pointer :: litr2c_to_soil2c(:)
   real(r8), pointer :: litr3_hr(:)
   real(r8), pointer :: litr3c_to_soil3c(:)
   real(r8), pointer :: soil1_hr(:)
   real(r8), pointer :: soil1c_to_soil2c(:)
   real(r8), pointer :: soil2_hr(:)
   real(r8), pointer :: soil2c_to_soil3c(:)
   real(r8), pointer :: soil3_hr(:)
   real(r8), pointer :: soil3c_to_soil4c(:)
   real(r8), pointer :: soil4_hr(:)
   real(r8), pointer :: col_ctrunc(:)    ! (gC/m2) column-level sink for C truncation
   integer , pointer :: ivt(:)           ! pft vegetation type
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)
   real(r8), pointer :: leafc_xfer_to_leafc(:)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:)
   real(r8), pointer :: cpool_to_xsmrpool(:)
   real(r8), pointer :: cpool_to_deadcrootc(:)
   real(r8), pointer :: cpool_to_deadcrootc_storage(:)
   real(r8), pointer :: cpool_to_deadstemc(:)
   real(r8), pointer :: cpool_to_deadstemc_storage(:)
   real(r8), pointer :: cpool_to_frootc(:)
   real(r8), pointer :: cpool_to_frootc_storage(:)
   real(r8), pointer :: cpool_to_gresp_storage(:)
   real(r8), pointer :: cpool_to_leafc(:)
   real(r8), pointer :: cpool_to_leafc_storage(:)
   real(r8), pointer :: cpool_to_livecrootc(:)
   real(r8), pointer :: cpool_to_livecrootc_storage(:)
   real(r8), pointer :: cpool_to_livestemc(:)
   real(r8), pointer :: cpool_to_livestemc_storage(:)
   real(r8), pointer :: deadcrootc_storage_to_xfer(:)
   real(r8), pointer :: deadstemc_storage_to_xfer(:)
   real(r8), pointer :: frootc_storage_to_xfer(:)
   real(r8), pointer :: frootc_to_litter(:)
   real(r8), pointer :: gresp_storage_to_xfer(:)
   real(r8), pointer :: leafc_storage_to_xfer(:)
   real(r8), pointer :: leafc_to_litter(:)
   real(r8), pointer :: livecrootc_storage_to_xfer(:)
   real(r8), pointer :: livecrootc_to_deadcrootc(:)
   real(r8), pointer :: livestemc_storage_to_xfer(:)
   real(r8), pointer :: livestemc_to_deadstemc(:)
   real(r8), pointer :: livestem_mr(:)
   real(r8), pointer :: froot_mr(:)
   real(r8), pointer :: leaf_mr(:)
   real(r8), pointer :: livecroot_mr(:)
   real(r8), pointer :: livestem_curmr(:)
   real(r8), pointer :: froot_curmr(:)
   real(r8), pointer :: leaf_curmr(:)
   real(r8), pointer :: livecroot_curmr(:)
   real(r8), pointer :: livestem_xsmr(:)
   real(r8), pointer :: froot_xsmr(:)
   real(r8), pointer :: leaf_xsmr(:)
   real(r8), pointer :: livecroot_xsmr(:)
   real(r8), pointer :: cpool_deadcroot_gr(:)
   real(r8), pointer :: cpool_deadcroot_storage_gr(:)
   real(r8), pointer :: cpool_deadstem_gr(:)
   real(r8), pointer :: cpool_deadstem_storage_gr(:)
   real(r8), pointer :: cpool_froot_gr(:)
   real(r8), pointer :: cpool_froot_storage_gr(:)
   real(r8), pointer :: cpool_leaf_gr(:)
   real(r8), pointer :: cpool_leaf_storage_gr(:)
   real(r8), pointer :: cpool_livecroot_gr(:)
   real(r8), pointer :: cpool_livecroot_storage_gr(:)
   real(r8), pointer :: cpool_livestem_gr(:)
   real(r8), pointer :: cpool_livestem_storage_gr(:)
   real(r8), pointer :: transfer_deadcroot_gr(:)
   real(r8), pointer :: transfer_deadstem_gr(:)
   real(r8), pointer :: transfer_froot_gr(:)
   real(r8), pointer :: transfer_leaf_gr(:)
   real(r8), pointer :: transfer_livecroot_gr(:)
   real(r8), pointer :: transfer_livestem_gr(:)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: cwdc(:)          ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)        ! (gC/m2) litter labile C
   real(r8), pointer :: litr2c(:)        ! (gC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)        ! (gC/m2) litter lignin C
   real(r8), pointer :: soil1c(:)        ! (gC/m2) soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)        ! (gC/m2) soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)        ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)        ! (gC/m2) soil organic matter C (slowest pool)
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: xsmrpool(:)           ! (gC/m2) execss maint resp C pool
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

! local pointers for dynamic landcover fluxes and states
   real(r8), pointer :: dwt_seedc_to_leaf(:)
   real(r8), pointer :: dwt_seedc_to_deadstem(:)
   real(r8), pointer :: dwt_frootc_to_litr1c(:)
   real(r8), pointer :: dwt_frootc_to_litr2c(:)
   real(r8), pointer :: dwt_frootc_to_litr3c(:)
   real(r8), pointer :: dwt_livecrootc_to_cwdc(:)
   real(r8), pointer :: dwt_deadcrootc_to_cwdc(:)
   real(r8), pointer :: seedc(:)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p     ! indices
   integer :: fp,fc   ! lake filter indices
   real(r8):: dt      ! radiation time step (seconds)
!
!EOP
!-----------------------------------------------------------------------
    ! assign local pointers
    woody                          => pftcon%woody

    ! assign local pointers at the column level
    cwdc_to_litr2c                 => cc13f%cwdc_to_litr2c
    cwdc_to_litr3c                 => cc13f%cwdc_to_litr3c
    frootc_to_litr1c               => cc13f%frootc_to_litr1c
    frootc_to_litr2c               => cc13f%frootc_to_litr2c
    frootc_to_litr3c               => cc13f%frootc_to_litr3c
    leafc_to_litr1c                => cc13f%leafc_to_litr1c
    leafc_to_litr2c                => cc13f%leafc_to_litr2c
    leafc_to_litr3c                => cc13f%leafc_to_litr3c
    litr1_hr                       => cc13f%litr1_hr
    litr1c_to_soil1c               => cc13f%litr1c_to_soil1c
    litr2_hr                       => cc13f%litr2_hr
    litr2c_to_soil2c               => cc13f%litr2c_to_soil2c
    litr3_hr                       => cc13f%litr3_hr
    litr3c_to_soil3c               => cc13f%litr3c_to_soil3c
    soil1_hr                       => cc13f%soil1_hr
    soil1c_to_soil2c               => cc13f%soil1c_to_soil2c
    soil2_hr                       => cc13f%soil2_hr
    soil2c_to_soil3c               => cc13f%soil2c_to_soil3c
    soil3_hr                       => cc13f%soil3_hr
    soil3c_to_soil4c               => cc13f%soil3c_to_soil4c
    soil4_hr                       => cc13f%soil4_hr
    col_ctrunc                     => cc13s%col_ctrunc
    cwdc                           => cc13s%cwdc
    litr1c                         => cc13s%litr1c
    litr2c                         => cc13s%litr2c
    litr3c                         => cc13s%litr3c
    soil1c                         => cc13s%soil1c
    soil2c                         => cc13s%soil2c
    soil3c                         => cc13s%soil3c
    soil4c                         => cc13s%soil4c
    ! new pointers for dynamic landcover
    dwt_seedc_to_leaf              => cc13f%dwt_seedc_to_leaf
    dwt_seedc_to_deadstem          => cc13f%dwt_seedc_to_deadstem
    dwt_frootc_to_litr1c           => cc13f%dwt_frootc_to_litr1c
    dwt_frootc_to_litr2c           => cc13f%dwt_frootc_to_litr2c
    dwt_frootc_to_litr3c           => cc13f%dwt_frootc_to_litr3c
    dwt_livecrootc_to_cwdc         => cc13f%dwt_livecrootc_to_cwdc
    dwt_deadcrootc_to_cwdc         => cc13f%dwt_deadcrootc_to_cwdc
    seedc                          => cc13s%seedc

    ! assign local pointers at the pft level
    ivt                            => pft%itype
    cpool_deadcroot_gr             => pc13f%cpool_deadcroot_gr
    cpool_deadcroot_storage_gr     => pc13f%cpool_deadcroot_storage_gr
    cpool_deadstem_gr              => pc13f%cpool_deadstem_gr
    cpool_deadstem_storage_gr      => pc13f%cpool_deadstem_storage_gr
    cpool_froot_gr                 => pc13f%cpool_froot_gr
    cpool_froot_storage_gr         => pc13f%cpool_froot_storage_gr
    cpool_leaf_gr                  => pc13f%cpool_leaf_gr
    cpool_leaf_storage_gr          => pc13f%cpool_leaf_storage_gr
    cpool_livecroot_gr             => pc13f%cpool_livecroot_gr
    cpool_livecroot_storage_gr     => pc13f%cpool_livecroot_storage_gr
    cpool_livestem_gr              => pc13f%cpool_livestem_gr
    cpool_livestem_storage_gr      => pc13f%cpool_livestem_storage_gr
    cpool_to_xsmrpool              => pc13f%cpool_to_xsmrpool
    cpool_to_deadcrootc            => pc13f%cpool_to_deadcrootc
    cpool_to_deadcrootc_storage    => pc13f%cpool_to_deadcrootc_storage
    cpool_to_deadstemc             => pc13f%cpool_to_deadstemc
    cpool_to_deadstemc_storage     => pc13f%cpool_to_deadstemc_storage
    cpool_to_frootc                => pc13f%cpool_to_frootc
    cpool_to_frootc_storage        => pc13f%cpool_to_frootc_storage
    cpool_to_gresp_storage         => pc13f%cpool_to_gresp_storage
    cpool_to_leafc                 => pc13f%cpool_to_leafc
    cpool_to_leafc_storage         => pc13f%cpool_to_leafc_storage
    cpool_to_livecrootc            => pc13f%cpool_to_livecrootc
    cpool_to_livecrootc_storage    => pc13f%cpool_to_livecrootc_storage
    cpool_to_livestemc             => pc13f%cpool_to_livestemc
    cpool_to_livestemc_storage     => pc13f%cpool_to_livestemc_storage
    deadcrootc_storage_to_xfer     => pc13f%deadcrootc_storage_to_xfer
    deadcrootc_xfer_to_deadcrootc  => pc13f%deadcrootc_xfer_to_deadcrootc
    deadstemc_storage_to_xfer      => pc13f%deadstemc_storage_to_xfer
    deadstemc_xfer_to_deadstemc    => pc13f%deadstemc_xfer_to_deadstemc
    froot_mr                       => pc13f%froot_mr
    froot_curmr                    => pc13f%froot_curmr
    froot_xsmr                     => pc13f%froot_xsmr
    frootc_storage_to_xfer         => pc13f%frootc_storage_to_xfer
    frootc_to_litter               => pc13f%frootc_to_litter
    frootc_xfer_to_frootc          => pc13f%frootc_xfer_to_frootc
    gresp_storage_to_xfer          => pc13f%gresp_storage_to_xfer
    leaf_mr                        => pc13f%leaf_mr
    leaf_curmr                        => pc13f%leaf_curmr
    leaf_xsmr                        => pc13f%leaf_xsmr
    leafc_storage_to_xfer          => pc13f%leafc_storage_to_xfer
    leafc_to_litter                => pc13f%leafc_to_litter
    leafc_xfer_to_leafc            => pc13f%leafc_xfer_to_leafc
    livecroot_mr                   => pc13f%livecroot_mr
    livecroot_curmr                   => pc13f%livecroot_curmr
    livecroot_xsmr                   => pc13f%livecroot_xsmr
    livecrootc_storage_to_xfer     => pc13f%livecrootc_storage_to_xfer
    livecrootc_to_deadcrootc       => pc13f%livecrootc_to_deadcrootc
    livecrootc_xfer_to_livecrootc  => pc13f%livecrootc_xfer_to_livecrootc
    livestem_mr                    => pc13f%livestem_mr
    livestem_curmr                    => pc13f%livestem_curmr
    livestem_xsmr                    => pc13f%livestem_xsmr
    livestemc_storage_to_xfer      => pc13f%livestemc_storage_to_xfer
    livestemc_to_deadstemc         => pc13f%livestemc_to_deadstemc
    livestemc_xfer_to_livestemc    => pc13f%livestemc_xfer_to_livestemc
    transfer_deadcroot_gr          => pc13f%transfer_deadcroot_gr
    transfer_deadstem_gr           => pc13f%transfer_deadstem_gr
    transfer_froot_gr              => pc13f%transfer_froot_gr
    transfer_leaf_gr               => pc13f%transfer_leaf_gr
    transfer_livecroot_gr          => pc13f%transfer_livecroot_gr
    transfer_livestem_gr           => pc13f%transfer_livestem_gr
    cpool                          => pc13s%cpool
    xsmrpool                       => pc13s%xsmrpool
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
 
       ! column level fluxes
 
       ! plant to litter fluxes
       ! leaf litter
       litr1c(c) = litr1c(c) + leafc_to_litr1c(c)*dt
       litr2c(c) = litr2c(c) + leafc_to_litr2c(c)*dt
       litr3c(c) = litr3c(c) + leafc_to_litr3c(c)*dt
       ! fine root litter
       litr1c(c) = litr1c(c) + frootc_to_litr1c(c)*dt
       litr2c(c) = litr2c(c) + frootc_to_litr2c(c)*dt
       litr3c(c) = litr3c(c) + frootc_to_litr3c(c)*dt
 
       ! seeding fluxes, from dynamic landcover
	   seedc(c) = seedc(c) - dwt_seedc_to_leaf(c) * dt
	   seedc(c) = seedc(c) - dwt_seedc_to_deadstem(c) * dt
	   
       ! fluxes into litter and CWD, from dynamic landcover
       litr1c(c) = litr1c(c) + dwt_frootc_to_litr1c(c)*dt
       litr2c(c) = litr2c(c) + dwt_frootc_to_litr2c(c)*dt
       litr3c(c) = litr3c(c) + dwt_frootc_to_litr3c(c)*dt
       cwdc(c)   = cwdc(c)   + dwt_livecrootc_to_cwdc(c)*dt
       cwdc(c)   = cwdc(c)   + dwt_deadcrootc_to_cwdc(c)*dt
       
       ! litter and SOM HR fluxes
       litr1c(c) = litr1c(c) - litr1_hr(c)*dt
       litr2c(c) = litr2c(c) - litr2_hr(c)*dt
       litr3c(c) = litr3c(c) - litr3_hr(c)*dt
       soil1c(c) = soil1c(c) - soil1_hr(c)*dt
       soil2c(c) = soil2c(c) - soil2_hr(c)*dt
       soil3c(c) = soil3c(c) - soil3_hr(c)*dt
       soil4c(c) = soil4c(c) - soil4_hr(c)*dt
 
       ! CWD to litter fluxes
       cwdc(c)   = cwdc(c)   - cwdc_to_litr2c(c)*dt
       litr2c(c) = litr2c(c) + cwdc_to_litr2c(c)*dt
       cwdc(c)   = cwdc(c)   - cwdc_to_litr3c(c)*dt
       litr3c(c) = litr3c(c) + cwdc_to_litr3c(c)*dt
 
       ! litter to SOM fluxes
       litr1c(c) = litr1c(c) - litr1c_to_soil1c(c)*dt
       soil1c(c) = soil1c(c) + litr1c_to_soil1c(c)*dt
       litr2c(c) = litr2c(c) - litr2c_to_soil2c(c)*dt
       soil2c(c) = soil2c(c) + litr2c_to_soil2c(c)*dt
       litr3c(c) = litr3c(c) - litr3c_to_soil3c(c)*dt
       soil3c(c) = soil3c(c) + litr3c_to_soil3c(c)*dt
 
       ! SOM to SOM fluxes
       soil1c(c) = soil1c(c) - soil1c_to_soil2c(c)*dt
       soil2c(c) = soil2c(c) + soil1c_to_soil2c(c)*dt
       soil2c(c) = soil2c(c) - soil2c_to_soil3c(c)*dt
       soil3c(c) = soil3c(c) + soil2c_to_soil3c(c)*dt
       soil3c(c) = soil3c(c) - soil3c_to_soil4c(c)*dt
       soil4c(c) = soil4c(c) + soil3c_to_soil4c(c)*dt
 
    end do ! end of columns loop
 
    ! pft loop
    do fp = 1,num_soilp
       p = filter_soilp(fp)
 
       ! phenology: transfer growth fluxes
       leafc(p)           = leafc(p)       + leafc_xfer_to_leafc(p)*dt
       leafc_xfer(p)      = leafc_xfer(p)  - leafc_xfer_to_leafc(p)*dt
       frootc(p)          = frootc(p)      + frootc_xfer_to_frootc(p)*dt
       frootc_xfer(p)     = frootc_xfer(p) - frootc_xfer_to_frootc(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           livestemc(p)       = livestemc(p)           + livestemc_xfer_to_livestemc(p)*dt
           livestemc_xfer(p)  = livestemc_xfer(p)  - livestemc_xfer_to_livestemc(p)*dt
           deadstemc(p)       = deadstemc(p)           + deadstemc_xfer_to_deadstemc(p)*dt
           deadstemc_xfer(p)  = deadstemc_xfer(p)  - deadstemc_xfer_to_deadstemc(p)*dt
           livecrootc(p)      = livecrootc(p)          + livecrootc_xfer_to_livecrootc(p)*dt
           livecrootc_xfer(p) = livecrootc_xfer(p) - livecrootc_xfer_to_livecrootc(p)*dt
           deadcrootc(p)      = deadcrootc(p)          + deadcrootc_xfer_to_deadcrootc(p)*dt
           deadcrootc_xfer(p) = deadcrootc_xfer(p) - deadcrootc_xfer_to_deadcrootc(p)*dt
       end if
 
       ! phenology: litterfall fluxes
       leafc(p) = leafc(p) - leafc_to_litter(p)*dt
       frootc(p) = frootc(p) - frootc_to_litter(p)*dt
 
       ! livewood turnover fluxes
       if (woody(ivt(p)) == 1._r8) then
           livestemc(p)  = livestemc(p)  - livestemc_to_deadstemc(p)*dt
           deadstemc(p)  = deadstemc(p)  + livestemc_to_deadstemc(p)*dt
           livecrootc(p) = livecrootc(p) - livecrootc_to_deadcrootc(p)*dt
           deadcrootc(p) = deadcrootc(p) + livecrootc_to_deadcrootc(p)*dt
       end if
 
       ! maintenance respiration fluxes
       cpool(p) = cpool(p) - cpool_to_xsmrpool(p)*dt
       cpool(p) = cpool(p) - leaf_curmr(p)*dt
       cpool(p) = cpool(p) - froot_curmr(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           cpool(p) = cpool(p) - livestem_curmr(p)*dt
           cpool(p) = cpool(p) - livecroot_curmr(p)*dt
       end if
 
       ! maintenance respiration fluxes
       xsmrpool(p) = xsmrpool(p) + cpool_to_xsmrpool(p)*dt
       xsmrpool(p) = xsmrpool(p) - leaf_xsmr(p)*dt
       xsmrpool(p) = xsmrpool(p) - froot_xsmr(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           xsmrpool(p) = xsmrpool(p) - livestem_xsmr(p)*dt
           xsmrpool(p) = xsmrpool(p) - livecroot_xsmr(p)*dt
       end if
 
       ! allocation fluxes
       cpool(p)           = cpool(p)          - cpool_to_leafc(p)*dt
       leafc(p)           = leafc(p)          + cpool_to_leafc(p)*dt
       cpool(p)           = cpool(p)          - cpool_to_leafc_storage(p)*dt
       leafc_storage(p)   = leafc_storage(p)  + cpool_to_leafc_storage(p)*dt
       cpool(p)           = cpool(p)          - cpool_to_frootc(p)*dt
       frootc(p)          = frootc(p)         + cpool_to_frootc(p)*dt
       cpool(p)           = cpool(p)          - cpool_to_frootc_storage(p)*dt
       frootc_storage(p)  = frootc_storage(p) + cpool_to_frootc_storage(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           cpool(p)               = cpool(p)              - cpool_to_livestemc(p)*dt
           livestemc(p)           = livestemc(p)          + cpool_to_livestemc(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_livestemc_storage(p)*dt
           livestemc_storage(p)   = livestemc_storage(p)  + cpool_to_livestemc_storage(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_deadstemc(p)*dt
           deadstemc(p)           = deadstemc(p)          + cpool_to_deadstemc(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_deadstemc_storage(p)*dt
           deadstemc_storage(p)   = deadstemc_storage(p)  + cpool_to_deadstemc_storage(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_livecrootc(p)*dt
           livecrootc(p)          = livecrootc(p)         + cpool_to_livecrootc(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_livecrootc_storage(p)*dt
           livecrootc_storage(p)  = livecrootc_storage(p) + cpool_to_livecrootc_storage(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_deadcrootc(p)*dt
           deadcrootc(p)          = deadcrootc(p)         + cpool_to_deadcrootc(p)*dt
           cpool(p)               = cpool(p)              - cpool_to_deadcrootc_storage(p)*dt
           deadcrootc_storage(p)  = deadcrootc_storage(p) + cpool_to_deadcrootc_storage(p)*dt
       end if
 
       ! growth respiration fluxes for current growth
       cpool(p) = cpool(p) - cpool_leaf_gr(p)*dt
       cpool(p) = cpool(p) - cpool_froot_gr(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           cpool(p) = cpool(p) - cpool_livestem_gr(p)*dt
           cpool(p) = cpool(p) - cpool_deadstem_gr(p)*dt
           cpool(p) = cpool(p) - cpool_livecroot_gr(p)*dt
           cpool(p) = cpool(p) - cpool_deadcroot_gr(p)*dt
       end if
 
       ! growth respiration for transfer growth
       gresp_xfer(p) = gresp_xfer(p) - transfer_leaf_gr(p)*dt
       gresp_xfer(p) = gresp_xfer(p) - transfer_froot_gr(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           gresp_xfer(p) = gresp_xfer(p) - transfer_livestem_gr(p)*dt
           gresp_xfer(p) = gresp_xfer(p) - transfer_deadstem_gr(p)*dt
           gresp_xfer(p) = gresp_xfer(p) - transfer_livecroot_gr(p)*dt
           gresp_xfer(p) = gresp_xfer(p) - transfer_deadcroot_gr(p)*dt
       end if
 
       ! growth respiration at time of storage
       cpool(p) = cpool(p) - cpool_leaf_storage_gr(p)*dt
       cpool(p) = cpool(p) - cpool_froot_storage_gr(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           cpool(p) = cpool(p) - cpool_livestem_storage_gr(p)*dt
           cpool(p) = cpool(p) - cpool_deadstem_storage_gr(p)*dt
           cpool(p) = cpool(p) - cpool_livecroot_storage_gr(p)*dt
           cpool(p) = cpool(p) - cpool_deadcroot_storage_gr(p)*dt
       end if
 
       ! growth respiration stored for release during transfer growth
       cpool(p)         = cpool(p)         - cpool_to_gresp_storage(p)*dt
       gresp_storage(p) = gresp_storage(p) + cpool_to_gresp_storage(p)*dt
 
       ! move storage pools into transfer pools
       leafc_storage(p)   = leafc_storage(p)   - leafc_storage_to_xfer(p)*dt
       leafc_xfer(p)  = leafc_xfer(p)  + leafc_storage_to_xfer(p)*dt
       frootc_storage(p)  = frootc_storage(p)  - frootc_storage_to_xfer(p)*dt
       frootc_xfer(p) = frootc_xfer(p) + frootc_storage_to_xfer(p)*dt
       if (woody(ivt(p)) == 1._r8) then
           livestemc_storage(p)  = livestemc_storage(p)   - livestemc_storage_to_xfer(p)*dt
           livestemc_xfer(p)     = livestemc_xfer(p)  + livestemc_storage_to_xfer(p)*dt
           deadstemc_storage(p)  = deadstemc_storage(p)   - deadstemc_storage_to_xfer(p)*dt
           deadstemc_xfer(p)     = deadstemc_xfer(p)  + deadstemc_storage_to_xfer(p)*dt
           livecrootc_storage(p) = livecrootc_storage(p)  - livecrootc_storage_to_xfer(p)*dt
           livecrootc_xfer(p)    = livecrootc_xfer(p) + livecrootc_storage_to_xfer(p)*dt
           deadcrootc_storage(p) = deadcrootc_storage(p)  - deadcrootc_storage_to_xfer(p)*dt
           deadcrootc_xfer(p)    = deadcrootc_xfer(p) + deadcrootc_storage_to_xfer(p)*dt
           gresp_storage(p)      = gresp_storage(p)       - gresp_storage_to_xfer(p)*dt
           gresp_xfer(p)         = gresp_xfer(p)      + gresp_storage_to_xfer(p)*dt
       end if
 
    end do ! end of pft loop

end subroutine C13StateUpdate1
!-----------------------------------------------------------------------

end module CNC13StateUpdate1Mod
