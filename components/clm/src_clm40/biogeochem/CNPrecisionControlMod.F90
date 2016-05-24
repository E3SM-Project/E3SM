module CNPrecisionControlMod

!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: CNPrecisionControlMod
! 
! !DESCRIPTION:
! controls on very low values in critical state variables 
! 
! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    save
    private
! !PUBLIC MEMBER FUNCTIONS:
    public:: CNPrecisionControl
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
! !IROUTINE: CNPrecisionControl
!
! !INTERFACE:
subroutine CNPrecisionControl(num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION: 
! On the radiation time step, force leaf and deadstem c and n to 0 if
! they get too small.
!
! !USES:
   use clmtype
   use abortutils,   only: endrun
   use clm_varctl,   only: iulog, use_c13
   use pftvarcon,    only: nc3crop
   use surfrdMod,    only: crop_prog
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
   real(r8), pointer :: col_ctrunc(:)         ! (gC/m2) column-level sink for C truncation
   real(r8), pointer :: cwdc(:)               ! (gC/m2) coarse woody debris C
   real(r8), pointer :: litr1c(:)             ! (gC/m2) litter labile C
   real(r8), pointer :: litr2c(:)             ! (gC/m2) litter cellulose C
   real(r8), pointer :: litr3c(:)             ! (gC/m2) litter lignin C
   real(r8), pointer :: soil1c(:)             ! (gC/m2) soil organic matter C (fast pool)
   real(r8), pointer :: soil2c(:)             ! (gC/m2) soil organic matter C (medium pool)
   real(r8), pointer :: soil3c(:)             ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: soil4c(:)             ! (gC/m2) soil organic matter C (slowest pool)

   real(r8), pointer :: c13_col_ctrunc(:)     ! (gC/m2) column-level sink for C truncation
   real(r8), pointer :: c13_cwdc(:)           ! (gC/m2) coarse woody debris C
   real(r8), pointer :: c13_litr1c(:)         ! (gC/m2) litter labile C
   real(r8), pointer :: c13_litr2c(:)         ! (gC/m2) litter cellulose C
   real(r8), pointer :: c13_litr3c(:)         ! (gC/m2) litter lignin C
   real(r8), pointer :: c13_soil1c(:)         ! (gC/m2) soil organic matter C (fast pool)
   real(r8), pointer :: c13_soil2c(:)         ! (gC/m2) soil organic matter C (medium pool)
   real(r8), pointer :: c13_soil3c(:)         ! (gC/m2) soil organic matter C (slow pool)
   real(r8), pointer :: c13_soil4c(:)         ! (gC/m2) soil organic matter C (slowest pool)

   real(r8), pointer :: col_ntrunc(:)         ! (gN/m2) column-level sink for N truncation
   real(r8), pointer :: cwdn(:)               ! (gN/m2) coarse woody debris N
   real(r8), pointer :: litr1n(:)             ! (gN/m2) litter labile N
   real(r8), pointer :: litr2n(:)             ! (gN/m2) litter cellulose N
   real(r8), pointer :: litr3n(:)             ! (gN/m2) litter lignin N
   real(r8), pointer :: soil1n(:)             ! (gN/m2) soil organic matter N (fast pool)
   real(r8), pointer :: soil2n(:)             ! (gN/m2) soil organic matter N (medium pool)
   real(r8), pointer :: soil3n(:)             ! (gN/m2) soil orgainc matter N (slow pool)
   real(r8), pointer :: soil4n(:)             ! (gN/m2) soil orgainc matter N (slowest pool)
   real(r8), pointer :: cpool(:)              ! (gC/m2) temporary photosynthate C pool
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
   real(r8), pointer :: pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation
   real(r8), pointer :: xsmrpool(:)           ! (gC/m2) execss maint resp C pool
   real(r8), pointer :: grainc(:)             ! (gC/m2) grain C
   real(r8), pointer :: grainc_storage(:)     ! (gC/m2) grain C storage
   real(r8), pointer :: grainc_xfer(:)        ! (gC/m2) grain C transfer

   real(r8), pointer :: c13_cpool(:)              ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: c13_deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: c13_deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: c13_deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: c13_deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: c13_deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: c13_deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: c13_frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: c13_frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: c13_frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: c13_gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: c13_gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: c13_leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: c13_leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: c13_leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: c13_livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: c13_livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: c13_livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: c13_livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: c13_livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: c13_livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: c13_pft_ctrunc(:)         ! (gC/m2) pft-level sink for C truncation

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
   real(r8), pointer :: grainn(:)             ! (gC/m2) grain N
   real(r8), pointer :: grainn_storage(:)     ! (gC/m2) grain N storage
   real(r8), pointer :: grainn_xfer(:)        ! (gC/m2) grain N transfer
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: npool(:)              ! (gN/m2) temporary plant N pool
   real(r8), pointer :: pft_ntrunc(:)         ! (gN/m2) pft-level sink for N truncation
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   integer , pointer :: ivt(:)                ! pft vegetation type
!
! local pointers to implicit in/out scalars
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p      ! indices
   integer :: fp,fc    ! lake filter indices
   real(r8):: pc,pn    ! truncation terms for pft-level corrections
   real(r8):: cc,cn    ! truncation terms for column-level corrections

   real(r8):: pc13     ! truncation terms for pft-level corrections
   real(r8):: cc13     ! truncation terms for column-level corrections

   real(r8):: ccrit    ! critical carbon state value for truncation
   real(r8):: ncrit    ! critical nitrogen state value for truncation
    
!EOP
!-----------------------------------------------------------------------
    ! assign local pointers at the column level
    col_ctrunc                     => ccs%col_ctrunc
    cwdc                           => ccs%cwdc
    litr1c                         => ccs%litr1c
    litr2c                         => ccs%litr2c
    litr3c                         => ccs%litr3c
    soil1c                         => ccs%soil1c
    soil2c                         => ccs%soil2c
    soil3c                         => ccs%soil3c
    soil4c                         => ccs%soil4c

    c13_col_ctrunc                     => cc13s%col_ctrunc
    c13_cwdc                           => cc13s%cwdc
    c13_litr1c                         => cc13s%litr1c
    c13_litr2c                         => cc13s%litr2c
    c13_litr3c                         => cc13s%litr3c
    c13_soil1c                         => cc13s%soil1c
    c13_soil2c                         => cc13s%soil2c
    c13_soil3c                         => cc13s%soil3c
    c13_soil4c                         => cc13s%soil4c

    col_ntrunc                     => cns%col_ntrunc
    cwdn                           => cns%cwdn
    litr1n                         => cns%litr1n
    litr2n                         => cns%litr2n
    litr3n                         => cns%litr3n
    soil1n                         => cns%soil1n
    soil2n                         => cns%soil2n
    soil3n                         => cns%soil3n
    soil4n                         => cns%soil4n

    ! assign local pointers at the pft level
    ivt                            => pft%itype
    cpool                          => pcs%cpool
    deadcrootc                     => pcs%deadcrootc
    deadcrootc_storage             => pcs%deadcrootc_storage
    deadcrootc_xfer                => pcs%deadcrootc_xfer
    deadstemc                      => pcs%deadstemc
    deadstemc_storage              => pcs%deadstemc_storage
    deadstemc_xfer                 => pcs%deadstemc_xfer
    frootc                         => pcs%frootc
    frootc_storage                 => pcs%frootc_storage
    frootc_xfer                    => pcs%frootc_xfer
    gresp_storage                  => pcs%gresp_storage
    gresp_xfer                     => pcs%gresp_xfer
    leafc                          => pcs%leafc
    leafc_storage                  => pcs%leafc_storage
    leafc_xfer                     => pcs%leafc_xfer
    livecrootc                     => pcs%livecrootc
    livecrootc_storage             => pcs%livecrootc_storage
    livecrootc_xfer                => pcs%livecrootc_xfer
    livestemc                      => pcs%livestemc
    livestemc_storage              => pcs%livestemc_storage
    livestemc_xfer                 => pcs%livestemc_xfer
    pft_ctrunc                     => pcs%pft_ctrunc
    xsmrpool                       => pcs%xsmrpool
    grainc                         => pcs%grainc
    grainc_storage                 => pcs%grainc_storage
    grainc_xfer                    => pcs%grainc_xfer

    c13_cpool                          => pc13s%cpool
    c13_deadcrootc                     => pc13s%deadcrootc
    c13_deadcrootc_storage             => pc13s%deadcrootc_storage
    c13_deadcrootc_xfer                => pc13s%deadcrootc_xfer
    c13_deadstemc                      => pc13s%deadstemc
    c13_deadstemc_storage              => pc13s%deadstemc_storage
    c13_deadstemc_xfer                 => pc13s%deadstemc_xfer
    c13_frootc                         => pc13s%frootc
    c13_frootc_storage                 => pc13s%frootc_storage
    c13_frootc_xfer                    => pc13s%frootc_xfer
    c13_gresp_storage                  => pc13s%gresp_storage
    c13_gresp_xfer                     => pc13s%gresp_xfer
    c13_leafc                          => pc13s%leafc
    c13_leafc_storage                  => pc13s%leafc_storage
    c13_leafc_xfer                     => pc13s%leafc_xfer
    c13_livecrootc                     => pc13s%livecrootc
    c13_livecrootc_storage             => pc13s%livecrootc_storage
    c13_livecrootc_xfer                => pc13s%livecrootc_xfer
    c13_livestemc                      => pc13s%livestemc
    c13_livestemc_storage              => pc13s%livestemc_storage
    c13_livestemc_xfer                 => pc13s%livestemc_xfer
    c13_pft_ctrunc                     => pc13s%pft_ctrunc

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
    grainn                         => pns%grainn
    grainn_storage                 => pns%grainn_storage
    grainn_xfer                    => pns%grainn_xfer
    livestemn                      => pns%livestemn
    livestemn_storage              => pns%livestemn_storage
    livestemn_xfer                 => pns%livestemn_xfer
    npool                          => pns%npool
    pft_ntrunc                     => pns%pft_ntrunc
    retransn                       => pns%retransn
   
   ! set the critical carbon state value for truncation (gC/m2)
   ccrit = 1.e-8_r8
   ! set the critical nitrogen state value for truncation (gN/m2)
   ncrit = 1.e-8_r8
   
   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      
      ! initialize the pft-level C and N truncation terms
      pc = 0._r8
      if (use_c13) then
         pc13 = 0._r8
      end if
      pn = 0._r8
      
      ! do tests on state variables for precision control
      ! for linked C-N state variables, perform precision test on
      ! the C component, but truncate C, C13, and N components
      
      ! leaf C and N
      if (abs(leafc(p)) < ccrit) then
          pc = pc + leafc(p)
          leafc(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_leafc(p)
             c13_leafc(p) = 0._r8
          endif
          pn = pn + leafn(p)
          leafn(p) = 0._r8
      end if

      ! leaf storage C and N
      if (abs(leafc_storage(p)) < ccrit) then
          pc = pc + leafc_storage(p)
          leafc_storage(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_leafc_storage(p)
             c13_leafc_storage(p) = 0._r8
          endif
          pn = pn + leafn_storage(p)
          leafn_storage(p) = 0._r8
      end if
          
      ! leaf transfer C and N
      if (abs(leafc_xfer(p)) < ccrit) then
          pc = pc + leafc_xfer(p)
          leafc_xfer(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_leafc_xfer(p)
             c13_leafc_xfer(p) = 0._r8
          endif
          pn = pn + leafn_xfer(p)
          leafn_xfer(p) = 0._r8
      end if
          
      ! froot C and N
      if (abs(frootc(p)) < ccrit) then
          pc = pc + frootc(p)
          frootc(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_frootc(p)
             c13_frootc(p) = 0._r8
          endif
          pn = pn + frootn(p)
          frootn(p) = 0._r8
      end if

      ! froot storage C and N
      if (abs(frootc_storage(p)) < ccrit) then
          pc = pc + frootc_storage(p)
          frootc_storage(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_frootc_storage(p)
             c13_frootc_storage(p) = 0._r8
          endif
          pn = pn + frootn_storage(p)
          frootn_storage(p) = 0._r8
      end if
          
      ! froot transfer C and N
      if (abs(frootc_xfer(p)) < ccrit) then
          pc = pc + frootc_xfer(p)
          frootc_xfer(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_frootc_xfer(p)
             c13_frootc_xfer(p) = 0._r8
          endif
          pn = pn + frootn_xfer(p)
          frootn_xfer(p) = 0._r8
      end if
          
      if ( crop_prog .and. ivt(p) >= nc3crop )then
         ! grain C and N
         if (abs(grainc(p)) < ccrit) then
             pc = pc + grainc(p)
             grainc(p) = 0._r8
             pn = pn + grainn(p)
             grainn(p) = 0._r8
         end if
   
         ! grain storage C and N
         if (abs(grainc_storage(p)) < ccrit) then
             pc = pc + grainc_storage(p)
             grainc_storage(p) = 0._r8
             pn = pn + grainn_storage(p)
             grainn_storage(p) = 0._r8
         end if
             
         ! grain transfer C and N
         if (abs(grainc_xfer(p)) < ccrit) then
             pc = pc + grainc_xfer(p)
             grainc_xfer(p) = 0._r8
             pn = pn + grainn_xfer(p)
             grainn_xfer(p) = 0._r8
         end if
      end if
          
      ! livestem C and N
      if (abs(livestemc(p)) < ccrit) then
          pc = pc + livestemc(p)
          livestemc(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_livestemc(p)
             c13_livestemc(p) = 0._r8
          endif
          pn = pn + livestemn(p)
          livestemn(p) = 0._r8
      end if

      ! livestem storage C and N
      if (abs(livestemc_storage(p)) < ccrit) then
          pc = pc + livestemc_storage(p)
          livestemc_storage(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_livestemc_storage(p)
             c13_livestemc_storage(p) = 0._r8
          endif
          pn = pn + livestemn_storage(p)
          livestemn_storage(p) = 0._r8
      end if
          
      ! livestem transfer C and N
      if (abs(livestemc_xfer(p)) < ccrit) then
          pc = pc + livestemc_xfer(p)
          livestemc_xfer(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_livestemc_xfer(p)
             c13_livestemc_xfer(p) = 0._r8
          endif
          pn = pn + livestemn_xfer(p)
          livestemn_xfer(p) = 0._r8
      end if
          
      ! deadstem C and N
      if (abs(deadstemc(p)) < ccrit) then
          pc = pc + deadstemc(p)
          deadstemc(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_deadstemc(p)
             c13_deadstemc(p) = 0._r8
          endif
          pn = pn + deadstemn(p)
          deadstemn(p) = 0._r8
      end if

      ! deadstem storage C and N
      if (abs(deadstemc_storage(p)) < ccrit) then
          pc = pc + deadstemc_storage(p)
          deadstemc_storage(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_deadstemc_storage(p)
             c13_deadstemc_storage(p) = 0._r8
          endif
          pn = pn + deadstemn_storage(p)
          deadstemn_storage(p) = 0._r8
      end if
          
      ! deadstem transfer C and N
      if (abs(deadstemc_xfer(p)) < ccrit) then
          pc = pc + deadstemc_xfer(p)
          deadstemc_xfer(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_deadstemc_xfer(p)
             c13_deadstemc_xfer(p) = 0._r8
          endif
          pn = pn + deadstemn_xfer(p)
          deadstemn_xfer(p) = 0._r8
      end if
          
      ! livecroot C and N
      if (abs(livecrootc(p)) < ccrit) then
          pc = pc + livecrootc(p)
          livecrootc(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_livecrootc(p)
             c13_livecrootc(p) = 0._r8
          endif
          pn = pn + livecrootn(p)
          livecrootn(p) = 0._r8
      end if

      ! livecroot storage C and N
      if (abs(livecrootc_storage(p)) < ccrit) then
          pc = pc + livecrootc_storage(p)
          livecrootc_storage(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_livecrootc_storage(p)
             c13_livecrootc_storage(p) = 0._r8
          endif
          pn = pn + livecrootn_storage(p)
          livecrootn_storage(p) = 0._r8
      end if
          
      ! livecroot transfer C and N
      if (abs(livecrootc_xfer(p)) < ccrit) then
          pc = pc + livecrootc_xfer(p)
          livecrootc_xfer(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_livecrootc_xfer(p)
             c13_livecrootc_xfer(p) = 0._r8
          endif
          pn = pn + livecrootn_xfer(p)
          livecrootn_xfer(p) = 0._r8
      end if
          
      ! deadcroot C and N
      if (abs(deadcrootc(p)) < ccrit) then
          pc = pc + deadcrootc(p)
          deadcrootc(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_deadcrootc(p)
             c13_deadcrootc(p) = 0._r8
          endif
          pn = pn + deadcrootn(p)
          deadcrootn(p) = 0._r8
      end if

      ! deadcroot storage C and N
      if (abs(deadcrootc_storage(p)) < ccrit) then
          pc = pc + deadcrootc_storage(p)
          deadcrootc_storage(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_deadcrootc_storage(p)
             c13_deadcrootc_storage(p) = 0._r8
          endif
          pn = pn + deadcrootn_storage(p)
          deadcrootn_storage(p) = 0._r8
      end if
          
      ! deadcroot transfer C and N
      if (abs(deadcrootc_xfer(p)) < ccrit) then
          pc = pc + deadcrootc_xfer(p)
          deadcrootc_xfer(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_deadcrootc_xfer(p)
             c13_deadcrootc_xfer(p) = 0._r8
          endif
          pn = pn + deadcrootn_xfer(p)
          deadcrootn_xfer(p) = 0._r8
      end if
          
      ! gresp_storage (C only)
      if (abs(gresp_storage(p)) < ccrit) then
          pc = pc + gresp_storage(p)
          gresp_storage(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_gresp_storage(p)
             c13_gresp_storage(p) = 0._r8
          endif
      end if

      ! gresp_xfer (C only)
      if (abs(gresp_xfer(p)) < ccrit) then
          pc = pc + gresp_xfer(p)
          gresp_xfer(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_gresp_xfer(p)
             c13_gresp_xfer(p) = 0._r8
          endif
      end if
          
      ! cpool (C only)
      if (abs(cpool(p)) < ccrit) then
          pc = pc + cpool(p)
          cpool(p) = 0._r8
          if (use_c13) then
             pc13 = pc13 + c13_cpool(p)
             c13_cpool(p) = 0._r8
          endif
      end if
          
      if ( crop_prog .and. ivt(p) >= nc3crop )then
         ! xsmrpool (C only)
         if (abs(xsmrpool(p)) < ccrit) then
             pc = pc + xsmrpool(p)
             xsmrpool(p) = 0._r8
         end if
      end if
          
      ! retransn (N only)
      if (abs(retransn(p)) < ncrit) then
          pn = pn + retransn(p)
          retransn(p) = 0._r8
      end if
          
      ! npool (N only)
      if (abs(npool(p)) < ncrit) then
          pn = pn + npool(p)
          npool(p) = 0._r8
      end if
      
      pft_ctrunc(p) = pft_ctrunc(p) + pc
      if (use_c13) then
         c13_pft_ctrunc(p) = c13_pft_ctrunc(p) + pc13
      endif
      pft_ntrunc(p) = pft_ntrunc(p) + pn
          
   end do ! end of pft loop

   ! column loop
   do fc = 1,num_soilc
      c = filter_soilc(fc)
      
      ! initialize the column-level C and N truncation terms
      cc = 0._r8
      if (use_c13) then
         cc13 = 0._r8
      endif
      cn = 0._r8
      
      ! do tests on state variables for precision control
      ! for linked C-N state variables, perform precision test on
      ! the C component, but truncate both C and N components
      
      ! coarse woody debris C and N
      if (abs(cwdc(c)) < ccrit) then
          cc = cc + cwdc(c)
          cwdc(c) = 0._r8
          if (use_c13) then
             cc13 = cc13 + c13_cwdc(c)
             c13_cwdc(c) = 0._r8
          endif
          cn = cn + cwdn(c)
          cwdn(c) = 0._r8
      end if

      ! litr1 C and N
      if (abs(litr1c(c)) < ccrit) then
          cc = cc + litr1c(c)
          litr1c(c) = 0._r8
          if (use_c13) then
             cc13 = cc13 + c13_litr1c(c)
             c13_litr1c(c) = 0._r8
          endif
          cn = cn + litr1n(c)
          litr1n(c) = 0._r8
      end if

      ! litr2 C and N
      if (abs(litr2c(c)) < ccrit) then
          cc = cc + litr2c(c)
          litr2c(c) = 0._r8
          if (use_c13) then
             cc13 = cc13 + c13_litr2c(c)
             c13_litr2c(c) = 0._r8
          endif
          cn = cn + litr2n(c)
          litr2n(c) = 0._r8
      end if

      ! litr3 C and N
      if (abs(litr3c(c)) < ccrit) then
          cc = cc + litr3c(c)
          litr3c(c) = 0._r8
          if (use_c13) then
             cc13 = cc13 + c13_litr3c(c)
             c13_litr3c(c) = 0._r8
          endif
          cn = cn + litr3n(c)
          litr3n(c) = 0._r8
      end if

      ! soil1 C and N
      if (abs(soil1c(c)) < ccrit) then
          cc = cc + soil1c(c)
          soil1c(c) = 0._r8
          if (use_c13) then
             cc13 = cc13 + c13_soil1c(c)
             c13_soil1c(c) = 0._r8
          endif
          cn = cn + soil1n(c)
          soil1n(c) = 0._r8
      end if

      ! soil2 C and N
      if (abs(soil2c(c)) < ccrit) then
          cc = cc + soil2c(c)
          soil2c(c) = 0._r8
          if (use_c13) then
             cc13 = cc13 + c13_soil2c(c)
             c13_soil2c(c) = 0._r8
          endif
          cn = cn + soil2n(c)
          soil2n(c) = 0._r8
      end if

      ! soil3 C and N
      if (abs(soil3c(c)) < ccrit) then
          cc = cc + soil3c(c)
          soil3c(c) = 0._r8
          if (use_c13) then
             cc13 = cc13 + c13_soil3c(c)
             c13_soil3c(c) = 0._r8
          endif
          cn = cn + soil3n(c)
          soil3n(c) = 0._r8
      end if
      
      ! soil4 C and N
      if (abs(soil4c(c)) < ccrit) then
          cc = cc + soil4c(c)
          soil4c(c) = 0._r8
          if (use_c13) then
             cc13 = cc13 + c13_soil4c(c)
             c13_soil4c(c) = 0._r8
          endif
          cn = cn + soil4n(c)
          soil4n(c) = 0._r8
      end if
      
      ! not doing precision control on soil mineral N, since it will
      ! be getting the N truncation flux anyway.
      
      col_ctrunc(c) = col_ctrunc(c) + cc
      if (use_c13) then
         c13_col_ctrunc(c) = c13_col_ctrunc(c) + cc13
      endif
      col_ntrunc(c) = col_ntrunc(c) + cn
      
   end do   ! end of column loop

end subroutine CNPrecisionControl
!-----------------------------------------------------------------------

end module CNPrecisionControlMod
