module CNGapMortalityMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNGapMortalityMod
!
! !DESCRIPTION:
! Module holding routines used in gap mortality for coupled carbon
! nitrogen code.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  save
  private
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNGapMortality
!
! !REVISION HISTORY:
! 3/29/04: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNGapMortality
!
! !INTERFACE:
subroutine CNGapMortality (num_soilc, filter_soilc, num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Gap-phase mortality routine for coupled carbon-nitrogen code (CN)
!
! !USES:
   use clmtype
   use clm_time_manager, only: get_days_per_year
   use clm_varcon      , only: secspday
   use clm_varctl      , only: use_cndv 
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! column filter for soil points
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! pft filter for soil points
!
! !CALLED FROM:
! subroutine CNEcosystemDyn
!
! !REVISION HISTORY:
! 3/29/04: Created by Peter Thornton
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in arrays
   integer , pointer :: ivt(:)         ! pft vegetation type
   real(r8), pointer :: woody(:)       ! binary flag for woody lifeform
                                       ! (1=woody, 0=not woody)
   real(r8), pointer :: leafc(:)              ! (gC/m2) leaf C
   real(r8), pointer :: frootc(:)             ! (gC/m2) fine root C
   real(r8), pointer :: livestemc(:)          ! (gC/m2) live stem C
   real(r8), pointer :: deadstemc(:)          ! (gC/m2) dead stem C
   real(r8), pointer :: livecrootc(:)         ! (gC/m2) live coarse root C
   real(r8), pointer :: deadcrootc(:)         ! (gC/m2) dead coarse root C
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: leafc_xfer(:)         ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc_xfer(:)        ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc_xfer(:)     ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc_xfer(:)     ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc_xfer(:)    ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc_xfer(:)    ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: gresp_xfer(:)         ! (gC/m2) growth respiration transfer
   real(r8), pointer :: leafn(:)              ! (gN/m2) leaf N
   real(r8), pointer :: frootn(:)             ! (gN/m2) fine root N
   real(r8), pointer :: livestemn(:)          ! (gN/m2) live stem N
   real(r8), pointer :: deadstemn(:)          ! (gN/m2) dead stem N
   real(r8), pointer :: livecrootn(:)         ! (gN/m2) live coarse root N
   real(r8), pointer :: deadcrootn(:)         ! (gN/m2) dead coarse root N
   real(r8), pointer :: retransn(:)           ! (gN/m2) plant pool of retranslocated N
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: leafn_xfer(:)         ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn_xfer(:)        ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn_xfer(:)     ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn_xfer(:)     ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn_xfer(:)    ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn_xfer(:)    ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: greffic(:)
   real(r8), pointer :: heatstress(:)
!
! local pointers to implicit in/out arrays
!
! local pointers to implicit out arrays
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
   real(r8), pointer :: m_leafn_to_litter(:)
   real(r8), pointer :: m_frootn_to_litter(:)
   real(r8), pointer :: m_livestemn_to_litter(:)
   real(r8), pointer :: m_deadstemn_to_litter(:)
   real(r8), pointer :: m_livecrootn_to_litter(:)
   real(r8), pointer :: m_deadcrootn_to_litter(:)
   real(r8), pointer :: m_retransn_to_litter(:)
   real(r8), pointer :: m_leafn_storage_to_litter(:)
   real(r8), pointer :: m_frootn_storage_to_litter(:)
   real(r8), pointer :: m_livestemn_storage_to_litter(:)
   real(r8), pointer :: m_deadstemn_storage_to_litter(:)
   real(r8), pointer :: m_livecrootn_storage_to_litter(:)
   real(r8), pointer :: m_deadcrootn_storage_to_litter(:)
   real(r8), pointer :: m_leafn_xfer_to_litter(:)
   real(r8), pointer :: m_frootn_xfer_to_litter(:)
   real(r8), pointer :: m_livestemn_xfer_to_litter(:)
   real(r8), pointer :: m_deadstemn_xfer_to_litter(:)
   real(r8), pointer :: m_livecrootn_xfer_to_litter(:)
   real(r8), pointer :: m_deadcrootn_xfer_to_litter(:)
!
! !OTHER LOCAL VARIABLES:
   integer :: p                         ! pft index
   integer :: fp                        ! pft filter index
   real(r8):: am                        ! rate for fractional mortality (1/yr)
   real(r8):: m                         ! rate for fractional mortality (1/s)
   real(r8):: mort_max                  ! asymptotic max mortality rate (/yr)
   real(r8), parameter :: k_mort = 0.3  !coeff of growth efficiency in mortality equation
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers
   woody                          => pftcon%woody

   ! assign local pointers to pft-level arrays
   ivt                            => pft%itype
   leafc                          => pcs%leafc
   frootc                         => pcs%frootc
   livestemc                      => pcs%livestemc
   deadstemc                      => pcs%deadstemc
   livecrootc                     => pcs%livecrootc
   deadcrootc                     => pcs%deadcrootc
   leafc_storage                  => pcs%leafc_storage
   frootc_storage                 => pcs%frootc_storage
   livestemc_storage              => pcs%livestemc_storage
   deadstemc_storage              => pcs%deadstemc_storage
   livecrootc_storage             => pcs%livecrootc_storage
   deadcrootc_storage             => pcs%deadcrootc_storage
   gresp_storage                  => pcs%gresp_storage
   leafc_xfer                     => pcs%leafc_xfer
   frootc_xfer                    => pcs%frootc_xfer
   livestemc_xfer                 => pcs%livestemc_xfer
   deadstemc_xfer                 => pcs%deadstemc_xfer
   livecrootc_xfer                => pcs%livecrootc_xfer
   deadcrootc_xfer                => pcs%deadcrootc_xfer
   gresp_xfer                     => pcs%gresp_xfer
   leafn                          => pns%leafn
   frootn                         => pns%frootn
   livestemn                      => pns%livestemn
   deadstemn                      => pns%deadstemn
   livecrootn                     => pns%livecrootn
   deadcrootn                     => pns%deadcrootn
   retransn                       => pns%retransn
   leafn_storage                  => pns%leafn_storage
   frootn_storage                 => pns%frootn_storage
   livestemn_storage              => pns%livestemn_storage
   deadstemn_storage              => pns%deadstemn_storage
   livecrootn_storage             => pns%livecrootn_storage
   deadcrootn_storage             => pns%deadcrootn_storage
   leafn_xfer                     => pns%leafn_xfer
   frootn_xfer                    => pns%frootn_xfer
   livestemn_xfer                 => pns%livestemn_xfer
   deadstemn_xfer                 => pns%deadstemn_xfer
   livecrootn_xfer                => pns%livecrootn_xfer
   deadcrootn_xfer                => pns%deadcrootn_xfer
   m_leafc_to_litter              => pcf%m_leafc_to_litter
   m_frootc_to_litter             => pcf%m_frootc_to_litter
   m_livestemc_to_litter          => pcf%m_livestemc_to_litter
   m_deadstemc_to_litter          => pcf%m_deadstemc_to_litter
   m_livecrootc_to_litter         => pcf%m_livecrootc_to_litter
   m_deadcrootc_to_litter         => pcf%m_deadcrootc_to_litter
   m_leafc_storage_to_litter      => pcf%m_leafc_storage_to_litter
   m_frootc_storage_to_litter     => pcf%m_frootc_storage_to_litter
   m_livestemc_storage_to_litter  => pcf%m_livestemc_storage_to_litter
   m_deadstemc_storage_to_litter  => pcf%m_deadstemc_storage_to_litter
   m_livecrootc_storage_to_litter => pcf%m_livecrootc_storage_to_litter
   m_deadcrootc_storage_to_litter => pcf%m_deadcrootc_storage_to_litter
   m_gresp_storage_to_litter      => pcf%m_gresp_storage_to_litter
   m_leafc_xfer_to_litter         => pcf%m_leafc_xfer_to_litter
   m_frootc_xfer_to_litter        => pcf%m_frootc_xfer_to_litter
   m_livestemc_xfer_to_litter     => pcf%m_livestemc_xfer_to_litter
   m_deadstemc_xfer_to_litter     => pcf%m_deadstemc_xfer_to_litter
   m_livecrootc_xfer_to_litter    => pcf%m_livecrootc_xfer_to_litter
   m_deadcrootc_xfer_to_litter    => pcf%m_deadcrootc_xfer_to_litter
   m_gresp_xfer_to_litter         => pcf%m_gresp_xfer_to_litter
   m_leafn_to_litter              => pnf%m_leafn_to_litter
   m_frootn_to_litter             => pnf%m_frootn_to_litter
   m_livestemn_to_litter          => pnf%m_livestemn_to_litter
   m_deadstemn_to_litter          => pnf%m_deadstemn_to_litter
   m_livecrootn_to_litter         => pnf%m_livecrootn_to_litter
   m_deadcrootn_to_litter         => pnf%m_deadcrootn_to_litter
   m_retransn_to_litter           => pnf%m_retransn_to_litter
   m_leafn_storage_to_litter      => pnf%m_leafn_storage_to_litter
   m_frootn_storage_to_litter     => pnf%m_frootn_storage_to_litter
   m_livestemn_storage_to_litter  => pnf%m_livestemn_storage_to_litter
   m_deadstemn_storage_to_litter  => pnf%m_deadstemn_storage_to_litter
   m_livecrootn_storage_to_litter => pnf%m_livecrootn_storage_to_litter
   m_deadcrootn_storage_to_litter => pnf%m_deadcrootn_storage_to_litter
   m_leafn_xfer_to_litter         => pnf%m_leafn_xfer_to_litter
   m_frootn_xfer_to_litter        => pnf%m_frootn_xfer_to_litter
   m_livestemn_xfer_to_litter     => pnf%m_livestemn_xfer_to_litter
   m_deadstemn_xfer_to_litter     => pnf%m_deadstemn_xfer_to_litter
   m_livecrootn_xfer_to_litter    => pnf%m_livecrootn_xfer_to_litter
   m_deadcrootn_xfer_to_litter    => pnf%m_deadcrootn_xfer_to_litter
   greffic                        => pdgvs%greffic
   heatstress                     => pdgvs%heatstress

   ! set the mortality rate based on annual rate
   am = 0.02_r8

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      if (use_cndv) then
         ! Stress mortality from lpj's subr Mortality.

         if (woody(ivt(p)) == 1._r8) then
            
            if (ivt(p) == 8) then
               mort_max = 0.03_r8 ! BDT boreal
            else
               mort_max = 0.01_r8 ! original value for all pfts
            end if
            
            ! heatstress and greffic calculated in Establishment once/yr
            
            ! Mortality rate inversely related to growth efficiency
            ! (Prentice et al 1993)
            am = mort_max / (1._r8 + k_mort * greffic(p))
            
            am = min(1._r8, am + heatstress(p))
         else ! lpj didn't set this for grasses; cn does
            ! set the mortality rate based on annual rate
            am = 0.02_r8
         end if
      end if

      m  = am/(get_days_per_year() * secspday)

      ! pft-level gap mortality carbon fluxes
      ! displayed pools
      m_leafc_to_litter(p)               = leafc(p)               * m
      m_frootc_to_litter(p)              = frootc(p)              * m
      m_livestemc_to_litter(p)           = livestemc(p)           * m
      m_deadstemc_to_litter(p)           = deadstemc(p)           * m
      m_livecrootc_to_litter(p)          = livecrootc(p)          * m
      m_deadcrootc_to_litter(p)          = deadcrootc(p)          * m

      ! storage pools
      m_leafc_storage_to_litter(p)       = leafc_storage(p)       * m
      m_frootc_storage_to_litter(p)      = frootc_storage(p)      * m
      m_livestemc_storage_to_litter(p)   = livestemc_storage(p)   * m
      m_deadstemc_storage_to_litter(p)   = deadstemc_storage(p)   * m
      m_livecrootc_storage_to_litter(p)  = livecrootc_storage(p)  * m
      m_deadcrootc_storage_to_litter(p)  = deadcrootc_storage(p)  * m
      m_gresp_storage_to_litter(p)       = gresp_storage(p)       * m

      ! transfer pools
      m_leafc_xfer_to_litter(p)          = leafc_xfer(p)          * m
      m_frootc_xfer_to_litter(p)         = frootc_xfer(p)         * m
      m_livestemc_xfer_to_litter(p)      = livestemc_xfer(p)      * m
      m_deadstemc_xfer_to_litter(p)      = deadstemc_xfer(p)      * m
      m_livecrootc_xfer_to_litter(p)     = livecrootc_xfer(p)     * m
      m_deadcrootc_xfer_to_litter(p)     = deadcrootc_xfer(p)     * m
      m_gresp_xfer_to_litter(p)          = gresp_xfer(p)          * m

      ! pft-level gap mortality nitrogen fluxes
      ! displayed pools
      m_leafn_to_litter(p)               = leafn(p)               * m
      m_frootn_to_litter(p)              = frootn(p)              * m
      m_livestemn_to_litter(p)           = livestemn(p)           * m
      m_deadstemn_to_litter(p)           = deadstemn(p)           * m
      m_livecrootn_to_litter(p)          = livecrootn(p)          * m
      m_deadcrootn_to_litter(p)          = deadcrootn(p)          * m
      m_retransn_to_litter(p)            = retransn(p)            * m

      ! storage pools
      m_leafn_storage_to_litter(p)       = leafn_storage(p)       * m
      m_frootn_storage_to_litter(p)      = frootn_storage(p)      * m
      m_livestemn_storage_to_litter(p)   = livestemn_storage(p)   * m
      m_deadstemn_storage_to_litter(p)   = deadstemn_storage(p)   * m
      m_livecrootn_storage_to_litter(p)  = livecrootn_storage(p)  * m
      m_deadcrootn_storage_to_litter(p)  = deadcrootn_storage(p)  * m

      ! transfer pools
      m_leafn_xfer_to_litter(p)          = leafn_xfer(p)          * m
      m_frootn_xfer_to_litter(p)         = frootn_xfer(p)         * m
      m_livestemn_xfer_to_litter(p)      = livestemn_xfer(p)      * m
      m_deadstemn_xfer_to_litter(p)      = deadstemn_xfer(p)      * m
      m_livecrootn_xfer_to_litter(p)     = livecrootn_xfer(p)     * m
      m_deadcrootn_xfer_to_litter(p)     = deadcrootn_xfer(p)     * m

   end do ! end of pft loop

   ! gather all pft-level litterfall fluxes to the column
   ! for litter C and N inputs

   call CNGapPftToColumn(num_soilc, filter_soilc)

end subroutine CNGapMortality
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNGapPftToColumn
!
! !INTERFACE:
subroutine CNGapPftToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! called in the middle of CNGapMoratlity to gather all pft-level gap mortality fluxes
! to the column level and assign them to the three litter pools
!
! !USES:
  use clmtype
  use clm_varpar, only : maxpatch_pft
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
   real(r8), pointer :: m_leafn_to_litter(:)
   real(r8), pointer :: m_frootn_to_litter(:)
   real(r8), pointer :: m_livestemn_to_litter(:)
   real(r8), pointer :: m_deadstemn_to_litter(:)
   real(r8), pointer :: m_livecrootn_to_litter(:)
   real(r8), pointer :: m_deadcrootn_to_litter(:)
   real(r8), pointer :: m_retransn_to_litter(:)
   real(r8), pointer :: m_leafn_storage_to_litter(:)
   real(r8), pointer :: m_frootn_storage_to_litter(:)
   real(r8), pointer :: m_livestemn_storage_to_litter(:)
   real(r8), pointer :: m_deadstemn_storage_to_litter(:)
   real(r8), pointer :: m_livecrootn_storage_to_litter(:)
   real(r8), pointer :: m_deadcrootn_storage_to_litter(:)
   real(r8), pointer :: m_leafn_xfer_to_litter(:)
   real(r8), pointer :: m_frootn_xfer_to_litter(:)
   real(r8), pointer :: m_livestemn_xfer_to_litter(:)
   real(r8), pointer :: m_deadstemn_xfer_to_litter(:)
   real(r8), pointer :: m_livecrootn_xfer_to_litter(:)
   real(r8), pointer :: m_deadcrootn_xfer_to_litter(:)
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
   real(r8), pointer :: m_leafn_to_litr1n(:)
   real(r8), pointer :: m_leafn_to_litr2n(:)
   real(r8), pointer :: m_leafn_to_litr3n(:)
   real(r8), pointer :: m_frootn_to_litr1n(:)
   real(r8), pointer :: m_frootn_to_litr2n(:)
   real(r8), pointer :: m_frootn_to_litr3n(:)
   real(r8), pointer :: m_livestemn_to_cwdn(:)
   real(r8), pointer :: m_deadstemn_to_cwdn(:)
   real(r8), pointer :: m_livecrootn_to_cwdn(:)
   real(r8), pointer :: m_deadcrootn_to_cwdn(:)
   real(r8), pointer :: m_retransn_to_litr1n(:)
   real(r8), pointer :: m_leafn_storage_to_litr1n(:)
   real(r8), pointer :: m_frootn_storage_to_litr1n(:)
   real(r8), pointer :: m_livestemn_storage_to_litr1n(:)
   real(r8), pointer :: m_deadstemn_storage_to_litr1n(:)
   real(r8), pointer :: m_livecrootn_storage_to_litr1n(:)
   real(r8), pointer :: m_deadcrootn_storage_to_litr1n(:)
   real(r8), pointer :: m_leafn_xfer_to_litr1n(:)
   real(r8), pointer :: m_frootn_xfer_to_litr1n(:)
   real(r8), pointer :: m_livestemn_xfer_to_litr1n(:)
   real(r8), pointer :: m_deadstemn_xfer_to_litr1n(:)
   real(r8), pointer :: m_livecrootn_xfer_to_litr1n(:)
   real(r8), pointer :: m_deadcrootn_xfer_to_litr1n(:)
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
   m_leafc_to_litr1c              => ccf%m_leafc_to_litr1c
   m_leafc_to_litr2c              => ccf%m_leafc_to_litr2c
   m_leafc_to_litr3c              => ccf%m_leafc_to_litr3c
   m_frootc_to_litr1c             => ccf%m_frootc_to_litr1c
   m_frootc_to_litr2c             => ccf%m_frootc_to_litr2c
   m_frootc_to_litr3c             => ccf%m_frootc_to_litr3c
   m_livestemc_to_cwdc            => ccf%m_livestemc_to_cwdc
   m_deadstemc_to_cwdc            => ccf%m_deadstemc_to_cwdc
   m_livecrootc_to_cwdc           => ccf%m_livecrootc_to_cwdc
   m_deadcrootc_to_cwdc           => ccf%m_deadcrootc_to_cwdc
   m_leafc_storage_to_litr1c      => ccf%m_leafc_storage_to_litr1c
   m_frootc_storage_to_litr1c     => ccf%m_frootc_storage_to_litr1c
   m_livestemc_storage_to_litr1c  => ccf%m_livestemc_storage_to_litr1c
   m_deadstemc_storage_to_litr1c  => ccf%m_deadstemc_storage_to_litr1c
   m_livecrootc_storage_to_litr1c => ccf%m_livecrootc_storage_to_litr1c
   m_deadcrootc_storage_to_litr1c => ccf%m_deadcrootc_storage_to_litr1c
   m_gresp_storage_to_litr1c      => ccf%m_gresp_storage_to_litr1c
   m_leafc_xfer_to_litr1c         => ccf%m_leafc_xfer_to_litr1c
   m_frootc_xfer_to_litr1c        => ccf%m_frootc_xfer_to_litr1c
   m_livestemc_xfer_to_litr1c     => ccf%m_livestemc_xfer_to_litr1c
   m_deadstemc_xfer_to_litr1c     => ccf%m_deadstemc_xfer_to_litr1c
   m_livecrootc_xfer_to_litr1c    => ccf%m_livecrootc_xfer_to_litr1c
   m_deadcrootc_xfer_to_litr1c    => ccf%m_deadcrootc_xfer_to_litr1c
   m_gresp_xfer_to_litr1c         => ccf%m_gresp_xfer_to_litr1c
   m_leafn_to_litr1n              => cnf%m_leafn_to_litr1n
   m_leafn_to_litr2n              => cnf%m_leafn_to_litr2n
   m_leafn_to_litr3n              => cnf%m_leafn_to_litr3n
   m_frootn_to_litr1n             => cnf%m_frootn_to_litr1n
   m_frootn_to_litr2n             => cnf%m_frootn_to_litr2n
   m_frootn_to_litr3n             => cnf%m_frootn_to_litr3n
   m_livestemn_to_cwdn            => cnf%m_livestemn_to_cwdn
   m_deadstemn_to_cwdn            => cnf%m_deadstemn_to_cwdn
   m_livecrootn_to_cwdn           => cnf%m_livecrootn_to_cwdn
   m_deadcrootn_to_cwdn           => cnf%m_deadcrootn_to_cwdn
   m_retransn_to_litr1n           => cnf%m_retransn_to_litr1n
   m_leafn_storage_to_litr1n      => cnf%m_leafn_storage_to_litr1n
   m_frootn_storage_to_litr1n     => cnf%m_frootn_storage_to_litr1n
   m_livestemn_storage_to_litr1n  => cnf%m_livestemn_storage_to_litr1n
   m_deadstemn_storage_to_litr1n  => cnf%m_deadstemn_storage_to_litr1n
   m_livecrootn_storage_to_litr1n => cnf%m_livecrootn_storage_to_litr1n
   m_deadcrootn_storage_to_litr1n => cnf%m_deadcrootn_storage_to_litr1n
   m_leafn_xfer_to_litr1n         => cnf%m_leafn_xfer_to_litr1n
   m_frootn_xfer_to_litr1n        => cnf%m_frootn_xfer_to_litr1n
   m_livestemn_xfer_to_litr1n     => cnf%m_livestemn_xfer_to_litr1n
   m_deadstemn_xfer_to_litr1n     => cnf%m_deadstemn_xfer_to_litr1n
   m_livecrootn_xfer_to_litr1n    => cnf%m_livecrootn_xfer_to_litr1n
   m_deadcrootn_xfer_to_litr1n    => cnf%m_deadcrootn_xfer_to_litr1n

   ! assign local pointers to pft-level arrays
   ivt                            => pft%itype
   wtcol                          => pft%wtcol
   pwtgcell                       => pft%wtgcell  
   m_leafc_to_litter              => pcf%m_leafc_to_litter
   m_frootc_to_litter             => pcf%m_frootc_to_litter
   m_livestemc_to_litter          => pcf%m_livestemc_to_litter
   m_deadstemc_to_litter          => pcf%m_deadstemc_to_litter
   m_livecrootc_to_litter         => pcf%m_livecrootc_to_litter
   m_deadcrootc_to_litter         => pcf%m_deadcrootc_to_litter
   m_leafc_storage_to_litter      => pcf%m_leafc_storage_to_litter
   m_frootc_storage_to_litter     => pcf%m_frootc_storage_to_litter
   m_livestemc_storage_to_litter  => pcf%m_livestemc_storage_to_litter
   m_deadstemc_storage_to_litter  => pcf%m_deadstemc_storage_to_litter
   m_livecrootc_storage_to_litter => pcf%m_livecrootc_storage_to_litter
   m_deadcrootc_storage_to_litter => pcf%m_deadcrootc_storage_to_litter
   m_gresp_storage_to_litter      => pcf%m_gresp_storage_to_litter
   m_leafc_xfer_to_litter         => pcf%m_leafc_xfer_to_litter
   m_frootc_xfer_to_litter        => pcf%m_frootc_xfer_to_litter
   m_livestemc_xfer_to_litter     => pcf%m_livestemc_xfer_to_litter
   m_deadstemc_xfer_to_litter     => pcf%m_deadstemc_xfer_to_litter
   m_livecrootc_xfer_to_litter    => pcf%m_livecrootc_xfer_to_litter
   m_deadcrootc_xfer_to_litter    => pcf%m_deadcrootc_xfer_to_litter
   m_gresp_xfer_to_litter         => pcf%m_gresp_xfer_to_litter
   m_leafn_to_litter              => pnf%m_leafn_to_litter
   m_frootn_to_litter             => pnf%m_frootn_to_litter
   m_livestemn_to_litter          => pnf%m_livestemn_to_litter
   m_deadstemn_to_litter          => pnf%m_deadstemn_to_litter
   m_livecrootn_to_litter         => pnf%m_livecrootn_to_litter
   m_deadcrootn_to_litter         => pnf%m_deadcrootn_to_litter
   m_retransn_to_litter           => pnf%m_retransn_to_litter
   m_leafn_storage_to_litter      => pnf%m_leafn_storage_to_litter
   m_frootn_storage_to_litter     => pnf%m_frootn_storage_to_litter
   m_livestemn_storage_to_litter  => pnf%m_livestemn_storage_to_litter
   m_deadstemn_storage_to_litter  => pnf%m_deadstemn_storage_to_litter
   m_livecrootn_storage_to_litter => pnf%m_livecrootn_storage_to_litter
   m_deadcrootn_storage_to_litter => pnf%m_deadcrootn_storage_to_litter
   m_leafn_xfer_to_litter         => pnf%m_leafn_xfer_to_litter
   m_frootn_xfer_to_litter        => pnf%m_frootn_xfer_to_litter
   m_livestemn_xfer_to_litter     => pnf%m_livestemn_xfer_to_litter
   m_deadstemn_xfer_to_litter     => pnf%m_deadstemn_xfer_to_litter
   m_livecrootn_xfer_to_litter    => pnf%m_livecrootn_xfer_to_litter
   m_deadcrootn_xfer_to_litter    => pnf%m_deadcrootn_xfer_to_litter

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

               ! leaf gap mortality nitrogen fluxes
               m_leafn_to_litr1n(c) = m_leafn_to_litr1n(c) + &
                  m_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               m_leafn_to_litr2n(c) = m_leafn_to_litr2n(c) + &
                  m_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               m_leafn_to_litr3n(c) = m_leafn_to_litr3n(c) + &
                  m_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root litter nitrogen fluxes
               m_frootn_to_litr1n(c) = m_frootn_to_litr1n(c) + &
                  m_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               m_frootn_to_litr2n(c) = m_frootn_to_litr2n(c) + &
                  m_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               m_frootn_to_litr3n(c) = m_frootn_to_litr3n(c) + &
                  m_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! wood gap mortality nitrogen fluxes
               m_livestemn_to_cwdn(c)  = m_livestemn_to_cwdn(c)  + &
                  m_livestemn_to_litter(p)  * wtcol(p)
               m_deadstemn_to_cwdn(c)  = m_deadstemn_to_cwdn(c)  + &
                  m_deadstemn_to_litter(p)  * wtcol(p)
               m_livecrootn_to_cwdn(c) = m_livecrootn_to_cwdn(c) + &
                  m_livecrootn_to_litter(p) * wtcol(p)
               m_deadcrootn_to_cwdn(c) = m_deadcrootn_to_cwdn(c) + &
                  m_deadcrootn_to_litter(p) * wtcol(p)

               ! retranslocated N pool gap mortality fluxes
               m_retransn_to_litr1n(c) = m_retransn_to_litr1n(c) + &
                  m_retransn_to_litter(p) * wtcol(p)

               ! storage gap mortality nitrogen fluxes
               m_leafn_storage_to_litr1n(c)      = m_leafn_storage_to_litr1n(c)      + &
                  m_leafn_storage_to_litter(p)      * wtcol(p)
               m_frootn_storage_to_litr1n(c)     = m_frootn_storage_to_litr1n(c)     + &
                  m_frootn_storage_to_litter(p)     * wtcol(p)
               m_livestemn_storage_to_litr1n(c)  = m_livestemn_storage_to_litr1n(c)  + &
                  m_livestemn_storage_to_litter(p)  * wtcol(p)
               m_deadstemn_storage_to_litr1n(c)  = m_deadstemn_storage_to_litr1n(c)  + &
                  m_deadstemn_storage_to_litter(p)  * wtcol(p)
               m_livecrootn_storage_to_litr1n(c) = m_livecrootn_storage_to_litr1n(c) + &
                  m_livecrootn_storage_to_litter(p) * wtcol(p)
               m_deadcrootn_storage_to_litr1n(c) = m_deadcrootn_storage_to_litr1n(c) + &
                  m_deadcrootn_storage_to_litter(p) * wtcol(p)

               ! transfer gap mortality nitrogen fluxes
               m_leafn_xfer_to_litr1n(c)      = m_leafn_xfer_to_litr1n(c)      + &
                  m_leafn_xfer_to_litter(p)      * wtcol(p)
               m_frootn_xfer_to_litr1n(c)     = m_frootn_xfer_to_litr1n(c)     + &
                  m_frootn_xfer_to_litter(p)     * wtcol(p)
               m_livestemn_xfer_to_litr1n(c)  = m_livestemn_xfer_to_litr1n(c)  + &
                  m_livestemn_xfer_to_litter(p)  * wtcol(p)
               m_deadstemn_xfer_to_litr1n(c)  = m_deadstemn_xfer_to_litr1n(c)  + &
                  m_deadstemn_xfer_to_litter(p)  * wtcol(p)
               m_livecrootn_xfer_to_litr1n(c) = m_livecrootn_xfer_to_litr1n(c) + &
                  m_livecrootn_xfer_to_litter(p) * wtcol(p)
               m_deadcrootn_xfer_to_litr1n(c) = m_deadcrootn_xfer_to_litr1n(c) + &
                  m_deadcrootn_xfer_to_litter(p) * wtcol(p)

            end if
         end if

      end do

   end do

end subroutine CNGapPftToColumn
!-----------------------------------------------------------------------

end module CNGapMortalityMod
