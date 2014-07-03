
module CNGapMortalityMod
#ifdef CN

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
   use pftvarcon       , only: npcropmin
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
! F. Li and S. Levis (11/06/12)
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
#if (defined CNDV)
   real(r8), pointer :: greffic(:)
   real(r8), pointer :: heatstress(:)
   real(r8), pointer :: nind(:)         ! number of individuals (#/m2) added by F. Li and S. Levis
#endif
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
   ivt                            => clm3%g%l%c%p%itype
   leafc                          => clm3%g%l%c%p%pcs%leafc
   frootc                         => clm3%g%l%c%p%pcs%frootc
   livestemc                      => clm3%g%l%c%p%pcs%livestemc
   deadstemc                      => clm3%g%l%c%p%pcs%deadstemc
   livecrootc                     => clm3%g%l%c%p%pcs%livecrootc
   deadcrootc                     => clm3%g%l%c%p%pcs%deadcrootc
   leafc_storage                  => clm3%g%l%c%p%pcs%leafc_storage
   frootc_storage                 => clm3%g%l%c%p%pcs%frootc_storage
   livestemc_storage              => clm3%g%l%c%p%pcs%livestemc_storage
   deadstemc_storage              => clm3%g%l%c%p%pcs%deadstemc_storage
   livecrootc_storage             => clm3%g%l%c%p%pcs%livecrootc_storage
   deadcrootc_storage             => clm3%g%l%c%p%pcs%deadcrootc_storage
   gresp_storage                  => clm3%g%l%c%p%pcs%gresp_storage
   leafc_xfer                     => clm3%g%l%c%p%pcs%leafc_xfer
   frootc_xfer                    => clm3%g%l%c%p%pcs%frootc_xfer
   livestemc_xfer                 => clm3%g%l%c%p%pcs%livestemc_xfer
   deadstemc_xfer                 => clm3%g%l%c%p%pcs%deadstemc_xfer
   livecrootc_xfer                => clm3%g%l%c%p%pcs%livecrootc_xfer
   deadcrootc_xfer                => clm3%g%l%c%p%pcs%deadcrootc_xfer
   gresp_xfer                     => clm3%g%l%c%p%pcs%gresp_xfer
   leafn                          => clm3%g%l%c%p%pns%leafn
   frootn                         => clm3%g%l%c%p%pns%frootn
   livestemn                      => clm3%g%l%c%p%pns%livestemn
   deadstemn                      => clm3%g%l%c%p%pns%deadstemn
   livecrootn                     => clm3%g%l%c%p%pns%livecrootn
   deadcrootn                     => clm3%g%l%c%p%pns%deadcrootn
   retransn                       => clm3%g%l%c%p%pns%retransn
   leafn_storage                  => clm3%g%l%c%p%pns%leafn_storage
   frootn_storage                 => clm3%g%l%c%p%pns%frootn_storage
   livestemn_storage              => clm3%g%l%c%p%pns%livestemn_storage
   deadstemn_storage              => clm3%g%l%c%p%pns%deadstemn_storage
   livecrootn_storage             => clm3%g%l%c%p%pns%livecrootn_storage
   deadcrootn_storage             => clm3%g%l%c%p%pns%deadcrootn_storage
   leafn_xfer                     => clm3%g%l%c%p%pns%leafn_xfer
   frootn_xfer                    => clm3%g%l%c%p%pns%frootn_xfer
   livestemn_xfer                 => clm3%g%l%c%p%pns%livestemn_xfer
   deadstemn_xfer                 => clm3%g%l%c%p%pns%deadstemn_xfer
   livecrootn_xfer                => clm3%g%l%c%p%pns%livecrootn_xfer
   deadcrootn_xfer                => clm3%g%l%c%p%pns%deadcrootn_xfer
   m_leafc_to_litter              => clm3%g%l%c%p%pcf%m_leafc_to_litter
   m_frootc_to_litter             => clm3%g%l%c%p%pcf%m_frootc_to_litter
   m_livestemc_to_litter          => clm3%g%l%c%p%pcf%m_livestemc_to_litter
   m_deadstemc_to_litter          => clm3%g%l%c%p%pcf%m_deadstemc_to_litter
   m_livecrootc_to_litter         => clm3%g%l%c%p%pcf%m_livecrootc_to_litter
   m_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%m_deadcrootc_to_litter
   m_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%m_leafc_storage_to_litter
   m_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%m_frootc_storage_to_litter
   m_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter
   m_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter
   m_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter
   m_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter
   m_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%m_gresp_storage_to_litter
   m_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%m_leafc_xfer_to_litter
   m_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%m_frootc_xfer_to_litter
   m_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%m_livestemc_xfer_to_litter
   m_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_litter
   m_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_litter
   m_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_litter
   m_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%m_gresp_xfer_to_litter
   m_leafn_to_litter              => clm3%g%l%c%p%pnf%m_leafn_to_litter
   m_frootn_to_litter             => clm3%g%l%c%p%pnf%m_frootn_to_litter
   m_livestemn_to_litter          => clm3%g%l%c%p%pnf%m_livestemn_to_litter
   m_deadstemn_to_litter          => clm3%g%l%c%p%pnf%m_deadstemn_to_litter
   m_livecrootn_to_litter         => clm3%g%l%c%p%pnf%m_livecrootn_to_litter
   m_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%m_deadcrootn_to_litter
   m_retransn_to_litter           => clm3%g%l%c%p%pnf%m_retransn_to_litter
   m_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%m_leafn_storage_to_litter
   m_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%m_frootn_storage_to_litter
   m_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter
   m_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter
   m_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter
   m_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter
   m_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter
   m_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter
   m_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter
   m_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter
   m_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter
   m_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter
#if (defined CNDV)
   greffic                        => clm3%g%l%c%p%pdgvs%greffic
   heatstress                     => clm3%g%l%c%p%pdgvs%heatstress
   nind                           => clm3%g%l%c%p%pdgvs%nind     ! F. Li and S. Levis
#endif

   ! set the mortality rate based on annual rate
   am = 0.02_r8

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

#if (defined CNDV)
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
#endif

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
      if (ivt(p) < npcropmin) m_retransn_to_litter(p) = retransn(p) * m

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

! added by F. Li and S. Levis
#if (defined CNDV)
    if (woody(ivt(p)) == 1._r8)then
     if (livestemc(p)+deadstemc(p)> 0._r8)then
         nind(p)=nind(p)*(1._r8-m)
     else
         nind(p) = 0._r8 
     end if
    end if
#endif

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
  use clm_varpar, only : maxpatch_pft, nlevdecomp
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
   logical , pointer :: pactive(:)  ! true=>do computations on this pft (see reweightMod for details)
   integer , pointer :: ivt(:)      ! pft vegetation type
   real(r8), pointer :: wtcol(:)    ! pft weight relative to column (0-1)
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
   real(r8), pointer :: gap_mortality_c_to_litr_met_c(:,:)         ! C fluxes associated with gap mortality to litter metabolic pool (gC/m3/s)
   real(r8), pointer :: gap_mortality_c_to_litr_cel_c(:,:)         ! C fluxes associated with gap mortality to litter cellulose pool (gC/m3/s)
   real(r8), pointer :: gap_mortality_c_to_litr_lig_c(:,:)         ! C fluxes associated with gap mortality to litter lignin pool (gC/m3/s)
   real(r8), pointer :: gap_mortality_c_to_cwdc(:,:)               ! C fluxes associated with gap mortality to CWD pool (gC/m3/s)
   real(r8), pointer :: gap_mortality_n_to_litr_met_n(:,:)         ! N fluxes associated with gap mortality to litter metabolic pool (gN/m3/s)
   real(r8), pointer :: gap_mortality_n_to_litr_cel_n(:,:)         ! N fluxes associated with gap mortality to litter cellulose pool (gN/m3/s)
   real(r8), pointer :: gap_mortality_n_to_litr_lig_n(:,:)         ! N fluxes associated with gap mortality to litter lignin pool (gN/m3/s)
   real(r8), pointer :: gap_mortality_n_to_cwdn(:,:)               ! N fluxes associated with gap mortality to CWD pool (gN/m3/s)
!
! local pointers to implicit in/out arrays
   real(r8), pointer :: leaf_prof(:,:)          ! (1/m) profile of leaves
   real(r8), pointer :: froot_prof(:,:)         ! (1/m) profile of fine roots
   real(r8), pointer :: croot_prof(:,:)         ! (1/m) profile of coarse roots
   real(r8), pointer :: stem_prof(:,:)          ! (1/m) profile of stems

!
! local pointers to implicit out arrays
!
!
! !OTHER LOCAL VARIABLES:
   integer :: fc,c,pi,p,j               ! indices
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
   npfts                          => clm3%g%l%c%npfts
   pfti                           => clm3%g%l%c%pfti

   ! assign local pointers to pft-level arrays
   pactive                        => clm3%g%l%c%p%active
   ivt                            => clm3%g%l%c%p%itype
   wtcol                          => clm3%g%l%c%p%wtcol
   m_leafc_to_litter              => clm3%g%l%c%p%pcf%m_leafc_to_litter
   m_frootc_to_litter             => clm3%g%l%c%p%pcf%m_frootc_to_litter
   m_livestemc_to_litter          => clm3%g%l%c%p%pcf%m_livestemc_to_litter
   m_deadstemc_to_litter          => clm3%g%l%c%p%pcf%m_deadstemc_to_litter
   m_livecrootc_to_litter         => clm3%g%l%c%p%pcf%m_livecrootc_to_litter
   m_deadcrootc_to_litter         => clm3%g%l%c%p%pcf%m_deadcrootc_to_litter
   m_leafc_storage_to_litter      => clm3%g%l%c%p%pcf%m_leafc_storage_to_litter
   m_frootc_storage_to_litter     => clm3%g%l%c%p%pcf%m_frootc_storage_to_litter
   m_livestemc_storage_to_litter  => clm3%g%l%c%p%pcf%m_livestemc_storage_to_litter
   m_deadstemc_storage_to_litter  => clm3%g%l%c%p%pcf%m_deadstemc_storage_to_litter
   m_livecrootc_storage_to_litter => clm3%g%l%c%p%pcf%m_livecrootc_storage_to_litter
   m_deadcrootc_storage_to_litter => clm3%g%l%c%p%pcf%m_deadcrootc_storage_to_litter
   m_gresp_storage_to_litter      => clm3%g%l%c%p%pcf%m_gresp_storage_to_litter
   m_leafc_xfer_to_litter         => clm3%g%l%c%p%pcf%m_leafc_xfer_to_litter
   m_frootc_xfer_to_litter        => clm3%g%l%c%p%pcf%m_frootc_xfer_to_litter
   m_livestemc_xfer_to_litter     => clm3%g%l%c%p%pcf%m_livestemc_xfer_to_litter
   m_deadstemc_xfer_to_litter     => clm3%g%l%c%p%pcf%m_deadstemc_xfer_to_litter
   m_livecrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%m_livecrootc_xfer_to_litter
   m_deadcrootc_xfer_to_litter    => clm3%g%l%c%p%pcf%m_deadcrootc_xfer_to_litter
   m_gresp_xfer_to_litter         => clm3%g%l%c%p%pcf%m_gresp_xfer_to_litter
   m_leafn_to_litter              => clm3%g%l%c%p%pnf%m_leafn_to_litter
   m_frootn_to_litter             => clm3%g%l%c%p%pnf%m_frootn_to_litter
   m_livestemn_to_litter          => clm3%g%l%c%p%pnf%m_livestemn_to_litter
   m_deadstemn_to_litter          => clm3%g%l%c%p%pnf%m_deadstemn_to_litter
   m_livecrootn_to_litter         => clm3%g%l%c%p%pnf%m_livecrootn_to_litter
   m_deadcrootn_to_litter         => clm3%g%l%c%p%pnf%m_deadcrootn_to_litter
   m_retransn_to_litter           => clm3%g%l%c%p%pnf%m_retransn_to_litter
   m_leafn_storage_to_litter      => clm3%g%l%c%p%pnf%m_leafn_storage_to_litter
   m_frootn_storage_to_litter     => clm3%g%l%c%p%pnf%m_frootn_storage_to_litter
   m_livestemn_storage_to_litter  => clm3%g%l%c%p%pnf%m_livestemn_storage_to_litter
   m_deadstemn_storage_to_litter  => clm3%g%l%c%p%pnf%m_deadstemn_storage_to_litter
   m_livecrootn_storage_to_litter => clm3%g%l%c%p%pnf%m_livecrootn_storage_to_litter
   m_deadcrootn_storage_to_litter => clm3%g%l%c%p%pnf%m_deadcrootn_storage_to_litter
   m_leafn_xfer_to_litter         => clm3%g%l%c%p%pnf%m_leafn_xfer_to_litter
   m_frootn_xfer_to_litter        => clm3%g%l%c%p%pnf%m_frootn_xfer_to_litter
   m_livestemn_xfer_to_litter     => clm3%g%l%c%p%pnf%m_livestemn_xfer_to_litter
   m_deadstemn_xfer_to_litter     => clm3%g%l%c%p%pnf%m_deadstemn_xfer_to_litter
   m_livecrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%m_livecrootn_xfer_to_litter
   m_deadcrootn_xfer_to_litter    => clm3%g%l%c%p%pnf%m_deadcrootn_xfer_to_litter
   gap_mortality_c_to_litr_met_c  => clm3%g%l%c%ccf%gap_mortality_c_to_litr_met_c
   gap_mortality_c_to_litr_cel_c  => clm3%g%l%c%ccf%gap_mortality_c_to_litr_cel_c
   gap_mortality_c_to_litr_lig_c  => clm3%g%l%c%ccf%gap_mortality_c_to_litr_lig_c
   gap_mortality_c_to_cwdc        => clm3%g%l%c%ccf%gap_mortality_c_to_cwdc
   gap_mortality_n_to_litr_met_n  => clm3%g%l%c%cnf%gap_mortality_n_to_litr_met_n
   gap_mortality_n_to_litr_cel_n  => clm3%g%l%c%cnf%gap_mortality_n_to_litr_cel_n
   gap_mortality_n_to_litr_lig_n  => clm3%g%l%c%cnf%gap_mortality_n_to_litr_lig_n
   gap_mortality_n_to_cwdn        => clm3%g%l%c%cnf%gap_mortality_n_to_cwdn
   leaf_prof                      => clm3%g%l%c%p%pps%leaf_prof
   froot_prof                     => clm3%g%l%c%p%pps%froot_prof
   croot_prof                     => clm3%g%l%c%p%pps%croot_prof
   stem_prof                      => clm3%g%l%c%p%pps%stem_prof


   do j = 1,nlevdecomp
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
                       (m_livestemc_to_litter(p) + m_deadstemc_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_c_to_cwdc(c,j) = gap_mortality_c_to_cwdc(c,j) + &
                       (m_livecrootc_to_litter(p) + m_deadcrootc_to_litter(p)) * wtcol(p) * croot_prof(p,j)
                  
                  ! storage gap mortality carbon fluxes
                  gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                       (m_leafc_storage_to_litter(p) + m_gresp_storage_to_litter(p))      * wtcol(p) * leaf_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)     = gap_mortality_c_to_litr_met_c(c,j)     + &
                       m_frootc_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                       (m_livestemc_storage_to_litter(p) + m_deadstemc_storage_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                       (m_livecrootc_storage_to_litter(p) + m_deadcrootc_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j)
                  
                  ! transfer gap mortality carbon fluxes
                  gap_mortality_c_to_litr_met_c(c,j)      = gap_mortality_c_to_litr_met_c(c,j)      + &
                       (m_leafc_xfer_to_litter(p) + m_gresp_xfer_to_litter(p))     * wtcol(p) * leaf_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)     = gap_mortality_c_to_litr_met_c(c,j)     + &
                       m_frootc_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j)  = gap_mortality_c_to_litr_met_c(c,j)  + &
                       (m_livestemc_xfer_to_litter(p) + m_deadstemc_xfer_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_c_to_litr_met_c(c,j) = gap_mortality_c_to_litr_met_c(c,j) + &
                       (m_livecrootc_xfer_to_litter(p) + m_deadcrootc_xfer_to_litter(p)) * wtcol(p) * croot_prof(p,j)
                  
                  ! leaf gap mortality nitrogen fluxes
                  gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                       m_leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  gap_mortality_n_to_litr_cel_n(c,j) = gap_mortality_n_to_litr_cel_n(c,j) + &
                       m_leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  gap_mortality_n_to_litr_lig_n(c,j) = gap_mortality_n_to_litr_lig_n(c,j) + &
                       m_leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                  
                  ! fine root litter nitrogen fluxes
                  gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                       m_frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  gap_mortality_n_to_litr_cel_n(c,j) = gap_mortality_n_to_litr_cel_n(c,j) + &
                       m_frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  gap_mortality_n_to_litr_lig_n(c,j) = gap_mortality_n_to_litr_lig_n(c,j) + &
                       m_frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)
                  
                  ! wood gap mortality nitrogen fluxes
                  gap_mortality_n_to_cwdn(c,j)  = gap_mortality_n_to_cwdn(c,j)  + &
                       (m_livestemn_to_litter(p) + m_deadstemn_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_n_to_cwdn(c,j) = gap_mortality_n_to_cwdn(c,j) + &
                       (m_livecrootn_to_litter(p) + m_deadcrootn_to_litter(p)) * wtcol(p) * croot_prof(p,j)
                  
                  ! retranslocated N pool gap mortality fluxes
                  gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                       m_retransn_to_litter(p) * wtcol(p) * leaf_prof(p,j)
                  
                  ! storage gap mortality nitrogen fluxes
                  gap_mortality_n_to_litr_met_n(c,j)      = gap_mortality_n_to_litr_met_n(c,j)      + &
                       m_leafn_storage_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  gap_mortality_n_to_litr_met_n(c,j)     = gap_mortality_n_to_litr_met_n(c,j)     + &
                       m_frootn_storage_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  gap_mortality_n_to_litr_met_n(c,j)  = gap_mortality_n_to_litr_met_n(c,j)  + &
                       (m_livestemn_storage_to_litter(p) + m_deadstemn_storage_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                       (m_livecrootn_storage_to_litter(p) + m_deadcrootn_storage_to_litter(p)) * wtcol(p) * croot_prof(p,j)
                  
                  ! transfer gap mortality nitrogen fluxes
                  gap_mortality_n_to_litr_met_n(c,j)      = gap_mortality_n_to_litr_met_n(c,j)      + &
                       m_leafn_xfer_to_litter(p)      * wtcol(p) * leaf_prof(p,j)
                  gap_mortality_n_to_litr_met_n(c,j)     = gap_mortality_n_to_litr_met_n(c,j)     + &
                       m_frootn_xfer_to_litter(p)     * wtcol(p) * froot_prof(p,j)
                  gap_mortality_n_to_litr_met_n(c,j)  = gap_mortality_n_to_litr_met_n(c,j)  + &
                       (m_livestemn_xfer_to_litter(p) + m_deadstemn_xfer_to_litter(p))  * wtcol(p) * stem_prof(p,j)
                  gap_mortality_n_to_litr_met_n(c,j) = gap_mortality_n_to_litr_met_n(c,j) + &
                       (m_livecrootn_xfer_to_litter(p) + m_deadcrootn_xfer_to_litter(p)) * wtcol(p) * croot_prof(p,j)
                  
                  
               end if
            end if
            
         end do
      end do
   end do

end subroutine CNGapPftToColumn
!-----------------------------------------------------------------------

#endif

end module CNGapMortalityMod
