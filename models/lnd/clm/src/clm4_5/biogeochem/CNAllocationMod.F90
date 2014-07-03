module CNAllocationMod
#ifdef CN

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNAllocationMod
!
! !DESCRIPTION:
! Module holding routines used in allocation model for coupled carbon
! nitrogen code.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varcon, only: dzsoi_decomp
  use abortutils  , only: endrun
  use clm_varctl, only: use_c13, use_c14
  implicit none
  save
  private
! !PUBLIC MEMBER FUNCTIONS:
  public :: CNAllocationInit         ! Initialization
  public :: CNAllocation             ! run method
  public :: CNAllocation_Carbon_only ! Return Carbon_only status

! !PUBLIC DATA MEMBERS:
   character(len=*), parameter, public :: suplnAll=& ! Supplemental Nitrogen for all PFT's
                     'ALL'
   character(len=*), parameter, public :: suplnNon=& ! No supplemental Nitrogen
                     'NONE'
   character(len=15), public :: suplnitro = suplnNon ! Supplemental Nitrogen mode
   logical, public :: Carbon_only = .false.  ! Carbon only mode 
                                             ! (Nitrogen is prescribed NOT prognostic)

! !PRIVATE DATA MEMBERS:
   real(r8):: dt                            !decomp timestep (seconds)
   real(r8):: bdnr                          !bulk denitrification rate (1/s)
   real(r8):: dayscrecover                  !number of days to recover negative cpool
   real(r8), pointer :: arepr(:)            !reproduction allocation coefficient
   real(r8), pointer :: aroot(:)            !root allocation coefficient
   real(r8), pointer:: col_plant_ndemand(:) !column-level plant N demand
!  logical :: crop_supln  = .false.         ! Prognostic crop receives supplemental Nitrogen
!
! !REVISION HISTORY:
! 8/5/03: Created by Peter Thornton
!
!EOP
!-----------------------------------------------------------------------

contains


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNAllocation_Carbon_only
!
! !INTERFACE:
logical function CNAllocation_Carbon_only()
!
! !DESCRIPTION: Return Carbon_only flag.
!
!EOP
!-----------------------------------------------------------------------
  CNAllocation_Carbon_only = Carbon_only
end function CNAllocation_Carbon_only

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNAllocationInit
!
! !INTERFACE:
subroutine CNAllocationInit ( lbc, ubc, lbp, ubp )
!
! !DESCRIPTION:
!
! !USES:
   use clm_varcon      , only: secspday
   use clm_time_manager, only: get_step_size
   use surfrdMod       , only: crop_prog
   use clm_varctl      , only: iulog
   use nanMod          , only: nan
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column-index bounds
   integer, intent(in) :: lbp, ubp        ! pft-index bounds
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! 4/6/11: Created by Erik Kluzek
!
! !LOCAL VARIABLES:
    character(len=32) :: subname = 'CNAllocationInit'
!EOP
!-----------------------------------------------------------------------
   if ( crop_prog )then
      allocate(arepr(lbp:ubp))
      allocate(aroot(lbp:ubp))
      arepr(:) = nan
      aroot(:) = nan
   end if
   allocate(col_plant_ndemand(lbc:ubc))
   col_plant_ndemand(:) = nan

   ! set time steps
   dt = real( get_step_size(), r8 )

   ! set some space-and-time constant parameters 
   bdnr         = 0.5_r8 * (dt/secspday)
   dayscrecover = 30.0_r8

   ! Change namelist settings into private logical variables
   select case(suplnitro)
      case(suplnNon)
         Carbon_only = .false.
      case(suplnAll)
         Carbon_only = .true.
      case default
         write(iulog,*) 'Supplemental Nitrogen flag (suplnitro) can only be: ', &
                        suplnNon, ' or ', suplnAll
         call endrun( trim(subname)//'ERROR: supplemental Nitrogen flag is not correct' )
   end select

end subroutine CNAllocationInit

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNAllocation
!
! !INTERFACE:
subroutine CNAllocation (lbp, ubp, lbc, ubc, &
       num_soilc, filter_soilc, num_soilp, filter_soilp )
!
! !DESCRIPTION:
!
! !USES:
   use clmtype
   use clm_varctl, only: iulog
   use shr_sys_mod, only: shr_sys_flush
   use pft2colMod, only: p2c
   use clm_varpar   , only: nlevsoi, nlevdecomp
#ifdef NITRIF_DENITRIF
   use clm_varcon, only: nitrif_n2o_loss_frac
#endif

   use pftvarcon , only: npcropmin, declfact, bfact, aleaff, arootf, astemf, &
                         arooti, fleafi, allconsl, allconss, grperc, grpnow, &
                         nsoybean
   use clm_varcon, only: secspday, istsoil, istcrop
   use clm_varpar, only: max_pft_per_col
   use clm_time_manager , only : get_step_size
!
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbp, ubp        ! pft-index bounds
   integer, intent(in) :: lbc, ubc        ! column-index bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(ubc-lbc+1) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(ubp-lbp+1) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNdecompAlloc in module CNdecompMod.F90
!
! !REVISION HISTORY:
! 8/5/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: migrated to vector data structures
!
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
   ! pft level
   integer , pointer :: ivt(:)        ! pft vegetation type
   integer , pointer :: pcolumn(:)    ! pft's column index
   integer , pointer :: pfti(:)       ! initial pft index in landunit
   real(r8), pointer :: lgsf(:)       ! long growing season factor [0-1]
   real(r8), pointer :: xsmrpool(:)   ! (gC/m2) temporary photosynthate C pool
   real(r8), pointer :: retransn(:)   ! (gN/m2) plant pool of retranslocated N
   real(r8), pointer :: psnsun(:)     ! sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: psnsha(:)     ! shaded leaf-level photosynthesis (umol CO2 /m**2/ s)

   real(r8), pointer :: c13_psnsun(:) ! C13 sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: c13_psnsha(:) ! C13 shaded leaf-level photosynthesis (umol CO2 /m**2/ s)

   real(r8), pointer :: c14_psnsun(:) ! C14 sunlit leaf-level photosynthesis (umol CO2 /m**2/ s)
   real(r8), pointer :: c14_psnsha(:) ! C14 shaded leaf-level photosynthesis (umol CO2 /m**2/ s)

   real(r8), pointer :: laisun(:)     ! sunlit projected leaf area index
   real(r8), pointer :: laisha(:)     ! shaded projected leaf area index
   real(r8), pointer :: leafc(:)
   real(r8), pointer :: frootc(:)
   real(r8), pointer :: livestemc(:)
   real(r8), pointer :: leaf_mr(:)
   real(r8), pointer :: froot_mr(:)
   real(r8), pointer :: livestem_mr(:)
   real(r8), pointer :: livecroot_mr(:)
   real(r8), pointer :: grain_mr(:)
   real(r8), pointer :: leaf_curmr(:)
   real(r8), pointer :: froot_curmr(:)
   real(r8), pointer :: livestem_curmr(:)
   real(r8), pointer :: livecroot_curmr(:)
   real(r8), pointer :: grain_curmr(:)
   real(r8), pointer :: leaf_xsmr(:)
   real(r8), pointer :: froot_xsmr(:)
   real(r8), pointer :: livestem_xsmr(:)
   real(r8), pointer :: livecroot_xsmr(:)
   real(r8), pointer :: grain_xsmr(:)
   ! column level
   real(r8), pointer :: sminn_vr(:,:)      ! (gN/m3) soil mineral N
   ! ecophysiological constants
   real(r8), pointer :: woody(:)      ! binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: froot_leaf(:) ! allocation parameter: new fine root C per new leaf C (gC/gC)
   real(r8), pointer :: croot_stem(:) ! allocation parameter: new coarse root C per new stem C (gC/gC)
   real(r8), pointer :: stem_leaf(:)  ! allocation parameter: new stem c per new leaf C (gC/gC)
   real(r8), pointer :: flivewd(:)    ! allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
   real(r8), pointer :: leafcn(:)     ! leaf C:N (gC/gN)
   real(r8), pointer :: frootcn(:)    ! fine root C:N (gC/gN)
   real(r8), pointer :: livewdcn(:)   ! live wood (phloem and ray parenchyma) C:N (gC/gN)
   real(r8), pointer :: deadwdcn(:)   ! dead wood (xylem and heartwood) C:N (gC/gN)
   real(r8), pointer :: fcur2(:)      ! allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
   integer, pointer :: plandunit(:)   ! index into landunit level quantities
   integer, pointer :: clandunit(:)   ! index into landunit level quantities
   integer , pointer :: itypelun(:)   ! landunit type
   logical , pointer :: croplive(:)   ! flag, true if planted, not harvested
   integer , pointer :: peaklai(:)    ! 1: max allowed lai; 0: not at max
   real(r8), pointer :: gddmaturity(:)! gdd needed to harvest
   real(r8), pointer :: huileaf(:)    ! heat unit index needed from planting to leaf emergence
   real(r8), pointer :: huigrain(:)   ! same to reach vegetative maturity
   real(r8), pointer :: hui(:)        ! =gdd since planting (gddplant)
   real(r8), pointer :: leafout(:)    ! =gdd from top soil layer temperature
   real(r8), pointer :: aleafi(:)     ! saved allocation coefficient from phase 2
   real(r8), pointer :: astemi(:)     ! saved allocation coefficient from phase 2
   real(r8), pointer :: aleaf(:)      ! leaf allocation coefficient
   real(r8), pointer :: astem(:)      ! stem allocation coefficient
   real(r8), pointer :: graincn(:)    ! grain C:N (gC/gN)
   real(r8), pointer :: fleafcn(:)    ! leaf c:n during organ fill
   real(r8), pointer :: fstemcn(:)    ! stem c:n during organ fill
   real(r8), pointer :: ffrootcn(:)   ! froot c:n during organ fill
!
! local pointers to implicit in/out arrays
!
   ! pft level
   real(r8), pointer :: grain_flag(:)            ! 1: grain fill stage; 0: not
   real(r8), pointer :: gpp(:)                   ! GPP flux before downregulation (gC/m2/s)
   real(r8), pointer :: availc(:)                ! C flux available for allocation (gC/m2/s)
   real(r8), pointer :: xsmrpool_recover(:)         ! C flux assigned to recovery of negative cpool (gC/m2/s)
   real(r8), pointer :: c_allometry(:)           ! C allocation index (DIM)
   real(r8), pointer :: n_allometry(:)           ! N allocation index (DIM)
   real(r8), pointer :: plant_ndemand(:)         ! N flux required to support initial GPP (gN/m2/s)
   real(r8), pointer :: tempsum_potential_gpp(:) ! temporary annual sum of potential GPP 
   real(r8), pointer :: tempmax_retransn(:)      ! temporary annual max of retranslocated N pool (gN/m2)
   real(r8), pointer :: annsum_potential_gpp(:)  ! annual sum of potential GPP
   real(r8), pointer :: avail_retransn(:)        ! N flux available from retranslocation pool (gN/m2/s)
   real(r8), pointer :: annmax_retransn(:)       ! annual max of retranslocated N pool
   real(r8), pointer :: plant_nalloc(:)          ! total allocated N flux (gN/m2/s)
   real(r8), pointer :: plant_calloc(:)          ! total allocated C flux (gC/m2/s)
   real(r8), pointer :: excess_cflux(:)          ! C flux not allocated due to downregulation (gC/m2/s)
   real(r8), pointer :: downreg(:)               ! fractional reduction in GPP due to N limitation (DIM)
   real(r8), pointer :: annsum_npp(:)            ! annual sum of NPP, for wood allocation
   real(r8), pointer :: cpool_to_xsmrpool(:)
   real(r8), pointer :: psnsun_to_cpool(:)
   real(r8), pointer :: psnshade_to_cpool(:)

   real(r8), pointer :: c13_psnsun_to_cpool(:)
   real(r8), pointer :: c13_psnshade_to_cpool(:)

   real(r8), pointer :: c14_psnsun_to_cpool(:)
   real(r8), pointer :: c14_psnshade_to_cpool(:)

   real(r8), pointer :: cpool_to_leafc(:)
   real(r8), pointer :: cpool_to_leafc_storage(:)
   real(r8), pointer :: cpool_to_frootc(:)
   real(r8), pointer :: cpool_to_frootc_storage(:)
   real(r8), pointer :: cpool_to_livestemc(:)
   real(r8), pointer :: cpool_to_livestemc_storage(:)
   real(r8), pointer :: cpool_to_deadstemc(:)
   real(r8), pointer :: cpool_to_deadstemc_storage(:)
   real(r8), pointer :: cpool_to_livecrootc(:)
   real(r8), pointer :: cpool_to_livecrootc_storage(:)
   real(r8), pointer :: cpool_to_deadcrootc(:)
   real(r8), pointer :: cpool_to_deadcrootc_storage(:)
   real(r8), pointer :: cpool_to_gresp_storage(:)         ! allocation to growth respiration storage (gC/m2/s)
   real(r8), pointer :: retransn_to_npool(:)              ! deployment of retranslocated N (gN/m2/s)
   real(r8), pointer :: sminn_to_npool(:)                 ! deployment of soil mineral N uptake (gN/m2/s)
   real(r8), pointer :: cpool_to_grainc(:)                ! allocation to grain C (gC/m2/s)
   real(r8), pointer :: cpool_to_grainc_storage(:)        ! allocation to grain C storage (gC/m2/s)
   real(r8), pointer :: npool_to_grainn(:)                ! allocation to grain N (gN/m2/s)
   real(r8), pointer :: npool_to_grainn_storage(:)        ! allocation to grain N storage (gN/m2/s)
   real(r8), pointer :: npool_to_leafn(:)                 ! allocation to leaf N (gN/m2/s)
   real(r8), pointer :: npool_to_leafn_storage(:)         ! allocation to leaf N storage (gN/m2/s)
   real(r8), pointer :: npool_to_frootn(:)                ! allocation to fine root N (gN/m2/s)
   real(r8), pointer :: npool_to_frootn_storage(:)        ! allocation to fine root N storage (gN/m2/s)
   real(r8), pointer :: npool_to_livestemn(:)
   real(r8), pointer :: npool_to_livestemn_storage(:)
   real(r8), pointer :: npool_to_deadstemn(:)
   real(r8), pointer :: npool_to_deadstemn_storage(:)
   real(r8), pointer :: npool_to_livecrootn(:)
   real(r8), pointer :: npool_to_livecrootn_storage(:)
   real(r8), pointer :: npool_to_deadcrootn(:)
   real(r8), pointer :: npool_to_deadcrootn_storage(:)
   ! column level
   real(r8), pointer :: fpi(:)                          ! fraction of potential immobilization (no units)
   real(r8), pointer :: fpg(:)                          ! fraction of potential gpp (no units)
   real(r8), pointer :: potential_immob(:)
   real(r8), pointer :: actual_immob(:)
   real(r8), pointer :: sminn_to_plant(:)
   real(r8), pointer :: fpi_vr(:,:) ! fraction of potential immobilization (no units)
#ifndef NITRIF_DENITRIF
   real(r8), pointer :: sminn_to_denit_excess_vr(:,:)
#else
   real(r8), parameter :: compet_plant_no3   = 1.0 ! (unitless) relative compettiveness of plants for NO3
   real(r8), parameter :: compet_plant_nh4   = 1.0 ! (unitless) relative compettiveness of plants for NH4
   real(r8), parameter :: compet_decomp_no3  = 1.0 ! (unitless) relative competitiveness of immobilizers for NO3
   real(r8), parameter :: compet_decomp_nh4  = 1.0 ! (unitless) relative competitiveness of immobilizers for NH4
   real(r8), parameter :: compet_denit       = 1.0 ! (unitless) relative competitiveness of denitrifiers for NO3
   real(r8), parameter :: compet_nit         = 1.0 ! (unitless) relative competitiveness of nitrifiers for NH4
   real(r8) :: fpi_no3_vr(lbc:ubc,1:nlevdecomp)               ! fraction of potential immobilization supplied by no3(no units)
   real(r8) :: fpi_nh4_vr(lbc:ubc,1:nlevdecomp)               ! fraction of potential immobilization supplied by nh4 (no units)
   real(r8) :: sum_nh4_demand(lbc:ubc,1:nlevdecomp), sum_nh4_demand_scaled(lbc:ubc,1:nlevdecomp)
   real(r8) :: sum_no3_demand(lbc:ubc,1:nlevdecomp), sum_no3_demand_scaled(lbc:ubc,1:nlevdecomp)
   real(r8), pointer :: smin_no3_vr(:,:)           ! (gN/m3) soil mineral NO3
   real(r8), pointer :: smin_nh4_vr(:,:)           ! (gN/m3) soil mineral NH4
   real(r8), pointer :: f_nit_vr(:,:)              ! (gN/m3/s) soil nitrification flux
   real(r8), pointer :: f_denit_vr(:,:)            ! (gN/m3/s) soil denitrification flux
   real(r8), pointer :: pot_f_nit_vr(:,:)          ! (gN/m3/s) potential soil nitrification flux
   real(r8), pointer :: pot_f_denit_vr(:,:)        ! (gN/m3/s) potential soil denitrification flux
   real(r8), pointer :: actual_immob_no3_vr(:,:)
   real(r8), pointer :: actual_immob_nh4_vr(:,:)
   real(r8), pointer :: smin_no3_to_plant_vr(:,:)
   real(r8), pointer :: smin_nh4_to_plant_vr(:,:)
   real(r8), pointer :: n2_n2o_ratio_denit_vr(:,:)    ! ratio of N2 to N2O production by denitrification [gN/gN]
   real(r8), pointer :: f_n2o_denit_vr(:,:)           ! flux of N2O from denitrification [gN/m3/s]
   real(r8), pointer :: f_n2o_nit_vr(:,:)             ! flux of N2O from nitrification [gN/m3/s]
#endif
   real(r8), pointer :: sminn_to_plant_vr(:,:)
   real(r8), pointer :: supplement_to_sminn_vr(:,:)
   real(r8), pointer :: nfixation_prof(:,:)
   real(r8), pointer :: potential_immob_vr(:,:)
   real(r8), pointer :: actual_immob_vr(:,:)
!
! local pointers to implicit out arrays
!
   real(r8), pointer :: leafn_to_retransn(:)
   real(r8), pointer :: frootn_to_retransn(:)
   real(r8), pointer :: livestemn_to_retransn(:)
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p,l,pi             !indices
   integer :: fp                   !lake filter pft index
   integer :: fc                   !lake filter column index
   real(r8):: mr                   !maintenance respiration (gC/m2/s)
   real(r8):: f1,f2,f3,f4,g1,g2    !allocation parameters
   real(r8):: cnl,cnfr,cnlw,cndw   !C:N ratios for leaf, fine root, and wood
   real(r8):: fcur                 !fraction of current psn displayed as growth
   real(r8):: gresp_storage        !temporary variable for growth resp to storage
   real(r8):: nlc                  !temporary variable for total new leaf carbon allocation
   real(r8):: curmr, curmr_ratio   !xsmrpool temporary variables
   real(r8):: sum_ndemand_vr(lbc:ubc, 1:nlevdecomp)          !total column N demand (gN/m3/s) at a given level
   integer :: j                    ! level index
   real(r8):: nuptake_prof(lbc:ubc, 1:nlevdecomp)
   real(r8):: sminn_tot(lbc:ubc)
   real(r8) f5                     !grain allocation parameter
   real(r8) cng                    !C:N ratio for grain (= cnlw for now; slevis)
   real(r8) fleaf                  !fraction allocated to leaf
   real(r8) t1                     !temporary variable

#ifndef NITRIF_DENITRIF
   integer :: nlimit(lbc:ubc,0:nlevdecomp)               !flag for N limitation
   real(r8):: residual_sminn_vr(lbc:ubc, 1:nlevdecomp)
   real(r8):: residual_sminn(lbc:ubc)
#else
   integer :: nlimit_no3(lbc:ubc,0:nlevdecomp)               !flag for NO3 limitation
   integer :: nlimit_nh4(lbc:ubc,0:nlevdecomp)               !flag for NH4 limitation
   real(r8):: residual_smin_nh4_vr(lbc:ubc, 1:nlevdecomp)
   real(r8):: residual_smin_no3_vr(lbc:ubc, 1:nlevdecomp)
   real(r8):: residual_smin_nh4(lbc:ubc)
   real(r8):: residual_smin_no3(lbc:ubc)
#endif
   real(r8):: residual_plant_ndemand(lbc:ubc)


!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
   ivt                         => clm3%g%l%c%p%itype
   pcolumn                     => clm3%g%l%c%p%column
   plandunit                   => clm3%g%l%c%p%landunit
   clandunit                   => clm3%g%l%c%landunit
   pfti                        => clm3%g%l%c%pfti
   itypelun                    => clm3%g%l%itype
   lgsf                        => clm3%g%l%c%p%pepv%lgsf
   xsmrpool                    => clm3%g%l%c%p%pcs%xsmrpool
   retransn                    => clm3%g%l%c%p%pns%retransn
   psnsun                      => clm3%g%l%c%p%pcf%psnsun
   psnsha                      => clm3%g%l%c%p%pcf%psnsha
   if ( use_c13 ) then
      c13_psnsun                  => clm3%g%l%c%p%pc13f%psnsun
      c13_psnsha                  => clm3%g%l%c%p%pc13f%psnsha
      c13_psnsun_to_cpool         => clm3%g%l%c%p%pc13f%psnsun_to_cpool
      c13_psnshade_to_cpool       => clm3%g%l%c%p%pc13f%psnshade_to_cpool
   endif
   if ( use_c14 ) then
      c14_psnsun                  => clm3%g%l%c%p%pc14f%psnsun
      c14_psnsha                  => clm3%g%l%c%p%pc14f%psnsha
      c14_psnsun_to_cpool         => clm3%g%l%c%p%pc14f%psnsun_to_cpool
      c14_psnshade_to_cpool       => clm3%g%l%c%p%pc14f%psnshade_to_cpool
   endif
   laisun                      => clm3%g%l%c%p%pps%laisun
   laisha                      => clm3%g%l%c%p%pps%laisha
   leafc                       => clm3%g%l%c%p%pcs%leafc
   frootc                      => clm3%g%l%c%p%pcs%frootc
   livestemc                   => clm3%g%l%c%p%pcs%livestemc
   leaf_mr                     => clm3%g%l%c%p%pcf%leaf_mr
   froot_mr                    => clm3%g%l%c%p%pcf%froot_mr
   livestem_mr                 => clm3%g%l%c%p%pcf%livestem_mr
   livecroot_mr                => clm3%g%l%c%p%pcf%livecroot_mr
   grain_mr                    => clm3%g%l%c%p%pcf%grain_mr
   leaf_curmr                  => clm3%g%l%c%p%pcf%leaf_curmr
   froot_curmr                 => clm3%g%l%c%p%pcf%froot_curmr
   livestem_curmr              => clm3%g%l%c%p%pcf%livestem_curmr
   livecroot_curmr             => clm3%g%l%c%p%pcf%livecroot_curmr
   grain_curmr                 => clm3%g%l%c%p%pcf%grain_curmr
   leaf_xsmr                   => clm3%g%l%c%p%pcf%leaf_xsmr
   froot_xsmr                  => clm3%g%l%c%p%pcf%froot_xsmr
   livestem_xsmr               => clm3%g%l%c%p%pcf%livestem_xsmr
   livecroot_xsmr              => clm3%g%l%c%p%pcf%livecroot_xsmr
   grain_xsmr                  => clm3%g%l%c%p%pcf%grain_xsmr
   sminn_vr                       => clm3%g%l%c%cns%sminn_vr
   woody                       => pftcon%woody
   froot_leaf                  => pftcon%froot_leaf
   croot_stem                  => pftcon%croot_stem
   stem_leaf                   => pftcon%stem_leaf
   flivewd                     => pftcon%flivewd
   leafcn                      => pftcon%leafcn
   frootcn                     => pftcon%frootcn
   livewdcn                    => pftcon%livewdcn
   deadwdcn                    => pftcon%deadwdcn
   fcur2                       => pftcon%fcur
   gddmaturity                 => clm3%g%l%c%p%pps%gddmaturity
   huileaf                     => clm3%g%l%c%p%pps%huileaf
   huigrain                    => clm3%g%l%c%p%pps%huigrain
   hui                         => clm3%g%l%c%p%pps%gddplant
   leafout                     => clm3%g%l%c%p%pps%gddtsoi
   croplive                    => clm3%g%l%c%p%pps%croplive
   peaklai                     => clm3%g%l%c%p%pps%peaklai
   graincn                     => pftcon%graincn
   fleafcn                     => pftcon%fleafcn
   ffrootcn                    => pftcon%ffrootcn
   fstemcn                     => pftcon%fstemcn
   ! Assign local pointers to derived type arrays (out)
   grain_flag                  => clm3%g%l%c%p%pepv%grain_flag
   gpp                         => clm3%g%l%c%p%pepv%gpp
   availc                      => clm3%g%l%c%p%pepv%availc
   xsmrpool_recover            => clm3%g%l%c%p%pepv%xsmrpool_recover
   c_allometry                 => clm3%g%l%c%p%pepv%c_allometry
   n_allometry                 => clm3%g%l%c%p%pepv%n_allometry
   plant_ndemand               => clm3%g%l%c%p%pepv%plant_ndemand
   tempsum_potential_gpp       => clm3%g%l%c%p%pepv%tempsum_potential_gpp
   tempmax_retransn            => clm3%g%l%c%p%pepv%tempmax_retransn
   annsum_potential_gpp        => clm3%g%l%c%p%pepv%annsum_potential_gpp
   avail_retransn              => clm3%g%l%c%p%pepv%avail_retransn
   annmax_retransn             => clm3%g%l%c%p%pepv%annmax_retransn
   plant_nalloc                => clm3%g%l%c%p%pepv%plant_nalloc
   plant_calloc                => clm3%g%l%c%p%pepv%plant_calloc
   excess_cflux                => clm3%g%l%c%p%pepv%excess_cflux
   downreg                     => clm3%g%l%c%p%pepv%downreg
   annsum_npp                  => clm3%g%l%c%p%pepv%annsum_npp
   cpool_to_xsmrpool           => clm3%g%l%c%p%pcf%cpool_to_xsmrpool
   psnsun_to_cpool             => clm3%g%l%c%p%pcf%psnsun_to_cpool
   psnshade_to_cpool           => clm3%g%l%c%p%pcf%psnshade_to_cpool
   cpool_to_leafc              => clm3%g%l%c%p%pcf%cpool_to_leafc
   cpool_to_leafc_storage      => clm3%g%l%c%p%pcf%cpool_to_leafc_storage
   cpool_to_frootc             => clm3%g%l%c%p%pcf%cpool_to_frootc
   cpool_to_frootc_storage     => clm3%g%l%c%p%pcf%cpool_to_frootc_storage
   cpool_to_livestemc          => clm3%g%l%c%p%pcf%cpool_to_livestemc
   cpool_to_livestemc_storage  => clm3%g%l%c%p%pcf%cpool_to_livestemc_storage
   cpool_to_deadstemc          => clm3%g%l%c%p%pcf%cpool_to_deadstemc
   cpool_to_deadstemc_storage  => clm3%g%l%c%p%pcf%cpool_to_deadstemc_storage
   cpool_to_livecrootc         => clm3%g%l%c%p%pcf%cpool_to_livecrootc
   cpool_to_livecrootc_storage => clm3%g%l%c%p%pcf%cpool_to_livecrootc_storage
   cpool_to_deadcrootc         => clm3%g%l%c%p%pcf%cpool_to_deadcrootc
   cpool_to_deadcrootc_storage => clm3%g%l%c%p%pcf%cpool_to_deadcrootc_storage
   cpool_to_gresp_storage      => clm3%g%l%c%p%pcf%cpool_to_gresp_storage
   cpool_to_grainc             => clm3%g%l%c%p%pcf%cpool_to_grainc
   cpool_to_grainc_storage     => clm3%g%l%c%p%pcf%cpool_to_grainc_storage
   npool_to_grainn             => clm3%g%l%c%p%pnf%npool_to_grainn
   npool_to_grainn_storage     => clm3%g%l%c%p%pnf%npool_to_grainn_storage
   retransn_to_npool           => clm3%g%l%c%p%pnf%retransn_to_npool
   sminn_to_npool              => clm3%g%l%c%p%pnf%sminn_to_npool
   npool_to_leafn              => clm3%g%l%c%p%pnf%npool_to_leafn
   npool_to_leafn_storage      => clm3%g%l%c%p%pnf%npool_to_leafn_storage
   npool_to_frootn             => clm3%g%l%c%p%pnf%npool_to_frootn
   npool_to_frootn_storage     => clm3%g%l%c%p%pnf%npool_to_frootn_storage
   npool_to_livestemn          => clm3%g%l%c%p%pnf%npool_to_livestemn
   npool_to_livestemn_storage  => clm3%g%l%c%p%pnf%npool_to_livestemn_storage
   npool_to_deadstemn          => clm3%g%l%c%p%pnf%npool_to_deadstemn
   npool_to_deadstemn_storage  => clm3%g%l%c%p%pnf%npool_to_deadstemn_storage
   npool_to_livecrootn         => clm3%g%l%c%p%pnf%npool_to_livecrootn
   npool_to_livecrootn_storage => clm3%g%l%c%p%pnf%npool_to_livecrootn_storage
   npool_to_deadcrootn         => clm3%g%l%c%p%pnf%npool_to_deadcrootn
   npool_to_deadcrootn_storage => clm3%g%l%c%p%pnf%npool_to_deadcrootn_storage
   leafn_to_retransn           => clm3%g%l%c%p%pnf%leafn_to_retransn
   frootn_to_retransn          => clm3%g%l%c%p%pnf%frootn_to_retransn
   livestemn_to_retransn       => clm3%g%l%c%p%pnf%livestemn_to_retransn
   fpg                         => clm3%g%l%c%cps%fpg
   potential_immob             => clm3%g%l%c%cnf%potential_immob
   actual_immob                => clm3%g%l%c%cnf%actual_immob
   sminn_to_plant              => clm3%g%l%c%cnf%sminn_to_plant
   fpi                         => clm3%g%l%c%cps%fpi
   fpi_vr                         => clm3%g%l%c%cps%fpi_vr
#ifndef NITRIF_DENITRIF
   sminn_to_denit_excess_vr       => clm3%g%l%c%cnf%sminn_to_denit_excess_vr
#else
   smin_nh4_vr                  => clm3%g%l%c%cns%smin_nh4_vr
   smin_no3_vr                  => clm3%g%l%c%cns%smin_no3_vr
   pot_f_nit_vr                 => clm3%g%l%c%cnf%pot_f_nit_vr
   pot_f_denit_vr               => clm3%g%l%c%cnf%pot_f_denit_vr
   f_nit_vr                     => clm3%g%l%c%cnf%f_nit_vr
   f_denit_vr                   => clm3%g%l%c%cnf%f_denit_vr
   actual_immob_no3_vr          => clm3%g%l%c%cnf%actual_immob_no3_vr
   actual_immob_nh4_vr          => clm3%g%l%c%cnf%actual_immob_nh4_vr
   smin_no3_to_plant_vr         => clm3%g%l%c%cnf%smin_no3_to_plant_vr
   smin_nh4_to_plant_vr         => clm3%g%l%c%cnf%smin_nh4_to_plant_vr
   n2_n2o_ratio_denit_vr        => clm3%g%l%c%cnf%n2_n2o_ratio_denit_vr
   f_n2o_denit_vr               => clm3%g%l%c%cnf%f_n2o_denit_vr
   f_n2o_nit_vr                 => clm3%g%l%c%cnf%f_n2o_nit_vr
#endif
   supplement_to_sminn_vr         => clm3%g%l%c%cnf%supplement_to_sminn_vr
   sminn_to_plant_vr              => clm3%g%l%c%cnf%sminn_to_plant_vr
   nfixation_prof                 => clm3%g%l%c%cps%nfixation_prof
   potential_immob_vr             => clm3%g%l%c%cnf%potential_immob_vr
   actual_immob_vr                => clm3%g%l%c%cnf%actual_immob_vr
   aleafi                      => clm3%g%l%c%p%pps%aleafi
   astemi                      => clm3%g%l%c%p%pps%astemi
   aleaf                       => clm3%g%l%c%p%pps%aleaf
   astem                       => clm3%g%l%c%p%pps%astem


   ! set time steps
   dt = real( get_step_size(), r8 )

   ! loop over pfts to assess the total plant N demand
   do fp=1,num_soilp
      p = filter_soilp(fp)

      ! get the time step total gross photosynthesis
      ! this is coming from the canopy fluxes code, and is the
      ! gpp that is used to control stomatal conductance.
      ! For the nitrogen downregulation code, this is assumed
      ! to be the potential gpp, and the actual gpp will be
      ! reduced due to N limitation. 
      
      ! Convert psn from umol/m2/s -> gC/m2/s

      ! The input psn (psnsun and psnsha) are expressed per unit LAI
      ! in the sunlit and shaded canopy, respectively. These need to be
      ! scaled by laisun and laisha to get the total gpp for allocation

      psnsun_to_cpool(p) = psnsun(p) * laisun(p) * 12.011e-6_r8
      psnshade_to_cpool(p) = psnsha(p) * laisha(p) * 12.011e-6_r8

      if ( use_c13 ) then
         c13_psnsun_to_cpool(p) = c13_psnsun(p) * laisun(p) * 12.011e-6_r8
         c13_psnshade_to_cpool(p) = c13_psnsha(p) * laisha(p) * 12.011e-6_r8
      endif

      if ( use_c14 ) then
         c14_psnsun_to_cpool(p) = c14_psnsun(p) * laisun(p) * 12.011e-6_r8
         c14_psnshade_to_cpool(p) = c14_psnsha(p) * laisha(p) * 12.011e-6_r8
      endif
      
      gpp(p) = psnsun_to_cpool(p) + psnshade_to_cpool(p)

      ! get the time step total maintenance respiration
      ! These fluxes should already be in gC/m2/s

      mr = leaf_mr(p) + froot_mr(p)
      if (woody(ivt(p)) == 1.0_r8) then
         mr = mr + livestem_mr(p) + livecroot_mr(p)
      else if (ivt(p) >= npcropmin) then
         if (croplive(p)) mr = mr + livestem_mr(p) + grain_mr(p)
      end if

      ! carbon flux available for allocation
      availc(p) = gpp(p) - mr
      
      ! new code added for isotope calculations, 7/1/05, PET
      ! If mr > gpp, then some mr comes from gpp, the rest comes from
      ! cpool (xsmr)
      if (mr > 0._r8 .and. availc(p) < 0._r8) then
         curmr = gpp(p)
         curmr_ratio = curmr / mr
      else
         curmr_ratio = 1._r8
      end if
      leaf_curmr(p) = leaf_mr(p) * curmr_ratio
      leaf_xsmr(p) = leaf_mr(p) - leaf_curmr(p)
      froot_curmr(p) = froot_mr(p) * curmr_ratio
      froot_xsmr(p) = froot_mr(p) - froot_curmr(p)
      livestem_curmr(p) = livestem_mr(p) * curmr_ratio
      livestem_xsmr(p) = livestem_mr(p) - livestem_curmr(p)
      livecroot_curmr(p) = livecroot_mr(p) * curmr_ratio
      livecroot_xsmr(p) = livecroot_mr(p) - livecroot_curmr(p)
      grain_curmr(p) = grain_mr(p) * curmr_ratio
      grain_xsmr(p) = grain_mr(p) - grain_curmr(p)
      
      ! no allocation when available c is negative
      availc(p) = max(availc(p),0.0_r8)

      ! test for an xsmrpool deficit
      if (xsmrpool(p) < 0.0_r8) then
         ! Running a deficit in the xsmrpool, so the first priority is to let
         ! some availc from this timestep accumulate in xsmrpool.
         ! Determine rate of recovery for xsmrpool deficit

         xsmrpool_recover(p) = -xsmrpool(p)/(dayscrecover*secspday)
         if (xsmrpool_recover(p) < availc(p)) then
             ! available carbon reduced by amount for xsmrpool recovery
             availc(p) = availc(p) - xsmrpool_recover(p)
         else
             ! all of the available carbon goes to xsmrpool recovery
             xsmrpool_recover(p) = availc(p)
             availc(p) = 0.0_r8
         end if
         cpool_to_xsmrpool(p) = xsmrpool_recover(p)
      end if

      f1 = froot_leaf(ivt(p))
      f2 = croot_stem(ivt(p))
     
      ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
      ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
      ! This variable allocation is only for trees. Shrubs have a constant
      ! allocation as specified in the pft-physiology file.  The value is also used
      ! as a trigger here: -1.0 means to use the dynamic allocation (trees).
      if (stem_leaf(ivt(p)) == -1._r8) then
         f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
      else
         f3 = stem_leaf(ivt(p))
      end if
      
      f4 = flivewd(ivt(p))
      g1 = grperc(ivt(p))
      g2 = grpnow(ivt(p))
      cnl = leafcn(ivt(p))
      cnfr = frootcn(ivt(p))
      cnlw = livewdcn(ivt(p))
      cndw = deadwdcn(ivt(p))

      ! calculate f1 to f5 for prog crops following AgroIBIS subr phenocrop

      f5 = 0._r8 ! continued intializations from above

      if (ivt(p) >= npcropmin) then ! skip 2 generic crops

         if (croplive(p)) then
            ! same phases appear in subroutine CropPhenology

            ! Phase 1 completed:
            ! ==================
            ! if hui is less than the number of gdd needed for filling of grain
            ! leaf emergence also has to have taken place for lai changes to occur
            ! and carbon assimilation
            ! Next phase: leaf emergence to start of leaf decline

            if (leafout(p) >= huileaf(p) .and. hui(p) < huigrain(p)) then

               ! allocation rules for crops based on maturity and linear decrease
               ! of amount allocated to roots over course of the growing season
   
               if (peaklai(p) == 1) then ! lai at maximum allowed
                  arepr(p) = 0._r8
                  aleaf(p) = 1.e-5_r8
                  astem(p) = 0._r8
                  aroot(p) = 1._r8 - arepr(p) - aleaf(p) - astem(p)
               else
                  arepr(p) = 0._r8
                  aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) -   &
                                 (arooti(ivt(p)) - arootf(ivt(p))) *  &
                                 min(1._r8, hui(p)/gddmaturity(p))))
                  fleaf = fleafi(ivt(p)) * (exp(-bfact(ivt(p))) -         &
                                exp(-bfact(ivt(p))*hui(p)/huigrain(p))) / &
                                (exp(-bfact(ivt(p)))-1) ! fraction alloc to leaf (from J Norman alloc curve)
                  aleaf(p) = max(1.e-5_r8, (1._r8 - aroot(p)) * fleaf)
                  astem(p) = 1._r8 - arepr(p) - aleaf(p) - aroot(p)
               end if
   
               ! AgroIBIS included here an immediate adjustment to aleaf & astem if the 
               ! predicted lai from the above allocation coefficients exceeded laimx.
               ! We have decided to live with lais slightly higher than laimx by
               ! enforcing the cap in the following tstep through the peaklai logic above.

               astemi(p) = astem(p) ! save for use by equations after shift
               aleafi(p) = aleaf(p) ! to reproductive phenology stage begins
               grain_flag(p) = 0._r8 ! setting to 0 while in phase 2

               ! Phase 2 completed:
               ! ==================
               ! shift allocation either when enough gdd are accumulated or maximum number
               ! of days has elapsed since planting

            else if (hui(p) >= huigrain(p)) then
   
               aroot(p) = max(0._r8, min(1._r8, arooti(ivt(p)) - &
                         (arooti(ivt(p)) - arootf(ivt(p))) * min(1._r8, hui(p)/gddmaturity(p))))
               if (astemi(p) > astemf(ivt(p))) then
                  astem(p) = max(0._r8, max(astemf(ivt(p)), astem(p) * &
                                 (1._r8 - min((hui(p)-                 &
                                 huigrain(p))/((gddmaturity(p)*declfact(ivt(p)))- &
                                 huigrain(p)),1._r8)**allconss(ivt(p)) )))
               end if
               if (aleafi(p) > aleaff(ivt(p))) then
                  aleaf(p) = max(1.e-5_r8, max(aleaff(ivt(p)), aleaf(p) * &
                                 (1._r8 - min((hui(p)-                    &
                                 huigrain(p))/((gddmaturity(p)*declfact(ivt(p)))- &
                                 huigrain(p)),1._r8)**allconsl(ivt(p)) )))
               end if

               !Beth's retranslocation of leafn, stemn, rootn to organ
               !Filter excess plant N to retransn pool for organ N
               !Only do one time then hold grain_flag till onset next season

               ! slevis: Will astem ever = astemf exactly?
               ! Beth's response: ...looks like astem can equal astemf under the right circumstances. 
               !It might be worth a rewrite to capture what I was trying to do, but the retranslocation for 
               !corn and wheat begins at the beginning of the grain fill stage, but for soybean I was holding it 
               !until after the leaf and stem decline were complete. Looking at how astem is calculated, once the 
               !stem decline is near complete, astem should (usually) be set to astemf. The reason for holding off 
               !on soybean is that the retranslocation scheme begins at the beginning of the grain phase, when the 
               !leaf and stem are still growing, but declining. Since carbon is still getting allocated and now 
               !there is more nitrogen available, the nitrogen can be diverted from grain. For corn and wheat 
               !the impact was probably enough to boost productivity, but for soybean the nitrogen was better off 
               !fulfilling the grain fill. It seems that if the peak lai is reached for soybean though that this 
               !would be bypassed altogether, not the intended outcome. I checked several of my output files and 
               !they all seemed to be going through the retranslocation loop for soybean - good news.

               if (ivt(p) /= nsoybean .or. astem(p) == astemf(ivt(p))) then
                  if (grain_flag(p) == 0._r8) then
                     t1 = 1 / dt
                     leafn_to_retransn(p) = t1 * ((leafc(p) / leafcn(ivt(p))) - (leafc(p) / fleafcn(ivt(p))))
                     livestemn_to_retransn(p) = t1 * ((livestemc(p) / livewdcn(ivt(p))) - (livestemc(p) / fstemcn(ivt(p))))
                     frootn_to_retransn(p) = 0._r8
                     if (ffrootcn(ivt(p)) > 0._r8) then
                        frootn_to_retransn(p) = t1 * ((frootc(p) / frootcn(ivt(p))) - (frootc(p) / ffrootcn(ivt(p))))
                     end if
                     grain_flag(p) = 1._r8
                  end if
               end if

               arepr(p) = 1._r8 - aroot(p) - astem(p) - aleaf(p)

            else                   ! pre emergence
               aleaf(p) = 1.e-5_r8 ! allocation coefficients should be irrelevant
               astem(p) = 0._r8    ! because crops have no live carbon pools;
               aroot(p) = 0._r8    ! this applies to this "else" and to the "else"
               arepr(p) = 0._r8    ! a few lines down
            end if

            f1 = aroot(p) / aleaf(p)
            f3 = astem(p) / aleaf(p)
            f5 = arepr(p) / aleaf(p)
            g1 = 0.25_r8

         else   ! .not croplive
            f1 = 0._r8
            f3 = 0._r8
            f5 = 0._r8
            g1 = 0.25_r8
         end if
      end if

      ! based on available C, use constant allometric relationships to
      ! determine N requirements
      if (woody(ivt(p)) == 1.0_r8) then
         c_allometry(p) = (1._r8+g1)*(1._r8+f1+f3*(1._r8+f2))
         n_allometry(p) = 1._r8/cnl + f1/cnfr + (f3*f4*(1._r8+f2))/cnlw + &
                       (f3*(1._r8-f4)*(1._r8+f2))/cndw
      else if (ivt(p) >= npcropmin) then ! skip generic crops
         cng = graincn(ivt(p))
         c_allometry(p) = (1._r8+g1)*(1._r8+f1+f5+f3*(1._r8+f2))
         n_allometry(p) = 1._r8/cnl + f1/cnfr + f5/cng + (f3*f4*(1._r8+f2))/cnlw + &
                       (f3*(1._r8-f4)*(1._r8+f2))/cndw
      else
         c_allometry(p) = 1._r8+g1+f1+f1*g1
         n_allometry(p) = 1._r8/cnl + f1/cnfr
      end if
      plant_ndemand(p) = availc(p)*(n_allometry(p)/c_allometry(p))

      ! retranslocated N deployment depends on seasonal cycle of potential GPP
      ! (requires one year run to accumulate demand)

      tempsum_potential_gpp(p) = tempsum_potential_gpp(p) + gpp(p)

      ! Adding the following line to carry max retransn info to CN Annual Update
      tempmax_retransn(p) = max(tempmax_retransn(p),retransn(p))

      ! Beth's code: crops pull from retransn pool only during grain fill;
      !              retransn pool has N from leaves, stems, and roots for
      !              retranslocation

      if (ivt(p) >= npcropmin .and. grain_flag(p) == 1._r8) then
         avail_retransn(p) = plant_ndemand(p)
      else if (ivt(p) < npcropmin .and. annsum_potential_gpp(p) > 0._r8) then
         avail_retransn(p) = (annmax_retransn(p)/2._r8)*(gpp(p)/annsum_potential_gpp(p))/dt
      else
         avail_retransn(p) = 0.0_r8
      end if

      ! make sure available retrans N doesn't exceed storage
      avail_retransn(p) = min(avail_retransn(p), retransn(p)/dt)

      ! modify plant N demand according to the availability of
      ! retranslocated N
      ! take from retransn pool at most the flux required to meet
      ! plant ndemand

      if (plant_ndemand(p) > avail_retransn(p)) then
         retransn_to_npool(p) = avail_retransn(p)
      else
         retransn_to_npool(p) = plant_ndemand(p)
      end if
      plant_ndemand(p) = plant_ndemand(p) - retransn_to_npool(p)

   end do ! end pft loop

   ! now use the p2c routine to get the column-averaged plant_ndemand
   call p2c(num_soilc,filter_soilc,plant_ndemand,col_plant_ndemand)

#ifndef NITRIF_DENITRIF
   ! column loops to resolve plant/heterotroph competition for mineral N

   ! init sminn_tot
   do fc=1,num_soilc
      c = filter_soilc(fc)
      sminn_tot(c) = 0.
   end do

   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)
         sminn_tot(c) = sminn_tot(c) + sminn_vr(c,j) * dzsoi_decomp(j)
      end do
   end do
   
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)      
         if (sminn_tot(c) .gt. 0.) then
               nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
         else
               nuptake_prof(c,j) = nfixation_prof(c,j)
         endif
         
         sum_ndemand_vr(c,j) = col_plant_ndemand(c) * nuptake_prof(c,j) + potential_immob_vr(c,j)
      end do
   end do

   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)      
         l = clandunit(c)
         if (sum_ndemand_vr(c,j)*dt < sminn_vr(c,j)) then
            
            ! N availability is not limiting immobilization or plant
            ! uptake, and both can proceed at their potential rates
            nlimit(c,j) = 0
            fpi_vr(c,j) = 1.0_r8
            actual_immob_vr(c,j) = potential_immob_vr(c,j)
            sminn_to_plant_vr(c,j) = col_plant_ndemand(c) * nuptake_prof(c,j)
         else if ( Carbon_only) then !.or. &
!                (crop_supln .and. (itypelun(l) == istcrop) .and. &
!                (ivt(pfti(c)) >= npcropmin)) )then
            ! this code block controls the addition of N to sminn pool
            ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
            ! model behave essentially as a carbon-only model, but with the
            ! benefit of keeping track of the N additions needed to
            ! eliminate N limitations, so there is still a diagnostic quantity
            ! that describes the degree of N limitation at steady-state.

             nlimit(c,j) = 1
             fpi_vr(c,j) = 1.0_r8
             actual_immob_vr(c,j) = potential_immob_vr(c,j)
             sminn_to_plant_vr(c,j) =  col_plant_ndemand(c) * nuptake_prof(c,j)
             supplement_to_sminn_vr(c,j) = sum_ndemand_vr(c,j) - (sminn_vr(c,j)/dt)
         else
            ! N availability can not satisfy the sum of immobilization and
            ! plant growth demands, so these two demands compete for available
            ! soil mineral N resource.
            
            nlimit(c,j) = 1
            if (sum_ndemand_vr(c,j) > 0.0_r8) then
               actual_immob_vr(c,j) = (sminn_vr(c,j)/dt)*(potential_immob_vr(c,j) / sum_ndemand_vr(c,j))
            else
               actual_immob_vr(c,j) = 0.0_r8
            end if
            
            if (potential_immob_vr(c,j) > 0.0_r8) then
               fpi_vr(c,j) = actual_immob_vr(c,j) / potential_immob_vr(c,j)
            else
               fpi_vr(c,j) = 0.0_r8
            end if
            
            sminn_to_plant_vr(c,j) = (sminn_vr(c,j)/dt) - actual_immob_vr(c,j)
         end if
      end do
   end do

   ! sum up N fluxes to plant
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)    
         sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
      end do
   end do

   ! give plants a second pass to see if there is any mineral N left over with which to satisfy residual N demand.
   do fc=1,num_soilc
      c = filter_soilc(fc)    
      residual_sminn(c) = 0._r8
   end do
   
   ! sum up total N left over after initial plant and immobilization fluxes
   do fc=1,num_soilc
      c = filter_soilc(fc)    
      residual_plant_ndemand(c) = col_plant_ndemand(c) - sminn_to_plant(c)
   end do
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)    
         if (residual_plant_ndemand(c) .gt. 0._r8 ) then
            if (nlimit(c,j) .eq. 0) then
               residual_sminn_vr(c,j) = max(sminn_vr(c,j) - (actual_immob_vr(c,j) + sminn_to_plant_vr(c,j) ) * dt, 0._r8)
               residual_sminn(c) = residual_sminn(c) + residual_sminn_vr(c,j) * dzsoi_decomp(j)
            else
               residual_sminn_vr(c,j)  = 0._r8
            endif
         endif
      end do
   end do

   ! distribute residual N to plants
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)    
         if ( residual_plant_ndemand(c) .gt. 0._r8 .and. residual_sminn(c) .gt. 0._r8 .and. nlimit(c,j) .eq. 0) then
               sminn_to_plant_vr(c,j) = sminn_to_plant_vr(c,j) + residual_sminn_vr(c,j) * &
                    min(( residual_plant_ndemand(c) *  dt ) / residual_sminn(c), 1._r8) / dt
         endif
      end do
   end do

   ! re-sum up N fluxes to plant
   do fc=1,num_soilc
      c = filter_soilc(fc)    
      sminn_to_plant(c) = 0._r8
   end do
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)    
         sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
         sum_ndemand_vr(c,j) = potential_immob_vr(c,j) + sminn_to_plant_vr(c,j)
      end do
   end do

   ! under conditions of excess N, some proportion is assumed to
   ! be lost to denitrification, in addition to the constant
   ! proportion lost in the decomposition pathways
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)    
         if ((sminn_to_plant_vr(c,j) + actual_immob_vr(c,j))*dt < sminn_vr(c,j)) then
            sminn_to_denit_excess_vr(c,j) = max(bdnr*((sminn_vr(c,j)/dt) - sum_ndemand_vr(c,j)),0._r8)
         else
            sminn_to_denit_excess_vr(c,j) = 0._r8
         endif
      end do
   end do

   ! sum up N fluxes to immobilization
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)    
         actual_immob(c) = actual_immob(c) + actual_immob_vr(c,j) * dzsoi_decomp(j)
         potential_immob(c) = potential_immob(c) + potential_immob_vr(c,j) * dzsoi_decomp(j)
      end do
   end do

   do fc=1,num_soilc
      c = filter_soilc(fc)    
      ! calculate the fraction of potential growth that can be
      ! acheived with the N available to plants      
      if (col_plant_ndemand(c) > 0.0_r8) then
         fpg(c) = sminn_to_plant(c) / col_plant_ndemand(c)
      else
         fpg(c) = 1.0_r8
      end if
      
      ! calculate the fraction of immobilization realized (for diagnostic purposes)
      if (potential_immob(c) > 0.0_r8) then
         fpi(c) = actual_immob(c) / potential_immob(c)
      else
         fpi(c) = 1.0_r8
      end if
   end do

#else  !----------NITRIF_DENITRIF-------------!
   ! column loops to resolve plant/heterotroph/nitrifier/denitrifier competition for mineral N

   ! init total mineral N pools
   do fc=1,num_soilc
      c = filter_soilc(fc)
      sminn_tot(c) = 0.
   end do

   ! sum up total mineral N pools
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)
         sminn_tot(c) = sminn_tot(c) + (smin_no3_vr(c,j) + smin_nh4_vr(c,j)) * dzsoi_decomp(j)
      end do
   end do

   ! define N uptake profile for initial vertical distribution of plant N uptake, assuming plant seeks N from where it is most abundant
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)
         if (sminn_tot(c) .gt. 0.) then
               nuptake_prof(c,j) = sminn_vr(c,j) / sminn_tot(c)
         else
               nuptake_prof(c,j) = nfixation_prof(c,j)
         endif
      end do
   end do

   ! main column/vertical loop
   do j = 1, nlevdecomp  
      do fc=1,num_soilc
         c = filter_soilc(fc)
         l = clandunit(c)

         !  first compete for nh4
         sum_nh4_demand(c,j) = col_plant_ndemand(c) * nuptake_prof(c,j) + potential_immob_vr(c,j) + pot_f_nit_vr(c,j)
         sum_nh4_demand_scaled(c,j) = col_plant_ndemand(c)* nuptake_prof(c,j) * compet_plant_nh4 + &
              potential_immob_vr(c,j)*compet_decomp_nh4 + pot_f_nit_vr(c,j)*compet_nit
         if (sum_nh4_demand(c,j)*dt < smin_nh4_vr(c,j)) then
            ! NH4 availability is not limiting immobilization or plant
            ! uptake, and all can proceed at their potential rates
            nlimit_nh4(c,j) = 0
            fpi_nh4_vr(c,j) = 1.0_r8
            actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j)
            smin_nh4_to_plant_vr(c,j) = col_plant_ndemand(c) * nuptake_prof(c,j)
            
            f_nit_vr(c,j) = pot_f_nit_vr(c,j)
            
         else
            
            ! NH4 availability can not satisfy the sum of immobilization, nitrification, and
            ! plant growth demands, so these three demands compete for available
            ! soil mineral NH4 resource.
            nlimit_nh4(c,j) = 1
            if (sum_nh4_demand(c,j) > 0.0_r8) then
               actual_immob_nh4_vr(c,j) = min((smin_nh4_vr(c,j)/dt)*(potential_immob_vr(c,j)*compet_decomp_nh4 / sum_nh4_demand_scaled(c,j)), potential_immob_vr(c,j))
               smin_nh4_to_plant_vr(c,j) = min((smin_nh4_vr(c,j)/dt)*(col_plant_ndemand(c)*nuptake_prof(c,j)*compet_plant_nh4 / sum_nh4_demand_scaled(c,j)), col_plant_ndemand(c)*nuptake_prof(c,j))
               f_nit_vr(c,j) =  min((smin_nh4_vr(c,j)/dt)*(pot_f_nit_vr(c,j)*compet_nit / sum_nh4_demand_scaled(c,j)), pot_f_nit_vr(c,j))
            else
               actual_immob_nh4_vr(c,j) = 0.0_r8
               smin_nh4_to_plant_vr(c,j) = 0.0_r8
               f_nit_vr(c,j) = 0.0_r8
            end if
            
            if (potential_immob_vr(c,j) > 0.0_r8) then
               fpi_nh4_vr(c,j) = actual_immob_nh4_vr(c,j) / potential_immob_vr(c,j)
            else
               fpi_nh4_vr(c,j) = 0.0_r8
            end if
            
         end if

         ! next compete for no3
         sum_no3_demand(c,j) = (col_plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j)) + (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j)) + pot_f_denit_vr(c,j)
         sum_no3_demand_scaled(c,j) = (col_plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))*compet_plant_no3 + &
              (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))*compet_decomp_no3 + pot_f_denit_vr(c,j)*compet_denit

         if (sum_no3_demand(c,j)*dt < smin_no3_vr(c,j)) then
            ! NO3 availability is not limiting immobilization or plant
            ! uptake, and all can proceed at their potential rates
            nlimit_no3(c,j) = 1
            fpi_no3_vr(c,j) = 1.0_r8 -  fpi_nh4_vr(c,j)
            actual_immob_no3_vr(c,j) = (potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))
            smin_no3_to_plant_vr(c,j) = (col_plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))
            
            f_denit_vr(c,j) = pot_f_denit_vr(c,j)
            
         else 
            
            ! NO3 availability can not satisfy the sum of immobilization, denitrification, and
            ! plant growth demands, so these three demands compete for available
            ! soil mineral NO3 resource.
            nlimit_no3(c,j) = 1
            if (sum_no3_demand(c,j) > 0.0_r8) then
               actual_immob_no3_vr(c,j) = min((smin_no3_vr(c,j)/dt)*((potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))*compet_decomp_no3 / sum_no3_demand_scaled(c,j)), &
                    potential_immob_vr(c,j)-actual_immob_nh4_vr(c,j))
               smin_no3_to_plant_vr(c,j) = min((smin_no3_vr(c,j)/dt)*((col_plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))*compet_plant_no3 / sum_no3_demand_scaled(c,j)), &
                    col_plant_ndemand(c)*nuptake_prof(c,j)-smin_nh4_to_plant_vr(c,j))
               f_denit_vr(c,j) =  min((smin_no3_vr(c,j)/dt)*(pot_f_denit_vr(c,j)*compet_denit / sum_no3_demand_scaled(c,j)), pot_f_denit_vr(c,j))
            else
               actual_immob_no3_vr(c,j) = 0.0_r8
               smin_no3_to_plant_vr(c,j) = 0.0_r8
               f_denit_vr(c,j) = 0.0_r8
            end if
            
            if (potential_immob_vr(c,j) > 0.0_r8) then
               fpi_no3_vr(c,j) = actual_immob_no3_vr(c,j) / potential_immob_vr(c,j)
            else
               fpi_no3_vr(c,j) = 0.0_r8
            end if
            
         end if
         

         ! n2o emissions: n2o from nitr is const fraction, n2o from denitr is calculated in nitrif_denitrif
         f_n2o_nit_vr(c,j) = f_nit_vr(c,j) * nitrif_n2o_loss_frac
         f_n2o_denit_vr(c,j) = f_denit_vr(c,j) / (1._r8 + n2_n2o_ratio_denit_vr(c,j))
         
         
         ! this code block controls the addition of N to sminn pool
         ! to eliminate any N limitation, when Carbon_Only is set.  This lets the
         ! model behave essentially as a carbon-only model, but with the
         ! benefit of keeping track of the N additions needed to
         ! eliminate N limitations, so there is still a diagnostic quantity
         ! that describes the degree of N limitation at steady-state.
         
         if ( Carbon_only) then !.or. &
!              (crop_supln .and. (itypelun(l) == istcrop) .and. &
!               (ivt(pfti(c)) >= npcropmin)) ) then

            if ( fpi_no3_vr(c,j) + fpi_nh4_vr(c,j) .lt. 1._r8 ) then
               fpi_nh4_vr(c,j) = 1.0_r8 - fpi_no3_vr(c,j)
               supplement_to_sminn_vr(c,j) = (potential_immob_vr(c,j) -  actual_immob_no3_vr(c,j)) - actual_immob_nh4_vr(c,j) 
               ! update to new values that satisfy demand
               actual_immob_nh4_vr(c,j) = potential_immob_vr(c,j) -  actual_immob_no3_vr(c,j)   
            end if
            if ( smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j) .lt. col_plant_ndemand(c)*nuptake_prof(c,j) ) then
               supplement_to_sminn_vr(c,j) = supplement_to_sminn_vr(c,j) + &
                    (col_plant_ndemand(c)*nuptake_prof(c,j) - smin_no3_to_plant_vr(c,j)) - smin_nh4_to_plant_vr(c,j)  ! use old values
               smin_nh4_to_plant_vr(c,j) = col_plant_ndemand(c)*nuptake_prof(c,j) - smin_no3_to_plant_vr(c,j)
            end if
            sminn_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)
         end if

         ! sum up no3 and nh4 fluxes
         fpi_vr(c,j) = fpi_no3_vr(c,j) + fpi_nh4_vr(c,j)
         sminn_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + smin_nh4_to_plant_vr(c,j)
         actual_immob_vr(c,j) = actual_immob_no3_vr(c,j) + actual_immob_nh4_vr(c,j)
      end do
   end do
   
   do fc=1,num_soilc
      c = filter_soilc(fc)
      ! sum up N fluxes to plant after initial competition
      sminn_to_plant(c) = 0._r8
   end do   
   do j = 1, nlevdecomp  
      do fc=1,num_soilc
         c = filter_soilc(fc)
         sminn_to_plant(c) = sminn_to_plant(c) + sminn_to_plant_vr(c,j) * dzsoi_decomp(j)
      end do
   end do

   ! give plants a second pass to see if there is any mineral N left over with which to satisfy residual N demand.
   ! first take frm nh4 pool; then take from no3 pool
   do fc=1,num_soilc
      c = filter_soilc(fc)
      residual_plant_ndemand(c) = col_plant_ndemand(c) - sminn_to_plant(c)
      residual_smin_nh4(c) = 0._r8
   end do
   do j = 1, nlevdecomp  
      do fc=1,num_soilc
         c = filter_soilc(fc)
         if (residual_plant_ndemand(c) .gt. 0._r8 ) then
            if (nlimit_nh4(c,j) .eq. 0) then
               residual_smin_nh4_vr(c,j) = max(smin_nh4_vr(c,j) - (actual_immob_vr(c,j) + smin_nh4_to_plant_vr(c,j) ) * dt, 0._r8)
               residual_smin_nh4(c) = residual_smin_nh4(c) + residual_smin_nh4_vr(c,j) * dzsoi_decomp(j)
            else
               residual_smin_nh4_vr(c,j)  = 0._r8
            endif
            
            if ( residual_smin_nh4(c) .gt. 0._r8 .and. nlimit_nh4(c,j) .eq. 0 ) then
               smin_nh4_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + residual_smin_nh4_vr(c,j) * &
                    min(( residual_plant_ndemand(c) *  dt ) / residual_smin_nh4(c), 1._r8) / dt
            endif
         end if
      end do
   end do
   
   ! re-sum up N fluxes to plant after second pass for nh4
   do fc=1,num_soilc
      c = filter_soilc(fc)
      sminn_to_plant(c) = 0._r8
   end do
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)
         sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)
         sminn_to_plant(c) = sminn_to_plant(c) + (sminn_to_plant_vr(c,j)) * dzsoi_decomp(j)
      end do
   end do
   
   !
   ! and now do second pass for no3
   do fc=1,num_soilc
      c = filter_soilc(fc)
      residual_plant_ndemand(c) = col_plant_ndemand(c) - sminn_to_plant(c)
      residual_smin_no3(c) = 0._r8
   end do
   
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)
         if (residual_plant_ndemand(c) .gt. 0._r8 ) then
            if (nlimit_no3(c,j) .eq. 0) then
               residual_smin_no3_vr(c,j) = max(smin_no3_vr(c,j) - (actual_immob_vr(c,j) + smin_no3_to_plant_vr(c,j) ) * dt, 0._r8)
               residual_smin_no3(c) = residual_smin_no3(c) + residual_smin_no3_vr(c,j) * dzsoi_decomp(j)
            else
               residual_smin_no3_vr(c,j)  = 0._r8
            endif
            
            if ( residual_smin_no3(c) .gt. 0._r8 .and. nlimit_no3(c,j) .eq. 0) then
               smin_no3_to_plant_vr(c,j) = smin_no3_to_plant_vr(c,j) + residual_smin_no3_vr(c,j) * &
                    min(( residual_plant_ndemand(c) *  dt ) / residual_smin_no3(c), 1._r8) / dt
            endif
         endif
      end do
   end do
   
   ! re-sum up N fluxes to plant after second passes of both no3 and nh4
   do fc=1,num_soilc
      c = filter_soilc(fc)
      sminn_to_plant(c) = 0._r8
   end do
   do j = 1, nlevdecomp
      do fc=1,num_soilc
         c = filter_soilc(fc)
         sminn_to_plant_vr(c,j) = smin_nh4_to_plant_vr(c,j) + smin_no3_to_plant_vr(c,j)
         sminn_to_plant(c) = sminn_to_plant(c) + (sminn_to_plant_vr(c,j)) * dzsoi_decomp(j)
      end do
   end do
   
   ! sum up N fluxes to immobilization
   do fc=1,num_soilc
      c = filter_soilc(fc)
      actual_immob(c) = 0._r8
      potential_immob(c) = 0._r8
   end do
   do j = 1, nlevdecomp  
      do fc=1,num_soilc
         c = filter_soilc(fc)
         actual_immob(c) = actual_immob(c) + actual_immob_vr(c,j) * dzsoi_decomp(j)
         potential_immob(c) = potential_immob(c) + potential_immob_vr(c,j) * dzsoi_decomp(j)
      end do
   end do


   do fc=1,num_soilc
      c = filter_soilc(fc)   
      ! calculate the fraction of potential growth that can be
      ! acheived with the N available to plants
      if (col_plant_ndemand(c) > 0.0_r8) then
         fpg(c) = sminn_to_plant(c) / col_plant_ndemand(c)
      else
         fpg(c) = 1._r8
      end if
      
      ! calculate the fraction of immobilization realized (for diagnostic purposes)
      if (potential_immob(c) > 0.0_r8) then
         fpi(c) = actual_immob(c) / potential_immob(c)
      else
         fpi(c) = 1._r8
      end if
   end do ! end of column loops
#endif


   
   ! start new pft loop to distribute the available N between the
   ! competing pfts on the basis of relative demand, and allocate C and N to
   ! new growth and storage
   
   do fp=1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)

      ! set some local allocation variables
      f1 = froot_leaf(ivt(p))
      f2 = croot_stem(ivt(p))
      
      ! modified wood allocation to be 2.2 at npp=800 gC/m2/yr, 0.2 at npp=0,
      ! constrained so that it does not go lower than 0.2 (under negative annsum_npp)
      ! There was an error in this formula in previous version, where the coefficient
      ! was 0.004 instead of 0.0025.
      ! This variable allocation is only for trees. Shrubs have a constant
      ! allocation as specified in the pft-physiology file.  The value is also used
      ! as a trigger here: -1.0 means to use the dynamic allocation (trees).
      if (stem_leaf(ivt(p)) == -1._r8) then
        f3 = (2.7/(1.0+exp(-0.004*(annsum_npp(p) - 300.0)))) - 0.4
      else
        f3 = stem_leaf(ivt(p))
      end if
      
      f4 = flivewd(ivt(p))
      g1 = grperc(ivt(p))
      g2 = grpnow(ivt(p))
      cnl = leafcn(ivt(p))
      cnfr = frootcn(ivt(p))
      cnlw = livewdcn(ivt(p))
      cndw = deadwdcn(ivt(p))
      fcur = fcur2(ivt(p))

      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
         if (croplive(p)) then
            f1 = aroot(p) / aleaf(p)
            f3 = astem(p) / aleaf(p)
            f5 = arepr(p) / aleaf(p)
            g1 = 0.25_r8
         else
            f1 = 0._r8
            f3 = 0._r8
            f5 = 0._r8
            g1 = 0.25_r8
         end if
      end if

      ! increase fcur linearly with ndays_active, until fcur reaches 1.0 at
      ! ndays_active = days/year.  This prevents the continued storage of C and N.
      ! turning off this correction (PET, 12/11/03), instead using bgtr in
      ! phenology algorithm.
      !fcur = fcur + (1._r8 - fcur)*lgsf(p)
      sminn_to_npool(p) = plant_ndemand(p) * fpg(c)
      plant_nalloc(p) = sminn_to_npool(p) + retransn_to_npool(p)


      ! calculate the associated carbon allocation, and the excess
      ! carbon flux that must be accounted for through downregulation
      plant_calloc(p) = plant_nalloc(p) * (c_allometry(p)/n_allometry(p))
      excess_cflux(p) = availc(p) - plant_calloc(p)


      ! reduce gpp fluxes due to N limitation
      if (gpp(p) > 0.0_r8) then
         downreg(p) = excess_cflux(p)/gpp(p)
         psnsun_to_cpool(p) = psnsun_to_cpool(p)*(1._r8 - downreg(p))
         psnshade_to_cpool(p) = psnshade_to_cpool(p)*(1._r8 - downreg(p))
         if ( use_c13 ) then
            c13_psnsun_to_cpool(p) = c13_psnsun_to_cpool(p)*(1._r8 - downreg(p))
            c13_psnshade_to_cpool(p) = c13_psnshade_to_cpool(p)*(1._r8 - downreg(p))
         endif

         if ( use_c14 ) then
            c14_psnsun_to_cpool(p) = c14_psnsun_to_cpool(p)*(1._r8 - downreg(p))
            c14_psnshade_to_cpool(p) = c14_psnshade_to_cpool(p)*(1._r8 - downreg(p))
         endif
      end if

      ! calculate the amount of new leaf C dictated by these allocation
      ! decisions, and calculate the daily fluxes of C and N to current
      ! growth and storage pools

      ! fcur is the proportion of this day's growth that is displayed now,
      ! the remainder going into storage for display next year through the
      ! transfer pools

      nlc = plant_calloc(p) / c_allometry(p)
      cpool_to_leafc(p)          = nlc * fcur
      cpool_to_leafc_storage(p)  = nlc * (1._r8 - fcur)
      cpool_to_frootc(p)         = nlc * f1 * fcur
      cpool_to_frootc_storage(p) = nlc * f1 * (1._r8 - fcur)
      if (woody(ivt(p)) == 1._r8) then
         cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
         cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1._r8 - fcur)
         cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4) * fcur
         cpool_to_deadstemc_storage(p)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
         cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
         cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
         cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
         cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
         cpool_to_livestemc(p)          = nlc * f3 * f4 * fcur
         cpool_to_livestemc_storage(p)  = nlc * f3 * f4 * (1._r8 - fcur)
         cpool_to_deadstemc(p)          = nlc * f3 * (1._r8 - f4) * fcur
         cpool_to_deadstemc_storage(p)  = nlc * f3 * (1._r8 - f4) * (1._r8 - fcur)
         cpool_to_livecrootc(p)         = nlc * f2 * f3 * f4 * fcur
         cpool_to_livecrootc_storage(p) = nlc * f2 * f3 * f4 * (1._r8 - fcur)
         cpool_to_deadcrootc(p)         = nlc * f2 * f3 * (1._r8 - f4) * fcur
         cpool_to_deadcrootc_storage(p) = nlc * f2 * f3 * (1._r8 - f4) * (1._r8 - fcur)
         cpool_to_grainc(p)             = nlc * f5 * fcur
         cpool_to_grainc_storage(p)     = nlc * f5 * (1._r8 -fcur)
      end if

      ! corresponding N fluxes
      npool_to_leafn(p)          = (nlc / cnl) * fcur
      npool_to_leafn_storage(p)  = (nlc / cnl) * (1._r8 - fcur)
      npool_to_frootn(p)         = (nlc * f1 / cnfr) * fcur
      npool_to_frootn_storage(p) = (nlc * f1 / cnfr) * (1._r8 - fcur)
      if (woody(ivt(p)) == 1._r8) then
         npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw) * fcur
         npool_to_livestemn_storage(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
         npool_to_deadstemn(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
         npool_to_deadstemn_storage(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
         npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
         npool_to_livecrootn_storage(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
         npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
         npool_to_deadcrootn_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
         cng = graincn(ivt(p))
         npool_to_livestemn(p)          = (nlc * f3 * f4 / cnlw) * fcur
         npool_to_livestemn_storage(p)  = (nlc * f3 * f4 / cnlw) * (1._r8 - fcur)
         npool_to_deadstemn(p)          = (nlc * f3 * (1._r8 - f4) / cndw) * fcur
         npool_to_deadstemn_storage(p)  = (nlc * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
         npool_to_livecrootn(p)         = (nlc * f2 * f3 * f4 / cnlw) * fcur
         npool_to_livecrootn_storage(p) = (nlc * f2 * f3 * f4 / cnlw) * (1._r8 - fcur)
         npool_to_deadcrootn(p)         = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * fcur
         npool_to_deadcrootn_storage(p) = (nlc * f2 * f3 * (1._r8 - f4) / cndw) * (1._r8 - fcur)
         npool_to_grainn(p)             = (nlc * f5 / cng) * fcur
         npool_to_grainn_storage(p)     = (nlc * f5 / cng) * (1._r8 -fcur)
      end if

      ! Calculate the amount of carbon that needs to go into growth
      ! respiration storage to satisfy all of the storage growth demands.
      ! Allows for the fraction of growth respiration that is released at the
      ! time of fixation, versus the remaining fraction that is stored for
      ! release at the time of display. Note that all the growth respiration
      ! fluxes that get released on a given timestep are calculated in growth_resp(),
      ! but that the storage of C for growth resp during display of transferred
      ! growth is assigned here.

      gresp_storage = cpool_to_leafc_storage(p) + cpool_to_frootc_storage(p)
      if (woody(ivt(p)) == 1._r8) then
         gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
         gresp_storage = gresp_storage + cpool_to_deadstemc_storage(p)
         gresp_storage = gresp_storage + cpool_to_livecrootc_storage(p)
         gresp_storage = gresp_storage + cpool_to_deadcrootc_storage(p)
      end if
      if (ivt(p) >= npcropmin) then ! skip 2 generic crops
         gresp_storage = gresp_storage + cpool_to_livestemc_storage(p)
         gresp_storage = gresp_storage + cpool_to_grainc_storage(p)
      end if
      cpool_to_gresp_storage(p) = gresp_storage * g1 * (1._r8 - g2)

   end do ! end pft loop


end subroutine CNAllocation

#endif

end module CNAllocationMod
