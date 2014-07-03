
module CNDecompCascadeMod_CENTURY
#ifdef CN

#ifdef CENTURY_DECOMP


!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNDecompMod
!
! !DESCRIPTION:
! Module that sets the coeffiecients used in the decomposition cascade submodel.  This uses the CENTURY parameters
!
! !USES:
   use shr_kind_mod , only: r8 => shr_kind_r8
   use shr_const_mod, only: SHR_CONST_TKFRZ
   use clm_varpar   , only: nlevsoi, nlevgrnd, nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools, nsompools 
   use clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
   use clm_varctl   , only: iulog, spinup_state
   use clm_varcon   , only: zsoi
#ifdef LCH4
   use clm_varctl   , only: anoxia
   use ch4varcon    , only: mino2lim
#endif
   use abortutils,   only: endrun

   implicit none
   save
   private
! !PUBLIC MEMBER FUNCTIONS:
   public:: init_decompcascade, decomp_rate_constants
!
! !PUBLIC DATA MEMBERS:
#if (defined VERTSOILC)
   real(r8), public :: decomp_depth_efolding = 0.5_r8        ! (meters) e-folding depth for reduction in decomposition [set to large number for depth-independance]
#endif

   logical, public :: normalize_q10_to_century_tfunc = .true.! do we normalize the century decomp. rates so that they match the CLM Q10 at a given tep?
   real(r8), public :: normalization_tref = 15._r8           ! reference temperature for normalizaion (degrees C)
   logical, public :: use_century_tfunc = .false.
#ifdef LCH4
   logical,  public :: anoxia_wtsat = .false.                ! true ==> weight anoxia by inundated fraction
#endif
   real(r8), public :: froz_q10 = 1.5_r8                     ! separate q10 for frozen soil respiration rates.  default to same as above zero rates
   integer, public :: nlev_soildecomp_standard               ! used here and in ch4Mod

   !! parameters for AD spinup
   real(r8), public, parameter :: spinup_vector(nsompools) = (/ 1.0_r8, 15.0_r8, 675.0_r8 /) ! multipliers for soil decomp during accelerated spinup

!EOP
!-----------------------------------------------------------------------

contains


!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: init_decompcascade
!
! !INTERFACE:
subroutine init_decompcascade(begc, endc)

! !DESCRIPTION:
!
!  initialize rate constants and decomposition pathways following the decomposition cascade of the CENTURY model.
!  written by C. Koven 
!
! !USES:
   use clmtype
   use clm_time_manager    , only : get_step_size

! !ARGUMENTS:
   implicit none
!
! !CALLED FROM:
! 
!
! !REVISION HISTORY:
!
   ! column level
   integer  :: begc, endc       ! per-proc beginning and ending column indices

   !-- properties of each pathway along decomposition cascade 
   character(len=8), pointer :: cascade_step_name(:)       ! name of transition
   real(r8), pointer :: rf_decomp_cascade(:,:,:)           ! respired fraction in decomposition step (frac)
   integer,  pointer :: cascade_donor_pool(:)              ! which pool is C taken from for a given decomposition step
   integer,  pointer :: cascade_receiver_pool(:)           ! which pool is C added to for a given decomposition step
   real(r8), pointer :: pathfrac_decomp_cascade(:,:,:)     ! what fraction of C leaving a given pool passes through a given transition (frac)
   !-- properties of each decomposing pool
   logical,  pointer :: floating_cn_ratio_decomp_pools(:)  ! TRUE => pool has fixed C:N ratio
   character(len=8), pointer :: decomp_pool_name_restart(:)! name of pool for restart files
   character(len=8), pointer :: decomp_pool_name_history(:)! name of pool for history files
   character(len=20), pointer :: decomp_pool_name_long(:)  ! name of pool for netcdf long names
   character(len=8), pointer :: decomp_pool_name_short(:)  ! name of pool for netcdf short names
   logical, pointer :: is_litter(:)                        ! TRUE => pool is a litter pool
   logical, pointer :: is_soil(:)                          ! TRUE => pool is a soil pool
   logical, pointer :: is_cwd(:)                           ! TRUE => pool is a cwd pool
   real(r8), pointer :: initial_cn_ratio(:)                ! c:n ratio for initialization of pools
   real(r8), pointer :: initial_stock(:)                   ! initial concentration for seeding at spinup
   logical, pointer :: is_metabolic(:)                     ! TRUE => pool is metabolic material
   logical, pointer :: is_cellulose(:)                     ! TRUE => pool is cellulose
   logical, pointer :: is_lignin(:)                        ! TRUE => pool is lignin
   real(r8), pointer :: cellclay(:,:)                      ! column 3D clay
   real(r8), pointer :: cellsand(:,:)                      ! column 3D sand
   real(r8), pointer :: spinup_factor(:)                   ! factor for AD spinup associated with each pool
   real(r8) :: rf_l1s1
   real(r8) :: rf_l2s1
   real(r8) :: rf_l3s2
   real(r8) :: rf_s1s2(begc:endc,1:nlevdecomp)
   real(r8) :: rf_s1s3(begc:endc,1:nlevdecomp)
   real(r8) :: rf_s2s1
   real(r8) :: rf_s2s3
   real(r8) :: rf_s3s1
   real(r8) :: rf_cwdl2
   real(r8) :: rf_cwdl3
   real(r8):: cwd_fcel
   real(r8):: cwd_flig
   real(r8) :: cn_s1
   real(r8) :: cn_s2
   real(r8) :: cn_s3
   real(r8) :: cn_s4
   real(r8) :: f_s1s2(begc:endc,1:nlevdecomp)
   real(r8) :: f_s1s3(begc:endc,1:nlevdecomp)
   real(r8) :: f_s2s1
   real(r8) :: f_s2s3

   integer :: i_litr1
   integer :: i_litr2
   integer :: i_litr3
   integer :: i_soil1
   integer :: i_soil2
   integer :: i_soil3
   integer :: i_l1s1
   integer :: i_l2s1
   integer :: i_l3s2
   integer :: i_s1s2
   integer :: i_s1s3
   integer :: i_s2s1
   integer :: i_s2s3
   integer :: i_s3s1
   integer :: i_cwdl2
   integer :: i_cwdl3

   integer :: c, j     ! indices
   real(r8) :: t       ! temporary variable


   cascade_step_name                       => decomp_cascade_con%cascade_step_name
   rf_decomp_cascade                       => clm3%g%l%c%cps%rf_decomp_cascade
   cascade_donor_pool                      => decomp_cascade_con%cascade_donor_pool
   cascade_receiver_pool                   => decomp_cascade_con%cascade_receiver_pool
   pathfrac_decomp_cascade                 => clm3%g%l%c%cps%pathfrac_decomp_cascade
   floating_cn_ratio_decomp_pools          => decomp_cascade_con%floating_cn_ratio_decomp_pools
   decomp_pool_name_restart                => decomp_cascade_con%decomp_pool_name_restart
   decomp_pool_name_history                => decomp_cascade_con%decomp_pool_name_history
   decomp_pool_name_long                   => decomp_cascade_con%decomp_pool_name_long
   decomp_pool_name_short                  => decomp_cascade_con%decomp_pool_name_short
   is_litter                               => decomp_cascade_con%is_litter
   is_soil                                 => decomp_cascade_con%is_soil
   is_cwd                                  => decomp_cascade_con%is_cwd
   initial_cn_ratio                        => decomp_cascade_con%initial_cn_ratio
   initial_stock                           => decomp_cascade_con%initial_stock
   is_metabolic                            => decomp_cascade_con%is_metabolic
   is_cellulose                            => decomp_cascade_con%is_cellulose
   is_lignin                               => decomp_cascade_con%is_lignin
   cellclay                                => clm3%g%l%c%cps%cellclay
   cellsand                                => clm3%g%l%c%cps%cellsand
   spinup_factor                           => decomp_cascade_con%spinup_factor

   !------- time-constant coefficients ---------- !
   ! set soil organic matter compartment C:N ratios
   cn_s1 = 8.0_r8
   cn_s2 = 11.0_r8
   cn_s3 = 11.0_r8

   ! set respiration fractions for fluxes between compartments
   rf_l1s1 = 0.55_r8
   rf_l2s1 = 0.5_r8
   rf_l3s2 = 0.5_r8
   rf_s2s1 = 0.55_r8
   rf_s2s3 = 0.55_r8
   rf_s3s1 = 0.55_r8
   rf_cwdl2 = 0._r8
   rf_cwdl3 = 0._r8


   ! set the cellulose and lignin fractions for coarse woody debris
   cwd_fcel = 0.76_r8
   cwd_flig = 0.24_r8

   ! set path fractions
   f_s2s1 = 0.42_r8/(0.45_r8)
   f_s2s3 = 0.03_r8/(0.45_r8)

   ! some of these are dependent on the soil texture properties
   do c = begc, endc
      do j = 1, nlevdecomp
         t = 0.85_r8 - 0.68_r8 * 0.01_r8 * (100._r8 - cellsand(c,j))
         f_s1s2(c,j) = 1._r8 - .004_r8 / (1._r8 - t)
         f_s1s3(c,j) = .004_r8 / (1._r8 - t)
         rf_s1s2(c,j) = t
         rf_s1s3(c,j) = t
      end do
   end do


   !-------------------  list of pools and their attributes  ------------

   i_litr1 = i_met_lit
   floating_cn_ratio_decomp_pools(i_litr1) = .true.
   decomp_pool_name_restart(i_litr1) = 'litr1'
   decomp_pool_name_history(i_litr1) = 'LITR1'
   decomp_pool_name_long(i_litr1) = 'litter 1'
   decomp_pool_name_short(i_litr1) = 'L1'
   is_litter(i_litr1) = .true.
   is_soil(i_litr1) = .false.
   is_cwd(i_litr1) = .false.
   initial_cn_ratio(i_litr1) = 90._r8
   initial_stock(i_litr1) = 0._r8
   is_metabolic(i_litr1) = .true.
   is_cellulose(i_litr1) = .false.
   is_lignin(i_litr1) = .false.
   
   i_litr2 = i_cel_lit
   floating_cn_ratio_decomp_pools(i_litr2) = .true.
   decomp_pool_name_restart(i_litr2) = 'litr2'
   decomp_pool_name_history(i_litr2) = 'LITR2'
   decomp_pool_name_long(i_litr2) = 'litter 2'
   decomp_pool_name_short(i_litr2) = 'L2'
   is_litter(i_litr2) = .true.
   is_soil(i_litr2) = .false.
   is_cwd(i_litr2) = .false.
   initial_cn_ratio(i_litr2) = 90._r8
   initial_stock(i_litr2) = 0._r8
   is_metabolic(i_litr2) = .false.
   is_cellulose(i_litr2) = .true.
   is_lignin(i_litr2) = .false.

   i_litr3 = i_lig_lit
   floating_cn_ratio_decomp_pools(i_litr3) = .true.
   decomp_pool_name_restart(i_litr3) = 'litr3'
   decomp_pool_name_history(i_litr3) = 'LITR3'
   decomp_pool_name_long(i_litr3) = 'litter 3'
   decomp_pool_name_short(i_litr3) = 'L3'
   is_litter(i_litr3) = .true.
   is_soil(i_litr3) = .false.
   is_cwd(i_litr3) = .false.
   initial_cn_ratio(i_litr3) = 90._r8
   initial_stock(i_litr3) = 0._r8
   is_metabolic(i_litr3) = .false.
   is_cellulose(i_litr3) = .false.
   is_lignin(i_litr3) = .true.

   ! CWD
   floating_cn_ratio_decomp_pools(i_cwd) = .true.
   decomp_pool_name_restart(i_cwd) = 'cwd'
   decomp_pool_name_history(i_cwd) = 'CWD'
   decomp_pool_name_long(i_cwd) = 'coarse woody debris'
   decomp_pool_name_short(i_cwd) = 'CWD'
   is_litter(i_cwd) = .false.
   is_soil(i_cwd) = .false.
   is_cwd(i_cwd) = .true.
   initial_cn_ratio(i_cwd) = 90._r8
   initial_stock(i_cwd) = 0._r8
   is_metabolic(i_cwd) = .false.
   is_cellulose(i_cwd) = .false.
   is_lignin(i_cwd) = .false.

   i_soil1 = 5
   floating_cn_ratio_decomp_pools(i_soil1) = .false.
   decomp_pool_name_restart(i_soil1) = 'soil1'
   decomp_pool_name_history(i_soil1) = 'SOIL1'
   decomp_pool_name_long(i_soil1) = 'soil 1'
   decomp_pool_name_short(i_soil1) = 'S1'
   is_litter(i_soil1) = .false.
   is_soil(i_soil1) = .true.
   is_cwd(i_soil1) = .false.
   initial_cn_ratio(i_soil1) = cn_s1
   initial_stock(i_soil1) = 20._r8
   is_metabolic(i_soil1) = .false.
   is_cellulose(i_soil1) = .false.
   is_lignin(i_soil1) = .false.

   i_soil2 = 6
   floating_cn_ratio_decomp_pools(i_soil2) = .false.
   decomp_pool_name_restart(i_soil2) = 'soil2'
   decomp_pool_name_history(i_soil2) = 'SOIL2'
   decomp_pool_name_long(i_soil2) = 'soil 2'
   decomp_pool_name_short(i_soil2) = 'S2'
   is_litter(i_soil2) = .false.
   is_soil(i_soil2) = .true.
   is_cwd(i_soil2) = .false.
   initial_cn_ratio(i_soil2) = cn_s2
   initial_stock(i_soil2) = 20._r8
   is_metabolic(i_soil2) = .false.
   is_cellulose(i_soil2) = .false.
   is_lignin(i_soil2) = .false.

   i_soil3 = 7
   floating_cn_ratio_decomp_pools(i_soil3) = .false.
   decomp_pool_name_restart(i_soil3) = 'soil3'
   decomp_pool_name_history(i_soil3) = 'SOIL3'
   decomp_pool_name_long(i_soil3) = 'soil 3'
   decomp_pool_name_short(i_soil3) = 'S3'
   is_litter(i_soil3) = .false.
   is_soil(i_soil3) = .true.
   is_cwd(i_soil3) = .false.
   initial_cn_ratio(i_soil3) = cn_s3
   initial_stock(i_soil3) = 20._r8
   is_metabolic(i_soil3) = .false.
   is_cellulose(i_soil3) = .false.
   is_lignin(i_soil3) = .false.

   spinup_factor(i_litr1) = 1._r8
   spinup_factor(i_litr2) = 1._r8
   spinup_factor(i_litr3) = 1._r8
   spinup_factor(i_cwd) = 1._r8
   spinup_factor(i_soil1) = spinup_vector(1)
   spinup_factor(i_soil2) = spinup_vector(2)
   spinup_factor(i_soil3) = spinup_vector(3)

   !----------------  list of transitions and their time-independent coefficients  ---------------!
   i_l1s1 = 1
   cascade_step_name(i_l1s1) = 'L1S1'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_l1s1) = rf_l1s1
   cascade_donor_pool(i_l1s1) = i_litr1
   cascade_receiver_pool(i_l1s1) = i_soil1
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_l1s1) = 1.0_r8

   i_l2s1 = 2
   cascade_step_name(i_l2s1) = 'L2S1'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_l2s1) = rf_l2s1
   cascade_donor_pool(i_l2s1) = i_litr2
   cascade_receiver_pool(i_l2s1) = i_soil1
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_l2s1)= 1.0_r8

   i_l3s2 = 3
   cascade_step_name(i_l3s2) = 'L3S2'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_l3s2) = rf_l3s2
   cascade_donor_pool(i_l3s2) = i_litr3
   cascade_receiver_pool(i_l3s2) = i_soil2
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_l3s2) = 1.0_r8

   i_s1s2 = 4
   cascade_step_name(i_s1s2) = 'S1S2'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s1s2) = rf_s1s2(begc:endc,1:nlevdecomp)
   cascade_donor_pool(i_s1s2) = i_soil1
   cascade_receiver_pool(i_s1s2) = i_soil2
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s1s2) = f_s1s2(begc:endc,1:nlevdecomp)

   i_s1s3 = 5
   cascade_step_name(i_s1s3) = 'S1S3'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s1s3) = rf_s1s3(begc:endc,1:nlevdecomp)
   cascade_donor_pool(i_s1s3) = i_soil1
   cascade_receiver_pool(i_s1s3) = i_soil3
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s1s3) = f_s1s3(begc:endc,1:nlevdecomp)

   i_s2s1 = 6
   cascade_step_name(i_s2s1) = 'S2S1'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s2s1) = rf_s2s1
   cascade_donor_pool(i_s2s1) = i_soil2
   cascade_receiver_pool(i_s2s1) = i_soil1
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s2s1) = f_s2s1

   i_s2s3 = 7 
   cascade_step_name(i_s2s3) = 'S2S3'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s2s3) = rf_s2s3
   cascade_donor_pool(i_s2s3) = i_soil2
   cascade_receiver_pool(i_s2s3) = i_soil3
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s2s3) = f_s2s3

   i_s3s1 = 8
   cascade_step_name(i_s3s1) = 'S3S1'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s3s1) = rf_s3s1
   cascade_donor_pool(i_s3s1) = i_soil3
   cascade_receiver_pool(i_s3s1) = i_soil1
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s3s1) = 1.0_r8

   i_cwdl2 = 9
   cascade_step_name(i_cwdl2) = 'CWDL2'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl2) = rf_cwdl2
   cascade_donor_pool(i_cwdl2) = i_cwd
   cascade_receiver_pool(i_cwdl2) = i_litr2
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl2) = cwd_fcel

   i_cwdl3 = 10
   cascade_step_name(i_cwdl3) = 'CWDL3'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl3) = rf_cwdl3
   cascade_donor_pool(i_cwdl3) = i_cwd
   cascade_receiver_pool(i_cwdl3) = i_litr3
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl3) = cwd_flig


end subroutine init_decompcascade





!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: decomp_rate_constants
!
! !INTERFACE:
subroutine decomp_rate_constants(lbc, ubc, num_soilc, filter_soilc)
 !
! !DESCRIPTION:
!
!  calculate rate constants and decomposition pathways for teh CENTURY decomposition cascade model
!  written by C. Koven based on original CLM4 decomposition cascade
!
! !USES:
   use clmtype
   use shr_const_mod, only : SHR_CONST_PI
   use clm_varcon, only: secspday
   use clm_time_manager, only : get_days_per_year

   !
! !ARGUMENTS:
   implicit none
   integer, intent(in) :: lbc, ubc        ! column bounds
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! 
!
! !REVISION HISTORY:
!
   ! column level
   real(r8), pointer :: decomp_k(:,:,:)   ! rate constant for decomposition (1./sec)
   real(r8), pointer :: t_scalar(:,:)     ! soil temperature scalar for decomp
   real(r8), pointer :: w_scalar(:,:)     ! soil water scalar for decomp
   real(r8), pointer :: o_scalar(:,:)     ! fraction by which decomposition is limited by anoxia

   real(r8), pointer :: dz(:,:)           ! soil layer thickness (m)
   real(r8), pointer :: t_soisno(:,:)     ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   real(r8), pointer :: sucsat(:,:)      ! minimum soil suction (mm)
   real(r8), pointer :: soilpsi(:,:)      ! soil water potential in each soil layer (MPa)
#ifdef LCH4
   real(r8), pointer :: o2stress_unsat(:,:) ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   real(r8), pointer :: o2stress_sat(:,:) ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   real(r8), pointer :: finundated(:)     ! fractional inundated area
#endif
   integer, pointer :: alt_indx(:)        ! current depth of thaw

!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   real(r8):: frw(lbc:ubc)          ! rooting fraction weight
   real(r8), pointer:: fr(:,:)      ! column-level rooting fraction by soil depth
   real(r8):: minpsi, maxpsi        ! limits for soil water scalar for decomp
   real(r8):: psi                   ! temporary soilpsi for water scalar
   !   real(r8):: w_scalar(lbc:ubc,1:nlevdecomp)     !soil water scalar for decomp
   real(r8):: rate_scalar           ! combined rate scalar for decomp
   real(r8):: k_l1                  ! decomposition rate constant litter 1 (1/sec)
   real(r8):: k_l2_l3               ! decomposition rate constant litter 2 and litter 3 (1/sec)
   real(r8):: k_s1                  ! decomposition rate constant SOM 1 (1/sec)
   real(r8):: k_s2                  ! decomposition rate constant SOM 2 (1/sec)
   real(r8):: k_s3                  ! decomposition rate constant SOM 3 (1/sec)
   real(r8):: k_frag                ! fragmentation rate constant CWD (1/sec)
   real(r8):: tau_l1                ! turnover time of  litter 1 (yr)
   real(r8):: tau_l2_l3             ! turnover time of  litter 2 and litter 3 (yr)
   real(r8):: tau_l3                ! turnover time of  litter 3 (yr)
   real(r8):: tau_s1                ! turnover time of  SOM 1 (yr)
   real(r8):: tau_s2                ! turnover time of  SOM 2 (yr)
   real(r8):: tau_s3                ! turnover time of  SOM 3 (yr)
   real(r8):: tau_cwd               ! corrected fragmentation rate constant CWD
   real(r8):: cwd_fcel              ! cellulose fraction of coarse woody debris
   real(r8):: cwd_flig              ! lignin fraction of coarse woody debris
   real(r8):: cwdc_loss             ! fragmentation rate for CWD carbon (gC/m2/s)
   real(r8):: cwdn_loss             ! fragmentation rate for CWD nitrogen (gN/m2/s)

   integer :: i_litr1
   integer :: i_litr2
   integer :: i_litr3
   integer :: i_soil1
   integer :: i_soil2
   integer :: i_soil3
   integer :: c, fc, j, k, l
   real(r8) :: q10 = 1.5_r8
   real(r8) :: catanf    ! hyperbolic temperature function from CENTURY
   real(r8) :: catanf_30 ! reference rate at 30C
   real(r8) :: t1        ! temperature argument

   real(r8) :: normalization_factor ! factor by which to offset the decomposition rates frm century to a q10 formulation
   real(r8):: days_per_year         ! days per year


#if (defined VERTSOILC)
   real(r8) :: depth_scalar(lbc:ubc,1:nlevdecomp) 
#endif

   !----- CENTURY T response function
   catanf(t1) = 11.75_r8 +(29.7_r8 / SHR_CONST_PI) * atan( SHR_CONST_PI * 0.031_r8  * ( t1 - 15.4_r8 ))

   ! Assign local pointers to derived type arrays
   t_soisno              => clm3%g%l%c%ces%t_soisno
   sucsat            => clm3%g%l%c%cps%sucsat
   soilpsi               => clm3%g%l%c%cps%soilpsi
   dz                    => clm3%g%l%c%cps%dz
   t_scalar           => clm3%g%l%c%ccf%t_scalar
   w_scalar           => clm3%g%l%c%ccf%w_scalar
   o_scalar           => clm3%g%l%c%ccf%o_scalar
   decomp_k              => clm3%g%l%c%ccf%decomp_k
#ifdef LCH4
   o2stress_sat          => clm3%g%l%c%cch4%o2stress_sat
   o2stress_unsat        => clm3%g%l%c%cch4%o2stress_unsat
   finundated            => clm3%g%l%c%cws%finundated
#endif
   alt_indx              => clm3%g%l%c%cps%alt_indx

   if ( use_century_tfunc .and. normalize_q10_to_century_tfunc ) then
      call endrun( 'ERROR: cannot have both use_century_tfunc and normalize_q10_to_century_tfunc set as true' )
   endif
   
   days_per_year = get_days_per_year()

   ! tau (yrs) at reference temperature
   ! the aboveground parameters from century
   ! tau_l1 = 1./14.8
   ! tau_l2 = 1./3.9
   ! tau_s1 = 1./6.0
   ! tau_l2 = 1./0.2
   ! tau_l3 = 1./.0045

   ! the belowground parameters from century
   tau_l1 = 1./18.5
   tau_l2_l3 = 1./4.9
   tau_s1 = 1./7.3
   tau_s2 = 1./0.2
   tau_s3 = 1./.0045

   ! century leaves wood decomposition rates open, within range of 0 - 0.5 yr^-1
   tau_cwd  = 1./0.3

   ! translate to per-second time constant
   k_l1 = 1._r8 / (secspday * days_per_year * tau_l1)
   k_l2_l3 = 1._r8 / (secspday * days_per_year * tau_l2_l3)
   k_s1 = 1._r8 / (secspday * days_per_year * tau_s1)
   k_s2 = 1._r8 / (secspday * days_per_year * tau_s2)
   k_s3 = 1._r8 / (secspday * days_per_year * tau_s3)
   k_frag = 1._r8 / (secspday * days_per_year * tau_cwd)
   

   ! calc ref rate
   catanf_30 = catanf(30._r8)
   ! The following code implements the acceleration part of the AD spinup algorithm

if ( spinup_state .eq. 1 ) then
   k_s1 = k_s1 * spinup_vector(1)
   k_s2 = k_s2 * spinup_vector(2)
   k_s3 = k_s3 * spinup_vector(3)
endif

   i_litr1 = 1
   i_litr2 = 2
   i_litr3 = 3
   i_soil1 = 5
   i_soil2 = 6
   i_soil3 = 7


   !--- time dependent coefficients-----!
   if ( nlevdecomp .eq. 1 ) then
      
      ! calculate function to weight the temperature and water potential scalars
      ! for decomposition control.  
      
      
      ! the following normalizes values in fr so that they
      ! sum to 1.0 across top nlevdecomp levels on a column
      frw(lbc:ubc) = 0._r8
      nlev_soildecomp_standard=5
      allocate(fr(lbc:ubc,nlev_soildecomp_standard))
      do j=1,nlev_soildecomp_standard
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            frw(c) = frw(c) + dz(c,j)
         end do
      end do
      do j = 1,nlev_soildecomp_standard
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if (frw(c) /= 0._r8) then
               fr(c,j) = dz(c,j) / frw(c)
            else
               fr(c,j) = 0._r8
            end if
         end do
      end do

      if ( .not. use_century_tfunc ) then
         ! calculate rate constant scalar for soil temperature
         ! assuming that the base rate constants are assigned for non-moisture
         ! limiting conditions at 25 C. 

         do j = 1,nlev_soildecomp_standard
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if (j==1) t_scalar(c,:) = 0._r8
               !! use separate (possibly equal) t funcs above and below freezing point
               !! t_scalar(c,1)=t_scalar(c,1) + (q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))*fr(c,j) 
               if (t_soisno(c,j) .ge. SHR_CONST_TKFRZ) then
                  t_scalar(c,1)=t_scalar(c,1) + (q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))*fr(c,j)
               else
                  t_scalar(c,1)=t_scalar(c,1) + (q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))*fr(c,j)
               endif
            end do
         end do

      else
         ! original century uses an arctangent function to calculate the temperature dependence of decomposition
         do j = 1,nlev_soildecomp_standard
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               if (j==1) t_scalar(c,:) = 0._r8
               
               t_scalar(c,1)=t_scalar(c,1) +max(catanf(t_soisno(c,j)-SHR_CONST_TKFRZ)/catanf_30*fr(c,j),0.01_r8)
            end do
         end do

      endif
      
      ! calculate the rate constant scalar for soil water content.
      ! Uses the log relationship with water potential given in
      ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
      ! a comparison of models. Ecology, 68(5):1190-1200.
      ! and supported by data in
      ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
      ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.
      
      minpsi = -10.0_r8;

      do j = 1,nlev_soildecomp_standard
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if (j==1) w_scalar(c,:) = 0._r8
            maxpsi = sucsat(c,j) * (-9.8e-6_r8)
            psi = min(soilpsi(c,j),maxpsi)
            ! decomp only if soilpsi is higher than minpsi
            if (psi > minpsi) then
               w_scalar(c,1) = w_scalar(c,1) + (log(minpsi/psi)/log(minpsi/maxpsi))*fr(c,j)
            end if
         end do
      end do

#ifdef LCH4
      if (anoxia_wtsat) then ! Adjust for saturated fraction if unfrozen
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if (alt_indx(c) >= nlev_soildecomp_standard .and. t_soisno(c,1) > SHR_CONST_TKFRZ) then
               w_scalar(c,1) = w_scalar(c,1)*(1._r8 - finundated(c)) + finundated(c)
            end if
         end do
      end if
#endif 

#ifdef LCH4
      ! Calculate ANOXIA
      if (anoxia) then
      ! Check for anoxia w/o LCH4 now done in controlMod.

         do j = 1,nlev_soildecomp_standard
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (j==1) o_scalar(c,:) = 0._r8

               if (.not. anoxia_wtsat) then
                  o_scalar(c,1) = o_scalar(c,1) + fr(c,j) * max(o2stress_unsat(c,j), mino2lim)
               else
                  o_scalar(c,1) = o_scalar(c,1) + fr(c,j) * &
                       (max(o2stress_unsat(c,j), mino2lim)*(1._r8 - finundated(c)) + &
                       max(o2stress_sat(c,j), mino2lim)*finundated(c) )
               end if
            end do
         end do
      else
         o_scalar(lbc:ubc,1:nlevdecomp) = 1._r8
      end if
#else
      o_scalar(lbc:ubc,1:nlevdecomp) = 1._r8
#endif 

      deallocate(fr)

   else

      if ( .not. use_century_tfunc ) then
         ! calculate rate constant scalar for soil temperature
         ! assuming that the base rate constants are assigned for non-moisture
         ! limiting conditions at 25 C. 
         ! Peter Thornton: 3/13/09
         ! Replaced the Lloyd and Taylor function with a Q10 formula, with Q10 = 1.5
         ! as part of the modifications made to improve the seasonal cycle of 
         ! atmospheric CO2 concentration in global simulations. This does not impact
         ! the base rates at 25 C, which are calibrated from microcosm studies.
         
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               !! use separate (possibly equal) t funcs above and below freezing point
               !! t_scalar(c,j)= (q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))
               if (t_soisno(c,j) .ge. SHR_CONST_TKFRZ) then
                  t_scalar(c,j)= (q10**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))
               else
                  t_scalar(c,j)= (q10**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))
               endif
            end do
         end do
         
      else
         
         do j = 1, nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               t_scalar(c,j)= max(catanf(t_soisno(c,j)-SHR_CONST_TKFRZ)/catanf_30, 0.01_r8)
            end do
         end do

      endif

      ! calculate the rate constant scalar for soil water content.
      ! Uses the log relationship with water potential given in
      ! Andren, O., and K. Paustian, 1987. Barley straw decomposition in the field:
      ! a comparison of models. Ecology, 68(5):1190-1200.
      ! and supported by data in
      ! Orchard, V.A., and F.J. Cook, 1983. Relationship between soil respiration
      ! and soil moisture. Soil Biol. Biochem., 15(4):447-453.
      
      minpsi = -10.0_r8;
      do j = 1,nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            maxpsi = sucsat(c,j) * (-9.8e-6_r8)
            psi = min(soilpsi(c,j),maxpsi)
            ! decomp only if soilpsi is higher than minpsi
            if (psi > minpsi) then
               w_scalar(c,j) = (log(minpsi/psi)/log(minpsi/maxpsi))
            else
               w_scalar(c,j) = 0._r8
            end if
#ifdef LCH4
            if (anoxia_wtsat .and. t_soisno(c,j) > SHR_CONST_TKFRZ) then ! wet area will have w_scalar of 1 if unfrozen
               w_scalar(c,j) = w_scalar(c,j)*(1._r8 - finundated(c)) + finundated(c)
            end if
#endif
         end do
      end do

#ifdef LCH4
      ! Calculate ANOXIA
      ! Check for anoxia w/o LCH4 now done in controlMod.

      if (anoxia) then
         do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if (.not. anoxia_wtsat) then
                  o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim)
               else
                  o_scalar(c,j) = max(o2stress_unsat(c,j), mino2lim) * (1._r8 - finundated(c)) + &
                       max(o2stress_sat(c,j), mino2lim) * finundated(c)
               end if
            end do
         end do
      else
         o_scalar(lbc:ubc,1:nlevdecomp) = 1._r8
      end if
#else
      o_scalar(lbc:ubc,1:nlevdecomp) = 1._r8
#endif

   end if

   if ( normalize_q10_to_century_tfunc ) then
      ! scale all decomposition rates by a constant to compensate for offset between original CENTURY temp func and Q10
      normalization_factor = (catanf(normalization_tref)/catanf_30) / (q10**((normalization_tref-25._r8)/10._r8))
      do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            t_scalar(c,j) = t_scalar(c,j) * normalization_factor
         end do
      end do
   endif
   
#if (defined VERTSOILC)
   ! add a term to reduce decomposition rate at depth
   ! for now used a fixed e-folding depth
   do j = 1, nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         depth_scalar(c,j) = exp(-zsoi(j)/decomp_depth_efolding)
      end do
   end do
#endif
   
#if (defined VERTSOILC)
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_litr2) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_litr3) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_soil2) = k_s2 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_soil3) = k_s3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j)
      end do
   end do
#else
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_litr2) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_litr3) = k_l2_l3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_soil2) = k_s2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
         decomp_k(c,j,i_soil3) = k_s3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j)
      end do
   end do
#endif


end subroutine decomp_rate_constants

#endif

#endif

end module CNDecompCascadeMod_CENTURY

