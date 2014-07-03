module CNDecompCascadeMod_BGC
#ifdef CN

#ifndef CENTURY_DECOMP

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNDecompMod
!
! !DESCRIPTION:
! Module that sets the coeffiecients used in the decomposition cascade submodel.  This uses the BGC parameters as in CLMCN 4.0
!
! !USES:
   use shr_kind_mod , only: r8 => shr_kind_r8
   use shr_const_mod, only: SHR_CONST_TKFRZ
   use clm_varpar   , only: nlevsoi, nlevgrnd, nlevdecomp, ndecomp_cascade_transitions, ndecomp_pools, nsompools
   use clm_varpar   , only: i_met_lit, i_cel_lit, i_lig_lit, i_cwd
   use clm_varctl  , only: iulog, spinup_state
   use clm_varcon, only: zsoi
#ifdef LCH4
   use clm_varctl, only: anoxia
   use ch4varcon, only: mino2lim
#endif

   implicit none
   save
   private
! !PUBLIC MEMBER FUNCTIONS:
   public:: init_decompcascade, decomp_rate_constants
!
! !PUBLIC DATA MEMBERS:
#if (defined VERTSOILC)
   real(r8), public :: decomp_depth_efolding = 0.5_r8    ! (meters) e-folding depth for reduction in decomposition [set to large number for depth-independance]
#endif
   real(r8), public :: froz_q10 = 1.5_r8                 ! separate q10 for frozen soil respiration rates.  default to same as above zero rates
#ifdef LCH4
   logical,  public :: anoxia_wtsat = .false.            ! true ==> weight anoxia by inundated fraction
#endif
   integer, public :: nlev_soildecomp_standard           ! used here and in ch4Mod

   !! parameters for AD spinup
   real(r8), public, parameter :: spinup_vector(nsompools) = (/ 1.0_r8, 1.0_r8, 5.0_r8, 70.0_r8 /) ! multipliers for soil decomp during accelerated spinup

!EOP
!-----------------------------------------------------------------------

contains


!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: init_decompcascade
!
! !INTERFACE:
subroutine init_decompcascade(begc, endc)

! !DESCRIPTION:
!
!  initialize rate constants and decomposition pathways for the BGC model originally implemented in CLM-CN
!  written by C. Koven based on original CLM4 decomposition cascade by P. Thornton
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
   character(len=8), pointer :: cascade_step_name(:)        ! name of transition
   real(r8), pointer :: rf_decomp_cascade(:,:,:)            ! respired fraction in decomposition step (frac)
   integer,  pointer :: cascade_donor_pool(:)               ! which pool is C taken from for a given decomposition step
   integer,  pointer :: cascade_receiver_pool(:)            ! which pool is C added to for a given decomposition step
   real(r8), pointer :: pathfrac_decomp_cascade(:,:,:)      ! what fraction of C leaving a given pool passes through a given transition (frac)
   !-- properties of each decomposing pool
   logical,  pointer :: floating_cn_ratio_decomp_pools(:)   ! TRUE => pool has fixed C:N ratio
   character(len=8), pointer :: decomp_pool_name_restart(:) ! name of pool for restart files
   character(len=8), pointer :: decomp_pool_name_history(:) ! name of pool for history files
   character(len=20), pointer :: decomp_pool_name_long(:)   ! name of pool for netcdf long names
   character(len=8), pointer :: decomp_pool_name_short(:)   ! name of pool for netcdf short names
   logical, pointer :: is_litter(:)                         ! TRUE => pool is a litter pool
   logical, pointer :: is_soil(:)                           ! TRUE => pool is a soil pool
   logical, pointer :: is_cwd(:)                            ! TRUE => pool is a cwd pool
   real(r8), pointer :: initial_cn_ratio(:)                 ! c:n ratio for initialization of pools
   real(r8), pointer :: initial_stock(:)                    ! initial concentration for seeding at spinup
   logical, pointer :: is_metabolic(:)                      ! TRUE => pool is metabolic material
   logical, pointer :: is_cellulose(:)                      ! TRUE => pool is cellulose
   logical, pointer :: is_lignin(:)                         ! TRUE => pool is lignin
   real(r8), pointer :: spinup_factor(:)      ! factor for AD spinup associated with each pool
   real(r8):: rf_l1s1      !respiration fraction litter 1 -> SOM 1
   real(r8):: rf_l2s2      !respiration fraction litter 2 -> SOM 2
   real(r8):: rf_l3s3      !respiration fraction litter 3 -> SOM 3
   real(r8):: rf_s1s2      !respiration fraction SOM 1 -> SOM 2
   real(r8):: rf_s2s3      !respiration fraction SOM 2 -> SOM 3
   real(r8):: rf_s3s4      !respiration fraction SOM 3 -> SOM 4
   real(r8):: cwd_fcel
   real(r8):: cwd_flig
   real(r8) :: cn_s1
   real(r8) :: cn_s2
   real(r8) :: cn_s3
   real(r8) :: cn_s4

   integer :: i_litr1
   integer :: i_litr2
   integer :: i_litr3
   integer :: i_soil1
   integer :: i_soil2
   integer :: i_soil3
   integer :: i_soil4
   integer :: i_atm
   integer :: i_l1s1
   integer :: i_l2s2
   integer :: i_l3s3
   integer :: i_s1s2
   integer :: i_s2s3
   integer :: i_s3s4
   integer :: i_s4atm
   integer :: i_cwdl2
   integer :: i_cwdl3



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
   spinup_factor                           => decomp_cascade_con%spinup_factor


   !------- time-constant coefficients ---------- !
   ! set soil organic matter compartment C:N ratios (from Biome-BGC v4.2.0)
   cn_s1 = 12.0_r8
   cn_s2 = 12.0_r8
   cn_s3 = 10.0_r8
   cn_s4 = 10.0_r8

   ! set respiration fractions for fluxes between compartments
   ! (from Biome-BGC v4.2.0)
   rf_l1s1 = 0.39_r8
   rf_l2s2 = 0.55_r8
   rf_l3s3 = 0.29_r8
   rf_s1s2 = 0.28_r8
   rf_s2s3 = 0.46_r8
   rf_s3s4 = 0.55_r8

   ! set the cellulose and lignin fractions for coarse woody debris
   cwd_fcel = 0.76_r8
   cwd_flig = 0.24_r8

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

   floating_cn_ratio_decomp_pools(i_cwd) = .true.
   decomp_pool_name_restart(i_cwd) = 'cwd'
   decomp_pool_name_history(i_cwd) = 'CWD'
   decomp_pool_name_long(i_cwd) = 'coarse woody debris'
   decomp_pool_name_short(i_cwd) = 'CWD'
   is_litter(i_cwd) = .false.
   is_soil(i_cwd) = .false.
   is_cwd(i_cwd) = .true.
   initial_cn_ratio(i_cwd) = 500._r8
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
   initial_stock(i_soil1) = 0._r8
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
   initial_stock(i_soil2) = 0._r8
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
   initial_stock(i_soil3) = 0._r8
   is_metabolic(i_soil3) = .false.
   is_cellulose(i_soil3) = .false.
   is_lignin(i_soil3) = .false.

   i_soil4 = 8
   floating_cn_ratio_decomp_pools(i_soil4) = .false.
   decomp_pool_name_restart(i_soil4) = 'soil4'
   decomp_pool_name_history(i_soil4) = 'SOIL4'
   decomp_pool_name_long(i_soil4) = 'soil 4'
   decomp_pool_name_short(i_soil4) = 'S4'
   is_litter(i_soil4) = .false.
   is_soil(i_soil4) = .true.
   is_cwd(i_soil4) = .false.
   initial_cn_ratio(i_soil4) = cn_s4
   initial_stock(i_soil4) = 10._r8
   is_metabolic(i_soil4) = .false.
   is_cellulose(i_soil4) = .false.
   is_lignin(i_soil4) = .false.

   i_atm = 0  !! for terminal pools (i.e. 100% respiration)
   floating_cn_ratio_decomp_pools(i_atm) = .false.
   decomp_pool_name_restart(i_atm) = 'atmosphere'
   decomp_pool_name_history(i_atm) = 'atmosphere'
   decomp_pool_name_long(i_atm) = 'atmosphere'
   decomp_pool_name_short(i_atm) = ''
   is_litter(i_atm) = .true.
   is_soil(i_atm) = .false.
   is_cwd(i_atm) = .false.
   initial_cn_ratio(i_atm) = 0._r8
   initial_stock(i_atm) = 0._r8
   is_metabolic(i_atm) = .false.
   is_cellulose(i_atm) = .false.
   is_lignin(i_atm) = .false.


   spinup_factor(i_litr1) = 1._r8
   spinup_factor(i_litr2) = 1._r8
   spinup_factor(i_litr3) = 1._r8
   spinup_factor(i_cwd) = 1._r8
   spinup_factor(i_soil1) = spinup_vector(1)
   spinup_factor(i_soil2) = spinup_vector(2)
   spinup_factor(i_soil3) = spinup_vector(3)
   spinup_factor(i_soil4) = spinup_vector(4)


   
   !----------------  list of transitions and their time-independent coefficients  ---------------!
   i_l1s1 = 1
   cascade_step_name(i_l1s1) = 'L1S1'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_l1s1) = rf_l1s1
   cascade_donor_pool(i_l1s1) = i_litr1
   cascade_receiver_pool(i_l1s1) = i_soil1
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_l1s1) = 1.0_r8

   i_l2s2 = 2
   cascade_step_name(i_l2s2) = 'L2S2'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_l2s2) = rf_l2s2
   cascade_donor_pool(i_l2s2) = i_litr2
   cascade_receiver_pool(i_l2s2) = i_soil2
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_l2s2) = 1.0_r8

   i_l3s3 = 3
   cascade_step_name(i_l3s3) = 'L3S3'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_l3s3) = rf_l3s3
   cascade_donor_pool(i_l3s3) = i_litr3
   cascade_receiver_pool(i_l3s3) = i_soil3
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_l3s3) = 1.0_r8

   i_s1s2 = 4
   cascade_step_name(i_s1s2) = 'S1S2'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s1s2) = rf_s1s2
   cascade_donor_pool(i_s1s2) = i_soil1
   cascade_receiver_pool(i_s1s2) = i_soil2
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s1s2) = 1.0_r8

   i_s2s3 = 5
   cascade_step_name(i_s2s3) = 'S2S3'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s2s3) = rf_s2s3
   cascade_donor_pool(i_s2s3) = i_soil2
   cascade_receiver_pool(i_s2s3) = i_soil3
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s2s3) = 1.0_r8

   i_s3s4 = 6 
   cascade_step_name(i_s3s4) = 'S3S4'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s3s4) = rf_s3s4
   cascade_donor_pool(i_s3s4) = i_soil3
   cascade_receiver_pool(i_s3s4) = i_soil4
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s3s4) = 1.0_r8

   i_s4atm = 7
   cascade_step_name(i_s4atm) = 'S4'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_s4atm) = 1.
   cascade_donor_pool(i_s4atm) = i_soil4
   cascade_receiver_pool(i_s4atm) = i_atm
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_s4atm) = 1.0_r8

   i_cwdl2 = 8
   cascade_step_name(i_cwdl2) = 'CWDL2'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl2) = 0._r8
   cascade_donor_pool(i_cwdl2) = i_cwd
   cascade_receiver_pool(i_cwdl2) = i_litr2
   pathfrac_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl2) = cwd_fcel

   i_cwdl3 = 9
   cascade_step_name(i_cwdl3) = 'CWDL3'
   rf_decomp_cascade(begc:endc,1:nlevdecomp,i_cwdl3) = 0._r8
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
!  calculate rate constants and decomposition pathways for the BGC model originally implemented in CLM-CN
!  written by C. Koven based on original CLM4 decomposition cascade by P. Thornton
!
! !USES:
   use clmtype
   use clm_time_manager    , only : get_step_size
   use clm_varcon, only: secspday

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
   real(r8), pointer :: decomp_k(:,:,:)     ! rate constant for decomposition (1./sec)
   real(r8), pointer :: t_scalar(:,:)       ! soil temperature scalar for decomp
   real(r8), pointer :: w_scalar(:,:)       ! soil water scalar for decomp
   real(r8), pointer :: o_scalar(:,:)       ! fraction by which decomposition is limited by anoxia

   real(r8), pointer :: dz(:,:)             ! soil layer thickness (m)
   real(r8), pointer :: t_soisno(:,:)       ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: sucsat(:,:)        ! minimum soil suction (mm)
   real(r8), pointer :: soilpsi(:,:)        ! soil water potential in each soil layer (MPa)
#ifdef LCH4
   real(r8), pointer :: o2stress_unsat(:,:) ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   real(r8), pointer :: o2stress_sat(:,:)   ! Ratio of oxygen available to that demanded by roots, aerobes, & methanotrophs (nlevsoi)
   real(r8), pointer :: finundated(:)       ! fractional inundated area (excluding dedicated wetland columns)
#endif
   integer, pointer :: alt_indx(:)          ! current depth of thaw

   real(r8) :: dt                           ! decomp timestep (seconds)   
   real(r8):: dtd                           ! decomp timestep (days)
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   real(r8):: frw(lbc:ubc)                  ! rooting fraction weight
   real(r8), pointer:: fr(:,:)              ! column-level rooting fraction by soil depth
   real(r8):: minpsi, maxpsi                ! limits for soil water scalar for decomp
   real(r8):: psi                           ! temporary soilpsi for water scalar
   !   real(r8):: w_scalar(lbc:ubc,1:nlevdecomp)     !soil water scalar for decomp
   real(r8):: rate_scalar  ! combined rate scalar for decomp
   real(r8):: k_l1         ! decomposition rate constant litter 1
   real(r8):: k_l2         ! decomposition rate constant litter 2
   real(r8):: k_l3         ! decomposition rate constant litter 3
   real(r8):: k_s1         ! decomposition rate constant SOM 1
   real(r8):: k_s2         ! decomposition rate constant SOM 2
   real(r8):: k_s3         ! decomposition rate constant SOM 3
   real(r8):: k_s4         ! decomposition rate constant SOM 3
   real(r8):: k_frag       ! fragmentation rate constant CWD
   real(r8):: ck_l1        ! corrected decomposition rate constant litter 1
   real(r8):: ck_l2        ! corrected decomposition rate constant litter 2
   real(r8):: ck_l3        ! corrected decomposition rate constant litter 3
   real(r8):: ck_s1        ! corrected decomposition rate constant SOM 1
   real(r8):: ck_s2        ! corrected decomposition rate constant SOM 2
   real(r8):: ck_s3        ! corrected decomposition rate constant SOM 3
   real(r8):: ck_s4        ! corrected decomposition rate constant SOM 3
   real(r8):: ck_frag      ! corrected fragmentation rate constant CWD
   real(r8):: cwd_fcel     ! cellulose fraction of coarse woody debris
   real(r8):: cwd_flig     ! lignin fraction of coarse woody debris
   real(r8):: cwdc_loss    ! fragmentation rate for CWD carbon (gC/m2/s)
   real(r8):: cwdn_loss    ! fragmentation rate for CWD nitrogen (gN/m2/s)

   integer :: i_litr1
   integer :: i_litr2
   integer :: i_litr3
   integer :: i_soil1
   integer :: i_soil2
   integer :: i_soil3
   integer :: i_soil4
   integer :: c, fc, j, k, l

#if (defined VERTSOILC)
   real(r8) :: depth_scalar(lbc:ubc,1:nlevdecomp) 
#endif

   ! Assign local pointers to derived type arrays
   t_soisno              => clm3%g%l%c%ces%t_soisno
   sucsat                => clm3%g%l%c%cps%sucsat
   soilpsi               => clm3%g%l%c%cps%soilpsi
   dz                    => clm3%g%l%c%cps%dz
   t_scalar              => clm3%g%l%c%ccf%t_scalar
   w_scalar              => clm3%g%l%c%ccf%w_scalar
   o_scalar              => clm3%g%l%c%ccf%o_scalar
   decomp_k              => clm3%g%l%c%ccf%decomp_k
#ifdef LCH4
   o2stress_sat          => clm3%g%l%c%cch4%o2stress_sat
   o2stress_unsat        => clm3%g%l%c%cch4%o2stress_unsat
   finundated            => clm3%g%l%c%cws%finundated
#endif
   alt_indx              => clm3%g%l%c%cps%alt_indx


   ! set time steps
   dt = real( get_step_size(), r8 )
   dtd = dt/secspday


   ! set initial base rates for decomposition mass loss (1/day)
   ! (from Biome-BGC v4.2.0, using three SOM pools)
   ! Value inside log function is the discrete-time values for a
   ! daily time step model, and the result of the log function is
   ! the corresponding continuous-time decay rate (1/day), following
   ! Olson, 1963.
   k_l1 = -log(1.0_r8-0.7_r8)
   k_l2 = -log(1.0_r8-0.07_r8)
   k_l3 = -log(1.0_r8-0.014_r8)
   k_s1 = -log(1.0_r8-0.07_r8)
   k_s2 = -log(1.0_r8-0.014_r8)
   k_s3 = -log(1.0_r8-0.0014_r8)
   k_s4 = -log(1.0_r8-0.0001_r8)
   k_frag = -log(1.0_r8-0.001_r8)

   ! calculate the new discrete-time decay rate for model timestep
   k_l1 = 1.0_r8-exp(-k_l1*dtd)
   k_l2 = 1.0_r8-exp(-k_l2*dtd)
   k_l3 = 1.0_r8-exp(-k_l3*dtd)
   k_s1 = 1.0_r8-exp(-k_s1*dtd)
   k_s2 = 1.0_r8-exp(-k_s2*dtd)
   k_s3 = 1.0_r8-exp(-k_s3*dtd)
   k_s4 = 1.0_r8-exp(-k_s4*dtd)
   k_frag = 1.0_r8-exp(-k_frag*dtd)
   
   ! The following code implements the acceleration part of the AD spinup
   ! algorithm, by multiplying all of the SOM decomposition base rates by 10.0.

if ( spinup_state .eq. 1 ) then
   k_s1 = k_s1 * spinup_vector(1)
   k_s2 = k_s2 * spinup_vector(2)
   k_s3 = k_s3 * spinup_vector(3)
   k_s4 = k_s4 * spinup_vector(4)
endif

   i_litr1 = 1
   i_litr2 = 2
   i_litr3 = 3
   i_soil1 = 5
   i_soil2 = 6
   i_soil3 = 7
   i_soil4 = 8



   ! ! CWD fragmentation -> litter pools
   !  thse have now been put into the regular decomposition cascade
   ! cwdc_loss = cwdc_vr(c,j) * ck_frag / dt
   ! cwdc_to_litr2c_vr(c,j) = cwdc_loss * cwd_fcel
   ! cwdc_to_litr3c_vr(c,j) = cwdc_loss * cwd_flig
   ! cwdn_loss = cwdn_vr(c,j) * ck_frag / dt
   ! cwdn_to_litr2n_vr(c,j) = cwdn_loss * cwd_fcel
   ! cwdn_to_litr3n_vr(c,j) = cwdn_loss * cwd_flig


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
      
      ! calculate rate constant scalar for soil temperature
      ! assuming that the base rate constants are assigned for non-moisture
      ! limiting conditions at 25 C. 
      ! Peter Thornton: 3/13/09
      ! Replaced the Lloyd and Taylor function with a Q10 formula, with Q10 = 1.5
      ! as part of the modifications made to improve the seasonal cycle of 
      ! atmospheric CO2 concentration in global simulations. This does not impact
      ! the base rates at 25 C, which are calibrated from microcosm studies.
      do j = 1,nlev_soildecomp_standard
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            if (j==1) t_scalar(c,:) = 0._r8
            !! use separate (possibly equal) t funcs above and below freezing point
            !! t_scalar(c,1)=t_scalar(c,1) + (1.5**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))*fr(c,j)
            if (t_soisno(c,j) .ge. SHR_CONST_TKFRZ) then
               t_scalar(c,1)=t_scalar(c,1) + (1.5**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))*fr(c,j)
            else
               t_scalar(c,1)=t_scalar(c,1) + (1.5**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))*fr(c,j)
            endif
         end do
      end do
      
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
      if (anoxia_wtsat) then ! Adjust for saturated fraction if unfrozen.
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
               !! t_scalar(c,j)= (1.5**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))
               if (t_soisno(c,j) .ge. SHR_CONST_TKFRZ) then
                  t_scalar(c,j)= (1.5**((t_soisno(c,j)-(SHR_CONST_TKFRZ+25._r8))/10._r8))
               else
                  t_scalar(c,j)= (1.5**(-25._r8/10._r8))*(froz_q10**((t_soisno(c,j)-SHR_CONST_TKFRZ)/10._r8))
               endif
         end do
      end do

      
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

   end if

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
         decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_litr2) = k_l2 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_litr3) = k_l3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil2) = k_s2 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil3) = k_s3 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil4) = k_s4 * t_scalar(c,j) * w_scalar(c,j) * depth_scalar(c,j) * o_scalar(c,j) / dt
      end do
   end do
#else
   do j = 1,nlevdecomp
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         decomp_k(c,j,i_litr1) = k_l1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_litr2) = k_l2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_litr3) = k_l3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_cwd) = k_frag * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil1) = k_s1 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil2) = k_s2 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil3) = k_s3 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
         decomp_k(c,j,i_soil4) = k_s4 * t_scalar(c,j) * w_scalar(c,j) * o_scalar(c,j) / dt
      end do
   end do
#endif



end subroutine decomp_rate_constants
#endif

#endif


end module CNDecompCascadeMod_BGC

