module CNPhenologyMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: CNPhenologyMod
!
! !DESCRIPTION:
! Module holding routines used in phenology model for coupled carbon
! nitrogen code.
!
! !USES:
  use clmtype
  use shr_kind_mod, only: r8 => shr_kind_r8
  use clm_varcon  , only: tfrz
  use clm_varctl  , only: iulog, use_cndv
  use clm_varpar  , only: numpft
  use shr_sys_mod , only: shr_sys_flush
  use abortutils  , only: endrun
  implicit none
  save
  private

! !PUBLIC MEMBER FUNCTIONS:
  public :: CNPhenologyInit      ! Initialization
  public :: CNPhenology          ! Update
!
! !REVISION HISTORY:
! 8/1/03: Created by Peter Thornton
! 10/23/03, Peter Thornton: migrated all routines to vector data structures
! 2/4/08,  slevis: adding crop phenology from AgroIBIS

! !PRIVATE DATA MEMBERS:

  real(r8)           :: dt              ! radiation time step delta t (seconds)
  real(r8)           :: fracday         ! dtime as a fraction of day
  real(r8)           :: crit_dayl       ! critical daylength for offset (seconds)
  real(r8)           :: ndays_on        ! number of days to complete onset
  real(r8)           :: ndays_off       ! number of days to complete offset
  real(r8)           :: fstor2tran      ! fraction of storage to move to transfer on each onset
  real(r8)           :: crit_onset_fdd  ! critical number of freezing days
  real(r8)           :: crit_onset_swi  ! water stress days for offset trigger
  real(r8)           :: soilpsi_on      ! water potential for onset trigger (MPa)
  real(r8)           :: crit_offset_fdd ! critical number of freezing degree days
                                        ! to trigger offset
  real(r8)           :: crit_offset_swi ! water stress days for offset trigger
  real(r8)           :: soilpsi_off     ! water potential for offset trigger (MPa)
  real(r8)           :: lwtop           ! live wood turnover proportion (annual fraction)
  !
  ! CropPhenology variables and constants
  !
  real(r8)           :: p1d, p1v            ! photoperiod factor constants for crop vernalization
  real(r8)           :: hti                 ! cold hardening index threshold for vernalization
  real(r8)           :: tbase               ! base temperature for vernalization
  integer, parameter :: NOT_Planted   = 999 ! If not planted   yet in year
  integer, parameter :: NOT_Harvested = 999 ! If not harvested yet in year
  integer, parameter :: inNH       = 1      ! Northern Hemisphere
  integer, parameter :: inSH       = 2      ! Southern Hemisphere
  integer, pointer   :: inhemi(:)           ! Hemisphere that pft is in
  integer            :: minplantjday(0:numpft,inSH) ! minimum planting julian day
  integer            :: maxplantjday(0:numpft,inSH) ! maximum planting julian day
  integer            :: jdayyrstart(inSH)           ! julian day of start of year

!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNPhenology
!
! !INTERFACE:
subroutine CNPhenology (num_soilc, filter_soilc, num_soilp, filter_soilp, &
                        num_pcropp, filter_pcropp, doalb)
!
! !DESCRIPTION:
! Dynamic phenology routine for coupled carbon-nitrogen code (CN)
! 1. grass phenology
!
! !USES:
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! filter for soil columns
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
   integer, intent(in) :: num_pcropp      ! number of prog. crop pfts in filter
   integer, intent(in) :: filter_pcropp(:)! filter for prognostic crop pfts
   logical, intent(in) :: doalb           ! true if time for sfc albedo calc
!
! !CALLED FROM:
! subroutine CNEcosystemDyn in module CNEcosystemDynMod.F90
!
! !REVISION HISTORY:
! 7/28/03: Created by Peter Thornton
! 9/05/03, Peter Thornton: moved from call with (p) to call with (c)
! 10/3/03, Peter Thornton: added subroutine calls for different phenology types
! 11/7/03, Peter Thornton: moved phenology type tests into phenology type
!    routines, and moved onset, offset, background litfall routines into
!    main phenology call.
! !LOCAL VARIABLES:
! local pointers to implicit in arrays
!
! local pointers to implicit in/out scalars
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
!EOP
!-----------------------------------------------------------------------

   ! each of the following phenology type routines includes a filter
   ! to operate only on the relevant pfts

   call CNPhenologyClimate(num_soilp, filter_soilp, num_pcropp, filter_pcropp)
   
   call CNEvergreenPhenology(num_soilp, filter_soilp)

   call CNSeasonDecidPhenology(num_soilp, filter_soilp)

   call CNStressDecidPhenology(num_soilp, filter_soilp)

   if (doalb .and. num_pcropp > 0 ) call CropPhenology(num_pcropp, filter_pcropp)

   ! the same onset and offset routines are called regardless of
   ! phenology type - they depend only on onset_flag, offset_flag, bglfr, and bgtr

   call CNOnsetGrowth(num_soilp, filter_soilp)

   call CNOffsetLitterfall(num_soilp, filter_soilp)

   call CNBackgroundLitterfall(num_soilp, filter_soilp)

   call CNLivewoodTurnover(num_soilp, filter_soilp)

   ! gather all pft-level litterfall fluxes to the column
   ! for litter C and N inputs

   call CNLitterToColumn(num_soilc, filter_soilc)

end subroutine CNPhenology

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNPhenologyInit
!
! !INTERFACE:
subroutine CNPhenologyInit( begp, endp )
!
! !DESCRIPTION:
! Initialization of CNPhenology. Must be called after time-manager is
! initialized, and after pftcon file is read in.
!
! !USES:
   use clm_time_manager, only: get_step_size
   use surfrdMod       , only: crop_prog
   use clm_varcon      , only: secspday
!
! !ARGUMENTS:
   implicit none
   integer, intent(IN) :: begp, endp ! Beginning and ending PFT index
! !CALLED FROM:
! subroutine initialize2 in module clm_initializeMod.F90
!
! !REVISION HISTORY:
! 3/28/11: Created by Erik Kluzek
!
! !LOCAL VARIABLES:
!EOP
!------------------------------------------------------------------------

    !
    ! Get time-step and what fraction of a day it is
    !
    dt      = real( get_step_size(), r8 )
    fracday = dt/secspday

    ! set some local parameters - these will be moved into
    ! parameter file after testing

    ! -----------------------------------------
    ! Constants for CNSeasonDecidPhenology
    ! -----------------------------------------
    !
    ! critical daylength from Biome-BGC, v4.1.2
    crit_dayl = 39300._r8

    ! -----------------------------------------
    ! Constants for CNSeasonDecidPhenology and CNStressDecidPhenology
    ! -----------------------------------------
    ndays_on  = 30._r8
    ndays_off = 15._r8

    ! transfer parameters
    fstor2tran = 0.5_r8
    ! -----------------------------------------
    ! Constants for CNStressDecidPhenology
    ! -----------------------------------------

    ! onset parameters
    crit_onset_fdd = 15.0_r8
    ! critical onset gdd now being calculated as a function of annual
    ! average 2m temp.
    ! crit_onset_gdd = 150.0 ! c3 grass value
    ! crit_onset_gdd = 1000.0   ! c4 grass value
    crit_onset_swi = 15.0_r8
    soilpsi_on     = -2.0_r8

    ! offset parameters
    crit_offset_fdd = 15.0_r8
    crit_offset_swi = 15.0_r8
    soilpsi_off     = -2.0_r8

    ! -----------------------------------------
    ! Constants for CNLivewoodTurnover
    ! -----------------------------------------

    ! set the global parameter for livewood turnover rate
    ! define as an annual fraction (0.7), and convert to fraction per second
    lwtop = 0.7_r8 / 31536000.0_r8

    ! -----------------------------------------
    ! Call any subroutine specific initialization routines
    ! -----------------------------------------

    if ( crop_prog ) call CropPhenologyInit( begp, endp )

end subroutine CNPhenologyInit
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNPhenologyClimate
!
! !INTERFACE:
subroutine CNPhenologyClimate (num_soilp, filter_soilp, num_pcropp, filter_pcropp)
!
! !DESCRIPTION:
! For coupled carbon-nitrogen code (CN).
!
! !USES:
   use clm_time_manager, only: get_days_per_year
   use clm_time_manager, only: get_curr_date, is_first_step
   use CropRestMod     , only: CropRestYear
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
   integer, intent(in) :: num_pcropp      ! number of prognostic crops in filter
   integer, intent(in) :: filter_pcropp(:)! filter for prognostic crop pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 3/13/07: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)             ! pft vegetation type
   ! ecophysiological constants
   real(r8), pointer :: t_ref2m(:)         ! 2m air temperature (K)
   real(r8), pointer :: tempavg_t2m(:)     ! temp. avg 2m air temperature (K)
   real(r8), pointer :: gdd0(:)            ! growing deg. days base 0 deg C (ddays)
   real(r8), pointer :: gdd8(:)            !    "     "    "    "   8  "  "    "
   real(r8), pointer :: gdd10(:)           !    "     "    "    "  10  "  "    "
   real(r8), pointer :: gdd020(:)          ! 20-yr mean of gdd0 (ddays)
   real(r8), pointer :: gdd820(:)          ! 20-yr mean of gdd8 (ddays)
   real(r8), pointer :: gdd1020(:)         ! 20-yr mean of gdd10 (ddays)
   integer , pointer :: pgridcell(:)       ! pft's gridcell index
!
! local pointers to implicit in/out scalars
!
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: p                    ! indices
   integer :: fp                   ! lake filter pft index
   integer :: nyrs                 ! number of years prognostic crop has run
   real(r8):: dayspyr              ! days per year (days)
   integer kyr                     ! current year
   integer kmo                     !         month of year  (1, ..., 12)
   integer kda                     !         day of month   (1, ..., 31)
   integer mcsec                   !         seconds of day (0, ..., seconds/day)
   real(r8), parameter :: yravg   = 20.0_r8      ! length of years to average for gdd
   real(r8), parameter :: yravgm1 = yravg-1.0_r8 ! minus 1 of above
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers to derived type arrays
   ivt                           => pft%itype
   t_ref2m                       => pes%t_ref2m
   tempavg_t2m                   => pepv%tempavg_t2m

   gdd0                          => pps%gdd0
   gdd8                          => pps%gdd8
   gdd10                         => pps%gdd10
   gdd020                        => pps%gdd020
   gdd820                        => pps%gdd820
   gdd1020                       => pps%gdd1020
   pgridcell                     => pft%gridcell

   ! set time steps

   dayspyr = get_days_per_year()

   do fp = 1,num_soilp
      p = filter_soilp(fp)
      tempavg_t2m(p) = tempavg_t2m(p) + t_ref2m(p) * (fracday/dayspyr)
   end do

   !
   ! The following crop related steps are done here rather than CropPhenology
   ! so that they will be completed each time-step rather than with doalb.
   !
   ! The following lines come from ibis's climate.f + stats.f
   ! gdd SUMMATIONS ARE RELATIVE TO THE PLANTING DATE (see subr. updateAccFlds)

   if (num_pcropp > 0) then
      ! get time-related info
      call get_curr_date(kyr, kmo, kda, mcsec)
      nyrs = CropRestYear()
   end if

   do fp = 1,num_pcropp
      p = filter_pcropp(fp)
      if (kmo == 1 .and. kda == 1 .and. nyrs == 0) then ! YR 1:
         gdd020(p)  = 0._r8                             ! set gdd..20 variables to 0
         gdd820(p)  = 0._r8                             ! and crops will not be planted
         gdd1020(p) = 0._r8
      end if
      if (kmo == 1 .and. kda == 1 .and. mcsec == 0) then        ! <-- END of EVERY YR:
         if (nyrs  == 1) then                                   ! <-- END of YR 1
            gdd020(p)  = gdd0(p)                                ! <-- END of YR 1
            gdd820(p)  = gdd8(p)                                ! <-- END of YR 1
            gdd1020(p) = gdd10(p)                               ! <-- END of YR 1
         end if                                                 ! <-- END of YR 1
         gdd020(p)  = (yravgm1* gdd020(p)  + gdd0(p))  / yravg  ! gdd..20 must be long term avgs
         gdd820(p)  = (yravgm1* gdd820(p)  + gdd8(p))  / yravg  ! so ignore results for yrs 1 & 2
         gdd1020(p) = (yravgm1* gdd1020(p) + gdd10(p)) / yravg 
      end if
   end do

end subroutine CNPhenologyClimate
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNEvergreenPhenology
!
! !INTERFACE:
subroutine CNEvergreenPhenology (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! For coupled carbon-nitrogen code (CN).
!
! !USES:
   use clm_varcon      , only: secspday
   use clm_time_manager, only: get_days_per_year
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/2/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)       ! pft vegetation type
   ! ecophysiological constants
   real(r8), pointer :: evergreen(:) ! binary flag for evergreen leaf habit (0 or 1)
   real(r8), pointer :: leaf_long(:) ! leaf longevity (yrs)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: bglfr(:)     ! background litterfall rate (1/s)
   real(r8), pointer :: bgtr(:)      ! background transfer growth rate (1/s)
   real(r8), pointer :: lgsf(:)      ! long growing season factor [0-1]
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   real(r8):: dayspyr                ! Days per year
   integer :: p                      ! indices
   integer :: fp                     ! lake filter pft index
!EOP
!-----------------------------------------------------------------------

   ! assign local pointers to derived type arrays
   ivt       => pft%itype
   evergreen => pftcon%evergreen
   leaf_long => pftcon%leaf_long
   bglfr     => pepv%bglfr
   bgtr      => pepv%bgtr
   lgsf      => pepv%lgsf
   dayspyr   = get_days_per_year()

   do fp = 1,num_soilp
      p = filter_soilp(fp)
      if (evergreen(ivt(p)) == 1._r8) then
          bglfr(p) = 1._r8/(leaf_long(ivt(p))*dayspyr*secspday)
          bgtr(p)  = 0._r8
          lgsf(p)  = 0._r8
      end if
   end do

end subroutine CNEvergreenPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNSeasonDecidPhenology
!
! !INTERFACE:
subroutine CNSeasonDecidPhenology (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! For coupled carbon-nitrogen code (CN).
! This routine handles the seasonal deciduous phenology code (temperate
! deciduous vegetation that has only one growing season per year).
!
! !USES:
   use shr_const_mod   , only: SHR_CONST_TKFRZ, SHR_CONST_PI
   use clm_varcon      , only: secspday
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/6/03: Created by Peter Thornton
! 10/24/03, Peter Thornton: migrated to vector data structures
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
   integer , pointer :: ivt(:)                ! pft vegetation type
   integer , pointer :: pcolumn(:)            ! pft's column index
   integer , pointer :: pgridcell(:)          ! pft's gridcell index
   real(r8), pointer :: latdeg(:)             ! latitude (radians)
   real(r8), pointer :: decl(:)               ! solar declination (radians)
   real(r8), pointer :: t_soisno(:,:)         ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   real(r8), pointer :: soilpsi(:,:)          ! soil water potential in each soil layer (MPa)
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   ! ecophysiological constants
   real(r8), pointer :: season_decid(:) ! binary flag for seasonal-deciduous leaf habit (0 or 1)
   real(r8), pointer :: woody(:)        ! binary flag for woody lifeform (1=woody, 0=not woody)
!
! local pointers to implicit in/out scalars
   real(r8), pointer :: dormant_flag(:)    ! dormancy flag
   real(r8), pointer :: days_active(:)     ! number of days since last dormancy
   real(r8), pointer :: onset_flag(:)      ! onset flag
   real(r8), pointer :: onset_counter(:)   ! onset counter (seconds)
   real(r8), pointer :: onset_gddflag(:)   ! onset freeze flag
   real(r8), pointer :: onset_gdd(:)       ! onset growing degree days
   real(r8), pointer :: offset_flag(:)     ! offset flag
   real(r8), pointer :: offset_counter(:)  ! offset counter (seconds)
   real(r8), pointer :: dayl(:)            ! daylength (seconds)
   real(r8), pointer :: prev_dayl(:)       ! daylength from previous albedo timestep (seconds)
   real(r8), pointer :: annavg_t2m(:)      ! annual average 2m air temperature (K)
   real(r8), pointer :: prev_leafc_to_litter(:)  ! previous timestep leaf C litterfall flux (gC/m2/s)
   real(r8), pointer :: prev_frootc_to_litter(:) ! previous timestep froot C litterfall flux (gC/m2/s)
   real(r8), pointer :: lgsf(:)            ! long growing season factor [0-1]
   real(r8), pointer :: bglfr(:)           ! background litterfall rate (1/s)
   real(r8), pointer :: bgtr(:)            ! background transfer growth rate (1/s)
   real(r8), pointer :: leafc_xfer_to_leafc(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: leafn_xfer_to_leafn(:)
   real(r8), pointer :: frootn_xfer_to_frootn(:)
   real(r8), pointer :: livestemn_xfer_to_livestemn(:)
   real(r8), pointer :: deadstemn_xfer_to_deadstemn(:)
   real(r8), pointer :: livecrootn_xfer_to_livecrootn(:)
   real(r8), pointer :: deadcrootn_xfer_to_deadcrootn(:)
   real(r8), pointer :: leafc_xfer(:)      ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc_xfer(:)     ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc_xfer(:)  ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc_xfer(:)  ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc_xfer(:) ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc_xfer(:) ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: leafn_xfer(:)      ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn_xfer(:)     ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn_xfer(:)  ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn_xfer(:)  ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn_xfer(:) ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn_xfer(:) ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: leafc_storage_to_xfer(:)
   real(r8), pointer :: frootc_storage_to_xfer(:)
   real(r8), pointer :: livestemc_storage_to_xfer(:)
   real(r8), pointer :: deadstemc_storage_to_xfer(:)
   real(r8), pointer :: livecrootc_storage_to_xfer(:)
   real(r8), pointer :: deadcrootc_storage_to_xfer(:)
   real(r8), pointer :: gresp_storage_to_xfer(:)
   real(r8), pointer :: leafn_storage_to_xfer(:)
   real(r8), pointer :: frootn_storage_to_xfer(:)
   real(r8), pointer :: livestemn_storage_to_xfer(:)
   real(r8), pointer :: deadstemn_storage_to_xfer(:)
   real(r8), pointer :: livecrootn_storage_to_xfer(:)
   real(r8), pointer :: deadcrootn_storage_to_xfer(:)
   logical , pointer :: pftmayexist(:)     ! exclude seasonal decid pfts from tropics
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: c,p            !indices
   integer :: fp             !lake filter pft index
   real(r8):: ws_flag        !winter-summer solstice flag (0 or 1)
   real(r8):: crit_onset_gdd !critical onset growing degree-day sum
   real(r8):: soilt
   real(r8):: lat            !latitude (radians)
   real(r8):: temp           !temporary variable for daylength calculation

!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
   ivt                           => pft%itype
   pcolumn                       => pft%column
   pgridcell                     => pft%gridcell
   latdeg                        => grc%latdeg
   decl                          => cps%decl
   t_soisno                      => ces%t_soisno
   leafc_storage                 => pcs%leafc_storage
   frootc_storage                => pcs%frootc_storage
   livestemc_storage             => pcs%livestemc_storage
   deadstemc_storage             => pcs%deadstemc_storage
   livecrootc_storage            => pcs%livecrootc_storage
   deadcrootc_storage            => pcs%deadcrootc_storage
   gresp_storage                 => pcs%gresp_storage
   leafn_storage                 => pns%leafn_storage
   frootn_storage                => pns%frootn_storage
   livestemn_storage             => pns%livestemn_storage
   deadstemn_storage             => pns%deadstemn_storage
   livecrootn_storage            => pns%livecrootn_storage
   deadcrootn_storage            => pns%deadcrootn_storage
   season_decid                  => pftcon%season_decid
   woody                         => pftcon%woody

   ! Assign local pointers to derived type arrays (out)
   dormant_flag                  => pepv%dormant_flag
   days_active                   => pepv%days_active
   onset_flag                    => pepv%onset_flag
   onset_counter                 => pepv%onset_counter
   onset_gddflag                 => pepv%onset_gddflag
   onset_gdd                     => pepv%onset_gdd
   offset_flag                   => pepv%offset_flag
   offset_counter                => pepv%offset_counter
   dayl                          => pepv%dayl
   prev_dayl                     => pepv%prev_dayl
   annavg_t2m                    => pepv%annavg_t2m
   prev_leafc_to_litter          => pepv%prev_leafc_to_litter
   prev_frootc_to_litter         => pepv%prev_frootc_to_litter
   bglfr                         => pepv%bglfr
   bgtr                          => pepv%bgtr
   lgsf                          => pepv%lgsf
   leafc_xfer_to_leafc           => pcf%leafc_xfer_to_leafc
   frootc_xfer_to_frootc         => pcf%frootc_xfer_to_frootc
   livestemc_xfer_to_livestemc   => pcf%livestemc_xfer_to_livestemc
   deadstemc_xfer_to_deadstemc   => pcf%deadstemc_xfer_to_deadstemc
   livecrootc_xfer_to_livecrootc => pcf%livecrootc_xfer_to_livecrootc
   deadcrootc_xfer_to_deadcrootc => pcf%deadcrootc_xfer_to_deadcrootc
   leafn_xfer_to_leafn           => pnf%leafn_xfer_to_leafn
   frootn_xfer_to_frootn         => pnf%frootn_xfer_to_frootn
   livestemn_xfer_to_livestemn   => pnf%livestemn_xfer_to_livestemn
   deadstemn_xfer_to_deadstemn   => pnf%deadstemn_xfer_to_deadstemn
   livecrootn_xfer_to_livecrootn => pnf%livecrootn_xfer_to_livecrootn
   deadcrootn_xfer_to_deadcrootn => pnf%deadcrootn_xfer_to_deadcrootn
   leafc_xfer                    => pcs%leafc_xfer
   frootc_xfer                   => pcs%frootc_xfer
   livestemc_xfer                => pcs%livestemc_xfer
   deadstemc_xfer                => pcs%deadstemc_xfer
   livecrootc_xfer               => pcs%livecrootc_xfer
   deadcrootc_xfer               => pcs%deadcrootc_xfer
   leafn_xfer                    => pns%leafn_xfer
   frootn_xfer                   => pns%frootn_xfer
   livestemn_xfer                => pns%livestemn_xfer
   deadstemn_xfer                => pns%deadstemn_xfer
   livecrootn_xfer               => pns%livecrootn_xfer
   deadcrootn_xfer               => pns%deadcrootn_xfer
   leafc_storage_to_xfer         => pcf%leafc_storage_to_xfer
   frootc_storage_to_xfer        => pcf%frootc_storage_to_xfer
   livestemc_storage_to_xfer     => pcf%livestemc_storage_to_xfer
   deadstemc_storage_to_xfer     => pcf%deadstemc_storage_to_xfer
   livecrootc_storage_to_xfer    => pcf%livecrootc_storage_to_xfer
   deadcrootc_storage_to_xfer    => pcf%deadcrootc_storage_to_xfer
   gresp_storage_to_xfer         => pcf%gresp_storage_to_xfer
   leafn_storage_to_xfer         => pnf%leafn_storage_to_xfer
   frootn_storage_to_xfer        => pnf%frootn_storage_to_xfer
   livestemn_storage_to_xfer     => pnf%livestemn_storage_to_xfer
   deadstemn_storage_to_xfer     => pnf%deadstemn_storage_to_xfer
   livecrootn_storage_to_xfer    => pnf%livecrootn_storage_to_xfer
   deadcrootn_storage_to_xfer    => pnf%deadcrootn_storage_to_xfer
   pftmayexist                   => pdgvs%pftmayexist

   ! start pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)

      if (season_decid(ivt(p)) == 1._r8) then

         ! set background litterfall rate, background transfer rate, and
         ! long growing season factor to 0 for seasonal deciduous types
         bglfr(p) = 0._r8
         bgtr(p) = 0._r8
         lgsf(p) = 0._r8

         ! onset gdd sum from Biome-BGC, v4.1.2
         crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_t2m(p) - SHR_CONST_TKFRZ))

         ! use solar declination information stored during Surface Albedo()
         ! and latitude from gps to calcluate daylength (convert latitude from degrees to radians)
         ! the constant 13750.9871 is the number of seconds per radian of hour-angle

         prev_dayl(p) = dayl(p)
         lat = (SHR_CONST_PI/180._r8)*latdeg(pgridcell(p))
         temp = -(sin(lat)*sin(decl(c)))/(cos(lat) * cos(decl(c)))
         temp = min(1._r8,max(-1._r8,temp))
         dayl(p) = 2.0_r8 * 13750.9871_r8 * acos(temp)

         ! set flag for solstice period (winter->summer = 1, summer->winter = 0)
         if (dayl(p) >= prev_dayl(p)) then
            ws_flag = 1._r8
         else
            ws_flag = 0._r8
         end if

         ! update offset_counter and test for the end of the offset period
         if (offset_flag(p) == 1.0_r8) then
            ! decrement counter for offset period
            offset_counter(p) = offset_counter(p) - dt

            ! if this is the end of the offset_period, reset phenology
            ! flags and indices
            if (offset_counter(p) == 0.0_r8) then
               ! this code block was originally handled by call cn_offset_cleanup(p)
               ! inlined during vectorization

               offset_flag(p) = 0._r8
               offset_counter(p) = 0._r8
               dormant_flag(p) = 1._r8
               days_active(p) = 0._r8
               if (use_cndv) then
                  pftmayexist(p) = .true.
               end if

               ! reset the previous timestep litterfall flux memory
               prev_leafc_to_litter(p) = 0._r8
               prev_frootc_to_litter(p) = 0._r8
            end if
         end if

         ! update onset_counter and test for the end of the onset period
         if (onset_flag(p) == 1.0_r8) then
            ! decrement counter for onset period
            onset_counter(p) = onset_counter(p) - dt

            ! if this is the end of the onset period, reset phenology
            ! flags and indices
            if (onset_counter(p) == 0.0_r8) then
               ! this code block was originally handled by call cn_onset_cleanup(p)
               ! inlined during vectorization

               onset_flag(p) = 0.0_r8
               onset_counter(p) = 0.0_r8
               ! set all transfer growth rates to 0.0
               leafc_xfer_to_leafc(p)   = 0.0_r8
               frootc_xfer_to_frootc(p) = 0.0_r8
               leafn_xfer_to_leafn(p)   = 0.0_r8
               frootn_xfer_to_frootn(p) = 0.0_r8
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_xfer_to_livestemc(p)   = 0.0_r8
                  deadstemc_xfer_to_deadstemc(p)   = 0.0_r8
                  livecrootc_xfer_to_livecrootc(p) = 0.0_r8
                  deadcrootc_xfer_to_deadcrootc(p) = 0.0_r8
                  livestemn_xfer_to_livestemn(p)   = 0.0_r8
                  deadstemn_xfer_to_deadstemn(p)   = 0.0_r8
                  livecrootn_xfer_to_livecrootn(p) = 0.0_r8
                  deadcrootn_xfer_to_deadcrootn(p) = 0.0_r8
               end if
               ! set transfer pools to 0.0
               leafc_xfer(p) = 0.0_r8
               leafn_xfer(p) = 0.0_r8
               frootc_xfer(p) = 0.0_r8
               frootn_xfer(p) = 0.0_r8
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_xfer(p) = 0.0_r8
                  livestemn_xfer(p) = 0.0_r8
                  deadstemc_xfer(p) = 0.0_r8
                  deadstemn_xfer(p) = 0.0_r8
                  livecrootc_xfer(p) = 0.0_r8
                  livecrootn_xfer(p) = 0.0_r8
                  deadcrootc_xfer(p) = 0.0_r8
                  deadcrootn_xfer(p) = 0.0_r8
               end if
            end if
         end if

         ! test for switching from dormant period to growth period
         if (dormant_flag(p) == 1.0_r8) then

            ! Test to turn on growing degree-day sum, if off.
            ! switch on the growing degree day sum on the winter solstice

            if (onset_gddflag(p) == 0._r8 .and. ws_flag == 1._r8) then
               onset_gddflag(p) = 1._r8
               onset_gdd(p) = 0._r8
            end if

            ! Test to turn off growing degree-day sum, if on.
            ! This test resets the growing degree day sum if it gets past
            ! the summer solstice without reaching the threshold value.
            ! In that case, it will take until the next winter solstice
            ! before the growing degree-day summation starts again.

            if (onset_gddflag(p) == 1._r8 .and. ws_flag == 0._r8) then
               onset_gddflag(p) = 0._r8
               onset_gdd(p) = 0._r8
            end if

            ! if the gdd flag is set, and if the soil is above freezing
            ! then accumulate growing degree days for onset trigger

            soilt = t_soisno(c,3)
            if (onset_gddflag(p) == 1.0_r8 .and. soilt > SHR_CONST_TKFRZ) then
               onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
            end if

            ! set onset_flag if critical growing degree-day sum is exceeded
            if (onset_gdd(p) > crit_onset_gdd) then
               onset_flag(p) = 1.0_r8
               dormant_flag(p) = 0.0_r8
               onset_gddflag(p) = 0.0_r8
               onset_gdd(p) = 0.0_r8
               onset_counter(p) = ndays_on * secspday

               ! move all the storage pools into transfer pools,
               ! where they will be transfered to displayed growth over the onset period.
               ! this code was originally handled with call cn_storage_to_xfer(p)
               ! inlined during vectorization

               ! set carbon fluxes for shifting storage pools to transfer pools
               leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
               frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
                  deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
                  livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
                  deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
                  gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
               end if

               ! set nitrogen fluxes for shifting storage pools to transfer pools
               leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
               frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
                  deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
                  livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
                  deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
               end if
            end if

         ! test for switching from growth period to offset period
         else if (offset_flag(p) == 0.0_r8) then
            if (use_cndv) then
               ! If days_active > 355, then remove pft in
               ! CNDVEstablishment at the end of the year.
               ! days_active > 355 is a symptom of seasonal decid. pfts occurring in
               ! gridcells where dayl never drops below crit_dayl.
               ! This results in TLAI>1e4 in a few gridcells.
               days_active(p) = days_active(p) + fracday
               if (days_active(p) > 355._r8) pftmayexist(p) = .false.
            end if

            ! only begin to test for offset daylength once past the summer sol
            if (ws_flag == 0._r8 .and. dayl(p) < crit_dayl) then
               offset_flag(p) = 1._r8
               offset_counter(p) = ndays_off * secspday
               prev_leafc_to_litter(p) = 0._r8
               prev_frootc_to_litter(p) = 0._r8
            end if
         end if

      end if ! end if seasonal deciduous

   end do ! end of pft loop

end subroutine CNSeasonDecidPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNStressDecidPhenology
!
! !INTERFACE:
subroutine CNStressDecidPhenology (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! This routine handles phenology for vegetation types, such as grasses and
! tropical drought deciduous trees, that respond to cold and drought stress
! signals and that can have multiple growing seasons in a given year.
! This routine allows for the possibility that leaves might persist year-round
! in the absence of a suitable stress trigger, by switching to an essentially
! evergreen habit, but maintaining a deciduous leaf longevity, while waiting
! for the next stress trigger.  This is in contrast to the seasonal deciduous
! algorithm (for temperate deciduous trees) that forces a single growing season
! per year.
!
! !USES:
   use clm_time_manager, only: get_days_per_year
   use clm_varcon      , only: secspday
   use shr_const_mod   , only: SHR_CONST_TKFRZ, SHR_CONST_PI
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/27/03: Created by Peter Thornton
! 01/29/04: Made onset_gdd critical sum a function of temperature, as in
!           seasonal deciduous algorithm.
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)                ! pft vegetation type
   integer , pointer :: pcolumn(:)            ! pft's column index
   integer , pointer :: pgridcell(:)          ! pft's gridcell index
   real(r8), pointer :: latdeg(:)             ! latitude (radians)
   real(r8), pointer :: decl(:)               ! solar declination (radians)
   real(r8), pointer :: leafc_storage(:)      ! (gC/m2) leaf C storage
   real(r8), pointer :: frootc_storage(:)     ! (gC/m2) fine root C storage
   real(r8), pointer :: livestemc_storage(:)  ! (gC/m2) live stem C storage
   real(r8), pointer :: deadstemc_storage(:)  ! (gC/m2) dead stem C storage
   real(r8), pointer :: livecrootc_storage(:) ! (gC/m2) live coarse root C storage
   real(r8), pointer :: deadcrootc_storage(:) ! (gC/m2) dead coarse root C storage
   real(r8), pointer :: gresp_storage(:)      ! (gC/m2) growth respiration storage
   real(r8), pointer :: leafn_storage(:)      ! (gN/m2) leaf N storage
   real(r8), pointer :: frootn_storage(:)     ! (gN/m2) fine root N storage
   real(r8), pointer :: livestemn_storage(:)  ! (gN/m2) live stem N storage
   real(r8), pointer :: deadstemn_storage(:)  ! (gN/m2) dead stem N storage
   real(r8), pointer :: livecrootn_storage(:) ! (gN/m2) live coarse root N storage
   real(r8), pointer :: deadcrootn_storage(:) ! (gN/m2) dead coarse root N storage
   real(r8), pointer :: t_soisno(:,:)         ! soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
   real(r8), pointer :: soilpsi(:,:)          ! soil water potential in each soil layer (MPa)
   real(r8), pointer :: leaf_long(:)          ! leaf longevity (yrs)
   real(r8), pointer :: stress_decid(:)       ! binary flag for stress-deciduous leaf habit (0 or 1)
   real(r8), pointer :: woody(:)              ! binary flag for woody lifeform (1=woody, 0=not woody)

!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: dormant_flag(:)    ! dormancy flag
   real(r8), pointer :: days_active(:)     ! number of days since last dormancy
   real(r8), pointer :: onset_flag(:)      ! onset flag
   real(r8), pointer :: onset_counter(:)   ! onset counter (seconds)
   real(r8), pointer :: onset_gddflag(:)   ! onset freeze flag
   real(r8), pointer :: onset_fdd(:)       ! onset freezing degree days counter
   real(r8), pointer :: onset_gdd(:)       ! onset growing degree days
   real(r8), pointer :: onset_swi(:)       ! onset soil water index
   real(r8), pointer :: offset_flag(:)     ! offset flag
   real(r8), pointer :: offset_counter(:)  ! offset counter (seconds)
   real(r8), pointer :: dayl(:)            ! daylength (seconds)
   real(r8), pointer :: offset_fdd(:)      ! offset freezing degree days counter
   real(r8), pointer :: offset_swi(:)      ! offset soil water index
   real(r8), pointer :: annavg_t2m(:)      ! annual average 2m air temperature (K)
   real(r8), pointer :: lgsf(:)            ! long growing season factor [0-1]
   real(r8), pointer :: bglfr(:)           ! background litterfall rate (1/s)
   real(r8), pointer :: bgtr(:)            ! background transfer growth rate (1/s)
   real(r8), pointer :: prev_leafc_to_litter(:)  ! previous timestep leaf C litterfall flux (gC/m2/s)
   real(r8), pointer :: prev_frootc_to_litter(:) ! previous timestep froot C litterfall flux (gC/m2/s)
   real(r8), pointer :: leafc_xfer_to_leafc(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: leafn_xfer_to_leafn(:)
   real(r8), pointer :: frootn_xfer_to_frootn(:)
   real(r8), pointer :: livestemn_xfer_to_livestemn(:)
   real(r8), pointer :: deadstemn_xfer_to_deadstemn(:)
   real(r8), pointer :: livecrootn_xfer_to_livecrootn(:)
   real(r8), pointer :: deadcrootn_xfer_to_deadcrootn(:)
   real(r8), pointer :: leafc_xfer(:)      ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc_xfer(:)     ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc_xfer(:)  ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc_xfer(:)  ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc_xfer(:) ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc_xfer(:) ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: leafn_xfer(:)      ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn_xfer(:)     ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn_xfer(:)  ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn_xfer(:)  ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn_xfer(:) ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn_xfer(:) ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: leafc_storage_to_xfer(:)
   real(r8), pointer :: frootc_storage_to_xfer(:)
   real(r8), pointer :: livestemc_storage_to_xfer(:)
   real(r8), pointer :: deadstemc_storage_to_xfer(:)
   real(r8), pointer :: livecrootc_storage_to_xfer(:)
   real(r8), pointer :: deadcrootc_storage_to_xfer(:)
   real(r8), pointer :: gresp_storage_to_xfer(:)
   real(r8), pointer :: leafn_storage_to_xfer(:)
   real(r8), pointer :: frootn_storage_to_xfer(:)
   real(r8), pointer :: livestemn_storage_to_xfer(:)
   real(r8), pointer :: deadstemn_storage_to_xfer(:)
   real(r8), pointer :: livecrootn_storage_to_xfer(:)
   real(r8), pointer :: deadcrootn_storage_to_xfer(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   real(r8),parameter :: secspqtrday = secspday / 4  ! seconds per quarter day
   integer :: c,p             ! indices
   integer :: fp              ! lake filter pft index
   real(r8):: dayspyr         ! days per year
   real(r8):: crit_onset_gdd  ! degree days for onset trigger
   real(r8):: soilt           ! temperature of top soil layer
   real(r8):: psi             ! water stress of top soil layer
   real(r8):: lat             !latitude (radians)
   real(r8):: temp            !temporary variable for daylength calculation
!EOP
!-----------------------------------------------------------------------
   ! Assign local pointers to derived type arrays (in)
    ivt                            => pft%itype
    pcolumn                        => pft%column
    pgridcell                      => pft%gridcell
    latdeg                         => grc%latdeg
    decl                           => cps%decl
    leafc_storage                  => pcs%leafc_storage
    frootc_storage                 => pcs%frootc_storage
    livestemc_storage              => pcs%livestemc_storage
    deadstemc_storage              => pcs%deadstemc_storage
    livecrootc_storage             => pcs%livecrootc_storage
    deadcrootc_storage             => pcs%deadcrootc_storage
    gresp_storage                  => pcs%gresp_storage
    leafn_storage                  => pns%leafn_storage
    frootn_storage                 => pns%frootn_storage
    livestemn_storage              => pns%livestemn_storage
    deadstemn_storage              => pns%deadstemn_storage
    livecrootn_storage             => pns%livecrootn_storage
    deadcrootn_storage             => pns%deadcrootn_storage
    soilpsi                        => cps%soilpsi
    t_soisno                       => ces%t_soisno
    leaf_long                      => pftcon%leaf_long
    woody                          => pftcon%woody
    stress_decid                   => pftcon%stress_decid

   ! Assign local pointers to derived type arrays (out)
    dormant_flag                   => pepv%dormant_flag
    days_active                    => pepv%days_active
    onset_flag                     => pepv%onset_flag
    onset_counter                  => pepv%onset_counter
    onset_gddflag                  => pepv%onset_gddflag
    onset_fdd                      => pepv%onset_fdd
    onset_gdd                      => pepv%onset_gdd
    onset_swi                      => pepv%onset_swi
    offset_flag                    => pepv%offset_flag
    offset_counter                 => pepv%offset_counter
    dayl                           => pepv%dayl
    offset_fdd                     => pepv%offset_fdd
    offset_swi                     => pepv%offset_swi
    annavg_t2m                     => pepv%annavg_t2m
    prev_leafc_to_litter           => pepv%prev_leafc_to_litter
    prev_frootc_to_litter          => pepv%prev_frootc_to_litter
    lgsf                           => pepv%lgsf
    bglfr                          => pepv%bglfr
    bgtr                           => pepv%bgtr
    leafc_xfer_to_leafc            => pcf%leafc_xfer_to_leafc
    frootc_xfer_to_frootc          => pcf%frootc_xfer_to_frootc
    livestemc_xfer_to_livestemc    => pcf%livestemc_xfer_to_livestemc
    deadstemc_xfer_to_deadstemc    => pcf%deadstemc_xfer_to_deadstemc
    livecrootc_xfer_to_livecrootc  => pcf%livecrootc_xfer_to_livecrootc
    deadcrootc_xfer_to_deadcrootc  => pcf%deadcrootc_xfer_to_deadcrootc
    leafn_xfer_to_leafn            => pnf%leafn_xfer_to_leafn
    frootn_xfer_to_frootn          => pnf%frootn_xfer_to_frootn
    livestemn_xfer_to_livestemn    => pnf%livestemn_xfer_to_livestemn
    deadstemn_xfer_to_deadstemn    => pnf%deadstemn_xfer_to_deadstemn
    livecrootn_xfer_to_livecrootn  => pnf%livecrootn_xfer_to_livecrootn
    deadcrootn_xfer_to_deadcrootn  => pnf%deadcrootn_xfer_to_deadcrootn
    leafc_xfer                     => pcs%leafc_xfer
    frootc_xfer                    => pcs%frootc_xfer
    livestemc_xfer                 => pcs%livestemc_xfer
    deadstemc_xfer                 => pcs%deadstemc_xfer
    livecrootc_xfer                => pcs%livecrootc_xfer
    deadcrootc_xfer                => pcs%deadcrootc_xfer
    leafn_xfer                     => pns%leafn_xfer
    frootn_xfer                    => pns%frootn_xfer
    livestemn_xfer                 => pns%livestemn_xfer
    deadstemn_xfer                 => pns%deadstemn_xfer
    livecrootn_xfer                => pns%livecrootn_xfer
    deadcrootn_xfer                => pns%deadcrootn_xfer
    leafc_storage_to_xfer          => pcf%leafc_storage_to_xfer
    frootc_storage_to_xfer         => pcf%frootc_storage_to_xfer
    livestemc_storage_to_xfer      => pcf%livestemc_storage_to_xfer
    deadstemc_storage_to_xfer      => pcf%deadstemc_storage_to_xfer
    livecrootc_storage_to_xfer     => pcf%livecrootc_storage_to_xfer
    deadcrootc_storage_to_xfer     => pcf%deadcrootc_storage_to_xfer
    gresp_storage_to_xfer          => pcf%gresp_storage_to_xfer
    leafn_storage_to_xfer          => pnf%leafn_storage_to_xfer
    frootn_storage_to_xfer         => pnf%frootn_storage_to_xfer
    livestemn_storage_to_xfer      => pnf%livestemn_storage_to_xfer
    deadstemn_storage_to_xfer      => pnf%deadstemn_storage_to_xfer
    livecrootn_storage_to_xfer     => pnf%livecrootn_storage_to_xfer
    deadcrootn_storage_to_xfer     => pnf%deadcrootn_storage_to_xfer

   ! set time steps
   dayspyr = get_days_per_year()

   do fp = 1,num_soilp
      p = filter_soilp(fp)
      c = pcolumn(p)

      if (stress_decid(ivt(p)) == 1._r8) then
         soilt = t_soisno(c,3)
         psi = soilpsi(c,3)

         ! use solar declination information stored during Surface Albedo()
         ! and latitude from gps to calcluate daylength (convert latitude from degrees to radians)
         ! the constant 13750.9871 is the number of seconds per radian of hour-angle

         lat = (SHR_CONST_PI/180._r8)*latdeg(pgridcell(p))
         temp = -(sin(lat)*sin(decl(c)))/(cos(lat) * cos(decl(c)))
         temp = min(1._r8,max(-1._r8,temp))
         dayl(p) = 2.0_r8 * 13750.9871_r8 * acos(temp)

         ! onset gdd sum from Biome-BGC, v4.1.2
         crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_t2m(p) - SHR_CONST_TKFRZ))


         ! update offset_counter and test for the end of the offset period
         if (offset_flag(p) == 1._r8) then
            ! decrement counter for offset period
            offset_counter(p) = offset_counter(p) - dt

            ! if this is the end of the offset_period, reset phenology
            ! flags and indices
            if (offset_counter(p) == 0._r8) then
               ! this code block was originally handled by call cn_offset_cleanup(p)
               ! inlined during vectorization
               offset_flag(p) = 0._r8
               offset_counter(p) = 0._r8
               dormant_flag(p) = 1._r8
               days_active(p) = 0._r8

               ! reset the previous timestep litterfall flux memory
               prev_leafc_to_litter(p) = 0._r8
               prev_frootc_to_litter(p) = 0._r8
            end if
         end if

         ! update onset_counter and test for the end of the onset period
         if (onset_flag(p) == 1.0_r8) then
            ! decrement counter for onset period
            onset_counter(p) = onset_counter(p) - dt

            ! if this is the end of the onset period, reset phenology
            ! flags and indices
            if (onset_counter(p) == 0.0_r8) then
               ! this code block was originally handled by call cn_onset_cleanup(p)
               ! inlined during vectorization
               onset_flag(p) = 0._r8
               onset_counter(p) = 0._r8
               ! set all transfer growth rates to 0.0
               leafc_xfer_to_leafc(p)   = 0._r8
               frootc_xfer_to_frootc(p) = 0._r8
               leafn_xfer_to_leafn(p)   = 0._r8
               frootn_xfer_to_frootn(p) = 0._r8
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_xfer_to_livestemc(p)   = 0._r8
                  deadstemc_xfer_to_deadstemc(p)   = 0._r8
                  livecrootc_xfer_to_livecrootc(p) = 0._r8
                  deadcrootc_xfer_to_deadcrootc(p) = 0._r8
                  livestemn_xfer_to_livestemn(p)   = 0._r8
                  deadstemn_xfer_to_deadstemn(p)   = 0._r8
                  livecrootn_xfer_to_livecrootn(p) = 0._r8
                  deadcrootn_xfer_to_deadcrootn(p) = 0._r8
               end if
               ! set transfer pools to 0.0
               leafc_xfer(p) = 0._r8
               leafn_xfer(p) = 0._r8
               frootc_xfer(p) = 0._r8
               frootn_xfer(p) = 0._r8
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_xfer(p) = 0._r8
                  livestemn_xfer(p) = 0._r8
                  deadstemc_xfer(p) = 0._r8
                  deadstemn_xfer(p) = 0._r8
                  livecrootc_xfer(p) = 0._r8
                  livecrootn_xfer(p) = 0._r8
                  deadcrootc_xfer(p) = 0._r8
                  deadcrootn_xfer(p) = 0._r8
               end if
            end if
         end if

         ! test for switching from dormant period to growth period
         if (dormant_flag(p) == 1._r8) then

            ! keep track of the number of freezing degree days in this
            ! dormancy period (only if the freeze flag has not previously been set
            ! for this dormancy period

            if (onset_gddflag(p) == 0._r8 .and. soilt < SHR_CONST_TKFRZ) onset_fdd(p) = onset_fdd(p) + fracday

            ! if the number of freezing degree days exceeds a critical value,
            ! then onset will require both wet soils and a critical soil
            ! temperature sum.  If this case is triggered, reset any previously
            ! accumulated value in onset_swi, so that onset now depends on
            ! the accumulated soil water index following the freeze trigger

            if (onset_fdd(p) > crit_onset_fdd) then
                onset_gddflag(p) = 1._r8
                onset_fdd(p) = 0._r8
                onset_swi(p) = 0._r8
            end if

            ! if the freeze flag is set, and if the soil is above freezing
            ! then accumulate growing degree days for onset trigger

            if (onset_gddflag(p) == 1._r8 .and. soilt > SHR_CONST_TKFRZ) then
               onset_gdd(p) = onset_gdd(p) + (soilt-SHR_CONST_TKFRZ)*fracday
            end if

            ! if soils are wet, accumulate soil water index for onset trigger
            if (psi >= soilpsi_on) onset_swi(p) = onset_swi(p) + fracday

            ! if critical soil water index is exceeded, set onset_flag, and
            ! then test for soil temperature criteria

            if (onset_swi(p) > crit_onset_swi) then
                onset_flag(p) = 1._r8

                ! only check soil temperature criteria if freeze flag set since
                ! beginning of last dormancy.  If freeze flag set and growing
                ! degree day sum (since freeze trigger) is lower than critical
                ! value, then override the onset_flag set from soil water.

                if (onset_gddflag(p) == 1._r8 .and. onset_gdd(p) < crit_onset_gdd) onset_flag(p) = 0._r8
            end if
            
            ! only allow onset if dayl > 6hrs
            if (onset_flag(p) == 1._r8 .and. dayl(p) <= secspqtrday) then
                onset_flag(p) = 0._r8
            end if

            ! if this is the beginning of the onset period
            ! then reset the phenology flags and indices

            if (onset_flag(p) == 1._r8) then
               dormant_flag(p) = 0._r8
               days_active(p) = 0._r8
               onset_gddflag(p) = 0._r8
               onset_fdd(p) = 0._r8
               onset_gdd(p) = 0._r8
               onset_swi(p) = 0._r8
               onset_counter(p) = ndays_on * secspday

               ! call subroutine to move all the storage pools into transfer pools,
               ! where they will be transfered to displayed growth over the onset period.
               ! this code was originally handled with call cn_storage_to_xfer(p)
               ! inlined during vectorization

               ! set carbon fluxes for shifting storage pools to transfer pools
               leafc_storage_to_xfer(p)  = fstor2tran * leafc_storage(p)/dt
               frootc_storage_to_xfer(p) = fstor2tran * frootc_storage(p)/dt
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemc_storage_to_xfer(p)  = fstor2tran * livestemc_storage(p)/dt
                  deadstemc_storage_to_xfer(p)  = fstor2tran * deadstemc_storage(p)/dt
                  livecrootc_storage_to_xfer(p) = fstor2tran * livecrootc_storage(p)/dt
                  deadcrootc_storage_to_xfer(p) = fstor2tran * deadcrootc_storage(p)/dt
                  gresp_storage_to_xfer(p)      = fstor2tran * gresp_storage(p)/dt
               end if

               ! set nitrogen fluxes for shifting storage pools to transfer pools
               leafn_storage_to_xfer(p)  = fstor2tran * leafn_storage(p)/dt
               frootn_storage_to_xfer(p) = fstor2tran * frootn_storage(p)/dt
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemn_storage_to_xfer(p)  = fstor2tran * livestemn_storage(p)/dt
                  deadstemn_storage_to_xfer(p)  = fstor2tran * deadstemn_storage(p)/dt
                  livecrootn_storage_to_xfer(p) = fstor2tran * livecrootn_storage(p)/dt
                  deadcrootn_storage_to_xfer(p) = fstor2tran * deadcrootn_storage(p)/dt
               end if
            end if

         ! test for switching from growth period to offset period
         else if (offset_flag(p) == 0._r8) then

            ! if soil water potential lower than critical value, accumulate
            ! as stress in offset soil water index

            if (psi <= soilpsi_off) then
               offset_swi(p) = offset_swi(p) + fracday

               ! if the offset soil water index exceeds critical value, and
               ! if this is not the middle of a previously initiated onset period,
               ! then set flag to start the offset period and reset index variables

               if (offset_swi(p) >= crit_offset_swi .and. onset_flag(p) == 0._r8) offset_flag(p) = 1._r8

            ! if soil water potential higher than critical value, reduce the
            ! offset water stress index.  By this mechanism, there must be a
            ! sustained period of water stress to initiate offset.

            else if (psi >= soilpsi_on) then
               offset_swi(p) = offset_swi(p) - fracday
               offset_swi(p) = max(offset_swi(p),0._r8)
            end if

            ! decrease freezing day accumulator for warm soil
            if (offset_fdd(p) > 0._r8 .and. soilt > SHR_CONST_TKFRZ) then
                offset_fdd(p) = offset_fdd(p) - fracday
                offset_fdd(p) = max(0._r8, offset_fdd(p))
            end if

            ! increase freezing day accumulator for cold soil
            if (soilt <= SHR_CONST_TKFRZ) then
               offset_fdd(p) = offset_fdd(p) + fracday

               ! if freezing degree day sum is greater than critical value, initiate offset
               if (offset_fdd(p) > crit_offset_fdd .and. onset_flag(p) == 0._r8) offset_flag(p) = 1._r8
            end if
            
            ! force offset if daylength is < 6 hrs
            if (dayl(p) <= secspqtrday) then
               offset_flag(p) = 1._r8
            end if

            ! if this is the beginning of the offset period
            ! then reset flags and indices
            if (offset_flag(p) == 1._r8) then
               offset_fdd(p) = 0._r8
               offset_swi(p) = 0._r8
               offset_counter(p) = ndays_off * secspday
               prev_leafc_to_litter(p) = 0._r8
               prev_frootc_to_litter(p) = 0._r8
            end if
         end if

         ! keep track of number of days since last dormancy for control on
         ! fraction of new growth to send to storage for next growing season

         if (dormant_flag(p) == 0.0_r8) then
             days_active(p) = days_active(p) + fracday
         end if

         ! calculate long growing season factor (lgsf)
         ! only begin to calculate a lgsf greater than 0.0 once the number
         ! of days active exceeds days/year.
         lgsf(p) = max(min((days_active(p)-dayspyr)/dayspyr, 1._r8),0._r8)

         ! set background litterfall rate, when not in the phenological offset period
         if (offset_flag(p) == 1._r8) then
            bglfr(p) = 0._r8
         else
            ! calculate the background litterfall rate (bglfr)
            ! in units 1/s, based on leaf longevity (yrs) and correction for long growing season

            bglfr(p) = (1._r8/(leaf_long(ivt(p))*dayspyr*secspday))*lgsf(p)
         end if

         ! set background transfer rate when active but not in the phenological onset period
         if (onset_flag(p) == 1._r8) then
            bgtr(p) = 0._r8
         else
            ! the background transfer rate is calculated as the rate that would result
            ! in complete turnover of the storage pools in one year at steady state,
            ! once lgsf has reached 1.0 (after 730 days active).

            bgtr(p) = (1._r8/(dayspyr*secspday))*lgsf(p)

            ! set carbon fluxes for shifting storage pools to transfer pools

            leafc_storage_to_xfer(p)  = leafc_storage(p) * bgtr(p)
            frootc_storage_to_xfer(p) = frootc_storage(p) * bgtr(p)
            if (woody(ivt(p)) == 1.0_r8) then
               livestemc_storage_to_xfer(p)  = livestemc_storage(p) * bgtr(p)
               deadstemc_storage_to_xfer(p)  = deadstemc_storage(p) * bgtr(p)
               livecrootc_storage_to_xfer(p) = livecrootc_storage(p) * bgtr(p)
               deadcrootc_storage_to_xfer(p) = deadcrootc_storage(p) * bgtr(p)
               gresp_storage_to_xfer(p)      = gresp_storage(p) * bgtr(p)
            end if

            ! set nitrogen fluxes for shifting storage pools to transfer pools
            leafn_storage_to_xfer(p)  = leafn_storage(p) * bgtr(p)
            frootn_storage_to_xfer(p) = frootn_storage(p) * bgtr(p)
            if (woody(ivt(p)) == 1.0_r8) then
               livestemn_storage_to_xfer(p)  = livestemn_storage(p) * bgtr(p)
               deadstemn_storage_to_xfer(p)  = deadstemn_storage(p) * bgtr(p)
               livecrootn_storage_to_xfer(p) = livecrootn_storage(p) * bgtr(p)
               deadcrootn_storage_to_xfer(p) = deadcrootn_storage(p) * bgtr(p)
            end if
         end if

      end if ! end if stress deciduous

   end do ! end of pft loop

end subroutine CNStressDecidPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropPhenology
!
! !INTERFACE:
subroutine CropPhenology(num_pcropp, filter_pcropp)

! !DESCRIPTION:
! Code from AgroIBIS to determine crop phenology and code from CN to
! handle CN fluxes during the phenological onset & offset periods.

! !USES:
  use clm_time_manager, only : get_curr_date, get_curr_calday, get_days_per_year
  use pftvarcon       , only : ncorn, nscereal, nwcereal, nsoybean, gddmin, hybgdd, &
                               lfemerg, grnfill, mxmat, minplanttemp, planttemp
  use clm_varcon      , only : spval, secspday

! !ARGUMENTS:
  integer, intent(in) :: num_pcropp       ! number of prog crop pfts in filter
  integer, intent(in) :: filter_pcropp(:) ! filter for prognostic crop pfts

! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 2/5/08:  slevis created according to AgroIBIS subroutines of Kucharik et al.
! 7/14/08: slevis adapted crop cycles to southern hemisphere
! 3/29/11: ekluzek simply logic using pftvarcon arrays

!EOP

! LOCAL VARAIBLES:
      integer kyr       ! current year
      integer kmo       !         month of year  (1, ..., 12)
      integer kda       !         day of month   (1, ..., 31)
      integer mcsec     !         seconds of day (0, ..., seconds/day)
      integer jday      ! julian day of the year
      integer fp,p      ! pft indices
      integer c         ! column indices
      integer g         ! gridcell indices
      integer h         ! hemisphere indices
      integer idpp      ! number of days past planting
      integer pmmin     ! earliest month to plant winter temperate cereal
      integer pdmin     ! earliest day in earliest month to plant
      integer pmmax     ! latest possible month (month) and
      integer pdmax     ! latest day in latest month to plant
      real(r8) dayspyr  ! days per year
      real(r8) crmcorn  ! comparitive relative maturity for corn

! local pointers to implicit in scalars

      integer , pointer :: pgridcell(:)! pft's gridcell index
      integer , pointer :: pcolumn(:)  ! pft's column index
      integer , pointer :: ivt(:)      ! pft
      real(r8), pointer :: hui(:)      ! =gdd since planting (gddplant)
      real(r8), pointer :: leafout(:)  ! =gdd from top soil layer temperature
      real(r8), pointer :: tlai(:)     ! one-sided leaf area index, no burying by snow
      real(r8), pointer :: gdd020(:)   ! 20 yr mean of gdd0
      real(r8), pointer :: gdd820(:)   ! 20 yr mean of gdd8
      real(r8), pointer :: gdd1020(:)  ! 20 yr mean of gdd10
      real(r8), pointer :: a5tmin(:)   ! 5-day running mean of min 2-m temperature
      real(r8), pointer :: a10tmin(:)  ! 10-day running mean of min 2-m temperature
      real(r8), pointer :: t10(:)      ! 10-day running mean of the 2 m temperature (K)
      real(r8), pointer :: t_ref2m_min(:)    !daily minimum of average 2 m height surface air temperature (K)
      real(r8), pointer :: bgtr(:)           ! background transfer growth rate (1/s)
      real(r8), pointer :: lgsf(:)           ! long growing season factor [0-1]
      real(r8), pointer :: offset_flag(:)    ! offset flag
      real(r8), pointer :: offset_counter(:) ! offset counter
      real(r8), pointer :: leaf_long(:)      ! leaf longevity (yrs)
      real(r8), pointer :: leafcn(:)         ! leaf C:N (gC/gN)
! local pointers to implicit out scalars
      integer , pointer :: idop(:)           ! date of planting
      integer , pointer :: harvdate(:)       ! harvest date
      logical , pointer :: croplive(:)       ! Flag, true if planted, not harvested
      logical , pointer :: cropplant(:)      ! Flag, true if crop may be planted
      real(r8), pointer :: cumvd(:)          ! cumulative vernalization d?ependence?
      real(r8), pointer :: hdidx(:)          ! cold hardening index?
      real(r8), pointer :: vf(:)             ! vernalization factor
      real(r8), pointer :: gddmaturity(:)    ! gdd needed to harvest
      real(r8), pointer :: bglfr(:)          ! background litterfall rate (1/s)
      real(r8), pointer :: huileaf(:)        ! heat unit index needed from planting to leaf emergence
      real(r8), pointer :: huigrain(:)       ! same to reach vegetative maturity
      real(r8), pointer :: onset_flag(:)     ! onset flag
      real(r8), pointer :: onset_counter(:)  ! onset counter
      real(r8), pointer :: leafc_xfer(:)     ! (gC/m2) leaf C transfer
      real(r8), pointer :: leafn_xfer(:)     ! (gN/m2) leaf N transfer
      real(r8), pointer :: dwt_seedc_to_leaf(:) ! (gC/m2/s) seed source to PFT-level
      real(r8), pointer :: dwt_seedn_to_leaf(:) ! (gN/m2/s) seed source to PFT-level
!------------------------------------------------------------------------

      pgridcell      => pft%gridcell
      pcolumn        => pft%column
      ivt            => pft%itype
      idop           => pps%idop
      harvdate       => pps%harvdate
      croplive       => pps%croplive
      cropplant      => pps%cropplant
      gddmaturity    => pps%gddmaturity
      huileaf        => pps%huileaf
      huigrain       => pps%huigrain
      hui            => pps%gddplant
      leafout        => pps%gddtsoi
      tlai           => pps%tlai
      gdd020         => pps%gdd020
      gdd820         => pps%gdd820
      gdd1020        => pps%gdd1020
      a5tmin         => pes%a5tmin
      a10tmin        => pes%a10tmin
      t10            => pes%t10
      cumvd          => pps%cumvd
      hdidx          => pps%hdidx
      vf             => pps%vf
      t_ref2m_min    => pes%t_ref2m_min
      bglfr          => pepv%bglfr
      bgtr           => pepv%bgtr
      lgsf           => pepv%lgsf
      onset_flag     => pepv%onset_flag
      offset_flag    => pepv%offset_flag
      onset_counter  => pepv%onset_counter
      offset_counter => pepv%offset_counter
      leafc_xfer     => pcs%leafc_xfer
      leafn_xfer     => pns%leafn_xfer
      leaf_long      => pftcon%leaf_long
      leafcn         => pftcon%leafcn
      dwt_seedc_to_leaf => ccf%dwt_seedc_to_leaf
      dwt_seedn_to_leaf => cnf%dwt_seedn_to_leaf
! ---------------------------------------

      ! get time info
      dayspyr = get_days_per_year()
      jday    = get_curr_calday()
      call get_curr_date(kyr, kmo, kda, mcsec)

      do fp = 1, num_pcropp
         p = filter_pcropp(fp)
         c = pcolumn(p)
         g = pgridcell(p)
         h = inhemi(p)

         ! background litterfall and transfer rates; long growing season factor

         bglfr(p) = 0._r8 ! this value changes later in a crop's life cycle
         bgtr(p)  = 0._r8
         lgsf(p)  = 0._r8

         ! ---------------------------------
         ! from AgroIBIS subroutine planting
         ! ---------------------------------

         ! in order to allow a crop to be planted only once each year
         ! initialize cropplant = .false., but hold it = .true. through the end of the year

         ! initialize other variables that are calculated for crops
         ! on an annual basis in cropresidue subroutine

         if ( jday == jdayyrstart(h) .and. mcsec == 0 )then

            ! make sure variables aren't changed at beginning of the year
            ! for a crop that is currently planted (e.g. winter temperate cereal)

            if (.not. croplive(p))  then
               cropplant(p) = .false.
               idop(p)      = NOT_Planted

               ! keep next for continuous, annual winter temperate cereal type crop;
               ! if we removed elseif,
               ! winter cereal grown continuously would amount to a cereal/fallow
               ! rotation because cereal would only be planted every other year

            else if (croplive(p) .and. ivt(p) == nwcereal) then
               cropplant(p) = .false.
!           else ! not possible to have croplive and ivt==cornORsoy? (slevis)
            end if

         end if

         if ( (.not. croplive(p)) .and. (.not. cropplant(p)) ) then

            ! gdd needed for * chosen crop and a likely hybrid (for that region) *
            ! to reach full physiological maturity

            ! based on accumulated seasonal average growing degree days from
            ! April 1 - Sept 30 (inclusive)
            ! for corn and soybeans in the United States -
            ! decided upon by what the typical average growing season length is
            ! and the gdd needed to reach maturity in those regions

            ! first choice is used for spring temperate cereal and/or soybeans and maize

            ! slevis: ibis reads xinpdate in io.f from control.crops.nc variable name 'plantdate'
            !         According to Chris Kucharik, the dataset of
            !         xinpdate was generated from a previous model run at 0.5 deg resolution

            ! winter temperate cereal : use gdd0 as a limit to plant winter cereal

            if (ivt(p) == nwcereal) then

               ! add check to only plant winter cereal after other crops (soybean, maize)
               ! have been harvested

               ! *** remember order of planting is crucial - in terms of which crops you want
               ! to be grown in what order ***

               ! in this case, corn or soybeans are assumed to be planted before
               ! cereal would be in any particular year that both pfts are allowed
               ! to grow in the same grid cell (e.g., double-cropping)

               ! slevis: harvdate below needs cropplant(p) above to be cropplant(p,ivt(p))
               !         where ivt(p) has rotated to winter cereal because
               !         cropplant through the end of the year for a harvested crop.
               !         Also harvdate(p) should be harvdate(p,ivt(p)) and should be
               !         updated on Jan 1st instead of at harvest (slevis)
               if (a5tmin(p)             /= spval                  .and. &
                   a5tmin(p)             <= minplanttemp(ivt(p))   .and. &
                   jday                  >= minplantjday(ivt(p),h) .and. &
                   (gdd020(p)            /= spval                  .and. &
                   gdd020(p)             >= gddmin(ivt(p)))) then

                  cumvd(p)       = 0._r8
                  hdidx(p)       = 0._r8
                  vf(p)          = 0._r8
                  croplive(p)    = .true.
                  cropplant(p)   = .true.
                  idop(p)        = jday
                  harvdate(p)    = NOT_Harvested
                  gddmaturity(p) = hybgdd(ivt(p))
                  leafc_xfer(p)  = 1._r8 ! initial seed at planting to appear
                  leafn_xfer(p)  = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  dwt_seedc_to_leaf(c) = dwt_seedc_to_leaf(c) + leafc_xfer(p)/dt
                  dwt_seedn_to_leaf(c) = dwt_seedn_to_leaf(c) + leafn_xfer(p)/dt

                  ! latest possible date to plant winter cereal and after all other 
                  ! crops were harvested for that year

               else if (jday       >=  maxplantjday(ivt(p),h) .and. &
                        gdd020(p)  /= spval                   .and. &
                        gdd020(p)  >= gddmin(ivt(p))) then

                  cumvd(p)       = 0._r8
                  hdidx(p)       = 0._r8
                  vf(p)          = 0._r8
                  croplive(p)    = .true.
                  cropplant(p)   = .true.
                  idop(p)        = jday
                  harvdate(p)    = NOT_Harvested
                  gddmaturity(p) = hybgdd(ivt(p))
                  leafc_xfer(p)  = 1._r8 ! initial seed at planting to appear
                  leafn_xfer(p)  = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  dwt_seedc_to_leaf(c) = dwt_seedc_to_leaf(c) + leafc_xfer(p)/dt
                  dwt_seedn_to_leaf(c) = dwt_seedn_to_leaf(c) + leafn_xfer(p)/dt
               else
                  gddmaturity(p) = 0._r8
               end if

            else ! not winter cereal... slevis: added distinction between NH and SH
               ! slevis: The idea is that jday will equal idop sooner or later in the year
               !         while the gdd part is either true or false for the year.
               if (t10(p) /= spval.and. a10tmin(p) /= spval   .and. &
                   t10(p)     > planttemp(ivt(p))             .and. &
                   a10tmin(p) > minplanttemp(ivt(p))          .and. &
                   jday       >= minplantjday(ivt(p),h)       .and. &
                   jday       <= maxplantjday(ivt(p),h)       .and. &
                   t10(p) /= spval .and. a10tmin(p) /= spval  .and. &
                   gdd820(p) /= spval                         .and. &
                   gdd820(p) >= gddmin(ivt(p))) then

                  ! impose limit on growing season length needed
                  ! for crop maturity - for cold weather constraints
                  croplive(p)  = .true.
                  cropplant(p) = .true.
                  idop(p)      = jday
                  harvdate(p)  = NOT_Harvested

                  ! go a specified amount of time before/after
                  ! climatological date
                  if (ivt(p)==nsoybean) gddmaturity(p)=min(gdd1020(p),hybgdd(ivt(p)))
                  if (ivt(p)==ncorn) then
                     gddmaturity(p)=max(950._r8, min(gdd820(p)*0.85_r8, hybgdd(ivt(p))))
                     gddmaturity(p)=max(950._r8, min(gddmaturity(p)+150._r8,1850._r8))
                  end if
                  if (ivt(p)==nscereal) gddmaturity(p)=min(gdd020(p),hybgdd(ivt(p)))

                  leafc_xfer(p) = 1._r8 ! initial seed at planting to appear
                  leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  dwt_seedc_to_leaf(c) = dwt_seedc_to_leaf(c) + leafc_xfer(p)/dt
                  dwt_seedn_to_leaf(c) = dwt_seedn_to_leaf(c) + leafn_xfer(p)/dt

               ! If hit the max planting julian day -- go ahead and plant
               else if (jday == maxplantjday(ivt(p),h) .and. gdd820(p) > 0._r8 .and. &
                        gdd820(p) /= spval ) then
                  croplive(p)  = .true.
                  cropplant(p) = .true.
                  idop(p)      = jday
                  harvdate(p)  = NOT_Harvested

                  if (ivt(p)==nsoybean) gddmaturity(p)=min(gdd1020(p),hybgdd(ivt(p)))
                  if (ivt(p)==ncorn) gddmaturity(p)=max(950._r8, min(gdd820(p)*0.85_r8, hybgdd(ivt(p))))
                  if (ivt(p)==nscereal) gddmaturity(p)=min(gdd020(p),hybgdd(ivt(p)))

                  leafc_xfer(p) = 1._r8 ! initial seed at planting to appear
                  leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  dwt_seedc_to_leaf(c) = dwt_seedc_to_leaf(c) + leafc_xfer(p)/dt
                  dwt_seedn_to_leaf(c) = dwt_seedn_to_leaf(c) + leafn_xfer(p)/dt

               else
                  gddmaturity(p) = 0._r8
               end if
            end if ! crop pft distinction

            ! crop phenology (gdd thresholds) controlled by gdd needed for
            ! maturity (physiological) which is based on the average gdd
            ! accumulation and hybrids in United States from April 1 - Sept 30

            ! calculate threshold from phase 1 to phase 2:
            ! threshold for attaining leaf emergence (based on fraction of
            ! gdd(i) -- climatological average)
            ! Hayhoe and Dwyer, 1990, Can. J. Soil Sci 70:493-497
            ! Carlson and Gage, 1989, Agric. For. Met., 45: 313-324
            ! J.T. Ritchie, 1991: Modeling Plant and Soil systems

            huileaf(p) = lfemerg(ivt(p)) * gddmaturity(p) ! 3-7% in cereal

            ! calculate threshhold from phase 2 to phase 3:
            ! from leaf emergence to beginning of grain-fill period
            ! this hypothetically occurs at the end of tassling, not the beginning
            ! tassel initiation typically begins at 0.5-0.55 * gddmaturity

            ! calculate linear relationship between huigrain fraction and relative
            ! maturity rating for maize

            if (ivt(p) == ncorn) then
               ! the following estimation of crmcorn from gddmaturity is based on a linear
               ! regression using data from Pioneer-brand corn hybrids (Kucharik, 2003,
               ! Earth Interactions 7:1-33: fig. 2)
               crmcorn = max(73._r8, min(135._r8, (gddmaturity(p)+ 53.683_r8)/13.882_r8))

               ! the following adjustment of grnfill based on crmcorn is based on a tuning
               ! of Agro-IBIS to give reasonable results for max LAI and the seasonal
               ! progression of LAI growth (pers. comm. C. Kucharik June 10, 2010)
               huigrain(p) = -0.002_r8  * (crmcorn - 73._r8) + grnfill(ivt(p))

               huigrain(p) = min(max(huigrain(p), grnfill(ivt(p))-0.1_r8), grnfill(ivt(p)))
               huigrain(p) = huigrain(p) * gddmaturity(p)     ! Cabelguenne et
            else
               huigrain(p) = grnfill(ivt(p)) * gddmaturity(p) ! al. 1999
            end if

         end if ! crop not live nor planted

         ! ----------------------------------
         ! from AgroIBIS subroutine phenocrop
         ! ----------------------------------

         ! all of the phenology changes are based on the total number of gdd needed
         ! to change to the next phase - based on fractions of the total gdd typical
         ! for  that region based on the April 1 - Sept 30 window of development

         ! crop phenology (gdd thresholds) controlled by gdd needed for
         ! maturity (physiological) which is based on the average gdd
         ! accumulation and hybrids in United States from April 1 - Sept 30
         
         ! Phase 1: Planting to leaf emergence (now in CNAllocation)
         ! Phase 2: Leaf emergence to beginning of grain fill (general LAI accumulation)
         ! Phase 3: Grain fill to physiological maturity and harvest (LAI decline)
         ! Harvest: if gdd past grain fill initiation exceeds limit
         ! or number of days past planting reaches a maximum, the crop has
         ! reached physiological maturity and plant is harvested;
         ! crop could be live or dead at this stage - these limits
         ! could lead to reaching physiological maturity or determining
         ! a harvest date for a crop killed by an early frost (see next comments)
         ! --- --- ---
         ! keeping comments without the code (slevis):
         ! if minimum temperature, t_ref2m_min <= freeze kill threshold, tkill
         ! for 3 consecutive days and lai is above a minimum,
         ! plant will be damaged/killed. This function is more for spring freeze events
         ! or for early fall freeze events

         ! spring temperate cereal is affected by this, winter cereal kill function
         ! is determined in crops.f - is a more elaborate function of
         ! cold hardening of the plant

         ! currently simulates too many grid cells killed by freezing temperatures

         ! removed on March 12 2002 - C. Kucharik
         ! until it can be a bit more refined, or used at a smaller scale.
         ! we really have no way of validating this routine
         ! too difficult to implement on 0.5 degree scale grid cells
         ! --- --- ---

         onset_flag(p)  = 0._r8 ! CN terminology to trigger certain
         offset_flag(p) = 0._r8 ! carbon and nitrogen transfers

         if (croplive(p)) then

            ! call vernalization if winter temperate cereal planted, living, and the
            ! vernalization factor is not 1;
            ! vf affects the calculation of gddtsoi & gddplant

            if (t_ref2m_min(p) < 1.e30_r8 .and. vf(p) /= 1._r8 .and. ivt(p) == nwcereal) then
               call vernalization(p)
            end if

            ! days past planting may determine harvest

            if (jday >= idop(p)) then
               idpp = jday - idop(p)
            else
               idpp = int(dayspyr) + jday - idop(p)
            end if

            ! onset_counter initialized to zero when .not. croplive
            ! offset_counter relevant only at time step of harvest

            onset_counter(p) = onset_counter(p) - dt

            ! enter phase 2 onset for one time step:
            ! transfer seed carbon to leaf emergence

            if (leafout(p) >= huileaf(p) .and. hui(p) < huigrain(p) .and. idpp < mxmat(ivt(p))) then
               if (abs(onset_counter(p)) > 1.e-6_r8) then
                  onset_flag(p)    = 1._r8
                  onset_counter(p) = dt
               else
                  onset_counter(p) = dt ! ensure no re-entry to onset of phase2
               end if

            ! enter harvest for one time step:
            ! - transfer live biomass to litter and to crop yield
            ! - send xsmrpool to the atmosphere
            ! if onset and harvest needed to last longer than one timestep
            ! the onset_counter would change from dt and you'd need to make
            ! changes to the offset subroutine below

            else if (hui(p) >= gddmaturity(p) .or. idpp >= mxmat(ivt(p))) then
               if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday
               croplive(p) = .false.     ! no re-entry in greater if-block
               if (tlai(p) > 0._r8) then ! plant had emerged before harvest
                  offset_flag(p) = 1._r8
                  offset_counter(p) = dt
               else                      ! plant never emerged from the ground
                  dwt_seedc_to_leaf(c) = dwt_seedc_to_leaf(c) - leafc_xfer(p)/dt
                  dwt_seedn_to_leaf(c) = dwt_seedn_to_leaf(c) - leafn_xfer(p)/dt
                  leafc_xfer(p) = 0._r8  ! revert planting transfers
                  leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
               end if

            ! enter phase 3 while previous criteria fail and next is true;
            ! in terms of order, phase 3 occurs before harvest, but when
            ! harvest *can* occur, we want it to have first priority.
            ! AgroIBIS uses a complex formula for lai decline.
            ! Use CN's simple formula at least as a place holder (slevis)

            else if (hui(p) >= huigrain(p)) then
               bglfr(p) = 1._r8/(leaf_long(ivt(p))*dayspyr*secspday)
            end if

         else   ! crop not live
            onset_counter(p) = 0._r8
            leafc_xfer(p) = 0._r8
            leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
         end if ! croplive

      end do ! prognostic crops loop

end subroutine CropPhenology
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropPhenologyInit
!
! !INTERFACE:
subroutine CropPhenologyInit( begp, endp )

! !DESCRIPTION:
! Initialization of CropPhenology. Must be called after time-manager is
! initialized, and after pftcon file is read in.
!
! !USES:
   use pftvarcon       , only: nwcereal, nsoybean, ncorn, nscereal,    &
                               npcropmin, npcropmax, mnNHplantdate,  &
                               mnSHplantdate, mxNHplantdate,         &
                               mxSHplantdate
   use clm_time_manager, only: get_calday
!
! !ARGUMENTS:
   implicit none
   integer, intent(IN) :: begp, endp ! Beginning and ending PFT index
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP

! LOCAL VARAIBLES:
   real(r8), pointer :: latdeg(:)                   ! latitude (radians)
   integer , pointer :: pgridcell(:)                ! pft's gridcell index
   integer           :: p,g,n,i                     ! indices
!------------------------------------------------------------------------
   latdeg         => grc%latdeg
   pgridcell      => pft%gridcell

   allocate( inhemi(begp:endp) )

   ! Julian day for the start of the year (mid-winter)
   jdayyrstart(inNH) =   1
   jdayyrstart(inSH) = 182

   ! Convert planting dates into julian day
   minplantjday(:,:) = huge(1)
   maxplantjday(:,:) = huge(1)
   do n = npcropmin, npcropmax
      minplantjday(n,inNH) = int( get_calday( mnNHplantdate(n), 0 ) )
      maxplantjday(n,inNH) = int( get_calday( mxNHplantdate(n), 0 ) )
   end do
   do n = npcropmin, npcropmax
      minplantjday(n,inSH) = int( get_calday( mnSHplantdate(n), 0 ) )
      maxplantjday(n,inSH) = int( get_calday( mxSHplantdate(n), 0 ) )
   end do

   ! Figure out what hemisphere each PFT is in
   do p = begp, endp
      g = pgridcell(p)
      ! Northern hemisphere
      if ( latdeg(g) > 0.0_r8 )then
         inhemi(p) = inNH
      else
         inhemi(p) = inSH
      end if
   end do

   !
   ! Constants for Crop vernalization
   !

   ! photoperiod factor calculation
   ! genetic constant - can be modified

   p1d = 0.004_r8  ! average for genotypes from Ritchey, 1991.
                   ! Modeling plant & soil systems: Wheat phasic developmt
   p1v = 0.003_r8  ! average for genotypes from Ritchey, 1991.

   hti   = 1._r8
   tbase = 0._r8

end subroutine CropPhenologyInit

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: vernalization
!
! !INTERFACE:
  subroutine vernalization(p)
!
! !DESCRIPTION:
!
! * * * only call for winter temperate cereal * * *
!
! subroutine calculates vernalization and photoperiod effects on
! gdd accumulation in winter temperate cereal varieties. Thermal time accumulation
! is reduced in 1st period until plant is fully vernalized. During this
! time of emergence to spikelet formation, photoperiod can also have a
! drastic effect on plant development.
!
! !ARGUMENTS:
      implicit none
      integer, intent(in) :: p    ! PFT index running over
!
! !REVISION HISTORY:
! Created by Sam Levis from AGROIBIS
!
!EOP

! LOCAL VARAIBLES:
      real(r8) tcrown                     ! ?
      real(r8) vd, vd1, vd2               ! vernalization dependence
      real(r8) tkil                       ! Freeze kill threshold
      integer  c,g                        ! indices
! local pointers to implicit in scalars
      integer , pointer :: pcolumn(:)     ! pft's column index
      logical , pointer :: croplive(:)    ! Flag, true if planted, not harvested
      real(r8), pointer :: tlai(:)        ! one-sided leaf area index, no burying by snow
      real(r8), pointer :: t_ref2m(:)     ! 2 m height surface air temperature (K)
      real(r8), pointer :: t_ref2m_min(:) !daily minimum of average 2 m height surface air temperature (K)
      real(r8), pointer :: t_ref2m_max(:) !daily maximum of average 2 m height surface air temperature (K)
      real(r8), pointer :: snowdp(:)      ! snow height (m)
! local pointers to implicit out scalars
      real(r8), pointer :: vf(:)          ! vernalization factor for cereal
      real(r8), pointer :: cumvd(:)       ! cumulative vernalization d?ependence?
      real(r8), pointer :: gddmaturity(:) ! gdd needed to harvest
      real(r8), pointer :: huigrain(:)    ! heat unit index needed to reach vegetative maturity
      real(r8), pointer :: hdidx(:)       ! cold hardening index?
!------------------------------------------------------------------------

        pcolumn     => pft%column
        croplive    => pps%croplive
        hdidx       => pps%hdidx
        cumvd       => pps%cumvd
        vf          => pps%vf
        gddmaturity => pps%gddmaturity
        huigrain    => pps%huigrain
        tlai        => pps%tlai
        t_ref2m     => pes%t_ref2m
        t_ref2m_min => pes%t_ref2m_min
        t_ref2m_max => pes%t_ref2m_max
        snowdp      => cps%snowdp

        c = pcolumn(p)

        ! for all equations - temperatures must be in degrees (C)
        ! calculate temperature of crown of crop (e.g., 3 cm soil temperature)
        ! snow depth in centimeters

        if (t_ref2m(p) < tfrz) then !slevis: t_ref2m inst of td=daily avg (K)
           tcrown = 2._r8 + (t_ref2m(p) - tfrz) * (0.4_r8 + 0.0018_r8 * &
                    (min(snowdp(c)*100._r8, 15._r8) - 15._r8)**2)
        else !slevis: snowdp inst of adsnod=daily average (m)
           tcrown = t_ref2m(p) - tfrz
        end if

        ! vernalization factor calculation
        ! if vf(p) = 1.  then plant is fully vernalized - and thermal time
        ! accumulation in phase 1 will be unaffected
        ! refers to gddtsoi & gddplant, defined in the accumulation routines (slevis)
        ! reset vf, cumvd, and hdidx to 0 at planting of crop (slevis)

        if (t_ref2m_max(p) > tfrz) then
           if (t_ref2m_min(p) <= tfrz+15._r8) then
             vd1      = 1.4_r8 - 0.0778_r8 * tcrown
             vd2      = 0.5_r8 + 13.44_r8 / ((t_ref2m_max(p)-t_ref2m_min(p)+3._r8)**2) * tcrown
             vd       = max(0._r8, min(1._r8, vd1, vd2))
             cumvd(p) = cumvd(p) + vd
           end if

           if (cumvd(p) < 10._r8 .and. t_ref2m_max(p) > tfrz+30._r8) then
             cumvd(p) = cumvd(p) - 0.5_r8 * (t_ref2m_max(p) - tfrz - 30._r8)
           end if
           cumvd(p) = max(0._r8, cumvd(p))       ! must be > 0

           vf(p) = 1._r8 - p1v * (50._r8 - cumvd(p))
           vf(p) = max(0._r8, min(vf(p), 1._r8)) ! must be between 0 - 1
        end if

        ! calculate cold hardening of plant
        ! determines for winter cereal varieties whether the plant has completed
        ! a period of cold hardening to protect it from freezing temperatures. If
        ! not, then exposure could result in death or killing of plants.

        ! there are two distinct phases of hardening

        if (t_ref2m_min(p) <= tfrz-3._r8 .or. hdidx(p) /= 0._r8) then
           if (hdidx(p) >= hti) then   ! done with phase 1
              hdidx(p) = hdidx(p) + 0.083_r8
              hdidx(p) = min(hdidx(p), hti*2._r8)
           end if

           if (t_ref2m_max(p) >= tbase + tfrz + 10._r8) then
              hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
              if (hdidx(p) > hti) hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
              hdidx(p) = max(0._r8, hdidx(p))
           end if

        else if (tcrown >= tbase-1._r8) then
           if (tcrown <= tbase+8._r8) then
              hdidx(p) = hdidx(p) + 0.1_r8 - (tcrown-tbase+3.5_r8)**2 / 506._r8
              if (hdidx(p) >= hti .and. tcrown <= tbase + 0._r8) then
                 hdidx(p) = hdidx(p) + 0.083_r8
                 hdidx(p) = min(hdidx(p), hti*2._r8)
              end if
           end if

           if (t_ref2m_max(p) >= tbase + tfrz + 10._r8) then
              hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
              if (hdidx(p) > hti) hdidx(p) = hdidx(p) - 0.02_r8 * (t_ref2m_max(p)-tbase-tfrz-10._r8)
              hdidx(p) = max(0._r8, hdidx(p))
           end if
        end if

        ! calculate what the cereal killing temperature
        ! there is a linear inverse relationship between
        ! hardening of the plant and the killing temperature or
        ! threshold that the plant can withstand
        ! when plant is fully-hardened (hdidx = 2), the killing threshold is -18 C

        ! will have to develop some type of relationship that reduces LAI and
        ! biomass pools in response to cold damaged crop

        if (t_ref2m_min(p) <= tfrz - 6._r8) then
           tkil = (tbase - 6._r8) - 6._r8 * hdidx(p)
           if (tkil >= tcrown) then
              if ((0.95_r8 - 0.02_r8 * (tcrown - tkil)**2) >= 0.02_r8) then
                 write (iulog,*)  'crop damaged by cold temperatures at p,c =', p,c
              else if (tlai(p) > 0._r8) then ! slevis: kill if past phase1
                 gddmaturity(p) = 0._r8      !         by forcing through
                 huigrain(p)    = 0._r8      !         harvest
                 write (iulog,*)  '95% of crop killed by cold temperatures at p,c =', p,c
              end if
           end if
        end if

  end subroutine vernalization

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNOnsetGrowth
!
! !INTERFACE:
subroutine CNOnsetGrowth (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Determines the flux of stored C and N from transfer pools to display
! pools during the phenological onset period.
!
! !USES:
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/27/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)             ! pft vegetation type
   real(r8), pointer :: onset_flag(:)      ! onset flag
   real(r8), pointer :: onset_counter(:)   ! onset days counter
   real(r8), pointer :: leafc_xfer(:)      ! (gC/m2) leaf C transfer
   real(r8), pointer :: frootc_xfer(:)     ! (gC/m2) fine root C transfer
   real(r8), pointer :: livestemc_xfer(:)  ! (gC/m2) live stem C transfer
   real(r8), pointer :: deadstemc_xfer(:)  ! (gC/m2) dead stem C transfer
   real(r8), pointer :: livecrootc_xfer(:) ! (gC/m2) live coarse root C transfer
   real(r8), pointer :: deadcrootc_xfer(:) ! (gC/m2) dead coarse root C transfer
   real(r8), pointer :: leafn_xfer(:)      ! (gN/m2) leaf N transfer
   real(r8), pointer :: frootn_xfer(:)     ! (gN/m2) fine root N transfer
   real(r8), pointer :: livestemn_xfer(:)  ! (gN/m2) live stem N transfer
   real(r8), pointer :: deadstemn_xfer(:)  ! (gN/m2) dead stem N transfer
   real(r8), pointer :: livecrootn_xfer(:) ! (gN/m2) live coarse root N transfer
   real(r8), pointer :: deadcrootn_xfer(:) ! (gN/m2) dead coarse root N transfer
   real(r8), pointer :: woody(:)           ! binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: bgtr(:)            ! background transfer growth rate (1/s)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: leafc_xfer_to_leafc(:)
   real(r8), pointer :: frootc_xfer_to_frootc(:)
   real(r8), pointer :: livestemc_xfer_to_livestemc(:)
   real(r8), pointer :: deadstemc_xfer_to_deadstemc(:)
   real(r8), pointer :: livecrootc_xfer_to_livecrootc(:)
   real(r8), pointer :: deadcrootc_xfer_to_deadcrootc(:)
   real(r8), pointer :: leafn_xfer_to_leafn(:)
   real(r8), pointer :: frootn_xfer_to_frootn(:)
   real(r8), pointer :: livestemn_xfer_to_livestemn(:)
   real(r8), pointer :: deadstemn_xfer_to_deadstemn(:)
   real(r8), pointer :: livecrootn_xfer_to_livecrootn(:)
   real(r8), pointer :: deadcrootn_xfer_to_deadcrootn(:)
!
! local pointers to implicit out scalars
!
! !OTHER LOCAL VARIABLES:
   integer :: p            ! indices
   integer :: fp           ! lake filter pft index
   real(r8):: t1           ! temporary variable

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => pft%itype
    onset_flag                     => pepv%onset_flag
    onset_counter                  => pepv%onset_counter
    leafc_xfer                     => pcs%leafc_xfer
    frootc_xfer                    => pcs%frootc_xfer
    livestemc_xfer                 => pcs%livestemc_xfer
    deadstemc_xfer                 => pcs%deadstemc_xfer
    livecrootc_xfer                => pcs%livecrootc_xfer
    deadcrootc_xfer                => pcs%deadcrootc_xfer
    leafn_xfer                     => pns%leafn_xfer
    frootn_xfer                    => pns%frootn_xfer
    livestemn_xfer                 => pns%livestemn_xfer
    deadstemn_xfer                 => pns%deadstemn_xfer
    livecrootn_xfer                => pns%livecrootn_xfer
    deadcrootn_xfer                => pns%deadcrootn_xfer
    bgtr                           => pepv%bgtr
    woody                          => pftcon%woody

   ! assign local pointers to derived type arrays (out)
    leafc_xfer_to_leafc            => pcf%leafc_xfer_to_leafc
    frootc_xfer_to_frootc          => pcf%frootc_xfer_to_frootc
    livestemc_xfer_to_livestemc    => pcf%livestemc_xfer_to_livestemc
    deadstemc_xfer_to_deadstemc    => pcf%deadstemc_xfer_to_deadstemc
    livecrootc_xfer_to_livecrootc  => pcf%livecrootc_xfer_to_livecrootc
    deadcrootc_xfer_to_deadcrootc  => pcf%deadcrootc_xfer_to_deadcrootc
    leafn_xfer_to_leafn            => pnf%leafn_xfer_to_leafn
    frootn_xfer_to_frootn          => pnf%frootn_xfer_to_frootn
    livestemn_xfer_to_livestemn    => pnf%livestemn_xfer_to_livestemn
    deadstemn_xfer_to_deadstemn    => pnf%deadstemn_xfer_to_deadstemn
    livecrootn_xfer_to_livecrootn  => pnf%livecrootn_xfer_to_livecrootn
    deadcrootn_xfer_to_deadcrootn  => pnf%deadcrootn_xfer_to_deadcrootn

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! only calculate these fluxes during onset period
      if (onset_flag(p) == 1._r8) then

         ! The transfer rate is a linearly decreasing function of time,
         ! going to zero on the last timestep of the onset period

         if (onset_counter(p) == dt) then
             t1 = 1.0_r8 / dt
         else
             t1 = 2.0_r8 / (onset_counter(p))
         end if
         leafc_xfer_to_leafc(p)   = t1 * leafc_xfer(p)
         frootc_xfer_to_frootc(p) = t1 * frootc_xfer(p)
         leafn_xfer_to_leafn(p)   = t1 * leafn_xfer(p)
         frootn_xfer_to_frootn(p) = t1 * frootn_xfer(p)
         if (woody(ivt(p)) == 1.0_r8) then
             livestemc_xfer_to_livestemc(p)   = t1 * livestemc_xfer(p)
             deadstemc_xfer_to_deadstemc(p)   = t1 * deadstemc_xfer(p)
             livecrootc_xfer_to_livecrootc(p) = t1 * livecrootc_xfer(p)
             deadcrootc_xfer_to_deadcrootc(p) = t1 * deadcrootc_xfer(p)
             livestemn_xfer_to_livestemn(p)   = t1 * livestemn_xfer(p)
             deadstemn_xfer_to_deadstemn(p)   = t1 * deadstemn_xfer(p)
             livecrootn_xfer_to_livecrootn(p) = t1 * livecrootn_xfer(p)
             deadcrootn_xfer_to_deadcrootn(p) = t1 * deadcrootn_xfer(p)
         end if

      end if ! end if onset period

      ! calculate the background rate of transfer growth (used for stress
      ! deciduous algorithm). In this case, all of the mass in the transfer
      ! pools should be moved to displayed growth in each timestep.

      if (bgtr(p) > 0._r8) then
         leafc_xfer_to_leafc(p)   = leafc_xfer(p) / dt
         frootc_xfer_to_frootc(p) = frootc_xfer(p) / dt
         leafn_xfer_to_leafn(p)   = leafn_xfer(p) / dt
         frootn_xfer_to_frootn(p) = frootn_xfer(p) / dt
         if (woody(ivt(p)) == 1.0_r8) then
             livestemc_xfer_to_livestemc(p)   = livestemc_xfer(p) / dt
             deadstemc_xfer_to_deadstemc(p)   = deadstemc_xfer(p) / dt
             livecrootc_xfer_to_livecrootc(p) = livecrootc_xfer(p) / dt
             deadcrootc_xfer_to_deadcrootc(p) = deadcrootc_xfer(p) / dt
             livestemn_xfer_to_livestemn(p)   = livestemn_xfer(p) / dt
             deadstemn_xfer_to_deadstemn(p)   = deadstemn_xfer(p) / dt
             livecrootn_xfer_to_livecrootn(p) = livecrootn_xfer(p) / dt
             deadcrootn_xfer_to_deadcrootn(p) = deadcrootn_xfer(p) / dt
         end if
      end if ! end if bgtr

   end do ! end pft loop

end subroutine CNOnsetGrowth
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNOffsetLitterfall
!
! !INTERFACE:
subroutine CNOffsetLitterfall (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Determines the flux of C and N from displayed pools to litter
! pools during the phenological offset period.
!
! !USES:
   use pftvarcon       , only: npcropmin
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/27/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)                   ! pft vegetation type
   real(r8), pointer :: offset_flag(:)           ! offset flag
   real(r8), pointer :: offset_counter(:)        ! offset days counter
   real(r8), pointer :: leafc(:)                 ! (gC/m2) leaf C
   real(r8), pointer :: frootc(:)                ! (gC/m2) fine root C
   real(r8), pointer :: cpool_to_leafc(:)        ! allocation to leaf C (gC/m2/s)
   real(r8), pointer :: cpool_to_frootc(:)       ! allocation to fine root C (gC/m2/s)
!  integer , pointer :: pcolumn(:)               ! pft's column index
   real(r8), pointer :: grainc(:)                ! (gC/m2) grain C
   real(r8), pointer :: livestemc(:)             ! (gC/m2) livestem C
   real(r8), pointer :: cpool_to_grainc(:)       ! allocation to grain C (gC/m2/s)
   real(r8), pointer :: cpool_to_livestemc(:)    ! allocation to live stem C (gC/m2/s)
   real(r8), pointer :: livewdcn(:)              ! live wood C:N (gC/gN)
   real(r8), pointer :: graincn(:)               ! grain C:N (gC/gN)
   real(r8), pointer :: leafcn(:)                ! leaf C:N (gC/gN)
   real(r8), pointer :: lflitcn(:)               ! leaf litter C:N (gC/gN)
   real(r8), pointer :: frootcn(:)               ! fine root C:N (gC/gN)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: prev_leafc_to_litter(:)  ! previous timestep leaf C litterfall flux (gC/m2/s)
   real(r8), pointer :: prev_frootc_to_litter(:) ! previous timestep froot C litterfall flux (gC/m2/s)
   real(r8), pointer :: leafc_to_litter(:)       ! leaf C litterfall (gC/m2/s)
   real(r8), pointer :: frootc_to_litter(:)      ! fine root C litterfall (gC/m2/s)
   real(r8), pointer :: leafn_to_litter(:)       ! leaf N litterfall (gN/m2/s)
   real(r8), pointer :: leafn_to_retransn(:)     ! leaf N to retranslocated N pool (gN/m2/s)
   real(r8), pointer :: frootn_to_litter(:)      ! fine root N litterfall (gN/m2/s)
   real(r8), pointer :: livestemc_to_litter(:)   ! live stem C litterfall (gC/m2/s)
   real(r8), pointer :: grainc_to_food(:)        ! grain C to food (gC/m2/s)
   real(r8), pointer :: livestemn_to_litter(:)   ! livestem N to litter (gN/m2/s)
   real(r8), pointer :: grainn_to_food(:)        ! grain N to food (gN/m2/s)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p, c         ! indices
   integer :: fp           ! lake filter pft index
   real(r8):: t1           ! temporary variable

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => pft%itype
    offset_flag                    => pepv%offset_flag
    offset_counter                 => pepv%offset_counter
    leafc                          => pcs%leafc
    frootc                         => pcs%frootc
    grainc                         => pcs%grainc
    livestemc                      => pcs%livestemc
    cpool_to_grainc                => pcf%cpool_to_grainc
    cpool_to_livestemc             => pcf%cpool_to_livestemc
    cpool_to_leafc                 => pcf%cpool_to_leafc
    cpool_to_frootc                => pcf%cpool_to_frootc
    leafcn                         => pftcon%leafcn
    lflitcn                        => pftcon%lflitcn
    frootcn                        => pftcon%frootcn
    livewdcn                       => pftcon%livewdcn
    graincn                        => pftcon%graincn

   ! assign local pointers to derived type arrays (out)
    prev_leafc_to_litter           => pepv%prev_leafc_to_litter
    prev_frootc_to_litter          => pepv%prev_frootc_to_litter
    leafc_to_litter                => pcf%leafc_to_litter
    frootc_to_litter               => pcf%frootc_to_litter
    livestemc_to_litter            => pcf%livestemc_to_litter
    grainc_to_food                 => pcf%grainc_to_food
    livestemn_to_litter            => pnf%livestemn_to_litter
    grainn_to_food                 => pnf%grainn_to_food
    leafn_to_litter                => pnf%leafn_to_litter
    leafn_to_retransn              => pnf%leafn_to_retransn
    frootn_to_litter               => pnf%frootn_to_litter

   ! The litterfall transfer rate starts at 0.0 and increases linearly
   ! over time, with displayed growth going to 0.0 on the last day of litterfall

   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! only calculate fluxes during offset period
      if (offset_flag(p) == 1._r8) then

         if (offset_counter(p) == dt) then
             t1 = 1.0_r8 / dt
             leafc_to_litter(p)  = t1 * leafc(p)  + cpool_to_leafc(p)
             frootc_to_litter(p) = t1 * frootc(p) + cpool_to_frootc(p)
             ! this assumes that offset_counter == dt for crops
             ! if this were ever changed, we'd need to add code to the "else"
             if (ivt(p) >= npcropmin) then
                grainc_to_food(p) = t1 * grainc(p)  + cpool_to_grainc(p) 
                livestemc_to_litter(p) = t1 * livestemc(p)  + cpool_to_livestemc(p)
             end if
         else
             t1 = dt * 2.0_r8 / (offset_counter(p) * offset_counter(p))
             leafc_to_litter(p)  = prev_leafc_to_litter(p)  + t1*(leafc(p)  - prev_leafc_to_litter(p)*offset_counter(p))
             frootc_to_litter(p) = prev_frootc_to_litter(p) + t1*(frootc(p) - prev_frootc_to_litter(p)*offset_counter(p))
         end if

         ! calculate the leaf N litterfall and retranslocation
         leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
         leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

         ! calculate fine root N litterfall (no retranslocation of fine root N)
         frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))

         if (ivt(p) >= npcropmin) then
            livestemn_to_litter(p) = livestemc_to_litter(p) / livewdcn(ivt(p))
            grainn_to_food(p) = grainc_to_food(p) / graincn(ivt(p))
         end if

         ! save the current litterfall fluxes
         prev_leafc_to_litter(p)  = leafc_to_litter(p)
         prev_frootc_to_litter(p) = frootc_to_litter(p)

      end if ! end if offset period

   end do ! end pft loop

end subroutine CNOffsetLitterfall
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNBackgroundLitterfall
!
! !INTERFACE:
subroutine CNBackgroundLitterfall (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Determines the flux of C and N from displayed pools to litter
! pools as the result of background litter fall.
!
! !USES:
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 10/2/03: Created by Peter Thornton
! 10/24/03, Peter Thornton: migrated to vector data structures
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! pft level
   integer , pointer :: ivt(:)       ! pft vegetation type
   real(r8), pointer :: bglfr(:)     ! background litterfall rate (1/s)
   real(r8), pointer :: leafc(:)     ! (gC/m2) leaf C
   real(r8), pointer :: frootc(:)    ! (gC/m2) fine root C
   ! ecophysiological constants
   real(r8), pointer :: leafcn(:)    ! leaf C:N (gC/gN)
   real(r8), pointer :: lflitcn(:)   ! leaf litter C:N (gC/gN)
   real(r8), pointer :: frootcn(:)   ! fine root C:N (gC/gN)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: leafc_to_litter(:)
   real(r8), pointer :: frootc_to_litter(:)
   real(r8), pointer :: leafn_to_litter(:)
   real(r8), pointer :: leafn_to_retransn(:)
   real(r8), pointer :: frootn_to_litter(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p            ! indices
   integer :: fp           ! lake filter pft index

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => pft%itype
    bglfr                          => pepv%bglfr
    leafc                          => pcs%leafc
    frootc                         => pcs%frootc
    leafcn                         => pftcon%leafcn
    lflitcn                        => pftcon%lflitcn
    frootcn                        => pftcon%frootcn

   ! assign local pointers to derived type arrays (out)
    leafc_to_litter                => pcf%leafc_to_litter
    frootc_to_litter               => pcf%frootc_to_litter
    leafn_to_litter                => pnf%leafn_to_litter
    leafn_to_retransn              => pnf%leafn_to_retransn
    frootn_to_litter               => pnf%frootn_to_litter

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! only calculate these fluxes if the background litterfall rate is non-zero
      if (bglfr(p) > 0._r8) then
         ! units for bglfr are already 1/s
         leafc_to_litter(p)  = bglfr(p) * leafc(p)
         frootc_to_litter(p) = bglfr(p) * frootc(p)

         ! calculate the leaf N litterfall and retranslocation
         leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
         leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

         ! calculate fine root N litterfall (no retranslocation of fine root N)
         frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))

      end if

   end do

end subroutine CNBackgroundLitterfall
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNLivewoodTurnover
!
! !INTERFACE:
subroutine CNLivewoodTurnover (num_soilp, filter_soilp)
!
! !DESCRIPTION:
! Determines the flux of C and N from live wood to
! dead wood pools, for stem and coarse root.
!
! !USES:
!
! !ARGUMENTS:
   integer, intent(in) :: num_soilp       ! number of soil pfts in filter
   integer, intent(in) :: filter_soilp(:) ! filter for soil pfts
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 12/5/03: created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   ! pft level
   integer , pointer :: ivt(:)         ! pft vegetation type
   real(r8), pointer :: livestemc(:)   ! (gC/m2) live stem C
   real(r8), pointer :: livecrootc(:)  ! (gC/m2) live coarse root C
   real(r8), pointer :: livestemn(:)   ! (gN/m2) live stem N
   real(r8), pointer :: livecrootn(:)  ! (gN/m2) live coarse root N
   ! ecophysiological constants
   real(r8), pointer :: woody(:)       ! binary flag for woody lifeform (1=woody, 0=not woody)
   real(r8), pointer :: livewdcn(:)    ! live wood (phloem and ray parenchyma) C:N (gC/gN)
   real(r8), pointer :: deadwdcn(:)    ! dead wood (xylem and heartwood) C:N (gC/gN)
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: livestemc_to_deadstemc(:)
   real(r8), pointer :: livecrootc_to_deadcrootc(:)
   real(r8), pointer :: livestemn_to_deadstemn(:)
   real(r8), pointer :: livestemn_to_retransn(:)
   real(r8), pointer :: livecrootn_to_deadcrootn(:)
   real(r8), pointer :: livecrootn_to_retransn(:)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
   integer :: p            ! indices
   integer :: fp           ! lake filter pft index
   real(r8):: ctovr        ! temporary variable for carbon turnover
   real(r8):: ntovr        ! temporary variable for nitrogen turnover

!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => pft%itype
    livestemc                      => pcs%livestemc
    livecrootc                     => pcs%livecrootc
    livestemn                      => pns%livestemn
    livecrootn                     => pns%livecrootn
    woody                          => pftcon%woody
    livewdcn                       => pftcon%livewdcn
    deadwdcn                       => pftcon%deadwdcn

   ! assign local pointers to derived type arrays (out)
    livestemc_to_deadstemc         => pcf%livestemc_to_deadstemc
    livecrootc_to_deadcrootc       => pcf%livecrootc_to_deadcrootc
    livestemn_to_deadstemn         => pnf%livestemn_to_deadstemn
    livestemn_to_retransn          => pnf%livestemn_to_retransn
    livecrootn_to_deadcrootn       => pnf%livecrootn_to_deadcrootn
    livecrootn_to_retransn         => pnf%livecrootn_to_retransn

   ! pft loop
   do fp = 1,num_soilp
      p = filter_soilp(fp)

      ! only calculate these fluxes for woody types
      if (woody(ivt(p)) > 0._r8) then

         ! live stem to dead stem turnover

         ctovr = livestemc(p) * lwtop
         ntovr = ctovr / livewdcn(ivt(p))
         livestemc_to_deadstemc(p) = ctovr
         livestemn_to_deadstemn(p) = ctovr / deadwdcn(ivt(p))
         livestemn_to_retransn(p)  = ntovr - livestemn_to_deadstemn(p)

         ! live coarse root to dead coarse root turnover

         ctovr = livecrootc(p) * lwtop
         ntovr = ctovr / livewdcn(ivt(p))
         livecrootc_to_deadcrootc(p) = ctovr
         livecrootn_to_deadcrootn(p) = ctovr / deadwdcn(ivt(p))
         livecrootn_to_retransn(p)  = ntovr - livecrootn_to_deadcrootn(p)

      end if

   end do

end subroutine CNLivewoodTurnover
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CNLitterToColumn
!
! !INTERFACE:
subroutine CNLitterToColumn (num_soilc, filter_soilc)
!
! !DESCRIPTION:
! called at the end of cn_phenology to gather all pft-level litterfall fluxes
! to the column level and assign them to the three litter pools
!
! !USES:
  use clm_varpar, only : max_pft_per_col
  use pftvarcon , only : npcropmin
!
! !ARGUMENTS:
  integer, intent(in) :: num_soilc       ! number of soil columns in filter
  integer, intent(in) :: filter_soilc(:) ! filter for soil columns
!
! !CALLED FROM:
! subroutine CNPhenology
!
! !REVISION HISTORY:
! 9/8/03: Created by Peter Thornton
!
! !LOCAL VARIABLES:
! local pointers to implicit in scalars
!
   integer , pointer :: ivt(:)          ! pft vegetation type
   real(r8), pointer :: wtcol(:)        ! weight (relative to column) for this pft (0-1)
   real(r8), pointer :: pwtgcell(:)     ! weight of pft relative to corresponding gridcell
   real(r8), pointer :: leafc_to_litter(:)     ! leaf C litterfall (gC/m2/s)
   real(r8), pointer :: frootc_to_litter(:)    ! fine root N litterfall (gN/m2/s)
   real(r8), pointer :: livestemc_to_litter(:) ! live stem C litterfall (gC/m2/s)
   real(r8), pointer :: grainc_to_food(:)      ! grain C to food (gC/m2/s)
   real(r8), pointer :: livestemn_to_litter(:) ! livestem N to litter (gN/m2/s)
   real(r8), pointer :: grainn_to_food(:)      ! grain N to food (gN/m2/s)
   real(r8), pointer :: leafn_to_litter(:)     ! leaf N litterfall (gN/m2/s)
   real(r8), pointer :: frootn_to_litter(:)    ! fine root N litterfall (gN/m2/s)
   real(r8), pointer :: lf_flab(:)      ! leaf litter labile fraction
   real(r8), pointer :: lf_fcel(:)      ! leaf litter cellulose fraction
   real(r8), pointer :: lf_flig(:)      ! leaf litter lignin fraction
   real(r8), pointer :: fr_flab(:)      ! fine root litter labile fraction
   real(r8), pointer :: fr_fcel(:)      ! fine root litter cellulose fraction
   real(r8), pointer :: fr_flig(:)      ! fine root litter lignin fraction
   integer , pointer :: npfts(:)        ! number of pfts for each column
   integer , pointer :: pfti(:)         ! beginning pft index for each column
!
! local pointers to implicit in/out scalars
!
   real(r8), pointer :: leafc_to_litr1c(:)     ! leaf C litterfall to litter 1 C (gC/m2/s)
   real(r8), pointer :: leafc_to_litr2c(:)     ! leaf C litterfall to litter 2 C (gC/m2/s)
   real(r8), pointer :: leafc_to_litr3c(:)     ! leaf C litterfall to litter 3 C (gC/m2/s)
   real(r8), pointer :: frootc_to_litr1c(:)    ! fine root C litterfall to litter 1 C (gC/m2/s)
   real(r8), pointer :: frootc_to_litr2c(:)    ! fine root C litterfall to litter 2 C (gC/m2/s)
   real(r8), pointer :: frootc_to_litr3c(:)    ! fine root C litterfall to litter 3 C (gC/m2/s)
   real(r8), pointer :: livestemc_to_litr1c(:) ! livestem C litterfall to litter 1 C (gC/m2/s)
   real(r8), pointer :: livestemc_to_litr2c(:) ! livestem C litterfall to litter 2 C (gC/m2/s)
   real(r8), pointer :: livestemc_to_litr3c(:) ! livestem C litterfall to litter 3 C (gC/m2/s)
   real(r8), pointer :: livestemn_to_litr1n(:) ! livestem N litterfall to litter 1 N (gN/m2/s)
   real(r8), pointer :: livestemn_to_litr2n(:) ! livestem N litterfall to litter 2 N (gN/m2/s)
   real(r8), pointer :: livestemn_to_litr3n(:) ! livestem N litterfall to litter 2 N (gN/m2/s)
   real(r8), pointer :: grainc_to_litr1c(:)    ! grain C litterfall to litter 1 C (gC/m2/s)
   real(r8), pointer :: grainc_to_litr2c(:)    ! grain C litterfall to litter 2 C (gC/m2/s)
   real(r8), pointer :: grainc_to_litr3c(:)    ! grain C litterfall to litter 3 C (gC/m2/s)
   real(r8), pointer :: grainn_to_litr1n(:)    ! grain N litterfall to litter 1 N (gN/m2/s)
   real(r8), pointer :: grainn_to_litr2n(:)    ! grain N litterfall to litter 2 N (gN/m2/s)
   real(r8), pointer :: grainn_to_litr3n(:)    ! grain N litterfall to litter 3 N (gN/m2/s)
   real(r8), pointer :: leafn_to_litr1n(:)     ! leaf N litterfall to litter 1 N (gN/m2/s)
   real(r8), pointer :: leafn_to_litr2n(:)     ! leaf N litterfall to litter 2 N (gN/m2/s)
   real(r8), pointer :: leafn_to_litr3n(:)     ! leaf N litterfall to litter 3 N (gN/m2/s)
   real(r8), pointer :: frootn_to_litr1n(:)    ! fine root N litterfall to litter 1 N (gN/m2/s)
   real(r8), pointer :: frootn_to_litr2n(:)    ! fine root N litterfall to litter 2 N (gN/m2/s)
   real(r8), pointer :: frootn_to_litr3n(:)    ! fine root N litterfall to litter 3 N (gN/m2/s)
!
! local pointers to implicit out scalars
!
!
! !OTHER LOCAL VARIABLES:
    integer :: fc,c,pi,p        ! indices
!EOP
!-----------------------------------------------------------------------
   ! assign local pointers to derived type arrays (in)
    ivt                            => pft%itype
    wtcol                          => pft%wtcol
    pwtgcell                       => pft%wtgcell  
    leafc_to_litter                => pcf%leafc_to_litter
    frootc_to_litter               => pcf%frootc_to_litter
    livestemc_to_litter            => pcf%livestemc_to_litter
    grainc_to_food                 => pcf%grainc_to_food
    livestemn_to_litter            => pnf%livestemn_to_litter
    grainn_to_food                 => pnf%grainn_to_food
    leafn_to_litter                => pnf%leafn_to_litter
    frootn_to_litter               => pnf%frootn_to_litter
    npfts                          => col%npfts
    pfti                           => col%pfti
    lf_flab                        => pftcon%lf_flab
    lf_fcel                        => pftcon%lf_fcel
    lf_flig                        => pftcon%lf_flig
    fr_flab                        => pftcon%fr_flab
    fr_fcel                        => pftcon%fr_fcel
    fr_flig                        => pftcon%fr_flig

   ! assign local pointers to derived type arrays (out)
    leafc_to_litr1c                => ccf%leafc_to_litr1c
    leafc_to_litr2c                => ccf%leafc_to_litr2c
    leafc_to_litr3c                => ccf%leafc_to_litr3c
    frootc_to_litr1c               => ccf%frootc_to_litr1c
    frootc_to_litr2c               => ccf%frootc_to_litr2c
    frootc_to_litr3c               => ccf%frootc_to_litr3c
    grainc_to_litr1c               => ccf%grainc_to_litr1c
    grainc_to_litr2c               => ccf%grainc_to_litr2c
    grainc_to_litr3c               => ccf%grainc_to_litr3c
    livestemc_to_litr1c            => ccf%livestemc_to_litr1c
    livestemc_to_litr2c            => ccf%livestemc_to_litr2c
    livestemc_to_litr3c            => ccf%livestemc_to_litr3c
    livestemn_to_litr1n            => cnf%livestemn_to_litr1n
    livestemn_to_litr2n            => cnf%livestemn_to_litr2n
    livestemn_to_litr3n            => cnf%livestemn_to_litr3n
    grainn_to_litr1n               => cnf%grainn_to_litr1n
    grainn_to_litr2n               => cnf%grainn_to_litr2n
    grainn_to_litr3n               => cnf%grainn_to_litr3n
    leafn_to_litr1n                => cnf%leafn_to_litr1n
    leafn_to_litr2n                => cnf%leafn_to_litr2n
    leafn_to_litr3n                => cnf%leafn_to_litr3n
    frootn_to_litr1n               => cnf%frootn_to_litr1n
    frootn_to_litr2n               => cnf%frootn_to_litr2n
    frootn_to_litr3n               => cnf%frootn_to_litr3n

   do pi = 1,max_pft_per_col
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if ( pi <=  npfts(c) ) then
            p = pfti(c) + pi - 1
            if (pwtgcell(p)>0._r8) then

               ! leaf litter carbon fluxes
               leafc_to_litr1c(c) = leafc_to_litr1c(c) + leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               leafc_to_litr2c(c) = leafc_to_litr2c(c) + leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               leafc_to_litr3c(c) = leafc_to_litr3c(c) + leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! leaf litter nitrogen fluxes
               leafn_to_litr1n(c) = leafn_to_litr1n(c) + leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
               leafn_to_litr2n(c) = leafn_to_litr2n(c) + leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
               leafn_to_litr3n(c) = leafn_to_litr3n(c) + leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

               ! fine root litter carbon fluxes
               frootc_to_litr1c(c) = frootc_to_litr1c(c) + frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               frootc_to_litr2c(c) = frootc_to_litr2c(c) + frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               frootc_to_litr3c(c) = frootc_to_litr3c(c) + frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)

               ! fine root litter nitrogen fluxes
               frootn_to_litr1n(c) = frootn_to_litr1n(c) + frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p)
               frootn_to_litr2n(c) = frootn_to_litr2n(c) + frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p)
               frootn_to_litr3n(c) = frootn_to_litr3n(c) + frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p)


               ! agroibis puts crop stem litter together with leaf litter
               ! so I've used the leaf lf_f* parameters instead of making
               ! new ones for now (slevis)
               ! also for simplicity I've put "food" into the litter pools
               if (ivt(p) >= npcropmin) then ! add livestemc to litter
                  ! stem litter carbon fluxes
                  livestemc_to_litr1c(c) = livestemc_to_litr1c(c) + livestemc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
                  livestemc_to_litr2c(c) = livestemc_to_litr2c(c) + livestemc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
                  livestemc_to_litr3c(c) = livestemc_to_litr3c(c) + livestemc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

                  ! stem litter nitrogen fluxes
                  livestemn_to_litr1n(c) = livestemn_to_litr1n(c) + livestemn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p)
                  livestemn_to_litr2n(c) = livestemn_to_litr2n(c) + livestemn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p)
                  livestemn_to_litr3n(c) = livestemn_to_litr3n(c) + livestemn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p)

                  ! grain litter carbon fluxes
                  grainc_to_litr1c(c) = grainc_to_litr1c(c) + grainc_to_food(p) * lf_flab(ivt(p)) * wtcol(p)
                  grainc_to_litr2c(c) = grainc_to_litr2c(c) + grainc_to_food(p) * lf_fcel(ivt(p)) * wtcol(p)
                  grainc_to_litr3c(c) = grainc_to_litr3c(c) + grainc_to_food(p) * lf_flig(ivt(p)) * wtcol(p)

                  ! grain litter nitrogen fluxes
                  grainn_to_litr1n(c) = grainn_to_litr1n(c) + grainn_to_food(p) * lf_flab(ivt(p)) * wtcol(p)
                  grainn_to_litr2n(c) = grainn_to_litr2n(c) + grainn_to_food(p) * lf_fcel(ivt(p)) * wtcol(p)
                  grainn_to_litr3n(c) = grainn_to_litr3n(c) + grainn_to_food(p) * lf_flig(ivt(p)) * wtcol(p)
               end if
            end if
         end if

      end do

   end do

end subroutine CNLitterToColumn
!-----------------------------------------------------------------------

end module CNPhenologyMod
