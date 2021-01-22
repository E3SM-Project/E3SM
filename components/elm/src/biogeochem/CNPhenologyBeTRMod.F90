module CNPhenologyBeTRMod
  !-----------------------------------------------------------------------
  ! !MODULE: CNPhenologyBeTRMod
  !
  ! !DESCRIPTION:
  ! Module holding routines used in phenology model for coupled carbon
  ! nitrogen code.

  !!Adding phosphorus -X.YANG
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_sys_mod         , only : shr_sys_flush
  use decompMod           , only : bounds_type
  use elm_varpar          , only : numpft
  use elm_varctl          , only : iulog
  use elm_varcon          , only : tfrz
  use abortutils          , only : endrun
  use CanopyStateType     , only : canopystate_type
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  use CNStateType         , only : cnstate_type
  use CropType            , only : crop_type
  use VegetationPropertiesType      , only : veg_vp
  use SoilStateType       , only : soilstate_type
  use TemperatureType     , only : temperature_type
  use WaterstateType      , only : waterstate_type
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use elm_varctl          , only : nu_com
  use CNBeTRIndicatorMod
  use GridcellType        , only : grc_pp
  use ColumnType          , only : col_pp
  use ColumnDataType      , only : col_es, col_ws, col_cf, col_nf, col_pf
  use TopounitDataType    , only : top_af, top_as
  use VegetationType      , only : veg_pp
  use VegetationDataType  , only : veg_es, veg_ef, veg_cs, veg_cf, veg_ns, veg_nf
  use VegetationDataType  , only : veg_ps, veg_pf
  !
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNPhenologyInit      ! Initialization
  public :: CNPhenology          ! Update
  public :: readCNPhenolBeTRParams   !
  !
  ! !PRIVATE DATA MEMBERS:
  type, private :: CNPnenolParamsType
     real(r8) :: crit_dayl       ! critical day length for senescence
     real(r8) :: ndays_on     	 ! number of days to complete leaf onset
     real(r8) :: ndays_off	 ! number of days to complete leaf offset
     real(r8) :: fstor2tran      ! fraction of storage to move to transfer for each onset
     real(r8) :: crit_onset_fdd  ! critical number of freezing days to set gdd counter
     real(r8) :: crit_onset_swi  ! critical number of days > soilpsi_on for onset
     real(r8) :: soilpsi_on      ! critical soil water potential for leaf onset
     real(r8) :: crit_offset_fdd ! critical number of freezing days to initiate offset
     real(r8) :: crit_offset_swi ! critical number of water stress days to initiate offset
     real(r8) :: soilpsi_off     ! critical soil water potential for leaf offset
     real(r8) :: lwtop   	 ! live wood turnover proportion (annual fraction)
  end type CNPnenolParamsType

  ! CNPhenolParamsInst is populated in readCNPhenolParams
  type(CNPnenolParamsType) ::  CNPhenolParamsInst

  real(r8) :: dt                            ! radiation time step delta t (seconds)
  real(r8) :: fracday                       ! dtime as a fraction of day
  real(r8) :: crit_dayl                     ! critical daylength for offset (seconds)
  real(r8) :: ndays_on                      ! number of days to complete onset
  real(r8) :: ndays_off                     ! number of days to complete offset
  real(r8) :: fstor2tran                    ! fraction of storage to move to transfer on each onset
  real(r8) :: crit_onset_fdd                ! critical number of freezing days
  real(r8) :: crit_onset_swi                ! water stress days for offset trigger
  real(r8) :: soilpsi_on                    ! water potential for onset trigger (MPa)
  real(r8) :: crit_offset_fdd               ! critical number of freezing degree days to trigger offset
  real(r8) :: crit_offset_swi               ! water stress days for offset trigger
  real(r8) :: soilpsi_off                   ! water potential for offset trigger (MPa)
  real(r8) :: lwtop                         ! live wood turnover proportion (annual fraction)

  ! CropPhenology variables and constants
  real(r8) :: p1d, p1v                      ! photoperiod factor constants for crop vernalization
  real(r8) :: hti                           ! cold hardening index threshold for vernalization
  real(r8) :: tbase                         ! base temperature for vernalization

  integer, parameter :: NOT_Planted   = 999 ! If not planted   yet in year
  integer, parameter :: NOT_Harvested = 999 ! If not harvested yet in year
  integer, parameter :: inNH       = 1      ! Northern Hemisphere
  integer, parameter :: inSH       = 2      ! Southern Hemisphere
  integer, pointer   :: inhemi(:)           ! Hemisphere that pft is in

  integer, allocatable :: minplantjday(:,:) ! minimum planting julian day
  integer, allocatable :: maxplantjday(:,:) ! maximum planting julian day
  integer              :: jdayyrstart(inSH) ! julian day of start of year
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine readCNPhenolBeTRParams ( ncid )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t,ncd_io

    ! !ARGUMENTS:
    implicit none
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNPhenolParamsType'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in parameter
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    !
    ! read in parameters
    !
    tString='crit_dayl'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%crit_dayl=tempr

    tString='ndays_on'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%ndays_on=tempr

    tString='ndays_off'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%ndays_off=tempr

    tString='fstor2tran'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%fstor2tran=tempr

    tString='crit_onset_fdd'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%crit_onset_fdd=tempr

    tString='crit_onset_swi'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%crit_onset_swi=tempr

    tString='soilpsi_on'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%soilpsi_on=tempr

    tString='crit_offset_fdd'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%crit_offset_fdd=tempr

    tString='crit_offset_swi'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%crit_offset_swi=tempr

    tString='soilpsi_off'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%soilpsi_off=tempr

    tString='lwtop_ann'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun( msg=trim(errCode)//trim(tString)//errMsg(__FILE__, __LINE__))
    CNPhenolParamsInst%lwtop=tempr

    call set_pheno_indicators
  end subroutine readCNPhenolBeTRParams

  !-----------------------------------------------------------------------
  subroutine CNPhenology (num_soilc, filter_soilc, num_soilp, filter_soilp, &
       num_pcropp, filter_pcropp, doalb, &
       waterstate_vars, temperature_vars, crop_vars, canopystate_vars, soilstate_vars, &
       cnstate_vars, carbonstate_vars, carbonflux_vars, &
       nitrogenstate_vars,nitrogenflux_vars,phosphorusstate_vars,phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! Dynamic phenology routine for coupled carbon-nitrogen code (CN)
    ! 1. grass phenology
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                  , intent(in)    :: num_pcropp      ! number of prog. crop patches in filter
    integer                  , intent(in)    :: filter_pcropp(:)! filter for prognostic crop patches
    logical                  , intent(in)    :: doalb           ! true if time for sfc albedo calc
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(temperature_type)   , intent(inout) :: temperature_vars
    type(crop_type)          , intent(inout) :: crop_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars

    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !-----------------------------------------------------------------------

    ! each of the following phenology type routines includes a filter
    ! to operate only on the relevant patches

    call CNPhenologyClimate(num_soilp, filter_soilp, num_pcropp, filter_pcropp, &
         temperature_vars, cnstate_vars,crop_vars)

    call CNEvergreenPhenology(num_soilp, filter_soilp, &
         cnstate_vars)

    call CNSeasonDecidPhenology(num_soilp, filter_soilp, &
         temperature_vars, cnstate_vars,  &
         carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
         phosphorusstate_vars,phosphorusflux_vars)

    call CNStressDecidPhenology(num_soilp, filter_soilp,   &
         soilstate_vars, temperature_vars, cnstate_vars, &
         carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
         phosphorusstate_vars,phosphorusflux_vars)

    if (num_pcropp > 0 ) then
       call CropPlantDate(num_soilp, filter_soilp, num_pcropp, filter_pcropp,&
             temperature_vars, cnstate_vars, crop_vars)
    end if

    if (doalb .and. num_pcropp > 0 ) then
       call CropPhenology(num_pcropp, filter_pcropp, &
            waterstate_vars, temperature_vars, crop_vars, canopystate_vars, cnstate_vars, &
            carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
            phosphorusstate_vars,phosphorusflux_vars)
    end if

    ! the same onset and offset routines are called regardless of
    ! phenology type - they depend only on onset_flag, offset_flag, bglfr, and bgtr

    call CNOnsetGrowth(num_soilp, filter_soilp, &
         cnstate_vars, &
         carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
         phosphorusstate_vars,phosphorusflux_vars)

   if (num_pcropp > 0 ) then
      call CNCropHarvest(num_pcropp, filter_pcropp, &
           num_soilc, filter_soilc, crop_vars, &
           cnstate_vars, carbonstate_vars, carbonflux_vars, nitrogenstate_vars, &
           nitrogenflux_vars, phosphorusstate_vars, phosphorusflux_vars)
   end if

    call CNOffsetLitterfall(num_soilp, filter_soilp, &
         cnstate_vars, carbonstate_vars, carbonflux_vars, nitrogenflux_vars,&
         phosphorusflux_vars, nitrogenstate_vars,phosphorusstate_vars)

    call CNBackgroundLitterfall(num_soilp, filter_soilp, &
         cnstate_vars, carbonstate_vars, carbonflux_vars, nitrogenflux_vars,&
         phosphorusflux_vars, nitrogenstate_vars, phosphorusstate_vars)

    call CNLivewoodTurnover(num_soilp, filter_soilp, &
         carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
         phosphorusstate_vars,phosphorusflux_vars)

    ! gather all patch-level litterfall fluxes to the column for litter C and N inputs

    call CNLitterToColumn(num_soilc, filter_soilc, &
         cnstate_vars, carbonflux_vars, nitrogenflux_vars,phosphorusflux_vars)

  end subroutine CNPhenology

  !-----------------------------------------------------------------------
  subroutine CNPhenologyInit(bounds)
    !
    ! !DESCRIPTION:
    ! Initialization of CNPhenology. Must be called after time-manager is
    ! initialized, and after ecophyscon file is read in.
    !
    ! !USES:
    use clm_time_manager, only: get_step_size
    use elm_varpar      , only: crop_prog
    use elm_varcon      , only: secspday
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    !------------------------------------------------------------------------

    !
    ! Get time-step and what fraction of a day it is
    !
    dt      = real( get_step_size(), r8 )
    fracday = dt/secspday

    ! set constants for CNSeasonDecidPhenology
    ! (critical daylength from Biome-BGC, v4.1.2)
    crit_dayl=CNPhenolParamsInst%crit_dayl

    ! Set constants for CNSeasonDecidPhenology and CNStressDecidPhenology
    ndays_on=CNPhenolParamsInst%ndays_on
    ndays_off=CNPhenolParamsInst%ndays_off

    ! set transfer parameters
    fstor2tran=CNPhenolParamsInst%fstor2tran

    ! -----------------------------------------
    ! Constants for CNStressDecidPhenology
    ! -----------------------------------------

    ! onset parameters
    crit_onset_fdd=CNPhenolParamsInst%crit_onset_fdd
    ! critical onset gdd now being calculated as a function of annual
    ! average 2m temp.
    ! crit_onset_gdd = 150.0 ! c3 grass value
    ! crit_onset_gdd = 1000.0   ! c4 grass value
    crit_onset_swi=CNPhenolParamsInst%crit_onset_swi
    soilpsi_on=CNPhenolParamsInst%soilpsi_on

    ! offset parameters
    crit_offset_fdd=CNPhenolParamsInst%crit_offset_fdd
    crit_offset_swi=CNPhenolParamsInst%crit_offset_swi
    soilpsi_off=CNPhenolParamsInst%soilpsi_off

    ! -----------------------------------------
    ! Constants for CNLivewoodTurnover
    ! -----------------------------------------

    ! set the global parameter for livewood turnover rate
    ! define as an annual fraction (0.7), and convert to fraction per second
    lwtop=CNPhenolParamsInst%lwtop/31536000.0_r8 !annual fraction converted to per second

    ! -----------------------------------------
    ! Call any subroutine specific initialization routines
    ! -----------------------------------------

    if ( crop_prog ) call CropPhenologyInit(bounds)

  end subroutine CNPhenologyInit

  !-----------------------------------------------------------------------
  subroutine CNPhenologyClimate (num_soilp, filter_soilp, num_pcropp, filter_pcropp, &
       temperature_vars, cnstate_vars, crop_vars)
    !
    ! !DESCRIPTION:
    ! For coupled carbon-nitrogen code (CN).
    !
    ! !USES:
    use clm_time_manager , only : get_days_per_year
    use clm_time_manager , only : get_curr_date, is_first_step
    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                , intent(in)    :: num_pcropp      ! number of prognostic crops in filter
    integer                , intent(in)    :: filter_pcropp(:)! filter for prognostic crop patches
    type(temperature_type) , intent(inout) :: temperature_vars
    type(cnstate_type)     , intent(inout) :: cnstate_vars
    type(crop_type)        , intent(inout) :: crop_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p                    ! indices
    integer :: fp                   ! lake filter pft index
    real(r8):: dayspyr              ! days per year (days)
    integer kyr                     ! current year
    integer kmo                     !         month of year  (1, ..., 12)
    integer kda                     !         day of month   (1, ..., 31)
    integer mcsec                   !         seconds of day (0, ..., seconds/day)
    real(r8), parameter :: yravg   = 20.0_r8      ! length of years to average for gdd
    real(r8), parameter :: yravgm1 = yravg-1.0_r8 ! minus 1 of above
    !-----------------------------------------------------------------------

    associate(                                                  &
         nyrs_crop_active => crop_vars%nyrs_crop_active_patch,   & ! InOut: [integer (:)  ]  number of years this crop patch has been active
         t_ref2m        => veg_es%t_ref2m     , & ! Input:  [real(r8) (:) ]  2m air temperature (K)
         gdd0           => veg_es%gdd0        , & ! Output: [real(r8) (:) ]  growing deg. days base 0 deg C (ddays)
         gdd8           => veg_es%gdd8        , & ! Output: [real(r8) (:) ]     "     "    "    "   8  "  "    "
         gdd10          => veg_es%gdd10       , & ! Output: [real(r8) (:) ]     "     "    "    "  10  "  "    "
         gdd020         => veg_es%gdd020      , & ! Output: [real(r8) (:) ]  20-yr mean of gdd0 (ddays)
         gdd820         => veg_es%gdd820      , & ! Output: [real(r8) (:) ]  20-yr mean of gdd8 (ddays)
         gdd1020        => veg_es%gdd1020     , & ! Output: [real(r8) (:) ]  20-yr mean of gdd10 (ddays)

         tempavg_t2m    => cnstate_vars%tempavg_t2m_patch       & ! Output: [real(r8) (:) ]  temp. avg 2m air temperature (K)
         )

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
      end if

      do fp = 1,num_pcropp
         p = filter_pcropp(fp)
         if (kmo == 1 .and. kda == 1 .and. nyrs_crop_active(p) == 0) then ! YR 1:
            gdd020(p)  = 0._r8                             ! set gdd..20 variables to 0
            gdd820(p)  = 0._r8                             ! and crops will not be planted
            gdd1020(p) = 0._r8
         end if
         if (kmo == 1 .and. kda == 1 .and. mcsec == 0) then        ! <-- END of EVERY YR:
            if (nyrs_crop_active(p)  == 1) then                    ! <-- END of YR 1
               gdd020(p)  = gdd0(p)                                ! <-- END of YR 1
               gdd820(p)  = gdd8(p)                                ! <-- END of YR 1
               gdd1020(p) = gdd10(p)                               ! <-- END of YR 1
            end if                                                 ! <-- END of YR 1
            gdd020(p)  = (yravgm1* gdd020(p)  + gdd0(p))  / yravg  ! gdd..20 must be long term avgs
            gdd820(p)  = (yravgm1* gdd820(p)  + gdd8(p))  / yravg  ! so ignore results for yrs 1 & 2
            gdd1020(p) = (yravgm1* gdd1020(p) + gdd10(p)) / yravg
         end if
      end do

    end associate

  end subroutine CNPhenologyClimate

  !-----------------------------------------------------------------------
  subroutine CNEvergreenPhenology (num_soilp, filter_soilp , &
       cnstate_vars)
    !
    ! !DESCRIPTION:
    ! For coupled carbon-nitrogen code (CN).
    !
    ! !USES:
    use elm_varcon       , only : secspday
    use clm_time_manager , only : get_days_per_year
    !
    ! !ARGUMENTS:
    integer           , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer           , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type), intent(inout) :: cnstate_vars
    !
    ! !LOCAL VARIABLES:
    real(r8):: dayspyr                ! Days per year
    integer :: p                      ! indices
    integer :: fp                     ! lake filter pft index
    !-----------------------------------------------------------------------

    associate(                                    &
         ivt        => veg_pp%itype                , & ! Input:  [integer  (:) ]  pft vegetation type

         evergreen  => veg_vp%evergreen     , & ! Input:  [real(r8) (:) ]  binary flag for evergreen leaf habit (0 or 1)
         leaf_long  => veg_vp%leaf_long     , & ! Input:  [real(r8) (:) ]  leaf longevity (yrs)

         bglfr      => cnstate_vars%bglfr_patch , & ! Output: [real(r8) (:) ]  background litterfall rate (1/s)
         bgtr       => cnstate_vars%bgtr_patch  , & ! Output: [real(r8) (:) ]  background transfer growth rate (1/s)
         lgsf       => cnstate_vars%lgsf_patch    & ! Output: [real(r8) (:) ]  long growing season factor [0-1]
         )

      dayspyr   = get_days_per_year()

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         if (evergreen(ivt(p)) == 1._r8) then
            bglfr(p) = 1._r8/(leaf_long(ivt(p)) * dayspyr * secspday)
            bgtr(p)  = 0._r8
            lgsf(p)  = 0._r8
         end if
      end do

    end associate

  end subroutine CNEvergreenPhenology

  !-----------------------------------------------------------------------
  subroutine CNSeasonDecidPhenology (num_soilp, filter_soilp       , &
       temperature_vars, cnstate_vars, &
       carbonstate_vars, nitrogenstate_vars, carbonflux_vars, nitrogenflux_vars,&
       phosphorusstate_vars, phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! For coupled carbon-nitrogen code (CN).
    ! This routine handles the seasonal deciduous phenology code (temperate
    ! deciduous vegetation that has only one growing season per year).
    !
    ! !USES:
    use shr_const_mod   , only: SHR_CONST_TKFRZ, SHR_CONST_PI
    use elm_varcon      , only: secspday
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: g,c,p          !indices
    integer :: fp             !lake filter pft index
    real(r8):: ws_flag        !winter-summer solstice flag (0 or 1)
    real(r8):: crit_onset_gdd !critical onset growing degree-day sum
    real(r8):: soilt
    !-----------------------------------------------------------------------

    associate(                                                                                             &
         ivt                                 =>    veg_pp%itype                                             , & ! Input:  [integer   (:)   ]  pft vegetation type
         dayl                                =>    grc_pp%dayl                                              , & ! Input:  [real(r8)  (:)   ]  daylength (s)
         prev_dayl                           =>    grc_pp%prev_dayl                                         , & ! Input:  [real(r8)  (:)   ]  daylength from previous time step (s)

         season_decid                        =>    veg_vp%season_decid                               , & ! Input:  [real(r8)  (:)   ]  binary flag for seasonal-deciduous leaf habit (0 or 1)
         woody                               =>    veg_vp%woody                                      , & ! Input:  [real(r8)  (:)   ]  binary flag for woody lifeform (1=woody, 0=not woody)

         t_soisno                            =>    col_es%t_soisno                         , & ! Input:  [real(r8)  (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

         annavg_t2m                          =>    cnstate_vars%annavg_t2m_patch                         , & ! Input:  [real(r8)  (:)   ]  annual average 2m air temperature (K)
         dormant_flag                        =>    cnstate_vars%dormant_flag_patch                       , & ! Output: [real(r8)  (:)   ]  dormancy flag
         days_active                         =>    cnstate_vars%days_active_patch                        , & ! Output: [real(r8)  (:)   ]  number of days since last dormancy
         onset_flag                          =>    cnstate_vars%onset_flag_patch                         , & ! Output: [real(r8)  (:)   ]  onset flag
         onset_counter                       =>    cnstate_vars%onset_counter_patch                      , & ! Output: [real(r8)  (:)   ]  onset counter (seconds)
         onset_gddflag                       =>    cnstate_vars%onset_gddflag_patch                      , & ! Output: [real(r8)  (:)   ]  onset freeze flag
         onset_gdd                           =>    cnstate_vars%onset_gdd_patch                          , & ! Output: [real(r8)  (:)   ]  onset growing degree days
         offset_flag                         =>    cnstate_vars%offset_flag_patch                        , & ! Output: [real(r8)  (:)   ]  offset flag
         offset_counter                      =>    cnstate_vars%offset_counter_patch                     , & ! Output: [real(r8)  (:)   ]  offset counter (seconds)
         bglfr                               =>    cnstate_vars%bglfr_patch                              , & ! Output: [real(r8)  (:)   ]  background litterfall rate (1/s)
         bgtr                                =>    cnstate_vars%bgtr_patch                               , & ! Output: [real(r8)  (:)   ]  background transfer growth rate (1/s)
         lgsf                                =>    cnstate_vars%lgsf_patch                               , & ! Output: [real(r8)  (:)   ]  long growing season factor [0-1]

         leafc_storage                       =>    veg_cs%leafc_storage                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) leaf C storage
         frootc_storage                      =>    veg_cs%frootc_storage                 , & ! Input:  [real(r8)  (:)   ]  (gC/m2) fine root C storage
         livestemc_storage                   =>    veg_cs%livestemc_storage              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live stem C storage
         deadstemc_storage                   =>    veg_cs%deadstemc_storage              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead stem C storage
         livecrootc_storage                  =>    veg_cs%livecrootc_storage             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live coarse root C storage
         deadcrootc_storage                  =>    veg_cs%deadcrootc_storage             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead coarse root C storage
         gresp_storage                       =>    veg_cs%gresp_storage                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) growth respiration storage
         leafc_xfer                          =>    veg_cs%leafc_xfer                     , & ! Output:  [real(r8) (:)   ]  (gC/m2) leaf C transfer
         frootc_xfer                         =>    veg_cs%frootc_xfer                    , & ! Output:  [real(r8) (:)   ]  (gC/m2) fine root C transfer
         livestemc_xfer                      =>    veg_cs%livestemc_xfer                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) live stem C transfer
         deadstemc_xfer                      =>    veg_cs%deadstemc_xfer                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead stem C transfer
         livecrootc_xfer                     =>    veg_cs%livecrootc_xfer                , & ! Output:  [real(r8) (:)   ]  (gC/m2) live coarse root C transfer
         deadcrootc_xfer                     =>    veg_cs%deadcrootc_xfer                , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead coarse root C transfer

         leafn_storage                       =>    veg_ns%leafn_storage                , & ! Input:  [real(r8)  (:)   ]  (gN/m2) leaf N storage
         frootn_storage                      =>    veg_ns%frootn_storage               , & ! Input:  [real(r8)  (:)   ]  (gN/m2) fine root N storage
         livestemn_storage                   =>    veg_ns%livestemn_storage            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live stem N storage
         deadstemn_storage                   =>    veg_ns%deadstemn_storage            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead stem N storage
         livecrootn_storage                  =>    veg_ns%livecrootn_storage           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live coarse root N storage
         deadcrootn_storage                  =>    veg_ns%deadcrootn_storage           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead coarse root N storage

         !! phosphorus
         leafp_storage                       =>    veg_ps%leafp_storage                , & ! Input:  [real(r8)  (:)   ]  (gP/m2) leaf P storage
         frootp_storage                      =>    veg_ps%frootp_storage               , & ! Input:  [real(r8)  (:)   ]  (gP/m2) fine root P storage
         livestemp_storage                   =>    veg_ps%livestemp_storage            , & ! Input:  [real(r8)  (:)   ]  (gP/m2) live stem P storage
         deadstemp_storage                   =>    veg_ps%deadstemp_storage            , & ! Input:  [real(r8)  (:)   ]  (gP/m2) dead stem P storage
         livecrootp_storage                  =>    veg_ps%livecrootp_storage           , & ! Input:  [real(r8)  (:)   ]  (gP/m2) live coarse root P storage
         deadcrootp_storage                  =>    veg_ps%deadcrootp_storage           , & ! Input:  [real(r8)  (:)   ]  (gP/m2) dead coarse root P storage

         prev_leafc_to_litter                =>    veg_cf%prev_leafc_to_litter            , & ! Output: [real(r8)  (:)   ]  previous timestep leaf C litterfall flux (gC/m2/s)
         prev_frootc_to_litter               =>    veg_cf%prev_frootc_to_litter           , & ! Output: [real(r8)  (:)   ]  previous timestep froot C litterfall flux (gC/m2/s)
         leafc_xfer_to_leafc                 =>    veg_cf%leafc_xfer_to_leafc             , & ! Output:  [real(r8) (:)   ]
         frootc_xfer_to_frootc               =>    veg_cf%frootc_xfer_to_frootc           , & ! Output:  [real(r8) (:)   ]
         livestemc_xfer_to_livestemc         =>    veg_cf%livestemc_xfer_to_livestemc     , & ! Output:  [real(r8) (:)   ]
         deadstemc_xfer_to_deadstemc         =>    veg_cf%deadstemc_xfer_to_deadstemc     , & ! Output:  [real(r8) (:)   ]
         livecrootc_xfer_to_livecrootc       =>    veg_cf%livecrootc_xfer_to_livecrootc   , & ! Output:  [real(r8) (:)   ]
         deadcrootc_xfer_to_deadcrootc       =>    veg_cf%deadcrootc_xfer_to_deadcrootc   , & ! Output:  [real(r8) (:)   ]
         leafc_storage_to_xfer               =>    veg_cf%leafc_storage_to_xfer           , & ! Output:  [real(r8) (:)   ]
         frootc_storage_to_xfer              =>    veg_cf%frootc_storage_to_xfer          , & ! Output:  [real(r8) (:)   ]
         livestemc_storage_to_xfer           =>    veg_cf%livestemc_storage_to_xfer       , & ! Output:  [real(r8) (:)   ]
         deadstemc_storage_to_xfer           =>    veg_cf%deadstemc_storage_to_xfer       , & ! Output:  [real(r8) (:)   ]
         livecrootc_storage_to_xfer          =>    veg_cf%livecrootc_storage_to_xfer      , & ! Output:  [real(r8) (:)   ]
         deadcrootc_storage_to_xfer          =>    veg_cf%deadcrootc_storage_to_xfer      , & ! Output:  [real(r8) (:)   ]
         gresp_storage_to_xfer               =>    veg_cf%gresp_storage_to_xfer           , & ! Output:  [real(r8) (:)   ]

         leafn_xfer_to_leafn                 =>    veg_nf%leafn_xfer_to_leafn           , & ! Output:  [real(r8) (:)   ]
         frootn_xfer_to_frootn               =>    veg_nf%frootn_xfer_to_frootn         , & ! Output:  [real(r8) (:)   ]
         livestemn_xfer_to_livestemn         =>    veg_nf%livestemn_xfer_to_livestemn   , & ! Output:  [real(r8) (:)   ]
         deadstemn_xfer_to_deadstemn         =>    veg_nf%deadstemn_xfer_to_deadstemn   , & ! Output:  [real(r8) (:)   ]
         livecrootn_xfer_to_livecrootn       =>    veg_nf%livecrootn_xfer_to_livecrootn , & ! Output:  [real(r8) (:)   ]
         deadcrootn_xfer_to_deadcrootn       =>    veg_nf%deadcrootn_xfer_to_deadcrootn , & ! Output:  [real(r8) (:)   ]
         leafn_xfer                          =>    veg_ns%leafn_xfer                   , & ! Output:  [real(r8) (:)   ]  (gN/m2) leaf N transfer
         frootn_xfer                         =>    veg_ns%frootn_xfer                  , & ! Output:  [real(r8) (:)   ]  (gN/m2) fine root N transfer
         livestemn_xfer                      =>    veg_ns%livestemn_xfer               , & ! Output:  [real(r8) (:)   ]  (gN/m2) live stem N transfer
         deadstemn_xfer                      =>    veg_ns%deadstemn_xfer               , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead stem N transfer
         livecrootn_xfer                     =>    veg_ns%livecrootn_xfer              , & ! Output:  [real(r8) (:)   ]  (gN/m2) live coarse root N transfer
         deadcrootn_xfer                     =>    veg_ns%deadcrootn_xfer              , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead coarse root N transfer
         leafn_storage_to_xfer               =>    veg_nf%leafn_storage_to_xfer         , & ! Output:  [real(r8) (:)   ]
         frootn_storage_to_xfer              =>    veg_nf%frootn_storage_to_xfer        , & ! Output:  [real(r8) (:)   ]
         livestemn_storage_to_xfer           =>    veg_nf%livestemn_storage_to_xfer     , & ! Output:  [real(r8) (:)   ]
         deadstemn_storage_to_xfer           =>    veg_nf%deadstemn_storage_to_xfer     , & ! Output:  [real(r8) (:)   ]
         livecrootn_storage_to_xfer          =>    veg_nf%livecrootn_storage_to_xfer    , & ! Output:  [real(r8) (:)   ]
         deadcrootn_storage_to_xfer          =>    veg_nf%deadcrootn_storage_to_xfer    , & ! Output:  [real(r8) (:)   ]

         leafp_xfer_to_leafp                 =>    veg_pf%leafp_xfer_to_leafp           , & ! Output:  [real(r8) (:)   ]
         frootp_xfer_to_frootp               =>    veg_pf%frootp_xfer_to_frootp         , & ! Output:  [real(r8) (:)   ]
         livestemp_xfer_to_livestemp         =>    veg_pf%livestemp_xfer_to_livestemp   , & ! Output:  [real(r8) (:)   ]
         deadstemp_xfer_to_deadstemp         =>    veg_pf%deadstemp_xfer_to_deadstemp   , & ! Output:  [real(r8) (:)   ]
         livecrootp_xfer_to_livecrootp       =>    veg_pf%livecrootp_xfer_to_livecrootp , & ! Output:  [real(r8) (:)   ]
         deadcrootp_xfer_to_deadcrootp       =>    veg_pf%deadcrootp_xfer_to_deadcrootp , & ! Output:  [real(r8) (:)   ]
         leafp_xfer                          =>    veg_ps%leafp_xfer                   , & ! Output:  [real(r8) (:)   ]  (gP/m2) leaf P transfer
         frootp_xfer                         =>    veg_ps%frootp_xfer                  , & ! Output:  [real(r8) (:)   ]  (gP/m2) fine root P transfer
         livestemp_xfer                      =>    veg_ps%livestemp_xfer               , & ! Output:  [real(r8) (:)   ]  (gP/m2) live stem P transfer
         deadstemp_xfer                      =>    veg_ps%deadstemp_xfer               , & ! Output:  [real(r8) (:)   ]  (gP/m2) dead stem P transfer
         livecrootp_xfer                     =>    veg_ps%livecrootp_xfer              , & ! Output:  [real(r8) (:)   ]  (gP/m2) live coarse root P transfer
         deadcrootp_xfer                     =>    veg_ps%deadcrootp_xfer              , & ! Output:  [real(r8) (:)   ]  (gP/m2) dead coarse root P transfer
         leafp_storage_to_xfer               =>    veg_pf%leafp_storage_to_xfer         , & ! Output:  [real(r8) (:)   ]
         frootp_storage_to_xfer              =>    veg_pf%frootp_storage_to_xfer        , & ! Output:  [real(r8) (:)   ]
         livestemp_storage_to_xfer           =>    veg_pf%livestemp_storage_to_xfer     , & ! Output:  [real(r8) (:)   ]
         deadstemp_storage_to_xfer           =>    veg_pf%deadstemp_storage_to_xfer     , & ! Output:  [real(r8) (:)   ]
         livecrootp_storage_to_xfer          =>    veg_pf%livecrootp_storage_to_xfer    , & ! Output:  [real(r8) (:)   ]
         deadcrootp_storage_to_xfer          =>    veg_pf%deadcrootp_storage_to_xfer      & ! Output:  [real(r8) (:)   ]

         )

      ! start pft loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = veg_pp%column(p)
         g = veg_pp%gridcell(p)

         if (season_decid(ivt(p)) == 1._r8) then

            ! set background litterfall rate, background transfer rate, and
            ! long growing season factor to 0 for seasonal deciduous types
            bglfr(p) = 0._r8
            bgtr(p) = 0._r8
            lgsf(p) = 0._r8

            ! onset gdd sum from Biome-BGC, v4.1.2
            crit_onset_gdd = exp(4.8_r8 + 0.13_r8*(annavg_t2m(p) - SHR_CONST_TKFRZ))

            ! set flag for solstice period (winter->summer = 1, summer->winter = 0)
            if (dayl(g) >= prev_dayl(g)) then
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
                  leafp_xfer_to_leafp(p)   = 0.0_r8
                  frootp_xfer_to_frootp(p) = 0.0_r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer_to_livestemc(p)   = 0.0_r8
                     deadstemc_xfer_to_deadstemc(p)   = 0.0_r8
                     livecrootc_xfer_to_livecrootc(p) = 0.0_r8
                     deadcrootc_xfer_to_deadcrootc(p) = 0.0_r8
                     livestemn_xfer_to_livestemn(p)   = 0.0_r8
                     deadstemn_xfer_to_deadstemn(p)   = 0.0_r8
                     livecrootn_xfer_to_livecrootn(p) = 0.0_r8
                     deadcrootn_xfer_to_deadcrootn(p) = 0.0_r8
                     livestemp_xfer_to_livestemp(p)   = 0.0_r8
                     deadstemp_xfer_to_deadstemp(p)   = 0.0_r8
                     livecrootp_xfer_to_livecrootp(p) = 0.0_r8
                     deadcrootp_xfer_to_deadcrootp(p) = 0.0_r8
                  end if
                  ! set transfer pools to 0.0
                  leafc_xfer(p) = 0.0_r8
                  leafn_xfer(p) = 0.0_r8
                  leafp_xfer(p) = 0.0_r8
                  frootc_xfer(p) = 0.0_r8
                  frootn_xfer(p) = 0.0_r8
                  frootp_xfer(p) = 0.0_r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer(p) = 0.0_r8
                     livestemn_xfer(p) = 0.0_r8
                     livestemp_xfer(p) = 0.0_r8
                     deadstemc_xfer(p) = 0.0_r8
                     deadstemn_xfer(p) = 0.0_r8
                     deadstemp_xfer(p) = 0.0_r8
                     livecrootc_xfer(p) = 0.0_r8
                     livecrootn_xfer(p) = 0.0_r8
                     livecrootp_xfer(p) = 0.0_r8
                     deadcrootc_xfer(p) = 0.0_r8
                     deadcrootn_xfer(p) = 0.0_r8
                     deadcrootp_xfer(p) = 0.0_r8
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

                  ! set phosphorus fluxes for shifting storage pools to transfer pools
                  leafp_storage_to_xfer(p)  = fstor2tran * leafp_storage(p)/dt
                  frootp_storage_to_xfer(p) = fstor2tran * frootp_storage(p)/dt
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemp_storage_to_xfer(p)  = fstor2tran * livestemp_storage(p)/dt
                     deadstemp_storage_to_xfer(p)  = fstor2tran * deadstemp_storage(p)/dt
                     livecrootp_storage_to_xfer(p) = fstor2tran * livecrootp_storage(p)/dt
                     deadcrootp_storage_to_xfer(p) = fstor2tran * deadcrootp_storage(p)/dt
                  end if
               end if

               ! test for switching from growth period to offset period
            else if (offset_flag(p) == 0.0_r8) then
               ! only begin to test for offset daylength once past the summer sol
               if (ws_flag == 0._r8 .and. dayl(g) < crit_dayl) then
                  offset_flag(p) = 1._r8
                  offset_counter(p) = ndays_off * secspday
                  prev_leafc_to_litter(p) = 0._r8
                  prev_frootc_to_litter(p) = 0._r8
               end if
            end if

         end if ! end if seasonal deciduous

      end do ! end of pft loop

    end associate

  end subroutine CNSeasonDecidPhenology

  !-----------------------------------------------------------------------
  subroutine CNStressDecidPhenology (num_soilp, filter_soilp , &
       soilstate_vars, temperature_vars, cnstate_vars        , &
       carbonstate_vars, nitrogenstate_vars, carbonflux_vars,nitrogenflux_vars,&
       phosphorusstate_vars,phosphorusflux_vars)
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
    use clm_time_manager , only : get_days_per_year
    use elm_varcon       , only : secspday
    use shr_const_mod    , only : SHR_CONST_TKFRZ, SHR_CONST_PI
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(soilstate_type)     , intent(in)    :: soilstate_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    real(r8),parameter :: secspqtrday = secspday / 4  ! seconds per quarter day
    integer :: g,c,p           ! indices
    integer :: fp              ! lake filter pft index
    real(r8):: dayspyr         ! days per year
    real(r8):: crit_onset_gdd  ! degree days for onset trigger
    real(r8):: soilt           ! temperature of top soil layer
    real(r8):: psi             ! water stress of top soil layer
    !-----------------------------------------------------------------------

    associate(                                                                                             &
         ivt                                 =>    veg_pp%itype                                             , & ! Input:  [integer   (:)   ]  pft vegetation type
         dayl                                =>    grc_pp%dayl                                              , & ! Input:  [real(r8)  (:)   ]  daylength (s)

         leaf_long                           =>    veg_vp%leaf_long                                  , & ! Input:  [real(r8)  (:)   ]  leaf longevity (yrs)
         woody                               =>    veg_vp%woody                                      , & ! Input:  [real(r8)  (:)   ]  binary flag for woody lifeform (1=woody, 0=not woody)
         stress_decid                        =>    veg_vp%stress_decid                               , & ! Input:  [real(r8)  (:)   ]  binary flag for stress-deciduous leaf habit (0 or 1)

         soilpsi                             =>    soilstate_vars%soilpsi_col                            , & ! Input:  [real(r8)  (:,:) ]  soil water potential in each soil layer (MPa)

         t_soisno                            =>    col_es%t_soisno                         , & ! Input:  [real(r8)  (:,:) ]  soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)

         dormant_flag                        =>    cnstate_vars%dormant_flag_patch                       , & ! Output:  [real(r8) (:)   ]  dormancy flag
         days_active                         =>    cnstate_vars%days_active_patch                        , & ! Output:  [real(r8) (:)   ]  number of days since last dormancy
         onset_flag                          =>    cnstate_vars%onset_flag_patch                         , & ! Output:  [real(r8) (:)   ]  onset flag
         onset_counter                       =>    cnstate_vars%onset_counter_patch                      , & ! Output:  [real(r8) (:)   ]  onset counter (seconds)
         onset_gddflag                       =>    cnstate_vars%onset_gddflag_patch                      , & ! Output:  [real(r8) (:)   ]  onset freeze flag
         onset_fdd                           =>    cnstate_vars%onset_fdd_patch                          , & ! Output:  [real(r8) (:)   ]  onset freezing degree days counter
         onset_gdd                           =>    cnstate_vars%onset_gdd_patch                          , & ! Output:  [real(r8) (:)   ]  onset growing degree days
         onset_swi                           =>    cnstate_vars%onset_swi_patch                          , & ! Output:  [real(r8) (:)   ]  onset soil water index
         offset_flag                         =>    cnstate_vars%offset_flag_patch                        , & ! Output:  [real(r8) (:)   ]  offset flag
         offset_counter                      =>    cnstate_vars%offset_counter_patch                     , & ! Output:  [real(r8) (:)   ]  offset counter (seconds)
         offset_fdd                          =>    cnstate_vars%offset_fdd_patch                         , & ! Output:  [real(r8) (:)   ]  offset freezing degree days counter
         offset_swi                          =>    cnstate_vars%offset_swi_patch                         , & ! Output:  [real(r8) (:)   ]  offset soil water index
         lgsf                                =>    cnstate_vars%lgsf_patch                               , & ! Output:  [real(r8) (:)   ]  long growing season factor [0-1]
         bglfr                               =>    cnstate_vars%bglfr_patch                              , & ! Output:  [real(r8) (:)   ]  background litterfall rate (1/s)
         bgtr                                =>    cnstate_vars%bgtr_patch                               , & ! Output:  [real(r8) (:)   ]  background transfer growth rate (1/s)
         annavg_t2m                          =>    cnstate_vars%annavg_t2m_patch                         , & ! Output:  [real(r8) (:)   ]  annual average 2m air temperature (K)

         leafc_storage                       =>    veg_cs%leafc_storage                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) leaf C storage
         frootc_storage                      =>    veg_cs%frootc_storage                 , & ! Input:  [real(r8)  (:)   ]  (gC/m2) fine root C storage
         livestemc_storage                   =>    veg_cs%livestemc_storage              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live stem C storage
         deadstemc_storage                   =>    veg_cs%deadstemc_storage              , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead stem C storage
         livecrootc_storage                  =>    veg_cs%livecrootc_storage             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) live coarse root C storage
         deadcrootc_storage                  =>    veg_cs%deadcrootc_storage             , & ! Input:  [real(r8)  (:)   ]  (gC/m2) dead coarse root C storage
         gresp_storage                       =>    veg_cs%gresp_storage                  , & ! Input:  [real(r8)  (:)   ]  (gC/m2) growth respiration storage
         leafc_xfer                          =>    veg_cs%leafc_xfer                     , & ! Output:  [real(r8) (:)   ]  (gC/m2) leaf C transfer
         frootc_xfer                         =>    veg_cs%frootc_xfer                    , & ! Output:  [real(r8) (:)   ]  (gC/m2) fine root C transfer
         livestemc_xfer                      =>    veg_cs%livestemc_xfer                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) live stem C transfer
         deadstemc_xfer                      =>    veg_cs%deadstemc_xfer                 , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead stem C transfer
         livecrootc_xfer                     =>    veg_cs%livecrootc_xfer                , & ! Output:  [real(r8) (:)   ]  (gC/m2) live coarse root C transfer
         deadcrootc_xfer                     =>    veg_cs%deadcrootc_xfer                , & ! Output:  [real(r8) (:)   ]  (gC/m2) dead coarse root C transfer

         leafn_storage                       =>    veg_ns%leafn_storage                , & ! Input:  [real(r8)  (:)   ]  (gN/m2) leaf N storage
         frootn_storage                      =>    veg_ns%frootn_storage               , & ! Input:  [real(r8)  (:)   ]  (gN/m2) fine root N storage
         livestemn_storage                   =>    veg_ns%livestemn_storage            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live stem N storage
         deadstemn_storage                   =>    veg_ns%deadstemn_storage            , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead stem N storage
         livecrootn_storage                  =>    veg_ns%livecrootn_storage           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) live coarse root N storage
         deadcrootn_storage                  =>    veg_ns%deadcrootn_storage           , & ! Input:  [real(r8)  (:)   ]  (gN/m2) dead coarse root N storage
         leafn_xfer                          =>    veg_ns%leafn_xfer                   , & ! Output:  [real(r8) (:)   ]  (gN/m2) leaf N transfer
         frootn_xfer                         =>    veg_ns%frootn_xfer                  , & ! Output:  [real(r8) (:)   ]  (gN/m2) fine root N transfer
         livestemn_xfer                      =>    veg_ns%livestemn_xfer               , & ! Output:  [real(r8) (:)   ]  (gN/m2) live stem N transfer
         deadstemn_xfer                      =>    veg_ns%deadstemn_xfer               , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead stem N transfer
         livecrootn_xfer                     =>    veg_ns%livecrootn_xfer              , & ! Output:  [real(r8) (:)   ]  (gN/m2) live coarse root N transfer
         deadcrootn_xfer                     =>    veg_ns%deadcrootn_xfer              , & ! Output:  [real(r8) (:)   ]  (gN/m2) dead coarse root N transfer


         leafp_storage                       =>    veg_ps%leafp_storage                , & ! Input:  [real(r8)  (:)   ]  (gP/m2) leaf P storage
         frootp_storage                      =>    veg_ps%frootp_storage               , & ! Input:  [real(r8)  (:)   ]  (gP/m2) fine root P storage
         livestemp_storage                   =>    veg_ps%livestemp_storage            , & ! Input:  [real(r8)  (:)   ]  (gP/m2) live stem P storage
         deadstemp_storage                   =>    veg_ps%deadstemp_storage            , & ! Input:  [real(r8)  (:)   ]  (gP/m2) dead stem P storage
         livecrootp_storage                  =>    veg_ps%livecrootp_storage           , & ! Input:  [real(r8)  (:)   ]  (gP/m2) live coarse root P storage
         deadcrootp_storage                  =>    veg_ps%deadcrootp_storage           , & ! Input:  [real(r8)  (:)   ]  (gP/m2) dead coarse root P storage
         leafp_xfer                          =>    veg_ps%leafp_xfer                   , & ! Output:  [real(r8) (:)   ]  (gP/m2) leaf P transfer
         frootp_xfer                         =>    veg_ps%frootp_xfer                  , & ! Output:  [real(r8) (:)   ]  (gP/m2) fine root P transfer
         livestemp_xfer                      =>    veg_ps%livestemp_xfer               , & ! Output:  [real(r8) (:)   ]  (gP/m2) live stem P transfer
         deadstemp_xfer                      =>    veg_ps%deadstemp_xfer               , & ! Output:  [real(r8) (:)   ]  (gP/m2) dead stem P transfer
         livecrootp_xfer                     =>    veg_ps%livecrootp_xfer              , & ! Output:  [real(r8) (:)   ]  (gP/m2) live coarse root P transfer
         deadcrootp_xfer                     =>    veg_ps%deadcrootp_xfer              , & ! Output:  [real(r8) (:)   ]  (gP/m2) dead coarse root P transfer

         prev_leafc_to_litter                =>    veg_cf%prev_leafc_to_litter            , & ! Output:  [real(r8) (:)   ]  previous timestep leaf C litterfall flux (gC/m2/s)
         prev_frootc_to_litter               =>    veg_cf%prev_frootc_to_litter           , & ! Output:  [real(r8) (:)   ]  previous timestep froot C litterfall flux (gC/m2/s)
         leafc_xfer_to_leafc                 =>    veg_cf%leafc_xfer_to_leafc             , & ! Output:  [real(r8) (:)   ]
         frootc_xfer_to_frootc               =>    veg_cf%frootc_xfer_to_frootc           , & ! Output:  [real(r8) (:)   ]
         livestemc_xfer_to_livestemc         =>    veg_cf%livestemc_xfer_to_livestemc     , & ! Output:  [real(r8) (:)   ]
         deadstemc_xfer_to_deadstemc         =>    veg_cf%deadstemc_xfer_to_deadstemc     , & ! Output:  [real(r8) (:)   ]
         livecrootc_xfer_to_livecrootc       =>    veg_cf%livecrootc_xfer_to_livecrootc   , & ! Output:  [real(r8) (:)   ]
         deadcrootc_xfer_to_deadcrootc       =>    veg_cf%deadcrootc_xfer_to_deadcrootc   , & ! Output:  [real(r8) (:)   ]
         leafc_storage_to_xfer               =>    veg_cf%leafc_storage_to_xfer           , & ! Output:  [real(r8) (:)   ]
         frootc_storage_to_xfer              =>    veg_cf%frootc_storage_to_xfer          , & ! Output:  [real(r8) (:)   ]
         livestemc_storage_to_xfer           =>    veg_cf%livestemc_storage_to_xfer       , & ! Output:  [real(r8) (:)   ]
         deadstemc_storage_to_xfer           =>    veg_cf%deadstemc_storage_to_xfer       , & ! Output:  [real(r8) (:)   ]
         livecrootc_storage_to_xfer          =>    veg_cf%livecrootc_storage_to_xfer      , & ! Output:  [real(r8) (:)   ]
         deadcrootc_storage_to_xfer          =>    veg_cf%deadcrootc_storage_to_xfer      , & ! Output:  [real(r8) (:)   ]
         gresp_storage_to_xfer               =>    veg_cf%gresp_storage_to_xfer           , & ! Output:  [real(r8) (:)   ]

         leafn_xfer_to_leafn                 =>    veg_nf%leafn_xfer_to_leafn           , & ! Output:  [real(r8) (:)   ]
         frootn_xfer_to_frootn               =>    veg_nf%frootn_xfer_to_frootn         , & ! Output:  [real(r8) (:)   ]
         livestemn_xfer_to_livestemn         =>    veg_nf%livestemn_xfer_to_livestemn   , & ! Output:  [real(r8) (:)   ]
         deadstemn_xfer_to_deadstemn         =>    veg_nf%deadstemn_xfer_to_deadstemn   , & ! Output:  [real(r8) (:)   ]
         livecrootn_xfer_to_livecrootn       =>    veg_nf%livecrootn_xfer_to_livecrootn , & ! Output:  [real(r8) (:)   ]
         deadcrootn_xfer_to_deadcrootn       =>    veg_nf%deadcrootn_xfer_to_deadcrootn , & ! Output:  [real(r8) (:)   ]
         leafn_storage_to_xfer               =>    veg_nf%leafn_storage_to_xfer         , & ! Output:  [real(r8) (:)   ]
         frootn_storage_to_xfer              =>    veg_nf%frootn_storage_to_xfer        , & ! Output:  [real(r8) (:)   ]
         livestemn_storage_to_xfer           =>    veg_nf%livestemn_storage_to_xfer     , & ! Output:  [real(r8) (:)   ]
         deadstemn_storage_to_xfer           =>    veg_nf%deadstemn_storage_to_xfer     , & ! Output:  [real(r8) (:)   ]
         livecrootn_storage_to_xfer          =>    veg_nf%livecrootn_storage_to_xfer    , & ! Output:  [real(r8) (:)   ]
         deadcrootn_storage_to_xfer          =>    veg_nf%deadcrootn_storage_to_xfer    ,  & ! Output:  [real(r8) (:)   ]

         leafp_xfer_to_leafp                 =>    veg_pf%leafp_xfer_to_leafp           , & ! Output:  [real(r8) (:)   ]
         frootp_xfer_to_frootp               =>    veg_pf%frootp_xfer_to_frootp         , & ! Output:  [real(r8) (:)   ]
         livestemp_xfer_to_livestemp         =>    veg_pf%livestemp_xfer_to_livestemp   , & ! Output:  [real(r8) (:)   ]
         deadstemp_xfer_to_deadstemp         =>    veg_pf%deadstemp_xfer_to_deadstemp   , & ! Output:  [real(r8) (:)   ]
         livecrootp_xfer_to_livecrootp       =>    veg_pf%livecrootp_xfer_to_livecrootp , & ! Output:  [real(r8) (:)   ]
         deadcrootp_xfer_to_deadcrootp       =>    veg_pf%deadcrootp_xfer_to_deadcrootp , & ! Output:  [real(r8) (:)   ]
         leafp_storage_to_xfer               =>    veg_pf%leafp_storage_to_xfer         , & ! Output:  [real(r8) (:)   ]
         frootp_storage_to_xfer              =>    veg_pf%frootp_storage_to_xfer        , & ! Output:  [real(r8) (:)   ]
         livestemp_storage_to_xfer           =>    veg_pf%livestemp_storage_to_xfer     , & ! Output:  [real(r8) (:)   ]
         deadstemp_storage_to_xfer           =>    veg_pf%deadstemp_storage_to_xfer     , & ! Output:  [real(r8) (:)   ]
         livecrootp_storage_to_xfer          =>    veg_pf%livecrootp_storage_to_xfer    , & ! Output:  [real(r8) (:)   ]
         deadcrootp_storage_to_xfer          =>    veg_pf%deadcrootp_storage_to_xfer      & ! Output:  [real(r8) (:)   ]

         )

      ! set time steps
      dayspyr = get_days_per_year()

      do fp = 1,num_soilp
         p = filter_soilp(fp)
         c = veg_pp%column(p)
         g = veg_pp%gridcell(p)

         if (stress_decid(ivt(p)) == 1._r8) then
            soilt = t_soisno(c,3)
            psi = soilpsi(c,3)

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
                  leafp_xfer_to_leafp(p)   = 0._r8
                  frootp_xfer_to_frootp(p) = 0._r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer_to_livestemc(p)   = 0._r8
                     deadstemc_xfer_to_deadstemc(p)   = 0._r8
                     livecrootc_xfer_to_livecrootc(p) = 0._r8
                     deadcrootc_xfer_to_deadcrootc(p) = 0._r8
                     livestemn_xfer_to_livestemn(p)   = 0._r8
                     deadstemn_xfer_to_deadstemn(p)   = 0._r8
                     livecrootn_xfer_to_livecrootn(p) = 0._r8
                     deadcrootn_xfer_to_deadcrootn(p) = 0._r8
                     livestemp_xfer_to_livestemp(p)   = 0._r8
                     deadstemp_xfer_to_deadstemp(p)   = 0._r8
                     livecrootp_xfer_to_livecrootp(p) = 0._r8
                     deadcrootp_xfer_to_deadcrootp(p) = 0._r8
                  end if
                  ! set transfer pools to 0.0
                  leafc_xfer(p) = 0._r8
                  leafn_xfer(p) = 0._r8
                  leafp_xfer(p) = 0._r8
                  frootc_xfer(p) = 0._r8
                  frootn_xfer(p) = 0._r8
                  frootp_xfer(p) = 0._r8
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemc_xfer(p) = 0._r8
                     livestemn_xfer(p) = 0._r8
                     livestemp_xfer(p) = 0._r8
                     deadstemc_xfer(p) = 0._r8
                     deadstemn_xfer(p) = 0._r8
                     deadstemp_xfer(p) = 0._r8
                     livecrootc_xfer(p) = 0._r8
                     livecrootn_xfer(p) = 0._r8
                     livecrootp_xfer(p) = 0._r8
                     deadcrootc_xfer(p) = 0._r8
                     deadcrootn_xfer(p) = 0._r8
                     deadcrootp_xfer(p) = 0._r8
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
               if (onset_flag(p) == 1._r8 .and. dayl(g) <= secspqtrday) then
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

                  ! set phosphorus fluxes for shifting storage pools to transfer pools
                  leafp_storage_to_xfer(p)  = fstor2tran * leafp_storage(p)/dt
                  frootp_storage_to_xfer(p) = fstor2tran * frootp_storage(p)/dt
                  if (woody(ivt(p)) == 1.0_r8) then
                     livestemp_storage_to_xfer(p)  = fstor2tran * livestemp_storage(p)/dt
                     deadstemp_storage_to_xfer(p)  = fstor2tran * deadstemp_storage(p)/dt
                     livecrootp_storage_to_xfer(p) = fstor2tran * livecrootp_storage(p)/dt
                     deadcrootp_storage_to_xfer(p) = fstor2tran * deadcrootp_storage(p)/dt
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
               if (dayl(g) <= secspqtrday) then
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


               ! set phosphorus fluxes for shifting storage pools to transfer pools
               leafp_storage_to_xfer(p)  = leafp_storage(p) * bgtr(p)
               frootp_storage_to_xfer(p) = frootp_storage(p) * bgtr(p)
               if (woody(ivt(p)) == 1.0_r8) then
                  livestemp_storage_to_xfer(p)  = livestemp_storage(p) * bgtr(p)
                  deadstemp_storage_to_xfer(p)  = deadstemp_storage(p) * bgtr(p)
                  livecrootp_storage_to_xfer(p) = livecrootp_storage(p) * bgtr(p)
                  deadcrootp_storage_to_xfer(p) = deadcrootp_storage(p) * bgtr(p)
               end if
            end if

         end if ! end if stress deciduous

      end do ! end of pft loop

    end associate

  end subroutine CNStressDecidPhenology

  !-----------------------------------------------------------------------
  subroutine CropPhenology(num_pcropp, filter_pcropp                     , &
       waterstate_vars, temperature_vars, crop_vars, canopystate_vars, cnstate_vars , &
       carbonstate_vars, nitrogenstate_vars,carbonflux_vars,nitrogenflux_vars,&
       phosphorusstate_vars, phosphorusflux_vars)

    ! !DESCRIPTION:
    ! Code from AgroIBIS to determine crop phenology and code from CN to
    ! handle CN fluxes during the phenological onset                       & offset periods.

    ! !USES:
    use clm_time_manager , only : get_curr_date, get_curr_calday, get_days_per_year
    use pftvarcon        , only : ncorn, nscereal, nwcereal, nsoybean, gddmin, hybgdd
    use pftvarcon        , only : nwcerealirrig, nsoybeanirrig, ncornirrig, nscerealirrig
    use pftvarcon        , only : lfemerg, grnfill, mxmat, minplanttemp, planttemp
    use elm_varcon       , only : spval, secspday
    use CropType         , only : tcvp, tcvt, cst
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_pcropp       ! number of prog crop patches in filter
    integer                  , intent(in)    :: filter_pcropp                                    (:) ! filter for prognostic crop patches
    type(waterstate_type)    , intent(in)    :: waterstate_vars
    type(temperature_type)   , intent(in)    :: temperature_vars
    type(crop_type)          , intent(inout) :: crop_vars
    type(canopystate_type)   , intent(in)    :: canopystate_vars
    type(cnstate_type)       , intent(inout) :: cnstate_vars
    type(carbonstate_type)   , intent(inout) :: carbonstate_vars
    type(nitrogenstate_type) , intent(inout) :: nitrogenstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusstate_type) , intent(inout) :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars

    !
    ! LOCAL VARAIBLES:
    integer kyr       ! current year
    integer kmo       !         month of year  (1, ..., 12)
    integer kda       !         day of month   (1, ..., 31)
    integer mcsec     !         seconds of day (0, ..., seconds/day)
    integer jday      ! julian day of the year
    integer fp,p      ! patch indices
    integer c         ! column indices
    integer g         ! gridcell indices
    integer t         ! topographic indices
    integer h         ! hemisphere indices
    integer idpp      ! number of days past planting
    real(r8) dayspyr  ! days per year
    real(r8) crmcorn  ! comparitive relative maturity for corn
    real(r8) ndays_on ! number of days to fertilize
    logical p_season  ! precipitation seasonal
    logical t_season  ! temperature seasonal
    logical no_season ! neither temperature or precipitation seasonal
    real(r8), parameter :: minrain = 0.1    ! minimum rainfall for planting
    real(r8), parameter :: minwet = 0.2     ! minimum fraction of saturation for planting
    real(r8), parameter :: maxwet = 0.8     ! maximum fraction of saturation for planting
    !------------------------------------------------------------------------

    associate(                                                                 &
         ivt                =>    veg_pp%itype                               , & ! Input:  [integer  (:) ]  pft vegetation type

         leaf_long          =>    veg_vp%leaf_long                           , & ! Input:  [real(r8) (:) ]  leaf longevity (yrs)
         leafcn             =>    veg_vp%leafcn                              , & ! Input:  [real(r8) (:) ]  leaf C:N (gC/gN)
         fertnitro          =>    veg_vp%fertnitro                           , & ! Input:  [real(r8) (:) ]  max fertilizer to be applied in total (kgN/m2)

         t_ref2m_min        =>    veg_es%t_ref2m_min         , & ! Input:  [real(r8) (:) ]  daily minimum of average 2 m height surface air temperature (K)
         t10                =>    veg_es%t_a10               , & ! Input:  [real(r8) (:) ]  10-day running mean of the 2 m temperature (K)
         a5tmin             =>    veg_es%t_a5min             , & ! Input:  [real(r8) (:) ]  5-day running mean of min 2-m temperature
         a10tmin            =>    veg_es%t_a10min            , & ! Input:  [real(r8) (:) ]  10-day running mean of min 2-m temperature
         gdd020             =>    veg_es%gdd020              , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd0
         gdd820             =>    veg_es%gdd820              , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd8
         gdd1020            =>    veg_es%gdd1020             , & ! Input:  [real(r8) (:) ]  20 yr mean of gdd10
         hui                =>    crop_vars%gddplant_patch                   , & ! Input:  [real(r8) (:) ]  gdd since planting (gddplant)
         leafout            =>    crop_vars%gddtsoi_patch                    , & ! Input:  [real(r8) (:) ]  gdd from top soil layer temperature

         tlai               =>    canopystate_vars%tlai_patch                , & ! Input:  [real(r8) (:) ]  one-sided leaf area index, no burying by snow

         idop               =>    cnstate_vars%idop_patch                    , & ! Output: [integer  (:) ]  date of planting
         harvdate           =>    crop_vars%harvdate_patch                   , & ! Output: [integer  (:) ]  harvest date
         croplive           =>    crop_vars%croplive_patch                   , & ! Output: [logical  (:) ]  Flag, true if planted, not harvested
         cropplant          =>    crop_vars%cropplant_patch                  , & ! Output: [logical  (:) ]  Flag, true if crop may be planted
         gddmaturity        =>    cnstate_vars%gddmaturity_patch             , & ! Output: [real(r8) (:) ]  gdd needed to harvest
         huileaf            =>    cnstate_vars%huileaf_patch                 , & ! Output: [real(r8) (:) ]  heat unit index needed from planting to leaf emergence
         huigrain           =>    cnstate_vars%huigrain_patch                , & ! Output: [real(r8) (:) ]  same to reach vegetative maturity
         cumvd              =>    cnstate_vars%cumvd_patch                   , & ! Output: [real(r8) (:) ]  cumulative vernalization d?ependence?
         hdidx              =>    cnstate_vars%hdidx_patch                   , & ! Output: [real(r8) (:) ]  cold hardening index?
         vf                 =>    crop_vars%vf_patch                         , & ! Output: [real(r8) (:) ]  vernalization factor
         bglfr              =>    cnstate_vars%bglfr_patch                   , & ! Output: [real(r8) (:) ]  background litterfall rate (1/s)
         bgtr               =>    cnstate_vars%bgtr_patch                    , & ! Output: [real(r8) (:) ]  background transfer growth rate (1/s)
         lgsf               =>    cnstate_vars%lgsf_patch                    , & ! Output: [real(r8) (:) ]  long growing season factor [0-1]
         onset_flag         =>    cnstate_vars%onset_flag_patch              , & ! Output: [real(r8) (:) ]  onset flag
         offset_flag        =>    cnstate_vars%offset_flag_patch             , & ! Output: [real(r8) (:) ]  offset flag
         onset_counter      =>    cnstate_vars%onset_counter_patch           , & ! Output: [real(r8) (:) ]  onset counter
         offset_counter     =>    cnstate_vars%offset_counter_patch          , & ! Output: [real(r8) (:) ]  offset counter

         leafc_xfer         =>    veg_cs%leafc_xfer          , & ! Output: [real(r8) (:) ]  (gC/m2)   leaf C transfer

         crop_seedc_to_leaf =>    veg_cf%crop_seedc_to_leaf   , & ! Output: [real(r8) (:) ]  (gC/m2/s) seed source to PFT-level

         fert_counter       =>    veg_nf%fert_counter       , & ! Output: [real(r8) (:) ]  >0 fertilize; <=0 not (seconds)
         leafn_xfer         =>    veg_ns%leafn_xfer        , & ! Output: [real(r8) (:) ]  (gN/m2)   leaf N transfer
         crop_seedn_to_leaf =>    veg_nf%crop_seedn_to_leaf , & ! Output: [real(r8) (:) ]  (gN/m2/s) seed source to PFT-level
         crpyld             =>    crop_vars%crpyld_patch                     , & ! Output:  [real(r8) ):)]  harvested crop (bu/acre)
         dmyield            =>    crop_vars%dmyield_patch                    , & ! Output:  [real(r8) ):)]  dry matter harvested crop (t/ha)
         leafcp             =>    veg_vp%leafcp                              , & ! Input:  [real(r8) (:) ]  leaf C:P (gC/gP)
         leafp_xfer         =>    veg_ps%leafp_xfer      , & ! Output: [real(r8) (:) ]  (gP/m2)   leaf P transfer
         crop_seedp_to_leaf =>    veg_pf%crop_seedp_to_leaf , & ! Output: [real(r8) (:) ]  (gP/m2/s) seed source to PFT-level
         fert               =>    veg_nf%fert               , & ! Output: [real(r8) (:) ]  (gN/m2/s) fertilizer applied each timestep
         cvt                =>    crop_vars%cvt_patch                     , & ! Output:  [real(r8) ):)]  exp weighted moving average CV precip
         cvp                =>    crop_vars%cvp_patch                     , & ! Output:  [real(r8) ):)]  exp weighted moving average average CV temp
         xt_bar             =>    crop_vars%xt_bar_patch                  , & ! Output:  [real(r8) ):)]  exp weighted moving average average monthly temp
         plantmonth         =>    crop_vars%plantmonth_patch              , & ! Output:  [real(r8) ):)]  plant month
         plantday           =>    crop_vars%plantday_patch                , & ! Output:  [real(r8) ):)]  plant day
         harvday            =>    crop_vars%harvday_patch                 , & ! Ouptut:  [real(r8) ):)]  harvest day
         forc_rain          =>    top_af%rain                             , & ! Input:   [real(r8) (:)]  rainfall rate
         wf2                =>    col_ws%wf2                                & ! Output:  [real(r8) (:)]  soil water as frac. of whc for top 0.17 m
         )

      ! get time info
      dayspyr = get_days_per_year()
      jday    = get_curr_calday()
      call get_curr_date(kyr, kmo, kda, mcsec)

      ndays_on = 20._r8 ! number of days to fertilize

      do fp = 1, num_pcropp
         p = filter_pcropp(fp)
         c = veg_pp%column(p)
         g = veg_pp%gridcell(p)
         t = veg_pp%topounit(p)
         h = inhemi(p)

         ! background litterfall and transfer rates; long growing season factor

         bglfr(p) = 0._r8 ! this value changes later in a crop's life cycle
         bgtr(p)  = 0._r8
         lgsf(p)  = 0._r8

         ! B.Drewniak - zero our yield calculator
         crpyld(p)  = 0._r8
         dmyield(p) = 0._r8

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
               plantday(p)  = NOT_Planted

               ! keep next for continuous, annual winter temperate cereal type crop;
               ! if we removed elseif,
               ! winter cereal grown continuously would amount to a cereal/fallow
               ! rotation because cereal would only be planted every other year

            else if (croplive(p) .and. (ivt(p) == nwcereal .or. ivt(p) == nwcerealirrig)) then
               cropplant(p) = .false.
               !           else ! not possible to have croplive and ivt==cornORsoy? (slevis)
               ! keep next for precip based annual crop that may be planted in
               ! January in the SH to prevent from cereal/fallow rotation
            else if (croplive(p) .and. cvp(p) > 0.4_r8) then
               cropplant(p) = .false.
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

            if (ivt(p) == nwcereal .or. ivt(p) == nwcerealirrig) then

               ! add check to only plant winter cereal after other crops (soybean, maize)
               ! have been harvested

               ! *** remember order of planting is crucial - in terms of which crops you want
               ! to be grown in what order ***

               ! in this case, corn or soybeans are assumed to be planted before
               ! cereal would be in any particular year that both patches are allowed
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
                  plantday(p)    = jday
                  harvday(p)     = NOT_Harvested
                  harvdate(p)    = NOT_Harvested
                  gddmaturity(p) = hybgdd(ivt(p))
                  leafc_xfer(p)  = 1._r8 ! initial seed at planting to appear
                  leafn_xfer(p)  = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  crop_seedc_to_leaf(p) = leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = leafn_xfer(p)/dt

                  leafp_xfer(p)  = leafc_xfer(p) / leafcp(ivt(p)) ! with onset
                  crop_seedp_to_leaf(p) = leafp_xfer(p)/dt

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
                  plantday(p)    = jday
                  harvday(p)     = NOT_Harvested
                  harvdate(p)    = NOT_Harvested
                  gddmaturity(p) = hybgdd(ivt(p))
                  leafc_xfer(p)  = 1._r8 ! initial seed at planting to appear
                  leafn_xfer(p)  = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  crop_seedc_to_leaf(p) = leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = leafn_xfer(p)/dt

                  leafp_xfer(p)  = leafc_xfer(p) / leafcp(ivt(p)) ! with onset
                  crop_seedp_to_leaf(p) = leafp_xfer(p)/dt
               else
                  gddmaturity(p) = 0._r8
               end if

            else ! not winter cereal... slevis: added distinction between NH and SH
               ! slevis: The idea is that jday will equal idop sooner or later in the year
               !         while the gdd part is either true or false for the year.
               ! Assign seasonality - either temperature or precipitation
               t_season = .false.
               p_season = .false.
               no_season = .false.
               if (cvp(p) > tcvp) then        ! precipitation CV indicator
                  if (cvt(p) >= tcvt) then    ! Both temperature and precip seasonality
                     if (minval(xt_bar(p,:)) .lt. cst)then ! cold season exists (i.e., min temp below 10 deg C)
                        t_season = .true.
                     else                            ! no cold season
                        p_season = .true.
                     end if
                  else                     ! precipitation seasonality
                     p_season = .true.
                  end if
               else if (cvt(p) >= tcvt) then   ! temperature seasonality
                  t_season = .true.
               else                        ! no seasonality at all - tropics
                  no_season = .true.
               end if

               if (kmo .eq. plantmonth(p)               .and. &
                    gdd820(p) /= spval                  .and. &
                    gdd820(p) >= gddmin(ivt(p))) then
                  ! also require reasonable soil moisture to plant
                  ! between 0.2 (dry) and 0.8 (wet) saturation in top soil
                  ! layers
                  if ( (wf2(c) .ge. minwet .and. wf2(c) .le. maxwet)  .and. &
                     (t_season                                  .and. &
                     t10(p) /= spval .and. a10tmin(p) /= spval  .and. &
                     t10(p)     > planttemp(ivt(p))             .and. &
                     a10tmin(p) > minplanttemp(ivt(p)))               &
                                                                .or.  &
                     (p_season                                  .and. &
                     forc_rain(t) .gt. minrain/dt)                    & ! rain threshold to trigger
                                                                .or.  & ! the beginnig of rain season
                     (no_season)) then


                     ! impose limit on growing season length needed
                     ! for crop maturity - for cold weather constraints
                     croplive(p)  = .true.
                     cropplant(p) = .true.
                     idop(p)      = jday
                     harvday(p)   = NOT_Harvested
                     harvdate(p)  = NOT_Harvested
                     harvdate(p)  = NOT_Harvested

                     ! go a specified amount of time before/after
                     ! climatological date
                     if (ivt(p)==nsoybean .or. ivt(p) == nsoybeanirrig) gddmaturity(p)=min(gdd1020(p),hybgdd(ivt(p)))
                     if (ivt(p)==ncorn .or. ivt(p)==ncornirrig) then
                        gddmaturity(p)=max(950._r8, min(gdd820(p)*0.85_r8, hybgdd(ivt(p))))
                        gddmaturity(p)=max(950._r8, min(gddmaturity(p)+150._r8,1850._r8))
                     end if
                     if (ivt(p)==nscereal .or. ivt(p) == nscerealirrig) gddmaturity(p)=min(gdd020(p),hybgdd(ivt(p)))

                     leafc_xfer(p) = 1._r8 ! initial seed at planting to appear
                     leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                     crop_seedc_to_leaf(p) = leafc_xfer(p)/dt
                     crop_seedn_to_leaf(p) = leafn_xfer(p)/dt

                     leafp_xfer(p) = leafc_xfer(p) / leafcp(ivt(p)) ! with onset
                     crop_seedp_to_leaf(p) = leafp_xfer(p)/dt
                  end if
                     ! If hit the max planting julian day -- go ahead and plant
               else if (kmo .eq. (plantmonth(p) + 1) .and. gdd820(p) > 0._r8 .and. &
                       gdd820(p) /= spval ) then
                  croplive(p)  = .true.
                  cropplant(p) = .true.
                  idop(p)      = jday
                  plantday(p)  = jday
                  harvday(p)   = NOT_Harvested
                  harvdate(p)  = NOT_Harvested

                  if (ivt(p)==nsoybean .or. ivt(p) == nsoybeanirrig) gddmaturity(p)=min(gdd1020(p),hybgdd(ivt(p)))
                  if (ivt(p)==ncorn .or. ivt(p)==ncornirrig) gddmaturity(p)=max(950._r8, min(gdd820(p)*0.85_r8, hybgdd(ivt(p))))
                  if (ivt(p)==nscereal .or. ivt(p) == nscerealirrig) gddmaturity(p)=min(gdd020(p),hybgdd(ivt(p)))

                  leafc_xfer(p) = 1._r8 ! initial seed at planting to appear
                  leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p)) ! with onset
                  crop_seedc_to_leaf(p) = leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = leafn_xfer(p)/dt

                  leafp_xfer(p) = leafc_xfer(p) / leafcp(ivt(p)) ! with onset
                  crop_seedp_to_leaf(p) = leafp_xfer(p)/dt
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

            if (ivt(p) == ncorn .or. ivt(p)==ncornirrig) then
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

            if (t_ref2m_min(p) < 1.e30_r8 .and. vf(p) /= 1._r8 .and. (ivt(p) == nwcereal .or. ivt(p) == nwcerealirrig)) then
               call vernalization(p, &
                    canopystate_vars, temperature_vars, waterstate_vars, cnstate_vars, crop_vars)
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
                  fert_counter(p)  = ndays_on * secspday
                  fert(p) = fertnitro(ivt(p)) * 1000._r8 / fert_counter(p)
               else
                  ! this ensures no re-entry to onset of phase2
                  ! b/c onset_counter(p) = onset_counter(p) - dt
                  ! at every time step

                  onset_counter(p) = dt
               end if

               ! enter harvest for one time step:
               ! - transfer live biomass to litter and to crop yield
               ! - send xsmrpool to the atmosphere
               ! if onset and harvest needed to last longer than one timestep
               ! the onset_counter would change from dt and you'd need to make
               ! changes to the offset subroutine below

            else if (hui(p) >= gddmaturity(p) .or. idpp >= mxmat(ivt(p))) then
               if (harvdate(p) >= NOT_Harvested) harvdate(p) = jday
               if (harvday(p) >= NOT_Harvested) harvday(p) = jday
               croplive(p) = .false.     ! no re-entry in greater if-block
               if (tlai(p) > 0._r8) then ! plant had emerged before harvest
                  offset_flag(p) = 1._r8
                  offset_counter(p) = dt
               else                      ! plant never emerged from the ground
                  crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
                  crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
                  crop_seedp_to_leaf(p) = crop_seedp_to_leaf(p) - leafp_xfer(p)/dt
                  leafc_xfer(p) = 0._r8  ! revert planting transfers
                  leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
                  leafp_xfer(p) = leafc_xfer(p) / leafcp(ivt(p))
               end if

               ! enter phase 3 while previous criteria fail and next is true;
               ! in terms of order, phase 3 occurs before harvest, but when
               ! harvest *can* occur, we want it to have first priority.
               ! AgroIBIS uses a complex formula for lai decline.
               ! Use CN's simple formula at least as a place holder (slevis)

            else if (hui(p) >= huigrain(p)) then
               bglfr(p) = 1._r8/(leaf_long(ivt(p))*dayspyr*secspday)
            end if

            ! continue fertilizer application while in phase 2;
            ! assumes that onset of phase 2 took one time step only

            if (fert_counter(p) <= 0._r8) then
               fert(p) = 0._r8
            else ! continue same fert application every timestep
               fert_counter(p) = fert_counter(p) - dt
            end if

         else   ! crop not live
            ! next 2 lines conserve mass if leaf*_xfer > 0 due to interpinic
            crop_seedc_to_leaf(p) = crop_seedc_to_leaf(p) - leafc_xfer(p)/dt
            crop_seedn_to_leaf(p) = crop_seedn_to_leaf(p) - leafn_xfer(p)/dt
            crop_seedp_to_leaf(p) = crop_seedp_to_leaf(p) - leafp_xfer(p)/dt
            onset_counter(p) = 0._r8
            leafc_xfer(p) = 0._r8
            leafn_xfer(p) = leafc_xfer(p) / leafcn(ivt(p))
            leafp_xfer(p) = leafc_xfer(p) / leafcp(ivt(p))
         end if ! croplive

      end do ! prognostic crops loop

    end associate

  end subroutine CropPhenology

  !-----------------------------------------------------------------------
  subroutine CropPhenologyInit(bounds)
    !
    ! !DESCRIPTION:
    ! Initialization of CropPhenology. Must be called after time-manager is
    ! initialized, and after ecophyscon file is read in.
    !
    ! !USES:
    use pftvarcon       , only: npcropmin, npcropmax, mnNHplantdate
    use pftvarcon       , only: mnSHplantdate, mxNHplantdate
    use pftvarcon       , only: mxSHplantdate
    use clm_time_manager, only: get_calday
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds
    !
    ! LOCAL VARAIBLES:
    integer           :: p,g,n,i                     ! indices
    !------------------------------------------------------------------------

    allocate( inhemi(bounds%begp:bounds%endp) )

    allocate( minplantjday(0:numpft,inSH)) ! minimum planting julian day
    allocate( maxplantjday(0:numpft,inSH)) ! minimum planting julian day

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
    do p = bounds%begp, bounds%endp
       g = veg_pp%gridcell(p)
       ! Northern hemisphere
       if ( grc_pp%latdeg(g) > 0.0_r8 )then
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
  subroutine vernalization(p, &
       canopystate_vars, temperature_vars, waterstate_vars, cnstate_vars, &
       crop_vars)
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
    integer                , intent(in) :: p    ! PATCH index running over
    type(canopystate_type) , intent(in) :: canopystate_vars
    type(temperature_type) , intent(in) :: temperature_vars
    type(waterstate_type)  , intent(in) :: waterstate_vars
    type(cnstate_type)     , intent(inout) :: cnstate_vars
    type(crop_type)        , intent(inout) :: crop_vars
    !
    ! LOCAL VARAIBLES:
    real(r8) tcrown                     ! ?
    real(r8) vd, vd1, vd2               ! vernalization dependence
    real(r8) tkil                       ! Freeze kill threshold
    integer  c,g                        ! indices
    !------------------------------------------------------------------------

    associate(                                               &
         tlai        => canopystate_vars%tlai_patch        , & ! Input:  [real(r8) (:) ]  one-sided leaf area index, no burying by snow

         t_ref2m     => veg_es%t_ref2m     , & ! Input:  [real(r8) (:) ]  2 m height surface air temperature (K)
         t_ref2m_min => veg_es%t_ref2m_min , & ! Input:  [real(r8) (:) ] daily minimum of average 2 m height surface air temperature (K)
         t_ref2m_max => veg_es%t_ref2m_max , & ! Input:  [real(r8) (:) ] daily maximum of average 2 m height surface air temperature (K)

         snow_depth  => col_ws%snow_depth     , & ! Input:  [real(r8) (:) ]  snow height (m)

         hdidx       => cnstate_vars%hdidx_patch           , & ! Output: [real(r8) (:) ]  cold hardening index?
         cumvd       => cnstate_vars%cumvd_patch           , & ! Output: [real(r8) (:) ]  cumulative vernalization d?ependence?
         vf          => crop_vars%vf_patch              , & ! Output: [real(r8) (:) ]  vernalization factor for cereal
         gddmaturity => cnstate_vars%gddmaturity_patch     , & ! Output: [real(r8) (:) ]  gdd needed to harvest
         huigrain    => cnstate_vars%huigrain_patch          & ! Output: [real(r8) (:) ]  heat unit index needed to reach vegetative maturity
         )

      c = veg_pp%column(p)

      ! for all equations - temperatures must be in degrees (C)
      ! calculate temperature of crown of crop (e.g., 3 cm soil temperature)
      ! snow depth in centimeters

      if (t_ref2m(p) < tfrz) then !slevis: t_ref2m inst of td=daily avg (K)
         tcrown = 2._r8 + (t_ref2m(p) - tfrz) * (0.4_r8 + 0.0018_r8 * &
              (min(snow_depth(c)*100._r8, 15._r8) - 15._r8)**2)
      else !slevis: snow_depth inst of adsnod=daily average (m)
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

    end associate

  end subroutine vernalization

  !-----------------------------------------------------------------------
  subroutine CropPlantDate (num_soilp, filter_soilp, num_pcropp, filter_pcropp, &
        temperature_vars, cnstate_vars, crop_vars)
    !
    ! !DESCRIPTION:
    ! For determining the plant month for crops, plant day is established in
    ! CropPhenologyMod
    !
    ! !USES:
    use clm_time_manager , only : get_curr_date
    use clm_time_manager , only : get_step_size
    use elm_varcon       , only : secspday
    use elm_varpar       , only : numpft
    use pftvarcon        , only : planttemp
    use CropMod          , only : calculate_eto, plant_month

    !
    ! !ARGUMENTS:
    integer                , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                , intent(in)    :: filter_soilp(:) ! filter for soil patches
    integer                , intent(in)    :: num_pcropp      ! number of prognostic crops in filter
    integer                , intent(in)    :: filter_pcropp(:)! filter for prognostic crop patches
    type(temperature_type) , intent(in)    :: temperature_vars
    type(cnstate_type)     , intent(inout) :: cnstate_vars
    type(crop_type)        , intent(inout) :: crop_vars

    !
    ! !LOCAL VARIABLES:
    integer :: p,c,m,t,n,h          ! indices
    integer :: fp                   ! lake filter pft index
    integer kyr                     ! current year
    integer kmo                     ! month of year  (1, ..., 12)
    integer kda                     ! day of month   (1, ..., 31)
    integer mcsec                   ! seconds of day (0, ..., seconds/day)
    integer , parameter :: nmon = 12        ! number months in year
    real(r8), parameter :: alpha = 0.05_r8  ! coefficient representing degree of weighting decrease
    real(r8), parameter :: mon = 12._r8     ! used for some calculations (number of days in month)
    real(r8) :: mu_p, mu_t, sigmasum_t, sigmasum_p, sigma_t, sigma_p
    real(r8) :: es, ETout
    integer, dimension(12) :: ndaypm= &
         (/31,28,31,30,31,30,31,31,30,31,30,31/) !days per month
    !-----------------------------------------------------------------------

    associate(                                                &
         ivt            => veg_pp%itype                     , &
         t_ref2m        => veg_es%t_ref2m                   , & ! Input: [real(r8) (:) ]  2m air temperature (K)
         forc_rain      => top_af%rain                      , & ! Input: [real(r8) (:) ]  rainfall rate
         forc_snow      => top_af%snow                      , & ! Input: [real(r8) (:) ]  downscaled snow
         forc_rh        => top_as%rhbot                     , & ! Input: [real(r8) (:) ]  relative humidity
         forc_wind      => top_as%windbot                   , & ! Input: [real(r8) (:) ]  atmospheric wind speed (m/s)
         forc_pbot      => top_as%pbot                      , & ! Input: [real(r8) (:) ]  downscaled surface pressure (Pa)
         eflx_soil_grnd => veg_ef%eflx_soil_grnd            , & ! Input: [real(r8) (:) ]  soil heat flux (W/m**2) [+ = into soil]
         netrad         => veg_ef%netrad                    , & ! Input: [real(r8) (:) ]  net radiation (positive downward) (W/m**2)
         nyrs_crop_active => crop_vars%nyrs_crop_active_patch,   & ! InOut:  [integer (:)  ]  number of years this crop patch has been active
         cvt            => crop_vars%cvt_patch              , & ! Output: [real(r8) (:) ]  coefficient of variance temperature
         cvp            => crop_vars%cvp_patch              , & ! Output: [real(r8) (:) ]     "     "    "    "    precipitation
         xt             => crop_vars%xt_patch               , & ! Output: [real(r8) (:) ]  monthly average temperature
         xp             => crop_vars%xp_patch               , & ! Output: [real(r8) (:) ]  monthly total precipitation
         xt_bar         => crop_vars%xt_bar_patch           , & ! Output: [real(r8) (:) ]  exp weighted moving ave of temp
         xp_bar         => crop_vars%xp_bar_patch           , & ! Output: [real(r8) (:) ]  exp weighted moving ave of precip
         prev_xt_bar    => crop_vars%prev_xt_bar_patch      , & ! Output: [real(r8) (:) ]  previous years exp weighted moving ave of temp
         prev_xp_bar    => crop_vars%prev_xp_bar_patch      , & ! Output: [real(r8) (:) ]  previous years exp weighted moving ave of prec
         p2ETo          => crop_vars%p2ETo_patch            , & ! Output: [real(r8) (:) ]  precipitation:evapotranspiration ratio
         p2ETo_bar      => crop_vars%p2ETo_bar_patch        , & ! Output: [real(r8) (:) ]  exp weighted moving ave of P:PET
         prev_p2ETo_bar => crop_vars%prev_p2ETo_bar_patch   , & ! Output: [real(r8) (:) ]  previous years exp weighted moving ave of P:PET
         P2E_rm         => crop_vars%P2E_rm_patch           , & ! Output: [real(r8) (:) ]  precipitation:evapotranspiration ratio 4-month sum
         ETo            => crop_vars%ETo_patch              , & ! Output: [real(r8) (:) ]  reference evapotranspiration (mm)
         plantmonth     => crop_vars%plantmonth_patch         & ! Output: [real(r8) (:) ]  month to plant crops
         )

      !
      ! The following crop related steps are done here rather than CropPhenology
      ! so that they will be completed each time-step rather than with doalb.
      !
      ! The following lines of code will determine the precipitation or
      ! temperature seasonality of the grid cell for determing the plant trigger
      ! The method is from Waha et al., 2012 which calculates the coefficient of
      ! variation of temperature and precipitation to determin which climate
      ! variable drives the growing season. Then thresholds for the warm season
      ! and wet season are determined from either temperature thresholds or the
      ! P:PET ratio

      dt      = real( get_step_size(), r8 )
      fracday = dt/secspday

      if (num_pcropp > 0) then
         ! get time-related info
         call get_curr_date(kyr, kmo, kda, mcsec)
      end if

      do fp = 1,num_pcropp
         p = filter_pcropp(fp)
         c = veg_pp%column(p)
         t = veg_pp%topounit(p)
         h = inhemi(p)
         if (kmo == 1 .and. kda == 1 .and. nyrs_crop_active(p) == 0) then ! YR 1:
            P2E_rm(p,:) = 0._r8
            if (mcsec == 3600) then
               xt(p,:) = 0._r8          ! set these values to zero the first of the year
               xp(p,:) = 0._r8
               ETo(p,:) = 0._r8
            end if
         end if
         if (kmo == 1 .and. kda == 1 .and. mcsec == 0) then        ! <-- END of EVERY YR:
            if (nyrs_crop_active(p) > 0) then
               ! grab the precipitation to evapotranspiration from previous year
               ! last months (oct-dec) before it is overwritten
               ! for the 4-month sum, zero out the remainder
               P2E_rm(p,1) = sum(prev_p2ETo_bar(p,10:12))
               P2E_rm(p,2) = sum(prev_p2ETo_bar(p,11:12))
               P2E_rm(p,3) = prev_p2ETo_bar(p,12)
               P2E_rm(p,4:12) = 0._r8
               ! calculate the annual mean temperature and precipitation
               mu_t = 0._r8
               mu_p = 0._r8
               do m=1,nmon
                  mu_t = mu_t + (xt_bar(p,m))
                  mu_p = mu_p + (xp_bar(p,m))
               end do
               mu_t = mu_t/mon
               mu_p = mu_p/mon
               ! calculate the standard deviation of temperature and
               ! precipitation
               sigmasum_t = 0._r8
               sigmasum_p = 0._r8
               do m = 1,nmon
                  sigmasum_t = sigmasum_t + (xt_bar(p,m) - mu_t)**2
                  sigmasum_p = sigmasum_p + (xp_bar(p,m) - mu_p)**2
               end do
               sigma_t = sqrt( (1._r8/(mon-1._r8)) * sigmasum_t)
               sigma_p = sqrt( (1._r8/(mon-1._r8)) * sigmasum_p)
               ! calculate the new coefficient of variance of temperature and
               ! precipitation
               cvt(p) = sigma_t/mu_t
               cvp(p) = sigma_p/mu_p
               ! calculate the P:ET0 ratio for precipitation seasonality
               ! now finish the P:PET calculation for the 4-month sum
               P2E_rm(p,1) = P2E_rm(p,1) + p2ETo_bar(p,1)
               P2E_rm(p,2) = P2E_rm(p,2) + p2ETo_bar(p,1) + p2ETo_bar(p,2)
               P2E_rm(p,3) = P2E_rm(p,3) + p2ETo_bar(p,1) + p2ETo_bar(p,2) + p2ETo_bar(p,3)
               do m=4,nmon
                  n = m-3
                  P2E_rm(p,m) = P2E_rm(p,m) + sum(p2ETo_bar(p,n:m))
               end do
            else
               cvt(p) = 0._r8
               cvp(p) = 0._r8
            end if
            ! the next call will determine next years plant month
            call plant_month(p, cvt(p), cvp(p), xt_bar(p,:), P2E_rm(p,:), minplantjday(ivt(p),h), plantmonth(p))
            prev_xt_bar(p,:) = xt_bar(p,:)   ! save last years values
            prev_xp_bar(p,:) = xp_bar(p,:)
            prev_p2ETo_bar(p,:) = p2ETo_bar(p,:)
            xt(p,:) = 0._r8          ! reset monthly averages
            xp(p,:) = 0._r8
            ETo(p,:) = 0._r8
         end if

         xt(p,kmo) = xt(p,kmo) + t_ref2m(p) * fracday/ndaypm(kmo) ! monthly average temperature
         xp(p,kmo) = xp(p,kmo) + (forc_rain(t)+forc_snow(t))*dt   ! monthly average precipitation
         ! calculate the potential evapotranspiration
         call calculate_eto(t_ref2m(p), netrad(p), eflx_soil_grnd(p), forc_pbot(t), forc_rh(t), forc_wind(t), es, dt, ETout)
         ! monthly ETo
         ETo(p,kmo) = ETo(p,kmo) + ETout
         ! calculate the P:PET for each month
         p2ETo(p,kmo) = xp(p,kmo)/ETo(p,kmo)

         if (nyrs_crop_active(p) == 0) then ! for the first year, use last years values
            prev_xt_bar(p,kmo) = xt(p,kmo)
            prev_xp_bar(p,kmo) = xp(p,kmo)
            prev_p2ETo_bar(p,kmo) = p2ETo(p,kmo)
         end if

         ! the following calculates the exponential weighted moving average of
         ! temperature, precipitation and P:PET for each month. This still  gives an
         ! average, but with more weight toward recent years to account for farmer
         ! memory, rather than a simple 20-year average with is used for growing
         ! season calculations (i.e., GDD). It is used in the plant_month call
         xt_bar(p,kmo) = alpha * xt(p,kmo) + (1-alpha) * prev_xt_bar(p,kmo)
         xp_bar(p,kmo) = alpha * xp(p,kmo) + (1-alpha) * prev_xp_bar(p,kmo)
         p2ETo_bar(p,kmo) = alpha * p2ETo(p,kmo) + (1-alpha) * prev_p2ETo_bar(p,kmo)

      end do

    end associate

  end subroutine CropPlantDate

  !-----------------------------------------------------------------------
  subroutine CNOnsetGrowth (num_soilp, filter_soilp, &
       cnstate_vars, &
       carbonstate_vars, nitrogenstate_vars, carbonflux_vars,nitrogenflux_vars,&
       phosphorusstate_vars,phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! Determines the flux of stored C and N from transfer pools to display
    ! pools during the phenological onset period.
    ! add flux for phosphorus - X.YANG
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)       , intent(in)    :: cnstate_vars
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusstate_type) , intent(in)    :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: fp           ! lake filter pft index
    real(r8):: t1           ! temporary variable
    !-----------------------------------------------------------------------

    associate(                                                                                             &
         ivt                                 =>    veg_pp%itype                                             , & ! Input:  [integer   (:) ]  pft vegetation type

         woody                               =>    veg_vp%woody                                      , & ! Input:  [real(r8)  (:) ]  binary flag for woody lifeform (1=woody, 0=not woody)

         onset_flag                          =>    cnstate_vars%onset_flag_patch                           , & ! Input:  [real(r8)  (:) ]  onset flag
         onset_counter                       =>    cnstate_vars%onset_counter_patch                        , & ! Input:  [real(r8)  (:) ]  onset days counter
         bgtr                                =>    cnstate_vars%bgtr_patch                                 , & ! Input:  [real(r8)  (:) ]  background transfer growth rate (1/s)

         leafc_xfer                          =>    veg_cs%leafc_xfer                     , & ! Input:  [real(r8)  (:) ]  (gC/m2) leaf C transfer
         frootc_xfer                         =>    veg_cs%frootc_xfer                    , & ! Input:  [real(r8)  (:) ]  (gC/m2) fine root C transfer
         livestemc_xfer                      =>    veg_cs%livestemc_xfer                 , & ! Input:  [real(r8)  (:) ]  (gC/m2) live stem C transfer
         deadstemc_xfer                      =>    veg_cs%deadstemc_xfer                 , & ! Input:  [real(r8)  (:) ]  (gC/m2) dead stem C transfer
         livecrootc_xfer                     =>    veg_cs%livecrootc_xfer                , & ! Input:  [real(r8)  (:) ]  (gC/m2) live coarse root C transfer
         deadcrootc_xfer                     =>    veg_cs%deadcrootc_xfer                , & ! Input:  [real(r8)  (:) ]  (gC/m2) dead coarse root C transfer

         leafn_xfer                          =>    veg_ns%leafn_xfer                   , & ! Input:  [real(r8)  (:) ]  (gN/m2) leaf N transfer
         frootn_xfer                         =>    veg_ns%frootn_xfer                  , & ! Input:  [real(r8)  (:) ]  (gN/m2) fine root N transfer
         livestemn_xfer                      =>    veg_ns%livestemn_xfer               , & ! Input:  [real(r8)  (:) ]  (gN/m2) live stem N transfer
         deadstemn_xfer                      =>    veg_ns%deadstemn_xfer               , & ! Input:  [real(r8)  (:) ]  (gN/m2) dead stem N transfer
         livecrootn_xfer                     =>    veg_ns%livecrootn_xfer              , & ! Input:  [real(r8)  (:) ]  (gN/m2) live coarse root N transfer
         deadcrootn_xfer                     =>    veg_ns%deadcrootn_xfer              , & ! Input:  [real(r8)  (:) ]  (gN/m2) dead coarse root N transfer


         leafp_xfer                          =>    veg_ps%leafp_xfer                   , & ! Input:  [real(r8)(:)   ]  (gP/m2) leaf P transfer
         frootp_xfer                         =>    veg_ps%frootp_xfer                  , & ! Input:  [real(r8)(:)   ]  (gP/m2) fine root P transfer
         livestemp_xfer                      =>    veg_ps%livestemp_xfer               , & ! Input:  [real(r8)(:)   ]  (gP/m2) live stem P transfer
         deadstemp_xfer                      =>    veg_ps%deadstemp_xfer               , & ! Input:  [real(r8)(:)   ]  (gP/m2) dead stem P transfer
         livecrootp_xfer                     =>    veg_ps%livecrootp_xfer              , & ! Input:  [real(r8)(:)   ]  (gP/m2) live coarse root P transfer
         deadcrootp_xfer                     =>    veg_ps%deadcrootp_xfer              , & ! Input:  [real(r8)(:)   ]  (gP/m2) dead coarse root P transfer


         leafc_xfer_to_leafc                 =>    veg_cf%leafc_xfer_to_leafc             , & ! Output:  [real(r8) (:) ]
         frootc_xfer_to_frootc               =>    veg_cf%frootc_xfer_to_frootc           , & ! Output:  [real(r8) (:) ]
         livestemc_xfer_to_livestemc         =>    veg_cf%livestemc_xfer_to_livestemc     , & ! Output:  [real(r8) (:) ]
         deadstemc_xfer_to_deadstemc         =>    veg_cf%deadstemc_xfer_to_deadstemc     , & ! Output:  [real(r8) (:) ]
         livecrootc_xfer_to_livecrootc       =>    veg_cf%livecrootc_xfer_to_livecrootc   , & ! Output:  [real(r8) (:) ]
         deadcrootc_xfer_to_deadcrootc       =>    veg_cf%deadcrootc_xfer_to_deadcrootc   , & ! Output:  [real(r8) (:) ]

         leafn_xfer_to_leafn                 =>    veg_nf%leafn_xfer_to_leafn           , & ! Output:  [real(r8) (:) ]
         frootn_xfer_to_frootn               =>    veg_nf%frootn_xfer_to_frootn         , & ! Output:  [real(r8) (:) ]
         livestemn_xfer_to_livestemn         =>    veg_nf%livestemn_xfer_to_livestemn   , & ! Output:  [real(r8) (:) ]
         deadstemn_xfer_to_deadstemn         =>    veg_nf%deadstemn_xfer_to_deadstemn   , & ! Output:  [real(r8) (:) ]
         livecrootn_xfer_to_livecrootn       =>    veg_nf%livecrootn_xfer_to_livecrootn , & ! Output:  [real(r8) (:) ]
         deadcrootn_xfer_to_deadcrootn       =>    veg_nf%deadcrootn_xfer_to_deadcrootn , & ! Output:  [real(r8) (:) ]

         leafp_xfer_to_leafp                 =>    veg_pf%leafp_xfer_to_leafp           , & ! Output:  [real(r8) (:) ]
         frootp_xfer_to_frootp               =>    veg_pf%frootp_xfer_to_frootp         , & ! Output:  [real(r8) (:) ]
         livestemp_xfer_to_livestemp         =>    veg_pf%livestemp_xfer_to_livestemp   , & ! Output:  [real(r8) (:) ]
         deadstemp_xfer_to_deadstemp         =>    veg_pf%deadstemp_xfer_to_deadstemp   , & ! Output:  [real(r8) (:) ]
         livecrootp_xfer_to_livecrootp       =>    veg_pf%livecrootp_xfer_to_livecrootp , & ! Output:  [real(r8) (:) ]
         deadcrootp_xfer_to_deadcrootp       =>    veg_pf%deadcrootp_xfer_to_deadcrootp   & ! Output:  [real(r8) (:) ]
         )

      ! patch loop
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
            leafp_xfer_to_leafp(p)   = t1 * leafp_xfer(p)
            frootp_xfer_to_frootp(p) = t1 * frootp_xfer(p)
            if (woody(ivt(p)) == 1.0_r8) then
               livestemc_xfer_to_livestemc(p)   = t1 * livestemc_xfer(p)
               deadstemc_xfer_to_deadstemc(p)   = t1 * deadstemc_xfer(p)
               livecrootc_xfer_to_livecrootc(p) = t1 * livecrootc_xfer(p)
               deadcrootc_xfer_to_deadcrootc(p) = t1 * deadcrootc_xfer(p)
               livestemn_xfer_to_livestemn(p)   = t1 * livestemn_xfer(p)
               deadstemn_xfer_to_deadstemn(p)   = t1 * deadstemn_xfer(p)
               livecrootn_xfer_to_livecrootn(p) = t1 * livecrootn_xfer(p)
               deadcrootn_xfer_to_deadcrootn(p) = t1 * deadcrootn_xfer(p)

               livestemp_xfer_to_livestemp(p)   = t1 * livestemp_xfer(p)
               deadstemp_xfer_to_deadstemp(p)   = t1 * deadstemp_xfer(p)
               livecrootp_xfer_to_livecrootp(p) = t1 * livecrootp_xfer(p)
               deadcrootp_xfer_to_deadcrootp(p) = t1 * deadcrootp_xfer(p)
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
            leafp_xfer_to_leafp(p)   = leafp_xfer(p) / dt
            frootp_xfer_to_frootp(p) = frootp_xfer(p) / dt
            if (woody(ivt(p)) == 1.0_r8) then
               livestemc_xfer_to_livestemc(p)   = livestemc_xfer(p) / dt
               deadstemc_xfer_to_deadstemc(p)   = deadstemc_xfer(p) / dt
               livecrootc_xfer_to_livecrootc(p) = livecrootc_xfer(p) / dt
               deadcrootc_xfer_to_deadcrootc(p) = deadcrootc_xfer(p) / dt
               livestemn_xfer_to_livestemn(p)   = livestemn_xfer(p) / dt
               deadstemn_xfer_to_deadstemn(p)   = deadstemn_xfer(p) / dt
               livecrootn_xfer_to_livecrootn(p) = livecrootn_xfer(p) / dt
               deadcrootn_xfer_to_deadcrootn(p) = deadcrootn_xfer(p) / dt
               livestemp_xfer_to_livestemp(p)   = livestemp_xfer(p) / dt
               deadstemp_xfer_to_deadstemp(p)   = deadstemp_xfer(p) / dt
               livecrootp_xfer_to_livecrootp(p) = livecrootp_xfer(p) / dt
               deadcrootp_xfer_to_deadcrootp(p) = deadcrootp_xfer(p) / dt
            end if
         end if ! end if bgtr

      end do ! end pft loop

    end associate

  end subroutine CNOnsetGrowth

 !----------------------------------------------------------------------
 subroutine CNCropHarvest (num_pcropp, filter_pcropp, num_soilc, filter_soilc, crop_vars, &
            cnstate_vars, carbonstate_vars, carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, &
            phosphorusstate_vars, phosphorusflux_vars)
   !
   ! !DESCRIPTION:
   ! This routine handles harvest for agriculture vegetation types, such as
   ! corn, soybean, and wheat. This routine allows harvest to be calculated
   ! instead of in the OffsetLitterfall subroutine. The harvest index is
   ! determined based on the LPJ model.
   !
   ! !ARGUMENTS:
   integer, intent(in) :: num_pcropp       ! number of prog crop pfts in filter
   integer, intent(in) :: filter_pcropp(:) ! filter for prognostic crop pfts
   integer, intent(in) :: num_soilc        ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:)  ! soil column filter

    type(crop_type)         , intent(inout) :: crop_vars
    type(cnstate_type)      , intent(inout) :: cnstate_vars
    type(carbonstate_type)  , intent(in)    :: carbonstate_vars
    type(carbonflux_type)   , intent(inout) :: carbonflux_vars
    type(nitrogenstate_type), intent(in)    :: nitrogenstate_vars
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars
    type(phosphorusstate_type),intent(inout):: phosphorusstate_vars
    type(phosphorusflux_type), intent(inout):: phosphorusflux_vars
   !
   ! !LOCAL VARIABLES:
   ! local pointers to implicit in scalars
   integer :: p                             ! indices
   integer :: fp                            ! lake filter pft index
   real(r8):: t1                            ! temporary variable
   real(r8):: cgrain                        ! amount of carbon in the grain
   !-------------------------------------------------------------------------
   associate(&
   ivt                   =>    veg_pp%itype                                   , & ! Input:  [integer (:)]  pft vegetation type
   offset_flag           =>    cnstate_vars%offset_flag_patch              , & ! Input:  [real(r8) (:) ]  offset flag
   offset_counter        =>    cnstate_vars%offset_counter_patch           , & ! Input:  [real(r8) (:) ]  offset days counter

   presharv              =>    veg_vp%presharv                         , & ! Input:  [real(r8) (:) ]  porportion of residue harvested
   fyield                =>    veg_vp%fyield                           , & ! Input:  [real(r8) (:) ]  fraction of grain actually harvested
   convfact              =>    veg_vp%convfact                         , & ! Input:  [real(r8) (:) ]  converstion factor for bu/acre

   leafc                 =>    veg_cs%leafc                , & ! Input:  [real(r8) (:) ]  (gC/m2) leaf C
   grainc                =>    veg_cs%grainc               , & ! Input:  [real(r8) (:) ]  (gC/m2) grain C
   livestemc             =>    veg_cs%livestemc            , & ! Input:  [real(r8) (:) ]  (gC/m2) livestem C
   leafn                 =>    veg_ns%leafn              , & ! Input:  [real(r8) (:) ]  (gN/m2) leaf N
   grainn                =>    veg_ns%grainn             , & ! Input:  [real(r8) (:) ]  (gN/m2) grain N
   livestemn             =>    veg_ns%livestemn          , & ! Input:  [real(r8) (:) ]  (gN/m2) livestem N
   leafp                 =>    veg_ps%leafp            , & ! Input:  [real(r8) (:) ]  (gP/m2) leaf P
   grainp                =>    veg_ps%grainp           , & ! Input:  [real(r8) (:) ]  (gP/m2) grain P
   livestemp             =>    veg_ps%livestemp        , & ! Input:  [real(r8) (:) ]  (gP/m2) livestem P
   cpool_to_grainc       =>    veg_cf%cpool_to_grainc       , & ! Input:  [real(r8) (:) ]  allocation to grain C (gC/m2/s)
   cpool_to_livestemc    =>    veg_cf%cpool_to_livestemc    , & ! Input:  [real(r8) (:) ]  allocation to live stem C (gC/m2/s)
   cpool_to_leafc        =>    veg_cf%cpool_to_leafc        , & ! Input:  [real(r8) (:) ]  allocation to leaf C (gC/m2/s)
   npool_to_leafn        =>    veg_nf%npool_to_leafn      , & ! Input:  [real(r8) (:)]  allocation to grain N (gN/m2/s)
   npool_to_livestemn    =>    veg_nf%npool_to_livestemn  , & ! Input:  [real(r8) (:)]  allocation to grain N (gN/m2/s)
   npool_to_grainn       =>    veg_nf%npool_to_grainn     , & ! Input:  [real(r8) (:)]  allocation to grain N (gN/m2/s)
   ppool_to_leafp        =>    veg_pf%ppool_to_leafp    , & ! Input:  [real(r8) (:)]  allocation to grain P (gP/m2/s)
   ppool_to_livestemp    =>    veg_pf%ppool_to_livestemp, & ! Input:  [real(r8) (:)]  allocation to grain P (gP/m2/s)
   ppool_to_grainp       =>    veg_pf%ppool_to_grainp   , & ! Input:  [real(r8) (:)]  allocation to grain P (gP/m2/s)
   hrv_leafc_to_prod1c   =>    veg_cf%hrv_leafc_to_prod1c   , & ! Input:  [real(r8) (:)] crop leafc harvested
   hrv_livestemc_to_prod1c  => veg_cf%hrv_livestemc_to_prod1c, & ! Input:  [real(r8) (:)] crop stemc harvested
   hrv_grainc_to_prod1c  =>    veg_cf%hrv_grainc_to_prod1c  , & ! Input:  [real(r8) (:)] crop grainc harvested
   hrv_leafn_to_prod1n   =>    veg_nf%hrv_leafn_to_prod1n , & ! Input:  [real(r8) (:)] crop leafn harvested
   hrv_livestemn_to_prod1n  => veg_nf%hrv_livestemn_to_prod1n, & ! Input:  [real(r8) (:)] crop stemn harvested
   hrv_grainn_to_prod1n  =>    veg_nf%hrv_grainn_to_prod1n, & ! Input:  [real(r8) (:)] crop grainn harvested
   hrv_leafp_to_prod1p   =>    veg_pf%hrv_leafp_to_prod1p , & ! Input:  [real(r8) (:)] crop leafp harvested
   hrv_livestemp_to_prod1p  => veg_pf%hrv_livestemp_to_prod1p, & ! Input:  [real(r8) (:)] crop stemp harvested
   hrv_grainp_to_prod1p  =>    veg_pf%hrv_grainp_to_prod1p, & ! Input:  [real(r8) (:)] crop grainp harvested
   crpyld                =>    crop_vars%crpyld_patch                      , & ! InOut:  [real(r8) ):)]  harvested crop (bu/acre)
   dmyield               =>    crop_vars%dmyield_patch                       & ! InOut:  [real(r8) ):)]  dry matter harvested crop (t/ha)
   )

   cgrain = 0.50_r8
   do fp = 1,num_pcropp
      p = filter_pcropp(fp)
      ! only calculate during the offset period
      if (offset_flag(p) == 1._r8) then

         if (offset_counter(p) == dt) then
         t1 = 1._r8 / dt
              !calculate yield (crpyld = bu/acre and dmyield = t/ha)
              crpyld(p)    = (grainc(p)+cpool_to_grainc(p)*dt) * fyield(ivt(p)) * convfact(ivt(p)) / (cgrain * 1000)
              dmyield(p)   = (grainc(p)+cpool_to_grainc(p)*dt) * fyield(ivt(p)) * 0.01 / cgrain

              !calculate harvested carbon and nitrogen; remaining goes into litterpool
              !except for grain which goes into next years availc for growth after
              !planting
              hrv_leafc_to_prod1c(p)  = presharv(ivt(p)) * ((t1 * leafc(p)) + cpool_to_leafc(p))
              hrv_livestemc_to_prod1c(p)  =  presharv(ivt(p)) * ((t1 * livestemc(p)) + cpool_to_livestemc(p))
              hrv_grainc_to_prod1c(p) = t1 * grainc(p) + cpool_to_grainc(p)

              ! Do the same for Nitrogen
              hrv_leafn_to_prod1n(p) = presharv(ivt(p)) * ((t1 * leafn(p)) + npool_to_leafn(p))
              hrv_livestemn_to_prod1n(p) = presharv(ivt(p)) * ((t1 * livestemn(p)) + npool_to_livestemn(p))
              hrv_grainn_to_prod1n(p) = t1 * grainn(p) + npool_to_grainn(p)

              ! Do the same for Phosphorus
              hrv_leafp_to_prod1p(p) = presharv(ivt(p)) * ((t1 * leafp(p)) + ppool_to_leafp(p))
              hrv_livestemp_to_prod1p(p) = presharv(ivt(p)) * ((t1 * livestemp(p)) + ppool_to_livestemp(p))
              hrv_grainp_to_prod1p(p) = t1 * grainp(p) + ppool_to_grainp(p)

         end if ! offseddt_counter

      end if ! offset_flag
   end do

   ! gather all pft-level fluxes from harvest to the column
   ! for C and N inputs

   call CNCropHarvestPftToColumn(num_soilc, filter_soilc,cnstate_vars, &
                   carbonflux_vars, nitrogenflux_vars, phosphorusflux_vars)
    end associate
 end subroutine CNCropHarvest

  !-----------------------------------------------------------------------
  subroutine CNOffsetLitterfall (num_soilp, filter_soilp, &
       cnstate_vars, carbonstate_vars, carbonflux_vars, nitrogenflux_vars,&
       phosphorusflux_vars, nitrogenstate_vars,phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! Determines the flux of C and N from displayed pools to litter
    ! pools during the phenological offset period.
    !
    ! !USES:
    use pftvarcon , only : npcropmin
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                 , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)      , intent(inout) :: cnstate_vars
    type(carbonstate_type)  , intent(in)    :: carbonstate_vars
    type(carbonflux_type)   , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type) , intent(inout) :: phosphorusflux_vars
    type(nitrogenstate_type)  , intent(in)   :: nitrogenstate_vars
    type(phosphorusstate_type), intent(in)   :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p, c         ! indices
    integer :: fp           ! lake filter pft index
    real(r8):: t1           ! temporary variable
    !-----------------------------------------------------------------------

    associate(                                                                     &
         ivt                   =>    veg_pp%itype                                   , & ! Input:  [integer  (:) ]  pft vegetation type

         leafcn                =>    veg_vp%leafcn                           , & ! Input:  [real(r8) (:) ]  leaf C:N (gC/gN)
         lflitcn               =>    veg_vp%lflitcn                          , & ! Input:  [real(r8) (:) ]  leaf litter C:N (gC/gN)
         frootcn               =>    veg_vp%frootcn                          , & ! Input:  [real(r8) (:) ]  fine root C:N (gC/gN)
         livewdcn              =>    veg_vp%livewdcn                         , & ! Input:  [real(r8) (:) ]  live wood C:N (gC/gN)
         graincn               =>    veg_vp%graincn                          , & ! Input:  [real(r8) (:) ]  grain C:N (gC/gN)
         presharv              =>    veg_vp%presharv                         , & ! Input:  [real(r8) (:) ]  porportion of residue harvested

         leafcp                =>    veg_vp%leafcp                           , & ! Input:  [real(r8) (:) ]  leaf C:P (gC/gP)
         lflitcp               =>    veg_vp%lflitcp                          , & ! Input:  [real(r8) (:) ]  leaf litter C:P (gC/gP)
         frootcp               =>    veg_vp%frootcp                          , & ! Input:  [real(r8) (:) ]  fine root C:P (gC/gP)
         livewdcp              =>    veg_vp%livewdcp                         , & ! Input:  [real(r8) (:) ]  live wood C:P (gC/gP)
         graincp               =>    veg_vp%graincp                          , & ! Input:  [real(r8) (:) ]  grain C:P (gC/gP)

         offset_flag           =>    cnstate_vars%offset_flag_patch              , & ! Input:  [real(r8) (:) ]  offset flag
         offset_counter        =>    cnstate_vars%offset_counter_patch           , & ! Input:  [real(r8) (:) ]  offset days counter

         leafc                 =>    veg_cs%leafc                , & ! Input:  [real(r8) (:) ]  (gC/m2) leaf C
         frootc                =>    veg_cs%frootc               , & ! Input:  [real(r8) (:) ]  (gC/m2) fine root C
         grainc                =>    veg_cs%grainc               , & ! Input:  [real(r8) (:) ]  (gC/m2) grain C
         livestemc             =>    veg_cs%livestemc            , & ! Input:  [real(r8) (:) ]  (gC/m2) livestem C

         cpool_to_grainc       =>    veg_cf%cpool_to_grainc       , & ! Input:  [real(r8) (:) ]  allocation to grain C (gC/m2/s)
         cpool_to_livestemc    =>    veg_cf%cpool_to_livestemc    , & ! Input:  [real(r8) (:) ]  allocation to live stem C (gC/m2/s)
         cpool_to_leafc        =>    veg_cf%cpool_to_leafc        , & ! Input:  [real(r8) (:) ]  allocation to leaf C (gC/m2/s)
         cpool_to_frootc       =>    veg_cf%cpool_to_frootc       , & ! Input:  [real(r8) (:) ]  allocation to fine root C (gC/m2/s)
         prev_leafc_to_litter  =>    veg_cf%prev_leafc_to_litter  , & ! Output: [real(r8) (:) ]  previous timestep leaf C litterfall flux (gC/m2/s)
         prev_frootc_to_litter =>    veg_cf%prev_frootc_to_litter , & ! Output: [real(r8) (:) ]  previous timestep froot C litterfall flux (gC/m2/s)
         leafc_to_litter       =>    veg_cf%leafc_to_litter       , & ! Output: [real(r8) (:) ]  leaf C litterfall (gC/m2/s)
         frootc_to_litter      =>    veg_cf%frootc_to_litter      , & ! Output: [real(r8) (:) ]  fine root C litterfall (gC/m2/s)
         livestemc_to_litter   =>    veg_cf%livestemc_to_litter   , & ! Output: [real(r8) (:) ]  live stem C litterfall (gC/m2/s)
         grainc_to_food        =>    veg_cf%grainc_to_food        , & ! Output: [real(r8) (:) ]  grain C to food (gC/m2/s)

         livestemn_to_litter   =>    veg_nf%livestemn_to_litter , & ! Output: [real(r8) (:) ]  livestem N to litter (gN/m2/s)
         grainn_to_food        =>    veg_nf%grainn_to_food      , & ! Output: [real(r8) (:) ]  grain N to food (gN/m2/s)
         leafn_to_litter       =>    veg_nf%leafn_to_litter     , & ! Output: [real(r8) (:) ]  leaf N litterfall (gN/m2/s)
         leafn_to_retransn     =>    veg_nf%leafn_to_retransn   , & ! Output: [real(r8) (:) ]  leaf N to retranslocated N pool (gN/m2/s)
         frootn_to_litter      =>    veg_nf%frootn_to_litter    , & ! Output: [real(r8) (:) ]  fine root N litterfall (gN/m2/s)

         livestemp_to_litter   =>    veg_pf%livestemp_to_litter , & ! Output: [real(r8) (:) ]  livestem P to litter (gP/m2/s)
         grainp_to_food        =>    veg_pf%grainp_to_food      , & ! Output: [real(r8) (:) ]  grain P to food (gP/m2/s)
         leafp_to_litter       =>    veg_pf%leafp_to_litter     , & ! Output: [real(r8) (:) ]  leaf P litterfall (gP/m2/s)
         leafp_to_retransp     =>    veg_pf%leafp_to_retransp   , & ! Output: [real(r8) (:) ]  leaf P to retranslocated P pool (gP/m2/s)
         frootp_to_litter      =>    veg_pf%frootp_to_litter    , & ! Output: [real(r8) (:) ]  fine root P litterfall (gP/m2/s)

         prev_leafn_to_litter  =>    veg_nf%prev_leafn_to_litter    , & ! Output: [real(r8) (:) ]  previous timestep leaf N litterfall flux (gN/m2/s)
         prev_frootn_to_litter =>    veg_nf%prev_frootn_to_litter   , & ! Output: [real(r8) (:) ]  previous timestep froot N litterfall flux (gN/m2/s)
         prev_leafp_to_litter  =>    veg_pf%prev_leafp_to_litter  , & ! Output: [real(r8) (:) ]  previous timestep leaf P litterfall flux (gP/m2/s)
         prev_frootp_to_litter =>    veg_pf%prev_frootp_to_litter , & ! Output: [real(r8) (:) ]  previous timestep froot P litterfall flux (gP/m2/s)
         leafn                 =>    veg_ns%leafn                  , & ! Input:  [real(r8) (:) ]  (gN/m2) leaf N
         frootn                =>    veg_ns%frootn                 , & ! Input:  [real(r8) (:) ]  (gN/m2) fine root N
         livestemn             =>    veg_ns%livestemn              , & ! Input:  [real(r8) (:) ]  (gN/m2) livestem N
         leafp                 =>    veg_ps%leafp                , & ! Input:  [real(r8) (:) ]  (gP/m2) leaf P
         frootp                =>    veg_ps%frootp               , & ! Input:  [real(r8) (:) ]  (gP/m2) fine root P
         livestemp             =>    veg_ps%livestemp            , & ! Input:  [real(r8) (:) ]  (gP/m2) livestem P
         npool_to_leafn        =>    veg_nf%npool_to_leafn          , &
         npool_to_frootn       =>    veg_nf%npool_to_frootn         , &
         npool_to_livestemn    =>    veg_nf%npool_to_livestemn      , &
         ppool_to_leafp        =>    veg_pf%ppool_to_leafp        , &
         ppool_to_frootp       =>    veg_pf%ppool_to_frootp       , &
         ppool_to_livestemp    =>    veg_pf%ppool_to_livestemp    , &
         hrv_leafc_to_prod1c   =>    veg_cf%hrv_leafc_to_prod1c       , & ! Input:  [real(r8) (:)] crop leafc harvested
         hrv_livestemc_to_prod1c  => veg_cf%hrv_livestemc_to_prod1c   , & ! Input:  [real(r8) (:)] crop stemc harvested
         hrv_leafn_to_prod1n   =>    veg_nf%hrv_leafn_to_prod1n     , & ! Input:  [real(r8) (:)] crop leafn harvested
         hrv_livestemn_to_prod1n  => veg_nf%hrv_livestemn_to_prod1n , & ! Input:  [real(r8) (:)] crop stemn harvested
         hrv_leafp_to_prod1p   =>    veg_pf%hrv_leafp_to_prod1p   , & ! Input:  [real(r8) (:)] crop leafp harvested
         hrv_livestemp_to_prod1p  => veg_pf%hrv_livestemp_to_prod1p & ! Input:  [real(r8) (:)] crop stemp harvested
         )

      ! The litterfall transfer rate starts at 0.0 and increases linearly
      ! over time, with displayed growth going to 0.0 on the last day of litterfall

      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! only calculate fluxes during offset period
         if (offset_flag(p) == 1._r8) then

            if (offset_counter(p) == dt) then
               t1 = 1.0_r8 / dt
               if (ivt(p) >= npcropmin) then
               ! this assumes that offset_counter == dt for crops
               ! if this were ever changed, we'd need to add code to the "else"
                  leafc_to_litter(p) = (1.0_r8 - presharv(ivt(p))) * ((t1 * leafc(p)) + cpool_to_leafc(p))
                  frootc_to_litter(p) = t1 * frootc(p) + cpool_to_frootc(p)
                  livestemc_to_litter(p) = (1.0_r8 - presharv(ivt(p))) * ((t1 * livestemc(p)) + cpool_to_livestemc(p))
               else
                  leafc_to_litter(p)  = t1 * leafc(p)  + cpool_to_leafc(p)
                  frootc_to_litter(p) = t1 * frootc(p) + cpool_to_frootc(p)
               end if
            else
               t1 = dt * 2.0_r8 / (offset_counter(p) * offset_counter(p))
               leafc_to_litter(p)  = prev_leafc_to_litter(p)  + t1*(leafc(p)  - prev_leafc_to_litter(p)*offset_counter(p))
               frootc_to_litter(p) = prev_frootc_to_litter(p) + t1*(frootc(p) - prev_frootc_to_litter(p)*offset_counter(p))
            end if

            if ( nu_com .eq. 'RD') then
               if (ivt(p) >= npcropmin) then
                  if (offset_counter(p) == dt) then
                      t1 = 1.0_r8 / dt

                     ! this assumes that offset_counter == dt for crops
                     ! if this were ever changed, we'd need to add code to the
                     ! "else"
                     leafn_to_litter(p) = (t1 * leafn(p) + npool_to_leafn(p)) - hrv_leafn_to_prod1n(p)
                     leafp_to_litter(p) = (t1 * leafp(p) + ppool_to_leafp(p)) - hrv_leafp_to_prod1p(p)

                     frootn_to_litter(p) = t1 * frootn(p) + npool_to_frootn(p)
                     frootp_to_litter(p) = t1 * frootp(p) + ppool_to_frootp(p)

                     livestemn_to_litter(p) = (t1 * livestemn(p) + npool_to_livestemn(p)) - hrv_livestemn_to_prod1n(p)
                     livestemp_to_litter(p) = (t1 * livestemp(p) + ppool_to_livestemp(p)) - hrv_livestemp_to_prod1p(p)
                  end if
               else
                  ! calculate the leaf N litterfall and retranslocation
                  leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
                  leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

                  ! calculate fine root N litterfall (no retranslocation of fine root N)
                  frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))

                  ! calculate the leaf P litterfall and retranslocation
                  leafp_to_litter(p)   = leafc_to_litter(p)  / lflitcp(ivt(p))
                  leafp_to_retransp(p) = (leafc_to_litter(p) / leafcp(ivt(p))) - leafp_to_litter(p)

                  ! calculate fine root P litterfall (no retranslocation of fine root N)
                  frootp_to_litter(p) = frootc_to_litter(p) / frootcp(ivt(p))
               end if
            else
               if (offset_counter(p) == dt) then
                  t1 = 1.0_r8 / dt
                  if (ivt(p) >= npcropmin) then
                     ! this assumes that offset_counter == dt for crops
                     ! if this were ever changed, we'd need to add code to the "else"
                     leafn_to_litter(p) = (1.0_r8 - presharv(ivt(p))) * ((t1 * leafn(p)) + npool_to_leafn(p))
                     leafp_to_litter(p) = (1.0_r8 - presharv(ivt(p))) * ((t1 * leafp(p)) + ppool_to_leafp(p))

                     frootn_to_litter(p) = t1 * frootn(p) + npool_to_frootn(p)
                     frootp_to_litter(p) = t1 * frootp(p) + ppool_to_frootp(p)

                     livestemn_to_litter(p) = (1.0_r8 - presharv(ivt(p))) * ((t1 * livestemn(p)) + npool_to_livestemn(p))
                     livestemp_to_litter(p) = (1.0_r8 - presharv(ivt(p))) * ((t1 * livestemp(p)) + ppool_to_livestemp(p))

                  else
                     leafn_to_litter(p)   = (t1 * leafn(p) + npool_to_leafn(p))*0.38_r8
                     leafn_to_retransn(p) = (t1 * leafn(p) + npool_to_leafn(p))*0.62_r8
                     frootn_to_litter(p)  = t1 * frootn(p) + npool_to_frootn(p)

                     leafp_to_litter(p)   = (t1 * leafp(p) + ppool_to_leafp(p))*0.35_r8
                     leafp_to_retransp(p) = (t1 * leafp(p) + ppool_to_leafp(p))*0.65_r8
                     frootp_to_litter(p)  = t1 * frootp(p) + ppool_to_frootp(p)
                  end if
               else
                  leafn_to_litter(p)   = leafc_to_litter(p) / max(leafc(p), 1.e-20_r8) * leafn(p) * 0.38_r8
                  leafn_to_retransn(p) = leafc_to_litter(p) / max(leafc(p), 1.e-20_r8) * leafn(p) * 0.62_r8
                  frootn_to_litter(p)  = frootc_to_litter(p)/ max(frootc(p), 1.e-20_r8) * frootn(p)

                  leafp_to_litter(p)   = leafc_to_litter(p) / max(leafc(p), 1.e-20_r8) * leafp(p) * 0.35_r8
                  leafp_to_retransp(p) = leafc_to_litter(p) / max(leafc(p), 1.e-20_r8) * leafp(p) * 0.65_r8
                  frootp_to_litter(p)  = frootc_to_litter(p)/ max(frootc(p), 1.e-20_r8) * frootp(p)
               end if
            end if
            leafn_to_litter(p)     = leafn_to_litter(p) * pheno_indicator(pid_leafn_to_litter)
            frootn_to_litter(p)    = frootn_to_litter(p) * pheno_indicator(pid_frootn_to_litter)
            livestemn_to_litter(p) = livestemn_to_litter(p) * pheno_indicator(pid_livestemn_to_litter)

            ! save the current litterfall fluxes
            prev_leafc_to_litter(p)  = leafc_to_litter(p)
            prev_frootc_to_litter(p) = frootc_to_litter(p)

         end if ! end if offset period

      end do ! end pft loop

    end associate

  end subroutine CNOffsetLitterfall

  !-----------------------------------------------------------------------
  subroutine CNBackgroundLitterfall (num_soilp, filter_soilp, &
       cnstate_vars, carbonstate_vars, carbonflux_vars, nitrogenflux_vars,&
       phosphorusflux_vars, nitrogenstate_vars, phosphorusstate_vars)
    !
    ! !DESCRIPTION:
    ! Determines the flux of C and N from displayed pools to litter
    ! pools as the result of background litter fall.
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                 , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(cnstate_type)      , intent(in)    :: cnstate_vars
    type(carbonstate_type)  , intent(in)    :: carbonstate_vars
    type(carbonflux_type)   , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type) , intent(inout) :: phosphorusflux_vars
    type(nitrogenstate_type)  , intent(in)    :: nitrogenstate_vars
    type(phosphorusstate_type), intent(in)    :: phosphorusstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: fp           ! lake filter pft index
    !-----------------------------------------------------------------------

    associate(                                                               &
         ivt               =>    veg_pp%itype                                 , & ! Input:  [integer  (:) ]  pft vegetation type

         leafcn            =>    veg_vp%leafcn                         , & ! Input:  [real(r8) (:) ]  leaf C:N (gC/gN)
         lflitcn           =>    veg_vp%lflitcn                        , & ! Input:  [real(r8) (:) ]  leaf litter C:N (gC/gN)
         frootcn           =>    veg_vp%frootcn                        , & ! Input:  [real(r8) (:) ]  fine root C:N (gC/gN)

         leafcp            =>    veg_vp%leafcp                         , & ! Input:  [real(r8) (:) ]  leaf C:P (gC/gP)
         lflitcp           =>    veg_vp%lflitcp                        , & ! Input:  [real(r8) (:) ]  leaf litter C:P (gC/gP)
         frootcp           =>    veg_vp%frootcp                        , & ! Input:  [real(r8) (:) ]  fine root C:P (gC/gP)

         bglfr             =>    cnstate_vars%bglfr_patch                  , & ! Input:  [real(r8) (:) ]  background litterfall rate (1/s)

         leafc             =>    veg_cs%leafc              , & ! Input:  [real(r8) (:) ]  (gC/m2) leaf C
         frootc            =>    veg_cs%frootc             , & ! Input:  [real(r8) (:) ]  (gC/m2) fine root C

         leafc_to_litter   =>    veg_cf%leafc_to_litter     , & ! Output: [real(r8) (:) ]
         frootc_to_litter  =>    veg_cf%frootc_to_litter    , & ! Output: [real(r8) (:) ]

         leafn_to_litter   =>    veg_nf%leafn_to_litter   , & ! Output: [real(r8) (:) ]
         leafn_to_retransn =>    veg_nf%leafn_to_retransn , & ! Output: [real(r8) (:) ]
         frootn_to_litter  =>    veg_nf%frootn_to_litter  , & ! Output: [real(r8) (:) ]

         leafp_to_litter   =>    veg_pf%leafp_to_litter   , & ! Output: [real(r8) (:) ]
         leafp_to_retransp =>    veg_pf%leafp_to_retransp , & ! Output: [real(r8) (:) ]
         frootp_to_litter  =>    veg_pf%frootp_to_litter  , & ! Output: [real(r8) (:) ]

         leafn             =>    veg_ns%leafn              , &
         frootn            =>    veg_ns%frootn             , &
         leafp             =>    veg_ps%leafp            , &
         frootp            =>    veg_ps%frootp             &
         )

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! only calculate these fluxes if the background litterfall rate is non-zero
         if (bglfr(p) > 0._r8) then
            ! units for bglfr are already 1/s
            leafc_to_litter(p)  = bglfr(p) * leafc(p)
            frootc_to_litter(p) = bglfr(p) * frootc(p)

            if ( nu_com .eq. 'RD') then
               ! calculate the leaf N litterfall and retranslocation
               leafn_to_litter(p)   = leafc_to_litter(p)  / lflitcn(ivt(p))
               leafn_to_litter(p)   = leafn_to_litter(p) * pheno_indicator(pid_leafn_to_litter)
               leafn_to_retransn(p) = (leafc_to_litter(p) / leafcn(ivt(p))) - leafn_to_litter(p)

               ! calculate fine root N litterfall (no retranslocation of fine root N)
               frootn_to_litter(p) = frootc_to_litter(p) / frootcn(ivt(p))

               ! calculate the leaf P litterfall and retranslocation
               leafp_to_litter(p)   = leafc_to_litter(p)  / lflitcp(ivt(p))
               leafp_to_retransp(p) = (leafc_to_litter(p) / leafcp(ivt(p))) - leafp_to_litter(p)

               ! calculate fine root P litterfall (no retranslocation of fine root P)
               frootp_to_litter(p) = frootc_to_litter(p) / frootcp(ivt(p))
            else
               ! calculate the leaf N litterfall and retranslocation
               leafn_to_litter(p)   = bglfr(p) * leafn(p) * 0.38_r8 ! 62% N resorption rate; LEONARDUS VERGUTZ 2012 Ecological Monographs 82(2) 205-220.
               leafn_to_litter(p)   = leafn_to_litter(p) * pheno_indicator(pid_leafn_to_litter)
               leafn_to_retransn(p) = bglfr(p) * leafn(p) - leafn_to_litter(p)

               ! calculate fine root N litterfall (no retranslocation of fine root N)
               frootn_to_litter(p) = bglfr(p) * frootn(p)

               ! calculate the leaf P litterfall and retranslocation
               leafp_to_litter(p)   = bglfr(p) * leafp(p) * 0.35_r8 ! 65% P resorption rate; LEONARDUS VERGUTZ 2012 Ecological Monographs 82(2) 205-220.
               leafp_to_retransp(p) = bglfr(p) * leafp(p) - leafp_to_litter(p)

               ! calculate fine root P litterfall (no retranslocation of fine root P)
               frootp_to_litter(p) = bglfr(p) * frootp(p) ! fine root P retranslocation occur (but not N retranslocation), why not include it here
            end if
         end if
         frootn_to_litter(p) = frootn_to_litter(p) * pheno_indicator(pid_frootn_to_litter)
      end do

    end associate

  end subroutine CNBackgroundLitterfall

  !-----------------------------------------------------------------------
  subroutine CNLivewoodTurnover (num_soilp, filter_soilp, &
       carbonstate_vars, nitrogenstate_vars, carbonflux_vars,nitrogenflux_vars,&
       phosphorusstate_vars,phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! Determines the flux of C and N from live wood to
    ! dead wood pools, for stem and coarse root.
    ! add phosphorus flux - X.YANG
    !
    ! !ARGUMENTS:
    integer                  , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                  , intent(in)    :: filter_soilp(:) ! filter for soil patches
    type(carbonstate_type)   , intent(in)    :: carbonstate_vars
    type(nitrogenstate_type) , intent(in)    :: nitrogenstate_vars
    type(carbonflux_type)    , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
    type(phosphorusstate_type) , intent(in)    :: phosphorusstate_vars
    type(phosphorusflux_type)  , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p            ! indices
    integer :: fp           ! lake filter pft index
    real(r8):: ctovr        ! temporary variable for carbon turnover
    real(r8):: ntovr        ! temporary variable for nitrogen turnover
    real(r8):: ptovr        ! temporary variable for phosphorus turnover
    !-----------------------------------------------------------------------

    associate(                                                                             &
         ivt                      =>    veg_pp%itype                                        , & ! Input:  [integer  (:) ]  pft vegetation type

         woody                    =>    veg_vp%woody                                 , & ! Input:  [real(r8) (:) ]  binary flag for woody lifeform (1=woody, 0=not woody)
         livewdcn                 =>    veg_vp%livewdcn                              , & ! Input:  [real(r8) (:) ]  live wood (phloem and ray parenchyma) C:N (gC/gN)
         deadwdcn                 =>    veg_vp%deadwdcn                              , & ! Input:  [real(r8) (:) ]  dead wood (xylem and heartwood) C:N (gC/gN)

         livestemc                =>    veg_cs%livestemc                 , & ! Input:  [real(r8) (:) ]  (gC/m2) live stem C
         livecrootc               =>    veg_cs%livecrootc                , & ! Input:  [real(r8) (:) ]  (gC/m2) live coarse root C

         livestemn                =>    veg_ns%livestemn               , & ! Input:  [real(r8) (:) ]  (gN/m2) live stem N
         livecrootn               =>    veg_ns%livecrootn              , & ! Input:  [real(r8) (:) ]  (gN/m2) live coarse root N

         livestemp                =>    veg_ps%livestemp               , & ! Input:  [real(r8) (:) ]  (gN/m2) live stem N
         livecrootp               =>    veg_ps%livecrootp              , & ! Input:  [real(r8) (:) ]  (gN/m2) live coarse root N

         livestemc_to_deadstemc   =>    veg_cf%livestemc_to_deadstemc     , & ! Output: [real(r8) (:) ]
         livecrootc_to_deadcrootc =>    veg_cf%livecrootc_to_deadcrootc   , & ! Output: [real(r8) (:) ]

         livestemn_to_deadstemn   =>    veg_nf%livestemn_to_deadstemn   , & ! Output: [real(r8) (:) ]
         livestemn_to_retransn    =>    veg_nf%livestemn_to_retransn    , & ! Output: [real(r8) (:) ]
         livecrootn_to_deadcrootn =>    veg_nf%livecrootn_to_deadcrootn , & ! Output: [real(r8) (:) ]
         livecrootn_to_retransn   =>    veg_nf%livecrootn_to_retransn   , & ! Output: [real(r8) (:) ]

         livewdcp                 =>    veg_vp%livewdcp                              , & ! Input:  [real(r8) (:) ]  live wood (phloem and ray parenchyma) C:P (gC/gP)
         deadwdcp                 =>    veg_vp%deadwdcp                              , & ! Input:  [real(r8) (:) ]  dead wood (xylem and heartwood) C:P (gC/gP)
         livestemp_to_deadstemp   =>    veg_pf%livestemp_to_deadstemp   , & ! Output: [real(r8) (:) ]
         livestemp_to_retransp    =>    veg_pf%livestemp_to_retransp    , & ! Output: [real(r8) (:) ]
         livecrootp_to_deadcrootp =>    veg_pf%livecrootp_to_deadcrootp , & ! Output: [real(r8) (:) ]
         livecrootp_to_retransp   =>    veg_pf%livecrootp_to_retransp     & ! Output: [real(r8) (:) ]
         )

      ! patch loop
      do fp = 1,num_soilp
         p = filter_soilp(fp)

         ! only calculate these fluxes for woody types
         if (woody(ivt(p)) > 0._r8) then
            if ( nu_com .eq. 'RD') then
               ! live stem to dead stem turnover

               ctovr = livestemc(p) * lwtop
               ntovr = ctovr / livewdcn(ivt(p))
               ptovr = ctovr / livewdcp(ivt(p))

               livestemc_to_deadstemc(p) = ctovr
               livestemn_to_deadstemn(p) = ctovr / deadwdcn(ivt(p))
               livestemn_to_retransn(p)  = ntovr - livestemn_to_deadstemn(p)

               livestemp_to_deadstemp(p) = ctovr / deadwdcp(ivt(p))
               livestemp_to_retransp(p)  = ptovr - livestemp_to_deadstemp(p)
               ! live coarse root to dead coarse root turnover

               ctovr = livecrootc(p) * lwtop
               ntovr = ctovr / livewdcn(ivt(p))
               ptovr = ctovr / livewdcp(ivt(p))

               livecrootc_to_deadcrootc(p) = ctovr
               livecrootn_to_deadcrootn(p) = ctovr / deadwdcn(ivt(p))
               livecrootn_to_retransn(p)  = ntovr - livecrootn_to_deadcrootn(p)

               livecrootp_to_deadcrootp(p) = ctovr / deadwdcp(ivt(p))
               livecrootp_to_retransp(p)  = ptovr - livecrootp_to_deadcrootp(p)
            else
               ! live stem to dead stem turnover

               ctovr = livestemc(p) * lwtop
               ntovr = livestemn(p) * lwtop
               ptovr = livestemp(p) * lwtop

               livestemc_to_deadstemc(p) = ctovr
               livestemn_to_deadstemn(p) = ntovr * livewdcn(ivt(p))/deadwdcn(ivt(p)) ! N retranslocation
               livestemn_to_retransn(p)  = ntovr - livestemn_to_deadstemn(p)

               livestemp_to_deadstemp(p) = ptovr*  livewdcp(ivt(p))/deadwdcp(ivt(p)) ! P retranslocation
               livestemp_to_retransp(p)  = ptovr - livestemp_to_deadstemp(p)
               ! live coarse root to dead coarse root turnover

               ctovr = livecrootc(p) * lwtop
               ntovr = livecrootn(p) * lwtop
               ptovr = livecrootp(p) * lwtop

               livecrootc_to_deadcrootc(p) = ctovr
               livecrootn_to_deadcrootn(p) = ntovr * livewdcn(ivt(p))/deadwdcn(ivt(p)) ! N retranslocation
               livecrootn_to_retransn(p)  = ntovr - livecrootn_to_deadcrootn(p)

               livecrootp_to_deadcrootp(p) = ptovr *  livewdcp(ivt(p))/deadwdcp(ivt(p)) ! P retranslocation
               livecrootp_to_retransp(p)  = ptovr - livecrootp_to_deadcrootp(p)
            end if

         end if

      end do

    end associate

  end subroutine CNLivewoodTurnover

  !-----------------------------------------------------------------------
  subroutine CNLitterToColumn (num_soilc, filter_soilc, &
       cnstate_vars, carbonflux_vars, nitrogenflux_vars,phosphorusflux_vars)
    !
    ! !DESCRIPTION:
    ! called at the end of cn_phenology to gather all pft-level litterfall fluxes
    ! to the column level and assign them to the three litter pools
    !
    ! !USES:
    use elm_varpar , only : max_patch_per_col, nlevdecomp
    use pftvarcon  , only : npcropmin
    !
    ! !ARGUMENTS:
    integer                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(cnstate_type)      , intent(in)    :: cnstate_vars
    type(carbonflux_type)   , intent(inout) :: carbonflux_vars
    type(nitrogenflux_type) , intent(inout) :: nitrogenflux_vars
    type(phosphorusflux_type) , intent(inout) :: phosphorusflux_vars
    !
    ! !LOCAL VARIABLES:
    integer :: fc,c,pi,p,j       ! indices
    !-----------------------------------------------------------------------

    associate(                                                                                       &
         ivt                                 =>    veg_pp%itype                                       , & ! Input:  [integer  (:)   ]  pft vegetation type
         wtcol                               =>    veg_pp%wtcol                                       , & ! Input:  [real(r8) (:)   ]  weight (relative to column) for this pft (0-1)

         lf_flab                             =>    veg_vp%lf_flab                              , & ! Input:  [real(r8) (:)   ]  leaf litter labile fraction
         lf_fcel                             =>    veg_vp%lf_fcel                              , & ! Input:  [real(r8) (:)   ]  leaf litter cellulose fraction
         lf_flig                             =>    veg_vp%lf_flig                              , & ! Input:  [real(r8) (:)   ]  leaf litter lignin fraction
         fr_flab                             =>    veg_vp%fr_flab                              , & ! Input:  [real(r8) (:)   ]  fine root litter labile fraction
         fr_fcel                             =>    veg_vp%fr_fcel                              , & ! Input:  [real(r8) (:)   ]  fine root litter cellulose fraction
         fr_flig                             =>    veg_vp%fr_flig                              , & ! Input:  [real(r8) (:)   ]  fine root litter lignin fraction

         leaf_prof                           =>    cnstate_vars%leaf_prof_patch                    , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of leaves
         froot_prof                          =>    cnstate_vars%froot_prof_patch                   , & ! Input:  [real(r8) (:,:) ]  (1/m) profile of fine roots

         leafc_to_litter                     =>    veg_cf%leafc_to_litter           , & ! Input:  [real(r8) (:)   ]  leaf C litterfall (gC/m2/s)
         frootc_to_litter                    =>    veg_cf%frootc_to_litter          , & ! Input:  [real(r8) (:)   ]  fine root N litterfall (gN/m2/s)
         livestemc_to_litter                 =>    veg_cf%livestemc_to_litter       , & ! Input:  [real(r8) (:)   ]  live stem C litterfall (gC/m2/s)
!         grainc_to_food                      =>    veg_cf%grainc_to_food            , & ! Input:  [real(r8) (:)   ]  grain C to food (gC/m2/s)
         phenology_c_to_litr_met_c           =>    col_cf%phenology_c_to_litr_met_c   , & ! Output: [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gC/m3/s)
         phenology_c_to_litr_cel_c           =>    col_cf%phenology_c_to_litr_cel_c   , & ! Output: [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gC/m3/s)
         phenology_c_to_litr_lig_c           =>    col_cf%phenology_c_to_litr_lig_c   , & ! Output: [real(r8) (:,:) ]  C fluxes associated with phenology (litterfall and crop) to litter lignin pool (gC/m3/s)

         livestemn_to_litter                 =>    veg_nf%livestemn_to_litter     , & ! Input:  [real(r8) (:)   ]  livestem N to litter (gN/m2/s)
!         grainn_to_food                      =>    veg_nf%grainn_to_food          , & ! Input:  [real(r8) (:)   ]  grain N to food (gN/m2/s)
         leafn_to_litter                     =>    veg_nf%leafn_to_litter         , & ! Input:  [real(r8) (:)   ]  leaf N litterfall (gN/m2/s)
         frootn_to_litter                    =>    veg_nf%frootn_to_litter        , & ! Input:  [real(r8) (:)   ]  fine root N litterfall (gN/m2/s)
         phenology_n_to_litr_met_n           =>    col_nf%phenology_n_to_litr_met_n , & ! Output: [real(r8) (:,:) ]  N fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gN/m3/s)
         phenology_n_to_litr_cel_n           =>    col_nf%phenology_n_to_litr_cel_n , & ! Output: [real(r8) (:,:) ]  N fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gN/m3/s)
         phenology_n_to_litr_lig_n           =>    col_nf%phenology_n_to_litr_lig_n , & ! Output: [real(r8) (:,:) ]  N fluxes associated with phenology (litterfall and crop) to litter lignin pool (gN/m3/s)

         livestemp_to_litter                 =>    veg_pf%livestemp_to_litter     , & ! Input:  [real(r8) (:)   ]  livestem P to litter (gP/m2/s)
!         grainp_to_food                      =>    veg_pf%grainp_to_food          , & ! Input:  [real(r8) (:)   ]  grain P to food (gP/m2/s)
         leafp_to_litter                     =>    veg_pf%leafp_to_litter         , & ! Input:  [real(r8) (:)   ]  leaf P litterfall (gP/m2/s)
         frootp_to_litter                    =>    veg_pf%frootp_to_litter        , & ! Input:  [real(r8) (:)   ]  fine root P litterfall (gP/m2/s)
         phenology_p_to_litr_met_p           =>    col_pf%phenology_p_to_litr_met_p , & ! Output: [real(r8) (:,:) ]  P fluxes associated with phenology (litterfall and crop) to litter metabolic pool (gP/m3/s)
         phenology_p_to_litr_cel_p           =>    col_pf%phenology_p_to_litr_cel_p , & ! Output: [real(r8) (:,:) ]  P fluxes associated with phenology (litterfall and crop) to litter cellulose pool (gP/m3/s)
         phenology_p_to_litr_lig_p           =>    col_pf%phenology_p_to_litr_lig_p   & ! Output: [real(r8) (:,:) ]  P fluxes associated with phenology (litterfall and crop) to litter lignin pool (gP/m3/s)
         )

         if(.false.)then
         do pi = 1,max_patch_per_col
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if ( pi <=  col_pp%npfts(c) ) then
                  p = col_pp%pfti(c) + pi - 1
                  if (veg_pp%active(p)) then
                     write(*,*)'=========================================='
                     write(*,*)'01 leafn_to_litter    =',p,leafn_to_litter(p)
                     write(*,*)'02 frootn_to_litter   =',p,frootn_to_litter(p)
                     !write(*,*)'03 livestemn_to_litter=',p,livestemn_to_litter(p)
                  endif
                endif
             enddo
          enddo
        endif
      do j = 1, nlevdecomp
         do pi = 1,max_patch_per_col
            do fc = 1,num_soilc
               c = filter_soilc(fc)

               if ( pi <=  col_pp%npfts(c) ) then
                  p = col_pp%pfti(c) + pi - 1
                  if (veg_pp%active(p)) then

                     ! leaf litter carbon fluxes
                     phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                          + leafc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                          + leafc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                          + leafc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     ! leaf litter nitrogen fluxes
                     phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                          + leafn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                          + leafn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                          + leafn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     ! leaf litter phosphorus fluxes
                     phenology_p_to_litr_met_p(c,j) = phenology_p_to_litr_met_p(c,j) &
                          + leafp_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_p_to_litr_cel_p(c,j) = phenology_p_to_litr_cel_p(c,j) &
                          + leafp_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                     phenology_p_to_litr_lig_p(c,j) = phenology_p_to_litr_lig_p(c,j) &
                          + leafp_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     ! fine root litter carbon fluxes
                     phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                          + frootc_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                          + frootc_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                          + frootc_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                     ! fine root litter nitrogen fluxes
                     phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                          + frootn_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                          + frootn_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                          + frootn_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                     ! fine root litter phosphorus fluxes
                     phenology_p_to_litr_met_p(c,j) = phenology_p_to_litr_met_p(c,j) &
                          + frootp_to_litter(p) * fr_flab(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_p_to_litr_cel_p(c,j) = phenology_p_to_litr_cel_p(c,j) &
                          + frootp_to_litter(p) * fr_fcel(ivt(p)) * wtcol(p) * froot_prof(p,j)
                     phenology_p_to_litr_lig_p(c,j) = phenology_p_to_litr_lig_p(c,j) &
                          + frootp_to_litter(p) * fr_flig(ivt(p)) * wtcol(p) * froot_prof(p,j)

                     ! agroibis puts crop stem litter together with leaf litter
                     ! so I've used the leaf lf_f* parameters instead of making
                     ! new ones for now (slevis)
                     ! The food is now directed to the product pools (BDrewniak)

                     if (ivt(p) >= npcropmin) then ! add livestemc to litter
                        ! stem litter carbon fluxes
                        phenology_c_to_litr_met_c(c,j) = phenology_c_to_litr_met_c(c,j) &
                             + livestemc_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        phenology_c_to_litr_cel_c(c,j) = phenology_c_to_litr_cel_c(c,j) &
                             + livestemc_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        phenology_c_to_litr_lig_c(c,j) = phenology_c_to_litr_lig_c(c,j) &
                             + livestemc_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                        ! stem litter nitrogen fluxes
                        phenology_n_to_litr_met_n(c,j) = phenology_n_to_litr_met_n(c,j) &
                             + livestemn_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        phenology_n_to_litr_cel_n(c,j) = phenology_n_to_litr_cel_n(c,j) &
                             + livestemn_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        phenology_n_to_litr_lig_n(c,j) = phenology_n_to_litr_lig_n(c,j) &
                             + livestemn_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                        ! stem litter phosphorus fluxes
                        phenology_p_to_litr_met_p(c,j) = phenology_p_to_litr_met_p(c,j) &
                             + livestemp_to_litter(p) * lf_flab(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        phenology_p_to_litr_cel_p(c,j) = phenology_p_to_litr_cel_p(c,j) &
                             + livestemp_to_litter(p) * lf_fcel(ivt(p)) * wtcol(p) * leaf_prof(p,j)
                        phenology_p_to_litr_lig_p(c,j) = phenology_p_to_litr_lig_p(c,j) &
                             + livestemp_to_litter(p) * lf_flig(ivt(p)) * wtcol(p) * leaf_prof(p,j)

                     end if

                  end if
               end if

            end do

         end do
      end do

    end associate

  end subroutine CNLitterToColumn

 !-----------------------------------------------------------------------
 subroutine CNCropHarvestPftToColumn (num_soilc, filter_soilc, &
            cnstate_vars, carbonflux_vars, nitrogenflux_vars, phosphorusflux_vars)
   !
   ! !DESCRIPTION:
   ! called at the end of CNCropHarvest to gather all pft-level harvest fluxes
   ! to the column level and assign them to a product pools
   !
   ! !USES:
   use elm_varpar, only : maxpatch_pft
   type(cnstate_type)       , intent(in)    :: cnstate_vars
   type(carbonflux_type)    , intent(inout) :: carbonflux_vars
   type(nitrogenflux_type)  , intent(inout) :: nitrogenflux_vars
   type(phosphorusflux_type), intent(inout) :: phosphorusflux_vars
   !
   ! !ARGUMENTS:
   integer, intent(in) :: num_soilc       ! number of soil columns in filter
   integer, intent(in) :: filter_soilc(:) ! soil column filter
   ! !LOCAL VARIABLES:
   integer :: fc,c,pi,p                   ! indices
   !-----------------------------------------------------------------------

   associate(&
   ivt                                 =>   veg_pp%itype                                    , & ! Input:  [integer (:)]  pft vegetation type
   wtcol                               =>   veg_pp%wtcol                                    , & ! Input:  [real(r8) (:)]  pft weight relative to column (0-1)
   phrv_leafc_to_prod1c                =>   veg_cf%hrv_leafc_to_prod1c    , & ! Input:  [real(r8) (:)] crop leafc harvested
   phrv_livestemc_to_prod1c            =>   veg_cf%hrv_livestemc_to_prod1c, & ! Input:  [real(r8) (:)] crop stemc harvested
   phrv_grainc_to_prod1c               =>   veg_cf%hrv_grainc_to_prod1c   , & ! Input:  [real(r8) (:)] crop grainc harvested
   phrv_cropc_to_prod1c                =>   veg_cf%hrv_cropc_to_prod1c    , & ! InOut:  [real(r8) (:)] crop carbon harvested
   phrv_leafn_to_prod1n                =>   veg_nf%hrv_leafn_to_prod1n    , & ! Input:  [real(r8) (:)] crop leafn harvested
   phrv_livestemn_to_prod1n            =>   veg_nf%hrv_livestemn_to_prod1n, & ! Input:  [real(r8) (:)] crop stemn harvested
   phrv_grainn_to_prod1n               =>   veg_nf%hrv_grainn_to_prod1n   , & ! Input:  [real(r8) (:)] crop grainn harvested
   phrv_cropn_to_prod1n                =>   veg_nf%hrv_cropn_to_prod1n    , & ! InOut:  [real(r8) (:)] crop grainn harvested
   phrv_leafp_to_prod1p                =>   veg_pf%hrv_leafp_to_prod1p  , & ! InOut:  [real(r8) (:)] crop grainp harvested
   phrv_livestemp_to_prod1p            =>   veg_pf%hrv_livestemp_to_prod1p, & ! InOut:  [real(r8) (:)] column level crop carbon harvested
   phrv_grainp_to_prod1p               =>   veg_pf%hrv_grainp_to_prod1p , & ! InOut:  [real(r8) (:)] column level crop nitrogen harvested
   phrv_cropp_to_prod1p                =>   veg_pf%hrv_cropp_to_prod1p  , & ! InOut:  [real(r8) (:)] column level crop phosphorus harvested
   chrv_cropc_to_prod1c                =>   col_cf%hrv_cropc_to_prod1c        , & ! InOut:  [real(r8) (:)] column level crop carbon harvested
   chrv_cropn_to_prod1n                =>   col_nf%hrv_cropn_to_prod1n      , & ! InOut:  [real(r8) (:)] column level crop nitrogen harvested
   chrv_cropp_to_prod1p                =>   col_pf%hrv_cropp_to_prod1p      & ! InOut:  [real(r8) (:)] column level crop phosphorus harvested
   )

   do pi = 1,maxpatch_pft
      do fc = 1,num_soilc
         c = filter_soilc(fc)

         if (pi <=  col_pp%npfts(c)) then
            p = col_pp%pfti(c) + pi - 1

            if (veg_pp%active(p)) then

                phrv_cropc_to_prod1c(p) = phrv_leafc_to_prod1c(p) + phrv_livestemc_to_prod1c(p) + &
                                         phrv_grainc_to_prod1c(p)

                chrv_cropc_to_prod1c(c) = chrv_cropc_to_prod1c(c) + phrv_cropc_to_prod1c(p) * wtcol(p)

                phrv_cropn_to_prod1n(p) = phrv_leafn_to_prod1n(p) + phrv_livestemn_to_prod1n(p) + &
                                         phrv_grainn_to_prod1n(p)

                chrv_cropn_to_prod1n(c) = chrv_cropn_to_prod1n(c) + phrv_cropn_to_prod1n(p) * wtcol(p)

                phrv_cropp_to_prod1p(p) = phrv_leafp_to_prod1p(p) + phrv_livestemp_to_prod1p(p) + &
                                         phrv_grainp_to_prod1p(p)

                chrv_cropp_to_prod1p(c) = chrv_cropp_to_prod1p(c) + phrv_cropp_to_prod1p(p) * wtcol(p)

            end if
         end if

      end do

   end do
 end associate
end subroutine CNCropHarvestPftToColumn

end module CNPhenologyBeTRMod
