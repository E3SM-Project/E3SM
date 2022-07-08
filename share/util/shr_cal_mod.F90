! !MODULE: shr_cal_mod -- calendar module, relates elapsed days to calendar date.
!
! !DESCRIPTION:
!   These calendar routines do conversions between...
!   \begin{itemize}
!   \item the integer number of elapsed days
!   \item the integers year, month, day (three inter-related integers)
!   \item the integer coded calendar date (yyyymmdd)
!   \end{itemize}
!   Possible uses include: a calling routine can increment the elapsed days
!   integer and use this module to determine what the corresponding calendar
!   date is;  this module can be used to determine how many days apart two
!   arbitrary calendar dates are.
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - created initial version, taken from cpl5
!     2012-feb-10 - T. Craig - update to esmf 5.2.0rp1, some noleap functions
!       associated with elapsed days must still be done internally
!
! !REMARKS:
!   Following are some internal assumptions.  These assumptions are somewhat
!   arbitrary -- they were chosen because they result in the simplest code given
!   the requirements of this module.  These assumptions can be relaxed as
!   necessary:
!   o the valid range of years is [-999,999999]
!   o elapsed days = 0 <=> January 1st, year 0000 for noleap
!
! !INTERFACE: ------------------------------------------------------------------

module shr_cal_mod

  ! !USES:

  use shr_kind_mod   ! kinds
  use shr_const_mod  ! constants
  use shr_sys_mod    ! system
  use shr_string_mod, only: shr_string_toLower
  use shr_log_mod, only: s_loglev  => shr_log_Level
  use shr_log_mod, only: s_logunit => shr_log_Unit
  use esmf

  implicit none

  private ! except

  ! !PUBLIC TYPES:

  type, public :: calParamType
     ! parameters to replace use of numbers for dates in code

     ! We'd like to have these be parameters, but then they
     ! can't be part of a derived type
     integer :: january     = 1
     integer :: february    = 2
     integer :: march       = 3
     integer :: april       = 4
     integer :: may         = 5
     integer :: june        = 6
     integer :: july        = 7
     integer :: august      = 8
     integer :: september   = 9
     integer :: october     = 10
     integer :: november    = 11
     integer :: december    = 12
     integer :: firstDayOfMonth = 1
  end type calParamType

  ! Instance of calender parameters as a protected type so can't change outside this module
  type(calParamType), public, protected :: calParams

  ! !PUBLIC MEMBER FUNCTIONS:

  public :: shr_cal_calendarName       ! checks and returns valid cal name
  public :: shr_cal_numDaysinMonth     ! number of days in a month
  public :: shr_cal_numDaysinYear      ! number of days in a year
  public :: shr_cal_elapsDaysStrtMonth ! elapsed days on start of month
  public :: shr_cal_timeset     ! set ESMF Time from ymd, s, shr_cal calendar
  public :: shr_cal_date2ymd   ! converts coded-date   to yr,month,day
  public :: shr_cal_date2julian! converts coded-date,sec to julian days
  public :: shr_cal_ymd2julian ! converts yr,month,day,sec to julian days
  public :: shr_cal_ymd2date   ! converts yr,month,day to coded-date
  public :: shr_cal_advDate    ! advance date/secs real seconds
  public :: shr_cal_advDateInt ! advance date/secs integer seconds
  public :: shr_cal_validDate  ! logical function: is coded-date valid?
  public :: shr_cal_validYMD   ! logical function: are yr,month,day valid?
  public :: shr_cal_validHMS   ! logical function: are hr, min, sec valid?
  public :: shr_cal_getDebug   ! get internal debug level
  public :: shr_cal_setDebug   ! set internal debug level
  public :: shr_cal_ymdtod2string ! translate ymdtod to string for filenames
  public :: shr_cal_datetod2string ! translate date to string for filenames
  public :: shr_cal_ymds2rday_offset ! translate yr,month,day,sec offset to a fractional day offset

  ! !PUBLIC DATA MEMBERS:

  ! none

  !EOP

  ! ! subroutine interfaces
  interface shr_cal_timeSet
     module procedure shr_cal_timeSet_int
     module procedure shr_cal_timeSet_long
  end interface shr_cal_timeSet

  interface shr_cal_date2ymd
     module procedure shr_cal_date2ymd_int
     module procedure shr_cal_date2ymd_long
  end interface shr_cal_date2ymd

  interface shr_cal_date2julian
     module procedure shr_cal_date2julian_int
     module procedure shr_cal_date2julian_long
  end interface shr_cal_date2julian

  interface shr_cal_ymd2date
     module procedure shr_cal_ymd2date_int
     module procedure shr_cal_ymd2date_long
  end interface shr_cal_ymd2date

  interface shr_cal_advDate
     module procedure shr_cal_advDate_int
     module procedure shr_cal_advDate_long
  end interface shr_cal_advDate

  interface shr_cal_advDateInt
     module procedure shr_cal_advDateInt_int
     module procedure shr_cal_advDateInt_long
  end interface shr_cal_advDateInt

  interface shr_cal_validdate
     module procedure shr_cal_validdate_int
     module procedure shr_cal_validdate_long
  end interface shr_cal_validdate

  interface  shr_cal_datetod2string
     module procedure shr_cal_datetod2string_int
     module procedure shr_cal_datetod2string_long
  end interface shr_cal_datetod2string

  integer(SHR_KIND_IN),parameter,public :: shr_cal_calMaxLen = 64
  character(len=*),parameter,public :: &
       shr_cal_noleap    = 'NO_LEAP', &
       shr_cal_gregorian = 'GREGORIAN'

  !--- trigger internal debug output ---
  integer(SHR_KIND_IN) :: debug = 0

  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_calendarName - check calendar name and translate
  !
  ! !DESCRIPTION:
  !    Check the validity of the calendar name and translate to standard naming
  !
  ! !REVISION HISTORY:
  !     2010-oct-7 - T Craig - initial version.
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  function shr_cal_calendarName(calendar,trap)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    character(*)         ,intent(in) :: calendar  ! calendar type
    logical     ,optional,intent(in) :: trap
    character(len=shr_cal_calMaxLen) :: shr_cal_calendarName

    !EOP

    character(len=shr_cal_calMaxLen) :: lcal
    logical :: ltrap
    character(len=shr_cal_calMaxLen) :: lowercalendar
    character(*),parameter :: subName = "(shr_cal_calendarName)"

    ltrap = .true.
    if (present(trap)) then
       ltrap = trap
    endif

    lcal = ' '
    lowercalendar = trim(shr_string_toLower(trim(calendar)))
    lcal = trim(calendar)

    selectcase(trim(lowercalendar))

    case ('noleap')
       lcal = trim(shr_cal_noleap)
    case ('no_leap')
       lcal = trim(shr_cal_noleap)
    case ('365_day')
       lcal = trim(shr_cal_noleap)
    case ('365day')
       lcal = trim(shr_cal_noleap)
    case ('gregorian')
       lcal = trim(shr_cal_gregorian)
    case ('standard')
       lcal = trim(shr_cal_gregorian)
    case ('proleptic_gregorian')
       lcal = trim(shr_cal_gregorian)
    case default
       if (ltrap) then
          write(s_logunit,*) trim(subname),' : ERROR calendar not supported : ',trim(calendar)
          call shr_sys_abort(trim(subname)//' : calendar not supported : '//trim(calendar))
       endif
    end select

    shr_cal_calendarName = trim(lcal)

  end function shr_cal_calendarName

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_timeSet - create and esmf_time type for ymd, sec and calendar
  !
  ! !DESCRIPTION:
  !    Create ESMF_Time type from ymd, sec, calendar
  !
  ! !REVISION HISTORY:
  !     2012-feb-10 - T. Craig - initial version.
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  subroutine shr_cal_timeSet_int(etime,ymd,sec,calendar, rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Time),intent(out) :: etime
    integer(SHR_KIND_IN),intent(in ) :: ymd,sec   ! ymd, sec
    character(*)        ,intent(in)  :: calendar  ! calendar type
    integer, intent(out), optional   :: rc
    !EOP

    integer(SHR_KIND_IN) :: year,month,day
    type(ESMF_CALKIND_FLAG) :: calkind
    character(len=shr_cal_calMaxLen) :: lcalendar
    integer :: lrc
    character(*),parameter :: subName = "(shr_cal_timeSet)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------


    lcalendar = shr_cal_calendarName(calendar)

    calkind = ESMF_CALKIND_NOLEAP
    if (trim(lcalendar) == trim(shr_cal_gregorian)) then
       calkind = ESMF_CALKIND_GREGORIAN
    endif

    call shr_cal_date2ymd(ymd,year,month,day)
    call ESMF_TimeSet(etime,yy=year,mm=month,dd=day,s=sec,calkindflag=calkind,rc=lrc)
    if (present(rc)) then
       rc = lrc
    else if(lrc /= ESMF_SUCCESS) then
       call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    endif
  end subroutine shr_cal_timeSet_int

  subroutine shr_cal_timeSet_long(etime,ymd,sec,calendar, rc)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Time),intent(out) :: etime
    integer(SHR_KIND_I8),intent(in ) :: ymd ! ymd
    integer(SHR_KIND_IN), intent(in) :: sec   ! ymd
    character(*)        ,intent(in)  :: calendar  ! calendar type
    integer, intent(out), optional   :: rc
    !EOP

    integer(SHR_KIND_IN) :: year,month,day
    type(ESMF_CALKIND_FLAG) :: calkind
    character(len=shr_cal_calMaxLen) :: lcalendar
    integer :: lrc
    character(*),parameter :: subName = "(shr_cal_timeSet_long)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------


    lcalendar = shr_cal_calendarName(calendar)

    calkind = ESMF_CALKIND_NOLEAP
    if (trim(lcalendar) == trim(shr_cal_gregorian)) then
       calkind = ESMF_CALKIND_GREGORIAN
    endif

    call shr_cal_date2ymd(ymd,year,month,day)
    call ESMF_TimeSet(etime,yy=year,mm=month,dd=day,s=sec,calkindflag=calkind,rc=lrc)
    if(present(rc)) then
       rc = lrc
    else if(rc /= ESMF_SUCCESS) then
       call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    endif
  end subroutine shr_cal_timeSet_long
  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_numDaysInMonth - return the number of days in a month.
  !
  ! !DESCRIPTION:
  !    Deturn the number of days in a month.
  !
  ! !REVISION HISTORY:
  !     2002-sep-18 - B. Kauffman - initial version.
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  integer function shr_cal_numDaysInMonth(year,month,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in ) :: year,month  ! calendar year,month
    character(*)        ,intent(in)  :: calendar  ! calendar type

    !EOP

    type(ESMF_time) :: time1,time2
    type(ESMF_timeInterval) :: timeint
    type(ESMF_CALKIND_FLAG) :: calkind
    integer(SHR_KIND_IN) :: eday
    character(len=shr_cal_calMaxLen) :: lcalendar
    integer :: rc
    character(*),parameter :: subName = "(shr_cal_numDaysInMonth)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lcalendar = shr_cal_calendarName(calendar)

    calkind = ESMF_CALKIND_NOLEAP
    if (trim(lcalendar) == trim(shr_cal_gregorian)) then
       calkind = ESMF_CALKIND_GREGORIAN
    endif

    call ESMF_TimeSet(time1,yy=year,mm=month,dd=1,s=0,calkindflag=calkind,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    if (month < 12) then
       call ESMF_TimeSet(time2,yy=year,mm=month+1,dd=1,s=0,calkindflag=calkind,rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    else
       call ESMF_TimeSet(time2,yy=year+1,mm=1,dd=1,s=0,calkindflag=calkind,rc=rc)
       if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    endif
    timeint = time2-time1
    call ESMF_TimeIntervalGet(timeint,StartTimeIn=time1,d=eday)
    shr_cal_numDaysInMonth = eday

  end function shr_cal_numDaysInMonth

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_numDaysInYear - return the number of days in a year.
  !
  ! !DESCRIPTION:
  !    Deturn the number of days in a year.
  !
  ! !REVISION HISTORY:
  !     2002-sep-18 - B. Kauffman - initial version.
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  integer function shr_cal_numDaysInYear(year,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in ) :: year      ! calendar year
    character(*)        ,intent(in)  :: calendar  ! calendar type

    !EOP

    type(ESMF_time) :: time1,time2
    type(ESMF_timeInterval) :: timeint
    type(ESMF_CALKIND_FLAG) :: calkind
    integer(SHR_KIND_IN) :: eday
    character(len=shr_cal_calMaxLen) :: lcalendar
    integer :: rc
    character(*),parameter :: subName = "(shr_cal_numDaysInYear)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    lcalendar = shr_cal_calendarName(calendar)

    calkind = ESMF_CALKIND_NOLEAP
    if (trim(lcalendar) == trim(shr_cal_gregorian)) then
       calkind = ESMF_CALKIND_GREGORIAN
    endif

    call ESMF_TimeSet(time1,yy=year,mm=1,dd=1,s=0,calkindflag=calkind,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    call ESMF_TimeSet(time2,yy=year+1,mm=1,dd=1,s=0,calkindflag=calkind,rc=rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
    timeint = time2-time1
    call ESMF_TimeIntervalGet(timeint,StartTimeIn=time1,d=eday)
    shr_cal_numDaysInYear = eday

  end function shr_cal_numDaysInYear

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_elapsDaysStrtMonth - return the number of elapsed days
  !            at start of month
  !
  ! !DESCRIPTION:
  !    Return the number of elapsed days at start of a month.
  !
  ! !REVISION HISTORY:
  !     2002-Oct-29 - R. Jacob - initial version
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  integer function shr_cal_elapsDaysStrtMonth(year,month,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in ) :: year,month  ! calendar year,month
    character(*)        ,intent(in)  :: calendar  ! calendar type

    !EOP

    integer :: k
    character(*),parameter :: subName = "(shr_cal_elapsDaysStrtMonth)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    shr_cal_elapsDaysStrtMonth = 0
    do k = 1,month-1
       shr_cal_elapsDaysStrtMonth = shr_cal_elapsDaysStrtMonth + &
            shr_cal_numDaysInMonth(year,k,calendar)
    enddo

  end function shr_cal_elapsDaysStrtMonth

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_date2ymd - converts coded-date to year/month/day.
  !
  ! !DESCRIPTION:
  !     Converts coded-date (yyyymmdd) to year/month/day.
  !
  ! !REVISION HISTORY:
  !     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  subroutine shr_cal_date2ymd_int (date,year,month,day)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in)  :: date             ! coded-date (yyyymmdd)
    integer(SHR_KIND_IN),intent(out) :: year,month,day   ! calendar year,month,day

    !EOP

    integer(SHR_KIND_IN) :: tdate   ! temporary date
    character(*),parameter :: subName = "(shr_cal_date2ymd)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    if (debug > 1) write(s_logunit,*) trim(subname),'_a ',date

    tdate = abs(date)
    year =int(     tdate       /10000)
    if (date < 0) year = -year
    month=int( mod(tdate,10000)/  100)
    day  =     mod(tdate,  100)

    if (debug > 1) write(s_logunit,*) trim(subname),'_b ',year,month,day

  end subroutine shr_cal_date2ymd_int

  subroutine shr_cal_date2ymd_long (date,year,month,day)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_I8),intent(in)  :: date             ! coded-date ([yy]yyyymmdd)
    integer(SHR_KIND_IN),intent(out) :: year,month,day   ! calendar year,month,day

    !EOP

    integer(SHR_KIND_I8) :: tdate   ! temporary date
    character(*),parameter :: subName = "(shr_cal_date2ymd_long)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    if (debug > 1) write(s_logunit,*) trim(subname),'_a ',date

    tdate = abs(date)
    year =int(     tdate       /10000)
    if (date < 0) year = -year
    month=int( mod(tdate,10000_SHR_KIND_I8)/  100)
    day  =     mod(tdate,  100_SHR_KIND_I8)

    if (debug > 1) write(s_logunit,*) trim(subname),'_b ',year,month,day

  end subroutine shr_cal_date2ymd_long

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_date2julian - converts coded-date to julian day of year
  !
  ! !DESCRIPTION:
  !     Converts coded-date to julian day of year
  !
  ! !REVISION HISTORY:
  !     2009-Oct-23 - T. Craig - initial version
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  subroutine shr_cal_date2julian_int(date,sec,jday,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in ) :: date            ! coded (yyyymmdd) calendar date
    integer(SHR_KIND_IN),intent(in ) :: sec             ! seconds
    real   (SHR_KIND_R8),intent(out) :: jday            ! julian day of year
    character(*)        ,intent(in)  :: calendar        ! calendar type

    !EOP

    !--- local ---
    integer(SHR_KIND_IN) :: year,month,day
    character(*),parameter :: subName = "(shr_cal_date2julian)"

    !-------------------------------------------------------------------------------
    ! NOTE:
    !   julian day of year since yy-01-01-00000
    !-------------------------------------------------------------------------------

    if (debug > 1) write(s_logunit,*) trim(subname),'_a ',date,sec

    call shr_cal_date2ymd(date,year,month,day)
    call shr_cal_ymd2julian(year,month,day,sec,jday,calendar)

    if (debug > 1) write(s_logunit,*) trim(subname),'_b ',jday

  end subroutine shr_cal_date2julian_int

  subroutine shr_cal_date2julian_long(date,sec,jday,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_I8),intent(in ) :: date            ! coded ([yy]yyyymmdd) calendar date
    integer(SHR_KIND_IN),intent(in ) :: sec             ! seconds
    real   (SHR_KIND_R8),intent(out) :: jday            ! julian day of year
    character(*)        ,intent(in)  :: calendar        ! calendar type

    !EOP

    !--- local ---
    integer(SHR_KIND_IN) :: year,month,day
    character(*),parameter :: subName = "(shr_cal_date2julian)"

    !-------------------------------------------------------------------------------
    ! NOTE:
    !   julian day of year since yy-01-01-00000
    !-------------------------------------------------------------------------------

    if (debug > 1) write(s_logunit,*) trim(subname),'_a ',date,sec

    call shr_cal_date2ymd(date,year,month,day)
    call shr_cal_ymd2julian(year,month,day,sec,jday,calendar)

    if (debug > 1) write(s_logunit,*) trim(subname),'_b ',jday

  end subroutine shr_cal_date2julian_long

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_ymd2julian - converts y,m,d,s to julian day of year
  !
  ! !DESCRIPTION:
  !     Converts y,m,d,s to julian day of year
  !
  ! !REVISION HISTORY:
  !     2009-Oct-23 - T. Craig - initial version
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  subroutine shr_cal_ymd2julian(year,month,day,sec,jday,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in ) :: year            ! year
    integer(SHR_KIND_IN),intent(in ) :: month           ! month
    integer(SHR_KIND_IN),intent(in ) :: day             ! day
    integer(SHR_KIND_IN),intent(in ) :: sec             ! seconds
    real   (SHR_KIND_R8),intent(out) :: jday            ! julian day of year
    character(*)        ,intent(in)  :: calendar        ! calendar type

    !EOP

    !--- local ---
    character(*),parameter :: subName = "(shr_cal_ymd2julian)"

    !-------------------------------------------------------------------------------
    ! NOTE:
    !   julian day of year since yy-01-01-000000
    !-------------------------------------------------------------------------------

    if (debug > 1) write(s_logunit,*) trim(subname),'_a ',year,month,day,sec

    if (.not. shr_cal_validYMD(year,month,day,calendar)) then
       write(s_logunit,*) trim(subname),' ERROR: invalid ymd',year,month,day
       call shr_sys_abort(trim(subname)//' ERROR: invalid ymd')
    endif

    jday = shr_cal_elapsDaysStrtMonth(year,month,calendar) + day + sec/SHR_CONST_CDAY

    if (debug > 1) write(s_logunit,*) trim(subname),'_b ',jday

  end subroutine shr_cal_ymd2julian

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_ymd2date - converts year, month, day to coded-date
  !
  ! !DESCRIPTION:
  !     Converts  year, month, day to coded-date
  !
  ! !REVISION HISTORY:
  !     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  subroutine shr_cal_ymd2date_int(year,month,day,date)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in ) :: year,month,day  ! calendar year,month,day
    integer(SHR_KIND_IN),intent(out) :: date            ! coded (yyyymmdd) calendar date

    !EOP

    !--- local ---
    character(*),parameter :: subName = "(shr_cal_ymd2date)"

    !-------------------------------------------------------------------------------
    ! NOTE:
    !   this calendar has a year zero (but no day or month zero)
    !-------------------------------------------------------------------------------

    if (debug > 1) write(s_logunit,*) trim(subname),'_a ',year,month,day

    date = abs(year)*10000 + month*100 + day  ! coded calendar date
    if (year < 0) date = -date

    if (debug > 1) write(s_logunit,*) trim(subname),'_b ',date

  end subroutine shr_cal_ymd2date_int

  subroutine shr_cal_ymd2date_long(year,month,day,date)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in ) :: year,month,day  ! calendar year,month,day
    integer(SHR_KIND_I8),intent(out) :: date            ! coded ([yy]yyyymmdd) calendar date

    !EOP

    !--- local ---
    character(*),parameter :: subName = "(shr_cal_ymd2date)"

    !-------------------------------------------------------------------------------
    ! NOTE:
    !   this calendar has a year zero (but no day or month zero)
    !-------------------------------------------------------------------------------

    if (debug > 1) write(s_logunit,*) trim(subname),'_a ',year,month,day

    date = abs(year)*10000_SHR_KIND_I8 + month*100 + day  ! coded calendar date
    if (year < 0) date = -date

    if (debug > 1) write(s_logunit,*) trim(subname),'_b ',date

  end subroutine shr_cal_ymd2date_long

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_advDate - advances a date and seconds with a delta time
  !
  ! !DESCRIPTION:
  !     Advances a date and seconds with a delta time
  !     Accuracy only good to nearest second
  !
  ! !REVISION HISTORY:
  !     2009-Jun-09 - T. Craig - allows delta < 0
  !     2005-Jun-10 - B. Kauffman - bug fix, simplified algorithm
  !     2005-May-15 - T. Craig - initial version
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  subroutine shr_cal_advDate_int(delta,units,dateIN,secIN,dateOUT,secOUT,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    real   (SHR_KIND_R8) ,intent(in)  :: delta     ! time increment
    character(*)         ,intent(in)  :: units     ! units of increment
    integer(SHR_KIND_IN) ,intent(in)  :: dateIN    ! base date, yyyymmdd
    real   (SHR_KIND_R8) ,intent(in)  :: secIN     ! base seconds
    integer(SHR_KIND_IN) ,intent(out) :: dateOUT   ! new date, yyyymmdd
    real   (SHR_KIND_R8) ,intent(out) :: secOUT    ! new seconds
    character(*)         ,intent(in)  :: calendar  ! calendar type

    !EOP

    !--- local ---
    type(ESMF_time) :: timeIn, timeOut
    type(ESMF_timeInterval) :: dt
    real   (SHR_KIND_R8)   :: dSec    ! delta-sec: advance date this many seconds
    integer(SHR_KIND_I8)   :: i8dsec, i8dday ! delta sec and day in i8
    integer(SHR_KIND_I8)   :: spd     ! seconds per day in i8
    integer(SHR_KIND_I4)   :: idday, idsec   ! delta sec and dat in i4
    integer(SHR_KIND_I4)   :: year, month, day, sec  ! calendar stuff
    character(len=shr_cal_calMaxLen) :: lcalendar
    type(ESMF_CALKIND_FLAG) :: calkind

    !--- formats ---
    character(*),parameter :: subName = "(shr_cal_advDate)"
    character(*),parameter :: F00 = "('(shr_cal_advDate) ',a,i5)"
    character(*),parameter :: F02 = "('(shr_cal_advDate) ',a,i8.8,f10.3)"

    !-------------------------------------------------------------------------------
    ! NOTE:
    !-------------------------------------------------------------------------------
    dSec = 0.0_SHR_KIND_R8
    !--- calculate delta-time in seconds ---
    if     (trim(units) == 'days'   ) then
       dSec = delta * SHR_CONST_CDAY
    elseif (trim(units) == 'hours'  ) then
       dSec = delta * 3600.0_SHR_KIND_R8
    elseif (trim(units) == 'minutes') then
       dSec = delta *   60.0_SHR_KIND_R8
    elseif (trim(units) == 'seconds') then
       dSec = delta *    1.0_SHR_KIND_R8
    else
       call shr_sys_abort(trim(subname)//' ERROR: unrecognized time units '//trim(units))
    endif

    ! take secIn into account here since it's real
    dSec = dSec + secIn

    ! i8 math, convert reals to nearest second
    i8dSec = nint(dSec,SHR_KIND_I8)
    spd = nint(SHR_CONST_CDAY)
    i8dday = i8dsec/spd
    i8dsec = i8dsec - i8dday*spd

    ! convert to i4
    idday = i8dday
    idsec = i8dsec

    calkind = ESMF_CALKIND_NOLEAP
    lcalendar = shr_cal_calendarName(calendar)

    if (trim(lcalendar) == trim(shr_cal_gregorian)) then
       calkind = ESMF_CALKIND_GREGORIAN
    endif

    call shr_cal_date2ymd(dateIn,year,month,day)
    call ESMF_TimeSet(timeIN,yy=year,mm=month,dd=day,calkindflag=calkind)
    call ESMF_TimeIntervalSet(dt,d=idday,s=idsec)

    timeOut = timeIn + dt

    call ESMF_TimeGet(timeOut,yy=year,mm=month,dd=day,s=sec)
    call shr_cal_ymd2date(year,month,day,dateOut)
    secOut = sec

    if (debug>0) then
       if (s_loglev > 0) write(s_logunit,*) subName," units,delta,calendar=",trim(units),delta,' ',trim(calendar)
       if (s_loglev > 0) write(s_logunit,F02) "dateIN ,secIN =",dateIN ,secIN
       if (s_loglev > 0) write(s_logunit,F02) "dateOUT,secOUT=",dateOUT,secOUT
    end if

  end subroutine shr_cal_advDate_int

  subroutine shr_cal_advDate_long(delta,units,dateIN,secIN,dateOUT,secOUT,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    real   (SHR_KIND_R8) ,intent(in)  :: delta     ! time increment
    character(*)         ,intent(in)  :: units     ! units of increment
    integer(SHR_KIND_I8) ,intent(in)  :: dateIN    ! base date, [yy]yyyymmdd
    real   (SHR_KIND_R8) ,intent(in)  :: secIN     ! base seconds
    integer(SHR_KIND_I8) ,intent(out) :: dateOUT   ! new date, [yy]yyyymmdd
    real   (SHR_KIND_R8) ,intent(out) :: secOUT    ! new seconds
    character(*)         ,intent(in)  :: calendar  ! calendar type

    !EOP

    !--- local ---
    type(ESMF_time) :: timeIn, timeOut
    type(ESMF_timeInterval) :: dt
    real   (SHR_KIND_R8)   :: dSec    ! delta-sec: advance date this many seconds
    integer(SHR_KIND_I8)   :: i8dsec, i8dday ! delta sec and day in i8
    integer(SHR_KIND_I8)   :: spd     ! seconds per day in i8
    integer(SHR_KIND_I4)   :: idday, idsec   ! delta sec and dat in i4
    integer(SHR_KIND_I4)   :: year, month, day, sec  ! calendar stuff
    character(len=shr_cal_calMaxLen) :: lcalendar
    type(ESMF_CALKIND_FLAG) :: calkind

    !--- formats ---
    character(*),parameter :: subName = "(shr_cal_advDate)"
    character(*),parameter :: F00 = "('(shr_cal_advDate) ',a,i5)"
    character(*),parameter :: F02 = "('(shr_cal_advDate) ',a,i10.8,f10.3)"

    !-------------------------------------------------------------------------------
    ! NOTE:
    !-------------------------------------------------------------------------------
    dSec = 0.0_SHR_KIND_R8
    !--- calculate delta-time in seconds ---
    if     (trim(units) == 'days'   ) then
       dSec = delta * SHR_CONST_CDAY
    elseif (trim(units) == 'hours'  ) then
       dSec = delta * 3600.0_SHR_KIND_R8
    elseif (trim(units) == 'minutes') then
       dSec = delta *   60.0_SHR_KIND_R8
    elseif (trim(units) == 'seconds') then
       dSec = delta *    1.0_SHR_KIND_R8
    else
       call shr_sys_abort(trim(subname)//' ERROR: unrecognized time units '//trim(units))
    endif

    ! take secIn into account here since it's real
    dSec = dSec + secIn

    ! i8 math, convert reals to nearest second
    i8dSec = nint(dSec,SHR_KIND_I8)
    spd = nint(SHR_CONST_CDAY)
    i8dday = i8dsec/spd
    i8dsec = i8dsec - i8dday*spd

    ! convert to i4
    idday = i8dday
    idsec = i8dsec

    calkind = ESMF_CALKIND_NOLEAP
    lcalendar = shr_cal_calendarName(calendar)

    if (trim(lcalendar) == trim(shr_cal_gregorian)) then
       calkind = ESMF_CALKIND_GREGORIAN
    endif

    call shr_cal_date2ymd(dateIn,year,month,day)
    call ESMF_TimeSet(timeIN,yy=year,mm=month,dd=day,calkindflag=calkind)
    call ESMF_TimeIntervalSet(dt,d=idday,s=idsec)

    timeOut = timeIn + dt

    call ESMF_TimeGet(timeOut,yy=year,mm=month,dd=day,s=sec)
    call shr_cal_ymd2date(year,month,day,dateOut)
    secOut = sec

    if (debug>0) then
       if (s_loglev > 0) write(s_logunit,*) subName," units,delta,calendar=",trim(units),delta,' ',trim(calendar)
       if (s_loglev > 0) write(s_logunit,F02) "dateIN ,secIN =",dateIN ,secIN
       if (s_loglev > 0) write(s_logunit,F02) "dateOUT,secOUT=",dateOUT,secOUT
    end if

  end subroutine shr_cal_advDate_long

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_advDateInt - advances a date and seconds with a delta time
  !
  ! !DESCRIPTION:
  !     Advances a date and seconds with a delta time
  !
  ! !REVISION HISTORY:
  !     2009-???-?? - ?? - replicated from shr_cal_advDate()
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  subroutine shr_cal_advDateInt_int(delta,units,dateIN,secIN,dateOUT,secOUT,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN) ,intent(in)  :: delta     ! time increment
    character(*)         ,intent(in)  :: units     ! units of increment
    integer(SHR_KIND_IN) ,intent(in)  :: dateIN    ! base date, yyyymmdd
    integer(SHR_KIND_IN) ,intent(in)  :: secIN     ! base seconds
    integer(SHR_KIND_IN) ,intent(out) :: dateOUT   ! new date, yyyymmdd
    integer(SHR_KIND_IN) ,intent(out) :: secOUT    ! new seconds
    character(*)         ,intent(in)  :: calendar  ! calendar type

    !EOP

    !--- local ---
    type(ESMF_time) :: timeIn, timeOut
    type(ESMF_timeInterval) :: dt
    integer(SHR_KIND_IN) :: year,month,day
    character(len=shr_cal_calMaxLen) :: lcalendar
    type(ESMF_CALKIND_FLAG) :: calkind

    !--- formats ---
    character(*),parameter :: subName = "(shr_cal_advDateInt)"
    character(*),parameter :: F00 = "('(shr_cal_advDateInt) ',a,i5)"
    character(*),parameter :: F02 = "('(shr_cal_advDateInt) ',a,i8.8,f10.3)"

    !-------------------------------------------------------------------------------
    ! NOTE:
    !-------------------------------------------------------------------------------
    lcalendar = shr_cal_calendarName(calendar)
    calkind = ESMF_CALKIND_NOLEAP
    if (trim(lcalendar) == trim(shr_cal_gregorian)) then
       calkind = ESMF_CALKIND_GREGORIAN
    endif

    call shr_cal_date2ymd(dateIn,year,month,day)
    call ESMF_TimeSet(timeIN,yy=year,mm=month,dd=day,s=secIn,calkindflag=calkind)

    !--- calculate delta-time in seconds ---
    if     (trim(units) == 'days'   ) then
       call ESMF_TimeIntervalSet(dt,d=delta)
    elseif (trim(units) == 'hours'  ) then
       call ESMF_TimeIntervalSet(dt,h=delta)
    elseif (trim(units) == 'minutes') then
       call ESMF_TimeIntervalSet(dt,m=delta)
    elseif (trim(units) == 'seconds') then
       call ESMF_TimeIntervalSet(dt,s=delta)
    else
       call shr_sys_abort(trim(subname)//' ERROR: unrecognized time units '//trim(units))
    endif

    timeOut = timeIn + dt

    call ESMF_TimeGet(timeOut,yy=year,mm=month,dd=day,s=secOut)
    call shr_cal_ymd2date(year,month,day,dateOut)

    if (debug>0) then
       if (s_loglev > 0) write(s_logunit,*) subName," units,delta,calendar=",trim(units),delta,' ',trim(calendar)
       if (s_loglev > 0) write(s_logunit,F02) "dateIN ,secIN =",dateIN ,secIN
       if (s_loglev > 0) write(s_logunit,F02) "dateOUT,secOUT=",dateOUT,secOUT
    end if

  end subroutine shr_cal_advDateInt_int

  subroutine shr_cal_advDateInt_long(delta,units,dateIN,secIN,dateOUT,secOUT,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN) ,intent(in)  :: delta     ! time increment
    character(*)         ,intent(in)  :: units     ! units of increment
    integer(SHR_KIND_I8) ,intent(in)  :: dateIN    ! base date, yyyymmdd
    integer(SHR_KIND_IN) ,intent(in)  :: secIN     ! base seconds
    integer(SHR_KIND_I8) ,intent(out) :: dateOUT   ! new date, yyyymmdd
    integer(SHR_KIND_IN) ,intent(out) :: secOUT    ! new seconds
    character(*)         ,intent(in)  :: calendar  ! calendar type

    !EOP

    !--- local ---
    type(ESMF_time) :: timeIn, timeOut
    type(ESMF_timeInterval) :: dt
    integer(SHR_KIND_IN) :: year,month,day
    character(len=shr_cal_calMaxLen) :: lcalendar
    type(ESMF_CALKIND_FLAG) :: calkind

    !--- formats ---
    character(*),parameter :: subName = "(shr_cal_advDateInt)"
    character(*),parameter :: F00 = "('(shr_cal_advDateInt) ',a,i5)"
    character(*),parameter :: F02 = "('(shr_cal_advDateInt) ',a,i10.8,f10.3)"

    !-------------------------------------------------------------------------------
    ! NOTE:
    !-------------------------------------------------------------------------------
    lcalendar = shr_cal_calendarName(calendar)
    calkind = ESMF_CALKIND_NOLEAP
    if (trim(lcalendar) == trim(shr_cal_gregorian)) then
       calkind = ESMF_CALKIND_GREGORIAN
    endif

    call shr_cal_date2ymd(dateIn,year,month,day)
    call ESMF_TimeSet(timeIN,yy=year,mm=month,dd=day,s=secIn,calkindflag=calkind)

    !--- calculate delta-time in seconds ---
    if     (trim(units) == 'days'   ) then
       call ESMF_TimeIntervalSet(dt,d=delta)
    elseif (trim(units) == 'hours'  ) then
       call ESMF_TimeIntervalSet(dt,h=delta)
    elseif (trim(units) == 'minutes') then
       call ESMF_TimeIntervalSet(dt,m=delta)
    elseif (trim(units) == 'seconds') then
       call ESMF_TimeIntervalSet(dt,s=delta)
    else
       call shr_sys_abort(trim(subname)//' ERROR: unrecognized time units '//trim(units))
    endif

    timeOut = timeIn + dt

    call ESMF_TimeGet(timeOut,yy=year,mm=month,dd=day,s=secOut)
    call shr_cal_ymd2date(year,month,day,dateOut)

    if (debug>0) then
       if (s_loglev > 0) write(s_logunit,*) subName," units,delta,calendar=",trim(units),delta,' ',trim(calendar)
       if (s_loglev > 0) write(s_logunit,F02) "dateIN ,secIN =",dateIN ,secIN
       if (s_loglev > 0) write(s_logunit,F02) "dateOUT,secOUT=",dateOUT,secOUT
    end if

  end subroutine shr_cal_advDateInt_long

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_validDate - determines if coded-date is a valid date
  !
  ! !DESCRIPTION:
  !    Determines if the given coded-date is a valid date.
  !
  ! !REVISION HISTORY:
  !     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  logical function shr_cal_validDate_int(date,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in ) :: date      ! coded (yyyymmdd) calendar date
    character(len=*)    ,intent(in ) :: calendar  ! calendar

    !EOP

    !--- local ---
    integer(SHR_KIND_IN) :: year,month,day
    character(*),parameter :: subName = "(shr_cal_validDate)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call shr_cal_date2ymd(date,year,month,day)
    shr_cal_validDate_int = shr_cal_validYMD(year,month,day,calendar)

  end function shr_cal_validDate_int


  logical function shr_cal_validDate_long(date,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_I8),intent(in ) :: date      ! coded ([yy]yyyymmdd) calendar date
    character(len=*)    ,intent(in ) :: calendar  ! calendar

    !EOP

    !--- local ---
    integer(SHR_KIND_IN) :: year,month,day

    character(*),parameter :: subName = "(shr_cal_validDate)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    call shr_cal_date2ymd(date,year,month,day)
    shr_cal_validDate_long = shr_cal_validYMD(year,month,day,calendar)

  end function shr_cal_validDate_long

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_validYMD - determines if year, month, day is a valid date
  !
  ! !DESCRIPTION:
  !    Determines if the given year, month, and day indicate a valid date.
  !
  ! !REVISION HISTORY:
  !     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  logical function shr_cal_validYMD(year,month,day,calendar)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in ) :: year,month,day  ! year,month,day
    character(len=*)    ,intent(in ) :: calendar        ! calendar

    !EOP

    !--- local ---

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    shr_cal_validYMD = .true.
    if (year  < -999) shr_cal_validYMD = .false.
    if (year  > 999999) shr_cal_validYMD = .false.
    if (month <    1) shr_cal_validYMD = .false.
    if (month >   12) shr_cal_validYMD = .false.
    if (day   <    1) shr_cal_validYMD = .false.
    if (day   > shr_cal_numDaysInMonth(year,month,calendar)) &
         shr_cal_validYMD = .false.

  end function shr_cal_validYMD

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_validHMS - determines if hour, min, sec is valid
  !
  ! !DESCRIPTION:
  !    Determines if the given hour, min, sec is valid
  !
  ! !REVISION HISTORY:
  !     2005-May-15 - T. Craig - initial version
  !
  ! !INTERFACE:  -----------------------------------------------------------------

  logical function shr_cal_validHMS(hr,min,sec)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer(SHR_KIND_IN),intent(in ) :: hr, min, sec   ! hour, minute, second

    !EOP

    !--- local ---
    character(*),parameter :: subName = "(shr_cal_validHMS)"

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    shr_cal_validHMS = .true.
    if (hr  <    0) shr_cal_validHMS = .false.
    if (hr  >   23) shr_cal_validHMS = .false.
    if (min <    0) shr_cal_validHMS = .false.
    if (min >   59) shr_cal_validHMS = .false.
    if (sec <    0) shr_cal_validHMS = .false.
    if (sec >   60) shr_cal_validHMS = .false.

  end function shr_cal_validHMS

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_setDebug -- Set local shr_cal debug level
  !
  ! !DESCRIPTION:
  !    Set local shr\_cal debug level, 0 = production
  !    \newline
  !    General Usage: call shr\_cal\_setDebug(2)
  !
  ! !REVISION HISTORY:
  !     2005-May-10  - B. Kauffman - first version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_cal_setDebug(level)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer,intent(in) :: level

    !EOP

    !--- formats ---
    character(*),parameter :: subName =   "(shr_cal_setDebug) "
    character(*),parameter :: F00     = "('(shr_cal_setDebug) ',a) "
    character(*),parameter :: F01     = "('(shr_cal_setDebug) ',a,i4) "

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    debug = level
    if (s_loglev > 0) write(s_logunit,F01) "debug level reset to ",level

  end subroutine shr_cal_setDebug

  !===============================================================================
  !BOP ===========================================================================
  !
  ! !IROUTINE: shr_cal_getDebug -- get shr_cal internal debug level
  !
  ! !DESCRIPTION:
  !    Get shr_cal internal debug level, 0 = production
  !    \newline
  !    General Usage: call shr\_cal\_getDebug(level)
  !
  ! !REVISION HISTORY:
  !     2005-May-10  - B. Kauffman - first version
  !
  ! !INTERFACE: ------------------------------------------------------------------

  subroutine shr_cal_getDebug(level)

    implicit none

    ! !INPUT/OUTPUT PARAMETERS:

    integer,intent(out) :: level

    !EOP

    !--- formats ---
    character(*),parameter :: subName =   "(shr_cal_getDebug) "
    character(*),parameter :: F00     = "('(shr_cal_getDebug) ',a) "

    !-------------------------------------------------------------------------------
    !
    !-------------------------------------------------------------------------------

    level = debug

  end subroutine shr_cal_getDebug

  !===============================================================================
  subroutine shr_cal_datetod2string_int(date_str, ymd, tod)
    ! Converts coded date (yyyymmdd) and optional time of day to a string like
    ! 'yyyy-mm-dd-ttttt' (if tod is present) or 'yyyy-mm-dd' (if tod is absent).
    !
    ! yyyy in the output string will have at least 4 but no more than 6 characters (with
    ! leading zeroes if necessary).
    character(len=*), intent(out) :: date_str
    integer(shr_kind_in), intent(in) :: ymd
    integer(shr_kind_in), intent(in), optional :: tod

    integer(shr_kind_in) :: yy, mm, dd

    call shr_cal_date2ymd(ymd, yy, mm, dd)
    call shr_cal_ymdtod2string(date_str, yy, mm, dd, tod)
  end subroutine shr_cal_datetod2string_int

  subroutine shr_cal_datetod2string_long(date_str, ymd, tod)
    ! Converts coded date (yyyymmdd) and optional time of day to a string like
    ! 'yyyy-mm-dd-ttttt' (if tod is present) or 'yyyy-mm-dd' (if tod is absent).
    !
    ! yyyy in the output string will have at least 4 but no more than 6 characters (with
    ! leading zeroes if necessary).
    character(len=*), intent(out) :: date_str
    integer(shr_kind_i8), intent(in) :: ymd
    integer(shr_kind_in), intent(in), optional :: tod

    integer(shr_kind_in) :: yy, mm, dd

    call shr_cal_date2ymd(ymd, yy, mm, dd)
    call shr_cal_ymdtod2string(date_str, yy, mm, dd, tod)
  end subroutine shr_cal_datetod2string_long

  subroutine shr_cal_ymdtod2string(date_str, yy, mm, dd, tod)
    ! Converts year, month, day and time of day to a string like 'yyyy-mm-dd-ttttt'.
    !
    ! The only required input is yy (year). Missing inputs will result in missing pieces in
    ! the output string. However, if tod is present, then mm and dd must be present; if dd
    ! is present, then mm must be present.
    !
    ! yyyy in the output string will have at least 4 but no more than 6 characters (with
    ! leading zeroes if necessary).
    integer(shr_kind_in), intent(in) :: yy
    integer(shr_kind_in), intent(in), optional :: mm, dd, tod
    character(len=*), intent(out) :: date_str

    character(len=6) :: year_str
    character(len=3) :: month_str
    character(len=3) :: day_str
    character(len=6) :: time_str
    character(len=*), parameter :: subname = 'shr_cal_ymdtod2string'

    if (yy > 999999) then
       write(s_logunit,*) trim(subname),' : ERROR: year too large: ', yy
       write(s_logunit,*) '(Currently, only years up to 999999 are supported)'
       call shr_sys_abort(trim(subname)//' : year too large (max of 999999)')
    end if
    write(year_str,'(i6.4)') yy
    year_str = adjustl(year_str)

    if (present(mm)) then
       write(month_str,'(a,i2.2)') '-',mm
    else
       month_str = ' '
    end if

    if (present(dd)) then
       if (.not. present(mm)) then
          call shr_sys_abort(trim(subname)//' : if dd is present, then mm must be present, too')
       end if
       write(day_str,'(a,i2.2)') '-',dd
    else
       day_str = ' '
    end if

    if (present(tod)) then
       if (.not. present(mm) .or. .not. present(dd)) then
          call shr_sys_abort(trim(subname)//' : if tod is present, then mm and dd must be present, too')
       end if
       write(time_str,'(a,i5.5)') '-',tod
    else
       time_str = ' '
    end if

    if (len_trim(year_str) + len_trim(month_str) + len_trim(day_str) + len_trim(time_str) > len(date_str)) then
       call shr_sys_abort(trim(subname//' : output string too short'))
    else
       date_str = trim(year_str) // trim(month_str) // trim(day_str) // trim(time_str)
    end if

  end subroutine shr_cal_ymdtod2string

  !===============================================================================
  subroutine shr_cal_ymds2rday_offset(etime, rdays_offset, &
       years_offset, months_offset, days_offset, seconds_offset)
    ! Given the current time (etime) and optional year, month, day and seconds offsets
    ! from the current time: Return an offset from the current time given in fractional
    ! days.
    !
    ! For example, if day_offset = -2 and seconds_offset = -21600, then rday_offset will
    ! be -2.25.

    ! !INPUT/OUTPUT PARAMETERS:

    type(ESMF_Time), intent(in) :: etime  ! current time
    real(SHR_KIND_R8), intent(out) :: rdays_offset  ! offset from current time in fractional days

    ! One or more of the following optional arguments should be provided:
    integer(SHR_KIND_IN), intent(in), optional :: years_offset   ! number of years offset from current time
    integer(SHR_KIND_IN), intent(in), optional :: months_offset  ! number of months offset from current time
    integer(SHR_KIND_IN), intent(in), optional :: days_offset    ! number of days offset from current time
    integer(SHR_KIND_IN), intent(in), optional :: seconds_offset ! number of seconds offset from current time

    !--- local ---
    type(ESMF_TimeInterval) :: timeinterval
    integer :: rc

    !-------------------------------------------------------------------------------

    call ESMF_TimeIntervalSet(timeinterval = timeinterval, &
         startTime = etime, &
         YY = years_offset, &
         MM = months_offset, &
         D  = days_offset, &
         S  = seconds_offset, &
         rc = rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)

    call ESMF_TimeIntervalGet(timeinterval = timeinterval, d_r8 = rdays_offset, rc = rc)
    if(rc /= ESMF_SUCCESS) call ESMF_Finalize(rc=rc, endflag=ESMF_END_ABORT)
  end subroutine shr_cal_ymds2rday_offset

  !===============================================================================
end module shr_cal_mod
