!===============================================================================
! SVN $Id: shr_cal_mod.F90 239 2006-02-08 19:02:33Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070423/shr/shr_cal_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
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
!
! !REMARKS:
!   Following are some internal assumptions.  These assumptions are somewhat
!   arbitrary -- they were chosen because they result in the simplest code given
!   the requirements of this module.  These assumptions can be relaxed as 
!   necessary: 
!   o the valid range of years is [0,9999]
!   o elapsed days = 0 <=> January 1st, year 0000
!   o all years have 365 days (no leap years)
!     This module is hard-coded to implement a 365-day calendar, ie. there
!     are no leap years.  This module can be modified to implement a calendar
!     with leap years if this becomes desireable.  This would make the internal
!     logic of this module more complex, but would not require any change to the
!     module API or the calling code because the module API hides these details
!     from all external routines.
!
! !INTERFACE: ------------------------------------------------------------------

module shr_cal_mod

! !USES:

   use shr_kind_mod   ! kinds
   use shr_const_mod  ! constants
   use shr_sys_mod    ! system

   implicit none

   private ! except

! !PUBLIC TYPES: 
 
   ! none

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_cal_set        ! converts elapsed days to coded-date
   public :: shr_cal_get        ! converts elapsed days to coded-date
   public :: shr_cal_eday2date  ! converts elapsed days to coded-date
   public :: shr_cal_eday2ymd   ! converts elapsed days to yr,month,day
   public :: shr_cal_date2ymd   ! converts coded-date   to yr,month,day
   public :: shr_cal_date2eday  ! converts coded-date   to elapsed days
   public :: shr_cal_ymd2date   ! converts yr,month,day to coded-date
   public :: shr_cal_ymd2eday   ! converts yr,month,day to elapsed days
   public :: shr_cal_advDate    ! advance date/secs
   public :: shr_cal_validDate  ! logical function: is coded-date valid?
   public :: shr_cal_validYMD   ! logical function: are yr,month,day valid?
   public :: shr_cal_validHMS   ! logical function: are hr, min, sec valid?
   public :: shr_cal_numDaysinMonth     ! number of days in a month
   public :: shr_cal_elapsDaysStrtMonth ! elapsed days on start of month
   public :: shr_cal_getDebug   ! get internal debug level
   public :: shr_cal_setDebug   ! set internal debug level

! !PUBLIC DATA MEMBERS:

   ! none

!EOP

   !----- local -----
   character(SHR_KIND_CL),save    :: calendar_type='noleap' ! calendar type
   integer(SHR_KIND_IN),parameter :: shr_cal_nvalidTypes = 2
   character(*),parameter :: &
     shr_cal_validTypes(shr_cal_nvalidTypes) = (/'noleap    ', &
                                                 '365_day   '/)

   !--- this is the noleap calendar ---
   integer(SHR_KIND_IN),parameter :: dsm(12) = &     ! elapsed Days on Start of Month
   &                    (/ 0,31,59, 90,120,151, 181,212,243, 273,304,334/)
   integer(SHR_KIND_IN),parameter :: dpm(12) = &     ! Days Per Month
   &                    (/31,28,31, 30, 31, 30,  31, 31, 30,  31, 30, 31/)


   !--- trigger internal debug output ---
   integer(SHR_KIND_IN) :: debug = 0

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_set - set the calenday type
!
! !DESCRIPTION:
!     Set the calendar type
!
! !REVISION HISTORY:
!     2005-May-28 - T. Craig - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_set(ctype)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(in)  :: ctype   ! calendar type

!EOP

   !--- local ---
   integer(SHR_KIND_IN)  :: n   ! counter
   logical :: found             ! check for valid type

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   found = .false.
   do n = 1,shr_cal_nvalidTypes
     if (trim(ctype) == trim(shr_cal_validTypes(n))) then
       calendar_type = trim(ctype)
       found = .true.
     endif
   enddo

   if (.not.found) call shr_sys_abort('shr_cal_set ERROR illegal calendar type')

end subroutine shr_cal_set

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_get - get the calenday type
!
! !DESCRIPTION:
!     Get the calendar type
!
! !REVISION HISTORY:
!     2005-May-28 - T. Craig - initial version
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_get(ctype)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   character(*),intent(out)  :: ctype   ! calendar type

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   ctype = calendar_type

end subroutine shr_cal_get

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_eday2date - converts elapsed days to coded-date
!
! !DESCRIPTION:
!     Converts elapsed days to coded-date.
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_eday2date(eday,date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in)  :: eday  ! number of elapsed days
   integer(SHR_KIND_IN),intent(out) :: date  ! coded (yyyymmdd) calendar date

!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: k,year,month,day

!-------------------------------------------------------------------------------
! ASSUMPTIONS:
!   this calendar has a year zero (but no day or month zero)
!-------------------------------------------------------------------------------

   year = eday/365       ! calendar year (note: Fortran truncation)
   day  = mod(eday,365)  ! elapsed days within current year
   do k=1,12
     IF (day >= dsm(k)) month=k   ! calendar month
   end do
   day = day-dsm(month) + 1         ! calendar day
  
   date = year*10000 + month*100 + day  ! coded calendar date

end subroutine shr_cal_eday2date

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_eday2ymd - converts elapsed days to year/month/day.
!
! !DESCRIPTION:
!     Converts elapsed days to year/month/day.
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_eday2ymd (eday,year,month,day)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in)  :: eday             ! elapsed days
   integer(SHR_KIND_IN),intent(out) :: year,month,day   ! calendar year,month,day

!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: k

!-------------------------------------------------------------------------------
! ASSUMPTIONS:
!   this calendar has a year zero (but no day or month zero)
!-------------------------------------------------------------------------------

   year = eday/365       ! calendar year (note: Fortran truncation)
   day  = mod(eday,365)  ! elapsed days within current year
   do k=1,12
     IF (day .ge. dsm(k)) month=k   ! calendar month
   end do
   day = day-dsm(month) + 1         ! calendar day

end subroutine shr_cal_eday2ymd 

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

subroutine shr_cal_date2ymd (date,year,month,day)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in)  :: date             ! coded-date (yyyymmdd)
   integer(SHR_KIND_IN),intent(out) :: year,month,day   ! calendar year,month,day

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (.not. shr_cal_validDate(date)) then
      write(6,*) "(cal_date2ymd) ERROR: invalid date = ",date
   endif

   year =int(     date       /10000)
   month=int( mod(date,10000)/  100)
   day  =     mod(date,  100) 

end subroutine shr_cal_date2ymd 

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_date2eday - converts coded-date to elapsed days
!
! !DESCRIPTION:
!     Converts coded-date to elapsed days
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_date2eday(date,eday)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: date            ! coded (yyyymmdd) calendar date
   integer(SHR_KIND_IN),intent(out) :: eday            ! number of elapsed days

!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: year,month,day

!-------------------------------------------------------------------------------
! NOTE:
!   elapsed days since yy-mm-dd = 00-01-01, with 0 elapsed seconds
!-------------------------------------------------------------------------------

   if (.not. shr_cal_validDate(date)) stop 

   year =int(     date       /10000)
   month=int( mod(date,10000)/  100)
   day  =     mod(date,  100) 

   eday = year*365 + dsm(month) + (day-1)

end subroutine shr_cal_date2eday

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

subroutine shr_cal_ymd2date(year,month,day,date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year,month,day  ! calendar year,month,day
   integer(SHR_KIND_IN),intent(out) :: date            ! coded (yyyymmdd) calendar date

!EOP

   !--- local ---

!-------------------------------------------------------------------------------
! NOTE:
!   this calendar has a year zero (but no day or month zero)
!-------------------------------------------------------------------------------

   if (.not. shr_cal_validYMD(year,month,day)) stop 

   date = year*10000 + month*100 + day  ! coded calendar date

end subroutine shr_cal_ymd2date

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_ymd2eday - converts year, month, day to elapsed days
!
! !DESCRIPTION:
!     Converts  year, month, day to elapsed days
!
! !REVISION HISTORY:
!     2001-dec-28 - B. Kauffman - initial version, taken from cpl5
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine  shr_cal_ymd2eday(year,month,day,eday)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year,month,day  ! calendar year,month,day
   integer(SHR_KIND_IN),intent(out) :: eday            ! number of elapsed days

!EOP

   !--- local ---

!-------------------------------------------------------------------------------
! NOTE:
!   elapsed days since yy-mm-dd = 00-01-01, with 0 elapsed seconds
!-------------------------------------------------------------------------------

   if (.not. shr_cal_validYMD(year,month,day)) stop 

   eday = year*365 + dsm(month) + (day-1)

end subroutine  shr_cal_ymd2eday

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_cal_advDate - advances a date and seconds with a delta time
!
! !DESCRIPTION:
!     Advances a date and seconds with a delta time
!
! !REVISION HISTORY:
!     2005-Jun-10 - B. Kauffman - bug fix, simplified algorithm
!     2005-May-15 - T. Craig - initial version 
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_cal_advDate(delta,units,dateIN,secIN,dateOUT,secOUT,calendar)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   real   (SHR_KIND_R8) ,intent(in)  :: delta     ! time increment
   character(*)         ,intent(in)  :: units     ! units of increment
   integer(SHR_KIND_IN) ,intent(in)  :: dateIN    ! base date, yyyymmdd
   real   (SHR_KIND_R8) ,intent(in)  :: secIN     ! base seconds
   integer(SHR_KIND_IN) ,intent(out) :: dateOUT   ! new date, yyyymmdd
   real   (SHR_KIND_R8) ,intent(out) :: secOUT    ! new seconds
   character(*),optional,intent(in)  :: calendar  ! calendar type

!EOP

   !--- local ---
   real   (SHR_KIND_R8)   :: dSec    ! delta-sec: advance date this many seconds
   integer(SHR_KIND_IN)   :: dDay    ! advance this many whole days
   real   (SHR_KIND_R8)   :: rSec    ! advance this many "remainder seconds"
   integer(SHR_KIND_IN)   :: eDay    ! date as elapsed dates from ref date
   character(SHR_KIND_CL) :: calOrig ! original calendar type

   !--- formats ---
   character(*),parameter :: subName = "(shr_cal_advDate)"
   character(*),parameter :: F00 = "('(shr_cal_advDate) ',a,i5)"
   character(*),parameter :: F02 = "('(shr_cal_advDate) ',a,i8.8,f10.3)"
   
!-------------------------------------------------------------------------------
! NOTE:
!-------------------------------------------------------------------------------

   call shr_cal_get(calOrig)

   !--- allow the temporary use of an alternate calendar ---
   if (present(calendar)) call shr_cal_set(calendar)

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
      call shr_sys_abort(' ERROR: unrecognized time units '//trim(units))
   endif

   !--- starting from dateIN and zero seconds...      ---
   !--- advance (secIN + dSec) seconds to arrive at output date ---

   call shr_cal_date2eday(dateIN,eDay) ! starting from this eDay...
   dSec = dSec + secIN                 ! advance this many seconds

   dDay = int(dSec/SHR_CONST_CDAY)     ! advance this many whole days...
   rSec = mod(dSec,SHR_CONST_CDAY)     ! ...plus this many secs (less than a day)

   call shr_cal_eday2date(eDay+dDay,dateOUT)
   secOUT = rSec    

   if (debug>0) then
      if (present(calendar)) then
	 write(6,*) subName," units,delta,calendar=",trim(units),' ',trim(calendar),delta
      else
	 write(6,*) subName," units,delta=",trim(units),delta
      endif
      write(6,F02) "dateIN ,secIN =",dateIN ,secIN
      write(6,F02) "dateOUT,secOUT=",dateOUT,secOUT
   end if

   call shr_cal_set(calOrig)

end subroutine shr_cal_advDate

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

logical function shr_cal_validDate(date) 

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: date            ! coded (yyyymmdd) calendar date

!EOP

   !--- local ---
   integer(SHR_KIND_IN) :: year,month,day

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   year =int(     date       /10000)
   month=int( mod(date,10000)/  100)
   day  =     mod(date,  100) 

   shr_cal_validDate = .true.
   if (year  <    0) shr_cal_validDate = .false.
   if (year  > 9999) shr_cal_validDate = .false.
   if (month <    1) shr_cal_validDate = .false.
   if (month >   12) shr_cal_validDate = .false.
   if (day   <    1) shr_cal_validDate = .false.
   if (shr_cal_validDate) then
      if (day > dpm(month)) shr_cal_validDate = .false.
   endif

end function shr_cal_validDate

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

logical function shr_cal_validYMD(year,month,day)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year,month,day  ! calendar year,month,day

!EOP

   !--- local ---

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   shr_cal_validYMD = .true.
   if (year  <    0) shr_cal_validYMD = .false.
   if (year  > 9999) shr_cal_validYMD = .false.
   if (month <    1) shr_cal_validYMD = .false.
   if (month >   12) shr_cal_validYMD = .false.
   if (day   <    1) shr_cal_validYMD = .false.
   if (shr_cal_validYMD) then
      if (day > dpm(month)) shr_cal_validYMD = .false.
   endif

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
! !IROUTINE: shr_cal_numDaysInMonth - return the number of days in a month.
!
! !DESCRIPTION:
!    Deturn the number of days in a month.
!
! !REVISION HISTORY:
!     2002-sep-18 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

integer function shr_cal_numDaysInMonth(year,month)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year,month  ! calendar year,month

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   shr_cal_numDaysInMonth = dpm(month)

end function shr_cal_numDaysInMonth

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

integer function shr_cal_elapsDaysStrtMonth(year,month)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in ) :: year,month  ! calendar year,month

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   shr_cal_elapsDaysStrtMonth = dsm(month)

end function shr_cal_elapsDaysStrtMonth

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
  write(6,F01) "debug level reset to ",level

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
!===============================================================================

end module shr_cal_mod
