!===============================================================================
! SVN $Id: shr_date_mod.F90 239 2006-02-08 19:02:33Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070423/shr/shr_date_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_date_mod -- date/time module, built upon a calendar module
!
! !DESCRIPTION:
!   Keeps track of model date, including elapsed seconds in current date. 
!
! !REVISION HISTORY:
!     2001-Sep-13 - B. Kauffman - created initial version
!
! !REMARKS:
!   This module is independant of a particular calendar, e.g. is ignorant of
!   whether the underlying calendar does or doesn't implement leap years.
!
! !INTERFACE: ------------------------------------------------------------------

module shr_date_mod

! !USES:

   use shr_cal_mod  ! underlying calendar
   use shr_sys_mod  ! system call wrappers
   use shr_kind_mod ! kinds

   implicit none

   private ! except

! !PUBLIC TYPES: 
 
   public :: shr_date
   type shr_date
      sequence               ! place in contiguous memory
      private                ! no public access to internal components
      integer(SHR_KIND_IN) :: y           ! calendar year
      integer(SHR_KIND_IN) :: m           ! calendar month
      integer(SHR_KIND_IN) :: d           ! calendar day
      integer(SHR_KIND_IN) :: s           ! elapsed seconds in current day
      integer(SHR_KIND_IN) :: cDate       ! coded calendar date (yymmdd)
      integer(SHR_KIND_IN) :: eDay        ! elsapsed days relative to calendar's reference date
      integer(SHR_KIND_IN) :: stepInDay   ! current time-step in current day
      integer(SHR_KIND_IN) :: stepsPerDay ! number of time-steps per day
   end type shr_date

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_date_adv1step   ! advance the date one time step
   public :: shr_date_advNextDay ! advance date to next day, 0 seconds
   public :: shr_date_initYMD    ! init date given YMD
   public :: shr_date_initEDay   ! init date given elapsed days
   public :: shr_date_initCDate  ! init date given coded date (yymmdd)
   public :: shr_date_getYMD     ! returns integers yy,mm,dd,sssss
   public :: shr_date_getEDay    ! returns elased day, sssss
   public :: shr_date_getCDate   ! returns coded date, ssssS
   public :: shr_date_getStepsPerDay ! returns steps per day
   public :: shr_date_getStepInDay ! returns step in the day
   public :: shr_date_getJulian  ! returns julian day number

   public :: assignment(=)       ! sets (date a) equal to (date b)
   public :: operator(==)        ! true iff (date a) == (date b)
   public :: operator(<)         ! true iff (date a) <  (date b)
   public :: operator(>)         ! true iff (date a) >  (date b)

! !PUBLIC DATA MEMBERS:

   real(SHR_KIND_R8),parameter,public :: shr_date_secsPerDay = 86400.0_SHR_KIND_R8 ! seconds per day

!EOP

   interface assignment(=)
      module procedure shr_date_assign
   end interface
   interface operator(==)
      module procedure shr_date_equals
   end interface
   interface operator(>)
      module procedure shr_date_greater
   end interface
   interface operator(<)
      module procedure shr_date_less
   end interface

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!===============================================================================

subroutine shr_date_assign(a,b)

   type(shr_date),intent(out) :: a
   type(shr_date),intent(in ) :: b

!-------------------------------------------------------------------------------
!  Make date a equal to date b 
!-------------------------------------------------------------------------------

   a%y           = b%y   
   a%m           = b%m   
   a%d           = b%d   
   a%s           = b%s   
   a%eday        = b%eday
   a%cDate       = b%cDate
   a%stepInDay   = b%stepInDay
   a%stepsPerDay = b%stepsPerDay

end subroutine shr_date_assign

!===============================================================================
!===============================================================================

function shr_date_equals(a,b)

   type(shr_date),intent(in) :: a,b
   logical                   :: shr_date_equals

!-------------------------------------------------------------------------------
!  Is date a equal to date b ??
!-------------------------------------------------------------------------------

   if (a%eday == b%eday .and. a%s == b%s) then
      shr_date_equals=.true.
   else
      shr_date_equals =.false.
   end if

end function shr_date_equals

!===============================================================================
!===============================================================================

function shr_date_greater(a,b)

   type(shr_date),intent(in) :: a,b
   logical                   :: shr_date_greater

!-------------------------------------------------------------------------------
!  Is date a greater than date b ??
!-------------------------------------------------------------------------------

   if      (a%eday < b%eday) then
      shr_date_greater = .false.
   else if (a%eday > b%eday) then
      shr_date_greater = .true.
   else if (a%s    > b%s   ) then
      shr_date_greater = .true.
   else 
      shr_date_greater = .false.
   end if

end function shr_date_greater

!===============================================================================
!===============================================================================

function shr_date_less(a,b)

   type(shr_date),intent(in) :: a,b
   logical                   :: shr_date_less

!-------------------------------------------------------------------------------
!  Is date a less than date b ??
!-------------------------------------------------------------------------------

   if      (a%eday < b%eday) then
      shr_date_less = .true.
   else if (a%eday > b%eday) then
      shr_date_less = .false.
   else if (a%s    < b%s   ) then
      shr_date_less = .true.
   else 
      shr_date_less = .false.
   end if

end function shr_date_less

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_initYMD - Initialize date given y,m,d, steps per day.
!
! !DESCRIPTION:
!     Initialize date given y,m,d, and steps per day.
!
! !REVISION HISTORY:
!     2001-Sep-13 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

type(shr_date) function shr_date_initYMD(y,m,d,ns)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: y     ! year
   integer(SHR_KIND_IN),intent(in) :: m     ! month
   integer(SHR_KIND_IN),intent(in) :: d     ! day
   integer(SHR_KIND_IN),intent(in) :: ns    ! number of steps per day.
 
!EOP

   !----- local -----
   type(shr_date)  :: date 

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (shr_cal_validYMD(y,m,d) ) then
      date%y = y
      date%m = m
      date%d = d
      date%s = 0
      date%stepInDay   = 0
      date%stepsPerDay = ns
      call shr_cal_ymd2date(y,m,d,date%cDate)
      call shr_cal_ymd2eDay(y,m,d,date%eDay )
   else
      write(6,*) 'ERROR: invalid y,m,d = ', y,m,d
      call shr_sys_abort('(shr_date_initYMD) invalid date')
   end if

   shr_date_initYMD = date

end function shr_date_initYMD

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_initEDay - Initialize date given an elapsed day
!
! !DESCRIPTION:
!     Initialize date given elapsed days and the number of steps per day.
!
! !REVISION HISTORY:
!     2003-May-27 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

type(shr_date) function shr_date_initEDay(eDay,ns)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: eDay  ! elpased days
   integer(SHR_KIND_IN),intent(in) :: ns    ! number of steps per day.
 
!EOP

   !----- local -----
   type(shr_date)       :: date   ! date to return
   integer(SHR_KIND_IN) :: y,m,d  ! year, month, day

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call shr_cal_eDay2ymd(eDay,y,m,d) ! convert eDay to y,m,d

   if ( shr_cal_validYMD(y,m,d) ) then
      date = shr_date_initYMD(y,m,d,ns)
   else
      write(6,*) 'ERROR: invalid eDay  = ', eDay
      write(6,*) 'ERROR: invalid y,m,d = ', y,m,d
      call shr_sys_abort('(shr_date_initEDay) invalid eDay')
   end if

   shr_date_initEDay = date

end function shr_date_initEDay

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_initCDate - Initialize date given coded date (yymmdd).
!
! !DESCRIPTION:
!     Initialize date given coded date (yymmdd) and the number of steps per day.
!
! !REVISION HISTORY:
!     2001-Sep-13 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

type(shr_date) function shr_date_initCDate(cDate,ns)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: cDate ! coded date (yymmdd)
   integer(SHR_KIND_IN),intent(in) :: ns    ! number of steps per day.
 
!EOP

   !----- local -----
   type(shr_date)       :: date   ! date to return
   integer(SHR_KIND_IN) :: y,m,d  ! year, month, day

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   if (shr_cal_validDate(cDate) ) then
      call shr_cal_date2ymd(cDate,y,m,d)
      date = shr_date_initYMD(y,m,d,ns)
   else
      write(6,*) 'ERROR: invalid cDate = ', cDate
      call shr_sys_abort('(shr_date_initCDate) invalid date')
   end if

   shr_date_initCDate = date

end function shr_date_initCDate

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_adv1step - advance the date by one time step
!
! !DESCRIPTION:
!     Advance the date by one time step.
!
! !REVISION HISTORY:
!     2001-Sep-13 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_date_adv1step(date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(inout) :: date  ! number of elapsed days

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   date%stepInDay = date%stepInDay + 1

   if (date%stepInDay < date%stepsPerDay) then
      date%s = nint((shr_date_secsPerDay*date%stepInDay)/date%stepsPerDay)
   else
      call shr_date_advNextDay(date)
   end if

end subroutine shr_date_adv1step

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_advNextDay - advance the date to start of next day.
!
! !DESCRIPTION:
!     Advance the date to the start of next day.
!
! !REVISION HISTORY:
!     2001-Sep-13 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_date_advNextDay(date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(inout) :: date  ! number of elapsed days

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   date%eDay      = date%eDay + 1
   date%stepInDay = 0
   date%s         = 0
   call shr_cal_eDay2ymd (date%eDay,date%y,date%m,date%d)
   call shr_cal_eDay2date(date%eDay,date%cDate)

end subroutine shr_date_advNextDay

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_getYMD - return yy,mm,dd,ss of date
!
! !DESCRIPTION:
!     return yy,mm,dd,ss of date.
!
! !REVISION HISTORY:
!     2001-Sep-13 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_date_getYMD(date,y,m,d,s)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(in) :: date  ! input date
   integer(SHR_KIND_IN),intent(out) :: y     ! year
   integer(SHR_KIND_IN),intent(out) :: m     ! month
   integer(SHR_KIND_IN),intent(out) :: d     ! day
   integer(SHR_KIND_IN),intent(out) :: s     ! elapsed seconds on date

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   y = date%y
   m = date%m
   d = date%d
   s = date%s

end subroutine shr_date_getYMD

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_getCDate - return coded date of a date
!
! !DESCRIPTION:
!     return coded date of a date
!
! !REVISION HISTORY:
!     2001-Sep-13 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_date_getCDATE(date,cDate,s)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(in) :: date  ! number of elapsed days
   integer(SHR_KIND_IN),intent(out) :: cDate ! coded date
   integer(SHR_KIND_IN),intent(out) :: s     ! elapsed seconds on date

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   cDate = date%cDate
   s     = date%s

end subroutine shr_date_getCDate

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_getEDay - return elapsed days of a date
!
! !DESCRIPTION:
!     return elapsed days of a date
!
! !REVISION HISTORY:
!     2001-Sep-13 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_date_getEDay(date,eDay,s)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(in) :: date  ! input date
   integer(SHR_KIND_IN),intent(out) :: eDay  ! elapsed days of date
   integer(SHR_KIND_IN),intent(out) :: s     ! elapsed seconds on date

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   eDay = date%eDay
   s    = date%s

end subroutine shr_date_getEDay

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_getStepsPerDay - return elapsed days of a date
!
! !DESCRIPTION:
!     return elapsed days of a date
!
! !REVISION HISTORY:
!     2001-Sep-13 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

integer function shr_date_getStepsPerDay(date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(in) :: date  ! input date

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   shr_date_getStepsPerDay = date%stepsPerDay

end function shr_date_getStepsPerDay

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_getStepInDay - return elapsed days of a date
!
! !DESCRIPTION:
!     return timestep in this day
!
! !REVISION HISTORY:
!     2001-Nov-13 - T. Craig - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

integer function shr_date_getStepInDay(date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(in) :: date  ! input date

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   shr_date_getStepInDay = date%StepInDay

end function shr_date_getStepInDay

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_date_getJulian - return julian day
!
! !DESCRIPTION:
!     return julian day
!
! !REVISION HISTORY:
!     2002-Oct-28 - R. Jacob -initial version
!
! !INTERFACE:  -----------------------------------------------------------------

real(SHR_KIND_R8) function shr_date_getJulian(date,shift)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date),intent(in)   :: date  ! input date
   integer(SHR_KIND_IN),intent(in),optional :: shift ! seconds to shift calculation

!EOP
  ! local
   integer(SHR_KIND_IN) :: totsec    ! total seconds

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   totsec = date%s

   if(present(shift)) totsec = totsec + shift

   shr_date_getJulian = shr_cal_elapsDaysStrtMonth(date%y,date%m) &
   + date%d + totsec/86400.0_SHR_KIND_R8

end function shr_date_getJulian

!===============================================================================
!===============================================================================

end module shr_date_mod
