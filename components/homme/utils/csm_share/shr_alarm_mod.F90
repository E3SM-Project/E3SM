!===============================================================================
! SVN $Id: shr_alarm_mod.F90 239 2006-02-08 19:02:33Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070423/shr/shr_alarm_mod.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: shr_alarm_mod -- date/time module, built upon calendar module
!
! !DESCRIPTION:
!   Keeps track of model date, including elapsed seconds in current date. 
!
! !REVISION HISTORY:
!     2002-Sep-17 - B. Kauffman - created initial version
!
! !REMARKS:
!   This module is independant of a particular calendar, e.g. is ignorant of
!   whether the underlying calendar does or doesn't implement leap years.
!
! !INTERFACE: ------------------------------------------------------------------

module shr_alarm_mod

! !USES:

   use shr_cal_mod  ! underlying calendar module
   use shr_date_mod ! underlying date module
   use shr_sys_mod  ! system call wrappers
   use shr_kind_mod ! kinds

   implicit none

   private ! except

! !PUBLIC TYPES: 
 
   public :: shr_alarm
   type shr_alarm
      private
      sequence
      integer(SHR_KIND_IN) :: type     ! type of alarm
      integer(SHR_KIND_IN) :: n        ! n wrt nDays or nMonths
      type(shr_date) :: oDate    ! offset date
      logical        :: isOn     ! true iff alarm is on
   end type shr_alarm

! !PUBLIC MEMBER FUNCTIONS:

   public :: shr_alarm_isOn         ! returns true iff alarm is on
   public :: shr_alarm_set          ! turn alarm on or off
   public :: shr_alarm_getType      ! returns alarm type
   public :: shr_alarm_initDate     ! initialize/return an alarm
   public :: shr_alarm_initMonthly  ! initialize/return an alarm
   public :: shr_alarm_initYearly   ! initialize/return an alarm
   public :: shr_alarm_initNDays    ! initialize/return an alarm
   public :: shr_alarm_initNMonths  ! initialize/return an alarm
   public :: shr_alarm_initNStep    ! initialize/return an alarm
   public :: shr_alarm_initifsec    ! initialize/return an alarm
   public :: shr_alarm_initifdays0  ! initialize/return an alarm
   public :: shr_alarm_initifday    ! initialize/return an alarm
   public :: shr_alarm_initifmon    ! initialize/return an alarm
   public :: shr_alarm_initifyear   ! initialize/return an alarm
   public :: shr_alarm_initNone     ! initialize/return an alarm
   public :: shr_alarm_dump         ! dump alarm internals for debugging

! !PUBLIC DATA MEMBERS:

   integer(SHR_KIND_IN),parameter,public :: shr_alarm_date    =  1 ! goes on on given date
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_monthly =  2 ! goes on start of month
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_yearly  =  3 ! goes on start of year 
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_nDays   =  4 ! periodic alarm
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_nMonths =  5 ! periodic alarm
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_nStep   =  6 ! periodic alarm
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_ifsec   =  7 ! goes on sec value
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_ifdays0 =  8 ! goes on day value
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_ifday   =  9 ! goes on day value
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_ifmon   = 10 ! goes on month value
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_ifyear  = 11 ! goes on year value
   integer(SHR_KIND_IN),parameter,public :: shr_alarm_none    =  1 ! goes on on given date

!EOP

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
contains
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_isOn - returns true iff alarm is on.
!
! !DESCRIPTION:
!    turns alarm on (when appropriate) and returns true iff alarm is on.
!
! !REVISION HISTORY:
!     2002-Sep-17 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

logical function shr_alarm_isOn(date,alarm)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date ),intent(in)    :: date            ! date
   type(shr_alarm),intent(inout) :: alarm           ! alarm

!EOP

   !----- local -----
   integer(SHR_KIND_IN) ::  y, m, d, s ! input  date's year,month,day,sec
   integer(SHR_KIND_IN) :: oy,om,od,os ! offset date's year,month,day,sec
   integer(SHR_KIND_IN) :: eDay        ! input  date's eDay
   integer(SHR_KIND_IN) :: SPDay       ! input  date's Steps Per Day
   integer(SHR_KIND_IN) :: SIDay       ! input  date's Step In Day
   integer(SHR_KIND_IN) :: oEDay       ! offset date's eDay

   !----- formats -----
   character(len=*),parameter :: F01 = "('(shr_alarm_isOn) ',a,i8)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   
   shr_alarm_isOn = .false.

   call shr_date_getYMD(date,y ,m ,d ,s )

   if ( alarm%isOn) then !--- alarm stays on until shr_alarm_set is called ---
      shr_alarm_isOn = .true.
   else 
      if      (alarm%type == shr_alarm_date    ) then
         if (s == 0) then
           call shr_date_getYMD(alarm%oDate,oy,om,od,os)
           if (y==oy .and. m==om .and. d==od) shr_alarm_isOn = .true.
         endif
      else if (alarm%type == shr_alarm_monthly ) then
         if (s == 0) then
           if (d == 1) shr_alarm_isOn = .true.
         endif
      else if (alarm%type == shr_alarm_yearly  ) then
         if (s == 0) then
           if (m==1 .and. d == 1) shr_alarm_isOn = .true.
         endif
      else if (alarm%type == shr_alarm_nDays   ) then
         if (s == 0) then
           call shr_date_getEday(      date, eDay , s)
           call shr_date_getEday(alarm%oDate,oEDay,os)
           if ( mod((eDay-oEDay),alarm%n) == 0 ) shr_alarm_isOn = .true.
         endif
      else if (alarm%type == shr_alarm_nMonths ) then
         if (s == 0) then
           call shr_date_getYMD(alarm%oDate,oy,om,od,os)
           if ( mod(((12*y+m)-(12*oy+om)),alarm%n) == 0 ) then
              !--- it's the right month, is it the right day? ---
              od = min(od,shr_cal_numDaysInMonth(y,m))
              if ( d == od ) shr_alarm_isOn = .true.
           end if
         endif
      else if (alarm%type == shr_alarm_nStep ) then
         call shr_date_getEday       ( date, eDay  , s)
         SPDay = shr_date_getStepsPerDay( date )
         SIDay = shr_date_getStepInDay  ( date )
         if ( mod((eDay*SPDay + SIDay),alarm%n) == 0) then
            shr_alarm_isOn = .true.
         endif
      else if (alarm%type == shr_alarm_ifsec ) then
         call shr_date_getYMD(date,oy,om,od,os)
         if (os == alarm%n) shr_alarm_isOn = .true.
      else if (alarm%type == shr_alarm_ifdays0 ) then
         if (s == 0) then
           call shr_date_getYMD(date,oy,om,od,os)
           if (od == alarm%n) shr_alarm_isOn = .true.
         endif
      else if (alarm%type == shr_alarm_ifday ) then
         call shr_date_getYMD(date,oy,om,od,os)
         if (od == alarm%n) shr_alarm_isOn = .true.
      else if (alarm%type == shr_alarm_ifmon ) then
         call shr_date_getYMD(date,oy,om,od,os)
         if (om == alarm%n) shr_alarm_isOn = .true.
      else if (alarm%type == shr_alarm_ifyear ) then
         call shr_date_getYMD(date,oy,om,od,os)
         if (oy == alarm%n) shr_alarm_isOn = .true.
      else if (alarm%type == shr_alarm_none ) then
         !--- do nothing ---
      else 
         write(6,F01) "ERROR: unrecognized alarm type = ",alarm%type
         call shr_sys_abort("(shr_alarm_isOn)")
      end if
   end if

end function shr_alarm_isOn

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_set - turn alarm on or off
!
! !DESCRIPTION:
!    Turns alarm on or off.
!
! !REVISION HISTORY:
!     2002-Sep-17 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_alarm_set(alarm,isOn)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_alarm),intent(inout) :: alarm  ! alarm
   logical        ,intent(in)    :: isOn

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   
    alarm%isOn = isOn

end subroutine shr_alarm_set

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_getType - return alarm type
!
! !DESCRIPTION:
!    Return alarm type.
!
! !REVISION HISTORY:
!     2002-Sep-17 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

integer function shr_alarm_getType(alarm)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_alarm),intent(inout) :: alarm  ! alarm

!EOP

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   
    shr_alarm_getType = alarm%type

end function shr_alarm_getType

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initDate(date) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm.
!
! !REVISION HISTORY:
!     2002-Sep-17 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initDate(date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_date ),intent(in)   :: date 
   type(shr_alarm)              :: shr_alarm_initDate

!EOP

   !----- local -----
   type(shr_alarm)              :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initDate) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_date
   alarm%n     = -999                       ! is irrelavant for this option
   alarm%oDate = date
   alarm%isOn = .false.

   shr_alarm_initDate = alarm

end function shr_alarm_initDate

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initMonthly(n) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm: goes off every n days relative to the given date.
!
! !REVISION HISTORY:
!     2002-Sep-17 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initMonthly()

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_alarm)           :: shr_alarm_initMonthly

!EOP

   !----- local -----
   type(shr_alarm)            :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initMonthly) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_monthly
   alarm%n     = -999                       ! is irrelavant for this option
   alarm%oDate = shr_date_initYMD(0,1,1,0)  ! is irrelavant for this option
   alarm%isOn = .false.

   shr_alarm_initMonthly = alarm

end function shr_alarm_initMonthly

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initYearly(n) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm: goes off every n days relative to the given date.
!
! !REVISION HISTORY:
!     2002-Sep-17 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initYearly()

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_alarm)           :: shr_alarm_initYearly

!EOP

   !----- local -----
   type(shr_alarm)            :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initYearly) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_yearly
   alarm%n     = -999                       ! is irrelavant for this option
   alarm%oDate = shr_date_initYMD(0,1,1,0)  ! is irrelavant for this option
   alarm%isOn = .false.

   shr_alarm_initYearly = alarm

end function shr_alarm_initYearly

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initNDays(n) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm: goes off every n days relative to the given date.
!
! !REVISION HISTORY:
!     2002-Sep-17 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initNDays(n,date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: n 
   type(shr_date),intent(in)       :: date
   type(shr_alarm)                 :: shr_alarm_initNDays

!EOP

   !----- local -----
   type(shr_alarm)            :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initNDays) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_nDays
   alarm%n     = n
   alarm%oDate = date
   alarm%isOn = .false.

   shr_alarm_initNDays = alarm

end function shr_alarm_initNDays

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initNMonths(n) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm: goes off every n months relative to the given date.
!
! !REVISION HISTORY:
!     2002-Sep-17 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initNMonths(n,date)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: n 
   type(shr_date),intent(in)       :: date
   type(shr_alarm)                 :: shr_alarm_initNMonths

!EOP

   !----- local -----
   type(shr_alarm)            :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initNMonths) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_nMonths
   alarm%n     = n
   alarm%oDate = date
   alarm%isOn = .false.

   shr_alarm_initNMonths = alarm

end function shr_alarm_initNMonths

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initNStep(n) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm: goes off every n steps relative to the given step 0.
!
! !REVISION HISTORY:
!     2002-Nov-17 - T. Craig - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initNStep(n)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: n 
   type(shr_alarm)                 :: shr_alarm_initNStep

!EOP

   !----- local -----
   type(shr_alarm)            :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initNStep) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_NStep
   alarm%n     = n
   alarm%oDate = shr_date_initYMD(0,1,1,0)  ! is irrelavant for this option
   alarm%isOn = .false.

   shr_alarm_initNStep = alarm

end function shr_alarm_initNStep

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initifsec(n) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm: goes off when sec is equal to n
!
! !REVISION HISTORY:
!     2002-Nov-26 - T. Craig - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initifsec(n)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: n 
   type(shr_alarm)                 :: shr_alarm_initifsec

!EOP

   !----- local -----
   type(shr_alarm)            :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initifsec) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_ifsec
   alarm%n     = n
   alarm%oDate = shr_date_initYMD(0,1,1,0)  ! is irrelavant for this option
   alarm%isOn = .false.

   shr_alarm_initifsec = alarm

end function shr_alarm_initifsec

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initifdays0(n) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm: goes off when day is equal to n and second is zero
!
! !REVISION HISTORY:
!     2003-Jun-15 - T. Craig - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initifdays0(n)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: n 
   type(shr_alarm)                 :: shr_alarm_initifdays0

!EOP

   !----- local -----
   type(shr_alarm)            :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initifdays0) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_ifdays0
   alarm%n     = n
   alarm%oDate = shr_date_initYMD(0,1,1,0)  ! is irrelavant for this option
   alarm%isOn = .false.

   shr_alarm_initifdays0 = alarm

end function shr_alarm_initifdays0

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initifday(n) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm: goes off when day is equal to n
!
! !REVISION HISTORY:
!     2002-Nov-26 - T. Craig - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initifday(n)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: n 
   type(shr_alarm)                 :: shr_alarm_initifday

!EOP

   !----- local -----
   type(shr_alarm)            :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initifday) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_ifday
   alarm%n     = n
   alarm%oDate = shr_date_initYMD(0,1,1,0)  ! is irrelavant for this option
   alarm%isOn = .false.

   shr_alarm_initifday = alarm

end function shr_alarm_initifday

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initifmon(n) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm: goes off when month is equal to n
!
! !REVISION HISTORY:
!     2002-Nov-26 - T. Craig - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initifmon(n)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: n 
   type(shr_alarm)                 :: shr_alarm_initifmon

!EOP

   !----- local -----
   type(shr_alarm)            :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initifmon) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_ifmon
   alarm%n     = n
   alarm%oDate = shr_date_initYMD(0,1,1,0)  ! is irrelavant for this option
   alarm%isOn = .false.

   shr_alarm_initifmon = alarm

end function shr_alarm_initifmon

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initifyear(n) -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm: goes off when year is equal to n
!
! !REVISION HISTORY:
!     2002-Nov-26 - T. Craig - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initifyear(n)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   integer(SHR_KIND_IN),intent(in) :: n 
   type(shr_alarm)                 :: shr_alarm_initifyear

!EOP

   !----- local -----
   type(shr_alarm)            :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initifyear) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_ifyear
   alarm%n     = n
   alarm%oDate = shr_date_initYMD(0,1,1,0)  ! is irrelavant for this option
   alarm%isOn = .false.

   shr_alarm_initifyear = alarm

end function shr_alarm_initifyear

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_initNone() -- initialize an alarm
!
! !DESCRIPTION:
!    Initialize an alarm.
!
! !REVISION HISTORY:
!     2002-Sep-17 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

function shr_alarm_initNone()

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_alarm)              :: shr_alarm_initNone

!EOP

   !----- local -----
   type(shr_alarm)              :: alarm

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_initNone) ',8a)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   alarm%type  = shr_alarm_none
   alarm%n     = -999                       ! is irrelavant for this option
   alarm%oDate = shr_date_initYMD(0,1,1,0)  ! is irrelavant for this option
   alarm%isOn = .false.

   shr_alarm_initNone = alarm

end function shr_alarm_initNone

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_alarm_dump -- dump alarm's internal data
!
! !DESCRIPTION:
!    Dump dump alarm's internal data for debugging.
!
! !REVISION HISTORY:
!     2003-Jan-06 - B. Kauffman - initial version.
!
! !INTERFACE:  -----------------------------------------------------------------

subroutine shr_alarm_dump(a)

   implicit none

! !INPUT/OUTPUT PARAMETERS:

   type(shr_alarm)  :: a  ! alarm in question

!EOP

   !----- local -----
   integer(SHR_KIND_IN) :: cDate ! coded date
   integer(SHR_KIND_IN) :: sec   ! seconds on date

   !----- formats -----
   character(len=*),parameter :: F00 = "('(shr_alarm_dump) ',8a)"
   character(len=*),parameter :: F01 = "('(shr_alarm_dump) ',a,i8)"
   character(len=*),parameter :: F02 = "('(shr_alarm_dump) ',a,i8.4,i6,'sec')"
   character(len=*),parameter :: F03 = "('(shr_alarm_dump) ',a,l7)"

!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

   call shr_date_getCDate(a%oDate,cDate,sec)

   write(6,F00) "dump alarm internal data..."
   write(6,F01) "* type  = ", a%type    ! type of alarm
   write(6,F01) "* n     = ", a%n       ! n wrt nDays or nMonths
   write(6,F02) "* oDate = ", cDate,sec ! offset date
   write(6,F03) "* isOn  = ", a%isOn    ! true iff alarm is on

end subroutine shr_alarm_dump

!===============================================================================
!===============================================================================

end module shr_alarm_mod

