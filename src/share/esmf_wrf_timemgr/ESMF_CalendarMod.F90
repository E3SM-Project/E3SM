! $Id$
!
! Earth System Modeling Framework
! Copyright 2002-2003, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the GPL.
!
!==============================================================================
!
!     ESMF Calendar Module
      module ESMF_CalendarMod
!
!==============================================================================
!
! This file contains the Calendar class definition and all Calendar class
! methods.
!
!------------------------------------------------------------------------------
! INCLUDES
#include <ESMF_TimeMgr.inc>

!==============================================================================
!BOPI
! !MODULE: ESMF_CalendarMod
!
! !DESCRIPTION:
! Part of Time Manager F90 API wrapper of C++ implemenation
!
! Defines F90 wrapper entry points for corresponding
! C++ class { \tt ESMC\_Calendar} implementation
!
! See {\tt ../include/ESMC\_Calendar.h} for complete description
!
!------------------------------------------------------------------------------
! !USES:
      ! inherit from ESMF base class
      use ESMF_BaseMod

      ! inherit from base time class
      use ESMF_BaseTimeMod

      implicit none
!
!------------------------------------------------------------------------------
! !PRIVATE TYPES:
      private
!------------------------------------------------------------------------------

      INTEGER, PARAMETER :: mday(MONTHS_PER_YEAR)   &
                          = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      INTEGER, PARAMETER :: mdayleap(MONTHS_PER_YEAR) &
                          = (/31,29,31,30,31,30,31,31,30,31,30,31/)
      INTEGER, DIMENSION(365) :: daym
      INTEGER, DIMENSION(366) :: daymleap
      INTEGER :: mdaycum(0:MONTHS_PER_YEAR)
      INTEGER :: mdayleapcum(0:MONTHS_PER_YEAR)
      TYPE(ESMF_BaseTime), TARGET :: monthbdys(0:MONTHS_PER_YEAR)
      TYPE(ESMF_BaseTime), TARGET :: monthbdysleap(0:MONTHS_PER_YEAR)
      TYPE(ESMF_BaseTime), TARGET :: monthedys(0:MONTHS_PER_YEAR)
      TYPE(ESMF_BaseTime), TARGET :: monthedysleap(0:MONTHS_PER_YEAR)


!------------------------------------------------------------------------------
!     ! ESMF_CalKind_Flag
!
!     ! F90 "enum" type to match C++ ESMC_CalKind_Flag enum

      type ESMF_CalKind_Flag
        integer :: caltype
      end type

      type(ESMF_CalKind_Flag), parameter :: &
                               ESMF_CALKIND_GREGORIAN =  ESMF_CalKind_Flag(1), &
                               ESMF_CALKIND_NOLEAP =     ESMF_CalKind_Flag(2)

!      type(ESMF_CalKind_Flag), parameter :: &
!                               ESMF_CALKIND_GREGORIAN =  ESMF_CalKind_Flag(1), &
!                               ESMF_CALKIND_JULIAN =     ESMF_CalKind_Flag(2), &
!                           ! like Gregorian, except Feb always has 28 days
!                               ESMF_CALKIND_NOLEAP =     ESMF_CalKind_Flag(3), &
!                           ! 12 months, 30 days each
!                               ESMF_CALKIND_360DAY =     ESMF_CalKind_Flag(4), &
!                           ! user defined
!                               ESMF_CALKIND_GENERIC =    ESMF_CalKind_Flag(5), &
!                           ! track base time seconds only
!                               ESMF_CALKIND_NOCALENDAR = ESMF_CalKind_Flag(6)

!------------------------------------------------------------------------------
!     ! ESMF_Calendar
!
!     ! F90 class type to match C++ Calendar class in size only;
!     !  all dereferencing within class is performed by C++ implementation
!
!------------------------------------------------------------------------------
!
!     ! ESMF_DaysPerYear
!
      type ESMF_DaysPerYear
        integer :: D = 0    ! whole days per year
        integer :: Dn = 0   ! fractional days per year numerator
        integer :: Dd = 1   ! fractional days per year denominator
      end type              ! e.g. for Venus, D=0, Dn=926, Dd=1000
!
!------------------------------------------------------------------------------
!     ! ESMF_Calendar
!
!
      type ESMF_Calendar
        type(ESMF_CalKind_Flag) :: Type
        logical :: Set = .false.
        integer, dimension(MONTHS_PER_YEAR) :: DaysPerMonth = 0
        integer :: SecondsPerDay = 0
        integer :: SecondsPerYear = 0
        type(ESMF_DaysPerYear) :: DaysPerYear
      end type
!------------------------------------------------------------------------------
! !PUBLIC DATA: added by Juanxiong He, in order to breakthe cycle call between
! ESMF_Stubs and ESMF_Time
   TYPE(ESMF_Calendar), public, save, pointer :: defaultCal   ! Default Calendar
   TYPE(ESMF_Calendar), public, save, pointer :: gregorianCal ! gregorian Calendar
   TYPE(ESMF_Calendar), public, save, pointer :: noleapCal    ! noleap Calendar

!
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
      public initdaym
!      public mday
!      public mdayleap
!      public monthbdys
!      public monthbdysleap
!      public monthedys
!      public monthedysleap
!      public daym
!      public daymleap
!      public mdaycum
!      public mdayleapcum
      public ndaysinmonth
      public nsecondsinmonth
      public ndaysinyear
      public nsecondsinyear
      public nmonthinyearsec
      public ndayinyearsec
      public nsecondsinyearmonth
      public isleap
      public ESMF_CalKind_Flag
      public ESMF_CALKIND_GREGORIAN, ESMF_CALKIND_NOLEAP
!             ESMF_CALKIND_360DAY, ESMF_CALKIND_NOCALENDAR
!      public ESMF_CAL_JULIAN
!      public ESMF_CAL_GENERIC
      public ESMF_Calendar
      public ESMF_DaysPerYear

!------------------------------------------------------------------------------
!
! !PUBLIC MEMBER FUNCTIONS:
      public ESMF_CalendarCreate

! Required inherited and overridden ESMF_Base class methods

      public ESMF_CalendarInitialized ! Only in this implementation, intended
                                      ! to be private within ESMF methods
!EOPI

!------------------------------------------------------------------------------
! The following line turns the CVS identifier string into a printable variable.
      character(*), parameter, private :: version = &
      '$Id$'

!==============================================================================

      contains


!==============================================================================
!BOP
! !IROUTINE: ESMF_CalendarCreate - Create a new ESMF Calendar of built-in type

! !INTERFACE:
      ! Private name; call using ESMF_CalendarCreate()
      function ESMF_CalendarCreate(name, calkindflag, rc)

! !RETURN VALUE:
      type(ESMF_Calendar) :: ESMF_CalendarCreate

! !ARGUMENTS:
      character (len=*),       intent(in),  optional :: name
      type(ESMF_CalKind_Flag), intent(in)            :: calkindflag
      integer,                 intent(out), optional :: rc

! !DESCRIPTION:
!     Creates and sets a {\tt calendar} to the given built-in
!     {\tt ESMF\_CalKind_Flag}.
!
!     This is a private method; invoke via the public overloaded entry point
!     {\tt ESMF\_CalendarCreate()}.
!
!     The arguments are:
!     \begin{description}
!     \item[{[name]}]
!          The name for the newly created calendar.  If not specified, a
!          default unique name will be generated: "CalendarNNN" where NNN
!          is a unique sequence number from 001 to 999.
!     \item[calkindflag]
!          The built-in {\tt ESMF\_CalKind_Flag}.  Valid values are:
!            {\tt ESMF\_CAL\_360DAY}, {\tt ESMF\_CAL\_GREGORIAN},
!            {\tt ESMF\_CAL\_JULIANDAY}, {\tt ESMF\_CAL\_NOCALENDAR}, and
!            {\tt ESMF\_CAL\_NOLEAP}.
!          See the "Time Manager Reference" document for a description of
!          each calendar type.
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
!EOP
! !REQUIREMENTS:
!     TMGn.n.n
      type(ESMF_DaysPerYear) :: dayspy

      if ( present(rc) ) rc = ESMF_FAILURE
! Calendar is hard-coded.  Use ESMF library if more flexibility is needed.
!      write(6,*) 'tcx ESMF_CalendarCreate ',calkindflag%caltype, ESMF_CALKIND_NOLEAP%caltype, ESMF_CALKIND_GREGORIAN%caltype
      if ( calkindflag%caltype  == ESMF_CALKIND_NOLEAP%caltype ) then
!         write(6,*) 'tcx ESMF_CalendarCreate: initialize noleap calendar '
         ESMF_CalendarCreate%Type = ESMF_CALKIND_NOLEAP
      elseif ( calkindflag%caltype  == ESMF_CALKIND_GREGORIAN%caltype ) then
!         write(6,*) 'tcx ESMF_CalendarCreate: initialize gregorian calendar '
         ESMF_CalendarCreate%Type = ESMF_CALKIND_GREGORIAN
      else
!         write(6,*) 'tcx ESMF_CalendarCreate: ERROR initialize invalid calendar'
         call wrf_error_fatal( "Error:: ESMF_CalendarCreate invalid calendar")
      endif

!$$$ This is a bug on some systems -- need initial value set by compiler at
!$$$ startup.
      ESMF_CalendarCreate%Set = .true.
      ESMF_CalendarCreate%SecondsPerDay = SECONDS_PER_DAY
! DaysPerYear and SecondsPerYear are incorrect for Gregorian calendars...
      dayspy%D = size(daym)
      dayspy%Dn = 0
      dayspy%Dd = 1
      ESMF_CalendarCreate%DaysPerYear = dayspy
      ESMF_CalendarCreate%SecondsPerYear = ESMF_CalendarCreate%SecondsPerDay &
                                       * dayspy%D
      ESMF_CalendarCreate%DaysPerMonth(:) = mday(:)

      if ( present(rc) ) rc = ESMF_SUCCESS

      end function ESMF_CalendarCreate


!==============================================================================
!BOP
! !IROUTINE: ESMF_CalendarInitialized - check if calendar was created

! !INTERFACE:
      function ESMF_CalendarInitialized(calendar)

! !RETURN VALUE:
      logical ESMF_CalendarInitialized

! !ARGUMENTS:
      type(ESMF_Calendar), intent(in)            :: calendar

! !DESCRIPTION:
!EOP
! !REQUIREMENTS:
!     TMGn.n.n
        ESMF_CalendarInitialized = calendar%set
        if ( calendar%SecondsPerDay == 0 ) &
              ESMF_CalendarInitialized = .false.

     end function ESMF_CalendarInitialized

!==============================================================================

SUBROUTINE initdaym
  IMPLICIT NONE
  INTEGER i,j,m

  m = 1
  mdaycum(0) = 0
!$$$ push this down into ESMF_BaseTime constructor
  monthbdys(0)%S  = 0
  monthbdys(0)%Sn = 0
  monthbdys(0)%Sd = 0
  DO i = 1,MONTHS_PER_YEAR
    DO j = 1,mday(i)
      daym(m) = i
      m = m + 1
    ENDDO
    mdaycum(i) = mdaycum(i-1) + mday(i)
!$$$ push this down into ESMF_BaseTime constructor
    monthbdys(i)%S  = SECONDS_PER_DAY * INT( mdaycum(i), ESMF_KIND_I8 )
    monthbdys(i)%Sn = 0
    monthbdys(i)%Sd = 0
  ENDDO
  ! End of month seconds, day before the beginning of next month
  DO i = 0,MONTHS_PER_YEAR
    j = i + 1
    if ( i == MONTHS_PER_YEAR ) j = 0
    monthedys(i)   = monthbdys(j)
    monthedys(i)%S = monthedys(i)%S - SECONDS_PER_DAY
  ENDDO

  m = 1
  mdayleapcum(0) = 0
!$$$ push this down into ESMF_BaseTime constructor
  monthbdysleap(0)%S  = 0
  monthbdysleap(0)%Sn = 0
  monthbdysleap(0)%Sd = 0
  DO i = 1,MONTHS_PER_YEAR
    DO j = 1,mdayleap(i)
      daymleap(m) = i
      m = m + 1
    ENDDO
    mdayleapcum(i) = mdayleapcum(i-1) + mdayleap(i)
!$$$ push this down into ESMF_BaseTime constructor
    monthbdysleap(i)%S  = SECONDS_PER_DAY * INT( mdayleapcum(i), ESMF_KIND_I8 )
    monthbdysleap(i)%Sn = 0
    monthbdysleap(i)%Sd = 0
  ENDDO
  ! End of month seconds, day before the beginning of next month
  DO i = 0,MONTHS_PER_YEAR
    j = i + 1
    if ( i == MONTHS_PER_YEAR ) j = 0
    monthedysleap(i)   = monthbdysleap(j)
    monthedysleap(i)%S = monthedysleap(i)%S - SECONDS_PER_DAY
  ENDDO

END SUBROUTINE initdaym

!==============================================================================

integer(esmf_kind_i8) FUNCTION nsecondsinyear ( year, calkindflag )
  ! Compute the number of seconds in the given year
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: year
  type(ESMF_CalKind_Flag),intent(in) :: calkindflag

  nsecondsinyear = SECONDS_PER_DAY * INT( ndaysinyear(year, calkindflag) , ESMF_KIND_I8 )

END FUNCTION nsecondsinyear

!==============================================================================

integer function ndaysinmonth( year,month,calkindflag)
  ! Compute number of days in month for year, month, cal
  IMPLICIT NONE
  INTEGER, INTENT(in) :: year,month
  type(ESMF_CalKind_Flag),intent(in) :: calkindflag
  ! locals

  IF ( ( MONTH < 1 ) .OR. ( MONTH > MONTHS_PER_YEAR ) ) THEN
    CALL wrf_error_fatal( 'ERROR ndaysinmonth:  MONTH out of range' )
  ENDIF

  IF ( isleap(year,calkindflag) ) THEN
    ndaysinmonth = mdayleap(month)
  ELSE
    ndaysinmonth = mday(month)
  ENDIF

END function ndaysinmonth
!==============================================================================

integer(esmf_kind_i8) function nsecondsinmonth( year,month,calkindflag)
  ! Compute number of days in month for year, month, cal
  IMPLICIT NONE
  INTEGER, INTENT(in) :: year,month
  type(ESMF_CalKind_Flag),intent(in) :: calkindflag
  ! locals

  nsecondsinmonth = ndaysinmonth(year,month,calkindflag)*SECONDS_PER_DAY

END function nsecondsinmonth

!==============================================================================

integer function nmonthinyearsec(year,basetime,calkindflag)
  ! Compute month for year, basetime, cal
  IMPLICIT NONE
  INTEGER, INTENT(in) :: year
  type(ESMF_BaseTime), intent(in) :: basetime
  type(ESMF_CalKind_Flag),intent(in) :: calkindflag
  ! locals
  TYPE(ESMF_BaseTime), pointer :: MMbdys(:)
  integer :: mm,i

  IF ( isleap(year,calkindflag) ) THEN
    MMbdys => monthbdysleap
  ELSE
    MMbdys => monthbdys
  ENDIF
  MM = -1
  DO i = 1,MONTHS_PER_YEAR
    IF ( ( basetime >= MMbdys(i-1) ) .AND. ( basetime < MMbdys(i) ) ) THEN
      MM = i
      EXIT
    ENDIF
  ENDDO
  IF ( MM == -1 ) THEN
    CALL wrf_error_fatal( 'nmonthinyearsec:  could not extract month of year from time' )
  ENDIF
  nmonthinyearsec = mm

END function nmonthinyearsec

!==============================================================================
integer function ndayinyearsec(year, basetime, calkindflag)
  ! Compute day of year for year, basetime, cal
  IMPLICIT NONE
  INTEGER, INTENT(in) :: year
  type(ESMF_BaseTime), intent(in) :: basetime
  type(ESMF_CalKind_Flag),intent(in) :: calkindflag
  ! locals
  TYPE(ESMF_BaseTime), pointer :: MMbdys(:)
  TYPE(ESMF_BaseTime) :: tmpbasetime
  integer :: mm

  mm = nmonthinyearsec(year, basetime, calkindflag)

  IF ( isleap(year,calkindflag) ) THEN
    MMbdys => monthbdysleap
  ELSE
    MMbdys => monthbdys
  ENDIF
  tmpbasetime = basetime - MMbdys(mm-1)
  ndayinyearsec = ( tmpbasetime%S / SECONDS_PER_DAY ) + 1

end function ndayinyearsec
!==============================================================================
integer(esmf_kind_i8) function nsecondsinyearmonth(year, month, calkindflag)
  ! Compute number of seconds from start of year for year, month, cal
  IMPLICIT NONE
  INTEGER, INTENT(in) :: year
  INTEGER, INTENT(in) :: month
  type(ESMF_CalKind_Flag),intent(in) :: calkindflag

  ! locals
  TYPE(ESMF_BaseTime), pointer :: MMbdys(:)

  IF ( ( MONTH < 1 ) .OR. ( MONTH > MONTHS_PER_YEAR ) ) THEN
    CALL wrf_error_fatal( 'ERROR nsecondsinyearmonth():  MONTH out of range' )
  ENDIF

  IF ( isleap(year, calkindflag) ) THEN
    MMbdys => monthbdysleap
  ELSE
    MMbdys => monthbdys
  ENDIF

  nsecondsinyearmonth = MMbdys(month-1)%s

end function nsecondsinyearmonth
!==============================================================================

integer FUNCTION ndaysinyear ( year,calkindflag )
  ! Compute the number of days in the given year
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: year
  type(ESMF_CalKind_Flag),intent(in) :: calkindflag

  IF ( isleap( year,calkindflag ) ) THEN
    ndaysinyear = 366
  ELSE
    ndaysinyear = 365
  ENDIF
END FUNCTION ndaysinyear

!==============================================================================

logical FUNCTION isleap ( year, calkindflag )
  ! Compute the number of days in February for the given year
  IMPLICIT NONE
  INTEGER,intent(in)  :: year
  type(ESMF_CalKind_Flag) :: calkindflag
  ! local
  INTEGER :: lyear

  lyear = abs(year)  ! make sure it handles negative years

  isleap = .false. ! By default, February has 28 days ...

  if (calkindflag%caltype == ESMF_CALKIND_GREGORIAN%caltype) then
     IF (MOD(lyear,4).eq.0) THEN
        isleap = .true.  ! But every four years, it has 29 days ...
        IF (MOD(lyear,100).eq.0) THEN
           isleap = .false.  ! Except every 100 years, when it has 28 days ...
           IF (MOD(lyear,400).eq.0) THEN
              isleap = .true.  ! Except every 400 years, when it has 29 days.
           END IF
        END IF
     END IF
  endif

END FUNCTION isleap

!==============================================================================
end module ESMF_CalendarMod
