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
!     ESMF Time Module
      module ESMF_TimeMod
!
!==============================================================================
!
! This file contains the Time class definition and all Time class methods.
!
!------------------------------------------------------------------------------
! INCLUDES
#include <ESMF_TimeMgr.inc>

!==============================================================================
!BOPI
! !MODULE: ESMF_TimeMod
!
! !DESCRIPTION:
! Part of Time Manager F90 API wrapper of C++ implemenation
!
! Defines F90 wrapper entry points for corresponding
! C++ class {\tt ESMC\_Time} implementation
!
! See {\tt ../include/ESMC\_Time.h} for complete description
!
!------------------------------------------------------------------------------
! !USES:
      ! inherit from ESMF base class
      use ESMF_BaseMod

      ! inherit from base time class
      use ESMF_BaseTimeMod

      ! associated derived types
      use ESMF_TimeIntervalMod
      use ESMF_CalendarMod
      use ESMF_ShrTimeMod, only : ESMF_Time
! added by Jhe
      use ESMF_Stubs

      implicit none
!
!------------------------------------------------------------------------------
! !PRIVATE TYPES:
      private
!------------------------------------------------------------------------------
!     ! ESMF_Time
!
!     ! F90 class type to match C++ Time class in size only;
!     !  all dereferencing within class is performed by C++ implementation

! move to ESMF_ShrTimeMod
!     type ESMF_Time
!       type(ESMF_BaseTime) :: basetime           ! inherit base class
!       ! time instant is expressed as year + basetime
!       integer :: YR
!       type(ESMF_Calendar), pointer :: calendar  ! associated calendar
!     end type
!------------------------------------------------------------------------------
! !PUBLIC DATA:

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
      public ESMF_Time
!------------------------------------------------------------------------------
!
! !PUBLIC MEMBER FUNCTIONS:
      public ESMF_TimeGet
      public ESMF_TimeSet
      public ESMF_TimePrint

! Required inherited and overridden ESMF_Base class methods

      public ESMF_TimeCopy
      public ESMF_SetYearWidth

! !PRIVATE MEMBER FUNCTIONS:

      private ESMF_TimeGetDayOfYear
      private ESMF_TimeGetDayOfYearInteger

! Inherited and overloaded from ESMF_BaseTime

      public operator(+)
      public ESMF_TimeInc

      public operator(-)
      private ESMF_TimeDec
      private ESMF_TimeDiff

      public operator(.EQ.)
      public ESMF_TimeEQ

      public operator(.NE.)
      public ESMF_TimeNE

      public operator(.LT.)
      public ESMF_TimeLT

      public operator(.GT.)
      public ESMF_TimeGT

      public operator(.LE.)
      public ESMF_TimeLE

      public operator(.GE.)
      public ESMF_TimeGE

!EOPI

!------------------------------------------------------------------------------
! The following line turns the CVS identifier string into a printable variable.
      character(*), parameter, private :: version = &
      '$Id$'

      integer :: yearWidth = 4

!==============================================================================
!
! INTERFACE BLOCKS
!
!==============================================================================
!BOP
! !INTERFACE:
      interface ESMF_TimeGetDayOfYear

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeGetDayOfYearInteger

! !DESCRIPTION:
!     This interface overloads the {\tt ESMF\_GetDayOfYear} method
!     for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(+)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeInc, ESMF_TimeInc2

! !DESCRIPTION:
!     This interface overloads the + operator for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface assignment (=)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeCopy

! !DESCRIPTION:
!     This interface overloads the = operator for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(-)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeDec, ESMF_TimeDec2

! !DESCRIPTION:
!     This interface overloads the - operator for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(-)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeDiff

! !DESCRIPTION:
!     This interface overloads the - operator for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.EQ.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeEQ

! !DESCRIPTION:
!     This interface overloads the .EQ. operator for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.NE.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeNE

! !DESCRIPTION:
!     This interface overloads the .NE. operator for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.LT.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeLT

! !DESCRIPTION:
!     This interface overloads the .LT. operator for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.GT.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeGT

! !DESCRIPTION:
!     This interface overloads the .GT. operator for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.LE.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeLE

! !DESCRIPTION:
!     This interface overloads the .LE. operator for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.GE.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeGE

! !DESCRIPTION:
!     This interface overloads the .GE. operator for the {\tt ESMF\_Time} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------

!==============================================================================

      contains

!==============================================================================
!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeGet - Get value in user-specified units

! !INTERFACE:
!      subroutine ESMF_TimeGet(time, YY, YRl, MM, DD, D, Dl, H, M, S, Sl, MS, &
!                              US, NS, d_, h_, m_, s_, ms_, us_, ns_, Sn, Sd, &
!                              dayOfYear, dayOfYear_r8, dayOfYear_intvl,      &
!                              timeString, rc)

recursive subroutine ESMF_TimeGet(time, YY, MM, DD, D, Dl, H, M, S, MS, &
                              Sn, Sd, &
                              dayOfYear, dayOfYear_r8, dayOfYear_intvl,      &
                              timeString, rc)

! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time
      integer, intent(out), optional :: YY
!      integer(ESMF_KIND_I8), intent(out), optional :: YRl
      integer, intent(out), optional :: MM
      integer, intent(out), optional :: DD
      integer, intent(out), optional :: D
      integer(ESMF_KIND_I8), intent(out), optional :: Dl
      integer, intent(out), optional :: H
      integer, intent(out), optional :: M
      integer, intent(out), optional :: S
!      integer(ESMF_KIND_I8), intent(out), optional :: Sl
      integer, intent(out), optional :: MS
!      integer, intent(out), optional :: US
!      integer, intent(out), optional :: NS
!      double precision, intent(out), optional :: d_
!      double precision, intent(out), optional :: h_
!      double precision, intent(out), optional :: m_
!      double precision, intent(out), optional :: s_
!      double precision, intent(out), optional :: ms_
!      double precision, intent(out), optional :: us_
!      double precision, intent(out), optional :: ns_
      integer, intent(out), optional :: Sn
      integer, intent(out), optional :: Sd
      integer, intent(out), optional :: dayOfYear
      real(ESMF_KIND_R8), intent(out), optional :: dayOfYear_r8
      character (len=*), intent(out), optional :: timeString
      type(ESMF_TimeInterval), intent(out), optional :: dayOfYear_intvl
      integer, intent(out), optional :: rc


! !DESCRIPTION:
!     Get the value of the {\tt ESMF\_Time} in units specified by the user
!     via F90 optional arguments.
!
!     Time manager represents and manipulates time internally with integers
!     to maintain precision. Hence, user-specified floating point values are
!     converted internally from integers.
!
!     See {\tt ../include/ESMC\_BaseTime.h and ../include/ESMC\_Time.h} for
!     complete description.
!
!     The arguments are:
!     \begin{description}
!     \item[time]
!          The object instance to query
!     \item[{[YY]}]
!          Integer year CCYR (>= 32-bit)
!     \item[{[YRl]}]
!          Integer year CCYR (large, >= 64-bit)
!     \item[{[MM]}]
!          Integer month 1-12
!     \item[{[DD]}]
!          Integer day of the month 1-31
!     \item[{[D]}]
!          Integer Julian days (>= 32-bit)
!     \item[{[Dl]}]
!          Integer Julian days (large, >= 64-bit)
!     \item[{[H]}]
!          Integer hours
!     \item[{[M]}]
!          Integer minutes
!     \item[{[S]}]
!          Integer seconds (>= 32-bit)
!     \item[{[Sl]}]
!          Integer seconds (large, >= 64-bit)
!     \item[{[MS]}]
!          Integer milliseconds
!     \item[{[US]}]
!          Integer microseconds
!     \item[{[NS]}]
!          Integer nanoseconds
!     \item[{[d\_]}]
!          Double precision days
!     \item[{[h\_]}]
!          Double precision hours
!     \item[{[m\_]}]
!          Double precision minutes
!     \item[{[s\_]}]
!          Double precision seconds
!     \item[{[ms\_]}]
!          Double precision milliseconds
!     \item[{[us\_]}]
!          Double precision microseconds
!     \item[{[ns\_]}]
!          Double precision nanoseconds
!     \item[{[Sn]}]
!          Integer fractional seconds - numerator
!     \item[{[Sd]}]
!          Integer fractional seconds - denominator
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMG2.1, TMG2.5.1, TMG2.5.6
!EOP
      type(ESMF_TimeInterval) :: day_step
      integer :: ierr
      TYPE(ESMF_Time) :: begofyear
      TYPE(ESMF_TimeInterval) :: difftobegofyear
      INTEGER :: year, month, dayofmonth, hour, minute, second
      INTEGER :: i
      INTEGER(ESMF_KIND_I8) :: cnt

      ierr = ESMF_SUCCESS

      IF ( PRESENT( YY ) ) THEN
        YY = time%YR
      ENDIF
      IF ( PRESENT( MM ) ) THEN
        CALL timegetmonth( time, MM )
      ENDIF
      IF ( PRESENT( DD ) ) THEN
        CALL timegetdayofmonth( time, DD )
      ENDIF

      if (present(d) .or. present(dl)) then
        cnt = 0
        do i = 0,time%yr-1
           cnt = cnt + ndaysinyear(i,time%calendar%type)
        enddo
        do i = time%yr,-1
           cnt = cnt - ndaysinyear(i,time%calendar%type)
        enddo
        call timegetmonth(time,month)
        do i = 1,month-1
           cnt = cnt + ndaysinmonth(time%yr,i,time%calendar%type)
        enddo
        call timegetdayofmonth( time, dayofmonth)
        cnt = cnt + dayofmonth
        if (present(d)) then
          d = cnt
        endif
        if (present(dl)) then
          dl = cnt
        endif
      endif
!
!$$$ push HMS down into ESMF_BaseTime
      IF ( PRESENT( H ) ) THEN
        H = mod( time%basetime%S, SECONDS_PER_DAY ) / SECONDS_PER_HOUR
      ENDIF
      IF ( PRESENT( M ) ) THEN
        M = mod( time%basetime%S, SECONDS_PER_HOUR) / SECONDS_PER_MINUTE
      ENDIF
      IF ( PRESENT( S ) ) THEN
        S = mod( time%basetime%S, SECONDS_PER_MINUTE )
      ENDIF

      IF ( PRESENT( S ) .AND. PRESENT( DD ) ) THEN
        IF ( ( .NOT. PRESENT( H ) ) .AND. ( .NOT. PRESENT( M ) ) ) THEN
          S = mod( time%basetime%S, SECONDS_PER_DAY )
        ENDIF
      ENDIF
      IF ( PRESENT( MS ) ) THEN
        IF ( time%basetime%Sd /= 0 ) THEN
          MS = NINT( ( time%basetime%Sn*1.0D0 / time%basetime%Sd*1.0D0 ) * 1000.0D0 )
        ELSE
          MS = 0
        ENDIF
      ENDIF
      IF ( PRESENT( Sd ) .AND. PRESENT( Sn ) ) THEN
        Sd = time%basetime%Sd
        Sn = time%basetime%Sn
      ENDIF
      IF ( PRESENT( dayOfYear ) ) THEN
        CALL ESMF_TimeGetDayOfYear( time, dayOfYear, rc=ierr )
      ENDIF
      IF ( PRESENT( timeString ) ) THEN
        ! This duplication for YMD is an optimization that avoids calling
        ! timegetmonth() and timegetdayofmonth() when it is not needed.
        year = time%YR
        CALL timegetmonth( time, month )
        CALL timegetdayofmonth( time, dayofmonth )
!$$$ push HMS down into ESMF_BaseTime
        hour = mod( time%basetime%S, SECONDS_PER_DAY ) / SECONDS_PER_HOUR
        minute = mod( time%basetime%S, SECONDS_PER_HOUR) / SECONDS_PER_MINUTE
        second = mod( time%basetime%S, SECONDS_PER_MINUTE )
        CALL ESMFold_TimeGetString( year, month, dayofmonth, &
                                    hour, minute, second, timeString )
      ENDIF
      IF ( PRESENT( dayOfYear_intvl ) ) THEN
        year = time%YR
        CALL ESMF_TimeSet( begofyear, yy=year, mm=1, dd=1, s=0, &
                           calendar=time%calendar, rc=ierr )
        IF ( ierr == ESMF_FAILURE)THEN
           rc = ierr
           RETURN
        END IF
        dayOfYear_intvl = time - begofyear
      ENDIF
      IF ( PRESENT( dayOfYear_r8) ) THEN
        year = time%YR
        CALL ESMF_TimeSet( begofyear, yy=year, mm=1, dd=1, s=0, &
                           calendar=time%calendar, rc=ierr )
        IF ( ierr == ESMF_FAILURE)THEN
           rc = ierr
           RETURN
        END IF
        CALL ESMF_TimeIntervalSet( day_step, d=1, s=0, rc=ierr )
        IF ( ierr == ESMF_FAILURE)THEN
           rc = ierr
           RETURN
        END IF
        difftobegofyear = time - begofyear + day_step
        CALL ESMF_TimeIntervalGet( difftobegofyear, d_r8=dayOfYear_r8, rc=ierr )
        IF ( ierr == ESMF_FAILURE)THEN
           rc = ierr
           RETURN
        END IF
      ENDIF

      IF ( PRESENT( rc ) ) THEN
        rc = ierr
      ENDIF

      end subroutine ESMF_TimeGet

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeSet - Initialize via user-specified unit set

! !INTERFACE:
!      subroutine ESMF_TimeSet(time, YY, YRl, MM, DD, D, Dl, H, M, S, Sl, &
!                              MS, US, NS, d_, h_, m_, s_, ms_, us_, ns_, &
!                              Sn, Sd, calendar, calkindflag, rc)

      subroutine ESMF_TimeSet(time, YY, MM, DD, D, Dl, H, M, S, &
                              MS, &
                              Sn, Sd, calendar, calkindflag, rc)

! !ARGUMENTS:
      type(ESMF_Time), intent(inout) :: time
      integer, intent(in), optional :: YY
!      integer(ESMF_KIND_I8), intent(in), optional :: YRl
      integer, intent(in), optional :: MM
      integer, intent(in), optional :: DD
      integer, intent(in), optional :: D
      integer(ESMF_KIND_I8), intent(in), optional :: Dl
      integer, intent(in), optional :: H
      integer, intent(in), optional :: M
      integer, intent(in), optional :: S
!      integer(ESMF_KIND_I8), intent(in), optional :: Sl
      integer, intent(in), optional :: MS
!      integer, intent(in), optional :: US
!      integer, intent(in), optional :: NS
!      double precision, intent(in), optional :: d_
!      double precision, intent(in), optional :: h_
!      double precision, intent(in), optional :: m_
!      double precision, intent(in), optional :: s_
!      double precision, intent(in), optional :: ms_
!      double precision, intent(in), optional :: us_
!      double precision, intent(in), optional :: ns_
      integer, intent(in), optional :: Sn
      integer, intent(in), optional :: Sd
      type(ESMF_Calendar), intent(in), target, optional :: calendar
      type(ESMF_CalKind_Flag), intent(in),  optional :: calkindflag
      integer, intent(out), optional :: rc

      ! locals
      INTEGER :: ierr
      logical :: dset

! !DESCRIPTION:
!     Initializes a {\tt ESMF\_Time} with a set of user-specified units
!     via F90 optional arguments.
!
!     Time manager represents and manipulates time internally with integers
!     to maintain precision. Hence, user-specified floating point values are
!     converted internally to integers.
!
!     See {\tt ../include/ESMC\_BaseTime.h and ../include/ESMC\_Time.h} for
!     complete description.
!
!     The arguments are:
!     \begin{description}
!     \item[time]
!          The object instance to initialize
!     \item[{[YY]}]
!          Integer year CCYR (>= 32-bit)
!     \item[{[YRl]}]
!          Integer year CCYR (large, >= 64-bit)
!     \item[{[MM]}]
!          Integer month 1-12
!     \item[{[DD]}]
!          Integer day of the month 1-31
!     \item[{[D]}]
!          Integer Julian days (>= 32-bit)
!     \item[{[Dl]}]
!          Integer Julian days (large, >= 64-bit)
!     \item[{[H]}]
!          Integer hours
!     \item[{[M]}]
!          Integer minutes
!     \item[{[S]}]
!          Integer seconds (>= 32-bit)
!     \item[{[Sl]}]
!          Integer seconds (large, >= 64-bit)
!     \item[{[MS]}]
!          Integer milliseconds
!     \item[{[US]}]
!          Integer microseconds
!     \item[{[NS]}]
!          Integer nanoseconds
!     \item[{[d\_]}]
!          Double precision days
!     \item[{[h\_]}]
!          Double precision hours
!     \item[{[m\_]}]
!          Double precision minutes
!     \item[{[s\_]}]
!          Double precision seconds
!     \item[{[ms\_]}]
!          Double precision milliseconds
!     \item[{[us\_]}]
!          Double precision microseconds
!     \item[{[ns\_]}]
!          Double precision nanoseconds
!     \item[{[Sn]}]
!          Integer fractional seconds - numerator
!     \item[{[Sd]}]
!          Integer fractional seconds - denominator
!     \item[{[cal]}]
!          Associated {\tt Calendar}
!     \item[{[tz]}]
!          Associated timezone (hours offset from GMT, e.g. EST = -5)
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMGn.n.n
!EOP
!  PRINT *,'DEBUG:  BEGIN ESMF_TimeSet()'
!$$$ push this down into ESMF_BaseTime constructor

      IF ( PRESENT( rc ) ) then
        rc = ESMF_FAILURE
      ENDIF

      time%YR = 0
      time%basetime%S  = 0
      time%basetime%Sn = 0
      time%basetime%Sd = 0

      IF ( PRESENT(calendar) )THEN
!  PRINT *,'DEBUG:  ESMF_TimeSet():  using passed-in calendar'
        IF ( .not. ESMF_CalendarInitialized( calendar ) )THEN
           call wrf_error_fatal( "Error:: ESMF_CalendarCreate not "// &
                                 "called on input Calendar")
        END IF
!        call flush(6)
!        write(6,*) 'tcx1 ESMF_TimeSet point to calendar'
!        call flush(6)
        time%Calendar => calendar
      ELSE
!  PRINT *,'DEBUG:  ESMF_TimeSet():  using default calendar'
! for the sake of WRF, check ESMF_IsInitialized, revised by Juanxiong He
        IF ( .not. ESMF_IsInitialized() )THEN
           call wrf_error_fatal( "Error:: ESMF_Initialize not called")
        END IF
!        IF ( .not. ESMF_CalendarInitialized( defaultCal ) )THEN
!           call wrf_error_fatal( "Error:: ESMF_Initialize not called")
!        END IF
        if (present(calkindflag)) then
!           write(6,*) 'tcx2 ESMF_TimeSet point to calendarkindflag',calkindflag%caltype
!           call flush(6)
           if (calkindflag%caltype == ESMF_CALKIND_GREGORIAN%caltype) then
              time%Calendar => gregorianCal
           elseif (calkindflag%caltype == ESMF_CALKIND_NOLEAP%caltype) then
              time%Calendar => noleapCal
           else
              call wrf_error_fatal( "Error:: ESMF_TimeSet invalid calkindflag")
           endif
        else
!           write(6,*) 'tcx3 ESMF_TimeSet point to defaultcal'
!           call flush(6)
           time%Calendar => defaultCal
        endif
      END IF
!      write(6,*) 'tcxn ESMF_TimeSet ',ESMF_CALKIND_NOLEAP%caltype
!      call flush(6)
!      write(6,*) 'tcxg ESMF_TimeSet ',ESMF_CALKIND_GREGORIAN%caltype
!      call flush(6)
!      write(6,*) 'tcxt ESMF_TimeSet ',time%calendar%type%caltype
!      call flush(6)

      dset = .false.
      if (present(D)) then
         if (present(Dl)) CALL wrf_error_fatal( 'ESMF_TimeSet:  D and Dl not both valid')
         time%basetime%s = SECONDS_PER_DAY * INT(D-1,ESMF_KIND_I8)
         dset=.true.
      elseif (present(Dl)) then
         time%basetime%s = SECONDS_PER_DAY * Dl-1_ESMF_KIND_I8
         dset=.true.
      endif

      IF ( PRESENT( YY ) ) THEN
!  PRINT *,'DEBUG:  ESMF_TimeSet():  YY = ',YY
        if (dset) CALL wrf_error_fatal( 'ESMF_TimeSet:  D or DL and YY,MM,DD not both valid')
        time%YR = YY
      ENDIF
      IF ( PRESENT( MM ) ) THEN
        if (dset) CALL wrf_error_fatal( 'ESMF_TimeSet:  D or DL and YY,MM,DD not both valid')
!  PRINT *,'DEBUG:  ESMF_TimeSet():  MM = ',MM
        CALL timeaddmonths( time, MM, ierr )
        IF ( ierr == ESMF_FAILURE ) THEN
          IF ( PRESENT( rc ) ) THEN
            rc = ESMF_FAILURE
            RETURN
          ELSE
            CALL wrf_error_fatal( 'ESMF_TimeSet:  MM out of range' )
          ENDIF
        ENDIF
!  PRINT *,'DEBUG:  ESMF_TimeSet():  back from timeaddmonths'
      ENDIF
      IF ( PRESENT( DD ) ) THEN
        if (dset) CALL wrf_error_fatal( 'ESMF_TimeSet:  D or DL and YY,MM,DD not both valid')
!$$$ no check for DD in range of days of month MM yet
!$$$ Must separate D and DD for correct interface!
!  PRINT *,'DEBUG:  ESMF_TimeSet():  DD = ',DD
        time%basetime%S = time%basetime%S + &
          ( SECONDS_PER_DAY * INT( (DD-1), ESMF_KIND_I8 ) )
      ENDIF
!$$$ push H,M,S,Sn,Sd,MS down into ESMF_BaseTime constructor
      IF ( PRESENT( H ) ) THEN
!  PRINT *,'DEBUG:  ESMF_TimeSet():  H = ',H
        time%basetime%S = time%basetime%S + &
          ( SECONDS_PER_HOUR * INT( H, ESMF_KIND_I8 ) )
      ENDIF
      IF ( PRESENT( M ) ) THEN
!  PRINT *,'DEBUG:  ESMF_TimeSet():  M = ',M
        time%basetime%S = time%basetime%S + &
          ( SECONDS_PER_MINUTE * INT( M, ESMF_KIND_I8 ) )
      ENDIF
      IF ( PRESENT( S ) ) THEN
!  PRINT *,'DEBUG:  ESMF_TimeSet():  S = ',S
        time%basetime%S = time%basetime%S + &
          INT( S, ESMF_KIND_I8 )
      ENDIF
      IF ( PRESENT( Sn ) .AND. ( .NOT. PRESENT( Sd ) ) ) THEN
        CALL wrf_error_fatal( &
          "ESMF_TimeSet:  Must specify Sd if Sn is specified")
      ENDIF
      IF ( PRESENT( Sd ) .AND. PRESENT( MS ) ) THEN
        CALL wrf_error_fatal( &
          "ESMF_TimeSet:  Must not specify both Sd and MS")
      ENDIF
      time%basetime%Sn = 0
      time%basetime%Sd = 0
      IF ( PRESENT( MS ) ) THEN
!  PRINT *,'DEBUG:  ESMF_TimeSet():  MS = ',MS
        time%basetime%Sn = MS
        time%basetime%Sd = 1000_ESMF_KIND_I8
      ELSE IF ( PRESENT( Sd ) ) THEN
!  PRINT *,'DEBUG:  ESMF_TimeSet():  Sd = ',Sd
        time%basetime%Sd = Sd
        IF ( PRESENT( Sn ) ) THEN
!  PRINT *,'DEBUG:  ESMF_TimeSet():  Sn = ',Sn
          time%basetime%Sn = Sn
        ENDIF
      ENDIF

!  PRINT *,'DEBUG:  ESMF_TimeSet():  calling normalize_time()'
!$$$DEBUG
!IF ( time%basetime%Sd > 0 ) THEN
!  PRINT *,'DEBUG ESMF_TimeSet() before normalize:  S,Sn,Sd = ', &
!    time%basetime%S, time%basetime%Sn, time%basetime%Sd
!ENDIF
!$$$END DEBUG
      CALL normalize_time( time )
!$$$DEBUG
!IF ( time%basetime%Sd > 0 ) THEN
!  PRINT *,'DEBUG ESMF_TimeSet() after normalize:  S,Sn,Sd = ', &
!    time%basetime%S, time%basetime%Sn, time%basetime%Sd
!ENDIF
!$$$END DEBUG

!  PRINT *,'DEBUG:  ESMF_TimeSet():  back from normalize_time()'
      IF ( PRESENT( rc ) ) THEN
        rc = ESMF_SUCCESS
      ENDIF

      end subroutine ESMF_TimeSet

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMFold_TimeGetString - Get time instant value in string format

! !INTERFACE:
      subroutine ESMFold_TimeGetString( year, month, dayofmonth, &
                                        hour, minute, second, TimeString )

! !ARGUMENTS:
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: dayofmonth
      integer, intent(in) :: hour
      integer, intent(in) :: minute
      integer, intent(in) :: second
      character*(*), intent(out) :: TimeString
      character*(256) :: TimeFormatString
! !DESCRIPTION:
!     Convert {\tt ESMF\_Time}'s value into ISO 8601 format YYYY-MM-DDThh:mm:ss
!
!     The arguments are:
!     \begin{description}
!     \item[time]
!          The object instance to convert
!     \item[TimeString]
!          The string to return
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMG2.4.7
!EOP

!PRINT *,'DEBUG:  ESMF_TimePrint():  YR,S,Sn,Sd = ',time%YR,time%basetime%S,time%basetime%Sn,time%basetime%Sd
!PRINT *,'DEBUG:  ESMF_TimePrint():  year = ',year
!PRINT *,'DEBUG:  ESMF_TimePrint():  month, dayofmonth = ',month,dayofmonth
!PRINT *,'DEBUG:  ESMF_TimePrint():  hour = ',hour
!PRINT *,'DEBUG:  ESMF_TimePrint():  minute = ',minute
!PRINT *,'DEBUG:  ESMF_TimePrint():  second = ',second

!$$$here...  add negative sign for YR<0
!$$$here...  add Sn, Sd ??
      write(TimeFormatString,FMT="(A,I4.4,A,I4.4,A)") &
           "(I", yearWidth, ".", yearWidth, ",'-',I2.2,'-',I2.2,'_',I2.2,':',I2.2,':',I2.2)"
      write(TimeString,FMT=TimeFormatString) year,month,dayofmonth,hour,minute,second

      end subroutine ESMFold_TimeGetString

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeGetDayOfYearInteger - Get time instant's day of the year as an integer value
!
! !INTERFACE:
      subroutine ESMF_TimeGetDayOfYearInteger(time, DayOfYear, rc)
!
! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time
      integer, intent(out) :: DayOfYear
      integer, intent(out), optional :: rc
!
! !DESCRIPTION:
!     Get the day of the year the given {\tt ESMF\_Time} instant falls on
!     (1-365).  Returned as an integer value
!
!     The arguments are:
!     \begin{description}
!     \item[time]
!          The object instance to query
!     \item[DayOfYear]
!          The {\tt ESMF\_Time} instant's day of the year (1-365)
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!EOP
      ! requires that time be normalized
!$$$ bug when Sn>0?  test
!$$$ add tests
      DayOfYear = ( time%basetime%S / SECONDS_PER_DAY ) + 1
      IF ( PRESENT( rc ) ) rc = ESMF_SUCCESS
      end subroutine ESMF_TimeGetDayOfYearInteger

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeInc - Increment time instant with a time interval
!
! !INTERFACE:
   function ESMF_TimeInc(time, timeinterval)
!
! !RETURN VALUE:
      type(ESMF_Time) :: ESMF_TimeInc
!
! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time
      type(ESMF_TimeInterval), intent(in) :: timeinterval
! !LOCAL:
      INTEGER :: year,month,day,sec,nmon,nyr,mpyi4
!
! !DESCRIPTION:
!     Increment {\tt ESMF\_Time} instant with a {\tt ESMF\_TimeInterval},
!     return resulting {\tt ESMF\_Time} instant
!
!     Maps overloaded (+) operator interface function to
!     {\tt ESMF\_BaseTime} base class
!
!     The arguments are:
!     \begin{description}
!     \item[time]
!          The given {\tt ESMF\_Time} to increment
!     \item[timeinterval]
!          The {\tt ESMF\_TimeInterval} to add to the given {\tt ESMF\_Time}
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.4, TMG2.4.4, TMG2.4.5, TMG2.4.6, TMG5.1, TMG5.2, TMG7.2
!EOP

      mpyi4 = MONTHS_PER_YEAR

      ! copy ESMF_Time specific properties (e.g. calendar, timezone)

      ESMF_TimeInc = time
!      write(6,*) 'tcx timeinc1 ',ESMF_TimeInc%yr,ESMF_TimeInc%basetime%s
      CALL normalize_time( ESMF_TimeInc )

!      write(6,*) 'tcx timeint ',timeinterval%yr,timeinterval%mm,timeinterval%basetime%s

      ! add years and months by manually forcing incremental years then adjusting the day of
      ! the month at the end if it's greater than the number of days in the month
      ! esmf seems to do exactly this based on testing

      nmon = timeinterval%mm
      nyr  = timeinterval%yr
      if (abs(nmon) > 0 .or. abs(nyr) > 0) then
         call ESMF_TimeGet(ESMF_TimeInc,yy=year,mm=month,dd=day,s=sec)
!         write(6,*) 'tcx timeinc mon1 ',year,month,day,sec,nyr,nmon
         year  = year  + nyr
         month = month + nmon
         do while (month > MONTHS_PER_YEAR)
            month = month - mpyi4
            year = year + 1
         enddo
         do while (month < 1)
            month = month + mpyi4
            year = year - 1
         enddo
!         write(6,*) 'tcx timeinc mon2 ',year,month,day,sec
         day = min(day,ndaysinmonth(year,month,ESMF_TimeInc%calendar%type))
         call ESMF_TimeSet(ESMF_TimeInc,yy=year,mm=month,dd=day,s=sec,calkindflag=time%calendar%type)
         call ESMF_TimeGet(ESMF_TimeInc,yy=year,mm=month,dd=day,s=sec)
!         write(6,*) 'tcx timeinc mon3 ',nmon,year,month,day,sec
      endif

      ! finally add seconds

!      write(6,*) 'tcx timeinc sec ',ESMF_TimeInc%basetime%s,timeinterval%basetime%s
      ESMF_TimeInc%basetime = ESMF_TimeInc%basetime + timeinterval%basetime

      ! and normalize

!      write(6,*) 'tcx timeinc2p ',ESMF_TimeInc%yr,ESMF_TimeInc%basetime%s

      CALL normalize_time( ESMF_TimeInc )

!      write(6,*) 'tcx timeinc2 ',ESMF_TimeInc%yr,ESMF_TimeInc%basetime%s

   end function ESMF_TimeInc

! this is added for certain compilers that don't deal with commutativity

   function ESMF_TimeInc2(timeinterval, time)
      type(ESMF_Time) :: ESMF_TimeInc2
      type(ESMF_Time), intent(in) :: time
      type(ESMF_TimeInterval), intent(in) :: timeinterval
      ESMF_TimeInc2 = ESMF_TimeInc( time, timeinterval )
   end function ESMF_TimeInc2

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeDec - Decrement time instant with a time interval
!
! !INTERFACE:
   function ESMF_TimeDec(time, timeinterval)
!
! !RETURN VALUE:
      type(ESMF_Time) :: ESMF_TimeDec
!
! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time
      type(ESMF_TimeInterval), intent(in) :: timeinterval
! !LOCAL:
      TYPE (ESMF_TimeInterval)  :: neginterval

! !DESCRIPTION:
!     Decrement {\tt ESMF\_Time} instant with a {\tt ESMF\_TimeInterval},
!     return resulting {\tt ESMF\_Time} instant
!
!     Maps overloaded (-) operator interface function to
!     {\tt ESMF\_BaseTime} base class
!
!     The arguments are:
!     \begin{description}
!     \item[time]
!          The given {\tt ESMF\_Time} to decrement
!     \item[timeinterval]
!          The {\tt ESMF\_TimeInterval} to subtract from the given
!          {\tt ESMF\_Time}
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.4, TMG2.4.4, TMG2.4.5, TMG2.4.6, TMG5.1, TMG5.2, TMG7.2
!EOP

      ESMF_TimeDec = time

      neginterval = timeinterval
!$$$push this down into a unary negation operator on TimeInterval
      neginterval%basetime%S = -neginterval%basetime%S
      neginterval%basetime%Sn = -neginterval%basetime%Sn
      neginterval%YR = -neginterval%YR
      neginterval%MM = -neginterval%MM
      ESMF_TimeDec = time + neginterval

   end function ESMF_TimeDec

!
! this is added for certain compilers that don't deal with commutativity
!
   function ESMF_TimeDec2(timeinterval, time)
      type(ESMF_Time) :: ESMF_TimeDec2
      type(ESMF_Time), intent(in) :: time
      type(ESMF_TimeInterval), intent(in) :: timeinterval
      ESMF_TimeDec2 = ESMF_TimeDec( time, timeinterval )
   end function ESMF_TimeDec2
!
!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeDiff - Return the difference between two time instants
!
! !INTERFACE:
   function ESMF_TimeDiff(time1, time2)
!
! !RETURN VALUE:
      type(ESMF_TimeInterval) :: ESMF_TimeDiff
!
! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time1
      type(ESMF_Time), intent(in) :: time2
! !LOCAL:
      TYPE(ESMF_BaseTime) :: cmptime, zerotime
      integer :: yr
      integer :: y1,m1,d1,s1,y2,m2,d2,s2
      integer :: rc

! !DESCRIPTION:
!     Return the {\tt ESMF\_TimeInterval} difference between two
!     {\tt ESMF\_Time} instants, time1 - time2
!
!     Maps overloaded (-) operator interface function to
!     {\tt ESMF\_BaseTime} base class
!
!     The arguments are:
!     \begin{description}
!     \item[time1]
!          The first {\tt ESMF\_Time} instant
!     \item[time2]
!          The second {\tt ESMF\_Time} instant
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.4, TMG2.4.4, TMG2.4.5, TMG2.4.6, TMG5.1, TMG5.2, TMG7.2
!EOP

      CALL ESMF_TimeIntervalSet( ESMF_TimeDiff, rc=rc )

      ESMF_TimeDiff%StartTime = time2
      ESMF_TimeDiff%StartTime_set = .true.

!      write(6,*) 'tcx timediff1 ',time2%yr,time2%basetime%s,time2%calendar%type%caltype
!      write(6,*) 'tcx timediff2 ',time1%yr,time1%basetime%s,time1%calendar%type%caltype

      call ESMF_TimeGet(time2,yy=y2,mm=m2,dd=d2,s=s2)
      call ESMF_TimeGet(time1,yy=y1,mm=m1,dd=d1,s=s1)

      ! Can either be yr/month based diff if diff is only in year and month
      ! or absolute seconds if diff in day/seconds as well

      if (d1 == d2 .and. s1 == s2) then
!         write(6,*) 'tcx timedifft ym'
         ESMF_TimeDiff%YR = y1 - y2
         ESMF_TimeDiff%MM = m1 - m2
         cmptime%S  = 0
         cmptime%Sn = 0
         cmptime%Sd = 0
         ESMF_TimeDiff%basetime = cmptime
      else
!         write(6,*) 'tcx timedifft sec'
         ESMF_TimeDiff%YR = 0
         ESMF_TimeDiff%MM = 0
         ESMF_TimeDiff%basetime = time1%basetime - time2%basetime
         IF ( time1%YR > time2%YR ) THEN
            DO yr = time2%YR, ( time1%YR - 1 )
!              write(6,*) 'tcx timediff3 ',yr,nsecondsinyear(yr,time2%calendar%type)
              ESMF_TimeDiff%basetime%S = ESMF_TimeDiff%basetime%S + nsecondsinyear(yr,time2%calendar%type)
            ENDDO
         ELSE IF ( time2%YR > time1%YR ) THEN
            DO yr = time1%YR, ( time2%YR - 1 )
!              write(6,*) 'tcx timediff4 ',yr,nsecondsinyear(yr,time2%calendar%type)
              ESMF_TimeDiff%basetime%S = ESMF_TimeDiff%basetime%S - nsecondsinyear(yr,time2%calendar%type)
            ENDDO
         ENDIF
      endif

!      write(6,*) 'tcx timediff5 ',ESMF_TimeDiff%YR, ESMF_TimeDiff%MM, ESMF_TimeDiff%basetime%s

      CALL normalize_timeint( ESMF_TimeDiff )

!      write(6,*) 'tcx timediff6 ',ESMF_TimeDiff%YR, ESMF_TimeDiff%MM, ESMF_TimeDiff%basetime%s

   end function ESMF_TimeDiff

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeEQ - Compare two times for equality
!
! !INTERFACE:
      function ESMF_TimeEQ(time1, time2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeEQ
!
! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time1
      type(ESMF_Time), intent(in) :: time2
!
! !DESCRIPTION:
!     Return true if both given {\tt ESMF\_Time} instants are equal, false
!     otherwise.  Maps overloaded (==) operator interface function to
!     {\tt ESMF\_BaseTime} base class.
!
!     The arguments are:
!     \begin{description}
!     \item[time1]
!          First time instant to compare
!     \item[time2]
!          Second time instant to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      integer :: res

      call timecmp(time1,time2,res)
      ESMF_TimeEQ = (res .EQ. 0)

      end function ESMF_TimeEQ

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeNE - Compare two times for non-equality
!
! !INTERFACE:
      function ESMF_TimeNE(time1, time2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeNE
!
! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time1
      type(ESMF_Time), intent(in) :: time2

! !DESCRIPTION:
!     Return true if both given {\tt ESMF\_Time} instants are not equal, false
!     otherwise.  Maps overloaded (/=) operator interface function to
!     {\tt ESMF\_BaseTime} base class.
!
!     The arguments are:
!     \begin{description}
!     \item[time1]
!          First time instant to compare
!     \item[time2]
!          Second time instant to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      integer :: res

      call timecmp(time1,time2,res)
      ESMF_TimeNE = (res .NE. 0)

      end function ESMF_TimeNE

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeLT - Time instant 1 less than time instant 2 ?
!
! !INTERFACE:
      function ESMF_TimeLT(time1, time2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeLT
!
! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time1
      type(ESMF_Time), intent(in) :: time2
!
! !DESCRIPTION:
!     Return true if first {\tt ESMF\_Time} instant is less than second
!     {\tt ESMF\_Time} instant, false otherwise.  Maps overloaded (<)
!     operator interface function to {\tt ESMF\_BaseTime} base class.
!
!     The arguments are:
!     \begin{description}
!     \item[time1]
!          First time instant to compare
!     \item[time2]
!          Second time instant to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      integer :: res

      call timecmp(time1,time2,res)
      ESMF_TimeLT = (res .LT. 0)

      end function ESMF_TimeLT

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeGT - Time instant 1 greater than time instant 2 ?
!
! !INTERFACE:
      function ESMF_TimeGT(time1, time2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeGT
!
! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time1
      type(ESMF_Time), intent(in) :: time2
!
! !DESCRIPTION:
!     Return true if first {\tt ESMF\_Time} instant is greater than second
!     {\tt ESMF\_Time} instant, false otherwise.  Maps overloaded (>) operator
!     interface function to {\tt ESMF\_BaseTime} base class.
!
!     The arguments are:
!     \begin{description}
!     \item[time1]
!          First time instant to compare
!     \item[time2]
!          Second time instant to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      integer :: res

      call timecmp(time1,time2,res)
      ESMF_TimeGT = (res .GT. 0)

      end function ESMF_TimeGT

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeLE - Time instant 1 less than or equal to time instant 2 ?
!
! !INTERFACE:
      function ESMF_TimeLE(time1, time2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeLE
!
! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time1
      type(ESMF_Time), intent(in) :: time2
!
! !DESCRIPTION:
!     Return true if first {\tt ESMF\_Time} instant is less than or equal to
!     second {\tt ESMF\_Time} instant, false otherwise.  Maps overloaded (<=)
!     operator interface function to {\tt ESMF\_BaseTime} base class.
!
!     The arguments are:
!     \begin{description}
!     \item[time1]
!          First time instant to compare
!     \item[time2]
!          Second time instant to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      integer :: res

      call timecmp(time1,time2,res)
      ESMF_TimeLE = (res .LE. 0)

      end function ESMF_TimeLE

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeGE - Time instant 1 greater than or equal to time instant 2 ?
!
! !INTERFACE:
      function ESMF_TimeGE(time1, time2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeGE
!
! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time1
      type(ESMF_Time), intent(in) :: time2
!
! !DESCRIPTION:
!     Return true if first {\tt ESMF\_Time} instant is greater than or equal to
!     second {\tt ESMF\_Time} instant, false otherwise.  Maps overloaded (>=)
!     operator interface function to {\tt ESMF\_BaseTime} base class.
!
!     The arguments are:
!     \begin{description}
!     \item[time1]
!          First time instant to compare
!     \item[time2]
!          Second time instant to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      integer :: res

      call timecmp(time1,time2,res)
      ESMF_TimeGE = (res .GE. 0)

      end function ESMF_TimeGE

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeCopy - Copy a time-instance

! !INTERFACE:
      subroutine ESMF_TimeCopy(timeout, timein)

! !ARGUMENTS:
      type(ESMF_Time), intent(out) :: timeout
      type(ESMF_Time), intent(in) :: timein

! !DESCRIPTION:
!     Copy a time-instance to a new instance.
!
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMGn.n.n
!EOP

      timeout%basetime = timein%basetime
      timeout%YR       = timein%YR
      timeout%Calendar => timein%Calendar
!tcx      timeout%Calendar = timein%Calendar
!      write(6,*) 'tcxa ESMF_TimeCopy'
!      call flush(6)
!      write(6,*) 'tcxb ESMF_TimeCopy',timein%calendar%type%caltype
!      call flush(6)
      timeout%Calendar = ESMF_CalendarCreate(calkindflag=timein%calendar%type)

      end subroutine ESMF_TimeCopy


!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimePrint - Print out a time instant's properties


! !INTERFACE:
      subroutine ESMF_TimePrint(time, options, rc)

! !ARGUMENTS:
      type(ESMF_Time), intent(in) :: time
      character (len=*), intent(in), optional :: options
      integer, intent(out), optional :: rc
      character (len=256) :: timestr

! !DESCRIPTION:
!     To support testing/debugging, print out a {\tt ESMF\_Time}'s
!     properties.
!
!     The arguments are:
!     \begin{description}
!     \item[time]
!          {\tt ESMF\_Time} instant to print out
!     \item[{[options]}]
!          Print options
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMGn.n.n
!EOP

      ! Quick hack to mimic ESMF 2.0.1
      ! Really should check value of options...
      IF ( PRESENT( options ) ) THEN
        CALL ESMF_TimeGet( time, timeString=timestr, rc=rc )
        timestr(11:11) = 'T'     ! ISO 8601 compatibility hack for debugging
        print *,' Time -----------------------------------'
        print *,' ',TRIM(timestr)
        print *,' end Time -------------------------------'
        print *
      ELSE
        call print_a_time (time)
      ENDIF

      end subroutine ESMF_TimePrint

!==============================================================================

SUBROUTINE print_a_time( time )
   IMPLICIT NONE
   type(ESMF_Time) time
   character*128 :: s
   integer rc
   CALL ESMF_TimeGet( time, timeString=s, rc=rc )
   print *,'Print a time|',TRIM(s),'|'
   write(0,*)'Print a time|',TRIM(s),'|'
   return
END SUBROUTINE print_a_time

!==============================================================================

SUBROUTINE timecmp(time1, time2, retval )
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: retval
!
! !ARGUMENTS:
  TYPE(ESMF_Time), INTENT(IN) :: time1
  TYPE(ESMF_Time), INTENT(IN) :: time2
  IF ( time1%YR .GT. time2%YR ) THEN ; retval = 1  ; RETURN ; ENDIF
  IF ( time1%YR .LT. time2%YR ) THEN ; retval = -1 ; RETURN ; ENDIF
  CALL seccmp( time1%basetime%S, time1%basetime%Sn, time1%basetime%Sd, &
               time2%basetime%S, time2%basetime%Sn, time2%basetime%Sd, &
	       retval )
END SUBROUTINE timecmp

!==============================================================================

SUBROUTINE normalize_time( time )
  ! A normalized time has time%basetime >= 0, time%basetime less than the current
  ! year expressed as a timeInterval, and time%YR can take any value
  IMPLICIT NONE
  TYPE(ESMF_Time), INTENT(INOUT) :: time
!  INTEGER(ESMF_KIND_I8) :: nsecondsinyear
  ! locals
  TYPE(ESMF_BaseTime) :: cmptime, zerotime
  INTEGER :: rc
  LOGICAL :: done

  ! first, normalize basetime
  ! this will force abs(Sn) < Sd and ensure that signs of S and Sn match

  CALL normalize_basetime( time%basetime )

  ! next, underflow negative seconds into YEARS
  ! time%basetime must end up non-negative

  zerotime%S  = 0
  zerotime%Sn = 0
  zerotime%Sd = 0
  DO WHILE ( time%basetime < zerotime )
    time%YR = time%YR - 1
    cmptime%S  = nsecondsinyear( time%YR, time%calendar%type )
    cmptime%Sn = 0
    cmptime%Sd = 0
    time%basetime = time%basetime + cmptime
  ENDDO

  ! next, overflow seconds into YEARS
  done = .FALSE.
  DO WHILE ( .NOT. done )
    cmptime%S  = nsecondsinyear( time%YR, time%calendar%type )
    cmptime%Sn = 0
    cmptime%Sd = 0
    IF ( time%basetime >= cmptime ) THEN
      time%basetime = time%basetime - cmptime
      time%YR = time%YR + 1
    ELSE
      done = .TRUE.
    ENDIF
  ENDDO

END SUBROUTINE normalize_time

!==============================================================================

SUBROUTINE timegetmonth( time, MM )
  IMPLICIT NONE
  TYPE(ESMF_Time), INTENT(IN) :: time
  INTEGER, INTENT(OUT) :: MM
  ! locals

  mm = nmonthinyearsec(time%yr,time%basetime,time%calendar%type)

END SUBROUTINE timegetmonth

!==============================================================================
SUBROUTINE timegetdayofmonth( time, DD )
  IMPLICIT NONE
  TYPE(ESMF_Time), INTENT(IN) :: time
  INTEGER, INTENT(OUT) :: DD
  ! locals

  dd = ndayinyearsec(time%yr, time%basetime, time%calendar%type)

END SUBROUTINE timegetdayofmonth

!==============================================================================

! Increment Time by number of seconds between start of year and start
! of month MM.
! 1 <= MM <= 12
! Time is NOT normalized.
SUBROUTINE timeaddmonths( time, MM, ierr )
  IMPLICIT NONE
  TYPE(ESMF_Time), INTENT(INOUT) :: time
  INTEGER, INTENT(IN) :: MM
  INTEGER, INTENT(OUT) :: ierr
  ! locals
  INTEGER(ESMF_KIND_I8) :: isec

  ierr = ESMF_SUCCESS
  IF ( ( MM < 1 ) .OR. ( MM > MONTHS_PER_YEAR ) ) THEN
    CALL wrf_message( 'ERROR timeaddmonths():  MM out of range' )
    ierr = ESMF_FAILURE
    return
  ENDIF

  isec = nsecondsinyearmonth(time%yr,MM,time%calendar%type)
  time%basetime%s = time%basetime%s + isec

END SUBROUTINE timeaddmonths

!==============================================================================

! Increment Time by number of seconds between start of year and start
! of month MM.
! 1 <= MM <= 12
! Time is NOT normalized.
SUBROUTINE ESMF_setYearWidth( yearWidthIn )

    integer, intent(in) :: yearWidthIn

    yearWidth = yearWidthIn

END SUBROUTINE ESMF_setYearWidth

!==============================================================================
!==============================================================================
end module ESMF_TimeMod
