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
!     ESMF TimeInterval Module

module ESMF_TimeIntervalMod

!
!==============================================================================
!
! This file contains the TimeInterval class definition and all TimeInterval
! class methods.
!
!------------------------------------------------------------------------------
! INCLUDES
#include <ESMF_TimeMgr.inc>
!
!===============================================================================
!BOPI
! !MODULE: ESMF_TimeIntervalMod
!
! !DESCRIPTION:
! Part of Time Manager F90 API wrapper of C++ implemenation
!
! Defines F90 wrapper entry points for corresponding
! C++ implementaion of class {\tt ESMC\_TimeInterval}
!
! See {\tt ../include/ESMC\_TimeInterval.h} for complete description
!
!------------------------------------------------------------------------------
! !USES:
      ! inherit from ESMF base class
      use ESMF_BaseMod

      ! inherit from base time class
      use ESMF_BaseTimeMod

      ! associated derived types
      use ESMF_FractionMod, only : ESMF_Fraction
      use ESMF_CalendarMod
      use ESMF_ShrTimeMod, only : ESMF_Time

      implicit none
!
!------------------------------------------------------------------------------
! !PRIVATE TYPES:
      private
!------------------------------------------------------------------------------
!     ! ESMF_TimeInterval
!
!     ! F90 class type to match C++ TimeInterval class in size only;
!     !  all dereferencing within class is performed by C++ implementation

      type ESMF_TimeInterval
        ! time interval is expressed as basetime
        type(ESMF_BaseTime) :: basetime  ! inherit base class
        ! Relative year and month fields support monthly or yearly time
        ! intervals.  Many operations are undefined when these fields are
        ! non-zero!
        INTEGER :: YR                    ! relative year
        INTEGER :: MM                    ! relative month
        logical :: starttime_set           ! reference time set
        type(ESMF_Time) :: starttime       ! reference time
      end type

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
      public ESMF_TimeInterval
!------------------------------------------------------------------------------
!
! for running WRF, add three subroutines or functions (WRFADDITION_TimeIntervalGet,
! ESMF_TimeIntervalDIVQuot, ESMF_TimeIntervalIsPositive), by jhe
! !PUBLIC MEMBER FUNCTIONS:
      public ESMF_TimeIntervalGet
      public ESMF_TimeIntervalSet
      public ESMF_TimeIntervalAbsValue
      public ESMF_TimeIntervalNegAbsValue
      public ESMF_TimeIntervalPrint
      public normalize_timeint

! Required inherited and overridden ESMF_Base class methods

!!!!!!!!! added by jhe
      public ESMF_TimeIntervalDIVQuot
      public ESMF_TimeIntervalIsPositive
!

! !PRIVATE MEMBER FUNCTIONS:

! overloaded operator functions

      public operator(/)
      private ESMF_TimeIntervalQuotI

      public operator(*)
      private ESMF_TimeIntervalProdI

! Inherited and overloaded from ESMF_BaseTime

      public operator(+)
      private ESMF_TimeIntervalSum

      public operator(-)
      private ESMF_TimeIntervalDiff

      public operator(.EQ.)
      private ESMF_TimeIntervalEQ

      public operator(.NE.)
      private ESMF_TimeIntervalNE

      public operator(.LT.)
      private ESMF_TimeIntervalLT

      public operator(.GT.)
      private ESMF_TimeIntervalGT

      public operator(.LE.)
      private ESMF_TimeIntervalLE

      public operator(.GE.)
      private ESMF_TimeIntervalGE
!EOPI

!------------------------------------------------------------------------------
! The following line turns the CVS identifier string into a printable variable.
      character(*), parameter, private :: version = &
      '$Id$'

!==============================================================================
!
! INTERFACE BLOCKS
!
!==============================================================================
!BOP
! !INTERFACE:
      interface operator(*)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeIntervalProdI

! !DESCRIPTION:
!     This interface overloads the * operator for the {\tt ESMF\_TimeInterval}
!     class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(/)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeIntervalQuotI

! !DESCRIPTION:
!     This interface overloads the / operator for the
!     {\tt ESMF\_TimeInterval} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(+)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeIntervalSum

! !DESCRIPTION:
!     This interface overloads the + operator for the
!     {\tt ESMF\_TimeInterval} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(-)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeIntervalDiff

! !DESCRIPTION:
!     This interface overloads the - operator for the
!     {\tt ESMF\_TimeInterval} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.EQ.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeIntervalEQ

! !DESCRIPTION:
!     This interface overloads the .EQ. operator for the
!     {\tt ESMF\_TimeInterval} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.NE.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeIntervalNE

! !DESCRIPTION:
!     This interface overloads the .NE. operator for the
!     {\tt ESMF\_TimeInterval} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.LT.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeIntervalLT

! !DESCRIPTION:
!     This interface overloads the .LT. operator for the
!     {\tt ESMF\_TimeInterval} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.GT.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeIntervalGT

! !DESCRIPTION:
!     This interface overloads the .GT. operator for the
!     {\tt ESMF\_TimeInterval} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.LE.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeIntervalLE

! !DESCRIPTION:
!     This interface overloads the .LE. operator for the
!     {\tt ESMF\_TimeInterval} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------
!BOP
! !INTERFACE:
      interface operator(.GE.)

! !PRIVATE MEMBER FUNCTIONS:
      module procedure ESMF_TimeIntervalGE

! !DESCRIPTION:
!     This interface overloads the .GE. operator for the
!     {\tt ESMF\_TimeInterval} class
!
!EOP
      end interface
!
!------------------------------------------------------------------------------

!==============================================================================

      contains

!==============================================================================
!
! Generic Get/Set routines which use F90 optional arguments
!
!---------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeIntervalGet - Get value in user-specified units

! !INTERFACE:
      subroutine ESMF_TimeIntervalGet(timeinterval, StartTimeIn, yy, mm, D, d_r8, S, S_i8, Sn, Sd, TimeString, rc )

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval
      type(ESMF_Time), optional, intent(in) :: StartTimeIn
      integer, intent(out), optional :: yy
      integer, intent(out), optional :: mm
      integer, intent(out), optional :: D
      real(ESMF_KIND_R8),   intent(out), optional :: d_r8
      integer(ESMF_KIND_I8),intent(out), optional :: S_i8
      integer, intent(out), optional :: S
      integer, intent(out), optional :: Sn
      integer, intent(out), optional :: Sd
      character*(*), optional, intent(out) :: TimeString
      integer, intent(out), optional :: rc


! !DESCRIPTION:
!     Get the value of the {\tt ESMF\_TimeInterval} in units specified by the
!     user via F90 optional arguments.
!
!     Time manager represents and manipulates time internally with integers
!     to maintain precision.  Hence, user-specified floating point values are
!     converted internally from integers.
!
!     See {\tt ../include/ESMC\_BaseTime.h} and
!     {\tt ../include/ESMC\_TimeInterval.h} for complete description.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval]
!          The object instance to query
!     \item[{[YY]}]
!          Integer years (>= 32-bit)
!     \item[{[YYl]}]
!          Integer years (large, >= 64-bit)
!     \item[{[MO]}]
!          Integer months (>= 32-bit)
!     \item[{[MOl]}]
!          Integer months (large, >= 64-bit)
!     \item[{[D]}]
!          Integer days (>= 32-bit)
!     \item[{[Dl]}]
!          Integer days (large, >= 64-bit)
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
!     TMG1.1
!EOP
      type(ESMF_Time) :: lstarttime
      logical :: lstarttime_set
      logical :: doyear
      INTEGER(ESMF_KIND_I8) :: seconds, secondsym, years
      INTEGER :: ierr
      INTEGER :: mpyi4, iyr,imo,mmon,nmon,mstart,ndays

      ierr = ESMF_FAILURE

      if (present(StartTimeIn)) then
         lstarttime_set = .true.
         lstarttime = StartTimeIn
      else
         lstarttime_set = timeinterval%StartTime_set
         lstarttime = timeinterval%StartTime
      endif


      CALL timeintchecknormalized( timeinterval,                &
                                   'ESMF_TimeIntervalGet arg1', &
                                    relative_interval=.true. )
      seconds = timeinterval%basetime%S
      years = timeinterval%YR

      secondsym = 0

      IF ( PRESENT( YY ) )THEN
        YY = years + timeinterval%MM / MONTHS_PER_YEAR
!        seconds = seconds - years * ( 365_ESMF_KIND_I8 * SECONDS_PER_DAY )
        IF ( PRESENT( MM ) )THEN
           mpyi4 = MONTHS_PER_YEAR
           MM = MOD( timeinterval%MM, mpyi4)
        else
           call wrf_error_fatal("ESMF_TimeIntervalGet: requires MM with YY")
        END IF
      ELSE IF ( PRESENT( MM ) )THEN
        MM = timeinterval%MM + years*12
      else if (lstarttime_set) then
        ! convert years and months to days carefully

        mpyi4 = MONTHS_PER_YEAR
        mmon = timeinterval%mm + timeinterval%yr*mpyi4
        mstart = nmonthinyearsec(lstarttime%yr,lstarttime%basetime,lstarttime%calendar%type)
!        write(6,*) 'tcxti1 ',mmon,lstarttime%yr,mstart,lstarttime%basetime%s

        iyr = lstarttime%yr
        if (mmon > 0) then
           imo = mstart-1  ! if adding months, start with this month after adding first +1
        else
           imo = mstart    ! if going backwards, start with last month after first -1
        endif
        nmon = 1
!        do nmon = 1,abs(mmon)
        do while (nmon <= abs(mmon))
           if (mmon > 0) then
              if (imo == 12 .and. (abs(mmon) - nmon) > 12) then
                 iyr = iyr + 1
                 nmon = nmon + 12
                 doyear = .true.
              else
                 imo = imo + 1
                 nmon = nmon + 1
                 doyear = .false.
              endif
           else
              if (imo == 1 .and. (abs(mmon) - nmon) > 12) then
                 iyr = iyr - 1
                 nmon = nmon + 12
                 doyear = .true.
              else
                 imo = imo - 1
                 nmon = nmon + 1
                 doyear = .false.
              endif
           endif

           do while (imo > 12)
              imo = imo - 12
              iyr = iyr + 1
           enddo
           do while (imo < 1)
              imo = imo + 12
              iyr = iyr - 1
           enddo

           if (doyear) then
              ndays = ndaysinyear(iyr,lstarttime%calendar%type)
           else
              ndays = ndaysinmonth(iyr,imo,lstarttime%calendar%type)
           endif
           secondsym = secondsym + (ndays * SECONDS_PER_DAY)
!           write(6,*) 'tcxti2 ',nmon,iyr,imo,ndays
        enddo
        if (mmon < 0) then
           secondsym = -secondsym
        endif
!        write(6,*) 'tcxti3 ',mmon,iyr,imo,secondsym
      elseif (PRESENT(D) .or. PRESENT(d_r8) .or. present(S) .or. present(S_i8)) then
        IF (timeinterval%MM /= 0) then
          CALL wrf_error_fatal("ESMF_TimeIntervalGet:  Need MM with D,d_r8,S,or S_i8")
        endif
        if (timeinterval%YR /= 0) then
          CALL wrf_error_fatal("ESMF_TimeIntervalGet:  Need YY or MM with D,d_r8,S,or S_i8")
        endif
      END IF

      seconds = seconds+secondsym

      IF ( PRESENT( D ) )THEN
        D = seconds / SECONDS_PER_DAY
        IF ( PRESENT(S) )   S    = mod( seconds, SECONDS_PER_DAY )
        IF ( PRESENT(S_i8)) S_i8 = mod( seconds, SECONDS_PER_DAY )
      ELSE
        IF ( PRESENT(S) )   S    = seconds
        IF ( PRESENT(S_i8)) S_i8 = seconds
      END IF

      IF ( PRESENT( d_r8 ) )THEN
        D_r8 = REAL( seconds, ESMF_KIND_R8 ) / &
               REAL( SECONDS_PER_DAY, ESMF_KIND_R8 )
      END IF

      ! If d_r8 present and sec present
      IF ( PRESENT( d_r8 ) )THEN
        IF ( PRESENT( S ) .or. present(s_i8) )THEN
          CALL wrf_error_fatal( &
            "ESMF_TimeIntervalGet:  Can not specify d_r8 and S S_i8 values" )
        END IF
      END IF

      ierr = ESMF_SUCCESS

      IF ( PRESENT( timeString ) ) THEN
        CALL ESMFold_TimeIntervalGetString( timeinterval, timeString, rc=ierr )
      ENDIF

      IF ( PRESENT(Sn) ) THEN
        Sn = timeinterval%basetime%Sn
      ENDIF
      IF ( PRESENT(Sd) ) THEN
        Sd = timeinterval%basetime%Sd
      ENDIF

      IF ( PRESENT(rc) ) rc = ierr

      end subroutine ESMF_TimeIntervalGet

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeIntervalSet - Initialize via user-specified unit set

! !INTERFACE:
!      subroutine ESMF_TimeIntervalSet(timeinterval, YY, YYl, MM, MOl, D, Dl, &
!                                      H, M, S, Sl, MS, US, NS, &
!                                      d_, d_r8, h_, m_, s_, ms_, us_, ns_, &
!                                      Sn, Sd, startTime, rc)
      subroutine ESMF_TimeIntervalSet(timeinterval, YY, MM, D, &
                                      H, M, S, S_i8, MS, &
                                      d_, d_r8, &
                                      Sn, Sd, startTime, rc)

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(out) :: timeinterval
      type(ESMF_Time), intent(in), optional :: StartTime
      integer, intent(in), optional :: YY
!      integer(ESMF_KIND_I8), intent(in), optional :: YYl
      integer, intent(in), optional :: MM
!      integer(ESMF_KIND_I8), intent(in), optional :: MOl
      integer, intent(in), optional :: D
!      integer(ESMF_KIND_I8), intent(in), optional :: Dl
      integer, intent(in), optional :: H
      integer, intent(in), optional :: M
      integer, intent(in), optional :: S
      integer(ESMF_KIND_I8), intent(in), optional :: S_i8
      integer, intent(in), optional :: MS
!      integer, intent(in), optional :: US
!      integer, intent(in), optional :: NS
      double precision, intent(in), optional :: d_
      double precision, intent(in), optional :: d_r8
!      double precision, intent(in), optional :: h_
!      double precision, intent(in), optional :: m_
!      double precision, intent(in), optional :: s_
!      double precision, intent(in), optional :: ms_
!      double precision, intent(in), optional :: us_
!      double precision, intent(in), optional :: ns_
      integer, intent(in), optional :: Sn
      integer, intent(in), optional :: Sd
      integer, intent(out), optional :: rc
      ! locals
      double precision :: din
      logical :: dinset

! !DESCRIPTION:
!     Set the value of the {\tt ESMF\_TimeInterval} in units specified by
!     the user via F90 optional arguments
!
!     Time manager represents and manipulates time internally with integers
!     to maintain precision.  Hence, user-specified floating point values are
!     converted internally to integers.
!
!     See {\tt ../include/ESMC\_BaseTime.h} and
!     {\tt ../include/ESMC\_TimeInterval.h} for complete description.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval]
!          The object instance to initialize
!     \item[{[YY]}]
!          Integer number of interval years (>= 32-bit)
!     \item[{[YYl]}]
!          Integer number of interval years (large, >= 64-bit)
!     \item[{[MM]}]
!          Integer number of interval months (>= 32-bit)
!     \item[{[MOl]}]
!          Integer number of interval months (large, >= 64-bit)
!     \item[{[D]}]
!          Integer number of interval days (>= 32-bit)
!     \item[{[Dl]}]
!          Integer number of interval days (large, >= 64-bit)
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
!     TMGn.n.n
!EOP

      IF ( PRESENT(rc) ) rc = ESMF_FAILURE

      timeinterval%startTime_set = .false.
      if (present(startTime)) then
         timeinterval%startTime = startTime
         timeinterval%startTime_set = .true.
      endif

      ! note that YR and MM are relative
      timeinterval%YR = 0
      IF ( PRESENT( YY ) ) THEN
        timeinterval%YR = YY
      ENDIF
      timeinterval%MM = 0
      IF ( PRESENT( MM ) ) THEN
        timeinterval%MM = MM
      ENDIF

      if (present(d_) .and. present(d_r8)) then
        CALL wrf_error_fatal( &
          "ESMF_TimeIntervalSet:  Cannot specify both d_r8 and d_")
      endif
      dinset = .false.
      if (present(d_))   then
        din = d_
        dinset = .true.
      endif
      if (present(d_r8)) then
        din = d_r8
        dinset = .true.
      endif
      IF ( dinset .AND. PRESENT( D ) ) THEN
        CALL wrf_error_fatal( &
          "ESMF_TimeIntervalSet:  Cannot specify both D and d_ or d_r8")
      ENDIF

      timeinterval%basetime%S = 0
      IF ( .NOT. dinset ) THEN
         IF ( PRESENT( D ) ) THEN
            timeinterval%basetime%S = timeinterval%basetime%S + &
                 ( SECONDS_PER_DAY * INT( D, ESMF_KIND_I8 ) )
         ENDIF
!$$$ push H,M,S,Sn,Sd,MS down into BaseTime constructor
         IF ( PRESENT( H ) ) THEN
            timeinterval%basetime%S = timeinterval%basetime%S + &
                 ( SECONDS_PER_HOUR * INT( H, ESMF_KIND_I8 ) )
         ENDIF
         IF ( PRESENT( M ) ) THEN
            timeinterval%basetime%S = timeinterval%basetime%S + &
                 ( SECONDS_PER_MINUTE * INT( M, ESMF_KIND_I8 ) )
         ENDIF
         IF ( PRESENT( S ) ) THEN
            timeinterval%basetime%S = timeinterval%basetime%S + &
                 INT( S, ESMF_KIND_I8 )
         ENDIF
         IF ( PRESENT( S_i8 ) ) THEN
            timeinterval%basetime%S = timeinterval%basetime%S + &
                 ( S_i8)
         ENDIF
      ELSE
         timeinterval%basetime%S = timeinterval%basetime%S + &
              INT( din * SECONDS_PER_DAY, ESMF_KIND_I8 )
      ENDIF
      IF ( PRESENT( Sn ) .AND. ( .NOT. PRESENT( Sd ) ) ) THEN
        CALL wrf_error_fatal( &
          "ESMF_TimeIntervalSet:  Must specify Sd if Sn is specified")
      ENDIF
      IF ( PRESENT( Sd ) .AND. PRESENT( MS ) ) THEN
        CALL wrf_error_fatal( &
          "ESMF_TimeIntervalSet:  Must not specify both Sd and MS")
      ENDIF
      timeinterval%basetime%Sn = 0
      timeinterval%basetime%Sd = 0
      IF ( PRESENT( MS ) ) THEN
        timeinterval%basetime%Sn = MS
        timeinterval%basetime%Sd = 1000_ESMF_KIND_I8
      ELSE IF ( PRESENT( Sd ) ) THEN
        timeinterval%basetime%Sd = Sd
        IF ( PRESENT( Sn ) ) THEN
          timeinterval%basetime%Sn = Sn
        ENDIF
      ENDIF
      CALL normalize_timeint( timeinterval )

      IF ( PRESENT(rc) ) rc = ESMF_SUCCESS

      end subroutine ESMF_TimeIntervalSet

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMFold_TimeIntervalGetString - Get time interval value in string format

! !INTERFACE:
      subroutine ESMFold_TimeIntervalGetString(timeinterval, TimeString, rc)

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval
      character*(*),  intent(out) :: TimeString
      integer, intent(out), optional :: rc
      ! locals
!      integer :: signnormtimeint
      LOGICAL :: negative
      INTEGER(ESMF_KIND_I8) :: iS, iSn, iSd, H, M, S, MM, D, YY
      character (len=1) :: signstr

! !DESCRIPTION:
!     Convert {\tt ESMF\_TimeInterval}'s value into string format
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval]
!          The object instance to convert
!     \item[TimeString]
!          The string to return
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.9
!EOP

! NOTE:  Sn, and Sd are not yet included in the returned string...
!PRINT *,'DEBUG ESMFold_TimeIntervalGetString():  YR,MM,S,Sn,Sd = ', &
!        timeinterval%YR, &
!        timeinterval%MM, &
!        timeinterval%basetime%S, &
!        timeinterval%basetime%Sn, &
!        timeinterval%basetime%Sd

      negative = ( signnormtimeint( timeInterval ) == -1 )
      IF ( negative ) THEN
        iS = -timeinterval%basetime%S
        iSn = -timeinterval%basetime%Sn
        signstr = '-'
      ELSE
        iS = timeinterval%basetime%S
        iSn = timeinterval%basetime%Sn
        signstr = ''
      ENDIF
      iSd = timeinterval%basetime%Sd

      D = iS / SECONDS_PER_DAY
      H = mod( iS, SECONDS_PER_DAY ) / SECONDS_PER_HOUR
      M = mod( iS, SECONDS_PER_HOUR) / SECONDS_PER_MINUTE
      S = mod( iS, SECONDS_PER_MINUTE )

!$$$here...  need to print Sn and Sd when they are used ???

      CALL timeintchecknormalized( timeinterval, 'ESMF_TimeIntervalGetString-arg1', &
                                   relative_interval=.true. )
      IF ( (timeinterval%MM == 0) .AND. (timeinterval%YR == 0) )THEN
         write(TimeString,FMT="(A,I10.10,'_',I3.3,':',I3.3,':',I3.3)") &
           TRIM(signstr), D, H, M, S
      ELSEif (timeinterval%YR == 0) then
         MM = timeinterval%MM
         write(TimeString,FMT="(I4.4, '_Months_',A,I10.10,'_',I3.3,':',I3.3,':',I3.3)") &
           MM, TRIM(signstr), D, H, M, S
      else
         YY = timeinterval%YR
         MM = timeinterval%MM
         write(TimeString,FMT="(I6.6,'_Years_',I4.4, '_Months_',A,I10.10,'_',I3.3,':',I3.3,':',I3.3)") &
           YY, MM, TRIM(signstr), D, H, M, S
      END IF

!write(0,*)'TimeIntervalGetString Sn ',timeinterval%basetime%Sn,' Sd ',timeinterval%basetime%Sd

      rc = ESMF_SUCCESS

      end subroutine ESMFold_TimeIntervalGetString

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalAbsValue - Get the absolute value of a time interval

! !INTERFACE:
      function ESMF_TimeIntervalAbsValue(timeinterval)

! !RETURN VALUE:
      type(ESMF_TimeInterval) :: ESMF_TimeIntervalAbsValue

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval
! !LOCAL:
      integer    :: rc

! !DESCRIPTION:
!     Return a {\tt ESMF\_TimeInterval}'s absolute value.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval]
!          The object instance to take the absolute value of.
!          Absolute value returned as value of function.
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.8
!EOP
      ESMF_TimeIntervalAbsValue = timeinterval
!$$$here...  move implementation into BaseTime
      ESMF_TimeIntervalAbsValue%basetime%S  = &
        abs(ESMF_TimeIntervalAbsValue%basetime%S)
      ESMF_TimeIntervalAbsValue%basetime%Sn = &
        abs(ESMF_TimeIntervalAbsValue%basetime%Sn )
      !
      ESMF_TimeIntervalAbsValue%MM = &
        abs(ESMF_TimeIntervalAbsValue%MM)
      ESMF_TimeIntervalAbsValue%YR = &
        abs(ESMF_TimeIntervalAbsValue%YR)

      end function ESMF_TimeIntervalAbsValue

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalNegAbsValue - Get the negative absolute value of a time interval

! !INTERFACE:
      function ESMF_TimeIntervalNegAbsValue(timeinterval)

! !RETURN VALUE:
      type(ESMF_TimeInterval) :: ESMF_TimeIntervalNegAbsValue

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval
! !LOCAL:
      integer    :: rc

! !DESCRIPTION:
!     Return a {\tt ESMF\_TimeInterval}'s negative absolute value.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval]
!          The object instance to take the negative absolute value of.
!          Negative absolute value returned as value of function.
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.8
!EOP
      ESMF_TimeIntervalNegAbsValue = timeinterval
!$$$here...  move implementation into BaseTime
      ESMF_TimeIntervalNegAbsValue%basetime%S  = &
        -abs(ESMF_TimeIntervalNegAbsValue%basetime%S)
      ESMF_TimeIntervalNegAbsValue%basetime%Sn = &
        -abs(ESMF_TimeIntervalNegAbsValue%basetime%Sn )
      !
      ESMF_TimeIntervalNegAbsValue%MM = &
        -abs(ESMF_TimeIntervalNegAbsValue%MM )
      ESMF_TimeIntervalNegAbsValue%YR = &
        -abs(ESMF_TimeIntervalNegAbsValue%YR )

      end function ESMF_TimeIntervalNegAbsValue

!------------------------------------------------------------------------------
!
! This section includes overloaded operators defined only for TimeInterval
! (not inherited from BaseTime)
! Note:  these functions do not have a return code, since F90 forbids more
! than 2 arguments for arithmetic overloaded operators
!
!------------------------------------------------------------------------------

! new WRF-specific function, Divide two time intervals and return the whole integer, without remainder
      function ESMF_TimeIntervalDIVQuot(timeinterval1, timeinterval2)

! !RETURN VALUE:
      INTEGER :: ESMF_TimeIntervalDIVQuot

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval1
      type(ESMF_TimeInterval), intent(in) :: timeinterval2

! !LOCAL
      INTEGER :: retval, isgn, rc
      type(ESMF_TimeInterval) :: zero, i1,i2

! !DESCRIPTION:
!     Returns timeinterval1 divided by timeinterval2 as a fraction quotient.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval1]
!          The dividend
!     \item[timeinterval2]
!          The divisor
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.5
!EOP

      CALL timeintchecknormalized( timeinterval1, 'ESMF_TimeIntervalDIVQuot arg1' )
      CALL timeintchecknormalized( timeinterval2, 'ESMF_TimeIntervalDIVQuot arg2' )

      call ESMF_TimeIntervalSet( zero, rc=rc )
      i1 = timeinterval1
      i2 = timeinterval2
      isgn = 1
      if ( i1 .LT. zero ) then
        i1 = WRFADDITION_TimeIntervalProdI(i1, -1)
        isgn = -isgn
      endif
      if ( i2 .LT. zero ) then
        i2 = WRFADDITION_TimeIntervalProdI(i2, -1)
        isgn = -isgn
      endif
! repeated subtraction
      retval = 0
      DO WHILE (  i1 .GE. i2 )
        i1 = i1 - i2
        retval = retval + 1
      ENDDO
      retval = retval * isgn

      ESMF_TimeIntervalDIVQuot = retval

      end function ESMF_TimeIntervalDIVQuot
! added by jhe
!------------------------------------------------------------------------------
!BOP
! !IROUTINE:   WRFADDITION_TimeIntervalProdI - Multiply a time interval by an
! integer

! !INTERFACE:
      function WRFADDITION_TimeIntervalProdI(timeinterval, multiplier)

! !RETURN VALUE:
      type(ESMF_TimeInterval) :: WRFADDITION_TimeIntervalProdI

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval
      integer, intent(in) :: multiplier
! !LOCAL:
      integer    :: rc

! !DESCRIPTION:
!     Multiply a {\tt ESMF\_TimeInterval} by an integer, return product
!     as a
!     {\tt ESMF\_TimeInterval}
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval]
!          The multiplicand
!     \item[mutliplier]
!          Integer multiplier
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.7, TMG7.2
!EOP
      CALL timeintchecknormalized( timeinterval, 'ESMF_TimeIntervalProdICarg1')

      CALL ESMF_TimeIntervalSet( WRFADDITION_TimeIntervalProdI, rc=rc )
!$$$move this into overloaded operator(*) in BaseTime
      WRFADDITION_TimeIntervalProdI%basetime%S  = &
        timeinterval%basetime%S * INT( multiplier, ESMF_KIND_I8 )
      WRFADDITION_TimeIntervalProdI%basetime%Sn = &
        timeinterval%basetime%Sn * INT( multiplier, ESMF_KIND_I8 )
      ! Don't multiply Sd
      WRFADDITION_TimeIntervalProdI%basetime%Sd = timeinterval%basetime%Sd
      CALL normalize_timeint( WRFADDITION_TimeIntervalProdI )

      end function WRFADDITION_TimeIntervalProdI

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalQuotI - Divide time interval by an integer, return time interval result

! !INTERFACE:
      function ESMF_TimeIntervalQuotI(timeinterval, divisor)

! !RETURN VALUE:
      type(ESMF_TimeInterval) :: ESMF_TimeIntervalQuotI

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval
      integer, intent(in) :: divisor

! !DESCRIPTION:
!     Divides a {\tt ESMF\_TimeInterval} by an integer divisor, returns
!     quotient as a {\tt ESMF\_TimeInterval}
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval]
!          The dividend
!     \item[divisor]
!          Integer divisor
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.6, TMG5.3, TMG7.2
!EOP

!PRINT *,'DEBUG ESMF_TimeIntervalQuotI() A:  S,Sn,Sd = ', &
!  timeinterval%basetime%S,timeinterval%basetime%Sn,timeinterval%basetime%Sd
!PRINT *,'DEBUG ESMF_TimeIntervalQuotI() A:  divisor = ', divisor

      CALL timeintchecknormalized( timeinterval, 'ESMF_TimeIntervalQuotI arg1' )

      IF ( divisor == 0 ) THEN
        CALL wrf_error_fatal( 'ESMF_TimeIntervalQuotI:  divide by zero' )
      ENDIF
      ESMF_TimeIntervalQuotI = timeinterval
!PRINT *,'DEBUG ESMF_TimeIntervalQuotI() B:  S,Sn,Sd = ', &
!  ESMF_TimeIntervalQuotI%basetime%S,ESMF_TimeIntervalQuotI%basetime%Sn,ESMF_TimeIntervalQuotI%basetime%Sd
      ESMF_TimeIntervalQuotI%basetime = timeinterval%basetime / divisor
!PRINT *,'DEBUG ESMF_TimeIntervalQuotI() C:  S,Sn,Sd = ', &
!  ESMF_TimeIntervalQuotI%basetime%S,ESMF_TimeIntervalQuotI%basetime%Sn,ESMF_TimeIntervalQuotI%basetime%Sd

      CALL normalize_timeint( ESMF_TimeIntervalQuotI )
!PRINT *,'DEBUG ESMF_TimeIntervalQuotI() D:  S,Sn,Sd = ', &
!  ESMF_TimeIntervalQuotI%basetime%S,ESMF_TimeIntervalQuotI%basetime%Sn,ESMF_TimeIntervalQuotI%basetime%Sd

      end function ESMF_TimeIntervalQuotI

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:   ESMF_TimeIntervalProdI - Multiply a time interval by an integer

! !INTERFACE:
      function ESMF_TimeIntervalProdI(timeinterval, multiplier)

! !RETURN VALUE:
      type(ESMF_TimeInterval) :: ESMF_TimeIntervalProdI

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval
      integer, intent(in) :: multiplier
! !LOCAL:
      integer    :: rc

! !DESCRIPTION:
!     Multiply a {\tt ESMF\_TimeInterval} by an integer, return product as a
!     {\tt ESMF\_TimeInterval}
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval]
!          The multiplicand
!     \item[mutliplier]
!          Integer multiplier
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.7, TMG7.2
!EOP
      CALL timeintchecknormalized( timeinterval, 'ESMF_TimeIntervalProdI arg1', &
                                   relative_interval=.true. )

      CALL ESMF_TimeIntervalSet( ESMF_TimeIntervalProdI, rc=rc )
!$$$move this into overloaded operator(*) in BaseTime
      ESMF_TimeIntervalProdI%basetime%S  = &
        timeinterval%basetime%S * INT( multiplier, ESMF_KIND_I8 )
      ESMF_TimeIntervalProdI%basetime%Sn = &
        timeinterval%basetime%Sn * INT( multiplier, ESMF_KIND_I8 )
      ! Don't multiply Sd
      ESMF_TimeIntervalProdI%basetime%Sd = timeinterval%basetime%Sd
      ESMF_TimeIntervalProdI%MM = timeinterval%MM * multiplier
      ESMF_TimeIntervalProdI%YR = timeinterval%YR * multiplier
      CALL normalize_timeint( ESMF_TimeIntervalProdI )

      end function ESMF_TimeIntervalProdI

!------------------------------------------------------------------------------
!
! This section includes the inherited ESMF_BaseTime class overloaded operators
!
!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalSum - Add two time intervals together

! !INTERFACE:
      function ESMF_TimeIntervalSum(timeinterval1, timeinterval2)

! !RETURN VALUE:
      type(ESMF_TimeInterval) :: ESMF_TimeIntervalSum

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval1
      type(ESMF_TimeInterval), intent(in) :: timeinterval2
! !LOCAL:
      integer                             :: rc
! !DESCRIPTION:
!     Add two {\tt ESMF\_TimeIntervals}, return sum as a
!     {\tt ESMF\_TimeInterval}.  Maps overloaded (+) operator interface
!     function to {\tt ESMF\_BaseTime} base class.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval1]
!          The augend
!     \item[timeinterval2]
!          The addend
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.4, TMG2.4.4, TMG2.4.5, TMG2.4.6, TMG5.1, TMG5.2,
!                 TMG7.2
!EOP
      CALL timeintchecknormalized( timeinterval1, 'ESMF_TimeIntervalSum arg1', &
                                   relative_interval=.true. )
      CALL timeintchecknormalized( timeinterval2, 'ESMF_TimeIntervalSum arg2', &
                                   relative_interval=.true. )

      ESMF_TimeIntervalSum = timeinterval1
      ESMF_TimeIntervalSum%basetime = ESMF_TimeIntervalSum%basetime + &
                                      timeinterval2%basetime
      ESMF_TimeIntervalSum%MM = ESMF_TimeIntervalSum%MM + &
                                      timeinterval2%MM
      ESMF_TimeIntervalSum%YR = ESMF_TimeIntervalSum%YR + &
                                      timeinterval2%YR

      CALL normalize_timeint( ESMF_TimeIntervalSum )

      end function ESMF_TimeIntervalSum

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalDiff - Subtract one time interval from another

! !INTERFACE:
      function ESMF_TimeIntervalDiff(timeinterval1, timeinterval2)

! !RETURN VALUE:
      type(ESMF_TimeInterval) :: ESMF_TimeIntervalDiff

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval1
      type(ESMF_TimeInterval), intent(in) :: timeinterval2
! !LOCAL:
      integer                             :: rc
! !DESCRIPTION:
!     Subtract timeinterval2 from timeinterval1, return remainder as a
!     {\tt ESMF\_TimeInterval}.
!     Map overloaded (-) operator interface function to {\tt ESMF\_BaseTime}
!     base class.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval1]
!          The minuend
!     \item[timeinterval2]
!          The subtrahend
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.4, TMG2.4.4, TMG2.4.5, TMG2.4.6, TMG5.1, TMG5.2, TMG7.2
!EOP
      CALL timeintchecknormalized( timeinterval1, 'ESMF_TimeIntervalDiff arg1', &
                                   relative_interval=.true. )
      CALL timeintchecknormalized( timeinterval2, 'ESMF_TimeIntervalDiff arg2', &
                                   relative_interval=.true. )

      ESMF_TimeIntervalDiff = timeinterval1
      ESMF_TimeIntervalDiff%basetime = ESMF_TimeIntervalDiff%basetime - &
                                       timeinterval2%basetime
      ESMF_TimeIntervalDiff%MM       = ESMF_TimeIntervalDiff%MM       - &
                                       timeinterval2%MM
      ESMF_TimeIntervalDiff%YR       = ESMF_TimeIntervalDiff%YR       - &
                                       timeinterval2%YR
      CALL normalize_timeint( ESMF_TimeIntervalDiff )

      end function ESMF_TimeIntervalDiff

!------------------------------------------------------------------------------
!BOP
! !IROUTINE: ESMF_TimeIntervalEQ - Compare two time intervals for equality

! !INTERFACE:
      function ESMF_TimeIntervalEQ(timeinterval1, timeinterval2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeIntervalEQ

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval1
      type(ESMF_TimeInterval), intent(in) :: timeinterval2

!DESCRIPTION:
!     Return true if both given time intervals are equal, false otherwise.
!     Maps overloaded (==) operator interface function to {\tt ESMF\_BaseTime}
!     base class.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval1]
!          First time interval to compare
!     \item[timeinterval2]
!          Second time interval to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      INTEGER :: res

      CALL timeintcmp(timeinterval1,timeinterval2,res)
      ESMF_TimeIntervalEQ = (res .EQ. 0)

      end function ESMF_TimeIntervalEQ

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalNE - Compare two time intervals for inequality

! !INTERFACE:
      function ESMF_TimeIntervalNE(timeinterval1, timeinterval2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeIntervalNE

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval1
      type(ESMF_TimeInterval), intent(in) :: timeinterval2

! !DESCRIPTION:
!     Return true if both given time intervals are not equal, false otherwise.
!     Maps overloaded (/=) operator interface function to {\tt ESMF\_BaseTime}
!     base class.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval1]
!          First time interval to compare
!     \item[timeinterval2]
!          Second time interval to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      INTEGER :: res

      CALL timeintcmp(timeinterval1,timeinterval2,res)
      ESMF_TimeIntervalNE = (res .NE. 0)

      end function ESMF_TimeIntervalNE

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalLT - Time interval 1 less than time interval 2 ?

! !INTERFACE:
      function ESMF_TimeIntervalLT(timeinterval1, timeinterval2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeIntervalLT

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval1
      type(ESMF_TimeInterval), intent(in) :: timeinterval2

! !DESCRIPTION:
!     Return true if first time interval is less than second time interval,
!     false otherwise. Maps overloaded (<) operator interface function to
!     {\tt ESMF\_BaseTime} base class.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval1]
!          First time interval to compare
!     \item[timeinterval2]
!          Second time interval to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      INTEGER :: res

      CALL timeintcmp(timeinterval1,timeinterval2,res)
      ESMF_TimeIntervalLT = (res .LT. 0)

      end function ESMF_TimeIntervalLT

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalGT - Time interval 1 greater than time interval 2?

! !INTERFACE:
      function ESMF_TimeIntervalGT(timeinterval1, timeinterval2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeIntervalGT

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval1
      type(ESMF_TimeInterval), intent(in) :: timeinterval2

! !DESCRIPTION:
!     Return true if first time interval is greater than second time interval,
!     false otherwise.  Maps overloaded (>) operator interface function to
!     {\tt ESMF\_BaseTime} base class.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval1]
!          First time interval to compare
!     \item[timeinterval2]
!          Second time interval to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      INTEGER :: res

      CALL timeintcmp(timeinterval1,timeinterval2,res)
      ESMF_TimeIntervalGT = (res .GT. 0)

      end function ESMF_TimeIntervalGT

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalLE - Time interval 1 less than or equal to time interval 2 ?

! !INTERFACE:
      function ESMF_TimeIntervalLE(timeinterval1, timeinterval2)

! !RETURN VALUE:
      logical :: ESMF_TimeIntervalLE

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval1
      type(ESMF_TimeInterval), intent(in) :: timeinterval2

! !DESCRIPTION:
!     Return true if first time interval is less than or equal to second time
!     interval, false otherwise.
!     Maps overloaded (<=) operator interface function to {\tt ESMF\_BaseTime}
!     base class.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval1]
!          First time interval to compare
!     \item[timeinterval2]
!          Second time interval to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      INTEGER :: res

      CALL timeintcmp(timeinterval1,timeinterval2,res)
      ESMF_TimeIntervalLE = (res .LE. 0)

      end function ESMF_TimeIntervalLE

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalGE - Time interval 1 greater than or equal to time interval 2 ?

! !INTERFACE:
      function ESMF_TimeIntervalGE(timeinterval1, timeinterval2)
!
! !RETURN VALUE:
      logical :: ESMF_TimeIntervalGE

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval1
      type(ESMF_TimeInterval), intent(in) :: timeinterval2

! !DESCRIPTION:
!     Return true if first time interval is greater than or equal to second
!     time interval, false otherwise. Maps overloaded (>=) operator interface
!     function to {\tt ESMF\_BaseTime} base class.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval1]
!          First time interval to compare
!     \item[timeinterval2]
!          Second time interval to compare
!     \end{description}
!
! !REQUIREMENTS:
!     TMG1.5.3, TMG2.4.3, TMG7.2
!EOP

      INTEGER :: res

      CALL timeintcmp(timeinterval1,timeinterval2,res)
      ESMF_TimeIntervalGE = (res .GE. 0)

      end function ESMF_TimeIntervalGE

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalIsPositive - Time interval greater than zero?

! !INTERFACE:
      function ESMF_TimeIntervalIsPositive(timeinterval)
!
! !RETURN VALUE:
      logical :: ESMF_TimeIntervalIsPositive

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval

! !LOCALS:
      type(ESMF_TimeInterval) :: zerotimeint
      integer :: rcint

! !DESCRIPTION:
!     Return true if time interval is greater than zero,
!     false otherwise.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval]
!          Time interval to compare
!     \end{description}
!EOP
      CALL timeintchecknormalized( timeinterval, &
                                   'ESMF_TimeIntervalIsPositive arg' )

      CALL ESMF_TimeIntervalSet ( zerotimeint, rc=rcint )
      IF ( rcint /= ESMF_SUCCESS ) THEN
        CALL wrf_error_fatal( &
          'ESMF_TimeIntervalIsPositive:  ESMF_TimeIntervalSet failed' )
      ENDIF
! hack for bug in PGI 5.1-x
!      ESMF_TimeIntervalIsPositive = timeinterval > zerotimeint
      ESMF_TimeIntervalIsPositive = ESMF_TimeIntervalGT( timeinterval, &
                                                         zerotimeint )
      end function ESMF_TimeIntervalIsPositive

!------------------------------------------------------------------------------
!BOP
! !IROUTINE:  ESMF_TimeIntervalPrint - Print out a time interval's properties

! !INTERFACE:
      subroutine ESMF_TimeIntervalPrint(timeinterval, opts, rc)

! !ARGUMENTS:
      type(ESMF_TimeInterval), intent(in) :: timeinterval
      character (len=*), intent(in), optional :: opts
      integer, intent(out), optional :: rc

! !DESCRIPTION:
!     To support testing/debugging, print out an {\tt ESMF\_TimeInterval}'s
!     properties.
!
!     The arguments are:
!     \begin{description}
!     \item[timeinterval]
!          Time interval to print out
!     \item[{[opts]}]
!          Print options
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
! !REQUIREMENTS:
!     TMGn.n.n
!EOP
      INTEGER :: ierr

      ierr = ESMF_SUCCESS
      call print_a_timeinterval( timeinterval )
      IF ( PRESENT(rc) ) rc = ierr

      end subroutine ESMF_TimeIntervalPrint

!------------------------------------------------------------------------------

! Exits with error message if timeInt is not normalized.
SUBROUTINE timeintchecknormalized( timeInt, msgstr, relative_interval )
  IMPLICIT NONE
  TYPE(ESMF_TimeInterval), INTENT(IN) :: timeInt
  CHARACTER(LEN=*), INTENT(IN) :: msgstr
  LOGICAL, INTENT(IN), optional :: relative_interval   ! If relative intervals are ok or not
  ! locals
  CHARACTER(LEN=256) :: outstr
  LOGICAL :: non_relative

  IF ( .NOT. PRESENT( relative_interval ) )THEN
     non_relative = .true.
  ELSE
     IF ( relative_interval )THEN
        non_relative = .false.
     ELSE
        non_relative = .true.
     END IF
  END IF
  IF ( non_relative )THEN
     IF ( ( timeInt%YR /= 0 ) .OR. &
          ( timeInt%MM /= 0 ) ) THEN
       outstr = 'un-normalized TimeInterval not allowed:  '//TRIM(msgstr)
       CALL wrf_error_fatal( outstr )
     ENDIF
  ELSE
     IF ( ( timeInt%YR /= 0 ) .OR. &
          ( timeInt%MM < -MONTHS_PER_YEAR) .OR. ( timeInt%MM > MONTHS_PER_YEAR ) ) THEN
! tcraig, don't require normalize TimeInterval for relative diffs
!       outstr = 'un-normalized TimeInterval not allowed:  '//TRIM(msgstr)
!       CALL wrf_error_fatal( outstr )
     ENDIF
  END IF
END SUBROUTINE timeintchecknormalized

!==============================================================================
SUBROUTINE print_a_timeinterval( time )
   IMPLICIT NONE
   type(ESMF_TimeInterval) time
   character*128 :: s
   integer rc
   CALL ESMFold_TimeIntervalGetString( time, s, rc )
   write(6,*)'Print a time interval|',time%yr, time%mm, time%basetime%s, time%starttime_set, time%starttime%calendar%type%caltype
   write(6,*)'Print a time interval|',TRIM(s),'|'
   return
END SUBROUTINE print_a_timeinterval

!==============================================================================

SUBROUTINE timeintcmp(timeint1in, timeint2in, retval )
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: retval
!
! !ARGUMENTS:
  TYPE(ESMF_TimeInterval), INTENT(IN) :: timeint1in
  TYPE(ESMF_TimeInterval), INTENT(IN) :: timeint2in

  TYPE(ESMF_TimeInterval) :: timeint1
  TYPE(ESMF_TimeInterval) :: timeint2

  timeint1 = timeint1in
  timeint2 = timeint2in
  call normalize_timeint(timeint1)
  call normalize_timeint(timeint2)

  IF ( (timeint1%MM /= timeint2%MM) .and. (timeint1%YR /= timeint2%YR) )THEN
    CALL wrf_error_fatal( &
      'timeintcmp:  Can not compare two intervals with different months and years' )
  END IF
  if (timeint1%YR .gt. timeint2%YR) then
     retval = 1
  elseif (timeint1%YR .lt. timeint2%YR) then
     retval = -1
  else
     if (timeint1%MM .gt. timeint2%MM) then
        retval = 1
     elseif (timeint1%MM .lt. timeint2%MM) then
        retval = 1
     else
        CALL seccmp( timeint1%basetime%S, timeint1%basetime%Sn, &
               timeint1%basetime%Sd,                      &
               timeint2%basetime%S, timeint2%basetime%Sn, &
               timeint2%basetime%Sd, retval )
     endif
  endif

END SUBROUTINE timeintcmp

!==============================================================================

SUBROUTINE normalize_timeint( timeInt )
  IMPLICIT NONE
  TYPE(ESMF_TimeInterval), INTENT(INOUT) :: timeInt
  INTEGER :: mpyi4

  ! normalize basetime
  ! this will force abs(Sn) < Sd and ensure that signs of S and Sn match
  ! YR and MM are ignored

  CALL normalize_basetime( timeInt%basetime )

  ! Rollover months to years

  mpyi4 = MONTHS_PER_YEAR
  IF      ( abs(timeInt%MM) .GE. MONTHS_PER_YEAR ) THEN
    timeInt%YR = timeInt%YR + timeInt%MM/MONTHS_PER_YEAR
    timeInt%MM = mod(timeInt%MM,mpyi4)
  ENDIF

  ! make sure yr and mm have same sign

  IF (timeInt%YR * timeInt%MM < 0) then
     if (timeInt%YR > 0) then
        timeInt%MM = timeInt%MM + MONTHS_PER_YEAR
        timeInt%YR = timeInt%YR - 1
     endif
     if (timeInt%YR < 0) then
        timeInt%MM = timeInt%MM - MONTHS_PER_YEAR
        timeInt%YR = timeInt%YR + 1
     endif
  endif

END SUBROUTINE normalize_timeint

!==============================================================================

integer FUNCTION signnormtimeint ( timeInt )
  ! Compute the sign of a time interval.
  ! YR and MM fields are *IGNORED*.
  ! returns 1, 0, or -1 or exits if timeInt fields have inconsistent signs.
  IMPLICIT NONE
  TYPE(ESMF_TimeInterval), INTENT(IN) :: timeInt
  LOGICAL :: positive, negative

  positive = .FALSE.
  negative = .FALSE.
  signnormtimeint = 0
  ! Note that Sd is required to be non-negative.  This is enforced in
  ! normalize_timeint().
  ! Note that Sn is required to be zero when Sd is zero.  This is enforced
  ! in normalize_timeint().
  IF ( ( timeInt%basetime%S > 0 ) .OR. &
       ( timeInt%basetime%Sn > 0 ) ) THEN
    positive = .TRUE.
  ENDIF
  IF ( ( timeInt%basetime%S < 0 ) .OR. &
       ( timeInt%basetime%Sn < 0 ) ) THEN
    negative = .TRUE.
  ENDIF
  IF ( positive .AND. negative ) THEN
    CALL wrf_error_fatal( &
      'signnormtimeint:  signs of fields cannot be mixed' )
  ELSE IF ( positive ) THEN
    signnormtimeint = 1
  ELSE IF ( negative ) THEN
    signnormtimeint = -1
  ENDIF
END FUNCTION signnormtimeint
!==============================================================================

end module ESMF_TimeIntervalMod

