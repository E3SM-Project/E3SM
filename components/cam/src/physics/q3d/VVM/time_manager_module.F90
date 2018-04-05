MODULE time_manager_module
! Defines the year and fractional Julian day of the current time step.

USE shr_kind_mod, only: r8 => shr_kind_r8

USE constld,  only: dt
USE timeinfo, only: rjday0,rjday,utc_time,iyr,imonth,iday

IMPLICIT NONE
PRIVATE

INTEGER, DIMENSION(12) :: N

PUBLIC :: time_manager

CONTAINS

! Local Subroutines:
!------------------------------------------------------------------------
! SUBROUTINE time_manager
! SUBROUTINE JULIAN : Calculates the julian day (1-365) given month, day, and year.
! SUBROUTINE MODAY  : Calculates the month and day given the julian day and year.
!------------------------------------------------------------------------

!========================================================================
   SUBROUTINE TIME_MANAGER(ITT)
!========================================================================
!  Defines the year and fractional Julian day of the current time step.
!
!  INPUT : DT, ITT
!  OUTPUT: RJDAY0,RJDAY,IYR
!-----------------------------------------------------------------------
      INTEGER, INTENT(IN) :: itt   ! time step count

      INTEGER :: iyear0,imonth0,iday0,jday0, ntspr
      REAL (KIND=r8) :: remder, hour0, secday
      DATA SECDAY/86400.0_r8/

! Specify year, month, day, and hour of initial time step

      iyear0  = 1
      imonth0 = 1
      iday0   = 1
      hour0   = 0.0_r8

!  Get Julian day
      CALL JULIAN(iyear0,imonth0,iday0,jday0)

!  Specify model start time in fractional Julian day
      rjday0 = REAL(jday0) + hour0/24.0_r8

!  Julian day of current time step
      rjday = rjday0 + (ITT-1) * DT / SECDAY
      iyr = iyear0

      IF (rjday .GE. 366.0_r8) THEN
        rjday = rjday - 365.0_r8
        iyr = iyr + 1
      ENDIF

!  Month, day, and UTC time of current time step
      CALL MODAY(iyr,imonth,iday,INT(rjday))
      utc_time = (rjday - INT(rjday)) * 24.0_r8

   END SUBROUTINE time_manager

!========================================================================      
   SUBROUTINE JULIAN(IYR,IMO,IDAY,rjday)
!========================================================================   
!  CALCULATES THE JULIAN DAY (1-365) GIVEN MONTH, DAY, AND YEAR.
!  IYR = YEAR (YYYY).
!  IMO = MONTH.
!  IDY = DAY.
!  RJDAY = JULIAN DAY.

      INTEGER, INTENT(IN) :: IYR
      INTEGER, INTENT(INOUT) :: IMO,IDAY
      INTEGER, INTENT(OUT) :: RJDAY
      INTEGER ::  M,I

      N = (/31,0,31,30,31,30,31,31,30,31,30,31/)
      N(2)=28

! Include this line if using leap years
!      IF(MOD(IYR,4).EQ.0) N(2)=29
      IF(IMO.LE.1) THEN
        RJDAY=IDAY
        RETURN
      ENDIF
      M=IMO-1
      RJDAY=0
      DO 20 I=1,M
   20 RJDAY=RJDAY+N(I)
      RJDAY=RJDAY+IDAY

   END SUBROUTINE julian

!========================================================================
   SUBROUTINE MODAY(IYR,IMO,IDAY,rjday)
!========================================================================
!  CALCULATES THE MONTH AND DAY GIVEN THE JULIAN DAY AND YEAR.

      INTEGER, INTENT(IN) :: IYR
      INTEGER, INTENT(INOUT) :: IMO,IDAY
      INTEGER, INTENT(IN) :: RJDAY

      INTEGER :: I,KDAY

      N = (/31,0,31,30,31,30,31,31,30,31,30,31/)
      N(2)=28

! Include this line if using leap years
!      IF(MOD(IYR,4).EQ.0) N(2)=29
      KDAY=RJDAY
      DO 30 I=1,12
        IF(KDAY.LE.N(I)) EXIT
        KDAY=KDAY-N(I)
   30 CONTINUE
      IMO=I
      IDAY=KDAY

   END SUBROUTINE moday

END MODULE time_manager_module
