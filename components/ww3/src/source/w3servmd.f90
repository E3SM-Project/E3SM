!/ ------------------------------------------------------------------- /
      MODULE W3SERVMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    For update log see individual subroutines.
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     In this module all WAVEWATCH specific service routines have
!     been gathered.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NDSTRC    Int.  Private  Data set number for output of STRACE
!                               (set in ITRACE).
!      NTRACE    Int.  Private  Maximum number of trace prints in
!                               strace (set in ITRACE).
!
!      PRFTB     Int.  Private  Base time for profiling.
!      FLPROF    Log.  Private  Flag for profiling initialization.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      ITRACE    Subr. Public   (Re-) Initialization for STRACE.
!      STRACE    Subr. Public   Enable subroutine tracing, usually
!                               activated with the !/S switch.
!      NEXTLN    Subr. Public   Get to next line in input command file.
!      W3S2XY    Subr. Public   Grid conversion routine.
!      EJ5P      R.F.  Public   Five parameter JONSWAP spectrum.
!      WWDATE    Subr. Public   Get system date.
!      WWTIME    Subr. Public   Get system time.
!      EXTCDE    Subr. Public   Abort program with exit code.
!      PRINIT    Subr. Public   Initialize profiling.
!      PRTIME    Subr. Public   Get profiling time.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!     None.
!
!  5. Remarks :
!
!  6. Switches
!
!       !/S    Enable subroutine tracing using STRACE in this module.
!
!       !/F90  FORTRAN 90 specific switches.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!
      INTEGER, PRIVATE        :: NDSTRC = 6, NTRACE = 0, PRFTB
      LOGICAL, PRIVATE        :: FLPROF = .FALSE.
!
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE ITRACE (NDS, NMAX)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         23-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    23-Nov-1999 : First version of routine.           ( version 2.00 )
!/
!  1. Purpose :
!
!     (Re-) initialization for module version of STRACE.
!
!  3. Parameter list
!     ----------------------------------------------------------------
!       NDS     Int.   I   Data set number ofr trace file.
!       NMAX    Int.   I   Maximum number of traces per routine.
!     ----------------------------------------------------------------
!
!     Private to module :
!     ----------------------------------------------------------------
!       NDSTRC  Int.  Output unit number for trace.     ( from NDS  )
!       NTRACE  Int.  Maximum number of trace prints.   ( from NMAX )
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None.
!
!  5. Called by :
!
!     Any program, multiple calls allowed.
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NDS, NMAX
!/
!/ ------------------------------------------------------------------- /
!/
      NTRACE = MAX ( 0 , NMAX )
      NDSTRC = NDS
!
      RETURN
!/
!/ End of ITRACE ----------------------------------------------------- /
!/
      END SUBROUTINE ITRACE
!/ ------------------------------------------------------------------- /
      SUBROUTINE STRACE (IENT, SNAME)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         25-Jan-2000 |
!/                  +-----------------------------------+
!/                                   Original version by N. Booij, DUT
!/
!/    30-Mar-1993 : Final FORTRAN 77                    ( version 1.18 )
!/    23-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    25-Jan-2000 : Force flushing of uniit.            ( version 2.00 )
!/
!  1. Purpose :
!
!     Keep track of entered subroutines.
!
!  3. Parameter list
!     ----------------------------------------------------------------
!       IENT    Int.  I/O  Number of times that STRACE has been
!                          called by the routine.
!       SNAME   Char.  I   Name of the subroutine (max. 6 characters)
!     ----------------------------------------------------------------
!
!     Private to module :
!     ----------------------------------------------------------------
!       NDSTRC  Int.  Output unit number for trace.
!       NTRACE  Int.  Maximum number of trace prints.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None.
!
!  5. Called by :
!
!     Any program, after private variables have been set by NTRACE.
!
!  9. Switches :
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(INOUT)  :: IENT
      CHARACTER, INTENT(IN)   :: SNAME*(*)
!/
!/ ------------------------------------------------------------------- /
!/
      IF (NTRACE.EQ.0 .OR. IENT.GE.NTRACE) RETURN
!
      IENT = IENT + 1
      IF (IENT.EQ.1) THEN
          WRITE (NDSTRC,10) SNAME
        ELSE
          WRITE (NDSTRC,11) SNAME, IENT
        END IF
!
      RETURN
!
! Formats
!
  10  FORMAT (' ---> TRACE SUBR : ',A6)
  11  FORMAT (' ---> TRACE SUBR : ',A6,'  ENTRY: ',I6)
!/
!/ End of STRACE ----------------------------------------------------- /
!/
      END SUBROUTINE STRACE
!/ ------------------------------------------------------------------- /
      SUBROUTINE NEXTLN ( CHCKC , NDSI , NDSE )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         18-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    15-Jan-1999 : Final FORTRAN 77                    ( version 1.18 )
!/    18-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Sets file pointer to next active line of input file, by skipping
!     lines starting with the character CHCKC.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       CHCKC   C*1   I  Check character for defining comment line.
!       NDSI    Int.  I  Input dataset number.
!       NDSE    Int.  I  Error output dataset number.
!                        (No output if NDSE < 0).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE ( !/S switch )
!
!  5. Called by :
!
!       Any routine.
!
!  6. Error messages :
!
!     - On EOF or error in input file.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NDSI, NDSE
      CHARACTER, INTENT(IN)   :: CHCKC*1
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IERR
      CHARACTER               :: TEST*1
!/
!/ ------------------------------------------------------------------- /
!/
!
  100 CONTINUE
      READ (NDSI,900,END=800,ERR=801,IOSTAT=IERR) TEST
      IF (TEST.EQ.CHCKC) THEN
          GOTO 100
        ELSE
          BACKSPACE (NDSI,ERR=802,IOSTAT=IERR)
        ENDIF
      RETURN
!
  800 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,910)
      CALL EXTCDE ( 1 )
!
  801 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,911) IERR
      CALL EXTCDE ( 2 )
!
  802 CONTINUE
      IF ( NDSE .GE. 0 ) WRITE (NDSE,911) IERR
      CALL EXTCDE ( 3 )
!
! Formats
!
  900 FORMAT (A)
  910 FORMAT (/' *** WAVEWATCH III ERROR IN NEXTLN : '/         &
               '     PREMATURE END OF INPUT FILE'/)
  911 FORMAT (/' *** WAVEWATCH III ERROR IN NEXTLN : '/         &
               '     ERROR IN READING FROM FILE'/               &
               '     IOSTAT =',I5/)
  912 FORMAT (/' *** WAVEWATCH III ERROR IN NEXTLN : '/         &
               '     ERROR ON BACKSPACE'/                       &
               '     IOSTAT =',I5/)
!/
!/ End of NEXTLN ----------------------------------------------------- /
!/
      END SUBROUTINE NEXTLN
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3S2XY ( NSEA, MSEA, MX, MY, S, MAPSF, XY )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III            NOAA/NMC |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         23-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    11-Dec-1996 : Final FORTRAN 77                    ( version 1.18 )
!/    23-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Convert a data array on the storage grid to a data array on the
!     full spatial grid. Land and ice points in the full grid are
!     not touched. Output array of conventional type XY(IX,IY).
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NSEA    Int.   I    Number of sea points.
!       MSEA, MX, MY
!               Int.   I    Array dimensions.
!       S       R.A.   I    Data on storage grid.
!       MAPSF   I.A.   I    Storage map for IX and IY, resp.
!       XY      R.A.   O    Data on XY grid.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None.
!
!  5. Called by :
!
!     Any WAVEWATCH III routine.
!
!  9. Switches :
!
!     None.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: MSEA, NSEA, MX, MY, MAPSF(MSEA,2)
      REAL, INTENT(IN)        :: S(MSEA)
      REAL, INTENT(OUT)       :: XY(MX,MY)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ISEA, IX, IY
!/
!/ ------------------------------------------------------------------- /
!/
      DO 100, ISEA=1, NSEA
        IX     = MAPSF(ISEA,1)
        IY     = MAPSF(ISEA,2)
        XY(IX,IY) = S(ISEA)
  100   CONTINUE
!/
!/ End of W3S2XY ----------------------------------------------------- /
!/
      END SUBROUTINE W3S2XY
!/ ------------------------------------------------------------------- /
      REAL FUNCTION EJ5P ( F, ALFA, FP, YLN, SIGA, SIGB )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         23-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    23-AMy-1985 : Original by G. Ph. van Vledder.
!/    23-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Computation of spectral density using a 5-parameter
!     JONSWAP-spectrum
!
!  2. Method
!
!     EJ5P(F) = A.EXP(B + LN(Y).EXP(C))
!
!     where: A = ALFA * 0.06175 * F**(-5)
!            B = -1.25*(FP/F)**4
!            C = -0.5 * ((F - FP)/(SIG * FP))**2
!     and
!            GRAV**2/(2.PI)**4 = 0.06175
!
!  3. Parameters :
!
!     Parameter list
!
!     ----------------------------------------------------------------
!       F       Real   I    Frequency in Hz
!       ALFA    Real   I    Energy scaling factor
!       FP      Real   I    Peak frequency in Hz
!       YLN     Real   I    Peak overshoot factor, given by LN-value
!       SIGA    Real   I    Spectral width, for F < FP
!       SIGB    Real   I    Spectral width, FOR F > FP
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     None.
!
!  5. Called by :
!
!     Any.
!
!  6. Error messages :
!
!  7. Remarks :
!
!     EXPMIN is a machine dependant constant such that
!     EXP(EXPMIN) can be successfully evaluated without
!     underflow by the compiler supllied EXP routine.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     None.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: F, ALFA, FP, YLN, SIGA, SIGB
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      REAL                    :: SIG, A, B, C
      REAL, SAVE              :: EPS=1.E-4, EXPMIN=-180.
!/
!/ ------------------------------------------------------------------- /
!/
      IF(F.LT.EPS) THEN
        EJ5P = 0.0
        RETURN
      END IF
!
      A = ALFA * 0.06175 / F**5
      B = -1.25 * (FP/F)**4
      B = MAX(B,EXPMIN)
!
      IF (YLN.LT.EPS) THEN
        EJ5P = A * EXP(B)
      ELSE
        IF( F.LE.FP) THEN
          SIG = SIGA
        ELSE
          SIG = SIGB
        END IF
        C = -0.5 * ((F - FP)/(SIG * FP))**2
        C = MAX(C,EXPMIN)
        EJ5P = A * EXP(B + EXP(C) * YLN)
      END IF
!
      RETURN
!/
!/ End of NEXTLN ----------------------------------------------------- /
!/
      END FUNCTION
!/ ------------------------------------------------------------------- /
      SUBROUTINE WWDATE (STRNG)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         08-May-2002 |
!/                  +-----------------------------------+
!/
!/    23-Dec-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    23-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    18-Sep-2000 : PGI switch added                    ( version 2.04 )
!/    13-Mar-2001 : LF95 switch added                   ( version 2.09 )
!/    08-May-2002 : Replace obsolete switches with F90  ( version 2.21 )
!/
!  1. Purpose :
!
!     Get date from machine dependent routine.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       STRNG   C*10   O   String with date in format YYYY/MM/DD
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     Machine dependent.
!
!  5. Called by :
!
!     Any routine.
!
!  9. Switches :
!
!     !/DUM  Dummy.
!     !/F90  FORTRAN 90 standard.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      CHARACTER, INTENT(OUT)  :: STRNG*10
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      CHARACTER*8             :: DATE
      CHARACTER*10            :: TIME
      CHARACTER*5             :: ZONE
      INTEGER                 :: VALUES(8)
!/
!/ ------------------------------------------------------------------- /
!/
! This is supposed to be standard F90
!
      STRNG = '----/--/--'
      CALL DATE_AND_TIME ( DATE, TIME, ZONE, VALUES )
      STRNG(1:4) = DATE(1:4)
      STRNG(6:7) = DATE(5:6)
      STRNG(9:10) = DATE(7:8)
!
! Dummy alternative
!
      RETURN
!/
!/ End of WWDATE ----------------------------------------------------- /
!/
      END SUBROUTINE WWDATE
!/ ------------------------------------------------------------------- /
      SUBROUTINE WWTIME (STRNG)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         08-May-2002 |
!/                  +-----------------------------------+
!/
!/    23-Dec-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    23-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    18-Sep-2000 : PGI switch added                    ( version 2.04 )
!/    13-Mar-2001 : LF95 switch added                   ( version 2.09 )
!/    08-May-2002 : Replace obsolete switches with F90  ( version 2.21 )
!/
!  1. Purpose :
!
!     Get time from machine dependent routine.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       STRNG   C*8    O   String with time in format hh:mm:ss
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     Machine dependent.
!
!  5. Called by :
!
!     Any routine.
!
!  9. Switches :
!
!     !/DUM  Dummy.
!     !/F90  FORTRAN 90 standard.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      CHARACTER, INTENT(OUT)  :: STRNG*8
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      CHARACTER*8             :: DATE
      CHARACTER*10            :: TIME
      CHARACTER*5             :: ZONE
      INTEGER                 :: VALUES(8)
!/
!/ ------------------------------------------------------------------- /
!/
! This is supposed to be standard F90
!
      STRNG = '--:--:--'
      CALL DATE_AND_TIME ( DATE, TIME, ZONE, VALUES )
      STRNG(1:2) = TIME(1:2)
      STRNG(4:5) = TIME(3:4)
      STRNG(7:8) = TIME(5:6)
!
! Dummy alternative
!
      RETURN
!/
!/ End of WWTIME ----------------------------------------------------- /
!/
      END SUBROUTINE WWTIME
!/ ------------------------------------------------------------------- /
      SUBROUTINE EXTCDE ( IEXIT )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         06-Jan-1999 |
!/                  +-----------------------------------+
!/
!/    06-Jan-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    23-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Perfor a program stop with an exit code.
!
!  2. Method :
!
!     Machine dependent.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IEXIT   Int.   I   Exit code to be used.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!  5. Called by :
!
!     Any.
!
!  9. Switches :
!
!     !/MPI  MPI finalize interface if active
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      use shr_sys_mod
      IMPLICIT NONE
!
      INCLUDE "mpif.h"
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IEXIT
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER                 :: IERR_MPI
      LOGICAL                 :: RUN
!/
!/ Test if MPI needs to be closed
!/
      CALL shr_sys_abort('WW3 EXTCDE abort')
      CALL MPI_INITIALIZED ( RUN, IERR_MPI )
      IF ( RUN ) THEN
          CALL MPI_BARRIER ( MPI_COMM_WORLD, IERR_MPI )
          CALL MPI_FINALIZE (IERR_MPI )
        END IF
!
      CALL EXIT ( IEXIT )
!/
!/ End of EXTCDE ----------------------------------------------------- /
!/
      END SUBROUTINE EXTCDE
!/ ------------------------------------------------------------------- /
      SUBROUTINE PRINIT
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         06-May-2005 !
!/                  +-----------------------------------+
!/
!/    06-May-2005 : Origination.                        ( version 3.07 )
!/
!  1. Purpose :
!
!     Initialize profilinf routine PRTIME.
!
!  2. Method :
!
!     FORTRAN 90 SYSTEM_CLOCK intrinsic routine.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      SYSTEM_CLOCK
!                Sur.    n/a    Get system time ( !/F90 )
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
!     !/F90  FORTRAN 90 specific calls.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
! -------------------------------------------------------------------- /
!
      CALL SYSTEM_CLOCK ( PRFTB )
!
      FLPROF = .TRUE.
!
      RETURN
!/
!/ End of PRINIT ----------------------------------------------------- /
!/
      END SUBROUTINE PRINIT
!/ ------------------------------------------------------------------- /
      SUBROUTINE PRTIME ( PTIME )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         06-May-2005 !
!/                  +-----------------------------------+
!/
!/    06-May-2005 : Origination.                        ( version 3.07 )
!/
!  1. Purpose :
!
!     Get wallclock time for profiling purposes.
!
!  2. Method :
!
!     FORTRAN 90 SYSTEM_CLOCK intrinsic routine.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       PTIME   Real   O   Time retrieced from system.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      SYSTEM_CLOCK
!                Sur.    n/a    Get system time ( !/F90 )
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!     Any, after PRINIT has been called.
!
!  6. Error messages :
!
!     - If no initialization, returned time equals -1.
!     - If no system clock, returned time equals -1.
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
!     !/F90  FORTRAN 90 specific calls.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(OUT)       :: PTIME
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: PRFTA, PRFINC, PRFMAX, PRFDT
!
! -------------------------------------------------------------------- /
!
      PTIME  = -1.
!
      IF ( .NOT. FLPROF ) RETURN
!
      CALL SYSTEM_CLOCK ( PRFTA, PRFINC, PRFMAX )
      IF ( PRFMAX .NE. 0 ) THEN
          IF ( PRFTA-PRFTB .GE. 0 ) THEN
              PRFDT  = PRFTA-PRFTB
            ELSE
              PRFDT  = PRFTA-PRFTB + PRFMAX
            END IF
          PTIME  = REAL(PRFDT)/REAL(PRFINC)
        END IF
!
      RETURN
!/
!/ End of PRTIME ----------------------------------------------------- /
!/
      END SUBROUTINE PRTIME
!/
!/ End of module W3SERVMD -------------------------------------------- /
!/
      END MODULE W3SERVMD
