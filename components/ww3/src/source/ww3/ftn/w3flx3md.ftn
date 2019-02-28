#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3FLX3MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         20-Apr-2010 |
!/                  +-----------------------------------+
!/
!/    05-Jul-2006 : Origination.                        ( version 3.09 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    20-Apr-2010 : Fix INTENT of UST.                ( version 3.14.1 )
!/
!/    Copyright 2009-2010 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     FLux/stress computations according Tolman and Chalikov (1996).
!     Cap on flux added compared to W3FLX2.
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3FLX3    Subr. Public   Stresses according to TC (1996).
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!     - Originally used with source term !/ST2.
!
!  6. Switches :
!
!     !/S  Enable subroutine tracing.
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLX3 ( ZWIND, DEPTH, FP, U, UDIR, UST, USTD, Z0, CD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Jan-2014 |
!/                  +-----------------------------------+
!/
!/    05-Jul-2006 : Origination.                        ( version 3.09 )
!/    20-Apr-2010 : Fix INTENT of UST.                ( version 3.14.1 )
!/    10-Jan-2014 : Add max on division by UST          ( version 4.18 )
!/                  (This was already done for W3FLX2 on 16 Sep 2011)
!/    10-Jan-2014 : Add a min value for FP              ( version 4.18 )
!/
!  1. Purpose :
!
!     FLux/stress computations according Tolman and Chalikov (1996).
!     Cap on flux added compared to W3FLX2.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       ZWIND   Real   I   Hight of wind.
!       DEPTH   Real   I   Depth.
!       FP      Real   I   Peak frequency.
!       U       Real   I   Wind speed.
!       UDIR    Real   I   Wind direction.
!       UST     Real   O   Friction velocity.
!       USTD    Real   0   Direction of friction velocity.
!       Z0      Real   O   z0 in profile law.
!       CD      Real   O   Drag coefficient.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SRCE    Subr. W3SRCEMD Source term integration.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3GDATMD, ONLY: NITTIN, CINXSI, CD_MAX, CAP_ID
      USE W3ODATMD, ONLY: NDSE, IAPROC, NAPERR
      USE W3SERVMD, ONLY: EXTCDE
!/S      USE W3SERVMD, ONLY: STRACE
      USE W3DISPMD, ONLY: DSIE, N1MAX, EWN1
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: ZWIND, DEPTH, FP, U, UDIR
      REAL, INTENT(INOUT)     :: UST
      REAL, INTENT(OUT)       :: USTD, Z0, CD
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I1, ITT
!/S      INTEGER, SAVE           :: IENT = 0
      REAL                    :: SQRTH, SIX, R1, WNP, CP, UNZ, ALPHA, &
                                 RDCH, AFP
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3FLX3')
!
! 1.  Peak phase velocity -------------------------------------------- *
!
! ----- start of inlined and reduced WAVNU1 -----
!
        AFP    = TPI * MAX ( FP, 0.001) 
!
        SQRTH  = SQRT ( DEPTH )
        SIX    = AFP * SQRTH
        I1     = INT ( SIX / DSIE )
        IF (I1.LE.N1MAX) THEN
            R1 = SIX/DSIE - REAL(I1)
            WNP    = ( (1.-R1)*EWN1(I1) + R1*EWN1(I1+1) ) / DEPTH
          ELSE
            WNP    = AFP * AFP / GRAV
          END IF
!
! -----  end of inlined and reduced WAVNU1  -----
!
      CP     = AFP / WNP
!
! 2.  Itterative stress computation ---------------------------------- *
!
      UNZ    = MAX ( 0.01 , U )
      USTD   = UDIR
!
      DO ITT=1, NITTIN
        ALPHA  = 0.57 / ( CP / MAX (UST,0.0001) )**(1.5)
        RDCH   = MAX ( 0. ,                                           &
            LOG ( ( ZWIND * GRAV) / ( CINXSI * SQRT(ALPHA) * UNZ**2) ) )
        CD     = 0.001 * ( 0.021 + 10.4 / (RDCH**1.23+1.85) )
        UST    = SQRT(CD) * UNZ
        Z0    = ZWIND * EXP ( -0.4 / SQRT(CD) )
        END DO
!
! 3.  Apply limit to drag coefficient -------------------------------- *
!
      IF ( CAP_ID .EQ. 0 ) THEN
          CD     = MIN ( CD_MAX, CD )
        ELSE
          CD     = CD_MAX * TANH ( CD / CD_MAX )
        END IF
!
      UST    = SQRT(CD) * UNZ
      Z0     = ZWIND * EXP ( -0.4 / SQRT(CD) )
!
      RETURN
!
! Formats
!
!/
!/ End of W3FLX3 ----------------------------------------------------- /
!/
      END SUBROUTINE W3FLX3
!/
!/ End of module W3FLX3MD -------------------------------------------- /
!/
      END MODULE W3FLX3MD
