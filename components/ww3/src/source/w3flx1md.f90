!/ ------------------------------------------------------------------- /
      MODULE W3FLX1MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    03-Jul-2006 : Origination.                        ( version 3.09 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Flux/stress computations according to Wu (1980)
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3FLX1    Subr. Public   Stresses according to Wu (1980).
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
!     - Originally used with source term !/ST1.
!
!  6. Switches :
!
!     !/S  Enable subroutine tracing.
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3FLX1 ( ZWND, U10, U10D, UST, USTD, Z0, CD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         03-Jul-2006 |
!/                  +-----------------------------------+
!/
!/    03-Jul-2006 : Origination.                        ( version 3.09 )
!/
!  1. Purpose :
!
!     FLux/stress computations according to Wu (1980)
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       ZWND    Real   I   Wind height.
!       U10     Real   I   Wind speed.
!       U10D    Real   I   Wind direction.
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
      USE W3ODATMD, ONLY: NDSE, IAPROC, NAPERR
      USE W3SERVMD, ONLY: EXTCDE
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: ZWND, U10, U10D
      REAL, INTENT(OUT)       :: UST, USTD, Z0, CD
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Tests ---------------------------------------------------------- *
!
      IF ( ABS(ZWND-10.) .GT. 0.01 ) THEN
          IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) ZWND
          CALL EXTCDE (1)
        END IF
!
! 2.  Computation ---------------------------------------------------- *
!
      CD     = 0.001 * (0.8+0.065*U10)
      Z0     = ZWND * EXP ( -0.4 / SQRT(CD) )
      UST    = U10 * SQRT(CD)
      USTD   = U10D
!
      RETURN
!
! Formats
!
 1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3STR1 : '/               &
               '     HIGHT OF WIND SHOULD BE 10m IN THIS APPRACH '/   &
               '     ZWND =',F8.2,'m'/)
!/
!/ End of W3FLX1 ----------------------------------------------------- /
!/
      END SUBROUTINE W3FLX1
!/
!/ End of module INFLX1MD -------------------------------------------- /
!/
      END MODULE W3FLX1MD
