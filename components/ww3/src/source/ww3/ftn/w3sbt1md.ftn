#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SBT1MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    For updates see W3SBT1 documentation.
!/
!  1. Purpose :
!
!     JONSWAP bottom friction routine.
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SBT1    Subr. Public   JONSWAP source term.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!     See subroutine documentation.
!
!  5. Remarks :
!
!  6. Switches :
!
!     See subroutine documentation.
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SBT1 (A, CG, WN, DEPTH, S, D)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    05-Dec-1996 : Final FORTRAN 77.                   ( version 1.18 )
!/    08-Dec-1999 : Upgrade to FORTRAN 90.              ( version 2.00 )
!/    20-Dec-2004 : Multiple model version.             ( version 3.06 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Bottom friction source term according to the empirical JONSWAP
!     formulation.
!
!  2. Method :
!
!              2 GAMMA   /    CG         \      SBTC1 /     \       .
!       Sbt = ---------- | ------- - 0.5 | E  = ----- | ... | E    (1)
!             GRAV DEPTH \  SI/WN        /      DEPTH \     /
!
!     Where GAMMA = -0.038 m2/s3 (JONSWAP)
!                 = -0.067 m2/s3 (Bouws and Komen 1983)
!
!     In the routine, the constant 2 GAMMA / GRAV = SBTC1.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D)
!       CG      R.A.  I   Group velocities.
!       WN      R.A.  I   Wavenumbers.
!       DEPTH   Real  I   Mean water depth.
!       S       R.A.  O   Source term (1-D version).
!       D       R.A.  O   Diagonal term of derivative (1-D version).
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing (!/S switch).
!      PRT2DS    Subr. W3ARRYMD Print plot output (!/T1 switch).
!      OUTMAT    Subr. W3ARRYMD Matrix output (!/T2 switch).
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SRCE    Subr. W3SRCEMD Source term integration.
!      W3EXPO    Subr.   N/A    Point output post-processor.
!      GXEXPO    Subr.   N/A    GrADS point output post-processor.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     None.
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
!     !/T   Enable general test output.
!     !/T0  2-D print plot of source term.
!     !/T1  Print arrays.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, MAPWN, SBTC1
!/T      USE W3ODATMD, ONLY: NDST
!/S      USE W3SERVMD, ONLY: STRACE
!/T0      USE W3ARRYMD, ONLY: PRT2DS
!/T1      USE W3ARRYMD, ONLY: OUTMAT
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: CG(NK), WN(NK), DEPTH, A(NSPEC)
      REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS, IK, NSCUT
!/S      INTEGER, SAVE           :: IENT = 0
!/T0      INTEGER                 :: ITH
      REAL                    :: FACTOR, CBETA(NK)
!/T0      REAL                    :: DOUT(NK,NTH)
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SBT1')
!
! 1.  Deep water ===================================================== *
!
      IF ( DEPTH*WN(1) .GT. 6 ) THEN
!
          D = 0.
          S = 0.
!
! 2.  Shallow water ================================================== *
!
        ELSE
!
! 2.a Set constant
!
          FACTOR = SBTC1 / DEPTH
!
!/T          WRITE (NDST,9000) FACTOR, DEPTH
!
! 2.b Wavenumber dependent part.
!
          DO IK=1, NK
            IF ( WN(IK)*DEPTH .GT. 6. ) EXIT
            CBETA(IK) = FACTOR *                                      &
               MAX(0., (CG(IK)*WN(IK)/SIG(IK)-0.5) )
            END DO
!
! 2.c Fill diagional matrix
!
          NSCUT  = (IK-1)*NTH
!
          DO IS=1, NSCUT
            D(IS) = CBETA(MAPWN(IS))
            END DO
!
          DO IS=NSCUT+1, NSPEC
            D(IS) = 0.
            END DO
!
          S = D * A
!
        END IF
!
! ... Test output of arrays
!
!/T0      DO IK=1, NK
!/T0        DO ITH=1, NTH
!/T0          DOUT(IK,ITH) = D(ITH+(IK-1)*NTH)
!/T0          END DO
!/T0        END DO
!
!/T0      CALL PRT2DS (NDST, NK, NK, NTH, DOUT, SIG(1:), '  ', 1.,    &
!/T0                         0.0, 0.001, 'Diag Sbt', ' ', 'NONAME')
!
!/T1      CALL OUTMAT (NDST, D, NTH, NTH, NK, 'diag Sbt')
!
      RETURN
!
! Formats
!
!/T 9000 FORMAT (' TEST W3SBT1 : FACTOR, DEPTH  : ',2E10.3)
!/
!/ End of W3SBT1 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SBT1
!/
!/ End of module W3SBT1MD -------------------------------------------- /
!/
      END MODULE W3SBT1MD
