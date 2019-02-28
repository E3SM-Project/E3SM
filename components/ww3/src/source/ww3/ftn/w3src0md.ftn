#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SRC0MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    05-Jul-2006 : Origination.                        ( version 3.09 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Mean wave parameter computation for case without input and
!     dissipation.
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SPR0    Subr. Public   Mean parameters from spectrum.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.            ( !/S )
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!       !/S      Enable subroutine tracing.
!       !/T      Test output, see subroutines.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SPR0 (A, CG, WN, EMEAN, FMEAN, WNMEAN, AMAX)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         05-Jul-2006 |
!/                  +-----------------------------------+
!/
!/    05-Jul-2006 : Origination.                        ( version 3.09 )
!/
!  1. Purpose :
!
!     Calculate mean wave parameters.
!
!  2. Method :
!
!     See source term routines.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action as a function of direction and
!                         wavenumber.
!       CG      R.A.  I   Group velocities.
!       WN      R.A.  I   Wavenumbers.
!       EMEAN   Real  O   Mean wave energy.
!       FMEAN   Real  O   Mean wave frequency.
!       WNMEAN  Real  O   Mean wavenumber.
!       AMAX    Real  O   Maximum action density in spectrum.
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
!     !/T  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, SIG, DDEN, FTE, FTF, FTWN
!/T      USE W3ODATMD, ONLY: NDST
!/S      USE W3SERVMD, ONLY: STRACE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NTH,NK), CG(NK), WN(NK)
      REAL, INTENT(OUT)       :: EMEAN, FMEAN, WNMEAN, AMAX
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IK, ITH
!/S      INTEGER, SAVE           :: IENT = 0
      REAL                    :: EB(NK), EBAND
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SPR0')
!
      EMEAN  = 0.
      FMEAN  = 0.
      WNMEAN = 0.
      AMAX   = 0.
!
! 1.  Integral over directions
!
      DO IK=1, NK
        EB(IK) = 0.
        DO ITH=1, NTH
          EB(IK) = EB(IK) + A(ITH,IK)
          AMAX   = MAX ( AMAX , A(ITH,IK) )
          END DO
        END DO
!
! 2.  Integrate over directions
!
      DO IK=1, NK
        EB(IK) = EB(IK) * DDEN(IK) / CG(IK)
        EMEAN  = EMEAN  + EB(IK)
        FMEAN  = FMEAN  + EB(IK) / SIG(IK)
        WNMEAN = WNMEAN + EB(IK) / SQRT(WN(IK))
        END DO
!
! 3.  Add tail beyond discrete spectrum
!     ( DTH * SIG absorbed in FTxx )
!
      EBAND  = EB(NK) / DDEN(NK)
      EMEAN  = EMEAN  + EBAND * FTE
      FMEAN  = FMEAN  + EBAND * FTF
      WNMEAN = WNMEAN + EBAND * FTWN
!
! 4.  Final processing
!
      FMEAN  = TPIINV * EMEAN / MAX ( 1.E-7 , FMEAN )
      WNMEAN = ( EMEAN / MAX ( 1.E-7 , WNMEAN ) )**2
!
!/T      WRITE (NDST,9000) EMEAN, FMEAN, WNMEAN
!
      RETURN
!
! Formats
!
!/T 9000 FORMAT (' TEST W3SPR0 : E,F,WN MEAN ',3E10.3)
!/
!/ End of W3SPR0 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SPR0
!/
!/ End of module W3SRC0MD -------------------------------------------- /
!/
      END MODULE W3SRC0MD
