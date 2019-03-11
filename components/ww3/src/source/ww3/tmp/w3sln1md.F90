#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SLN1MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    23-Jun-2006 : Origination.                        ( version 3.09 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Linear wind input according to Cavaleri and Melanotte-Rizzoli
!     (1982) filtered for low frequencies according to Tolman (1992).
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SLN1    Subr. Public   User supplied linear input.
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
!  6. Switches :
!
!     !/S  Enable subroutine tracing.
!     !/T  Test output.
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SLN1 (K, FHIGH, USTAR, USDIR, S)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         23-Jun-2006 |
!/                  +-----------------------------------+
!/
!/    23-Jun-2006 : Origination.                        ( version 3.09 )
!/
!  1. Purpose :
!
!     Linear wind input according to Cavaleri and Melanotte-Rizzoli
!     (1982) filtered for low frequencies according to Tolman (1992).
!
!  2. Method :
!
!     The expression of Cavaleri and Melanotte-Rizzoli, converted to
!     action spectra defined in terms of wavenumber and direction
!     becomes
!
!                       -1       /     /                \ \ 4
!       Sln  = SLNC1 * k   * max | 0., | U* cos(Dtheta) | |        (1)
!                                \     \                / /
!
!                             2     -2
!              SLNC1 = 80 RHOr  GRAV   FILT                        (2)
!
!     Where :
!
!        RHOr     Density of air dev. by density of water.
!        U*       Wind friction velocity.
!        Dtheta   Difference in wind and wave direction.
!        FILT     Filter based on PM and cut-off frequencies.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       K       R.A.  I   Wavenumber for entire spectrum.
!       FHIGH   R.A.  I   Cut-off frequency in integration (rad/s)
!       USTAR   Real  I   Friction velocity.
!       USDIR   Real  I   Direction of USTAR.
!       S       R.A.  O   Source term.
!     ----------------------------------------------------------------
!                         *) Stored as 1-D array with dimension NTH*NK
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
!     !/T  Test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3GDATMD, ONLY: NTH, NK, ECOS, ESIN, SIG, SLNC1, FSPM, FSHF
      USE W3ODATMD, ONLY: NDSE, NDST
      USE W3SERVMD, ONLY: EXTCDE
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: K(NK), FHIGH, USTAR, USDIR
      REAL, INTENT(OUT)       :: S(NTH,NK)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ITH, IK
      REAL                    :: COSU, SINU, DIRF(NTH), FAC, FF1, FF2, &
                                 FFILT, RFR, WNF(NK)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Set up factors ------------------------------------------------- *
!
      COSU   = COS(USDIR)
      SINU   = SIN(USDIR)
!
      DO ITH=1, NTH
        DIRF(ITH) = MAX ( 0. , (ECOS(ITH)*COSU+ESIN(ITH)*SINU) )**4
        END DO
!
      FAC    = SLNC1 * USTAR**4
      FF1    = FSPM * GRAV/(28.*USTAR)
      FF2    = FSHF * MIN(SIG(NK),FHIGH)
      FFILT  = MIN ( MAX(FF1,FF2) , 2.*SIG(NK) )
      DO IK=1, NK
        RFR    = SIG(IK) / FFILT
        IF ( RFR .LT. 0.5 ) THEN
            WNF(IK) = 0.
          ELSE
            WNF(IK) = FAC / K(IK) * EXP(-RFR**(-4))
          END IF
        END DO
!
! 2.  Compose source term -------------------------------------------- *
!
      DO IK=1, NK
        S(:,IK) = WNF(IK) * DIRF(:)
        END DO
!
      RETURN
!
! Formats
!
!/
!/ End of W3SLN1 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SLN1
!/
!/ End of module INSLN1MD -------------------------------------------- /
!/
      END MODULE W3SLN1MD
