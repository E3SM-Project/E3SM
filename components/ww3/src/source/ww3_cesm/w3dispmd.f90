!/ ------------------------------------------------------------------- /
      MODULE W3DISPMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    30-Nov-1999 : Fortran 90 version.                 ( version 2.00 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     A set of routines for solving the dispersion relation.
!
!  2. Variables and types :
!
!     All variables are retated to the interpolation tables. See
!     DISTAB for a more comprehensive description.
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NAR1D     I.P.  Public   Nmmer of elements in interpolation
!                               array.
!      DFAC      R.P.  Public   Value of KH at deep boundary.
!      EWN1      R.A.  Public   Wavenumber array.
!      ECG1      R.A.  Public   Group velocity array.
!      N1MAX     Int.  Public   Actual maximum position in array.
!      DSIE      Real  Public   SI step.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      WAVNU1    Subr. Public   Solve dispersion using lookup table.
!      WAVNU2    Subr. Public   Solve dispersion relation itteratively.
!      DISTAB    Subr. Public   Fill interpolation tables.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing            ( !/S )
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!       !/S   Enable subroutine tracing.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      PUBLIC
!/
!/ Set up of public interpolation table ------------------------------ /
!/
      INTEGER, PARAMETER      :: NAR1D  =  121
      REAL, PARAMETER         :: DFAC   =    6.
!/
      INTEGER                 :: N1MAX
      REAL                    :: ECG1(0:NAR1D), EWN1(0:NAR1D), DSIE
!/
!/ Set up of public subroutines -------------------------------------- /
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE WAVNU1 (SI,H,K,CG)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         30-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    04-Nov-1990 : Final FORTRAN 77                    ( version 1.18 )
!/    30-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Calculate wavenumber and group velocity from the interpolation
!     array filled by DISTAB from a given intrinsic frequency and the
!     waterdepth.
!
!  2. Method :
!
!     Linear interpolation from one-dimensional array.
!
!  3. Parameters used :
!
!     Parameter list
!     ----------------------------------------------------------------
!       SI      Real   I   Intrinsic frequency (moving frame)  (rad/s)
!       H       Real   I   Waterdepth                            (m)
!       K       Real   O   Wavenumber                          (rad/m)
!       CG      Real   O   Group velocity                       (m/s)
!     ----------------------------------------------------------------
!
!  4. Error messages :
!
!     - None.
!
!  5. Called by :
!
!     - Any main program
!
!  6. Subroutines used :
!
!     - None
!
!  7. Remarks :
!
!     - Calculated si* is always made positive without checks : check in
!       main program assumed !
!     - Depth is unlimited.
!
!  8. Structure :
!
!     +---------------------------------------------+
!     | calculate non-dimensional frequency         |
!     |---------------------------------------------|
!     | T            si* in range ?               F |
!     |----------------------|----------------------|
!     | calculate k* and cg* | deep water approx.   |
!     | calculate output     |                      |
!     |      parameters      |                      |
!     +---------------------------------------------+
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE W3CONSTANTS, ONLY : GRAV
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: SI, H
      REAL, INTENT(OUT)       :: K, CG
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I1, I2
      REAL                    :: SQRTH, SIX, R1, R2
!/
!/ ------------------------------------------------------------------- /
!/
!
      SQRTH  = SQRT(H)
      SIX    = SI * SQRTH
      I1     = INT(SIX/DSIE)
!
      IF (I1.LE.N1MAX) THEN
          I2 = I1 + 1
          R1 = SIX/DSIE - REAL(I1)
          R2 = 1. - R1
          K  = ( R2*EWN1(I1) + R1*EWN1(I2) ) / H
          CG = ( R2*ECG1(I1) + R1*ECG1(I2) ) * SQRTH
        ELSE
          K  = SI*SI/GRAV
          CG = 0.5 * GRAV / SI
        END IF
!
      RETURN
!/
!/ End of WAVNU1 ----------------------------------------------------- /
!/
      END SUBROUTINE WAVNU1
!/ ------------------------------------------------------------------- /
      SUBROUTINE WAVNU2 (W,H,K,CG,EPS,NMAX,ICON)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         30-Nov-1999 |
!/                  +-----------------------------------+
!/
!/    17-Jul-1990 : Final FORTRAN 77                    ( version 1.18 )
!/    30-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Calculation of wavenumber K from a given angular
!     frequency W and waterdepth H.
!
!  2. Method :
!
!     Used equation :
!                        2
!                       W  = G*K*TANH(K*H)
!
!     Because of the nature of the equation, K is calculated
!     with an itterative procedure.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       W       Real   I   Angular frequency
!       H       Real   I   Waterdepth
!       K       Real   O   Wavenumber ( same sign as W )
!       CG      Real   O   Group velocity (same sign as W)
!       EPS     Real   I   Wanted max. difference between K and Kold
!       NMAX    Int.   I   Max number of repetitions in calculation
!       ICON    Int.   O   Contol counter ( See error messages )
!     ----------------------------------------------------------------
!
!  9. Switches :
!
!     !/S  Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE W3CONSTANTS, ONLY : GRAV
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NMAX
      INTEGER, INTENT(OUT)    :: ICON
      REAL, INTENT(IN)        :: W, H, EPS
      REAL, INTENT(OUT)       :: CG, K
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I
      REAL                    :: F, W0, FD, DIF, RDIF, KOLD
!/
!/ ------------------------------------------------------------------- /
!/
!
!     Initialisations :
!
      CG   = 0
      KOLD = 0
      ICON = 0
      W0   = ABS(W)
!
!     1st approach :
!
      IF (W0.LT.SQRT(GRAV/H)) THEN
          K = W0/SQRT(GRAV*H)
        ELSE
          K = W0*W0/GRAV
        END IF
!
!     Refinement :
!
      DO I=1, NMAX
        DIF = ABS(K-KOLD)
        IF (K.NE.0) THEN
            RDIF = DIF/K
          ELSE
            RDIF = 0
          END IF
        IF (DIF .LT. EPS .AND. RDIF .LT. EPS) THEN
            ICON = 1
            GOTO 100
          ELSE
            KOLD = K
            F    = GRAV*KOLD*TANH(KOLD*H)-W0**2
            IF (KOLD*H.GT.25) THEN
                FD = GRAV*TANH(KOLD*H)
              ELSE
                FD = GRAV*TANH(KOLD*H) + GRAV*KOLD*H/((COSH(KOLD*H))**2)
              END IF
            K    = KOLD - F/FD
          END IF
        END DO
!
      DIF   = ABS(K-KOLD)
      RDIF  = DIF/K
      IF (DIF .LT. EPS .AND. RDIF .LT. EPS) ICON = 1
 100  CONTINUE
      IF (2*K*H.GT.25) THEN
          CG = W0/K * 0.5
        ELSE
          CG = W0/K * 0.5*(1+(2*K*H/SINH(2*K*H)))
        END IF
      IF (W.LT.0.0) THEN
          K  = (-1)*K
          CG = CG*(-1)
        END IF
!
      RETURN
!/
!/ End of WAVNU2 ----------------------------------------------------- /
!/
      END SUBROUTINE WAVNU2
!/ ------------------------------------------------------------------- /
      SUBROUTINE DISTAB
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         30-Nov-1990 |
!/                  +-----------------------------------+
!/
!/    04-Nov-1990 : Final FORTRAN 77                    ( version 1.18 )
!/    30-Nov-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/
!  1. Purpose :
!
!     Fill interpolation arrays for the calculation of wave parameters
!     according to the linear (Airy) wave theory given the intrinsic
!     frequency.
!
!  2. Method :
!
!     For a given set of non-dimensional frequencies the interpolation
!     arrays with non-dimensional depths and group velocity are filled.
!     The following non-dimensional parameters are used :
!
!       frequency   f*SQRT(h/g) = f*
!       depth       kh          = k*
!       group vel.  c/SQRT(gh)  = c*
!
!     Where k is the wavenumber, h the depth f the intrinsic frequency,
!     g the acceleration of gravity and c the group velocity.
!
!  3. Parameters :
!
!     See module documentation.
!
!  4. Error messages :
!
!     - None.
!
!  5. Called by :
!
!     - W3GRID
!     - Any main program.
!
!  6. Subroutines used :
!
!     - WAVNU2 (solve dispersion relation)
!
!  7. Remarks :
!
!     - In the filling of the arrays H = 1. is assumed and the factor
!       SQRT (g) is moved from the interpolation to the filling
!       procedure thus :
!
!         k* = k
!
!         c* = cg/SQRT(g)
!
!  8. Structure
!
!     -----------------------------------
!       include common block
!       calculate parameters
!       fill zero-th position of arrays
!       fill middle positions of arrays
!       fill last positions of arrays
!     -----------------------------------
!
!  9. Switches :
!
!       !/S   Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE W3CONSTANTS, ONLY : GRAV
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, ICON
      REAL                    :: DEPTH, CG, SIMAX, SI, K
!/
!/ ------------------------------------------------------------------- /
!/
!
! Calculate parameters ----------------------------------------------- *
!
      N1MAX  = NAR1D - 1
      DEPTH  = 1.
      SIMAX  = SQRT (GRAV * DFAC)
      DSIE   = SIMAX / REAL(N1MAX)
!
! Fill zero-th position of arrays ------------------------------------ *
!
      EWN1(0) = 0.
      ECG1(0) = SQRT(GRAV)
!
! Fill middle positions of arrays ------------------------------------ *
!
      DO I=1, N1MAX
        SI = REAL(I)*DSIE
        CALL WAVNU2 (SI,DEPTH,K,CG,1E-7,15,ICON)
        EWN1(I) = K
        ECG1(I) = CG
        END DO
!
! Fill last positions of arrays -------------------------------------- *
!
      I      = N1MAX+1
      SI     = REAL(I)*DSIE
      CALL WAVNU2 (SI,DEPTH,K,CG,1E-7,15,ICON)
      EWN1(I) = K
      ECG1(I) = CG
!
      RETURN
!/
!/ End of DISTAB ----------------------------------------------------- /
!/
      END SUBROUTINE DISTAB
!/
!/ End of module W3DISPMD -------------------------------------------- /
!/
      END MODULE W3DISPMD
