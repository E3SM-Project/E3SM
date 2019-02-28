#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SNL3MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Jul-2012 |
!/                  +-----------------------------------+
!/
!/    21-Jul-2008 : Origination as NLX option.          ( version 3.13 )
!/    03-Jan-2009 : Bug fixes INSNLX.                   ( version 3.13 )
!/                  See remarks section for module.
!/    25-Aug-2009 : Conversion to F(f,theta) form.      ( version 3.13 )
!/    13-Nov-2009 : Bug fix DELTH in initialization.    ( version 3.13 )
!/    01-Dec-2009 : Bug fix frequency filtering.        ( version 3.13 )
!/    13-Aug-2010 : Move to NL3.                        ( version 3.15 )
!/    13-Jul-2012 : Moved from version 3.15 to 4.08.    ( version 4.08 )
!/
!/    Copyright 2008-2012 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!
!  1. Purpose :
!
!     Generalized and optimized multiple DIA implementation.
!     Expressions in terms of original F(f,theta) spectrum.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NKD       I.P.  Private  Number of nondimensional depths in
!                               storage array.
!      KDMIN     R.P.  Private  Minimum relative depth in table.
!      KDMAX     R.P.  Private  Maximum relative depth in table.
!      LAMMAX    R.P.  Public   Maximum value for lambda or mu.
!      DELTHM    R.P.  Public   Maximum angle gap (degree).
!      SITMIN    Real  Private  Minimum nondimensional radian 
!                               frequency in table.
!      XSIT      Real  Private  Corresponding increment factor.
!     ----------------------------------------------------------------
!
!     See W3SNL3 and INSNL3 for documentation of variables in W3GDATMD
!     as used here.
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SNL3    Subr. Public   Multiple DIA for arbitrary depth.
!      EXPAND    Subr. W3SNL3   Expand spectrum for indirect address.
!      EXPND2    Subr. W3SNL3   Expand Snl and D contributions.
!      INSNL3    Subr. Public   Corresponding initialization routine.
!      MINLAM    R.F.  INSNL3   Minimum lambda for quadruplet.
!      MAXLAM    R.F.  INSNL3   Maximum lambda for quadruplet.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr. W2SERVMD Program abort.
!      WAVNU1    Subr. W3DISPMD Solve dispersion relation.
!      WAVNU2    Subr. W3DISPMD Solve dispersion relation.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!     - Filtering techniques for computation of quadruplet spectral
!       values and distribution in spectral space have been tested
!       but were not found worth the large coding effort involved.
!     - WAVNU1 is used in W3SNL3  for consistency with spectral grid 
!       description.
!     - WAVNU2 is used in INSNL3  for accuracy in the computation of
!       the layout of the quadruplets (higher computational cost is
!       not an issue with initialization routine).
!     - For large lambda or mu the original maximum kd = 10. still 
!       leads to significantly different quadruplet layout in 
!       secion 3. To remedy this, the orriginal settings of the
!       lookup tables
!
!     INTEGER, PRIVATE, PARAMETER :: NKD = 250
!     REAL, PRIVATE, PARAMETER    :: KDMIN = 0.025 ,  KDMAX = 10.
!
!       was reset to
!
!     INTEGER, PRIVATE, PARAMETER :: NKD = 275
!     REAL, PRIVATE, PARAMETER    :: KDMIN = 0.025 ,  KDMAX = 20.
!
!       for the bug fix of 03-Jan-2009. Note that with this, the 
!       estimate of NTHMAX in INSNL3 also is needed to guarantee 
!       consistent NTHMAX and NTHM2 for any lambda and mu.
!
!  6. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/Tn   Test output (see main subroutines).
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER, PRIVATE, PARAMETER :: NKD = 275
      REAL, PRIVATE, PARAMETER    :: KDMIN = 0.025 ,  KDMAX = 20.
      REAL, PUBLIC, PARAMETER     :: LAMMAX = 0.49999
      REAL, PUBLIC, PARAMETER     :: DELTHM = 90.
!
      REAL, PRIVATE               :: SITMIN, XSIT
!
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SNL3 (  A, CG, WN, DEPTH, S, D )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         01-Dec-2009 |
!/                  +-----------------------------------+
!/
!/    21-Jul-2008 : Origination as NLX option.          ( version 3.13 )
!/    25-Aug-2009 : Conversion to F(f,theta) form.      ( version 3.13 )
!/    01-Dec-2009 : Bug fix frequency filtering.        ( version 3.13 )
!/
!  1. Purpose :
!
!     Multiple Discrete Interaction Parameterization for arbitrary
!     depths with generalized quadruplet layout.
!
!  2. Method :
!
!     This is a direct implementation of the Discrete Interaction
!     Paramterization (DIA) with multiple representative quadruplets
!     (MDIA) for arbitrary water depths.
!
!     The outer loop of the code is over quadruplet realizations, 
!     which implies two realizations for a conventional quadruplet
!     definitions and four for extended definitions (with rescaling 
!     of the contants for consistency). Within this loop the compu-
!     tations are performed in two stages. First, interactions
!     contributions are computed for the entire spectral space, 
!     second all contributions are combined into the actual inter-
!     actions and diagonal contributions.
!
!     Arbitrary depths are addressed by generating a lookup table
!     for the relative depth. These tables are used for each discrete
!     frequency separately. Efficient memory usages requires relative
!     addressing to reduce the size of the lookup tables. To use this
!     the spectral space is expanded to higher and lower frequencies
!     as well as directional space is expanded/volded. This is done
!     for the input (pseudo-) spectrum (action spectrum devided by the
!     wavenumber) to determine spectral densities at the quadruplet
!     components, and the spectral space describing individual contri-
!     butions before they are combined into the actual interactions.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action spectrum A(ITH,IK) as a function of
!                         direction (rad)  and wavenumber.
!       CG      R.A.  I   Group velocities (dimension NK).
!       WN      R.A.  I   Wavenumbers (dimension NK).
!       DEPTH   Real  I   Water depth in meters.
!       S       R.A.  O   Source term.
!       D       R.A.  O   Diagonal term of derivative.
!     ----------------------------------------------------------------
!
!     Variables describing the expanded frequency space from the 
!     dynamic storage in w3gdatmd.
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NFR       Int.  Public   Number of frequencies or wavenumbers
!                               in discrete spectral space (NFR=>NK).
!      NFRMIN    Int.  Public   Minimum discrete frequency in the 
!                               expanded frequency space.
!      NFRMAX    Int.  Public   Idem maximum for first part.
!      NFRCUT    Int.  Public   Idem maximum for second part.
!      NTHMAX    Int.  Public   Extension of directional space.
!      NTHEXP    Int   Public   Number of bins in extended dir. space.
!      NSPMIN, NSPMAX, NSPMX2
!                Int.  Public   1D spectral space range.
!      FRQ       R.A.  Public   Expanded frequency range (Hz).
!      XSI       R.A.  Public   Expanded frequency range (rad/s).
!     ----------------------------------------------------------------
!
!     Variables describing lookup tables.
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NQA       Int.  Public   Number of actual quadruplets.
!      QST1      I.A.  Public   Spectral offsets for compuation of
!                               quadruplet spectral desnities.
!      QST2      R.A.  Public   Idem weights.
!      QST3      R.A.  Public   Norm. factors in product term and
!                               in diagonal strength.
!      QST4      I.A.  Public   Spectral offsets for combining of 
!                               interactions and diagonal.
!      QST5      R.A.  Public   Idem weights for interactions.
!      QST6      R.A.  Public   Idem weights for diagonal.
!     ----------------------------------------------------------------
!
!     Variables describing model setup.
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      SNLMSC    Real  Public   Tuning power 'deep' scaling.
!      SNLNSC    Real  Public   Tuning power 'shallow' scaling.
!      SNLSFD    Real  Public   'Deep' nondimensional filer freq.
!      SNLSFS    Real  Public   'Shallow' nondimensional filer freq.
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
!     - Note that this code uses explicit unroling of potential loop
!       structures for optimization purposes.
!     - Normalization with respect to the number of quadruplets is
!       included in the proportionality constant.
!     - Note that the outer loop in the routine considers one actual
!       quadruplet realization per loop cycle. For the traditional 
!       quadruplet layout two realizations occure, for the expanded 
!       four realizations occur. For consistency, strength of a
!       traditional layout is therefore doubled.
!     - 1D representation is used of 2D spectral space for optimization
!       purposes.
!     - Contributions are first computed in the convetional spectral
!       space and are then expancded "in place" into the expanded
!       spectral space in EXPND2.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S    Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3GDATMD, ONLY: NFR => NK, NTH, SIG, FACHFE, FACTI1, FACTI2,&
                          NFRMIN, NFRMAX, NFRCUT, NTHMAX, NTHEXP,     &
                          NSPMIN, NSPMAX, NSPMX2, FRQ, XSI, NQA,      &
                          QST1, QST2, QST3, QST4, QST5, QST6, SNLMSC, &
                          SNLNSC, SNLSFD, SNLSFS
      USE W3ODATMD, ONLY: NDSE, NDST
!
      USE W3SERVMD, ONLY: EXTCDE
      USE W3DISPMD, ONLY: WAVNU1
!/S      USE W3SERVMD, ONLY: STRACE
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)   :: A(NTH,NFR), CG(NFR), WN(NFR), DEPTH
      REAL, INTENT(OUT)  :: S(NTH,NFR), D(NTH,NFR)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER            :: IFR, IERR, IKD, JKD(NFRCUT), IQA, IF1MIN, &
                            IF1MAX, IF2MIN, IF2MAX, ISP0, ISPX0, ITH, &
                            ISP, ISPX
!/S      INTEGER, SAVE      :: IENT = 0
      INTEGER            :: LQST1(16), LQST4(16)
      REAL               :: XSITLN, SIT, FPROP, FQ1, FQ2, FQ3, FQ4,   &
                            AUX1, AUX2
      REAL               :: XWN(NFRMAX), XCG(NFRMAX), SCALE1(NFRCUT), &
                            SCALE2(NFRCUT), LQST2(16), FACT(6),       &
                            LQST5(16), LQST6(16)
      REAL               :: UE(NSPMIN:NSPMAX), DSB(NSPMIN:NSPMX2),    &
                            DD1(NSPMIN:NSPMX2), DD2(NSPMIN:NSPMX2),   &
                            DD3(NSPMIN:NSPMX2), DD4(NSPMIN:NSPMX2)
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SNL3')
!
! 1.  Initialization ------------------------------------------------- *
! 1.a Constants and arrays
!
      XSITLN = LOG(XSIT)
!
      S      = 0.
      D      = 0.
!     DSB    = 0.
!     DD1    = 0.
!     DD2    = 0.
!     DD3    = 0.
!     DD4    = 0.
!
! 1.a Extended frequency range
!
      XWN(1:NFR) = WN
      XCG(1:NFR) = CG
!
      DO IFR = NFR+1, NFRMAX
        CALL WAVNU1 ( XSI(IFR), DEPTH, XWN(IFR), XCG(IFR) )
        END DO
!
! 1.b Expanded pseudo spetrum
!
      CALL EXPAND ( UE )
!
! 1.c Set up scaling functions
!
      AUX1   = 1. / ( TPI**11 * GRAV**(4.-SNLMSC) )
      AUX2   = GRAV**2 / TPI**11
!
      DO IFR=1, NFRCUT
        SCALE1(IFR) = AUX1 * XWN(IFR)**(4.+SNLMSC) *                 &
                        XSI(IFR)**(13.-2.*SNLMSC) / XCG(IFR)**2
        SCALE2(IFR) = AUX2 * XWN(IFR)**11 *                          &
                        (XWN(IFR)*DEPTH)**SNLNSC / XCG(IFR)
        END DO
!
! 1.d Set up depth scaling counters
!
      DO IFR=1, NFRCUT
        SIT      = XSI(IFR) * SQRT(DEPTH/GRAV)
        IKD      = 1 + NINT ( ( LOG(SIT) - LOG(SITMIN) ) / XSITLN )
        JKD(IFR) = MAX ( 1 , MIN(IKD,NKD) )
        END DO
!
! 2.  Base loop over quadruplet realizations ------------------------- *
!
      DO IQA=1 , NQA
!
! 3.  Obtain quadruplet energies for all spectral bins --------------- *
! 3.a Set frequency ranges
!
        AUX1   = QST3(5,IQA,1)
        AUX2   = QST3(6,IQA,1)
!
        IF1MIN = 1
        IF1MAX = NFRCUT
        IF2MIN = 1
        IF2MAX = NFR
!
        IF ( AUX1 .LE. 0. .AND. AUX2 .LE. 0. ) THEN
!
            CYCLE
!
          ELSE IF ( AUX2 .LE. 0. ) THEN
!
            SIT    = SNLSFD * SQRT(GRAV/DEPTH)
            IFR    = NINT ( FACTI2 + FACTI1*LOG(SIT) )
            IF ( IFR .GT. NFR ) CYCLE
!
            IF ( IFR .GT. 1 ) THEN
                IF1MIN = MAX ( 1 , IFR )
                IF2MIN = MAX ( 1 , IF1MIN + NFRMIN )
                DSB(1:(IF1MIN-1)*NTH) = 0.
                DD1(1:(IF1MIN-1)*NTH) = 0.
                DD2(1:(IF1MIN-1)*NTH) = 0.
                DD3(1:(IF1MIN-1)*NTH) = 0.
                DD4(1:(IF1MIN-1)*NTH) = 0.
              END IF
!
          ELSE IF ( AUX1 .LE. 0. ) THEN
!
            SIT    = SNLSFS * SQRT(GRAV/DEPTH)
            IFR    = NINT ( FACTI2 + FACTI1*LOG(SIT) )
            IF ( IFR .LT. 1 ) CYCLE
!
            IF ( IFR .LT. NFRCUT ) THEN
                IF1MAX = MIN ( NFRCUT, IFR )
!               IF2MAX = NFR
                DSB(IF1MAX*NTH+1:NFRCUT*NTH) = 0.
                DD1(IF1MAX*NTH+1:NFRCUT*NTH) = 0.
                DD2(IF1MAX*NTH+1:NFRCUT*NTH) = 0.
                DD3(IF1MAX*NTH+1:NFRCUT*NTH) = 0.
                DD4(IF1MAX*NTH+1:NFRCUT*NTH) = 0.
              END IF
!
          END IF
!
! 3.b Loop over frequencies
!
        DO IFR=IF1MIN, IF1MAX
!
! 3.c Find discrete depths
!
          IKD    = JKD(IFR)
!
! 3.d Get offsets and weights
!
          LQST1  = QST1(:,IQA,IKD)
          LQST2  = QST2(:,IQA,IKD)
          FACT   = QST3(:,IQA,IKD)
          FACT(1:4) = FACT(1:4) * XCG(IFR) / ( XWN(IFR) *XSI(IFR) )
          FPROP  = SCALE1(IFR)*FACT(5) + SCALE2(IFR)*FACT(6)
!
! 3.e Loop over directions
!
          ISP0   = (IFR-1)*NTH
          ISPX0  = (IFR-1)*NTHEXP
!
          DO ITH=1, NTH
!
            ISP    = ISP0 + ITH
            ISPX   = ISPX0 + ITH
!
            FQ1    = ( UE(ISPX+LQST1( 1)) * LQST2( 1) +               &
                       UE(ISPX+LQST1( 2)) * LQST2( 2) +               &
                       UE(ISPX+LQST1( 3)) * LQST2( 3) +               &
                       UE(ISPX+LQST1( 4)) * LQST2( 4) ) * FACT(1)
            FQ2    = ( UE(ISPX+LQST1( 5)) * LQST2( 5) +               &
                       UE(ISPX+LQST1( 6)) * LQST2( 6) +               &
                       UE(ISPX+LQST1( 7)) * LQST2( 7) +               &
                       UE(ISPX+LQST1( 8)) * LQST2( 8) ) * FACT(2)
            FQ3    = ( UE(ISPX+LQST1( 9)) * LQST2( 9) +               &
                       UE(ISPX+LQST1(10)) * LQST2(10) +               &
                       UE(ISPX+LQST1(11)) * LQST2(11) +               &
                       UE(ISPX+LQST1(12)) * LQST2(12) ) * FACT(3)
            FQ4    = ( UE(ISPX+LQST1(13)) * LQST2(13) +               &
                       UE(ISPX+LQST1(14)) * LQST2(14) +               &
                       UE(ISPX+LQST1(15)) * LQST2(15) +               &
                       UE(ISPX+LQST1(16)) * LQST2(16) ) * FACT(4)
!
            AUX1   = FQ1 * FQ2 * ( FQ3 + FQ4 )
            AUX2   = FQ3 * FQ4 * ( FQ1 + FQ2 )
            DSB(ISP) = FPROP * ( AUX1 - AUX2 )
!
            AUX1   = FQ3 + FQ4
            AUX2   = FQ3 * FQ4     
            DD1(ISP) = FPROP * FACT(1) * ( FQ2 * AUX1 - AUX2 )
            DD2(ISP) = FPROP * FACT(2) * ( FQ1 * AUX1 - AUX2 )
!
            AUX1   = FQ1 + FQ2
            AUX2   = FQ1 * FQ2     
            DD3(ISP) = FPROP * FACT(3) * ( AUX2 - FQ4*AUX1 )
            DD4(ISP) = FPROP * FACT(4) * ( AUX2 - FQ3*AUX1 )
!
! ... End loop 3.e
!
            END DO
!
! ... End loop 3.b
!
          END DO
!
! 3.e Expand arrays
!
        CALL EXPND2 ( DSB(1:NTH*NFRCUT), DSB )
        CALL EXPND2 ( DD1(1:NTH*NFRCUT), DD1 )
        CALL EXPND2 ( DD2(1:NTH*NFRCUT), DD2 )
        CALL EXPND2 ( DD3(1:NTH*NFRCUT), DD3 )
        CALL EXPND2 ( DD4(1:NTH*NFRCUT), DD4 )
!
! 4.  Put it all together -------------------------------------------- *
! 4.a Loop over frequencies
!
        DO IFR=IF2MIN, IF2MAX
!
! 4.b Find discrete depths and storage
!
          IKD    = JKD(IFR)
!
! 4.c Get offsets and weights
!
          LQST4  = QST4(:,IQA,IKD)
          LQST5  = QST5(:,IQA,IKD)
          LQST6  = QST6(:,IQA,IKD)
!
! 4.d Loop over directions
!
          ISPX0  = (IFR-1)*NTHEXP
!
          DO ITH=1, NTH
!
            ISPX   = ISPX0 + ITH
!
            S(ITH,IFR) = S(ITH,IFR) + DSB(ISPX+LQST4( 1)) * LQST5( 1) &
                                    + DSB(ISPX+LQST4( 2)) * LQST5( 2) &
                                    + DSB(ISPX+LQST4( 3)) * LQST5( 3) &
                                    + DSB(ISPX+LQST4( 4)) * LQST5( 4) &
                                    + DSB(ISPX+LQST4( 5)) * LQST5( 5) &
                                    + DSB(ISPX+LQST4( 6)) * LQST5( 6) &
                                    + DSB(ISPX+LQST4( 7)) * LQST5( 7) &
                                    + DSB(ISPX+LQST4( 8)) * LQST5( 8) &
                                    + DSB(ISPX+LQST4( 9)) * LQST5( 9) &
                                    + DSB(ISPX+LQST4(10)) * LQST5(10) &
                                    + DSB(ISPX+LQST4(11)) * LQST5(11) &
                                    + DSB(ISPX+LQST4(12)) * LQST5(12) &
                                    + DSB(ISPX+LQST4(13)) * LQST5(13) &
                                    + DSB(ISPX+LQST4(14)) * LQST5(14) &
                                    + DSB(ISPX+LQST4(15)) * LQST5(15) &
                                    + DSB(ISPX+LQST4(16)) * LQST5(16) 
!
            D(ITH,IFR) = D(ITH,IFR) + DD1(ISPX+LQST4( 1)) * LQST6( 1) &
                                    + DD1(ISPX+LQST4( 2)) * LQST6( 2) &
                                    + DD1(ISPX+LQST4( 3)) * LQST6( 3) &
                                    + DD1(ISPX+LQST4( 4)) * LQST6( 4) &
                                    + DD2(ISPX+LQST4( 5)) * LQST6( 5) &
                                    + DD2(ISPX+LQST4( 6)) * LQST6( 6) &
                                    + DD2(ISPX+LQST4( 7)) * LQST6( 7) &
                                    + DD2(ISPX+LQST4( 8)) * LQST6( 8) &
                                    + DD3(ISPX+LQST4( 9)) * LQST6( 9) &
                                    + DD3(ISPX+LQST4(10)) * LQST6(10) &
                                    + DD3(ISPX+LQST4(11)) * LQST6(11) &
                                    + DD3(ISPX+LQST4(12)) * LQST6(12) &
                                    + DD4(ISPX+LQST4(13)) * LQST6(13) &
                                    + DD4(ISPX+LQST4(14)) * LQST6(14) &
                                    + DD4(ISPX+LQST4(15)) * LQST6(15) &
                                    + DD4(ISPX+LQST4(16)) * LQST6(16)
!
! ... End loop 4.d
!
            END DO
!
! ... End loop 4.a
!
          END DO
!
! ... End of loop 2.
!
        END DO
!
! 5.  Convert back to wave action ------------------------------------ *
!
      DO IFR=IF2MIN, IF2MAX
        S(:,IFR) = S(:,IFR) / XSI(IFR) * XCG(IFR) * TPIINV
        END DO
!
      RETURN
!/
!/ Embedded subroutines
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE EXPAND ( SPEC )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         21-Aug-2009 |
!/                  +-----------------------------------+
!/
!/    03-Jul-2008 : Origination.                        ( version 3.13 )
!/    21-Aug-2009 : Conversion to F(f,theta) form.      ( version 3.13 )
!/
!  1. Purpose :
!
!     Expand spectrum, subroutine used to simplify addressing.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       SPEC    R.A.  O   Expanded spectrum.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ Parameter list
!/
      REAL, INTENT(OUT)       :: SPEC(1-NTHMAX:NTH+NTHMAX,NFRMIN:NFRMAX)
!/
!/ Local parameters
!/
      INTEGER                 :: IFR, ITH
!/
!/ ------------------------------------------------------------------- /
!
      SPEC(:,NFRMIN:0) = 0.
!
      SPEC(1:NTH,1:NFR) = A * TPI
!
      DO IFR=1, NFR
        SPEC(1:NTH,IFR) = SPEC(1:NTH,IFR) * XSI(IFR) / XCG(IFR)
        END DO
!
      DO IFR=NFR+1, NFRMAX
        SPEC(1:NTH,IFR) = SPEC(1:NTH,IFR-1) * FACHFE
        END DO
!
      DO ITH=1, NTHMAX
        SPEC(NTH+ITH,1:NFRMAX) = SPEC(   ITH   ,1:NFRMAX)
        SPEC( 1 -ITH,1:NFRMAX) = SPEC(NTH+1-ITH,1:NFRMAX)
        END DO
!
      RETURN
!/
!/ End of EXPAND ----------------------------------------------------- /
!/
      END SUBROUTINE EXPAND
!/ ------------------------------------------------------------------- /
      SUBROUTINE EXPND2 ( ARIN, AROUT )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         16-Jul-2008 |
!/                  +-----------------------------------+
!/
!/    16-Jul-2008 : Origination.                        ( version 3.13 )
!/
!  1. Purpose :
!
!     Expand spectrum to simplify indirect addressing.
!     Done 'in place' with temporary array ( ARIN = AROUT )
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       SPIN    R.A.  I   Input array.
!       SPOUT   R.A.  I   Output array.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ Parameter list
!/
      REAL, INTENT(IN)      :: ARIN(NTH,NFRCUT)
      REAL, INTENT(OUT)     :: AROUT(1-NTHMAX:NTH+NTHMAX,NFRMIN:NFRCUT)
!/
!/ Local parameters
!/
      INTEGER               :: IFR, ITH
      REAL                  :: TEMP(NTH,NFRCUT)
!/
!/ ------------------------------------------------------------------- /
!
      TEMP   = ARIN
!
      AROUT(:,NFRMIN:0) = 0.
!
      AROUT(1:NTH,1:NFRCUT) = TEMP
!
      DO ITH=1, NTHMAX
        AROUT(NTH+ITH,1:NFRCUT) = AROUT(   ITH   ,1:NFRCUT)
        AROUT( 1 -ITH,1:NFRCUT) = AROUT(NTH+1-ITH,1:NFRCUT)
        END DO
!
      RETURN
!/
!/ End of EXPND2 ----------------------------------------------------- /
!/
      END SUBROUTINE EXPND2
!/
!/ End of W3SNL3 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SNL3
!/ ------------------------------------------------------------------- /
      SUBROUTINE INSNL3
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Nov-2009 |
!/                  +-----------------------------------+
!/
!/    21-Jul-2008 : Origination as NLX option.          ( version 3.13 )
!/    03-Jan-2009 : Bug fixes NTHMAX and NTHMX2.        ( version 3.13 )
!/    21-Aug-2009 : Conversion to F(f,theta) form.      ( version 3.13 )
!/    13-Nov-2009 : Harden DELTH computation.           ( version 3.13 )
!/
!  1. Purpose :
!
!     Initialization for generalized multiple DIA routine.
!
!  2. Method :
!
!     Fill storage aryays as described in the main subroutine with 
!     interpolation, model and distribution data.
!
!  3. Parameters :
!
!     Variables in W3GDATMD describing model setup.
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      SNLNQ     Int.  Public   Number of quadruplet definitions.
!      SNLL      R.A.  Public   Array with lambda for quadruplet.
!      SNLM      R.A.  Public   Array with mu for quadruplet.
!      SNLT      R.A.  Public   Array with Dtheta for quadruplet.
!      SNLCD     R.A.  Public   Array with Cd for quadruplet.
!      SNLCS     R.A.  Public   Array with Cs for quadruplet.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr. W3SERVMD Program abort.
!      WAVNU2    Subr. W3DISPMD Solve dispersion relation.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3IOGR    Subr. W3IOGRMD Process model definiton file.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     See error escape location.
!
!  8. Remarks :
!
!     - Allocation of arrays directly done in data structure, using
!       IGRID and resetting pointer of aliaases.
!     - In the 03-Jan-2009 bug fix !/T3 error output was fixed, and
!       NTHMAX is increased by 1 to assure that NTHMX2 .LE. NTHMAX
!       for any lambda and mu. With this, the label 810 test is 
!       changed from equality testing to .LE. testing.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/T    General test output.
!     !/T1   Filling of lookup table for quadruplet and interaction
!            strength.
!     !/T2   Filling of lookup table for combining interactions.
!     !/T3   Display raw lookup table of second type.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3ODATMD, ONLY: NDSE, NDST
      USE W3GDATMD,       NFR => NK
!
      USE W3DISPMD, ONLY: WAVNU2
      USE W3SERVMD, ONLY: EXTCDE
!/S      USE W3SERVMD, ONLY: STRACE
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IFRMIN, IFRMAX, IKD, IERR, IQ, NQD,  &
                                 NQS, J, IFR, IQA, JJ, JF, NTHMX2,    &
                                 JIQ, JOF, JQR, IST
      INTEGER                 :: JFR(4), JFR1(4), JTH(4), JTH1(4)
!/S      INTEGER, SAVE           :: IENT = 0
      INTEGER, ALLOCATABLE    :: AST1(:,:,:), AST2(:,:,:)
      REAL                    :: SITMAX, XFRLN
      REAL                    :: OFF12, OFF34, TH12, DEPTH,           &
                                 S0, S1, S2, S3, S4, AUXFR(4),        &
                                 WN0, WN1, WN2, WN3, WN4,             &
                                 CG0, CG1, CG2, CG3, CG4, AUXF,       &
                                 AA, BB, CC, DELTH(4), AUX1, AUX2,    &
                                 WFR(4), WFR1(4), WTH(4), WTH1(4),    &
                                 WFROFF, SIOFF, WF
!
      TYPE QST
        INTEGER               :: OFR(4), OFR1(4), OTH(4), OTH1(4)
        REAL                  :: HFR(4), HFR1(4), HTH(4), HTH1(4)
        REAL                  :: F1, F2, F3, F4, CQD, CQS
      END TYPE QST
!
      TYPE(QST), ALLOCATABLE  :: TSTORE(:,:)
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'INSNL3')
!
! 1.  Initialization ------------------------------------------------- *
! 1.a Checks
!
      XFRLN  = LOG(XFR)
!
      IF ( LAMMAX.LE.0. .OR. LAMMAX.GT.0.5 .OR. DELTHM.LT.0. ) GOTO 800
!
! 1.b Set up relative depths
!
      ALLOCATE ( TSTORE(SNLNQ*4,1:NKD) )
!
      DEPTH  = 1.
      SITMIN = SQRT ( KDMIN * TANH(KDMIN) )
      SITMAX = SQRT ( KDMAX * TANH(KDMAX) )
      XSIT   = (SITMAX/SITMIN)**(1./REAL(NKD-1))
!
!/T      WRITE (NDST,9010) NKD, KDMIN, KDMAX, XSIT
!
! 2.  Building quadruplet data base ---------------------------------- *
!     For quadruplet and interaction strength evaluation
!
      IFRMIN = 0
      IFRMAX = 0
      NTHMAX = 0
!
! 2.a Loop over relative depths
!
      S0     = SITMIN * SQRT ( GRAV / DEPTH )  / XSIT
!
      DO IKD=1, NKD
!
        S0     = S0 * XSIT
        CALL WAVNU2 ( S0, DEPTH, WN0, CG0, 1.E-6, 25, IERR)
!
! 2.b Loop over representative quadruplets
!
        NQA    = 0
        NQD    = 0
        NQS    = 0
!
        DO IQ=1, SNLNQ
!
!/T1          WRITE (NDST,9020) IKD, IQ, WN0*DEPTH, S0*TPIINV, DEPTH
!
          OFF12  = SNLM(IQ)
          OFF34  = SNLL(IQ)
          TH12   = SNLT(IQ) * DERA
          IF ( SNLCD(IQ) .GT. 0. ) NQD = NQD + 1
          IF ( SNLCS(IQ) .GT. 0. ) NQS = NQS + 1
!
          IF ( TH12 .LT. 0. ) THEN
              IF ( OFF12.LT.0. .OR. OFF12.GT.0.5 .OR.                 &
                   OFF34.LT.0. .OR. OFF34.GT.0.5 ) GOTO 801
            ELSE
              IF ( SNLT(IQ).GT.DELTHM .OR. OFF12.LT.0. .OR.           &
                   OFF12.GE.1.                                        &
                   .OR.  OFF34.LT.MINLAM(OFF12,SNLT(IQ)) .OR.         &
                         OFF34.GT.MAXLAM(OFF12,SNLT(IQ)) ) GOTO 802
            END IF
!
!/T1          WRITE (NDST,9021) SNLT(IQ), OFF12, OFF34,               &
!/T1                            SNLCD(IQ), SNLCS(IQ)
!
! 2.c Offset angles
!
          S1     = S0 * ( 1. + OFF12 )
          CALL WAVNU2 ( S1, DEPTH, WN1, CG1, 1.E-6, 25, IERR)
          S2     = S0 * ( 1. - OFF12 )
          CALL WAVNU2 ( S2, DEPTH, WN2, CG2, 1.E-6, 25, IERR)
          S3     = S0 * ( 1. + OFF34 )
          CALL WAVNU2 ( S3, DEPTH, WN3, CG3, 1.E-6, 25, IERR)
          S4     = S0 * ( 1. - OFF34 )
          CALL WAVNU2 ( S4, DEPTH, WN4, CG4, 1.E-6, 25, IERR)
!
          AUXFR(1) = S1 / S0
          AUXFR(2) = S2 / S0
          AUXFR(3) = S3 / S0
          AUXFR(4) = S4 / S0
!
          IF ( TH12 .LT. 0. ) THEN
              BB = 2. * WN0
            ELSE
              BB = WN1**2 + WN2**2 + 2.*WN1*WN2*COS(TH12)
              BB = SQRT ( MAX ( BB , 0. ) )
            END IF
!
          IF ( TH12.LT.0. .AND. ABS(OFF12).LE.1.E-4 ) THEN
              DELTH(1) = 0.
              DELTH(2) = 0.
            ELSE
              CC       = WN1
              AA       = WN2
              AUX1     = (CC**2+BB**2-AA**2) / (2.*BB*CC)
              AUX2     = (AA**2+BB**2-CC**2) / (2.*BB*AA)
              DELTH(1) = - ACOS( MAX ( 0. , MIN ( 1. , AUX1 ) ) )
              DELTH(2) =   ACOS( MAX ( 0. , MIN ( 1. , AUX2 ) ) )
            END IF
          CC       = WN3
          AA       = WN4
          AUX1     = (CC**2+BB**2-AA**2) / (2.*BB*CC)
          AUX2     = (AA**2+BB**2-CC**2) / (2.*BB*AA)
          DELTH(3) = - ACOS( MAX ( 0. , MIN ( 1. , AUX1 ) ) )
          DELTH(4) =   ACOS( MAX ( 0. , MIN ( 1. , AUX2 ) ) )
!
!/T1          WRITE (NDST,9022) DELTH(:) * RADE
!
! 2.d Frequency indices
! 
          DO J=1, 4
            JFR (J) = INT( LOG(AUXFR(J)) / XFRLN )
            JFR1(J) = JFR(J) + 1 * SIGN(1.,AUXFR(J)-1.)
            WFR (J) = (XFR**JFR1(J)-AUXFR(J))/(XFR**JFR1(J)-XFR**JFR(J))
            WFR1(J) = 1. - WFR(J)
            END DO
!
          IFRMIN = MIN ( IFRMIN , MINVAL(JFR1) )
          IFRMAX = MAX ( IFRMAX , MAXVAL(JFR1) )
!
!/T1            WRITE (NDST,9023) 1, JFR(1), JFR1(1), WFR(1), WFR1(1)
!/T1            DO, J=2, 4
!/T1              WRITE (NDST,9024) J, JFR(J), JFR1(J), WFR(J), WFR1(J)
!/T1              END DO
!
! 2.e Directional indices
!
          DO J=1, 4
            AUX1    = DELTH(J) / DTH
            JTH (J) = INT(AUX1)
            JTH1(J) = JTH(J) + 1 * SIGN(1.,DELTH(J))
            WTH1(J) = ABS(AUX1) - REAL(ABS(JTH(J)))
            WTH (J) = 1. - WTH1(J)
            END DO
!
          NTHMAX = MAX ( NTHMAX , MAXVAL(ABS(JTH1)) )
!
!/T1            WRITE (NDST,9025) 1, JTH(1), JTH1(1), WTH(1), WTH1(1)
!/T1            DO, J=2, 4
!/T1              WRITE (NDST,9024) J, JTH(J), JTH1(J), WTH(J), WTH1(J)
!/T1              END DO
!
! 2.f Temp storage of data
!
          IF ( SNLM(IQ).EQ.0. .AND. SNLT(IQ).LT.0. ) THEN
              JJ     = 2
            ELSE
              JJ     = 4
            END IF
!
          DO J=1, JJ
            SELECT CASE (J)
              CASE (2)
                JTH (3) = -JTH (3)
                JTH (4) = -JTH (4)
                JTH1(3) = -JTH1(3)
                JTH1(4) = -JTH1(4)
              CASE (3)
                JTH     = -JTH
                JTH1    = -JTH1
              CASE (4)
                JTH (3) = -JTH (3)
                JTH (4) = -JTH (4)
                JTH1(3) = -JTH1(3)
                JTH1(4) = -JTH1(4)
              CASE DEFAULT
            END SELECT
!
            NQA    = NQA + 1
            TSTORE(NQA,IKD)%OFR  = JFR
            TSTORE(NQA,IKD)%OFR1 = JFR1
            TSTORE(NQA,IKD)%HFR  = WFR
            TSTORE(NQA,IKD)%HFR1 = WFR1
            TSTORE(NQA,IKD)%OTH  = JTH
            TSTORE(NQA,IKD)%OTH1 = JTH1
            TSTORE(NQA,IKD)%HTH  = WTH
            TSTORE(NQA,IKD)%HTH1 = WTH1
            IF ( JJ .EQ. 2 ) THEN
                TSTORE(NQA,IKD)%CQD  = SNLCD(IQ) * 2.
                TSTORE(NQA,IKD)%CQS  = SNLCS(IQ) * 2.
              ELSE
                TSTORE(NQA,IKD)%CQD  = SNLCD(IQ)
                TSTORE(NQA,IKD)%CQS  = SNLCS(IQ)
              END IF
            AUXF                 = ( WN0 * S0 ) / CG0
            TSTORE(NQA,IKD)%F1   = AUXF * CG1 / ( WN1 * S1 )
            TSTORE(NQA,IKD)%F2   = AUXF * CG2 / ( WN2 * S2 )
            TSTORE(NQA,IKD)%F3   = AUXF * CG3 / ( WN3 * S3 )
            TSTORE(NQA,IKD)%F4   = AUXF * CG4 / ( WN4 * S4 )
!
            END DO
!
! ... End loop 2.b
!
          END DO
!
! ... End loop 2.a
!
        END DO
!
!/T1      WRITE (NDST,*)
!/T      WRITE (NDST,9026) NQA, SNLNQ*4, NQD, NQS
!
! 2.g Expanded spectral range
!
      NTHMAX = NTHMAX + 1
!
      NFRMIN =  1  + IFRMIN
      NFRMAX = NFR + IFRMAX - IFRMIN
      NFRCUT = NFR          - IFRMIN
      NTHEXP = NTH + 2*NTHMAX
!
      NSPMIN = 1 + (NFRMIN-1)*NTHEXP - NTHMAX
      NSPMAX = NFRMAX * NTHEXP - NTHMAX
      NSPMX2 = NFRCUT * NTHEXP - NTHMAX
!
!/T      WRITE (NDST,9027) NFR, NFRMIN, NFRMAX, NFRCUT, NTH,          &
!/T             1-NTHMAX, NTH+NTHMAX, NTHEXP
!
      ALLOCATE ( MPARS(IGRID)%SNLPS%FRQ(NFRMAX),                      &
                 MPARS(IGRID)%SNLPS%XSI(NFRMAX) )
      FRQ     => MPARS(IGRID)%SNLPS%FRQ
      XSI     => MPARS(IGRID)%SNLPS%XSI
!
      XSI(1:NFR) = SIG(1:NFR)
      DO IFR=NFR+1, NFRMAX
        XSI(IFR) = XSI(IFR-1) * XFR
        END DO
      FRQ    = XSI * TPIINV
!
! 2.h Final storage 
!
      ALLOCATE ( MPARS(IGRID)%SNLPS%QST1(16,NQA,NKD),                 &
                 MPARS(IGRID)%SNLPS%QST3(6,NQA,NKD),                  &
                 MPARS(IGRID)%SNLPS%QST2(16,NQA,NKD) )
      QST1   => MPARS(IGRID)%SNLPS%QST1
      QST2   => MPARS(IGRID)%SNLPS%QST2
      QST3   => MPARS(IGRID)%SNLPS%QST3
!
! 2.h.1 Basic data
!
      DO IKD=1, NKD
        DO IQA=1, NQA
!
          DO J=1, 4
!
            QST1((J-1)*4+1,IQA,IKD) = TSTORE(IQA,IKD)%OTH (J) +       &
                                      TSTORE(IQA,IKD)%OFR (J) * NTHEXP
            QST1((J-1)*4+2,IQA,IKD) = TSTORE(IQA,IKD)%OTH1(J) +       &
                                      TSTORE(IQA,IKD)%OFR (J) * NTHEXP
            QST1((J-1)*4+3,IQA,IKD) = TSTORE(IQA,IKD)%OTH (J) +       &
                                      TSTORE(IQA,IKD)%OFR1(J) * NTHEXP
            QST1((J-1)*4+4,IQA,IKD) = TSTORE(IQA,IKD)%OTH1(J) +       &
                                      TSTORE(IQA,IKD)%OFR1(J) * NTHEXP
!
            QST2((J-1)*4+1,IQA,IKD) = TSTORE(IQA,IKD)%HFR (J) *       &
                                      TSTORE(IQA,IKD)%HTH (J)
            QST2((J-1)*4+2,IQA,IKD) = TSTORE(IQA,IKD)%HFR (J) *       &
                                      TSTORE(IQA,IKD)%HTH1(J)
            QST2((J-1)*4+3,IQA,IKD) = TSTORE(IQA,IKD)%HFR1(J) *       &
                                      TSTORE(IQA,IKD)%HTH (J)
            QST2((J-1)*4+4,IQA,IKD) = TSTORE(IQA,IKD)%HFR1(J) *       &
                                      TSTORE(IQA,IKD)%HTH1(J)
!
            END DO
!
          QST3(1,IQA,IKD) = TSTORE(IQA,IKD)%F1
          QST3(2,IQA,IKD) = TSTORE(IQA,IKD)%F2
          QST3(3,IQA,IKD) = TSTORE(IQA,IKD)%F3
          QST3(4,IQA,IKD) = TSTORE(IQA,IKD)%F4
          QST3(5,IQA,IKD) = TSTORE(IQA,IKD)%CQD
          QST3(6,IQA,IKD) = TSTORE(IQA,IKD)%CQS
!
          END DO
        END DO
!
      IF ( NQD .GT. 0 ) QST3(5,:,:) = QST3(5,:,:) / REAL(NQD)
      IF ( NQS .GT. 0 ) QST3(6,:,:) = QST3(6,:,:) / REAL(NQS)
!
      DEALLOCATE ( TSTORE )
!
! 3.  Building quadruplet data base ---------------------------------- *
!     For constructing interactions and diagonal from contributions
!
      NTHMX2 = 0
      ALLOCATE ( MPARS(IGRID)%SNLPS%QST4(16,NQA,NKD),                 &
                 MPARS(IGRID)%SNLPS%QST5(16,NQA,NKD),                 &
                 MPARS(IGRID)%SNLPS%QST6(16,NQA,NKD) )
      QST4   => MPARS(IGRID)%SNLPS%QST4
      QST5   => MPARS(IGRID)%SNLPS%QST5
      QST6   => MPARS(IGRID)%SNLPS%QST6
      ALLOCATE ( AST1(16,NQA,NKD), AST2(16,NQA,NKD) )
!
! 3.a Loop over relative depths
!
      S0     = SITMIN * SQRT ( GRAV / DEPTH )  / XSIT
!
      DO IKD=1, NKD
!
        S0     = S0 * XSIT
        CALL WAVNU2 ( S0, DEPTH, WN0, CG0, 1.E-6, 25, IERR)
!
! 3.b Loop over representative quadruplets
!
        NQA    = 0
!
        DO IQ=1, SNLNQ
!
!/T2          WRITE (NDST,9030) IKD, IQ, WN0*DEPTH, S0*TPIINV, DEPTH
!
          OFF12  = SNLM(IQ)
          OFF34  = SNLL(IQ)
          TH12   = SNLT(IQ) * DERA
!
!/T2          WRITE (NDST,9031) SNLT(IQ), OFF12, OFF34
!
! 3.c Frequency indices
!
          AUXFR(1) = ( 1. + OFF12 )
          AUXFR(2) = ( 1. - OFF12 )
          AUXFR(3) = ( 1. + OFF34 )
          AUXFR(4) = ( 1. - OFF34 )
!
          DO J=1, 4
            JFR (J) = INT( LOG(AUXFR(J)) / XFRLN )
            JFR1(J) = JFR(J) + 1 * SIGN(1.,AUXFR(J)-1.)
            WFR (J) = (XFR**JFR1(J)-AUXFR(J))/(XFR**JFR1(J)-XFR**JFR(J))
            WFR1(J) = 1. - WFR(J)
            END DO
!
!/T2            WRITE (NDST,9032) 1, JFR(1), JFR1(1), WFR(1), WFR1(1)
!/T2            DO, J=2, 4
!/T2              WRITE (NDST,9033) J, JFR(J), JFR1(J), WFR(J), WFR1(J)
!/T2              END DO
!
! 3.d Loop over quadruplet components
!
          DO JIQ=1, 4  
!
            IF ( JIQ .LE. 2 ) THEN
                WF     = -1.
              ELSE
                WF     =  1.
              END IF
!
! 3.e Loop over frequency offsets, get directional offsets
!
            DO JOF=1, 2  
!
              IF ( JOF .EQ. 1 ) THEN
                  IFR    = -JFR(JIQ)
                  WFROFF =  WFR(JIQ)
                ELSE
                  IFR    = -JFR1(JIQ)
                  WFROFF =  WFR1(JIQ)
                END IF
!
              SIOFF  = S0 * XFR**IFR
              CALL WAVNU2 ( SIOFF, DEPTH, WN0, CG0, 1.E-6, 25, IERR)
              S1     = SIOFF * ( 1. + OFF12 )
              CALL WAVNU2 ( S1, DEPTH, WN1, CG1, 1.E-6, 25, IERR)
              S2     = SIOFF * ( 1. - OFF12 )
              CALL WAVNU2 ( S2, DEPTH, WN2, CG2, 1.E-6, 25, IERR)
              S3     = SIOFF * ( 1. + OFF34 )
              CALL WAVNU2 ( S3, DEPTH, WN3, CG3, 1.E-6, 25, IERR)
              S4     = SIOFF * ( 1. - OFF34 )
              CALL WAVNU2 ( S4, DEPTH, WN4, CG4, 1.E-6, 25, IERR)
!
!/T2              WRITE (NDST,9034) JIQ, JOF, IFR, WFROFF, SIOFF/S0
!
              IF ( TH12 .LT. 0. ) THEN
                  BB = 2. * WN0
                ELSE
                  BB = WN1**2 + WN2**2 + 2.*WN1*WN2*COS(TH12)
                  BB = SQRT ( MAX ( BB , 0. ) )
                END IF
!
              IF ( TH12.LT.0. .AND. ABS(OFF12).LE.1.E-4 ) THEN
                  DELTH(1) = 0.
                  DELTH(2) = 0.
                ELSE
                  CC       = WN1
                  AA       = WN2
                  AUX1     = (CC**2+BB**2-AA**2) / (2.*BB*CC)
                  AUX2     = (AA**2+BB**2-CC**2) / (2.*BB*AA)
                  DELTH(1) = - ACOS( MAX ( 0. , MIN ( 1. , AUX1 ) ) )
                  DELTH(2) =   ACOS( MAX ( 0. , MIN ( 1. , AUX2 ) ) )
                END IF
              CC       = WN3
              AA       = WN4
              AUX1     = (CC**2+BB**2-AA**2) / (2.*BB*CC)
              AUX2     = (AA**2+BB**2-CC**2) / (2.*BB*AA)
              DELTH(3) = - ACOS( MAX ( 0. , MIN ( 1. , AUX1 ) ) )
              DELTH(4) =   ACOS( MAX ( 0. , MIN ( 1. , AUX2 ) ) )
!
!/T2              WRITE (NDST,9035) DELTH(:) * RADE
!
              AUX1    = DELTH(JIQ) / DTH
              JTH (JIQ) = INT(AUX1)
              JTH1(JIQ) = JTH(JIQ) + 1 * SIGN(1.,DELTH(JIQ))
              WTH1(JIQ) = ABS(AUX1) - REAL(ABS(JTH(JIQ)))
              WTH (JIQ) = 1. - WTH1(JIQ)
!
              NTHMX2 = MAX ( NTHMX2 , ABS(JTH1(JIQ)) )
!
!/T2              WRITE (NDST,9036) JIQ, JTH(JIQ), JTH1(JIQ),         &
!/T2                                     WTH(JIQ), WTH1(JIQ)
!
! 3.f Loop over quadruplet realizations
!
              IF ( SNLM(IQ).EQ.0. .AND. SNLT(IQ).LT.0. ) THEN
                  JJ     = 2
                ELSE
                  JJ     = 4
                END IF
!
              DO JQR=1, JJ
!
                SELECT CASE (JQR)
                  CASE (2)
                    JTH (3) = -JTH (3)
                    JTH (4) = -JTH (4)
                    JTH1(3) = -JTH1(3)
                    JTH1(4) = -JTH1(4)
                  CASE (3)
                    JTH     = -JTH
                    JTH1    = -JTH1
                  CASE (4)
                    JTH (3) = -JTH (3)
                    JTH (4) = -JTH (4)
                    JTH1(3) = -JTH1(3)
                    JTH1(4) = -JTH1(4)
                  CASE DEFAULT
                    JTH     = -JTH
                    JTH1    = -JTH1
                END SELECT
!
                IST    = (JIQ-1)*4 + (JOF-1)*2 + 1
                AST1(IST,NQA+JQR,IKD) = IFR
                AST2(IST,NQA+JQR,IKD) = JTH(JIQ)
                QST5(IST,NQA+JQR,IKD) = WF * ( WFROFF * WTH(JIQ) )
                QST6(IST,NQA+JQR,IKD) = WF * ( WFROFF * WTH(JIQ) )**2
                IST    = IST + 1
                AST1(IST,NQA+JQR,IKD) = IFR
                AST2(IST,NQA+JQR,IKD) = JTH1(JIQ)
                QST5(IST,NQA+JQR,IKD) = WF * ( WFROFF * WTH1(JIQ) )
                QST6(IST,NQA+JQR,IKD) = WF * ( WFROFF * WTH1(JIQ) )**2
!
! ... End loop 3.f
!
                END DO
!
! ... End loop 3.e
!
              END DO
!
! ... End loop 3.d
!
            END DO
!
!/T3          DO JQR=1, JJ
!/T3            WRITE (NDST,9037) IKD, NQA+JQR
!/T3            DO IST=1, 16
!/T3              WRITE (NDST,9038) IST, AST1(IST,NQA+JQR,IKD),       &
!/T3                                     AST2(IST,NQA+JQR,IKD),       &
!/T3                                     QST5(IST,NQA+JQR,IKD),       &
!/T3                                     QST6(IST,NQA+JQR,IKD)
!/T3              END DO
!/T3            END DO
!
! ... End loop 3.b
!
          NQA    = NQA + JJ
!
          END DO
!
! ... End loop 3.a
!
        END DO
!
! 3.g Finalize storage 
!
      QST4 = AST1*NTHEXP + AST2
!
      IF ( NTHMAX .LT. NTHMX2 ) GOTO 810
      IF ( NQA .NE. SIZE(AST1(1,:,1)) ) GOTO 811
!
      DEALLOCATE ( AST1, AST2 )
!
      RETURN
!
! Error escape locations
!
  800 CONTINUE
      WRITE (NDSE,1000) LAMMAX, DELTHM
      CALL EXTCDE ( 1000 )
!
  801 CONTINUE
      WRITE (NDSE,1001) OFF12, OFF34
      CALL EXTCDE ( 1001 )
!
  802 CONTINUE
      WRITE (NDSE,1002) OFF12, OFF34, SNLT(IQ),                       &
                        MINLAM(OFF12,SNLT(IQ)), MAXLAM(OFF12,SNLT(IQ))
      CALL EXTCDE ( 1002 )
!
  810 CONTINUE
      WRITE (NDSE,1010) NTHMAX, NTHMX2
      CALL EXTCDE ( 1010 )
!
  811 CONTINUE
      WRITE (NDSE,1011) NQA, SIZE(AST1(1,:,1))
      CALL EXTCDE ( 1011 )
!
      RETURN
!
! Formats
!
 1000 FORMAT (/' *** WAVEWATCH-III ERROR IN INSNL3 :'/                &
               '     PARAMETER OUT OF RANGE '/                        &
               '     LAMMAX, DELTHM :', 2E12.4/)
 1001 FORMAT (/' *** WAVEWATCH-III ERROR IN INSNL3 :'/                &
               '     PARAMETER OUT OF RANGE '/                        &
               '     MU, LAMBDA :', 2E12.4/)
 1002 FORMAT (/' *** WAVEWATCH-III ERROR IN INSNL3 :'/                &
               '     PARAMETER OUT OF RANGE '/                        &
               '     MU, LAMBDA, TH12 :',3E12.4/                      &
               '     LAMBDA RANGE     :',2E12.4)
 1010 FORMAT (/' *** WAVEWATCH-III ERROR IN INSNL3 :'/                &
               '     NTHMAX LESS THAN NTHMX2 :', 2I8/)
 1011 FORMAT (/' *** WAVEWATCH-III ERROR IN INSNL3 :'/                &
               '     NQA INCONSISTENT :', 2I8/)
!
!/T 9010 FORMAT (/' TEST INSNL3: NKD, KDMIN/MAX/X : ',I8,3F10.4)
!
!/T1 9020 FORMAT (/' TEST INSNL3: IKD, IQ, KD, F, D: ',2I8,2F10.4,F10.2)
!/T1 9021 FORMAT (/' TEST INSNL3: TH12             : ',3X,F8.2/       &
!/T1               '              OFF12, OFF34     : ',3X,2F8.2/      &
!/T1               '              CD, CS           : ',3X,2E10.2)
!/T1 9022 FORMAT ( '              ANGLES (DEGR)    : ',1X,4F8.2)
!/T1 9023 FORMAT ( '              FREQUENCY IND.   : ',1X,3I4,2F6.2)
!/T1 9024 FORMAT ( '                               : ',1X,3I4,2F6.2)
!/T1 9025 FORMAT ( '              DIRECTIONAL IND. : ',1X,3I4,2F6.2)
!/T 9026 FORMAT ( ' TEST INSNL3: FILLING FIRST DATA TABLES :'/        &
!/T               '              NQA AND MAXIMUM  : ',2I8/            &
!/T               '              NQD AND NQS      : ',2I8)
!/T 9027 FORMAT ( '              NFR, MIN/MAX/CUT : ',4I8/            &
!/T               '              NTH, MIN/MAX/EXP : ',4I8)
!
!/T2 9030 FORMAT (/' TEST INSNL3: IKD, IQ, KD, F, D: ',2I8,2F10.4,F10.2)
!/T2 9031 FORMAT (/' TEST INSNL3: TH12             : ',3X,F8.2/       &
!/T2               '              OFF12, OFF34     : ',3X,2F8.2)
!/T2 9032 FORMAT ( '              FREQUENCY IND.   : ',1X,3I4,2F6.2)
!/T2 9033 FORMAT ( '                               : ',1X,3I4,2F6.2)
!/T2 9034 FORMAT ( '              J,J,J, W, SIn    : ',1X,3I4,2F6.2)
!/T2 9035 FORMAT ( '              ANGLES (DEGR)    : ',3X,4F8.2)
!/T2 9036 FORMAT ( '              DIRECTIONAL IND. : ',1X,3I4,2F6.2)
!/T3 9037 FORMAT (/' TEST INSNL3: STORAGE ARRAYS FOR IKD, IQA =',2I6)
!/T3 9038 FORMAT (23X,3I4,3F8.3)
!/
!/ Embedded subroutines
!/ 
      CONTAINS
!/ ------------------------------------------------------------------- /
      REAL FUNCTION MINLAM ( MU, THETA )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Jan-2004 |
!/                  +-----------------------------------+
!/
!/    28-Jan-2009 : Origination.
!/
!  1. Purpose :
!  
!     Calculate minimum allowed lambda for quadruplet configuration.
!     
!  3. Parameters :
!  
!     Parameter list
!     ----------------------------------------------------------------
!       MU, THETA  Real   Quadruplet parameters, theta in degree.
!     ----------------------------------------------------------------
!     
! 10. Source code :
! 
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ Parameter list
!/ 
      REAL, INTENT(IN)        :: MU, THETA
!/
!/ Local parameters
!/ 
      REAL                    :: MULOC, THETAR, BB, AUX
!/
!/ ------------------------------------------------------------------- /
!/ 
      IF ( THETA .LT. 0. ) THEN
          MINLAM = 0. 
        ELSE
          MULOC  = MAX ( 0. , MIN ( 1., MU ) )
          THETAR = THETA * ATAN(1.) / 45.
          BB     = (1.+MULOC)**4 + (1.-MULOC)**4 +                    &
                      2. * (1.+MULOC)**2 * (1.-MULOC)**2 * COS(THETAR)
          BB     = SQRT ( MAX ( BB , 0. ) )
          AUX    = MAX ( 0. , 0.5*BB-1. ) 
          MINLAM = SQRT ( AUX )
        END IF
!
      RETURN
!/
!/ End of MINLAM ----------------------------------------------------- /
!/
      END FUNCTION MINLAM
!/ ------------------------------------------------------------------- /
      REAL FUNCTION MAXLAM ( MU, THETA )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         28-Jan-2004 |
!/                  +-----------------------------------+
!/
!/    28-Jan-2009 : Origination.
!/
!  1. Purpose :
!
!     Calculate minimum allowed lambda for quadruplet configuration.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       MU, THETA  Real   Quadruplet parameters, theta in degree.
!     ----------------------------------------------------------------
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ Parameter list
!/
      REAL, INTENT(IN)        :: MU, THETA
!/
!/ Local parameters
!/
      REAL                    :: MULOC, THETAR, BB, AUX
!/
!/ ------------------------------------------------------------------- /
!/
      IF ( THETA .LT. 0. ) THEN
          MAXLAM = 0.5
        ELSE
          MULOC  = MAX ( 0. , MIN ( 1., MU ) )
          THETAR = THETA * ATAN(1.) / 45.
          BB     = (1.+MULOC)**4 + (1.-MULOC)**4 +                    &
                      2. * (1.+MULOC)**2 * (1.-MULOC)**2 * COS(THETAR)
          BB     = SQRT ( MAX ( BB , 0. ) )
          MAXLAM = 0.25 * BB
        END IF
!
      RETURN
!/
!/ End of MAXLAM ----------------------------------------------------- /
!/
      END FUNCTION MAXLAM
!/
!/ End of INSNL3 ----------------------------------------------------- /
!/
      END SUBROUTINE INSNL3
!/
!/ End of module W3SNL3MD -------------------------------------------- /
!/
      END MODULE W3SNL3MD
