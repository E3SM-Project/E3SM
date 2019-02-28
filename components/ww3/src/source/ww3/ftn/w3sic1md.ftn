#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SIC1MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           E. Rogers               |
!/                  |           S. Zieger               |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         11-Oct-2013 |
!/                  +-----------------------------------+
!/
!/    For updates see W3SIC1 documentation.
!/
!  1. Purpose :
!
!     Calculate ice source term S_{ice} according to simple methods.
!          Exponential decay rate is uniform in frequency, and 
!          specified directly by the user.  This method is, in effect, 
!          not sustantially different from handling sea ice via the 
!          "sub-grid" blocking approach, after improvements by 
!          Fabrice Ardhuin (in v4.00).
!
!  2. Variables and types :
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SIC1    Subr. Public   ice source term.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!     See subroutine documentation.
!
!  5. Remarks :
!
!     Reference:Rogers, W.E. and M.D. Orzech, 2013: Implementation and
!        Testing of Ice and Mud Source Functions in WAVEWATCH III(R), 
!        NRL/MR/7320--13-9462, 31pp.
!        available from http://www7320.nrlssc.navy.mil/pubs.php
!        Direct link: 
!        http://www7320.nrlssc.navy.mil/pubs/2013/rogers2-2013.pdf
!
!  6. Switches :
!
!     See subroutine documentation.
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      PUBLIC :: W3SIC1
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SIC1 (A, DEPTH, CG, IX, IY, S, D)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           E. Rogers               |
!/                  |           S. Zieger               |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         11-Oct-2013 |
!/                  +-----------------------------------+
!/
!/    16-Oct-2012 : Origination.                        ( version 4.04 )
!/                                                        (E. Rogers)
!/    09-Oct-2013 : W3SIC1 SUBTYPE=2 outsourced to W3SIC2 (S. Zieger)
!/
!/        FIXME   : Move field input to W3SRCE and provide
!/     (S.Zieger)   input parameter to W3SIC1 to make the subroutine
!/                : versatile for point output processors ww3_outp
!/                  and ww3_ounp.
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS. 
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     S_{ice} source term using 5 parameters read from input files.
!     These parameters are allowed to vary in space and time.
!     The parameters control the exponential decay rate k_i
!     Since there are 5 parameters, this permits description of
!     dependence of k_i on frequency or wavenumber.
!
!/ ------------------------------------------------------------------- /
!
!  2. Method :
!
!     Regarding i/o (general to all Sice modules): S_{ice} source term
!     is calculated using up to 5 parameters read from input files.
!     These parameters are allowed to vary in space and time.
!     The parameters control the exponential decay rate k_i
!     Since there are 5 parameters, this permits description of
!     dependence of k_i on frequency or wavenumber.
!
!     Sea ice affects the wavenumber k of wind-generated ocean waves.
!     The ice-modified wavenumber can be expressed as a complex number
!     k = k_r + i*k_i, with the real part k_r representing impact of
!     the sea ice on the physical wavelength and propagation speeds, 
!     producing something analogous to shoaling and refraction by 
!     bathymetry, whereas the imaginary part of the complex 
!     wavenumber, k_i, is an exponential decay coefficient 
!     k_i(x,y,t,sigma) (depending on location, time and frequency, 
!     respectively), representing wave attenuation, and can be 
!     introduced in a wave model such as WW3 as S_ice/E=-2*Cg*k_i, 
!     where S_ice is one of several dissipation mechanisms, along 
!     with whitecapping, for example, S_ds=S_wc+S_ice+⋯. The k_r - 
!     modified by ice would enter the model via the C calculations 
!     on the left-hand side of the governing equation.The fundamentals
!     are straightforward, e.g. Rogers and Holland (2009 and 
!     subsequent unpublished work) modified a similar model, SWAN 
!     (Booij et al. 1999) to include the effects of a viscous mud 
!     layer using the same approach (k = k_r + i*k_i) previously.
!
!     General approach is analogous to Rogers and Holland (2009) 
!         approach for mud.
!     See text near their eq. 1 :
!       k        = k_r  +  i * k_i
!       eta(x,t) = Real( a * exp( i * ( k * x - sigma * t ) ) )
!       a        = a0 * exp( -k_i * x )
!       S / E    = -2 * Cg * k_i (see also Komen et al. (1994, pg. 170)
!
!     Following W3SBT1 as a guide, equation 1 of W3SBT1 says:
!         S = D * E
!     However, the code of W3SBT1 has
!         S = D * A
!     This leads me to believe that the calling routine is 
!         expecting "S/sigma" not "S"
!     Thus we will use D = S/E = -2 * Cg * k_i
!
!     Notes regarding numerics:
!     Experiments with constant k_i values suggest that :
!       for dx=20.0 km, k_i should not exceed 3.5e-6
!      (assumes 2.7% Hs error in my particular test case is intolerable)
!       for dx=5.0 km,  k_i should not exceed 2.0e-5
!       for dx=2.5 km,  k_i should not exceed 5.0e-5
!       for dx=1.0 km,  k_i should not exceed 2.0e-4
!       for dx=0.35 km, error is less than 2.1% for all k_i tested
!       for dx=0.10 km, error is less than 1.3% for all k_i tested
!     "Ground truth" used for this is an exponential decay profile.
!
!      For reference, ACNFS is 1/12th deg, so delta_latitude=9.25 km.
!
!     {put more equations here}
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D)
!       DEPTH   Real  I   Local water depth
!       CG      R.A.  I   Group velocities.
!       IX,IY   I.S.  I   Grid indices.
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
!      W3EXPO    Subr.   N/A    ASCII Point output post-processor.
!      W3EXNC    Subr.   N/A    NetCDF Point output post-processor.
!      GXEXPO    Subr.   N/A    GrADS point output post-processor.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     None.
!
!  7. Remarks :
!
!     If ice parameter 1 is zero, no calculations are made.
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
      USE CONSTANTS, ONLY: TPI
      USE W3ODATMD, ONLY: NDSE
      USE W3SERVMD, ONLY: EXTCDE
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, MAPWN
      USE W3IDATMD, ONLY: ICEP1, ICEP2, ICEP3, ICEP4, ICEP5, INFLAGS2
!/T      USE W3ODATMD, ONLY: NDST
!/S      USE W3SERVMD, ONLY: STRACE
!/T0      USE W3ARRYMD, ONLY: PRT2DS
!/T1      USE W3ARRYMD, ONLY: OUTMAT
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
      REAL, INTENT(IN)        :: CG(NK),   A(NSPEC), DEPTH
      REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
      INTEGER, INTENT(IN)     :: IX, IY
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/S      INTEGER, SAVE           :: IENT = 0
!/T0      INTEGER                 :: ITH
!/T0      REAL                    :: DOUT(NK,NTH)
      INTEGER                 :: IKTH, IK
      REAL                    :: D1D(NK) !In SBT1: D1D was named "CBETA"
      REAL                    :: ICECOEF1, ICECOEF2, ICECOEF3, &
                                 ICECOEF4, ICECOEF5
      REAL, ALLOCATABLE       :: WN_I(:)  ! exponential decay rate for amplitude
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SIC1')
!
! 0.  Initializations ------------------------------------------------ *
!
      D        = 0.0
!
      ALLOCATE(WN_I(NK))
      WN_I     = 0.0
      ICECOEF1 = 0.0
      ICECOEF2 = 0.0
      ICECOEF3 = 0.0
      ICECOEF4 = 0.0
      ICECOEF5 = 0.0
!
      IF (.NOT.INFLAGS2(-7))THEN
         WRITE (NDSE,1001) 'ICE PARAMETER 1'
         CALL EXTCDE(2)
      ENDIF
!
      ICECOEF1 = ICEP1(IX,IY)

!
! 1.  No ice --------------------------------------------------------- /
!
      IF ( ICECOEF1==0. ) THEN
         D = 0.
!
! 2.  Ice ------------------------------------------------------------ /
      ELSE
!
! 2.a Set constant(s) and write test output -------------------------- /
!
!         (none)
!
!/T38        WRITE (NDST,9000) DEPTH,ICECOEF1,ICECOEF2,ICECOEF3,ICECOEF4
!
! 2.b Make calculations ---------------------------------------------- /
         WN_I = ICECOEF1 ! uniform in k

         DO IK=1, NK
!   SBT1 has: D1D(IK) = FACTOR *  MAX(0., (CG(IK)*WN(IK)/SIG(IK)-0.5) )
!             recall that D=S/E=-2*Cg*k_i
            D1D(IK) = -2. * CG(IK) * WN_I(IK)
         END DO
!
! 2.c Fill diagional matrix
!
         DO IKTH=1, NSPEC
            D(IKTH) = D1D(MAPWN(IKTH))
         END DO
!
      END IF
!
      S = D * A
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
!/T0                         0.0, 0.001, 'Diag Sice', ' ', 'NONAME')
!
!/T1      CALL OUTMAT (NDST, D, NTH, NTH, NK, 'diag Sice')
!
! Formats
!
 1001 FORMAT (/' *** WAVEWATCH III ERROR IN W3SIC1 : '/               &
               '     ',A,' REQUIRED BUT NOT SELECTED'/)
!
!/T 9000 FORMAT (' TEST W3SIC1 : DEPTH,ICECOEF1  : ',2E10.3)
!/
!/ End of W3SIC1 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SIC1
!/
!/ End of module W3SIC1MD -------------------------------------------- /
!/
      END MODULE W3SIC1MD
