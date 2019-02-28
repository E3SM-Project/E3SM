#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SNLSMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Jul-2012 |
!/                  +-----------------------------------+
!/
!/    04-Aug-2008 : Origination in research model.      ( version 3.13 )
!/    27-Sep-2010 : Added to svn repository.            ( version 3.15 )
!/    13-Jul-2012 : Moved from version 3.15 to 4.08.    ( version 4.08 )
!/
!/    Copyright 2009-2012 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Nonlinear interaction based `smoother' for high frequencies.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NKD       I.P.  Private  Number of nondimensional depths in
!                               storage array.
!      KDMIN     R.P.  Private  Minimum relative depth in table.
!      KDMAX     R.P.  Private  Maximum relative depth in table.
!      SITMIN    Real  Private  Minimum nondimensional radian
!                               frequency in table.
!      XSIT      Real  Private  Corresponding incremet factor.
!      ABMAX     R.P.  Public   Maximum value of a34, b3 and b4.
!     ----------------------------------------------------------------
!
!     Variables in W3GDATMD :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      CNLSA     Real  Public   a34 in quadruplet definition.
!      CNLSC     Real  Public   C in Snl definition.
!      CNLSFM    Real  Public   Maximum relative spectral change.
!      CNLSC1/3  Real  Public   Constant in frequency filter.
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SNLS    Subr. Public   Nonlinear 'smoother' algorithm.
!      EXPAND    Subr. W3SNLS   Expand spectrum for indirect address.
!      INSNLS    Subr. Public   Initialization routine.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WAVNU1    Subr. W3DISPMD Solve dispersion relation.
!      WAVNU2    Subr. W3DISPMD Solve dispersion relation.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr. W3SERVMD Program abort.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output.
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER, PRIVATE, PARAMETER :: NKD = 100
      REAL, PRIVATE, PARAMETER    :: KDMIN = 0.25 ,  KDMAX = 10.
      REAL, PRIVATE               :: SITMIN, XSIT
!
      REAL, PARAMETER             :: ABMAX = 0.25
!
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SNLS (  A, CG, WN, DEPTH, UABS, DT, SNL, AA )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Aug-2008 |
!/                  +-----------------------------------+
!/
!/    04-Aug-2008 : Origination.                        ( version 3.13 )
!/
!  1. Purpose :
!
!     High-frequeny filter based on the nonlinear interactions for 
!     an uresolved quadruplet.
!
!  2. Method :
!
!     Compute interactions for a quadruplet that is not resolved by
!     the discrete spectral rsolution, and then reduces to a simple
!     five-point stencil. Furthermore interactions are filtered by
!     frequency to allow for high-frequency impact only, and the
!     integration schem is embedded, and reduces to a filter technique
!     for large time steps or strong interactions.
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
!       UABS    Real  I   Wind speed (m/s).
!       DT      Real  I   Numerical time step (s).
!       SNL     R.A.  O   Nonlinear source term.                (Opt)
!       AA      R.A.  O   Averaged spectrum.                    (Opt)
!     ----------------------------------------------------------------
!       Note: A and AA may safely be same array/address.
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WAVNU1    Subr. W3DISPMD Solve dispersion relation.
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
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output.
!     !/T1   Test output frequency filter.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3GDATMD, ONLY: NFR => NK, NTH, SIG, XFR, FACHFA, DTH,      &
                          NTHX, NFRX, NSPL, NSPH, SNSST, CNLSC,       &
                          CNLSFM, CNLSC1, CNLSC2, CNLSC3
      USE W3ODATMD, ONLY: NDST, NDSE
!
      USE W3DISPMD, ONLY: WAVNU1
!/S      USE W3SERVMD, ONLY: STRACE
!/T2      USE W3ARRYMD, ONLY: PRT2DS
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)            :: A(NTH,NFR), CG(NFR), WN(NFR),    &
                                     DEPTH, UABS, DT
      REAL, INTENT(OUT), OPTIONAL :: SNL(NTH,NFR), AA(NTH,NFR)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IFR, IFRMIN, ITH, IFRMN2,            &
                                 IKD, JKD(0:NFR+2), ISPX0, ISPX
!/S      INTEGER, SAVE           :: IENT = 0
      REAL                    :: SIGP, CP, CM, XL, XH, EL, EH, DENOM, &
                                 SIT, XSITLN, MC, F3A,  F3B, F3C,     &
                                 F4A, F4B, F4C, F00, F31, F32, F41,   &
                                 F42, AUXB, AUX11, AUX21, AUX12,      &
                                 AUX22, FC1, FC2, FC3, FC4
      REAL                    :: XSI(NFR+2), XWN(NFR+2), XCG(NFR+2),  &
                                 UP(NSPL:NSPH), UN(NSPL:NSPH),        &
                                 E1(0:NFR+2), FILTFP(NFR+2),          &
                                 FPROP(NFR+2), DS1(NSPL:NSPH),        &
                                 DS2(NSPL:NSPH), DS3(NSPL:NSPH),      &
                                 DA1(NSPL:NSPH), DA2(NSPL:NSPH),      &
                                 DA3(NSPL:NSPH)
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'W3SNLS')
!
!/T      WRITE (NDST,9000) DEPTH, UABS, DT
!
! 1.  Initializations ------------------------------------------------ *
! 1.a Expanded frequency range
!
      XSI(1:NFR) = SIG(1:NFR)
      XWN(1:NFR) = WN
      XCG(1:NFR) = CG
!
      XSI(NFR+1) = XSI(NFR) * XFR
      CALL WAVNU1 ( XSI(NFR+1), DEPTH, XWN(NFR+1), XCG(NFR+1) )
      XSI(NFR+2) = XSI(NFR+1) * XFR
      CALL WAVNU1 ( XSI(NFR+2), DEPTH, XWN(NFR+2), XCG(NFR+2) )
!
! 1.b Expanded psuedo spectrum
!
      CALL EXPAND ( UP, UN )
!
! 1.c Get relevant spectral peak frequency
!
!/T1      E1     = -1.
      SIGP   = - TPI
      XL     = 1./XFR - 1.
      XH     =  XFR - 1.
!
! 1.c.1 Wind too weak
!
      IF ( UABS .LT. XSI(NFR)/XWN(NFR) ) THEN
          SIGP   = GRAV / MAX ( 0.01 , UABS )
        ELSE
!
! 1.c.2 Compute 1D spectrum
!
          E1(NFR+2) = SUM(A(:,NFR)) * FACHFA**2 * XSI(NFR+2)          &
                                           / XCG(NFR+2) * TPI * DTH 
          E1(NFR+1) = SUM(A(:,NFR)) * FACHFA    * XSI(NFR+1)          &
                                           / XCG(NFR+1) * TPI * DTH 
!
          DO IFR=NFR, 1, -1
            E1(IFR) = SUM(A(:,IFR)) * XSI(IFR) / XCG(IFR) * TPI * DTH 
!
! 1.c.3 Reached PM frequency
!
            IF ( UABS .LT. XSI(IFR)/XWN(IFR) ) THEN
                CP     = XSI(IFR)/XWN(IFR)
                CM     = XSI(IFR+1)/XWN(IFR+1)
                SIGP   = XSI( IFR ) * (UABS-CM)/(CP-CM) +             &
                         XSI(IFR+1) * (CP-UABS)/(CP-CM)
                EXIT
!
              ELSE IF ( E1(IFR) .LT. E1(IFR+1) ) THEN
!
! 1.c.4 Reached first peak
!
                EL     = E1(IFR  ) - E1(IFR+1)
                EH     = E1(IFR+2) - E1(IFR+1)
                DENOM  = XL*EH - XH*EL
                SIGP   = XSI(IFR+1) * (1.+0.5*(XL**2*EH-XH**2*EL)     &
                         / SIGN ( MAX(ABS(DENOM),1.E-15) , DENOM ) )
                EXIT
              ENDIF
!
! ... End loop 1.c.2
!
            END DO
!
! 1.c.5 Nothing found
!
          IF ( SIGP .LT. 0. ) THEN
!
! 1.c.5.a No energy there
!
              IF ( E1(1) .EQ. 0. ) THEN
                  SIGP   = 2. * SIG(NFR)
!
! 1.c.5.b Peak at low boundary
!
                ELSE
                  SIGP   = XSI(1)
                END IF
            END IF
!
        END IF
!
! 1.d Set up filter function etc.
!
      XSITLN = LOG(XSIT)
      IFRMIN =  1
      JKD    =  1
!/T1      FILTFP = -1.
!
      DO IFR=NFR+2, 1, -1
!
        FILTFP(IFR) = EXP(-CNLSC1/(XSI(IFR)/(CNLSC2*SIGP))**CNLSC3)
        FPROP (IFR) = FILTFP(IFR) * CNLSC * XWN(IFR)**8 *             &
                         XSI(IFR)**4 / TPI**9 / XCG(IFR)
        SIT      = XSI(IFR) * SQRT(DEPTH/GRAV)
        IKD      = 1 + NINT ( ( LOG(SIT) - LOG(SITMIN) ) / XSITLN )
        JKD(IFR) = MAX ( 1 , MIN(IKD,NKD) )
!
        IF ( FILTFP(IFR) .LT. 1.E-10 ) THEN
            IFRMIN = IFR
            EXIT
          END IF
!
        END DO
!
      IFRMN2 = MAX ( 1 , IFRMIN - 1 )
      SIT    = XSI(IFRMN2) * SQRT(DEPTH/GRAV)
      IKD    = 1 + NINT ( ( LOG(SIT) - LOG(SITMIN) ) / XSITLN )
      JKD(IFRMN2) = MAX ( 1 , MIN(IKD,NKD) )
!
!/T      WRITE (NDST,9010) IFRMIN, SIGP * TPIINV
!/T1      WRITE (NDST,9011)
!/T1      DO IFR=1, NFR
!/T1        WRITE (NDST,9012) IFR, XSI(IFR)/TPI, XSI(IFR)/XWN(IFR),   &
!/T1                          E1(IFR), FILTFP(IFR)
!/T1        END DO
!
! 1.e Initialize arrays
!
!
! 2.  Compute base interactions -------------------------------------- *
! 2.a Loop over frequencies
!
      DO IFR=IFRMIN, NFR+1
!
        ISPX0  = (IFR-1)*NTHX
        IKD    = JKD(IFR)
!
        MC     = SNSST( 1,IKD)
        F3A    = SNSST( 2,IKD)
        F3B    = SNSST( 3,IKD)
        F3C    = SNSST( 4,IKD)
        F4A    = SNSST( 5,IKD)
        F4B    = SNSST( 6,IKD)
        F4C    = F3C
!
! 2.b Loop over directions
!
        DO ITH=1, NTH
!
          ISPX   = ISPX0 + ITH
!
          F00    = UP(ISPX)
          F31    = UP(ISPX)*F3A + UP(ISPX+1)*F3B + UP(ISPX+NTHX)*F3C
          F41    = UP(ISPX)*F4A + UP(ISPX-1)*F4B + UP(ISPX-NTHX)*F4C
          F32    = UP(ISPX)*F3A + UP(ISPX-1)*F3B + UP(ISPX+NTHX)*F3C
          F42    = UP(ISPX)*F4A + UP(ISPX+1)*F4B + UP(ISPX-NTHX)*F4C
!
          DS1(ISPX) = FPROP(IFR) * (F00**2*(F31+F41)-2.*F00*F31*F41)
          DS2(ISPX) = FPROP(IFR) * (F00**2*(F32+F42)-2.*F00*F32*F42)
!
          AUX11  = DT * DS1(ISPX)
          AUX21  = DT * DS2(ISPX)
          AUXB   = CNLSFM * FILTFP(IFR) * MAX(1.E-10,UN(ISPX)) /   &
                      MAX ( 1.E-10 , ABS(AUX11)+ABS(AUX21) ) / MC 
          AUX12  = AUXB * ABS(AUX11)
          AUX22  = AUXB * ABS(AUX21)
!
! Expensive but more smooth limiter
!
!         DA1(ISPX) = AUX12 * TANH(AUX11/MAX(1.E-10,AUX12))
!         DA2(ISPX) = AUX22 * TANH(AUX21/MAX(1.E-10,AUX22))
!
! Crude but cheaper limiter
!
          DA1(ISPX) = MAX ( -AUX12 , MIN ( AUX11 , AUX12 ) )
          DA2(ISPX) = MAX ( -AUX22 , MIN ( AUX21 , AUX22 ) )
!
          END DO
!
! ... End loop 2.b
!
        END DO
!
! 2.c Complete expanded arrays
!
! ... End loop 2.a
!
! 3.  Compute source term if requested ------------------------------- *
! 3.a Check for request
!
      IF ( PRESENT(SNL) ) THEN
!/T          WRITE (NDST,9030) 'YES/--'
!
! 3.b Initializations
!
          SNL(:,1:IFRMN2-1) = 0.
!
          DS1(NSPL:IFRMN2*NTHX-1) = 0.
          DS2(NSPL:IFRMN2*NTHX-1) = 0.
          DS3(NSPL:IFRMN2*NTHX-1) = 0.
!
          ISPX  = IFRMN2*NTHX
          DS1(ISPX+NTH+1:NSPH:NTHX) = DS1(ISPX+ 1 :NSPH:NTHX)
          DS1(ISPX      :NSPH:NTHX) = DS1(ISPX+NTH:NSPH:NTHX)
          DS2(ISPX+NTH+1:NSPH:NTHX) = DS2(ISPX+ 1 :NSPH:NTHX)
          DS2(ISPX      :NSPH:NTHX) = DS2(ISPX+NTH:NSPH:NTHX)
          DS3(IFRMN2*NTHX:NSPH)     = DS1(IFRMN2*NTHX:NSPH)  +        &
                                      DS2(IFRMN2*NTHX:NSPH)
!
! 3.c Loop over frequencies
!
          DO IFR=IFRMN2, NFR
!
            ISPX0  = (IFR-1)*NTHX
            IKD    = JKD(IFR)
!
            FC1    = - SNSST(1,IKD)
            FC2    =   SNSST(4,IKD)
            FC3    =   SNSST(3,IKD)
            FC4    =   SNSST(6,IKD)
!
! 3.d Loop over directions
!
            DO ITH=1, NTH
              ISPX         = ISPX0 + ITH
              SNL(ITH,IFR) = FC1 *             DS3(   ISPX  )         &
                 +    FC2 * ( DS3(ISPX-NTHX) + DS3(ISPX+NTHX) )       &
                 +    FC3 * ( DS1(ISPX-  1 ) + DS2(ISPX+  1 ) )       &
                 +    FC4 * ( DS1(ISPX+  1 ) + DS2(ISPX-  1 ) )
!
              END DO
!
! ... End loop 3.d
!
            END DO
!
! ... End loop 3.c
!
!/T        ELSE
!/T          WRITE (NDST,9030) '---/NO'
        END IF
!
! 4.  Compute filtered spectrum if requested ------------------------- *
! 4.a Check for request
!
      IF ( PRESENT(AA) ) THEN
!/T          WRITE (NDST,9040) 'YES/--'
!
! 4.b Initializations
!
          AA(:,1:IFRMN2-1) = A(:,1:IFRMN2-1)
!
          DA1(NSPL:IFRMN2*NTHX-1) = 0.
          DA2(NSPL:IFRMN2*NTHX-1) = 0.
          DA3(NSPL:IFRMN2*NTHX-1) = 0.
!
          ISPX  = IFRMN2*NTHX
          DA1(ISPX+NTH+1:NSPH:NTHX) = DA1(ISPX+ 1 :NSPH:NTHX)
          DA1(ISPX      :NSPH:NTHX) = DA1(ISPX+NTH:NSPH:NTHX)
          DA2(ISPX+NTH+1:NSPH:NTHX) = DA2(ISPX+ 1 :NSPH:NTHX)
          DA2(ISPX      :NSPH:NTHX) = DA2(ISPX+NTH:NSPH:NTHX)
          DA3(IFRMN2*NTHX:NSPH)     = DA1(IFRMN2*NTHX:NSPH)  +        &
                                      DA2(IFRMN2*NTHX:NSPH)
!
! 4.c Loop over frequencies
!
          DO IFR=IFRMN2, NFR
!
            ISPX0  = (IFR-1)*NTHX
            IKD    = JKD(IFR)
!
            FC1    = - SNSST(1,IKD)
            FC2    =   SNSST(4,IKD)
            FC3    =   SNSST(3,IKD)
            FC4    =   SNSST(6,IKD)
!
! 4.d Loop over directions
!
            DO ITH=1, NTH
              ISPX         = ISPX0 + ITH
              AA(ITH,IFR) = MAX ( 0. , A(ITH,IFR) +                   &
                   FC1 *   DA3(ISPX)                                  &
                 + FC2 * ( DA3(ISPX-NTHX) + DA3(ISPX+NTHX) )          &
                 + FC3 * ( DA1(ISPX-  1 ) + DA2(ISPX+  1 ) )          &
                 + FC4 * ( DA1(ISPX+  1 ) + DA2(ISPX-  1 ) ) )
              END DO
!
! ... End loop 4.d
!
            END DO
!
! ... End loop 4.c
!
!/T        ELSE
!/T          WRITE (NDST,9040) '---/NO'
        END IF
!
!/T stop
      RETURN
!
! Formats
!
!/T 9000 FORMAT (/' TEST W3SNLS: DEPTH, UABS, DT :',F9.2,F7.2,F7.2)
!
!/T 9010 FORMAT ( '              IFRMIN, FP  :',I4,F8.4)
!/T1 9011 FORMAT ( ' TEST W3SNLS: IFR, FR, C, E1, FILT :')
!/T1 9012 FORMAT (13X,I4,F10.4,2F10.2,F10.4)
!
!/T 9030 FORMAT ( ' TEST W3SNLS: SOURCE TERM REQUESTED : ',A)
!/T 9040 FORMAT ( ' TEST W3SNLS: AVERAGING REQUESTED   : ',A)
!/
!/ Embedded subroutines
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE EXPAND ( PSPC, SPEC )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         23-Jul-2008 |
!/                  +-----------------------------------+
!/
!  1. Purpose :
!
!     Expand spectrum to simplify indirect addressing.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       PSPC    R.A.  O   Expanded spectrum.
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
      REAL, INTENT(OUT)       :: PSPC(0:NTH+1,0:NFR+2),               &
                                 SPEC(0:NTH+1,0:NFR+2)
!/
!/ Local parameters
!/ 
      INTEGER                 :: IFR, ITH
!/
!/ ------------------------------------------------------------------- /
!
      SPEC(:,0) = 0. 
!
      SPEC(1:NTH,1:NFR) = A
      SPEC(1:NTH,NFR+1) = SPEC(1:NTH,NFR) * FACHFA
      SPEC(1:NTH,NFR+2) = SPEC(1:NTH,NFR+1) * FACHFA
!
      SPEC(NTH+1,1:NFR+2) = SPEC( 1 ,1:NFR+2)
      SPEC(  0  ,1:NFR+2) = SPEC(NTH,1:NFR+2)
!
      DO IFR=1, NFR+2
        PSPC(:,IFR) = SPEC(:,IFR) / XWN(IFR)
        END DO 
!
      RETURN
!/
!/ End of EXPAND ----------------------------------------------------- /
!/
      END SUBROUTINE EXPAND
!/
!/ End of W3SNLS ----------------------------------------------------- /
!/
      END SUBROUTINE W3SNLS
!/ ------------------------------------------------------------------- /
      SUBROUTINE INSNLS
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-Aug-2008 |
!/                  +-----------------------------------+
!/
!/    04-Aug-2008 : Origination.                        ( version 3.13 )
!/
!  1. Purpose :
!
!     Initializations for the Snl / filter source term for high 
!     frequencies.
!
!  2. Method :
!
!     Precompute weight functions and store in array.
!
!  3. Parameters :
!
!     No parameter list.
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      WAVNU2    Subr. W3DISPMD Solve dispersion relation.
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr. W3SERVMD Program abort.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3IOGR    Subr. W3IOGRMD Process model definition file.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - Check a34, b4 and b5 against MAXAB to assure that the values
!       are consistent with a reduced 5-point stencil for unresolved
!       quadruplets. a34 is checked in ww3_grid, b3 and b4 are not.
!
!  7. Remarks :
!
!     - Small quadruplet compared to grid size reduces interactions
!       so that distribution of results is purely local. This results
!       in a much simpler model initialization than for the general
!       MDIA.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3ODATMD, ONLY: NDST, NDSE
      USE W3GDATMD,       NFR => NK, A34 => CNLSA
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
      INTEGER                 :: IKD, IERR
!/S      INTEGER, SAVE           :: IENT = 0
      REAL                    :: DEPTH, SITMAX, OFF, S0, WN0, CG0,    &
                                 S3, WN3, CG3, S4, WN4, CG4, WN12,    &
                                 DT3, DT4, B3, B4
!/
!/ ------------------------------------------------------------------- /
!/
!/S      CALL STRACE (IENT, 'INSNLS')
!
! 1.  Initializations ------------------------------------------------ *
! 1.a Set up relative depths
!
      DEPTH  = 1.
      SITMIN = SQRT ( KDMIN * TANH(KDMIN) )
      SITMAX = SQRT ( KDMAX * TANH(KDMAX) )
      XSIT   = (SITMAX/SITMIN)**(1./REAL(NKD-1))
!
!/T      WRITE (NDST,9010) NKD, KDMIN, KDMAX, XSIT
!
! 1.b Set up quadruplet
!
      OFF    = (XFR-1.) * A34
!
! 1.c Set up storage
!
      NTHX   = NTH + 2
      NFRX   = NFR + 2
      NSPL   = - NTHX
      NSPH   = NFRX*NTHX - 1
!
      ALLOCATE ( MPARS(IGRID)%SNLPS%SNSST(6,NKD) )
      SNSST  => MPARS(IGRID)%SNLPS%SNSST
!
! 2.  Building quadruplet data base ---------------------------------- *
!     For quadruplet and interaction strength evaluation
!
      S0     = SITMIN * SQRT ( GRAV / DEPTH )  / XSIT
!
! 2.a Loop over relative depths
!
      DO IKD=1, NKD
!
! 2.b Base quadruplet set up
!
        S0     = S0 * XSIT
        S3     = ( 1. + OFF ) * S0
        S4     = ( 1. - OFF ) * S0
!
        CALL WAVNU2 ( S0, DEPTH, WN0, CG0, 1.E-6, 25, IERR)
        CALL WAVNU2 ( S3, DEPTH, WN3, CG3, 1.E-6, 25, IERR)
        CALL WAVNU2 ( S4, DEPTH, WN4, CG4, 1.E-6, 25, IERR)
!
!/T        WRITE (NDST,9020) IKD, WN0*DEPTH, S0*TPIINV, DEPTH
!
! 2.c Offset angles
!
        WN12   = 2. * WN0
        DT3    = ACOS( (WN3**2+WN12**2-WN4**2) / (2.*WN12*WN3) )
        DT4    = ACOS( (WN4**2+WN12**2-WN3**2) / (2.*WN12*WN4) )
!
        B3     = DT3 / DTH
        B4     = DT4 / DTH
!
!/T        WRITE (NDST,9021) A34, B3, B4, DT3*RADE, DT4*RADE
!
        IF ( A34.GT.ABMAX .OR. B3.GT.ABMAX .OR. B4.GT.ABMAX .OR.     &
             A34.LT.0. .OR. B3.LT.0. .OR. B4.LT.0. ) GOTO 801
!
! 2.d Store weights
!
        SNSST( 1,IKD) = 2.*A34 + B3 + B4
        SNSST( 2,IKD) = 1. - A34 - B3
        SNSST( 3,IKD) = B3
        SNSST( 4,IKD) = A34
        SNSST( 5,IKD) = 1. - A34 - B4
        SNSST( 6,IKD) = B4
!
! ... End loop 2.a
!
        END DO
!
      RETURN
!
! Error escape locations
!
  801 CONTINUE
      WRITE (NDSE,1001) A34, B3, B4
      CALL EXTCDE (1001)
!
! Formats
!
 1001 FORMAT (/' *** WAVEWATCH-III ERROR IN INSNLS :'/                &
               '     PARAMETER FORCED OUT OF RANGE '/                 &
               '     A34, B3, B4 :', 3F10.4/)
!
!/T 9010 FORMAT (/' TEST INSNLS: NKD, KDMIN/MAX/X :',I5,3F10.4)
!/T 9020 FORMAT ( '              IKD, KD, F, D    :',I5,3F10.4)
!/T 9021 FORMAT ( '              A34, B3,B4, TH3/4:',3F7.3,2F6.2)
!/
! /End of INSNLS ------------------------------------------------------/
!/
      END SUBROUTINE INSNLS
!/
!/ End of module W3SNLSMD -------------------------------------------- /
!/
      END MODULE W3SNLSMD
