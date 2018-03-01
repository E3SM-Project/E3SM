!/ ------------------------------------------------------------------- /
      MODULE W3SNL1MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    04-Feb-2000 : Origination.                        ( version 2.00 )
!/    09-May-2002 : Switch clean up.                    ( version 2.21 )
!/    24-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Bundles routines calculate nonlinear wave-wave interactions
!     according to the Discrete Interaction Approximation (DIA) of
!     Hasselmann et al. (JPO, 1985).
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SNL1    Subr. Public   Calculate interactions.
!      INSNL1    Subr. Public   Initialization routine.
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
!       !/S      Enable subroutine tracing.
!       !/T(n)   Test output, see subroutines.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      PUBLIC
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SNL1 (A, CG, KDMEAN, S, D)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         24-Dec-2004 |
!/                  +-----------------------------------+
!/
!/    12-Jun-1996 : Final FORTRAN 77                    ( version 1.18 )
!/    04-Feb-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    09-May-2002 : Switch clean up.                    ( version 2.21 )
!/    24-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/
!  1. Purpose :
!
!     Calculate nonlinear interactions and the diagonal term of
!     its derivative.
!
!  2. Method :
!
!     Discrete interaction approximation. (Hasselmann and Hasselmann
!     1985; WAMDI group 1988)
!
!     The DIA is applied to the energy spectrum (instead of the action
!     spectrum), for which is was originally developped. Because the
!     frequency grid is invariant, the nonlinear interactions are
!     calculated for the frequency spectrum, as in WAM. This requires
!     only a single set of interpolation data which can be applied
!     throughout the spatial domain. For deep water this is idenitical
!     to a direct application to the wavenumber spectrum, for shallow
!     water it is not. As the shallow water correction is nothing but
!     a crude approximation, the choice between spectra is expected to
!     be irrelevant.
!
!     The nonlinear interactions are calculated for two "mirror image"
!     quadruplets as described in the manual. The central bin of these
!     quadruples is placed on the discrete complonents of the spectrum,
!     which requires interpolation to obtain other eneregy densities.
!     The figure below defines the diferent basic counters and weights
!     necessary for this interpolation.
!
!               IFRM1  IFRM
!                5        7    T |
!          ITHM1  +------+     H +
!                 |      |     E |      IFRP      IFRP1
!                 |   \  |     T |       3           1
!           ITHM  +------+     A +        +---------+  ITHP1
!                6       \8      |        |         |
!                                |        |  /      |
!                           \    +        +---------+  ITHP
!                                |      /4           2
!                              \ |  /
!          -+-----+------+-------#--------+---------+----------+
!                              / |  \        FREQ.
!                                |      \4           2
!                           /    +        +---------+  ITHP
!                                |        |  \      |
!                6       /8      |        |         |
!           ITHM  +------+       +        +---------+  ITHP1
!                 |   \  |       |       3           1
!                 |      |       |      IFRP      IFRP1
!          ITHM1  +------+       +
!                5        7      |
!
!     To create long vector loops and to efficiently deal with the
!     closed nature of the directional space, the relative counters
!     above are replaced by complete addresses stored in 32 arrays
!     (see section 3 and INSNL1). The interaction are furthermore
!     calucated for an extended spectrum, making it unnecessary to
!     introduce extra weight factors for low and high frequencies.
!     Therefore low and high frequencies are added to the local
!     (auxiliary) spectrum as illustraed below.
!
!              ^  +---+---------------------+---------+- NTH
!              |  |   :                     :         |
!                 |   :                     :         |
!              d  | 2 :  original spectrum  :    1    |
!              i  |   :                     :         |
!              r  |   :                     :         |
!                 +---+---------------------+---------+-  1
!                            Frequencies -->     ^
!         IFR =   0   1                    NFR   |  NFRHGH
!                                                |
!                                             NFRCHG
!
!     where : 1 : Extra tail added beyond NFR
!             2 : Empty bins at low frequencies
!
!             NFRHGH = NFR + IFRP1 - IFRM1
!             NFRCHG = NFR - IFRM1
!
!     All counters and arrays are set in INSNL1. See also section 3
!     and section 8.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action spectrum A(ISP) as a function of
!                         direction (rad)  and wavenumber.
!       CG      R.A.  I   Group velocities (dimension NK).
!       KDMEAN  Real  I   Mean relative depth.
!       S       R.A.  O   Source term.                           *)
!       D       R.A.  O   Diagonal term of derivative.           *)
!     ----------------------------------------------------------------
!                             *) 1-D array with dimension NTH*NK
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      PRT2DS    Subr. W3ARRYMD Print plot of spectra.
!      OUTMAT    Subr. W3WRRYMD Print out 2D matrix.
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
!       None.
!
!  8. Structure :
!
!     -------------------------------------------
!      1.  Calculate proportionality constant.
!      2.  Prepare auxiliary spectrum
!      3.  Calculate (unfolded) interactions
!        a Energy at interacting bins
!        b Contribution to interactions
!        c Fold interactions to side angles
!      4.  Put source and diagonal term together
!     -------------------------------------------
!
!  9. Switches :
!
!     !/S   Enable subroutine tracing.
!     !/T   Enable general test output.
!     !/T0  2-D print plot of source term.
!     !/T1  Print arrays.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, FACHFE,                &
                          KDCON, KDMN, SNLC1, SNLS1, SNLS2, SNLS3
      USE W3ADATMD, ONLY: NFR, NFRHGH, NFRCHG, NSPECX, NSPECY,        &
                    IP11, IP12, IP13, IP14, IM11, IM12, IM13, IM14,   &
                    IP21, IP22, IP23, IP24, IM21, IM22, IM23, IM24,   &
                    IC11, IC12, IC21, IC22, IC31, IC32, IC41, IC42,   &
                    IC51, IC52, IC61, IC62, IC71, IC72, IC81, IC82,   &
                    DAL1, DAL2, DAL3, AF11,                           &
                    AWG1, AWG2, AWG3, AWG4, AWG5, AWG6, AWG7, AWG8,   &
                    SWG1, SWG2, SWG3, SWG4, SWG5, SWG6, SWG7, SWG8
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NSPEC), CG(NK), KDMEAN
      REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ITH, IFR, ISP
      REAL                    :: X, X2, CONS, CONX, FACTOR,           &
                                 E00, EP1, EM1, EP2, EM2,             &
                                 SA1A, SA1B, SA2A, SA2B
      REAL               ::  UE  (1-NTH:NSPECY), SA1 (1-NTH:NSPECX),  &
                             SA2 (1-NTH:NSPECX), DA1C(1-NTH:NSPECX),  &
                             DA1P(1-NTH:NSPECX), DA1M(1-NTH:NSPECX),  &
                             DA2C(1-NTH:NSPECX), DA2P(1-NTH:NSPECX),  &
                             DA2M(1-NTH:NSPECX), CON (      NSPEC )
!/
!/ ------------------------------------------------------------------- /
!/
! initialisations
!
! 1.  Calculate prop. constant --------------------------------------- *
!
      X      = MAX ( KDCON*KDMEAN , KDMN )
      X2     = MAX ( -1.E15, SNLS3*X)
      CONS   = SNLC1 * ( 1. + SNLS1/X * (1.-SNLS2*X) * EXP(X2) )
!
! 2.  Prepare auxiliary spectrum and arrays -------------------------- *
!
      DO IFR=1, NFR
        CONX = TPIINV / SIG(IFR) * CG(IFR)
        DO ITH=1, NTH
          ISP       = ITH + (IFR-1)*NTH
          UE (ISP) = A(ISP) / CONX
          CON(ISP) = CONX
          END DO
        END DO
!
      DO IFR=NFR+1, NFRHGH
        DO ITH=1, NTH
          ISP      = ITH + (IFR-1)*NTH
          UE(ISP) = UE(ISP-NTH) * FACHFE
          END DO
        END DO
!
      DO ISP=1-NTH, 0
        UE  (ISP) = 0.
        SA1 (ISP) = 0.
        SA2 (ISP) = 0.
        DA1C(ISP) = 0.
        DA1P(ISP) = 0.
        DA1M(ISP) = 0.
        DA2C(ISP) = 0.
        DA2P(ISP) = 0.
        DA2M(ISP) = 0.
        END DO
!
! 3.  Calculate interactions for extended spectrum ------------------- *
!
      DO ISP=1, NSPECX
!
! 3.a Energy at interacting bins
!
        E00    =        UE(ISP)
        EP1    = AWG1 * UE(IP11(ISP)) + AWG2 * UE(IP12(ISP))        &
               + AWG3 * UE(IP13(ISP)) + AWG4 * UE(IP14(ISP))
        EM1    = AWG5 * UE(IM11(ISP)) + AWG6 * UE(IM12(ISP))        &
               + AWG7 * UE(IM13(ISP)) + AWG8 * UE(IM14(ISP))
        EP2    = AWG1 * UE(IP21(ISP)) + AWG2 * UE(IP22(ISP))        &
               + AWG3 * UE(IP23(ISP)) + AWG4 * UE(IP24(ISP))
        EM2    = AWG5 * UE(IM21(ISP)) + AWG6 * UE(IM22(ISP))        &
               + AWG7 * UE(IM23(ISP)) + AWG8 * UE(IM24(ISP))
!
! 3.b Contribution to interactions
!
        FACTOR = CONS * AF11(ISP) * E00
!
        SA1A   = E00 * ( EP1*DAL1 + EM1*DAL2 )
        SA1B   = SA1A - EP1*EM1*DAL3
        SA2A   = E00 * ( EP2*DAL1 + EM2*DAL2 )
        SA2B   = SA2A - EP2*EM2*DAL3
!
        SA1 (ISP) = FACTOR * SA1B
        SA2 (ISP) = FACTOR * SA2B
!
        DA1C(ISP) = CONS * AF11(ISP) * ( SA1A + SA1B )
        DA1P(ISP) = FACTOR * ( DAL1*E00 - DAL3*EM1 )
        DA1M(ISP) = FACTOR * ( DAL2*E00 - DAL3*EP1 )
!
        DA2C(ISP) = CONS * AF11(ISP) * ( SA2A + SA2B )
        DA2P(ISP) = FACTOR * ( DAL1*E00 - DAL3*EM2 )
        DA2M(ISP) = FACTOR * ( DAL2*E00 - DAL3*EP2 )
!
        END DO
!
! 4.  Put source and diagonal term together -------------------------- *
!
      DO ISP=1, NSPEC
!
        S(ISP) = CON(ISP) * ( - 2. * ( SA1(ISP) + SA2(ISP) )       &
                   + AWG1 * ( SA1(IC11(ISP)) + SA2(IC12(ISP)) )    &
                   + AWG2 * ( SA1(IC21(ISP)) + SA2(IC22(ISP)) )    &
                   + AWG3 * ( SA1(IC31(ISP)) + SA2(IC32(ISP)) )    &
                   + AWG4 * ( SA1(IC41(ISP)) + SA2(IC42(ISP)) )    &
                   + AWG5 * ( SA1(IC51(ISP)) + SA2(IC52(ISP)) )    &
                   + AWG6 * ( SA1(IC61(ISP)) + SA2(IC62(ISP)) )    &
                   + AWG7 * ( SA1(IC71(ISP)) + SA2(IC72(ISP)) )    &
                   + AWG8 * ( SA1(IC81(ISP)) + SA2(IC82(ISP)) ) )
!
        D(ISP) =  - 2. * ( DA1C(ISP) + DA2C(ISP) )                 &
                + SWG1 * ( DA1P(IC11(ISP)) + DA2P(IC12(ISP)) )     &
                + SWG2 * ( DA1P(IC21(ISP)) + DA2P(IC22(ISP)) )     &
                + SWG3 * ( DA1P(IC31(ISP)) + DA2P(IC32(ISP)) )     &
                + SWG4 * ( DA1P(IC41(ISP)) + DA2P(IC42(ISP)) )     &
                + SWG5 * ( DA1M(IC51(ISP)) + DA2M(IC52(ISP)) )     &
                + SWG6 * ( DA1M(IC61(ISP)) + DA2M(IC62(ISP)) )     &
                + SWG7 * ( DA1M(IC71(ISP)) + DA2M(IC72(ISP)) )     &
                + SWG8 * ( DA1M(IC81(ISP)) + DA2M(IC82(ISP)) )
!
        END DO
!
! ... Test output :
!
      RETURN
!
! Formats
!
!/
!/ End of W3SNL1 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SNL1
!/ ------------------------------------------------------------------- /
      SUBROUTINE INSNL1 ( IMOD )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         24-Dec-2004 |
!/                  +-----------------------------------+
!/
!/    19-Oct-1998 : Final FORTRAN 77                    ( version 1.18 )
!/    04-Feb-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    09-May-2002 : Switch clean up.                    ( version 2.21 )
!/    24-Dec-2004 : Multiple grid version.              ( version 3.06 )
!/
!  1. Purpose :
!
!     Preprocessing for nonlinear interactions (weights).
!
!  2. Method :
!
!     See W3SNL1.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.  I   Model number.
!     ----------------------------------------------------------------
!
!     Local variables
!     ----------------------------------------------------------------
!       ITHxn   Real  Directional indices.                 (relative)
!       IFRxn   Real  Frequency indices.                   (relative)
!       IT1     R.A.  Directional indices.                      (1-D)
!       IFn     R.A.  Frequency indices.                        (1-D)
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
!      W3IOGR    Subr. W3IOGRMD Model definition file processing.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - Check on array dimensions for local arrays in W3SNL.
!
!  7. Remarks :
!
!     - Test output is generated through W3IOGR.
!     - No testing of IMOD ir resetting of pointers.
!
!  8. Structure :
!
!     - See source code.
!
!  9. Switches :
!
!       !/S      Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, DTH, XFR, SIG, LAM
      USE W3ADATMD, ONLY: W3DMNL
      USE W3ADATMD, ONLY: NFR, NFRHGH, NFRCHG, NSPECX, NSPECY,        &
                    IP11, IP12, IP13, IP14, IM11, IM12, IM13, IM14,   &
                    IP21, IP22, IP23, IP24, IM21, IM22, IM23, IM24,   &
                    IC11, IC12, IC21, IC22, IC31, IC32, IC41, IC42,   &
                    IC51, IC52, IC61, IC62, IC71, IC72, IC81, IC82,   &
                    DAL1, DAL2, DAL3, AF11,                           &
                    AWG1, AWG2, AWG3, AWG4, AWG5, AWG6, AWG7, AWG8,   &
                    SWG1, SWG2, SWG3, SWG4, SWG5, SWG6, SWG7, SWG8
      USE W3ODATMD, ONLY: NDST, NDSE
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMOD
!/
!/ Local parameters
!/
      INTEGER                 :: IFR, ITH, ISP, ITHP, ITHP1, ITHM,    &
                                 ITHM1,IFRP, IFRP1, IFRM, IFRM1
      INTEGER, ALLOCATABLE    :: IF1(:), IF2(:), IF3(:), IF4(:),      &
                                 IF5(:), IF6(:), IF7(:), IF8(:),      &
                                 IT1(:), IT2(:), IT3(:), IT4(:),      &
                                 IT5(:), IT6(:), IT7(:), IT8(:)
      REAL                    :: DELTH3, DELTH4, LAMM2, LAMP2, CTHP,  &
                                 WTHP, WTHP1, CTHM, WTHM, WTHM1,      &
                                 XFRLN, WFRP, WFRP1, WFRM, WFRM1, FR, &
                                 AF11A
!/
!/ ------------------------------------------------------------------- /
!/
!
      NFR     = NK
!
! 1.  Internal angles of quadruplet.
!
      LAMM2  = (1.-LAM)**2
      LAMP2  = (1.+LAM)**2
      DELTH3 = ACOS( (LAMM2**2+4.-LAMP2**2) / (4.*LAMM2) )
      DELTH4 = ASIN(-SIN(DELTH3)*LAMM2/LAMP2)
!
! 2.  Lambda dependend weight factors.
!
      DAL1   = 1. / (1.+LAM)**4
      DAL2   = 1. / (1.-LAM)**4
      DAL3   = 2. * DAL1 * DAL2
!
! 3.  Directional indices.
!
      CTHP   = ABS(DELTH4/DTH)
      ITHP   = INT(CTHP)
      ITHP1  = ITHP + 1
      WTHP   = CTHP - REAL(ITHP)
      WTHP1  = 1.- WTHP
!
      CTHM   = ABS(DELTH3/DTH)
      ITHM   = INT(CTHM)
      ITHM1  = ITHM + 1
      WTHM   = CTHM - REAL(ITHM)
      WTHM1  = 1.- WTHM
!
! 4.  Frequency indices.
!
      XFRLN  = LOG(XFR)
!
      IFRP   = INT( LOG(1.+LAM) / XFRLN )
      IFRP1  = IFRP + 1
      WFRP   = (1.+LAM - XFR**IFRP) / (XFR**IFRP1 - XFR**IFRP)
      WFRP1  = 1. - WFRP
!
      IFRM   = INT( LOG(1.-LAM) / XFRLN )
      IFRM1  = IFRM - 1
      WFRM   = (XFR**IFRM -(1.-LAM)) / (XFR**IFRM - XFR**IFRM1)
      WFRM1  = 1. - WFRM
!
! 5.  Range of calculations
!
      NFRHGH = NFR + IFRP1 - IFRM1
      NFRCHG = NFR - IFRM1
      NSPECY = NFRHGH * NTH
      NSPECX = NFRCHG * NTH
!
! 6.  Allocate arrays or check array sizes
!
      CALL W3DMNL ( IMOD, NDSE, NDST, NSPEC, NSPECX )
!
      ALLOCATE ( IF1(NFRCHG), IF2(NFRCHG), IF3(NFRCHG), IF4(NFRCHG),  &
                 IF5(NFRCHG), IF6(NFRCHG), IF7(NFRCHG), IF8(NFRCHG),  &
                 IT1(NTH), IT2(NTH), IT3(NTH), IT4(NTH),              &
                 IT5(NTH), IT6(NTH), IT7(NTH), IT8(NTH) )
!
! 7.  Spectral addresses
!
      DO IFR=1, NFRCHG
        IF1(IFR) =           IFR+IFRP
        IF2(IFR) =           IFR+IFRP1
        IF3(IFR) = MAX ( 0 , IFR+IFRM  )
        IF4(IFR) = MAX ( 0 , IFR+IFRM1 )
        IF5(IFR) = MAX ( 0 , IFR-IFRP  )
        IF6(IFR) = MAX ( 0 , IFR-IFRP1 )
        IF7(IFR) =           IFR-IFRM
        IF8(IFR) =           IFR-IFRM1
        END DO
!
      DO ITH=1, NTH
        IT1(ITH) = ITH + ITHP
        IT2(ITH) = ITH + ITHP1
        IT3(ITH) = ITH + ITHM
        IT4(ITH) = ITH + ITHM1
        IT5(ITH) = ITH - ITHP
        IT6(ITH) = ITH - ITHP1
        IT7(ITH) = ITH - ITHM
        IT8(ITH) = ITH - ITHM1
        IF ( IT1(ITH).GT.NTH) IT1(ITH) = IT1(ITH) - NTH
        IF ( IT2(ITH).GT.NTH) IT2(ITH) = IT2(ITH) - NTH
        IF ( IT3(ITH).GT.NTH) IT3(ITH) = IT3(ITH) - NTH
        IF ( IT4(ITH).GT.NTH) IT4(ITH) = IT4(ITH) - NTH
        IF ( IT5(ITH).LT. 1 ) IT5(ITH) = IT5(ITH) + NTH
        IF ( IT6(ITH).LT. 1 ) IT6(ITH) = IT6(ITH) + NTH
        IF ( IT7(ITH).LT. 1 ) IT7(ITH) = IT7(ITH) + NTH
        IF ( IT8(ITH).LT. 1 ) IT8(ITH) = IT8(ITH) + NTH
        END DO
!
      DO ISP=1, NSPECX
        IFR       = 1 + (ISP-1)/NTH
        ITH       = 1 + MOD(ISP-1,NTH)
        IP11(ISP) = IT2(ITH) + (IF2(IFR)-1)*NTH
        IP12(ISP) = IT1(ITH) + (IF2(IFR)-1)*NTH
        IP13(ISP) = IT2(ITH) + (IF1(IFR)-1)*NTH
        IP14(ISP) = IT1(ITH) + (IF1(IFR)-1)*NTH
        IM11(ISP) = IT8(ITH) + (IF4(IFR)-1)*NTH
        IM12(ISP) = IT7(ITH) + (IF4(IFR)-1)*NTH
        IM13(ISP) = IT8(ITH) + (IF3(IFR)-1)*NTH
        IM14(ISP) = IT7(ITH) + (IF3(IFR)-1)*NTH
        IP21(ISP) = IT6(ITH) + (IF2(IFR)-1)*NTH
        IP22(ISP) = IT5(ITH) + (IF2(IFR)-1)*NTH
        IP23(ISP) = IT6(ITH) + (IF1(IFR)-1)*NTH
        IP24(ISP) = IT5(ITH) + (IF1(IFR)-1)*NTH
        IM21(ISP) = IT4(ITH) + (IF4(IFR)-1)*NTH
        IM22(ISP) = IT3(ITH) + (IF4(IFR)-1)*NTH
        IM23(ISP) = IT4(ITH) + (IF3(IFR)-1)*NTH
        IM24(ISP) = IT3(ITH) + (IF3(IFR)-1)*NTH
        END DO
!
      DO ISP=1, NSPEC
        IFR       = 1 + (ISP-1)/NTH
        ITH       = 1 + MOD(ISP-1,NTH)
        IC11(ISP) = IT6(ITH) + (IF6(IFR)-1)*NTH
        IC21(ISP) = IT5(ITH) + (IF6(IFR)-1)*NTH
        IC31(ISP) = IT6(ITH) + (IF5(IFR)-1)*NTH
        IC41(ISP) = IT5(ITH) + (IF5(IFR)-1)*NTH
        IC51(ISP) = IT4(ITH) + (IF8(IFR)-1)*NTH
        IC61(ISP) = IT3(ITH) + (IF8(IFR)-1)*NTH
        IC71(ISP) = IT4(ITH) + (IF7(IFR)-1)*NTH
        IC81(ISP) = IT3(ITH) + (IF7(IFR)-1)*NTH
        IC12(ISP) = IT2(ITH) + (IF6(IFR)-1)*NTH
        IC22(ISP) = IT1(ITH) + (IF6(IFR)-1)*NTH
        IC32(ISP) = IT2(ITH) + (IF5(IFR)-1)*NTH
        IC42(ISP) = IT1(ITH) + (IF5(IFR)-1)*NTH
        IC52(ISP) = IT8(ITH) + (IF8(IFR)-1)*NTH
        IC62(ISP) = IT7(ITH) + (IF8(IFR)-1)*NTH
        IC72(ISP) = IT8(ITH) + (IF7(IFR)-1)*NTH
        IC82(ISP) = IT7(ITH) + (IF7(IFR)-1)*NTH
        END DO
!
      DEALLOCATE ( IF1, IF2, IF3, IF4, IF5, IF6, IF7, IF8,  &
                   IT1, IT2, IT3, IT4, IT5, IT6, IT7, IT8 )
!
! 8.  Fill scaling array (f**11)
!
      DO IFR=1, NFR
        AF11A  = (SIG(IFR)*TPIINV)**11
        DO ITH=1, NTH
          AF11(ITH+(IFR-1)*NTH) = AF11A
          END DO
        END DO
!
      FR     = SIG(NFR)*TPIINV
      DO IFR=NFR+1, NFRCHG
        FR     = FR * XFR
        AF11A  = FR**11
        DO ITH=1, NTH
          AF11(ITH+(IFR-1)*NTH) = AF11A
          END DO
        END DO
!
! 9.  Interpolation weights
!
      AWG1   = WTHP  * WFRP
      AWG2   = WTHP1 * WFRP
      AWG3   = WTHP  * WFRP1
      AWG4   = WTHP1 * WFRP1
      AWG5   = WTHM  * WFRM
      AWG6   = WTHM1 * WFRM
      AWG7   = WTHM  * WFRM1
      AWG8   = WTHM1 * WFRM1
!
      SWG1   = AWG1**2
      SWG2   = AWG2**2
      SWG3   = AWG3**2
      SWG4   = AWG4**2
      SWG5   = AWG5**2
      SWG6   = AWG6**2
      SWG7   = AWG7**2
      SWG8   = AWG8**2
!
      RETURN
!
! Formats
!
!/
!/ End of INSNL1 ----------------------------------------------------- /
!/
      END SUBROUTINE INSNL1
!/
!/ End of module W3SNL1MD -------------------------------------------- /
!/
      END MODULE W3SNL1MD
