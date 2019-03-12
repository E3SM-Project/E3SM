#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3SRCEMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Dec-2015 |
!/                  +-----------------------------------+
!/
!/    For updates see subroutine.
!/
!  1. Purpose :
!
!     Source term integration routine.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      OFFSET    R.P.  Private  Offset in time integration scheme.
!                               0.5 in original WAM, now 1.0
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3SRCE    Subr. Public   Calculate and integrate source terms.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!     See corresponding documentation of W3SRCE.
!
!  5. Remarks :
!
!  6. Switches :
!
!       See section 9 of W3SRCE.
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
!/
      REAL, PARAMETER, PRIVATE:: OFFSET = 1.
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SRCE ( IT, IX, IY, IMOD, SPEC, ALPHA, WN1, CG1,    &
                          D_INP, U10ABS, U10DIR, AS, USTAR, USTDIR,   &
                          CX, CY,  ICE, ICEH, ICEF, ICEDMAX,          &
                          REFLEC, REFLED, DELX, DELY, DELA, TRNX,     &
                          TRNY, BERG, FPI, DTDYN, FCUT, DTG, TAUWX,   &
                          TAUWY, TAUOX, TAUOY, TAUWIX, TAUWIY, TAUWNX,&
                          TAUWNY, PHIAW, CHARN, PHIOC, WHITECAP, D50, &
                          PSIC, BEDFORM , PHIBBL, TAUBBL, TAUICE,     &
                          PHICE, COEF)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Dec-2015 |
!/                  +-----------------------------------+
!/
!/    06-Dec-1996 : Final FORTRAN 77                    ( version 1.18 )
!/    04-Feb-2000 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    14-Feb-2000 : Exact-NL added                      ( version 2.01 )
!/    04-May-2000 : Non-central integration             ( version 2.03 )
!/    02-Feb-2001 : Xnl version 3.0                     ( version 2.07 )
!/    09-May-2002 : Switch clean up.                    ( version 2.21 )
!/    13-Nov-2002 : Add stress vector.                  ( version 3.00 )
!/    27-Nov-2002 : First version of VDIA and MDIA.     ( version 3.01 )
!/    07-Oct-2003 : Output options for NN training.     ( version 3.05 )
!/    24-Dec-2004 : Multiple model version.             ( version 3.06 )
!/    23-Jun-2006 : Linear input added.                 ( version 3.09 )
!/    27-Jun-2006 : Adding file name preamble.          ( version 3.09 )
!/    04-Jul-2006 : Separation of stress computation.   ( version 3.09 )
!/    16-Apr-2007 : Miche style limiter added.          ( version 3.11 )
!/                  (J. H. Alves)
!/    25-Apr-2007 : Battjes-Janssen Sdb added.          ( version 3.11 )
!/                  (J. H. Alves)
!/    09-Oct-2007 : Adding WAM 4+ and SB1 options.      ( version 3.13 )
!/                  (F. Ardhuin)
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/    19-Aug-2010 : Making treatment of 0 water depth   ( version 3.14.6 )
!/                  consistent with the rest of the model.
!/    31-Mar-2010 : Adding ice conc. and reflections    ( version 3.14.4 )
!/    15-May-2010 : Adding transparencies               ( version 3.14.4 )
!/    01-Jun-2011 : Movable bed bottom friction in BT4  ( version 4.01 )
!/    01-Jul-2011 : Energy and momentum flux, friction  ( version 4.01 )
!/    24-Aug-2011 : Uses true depth for depth-induced   ( version 4.04 )
!/    16-Sep-2011 : Initialization of TAUWAX, TAUWAY    ( version 4.04 )
!/     1-Dec-2011 : Adding BYDRZ source term package    ( version 4.04 )
!/                  ST6 and optional Hwang (2011)
!/                  stresses FLX4.
!/    14-Mar-2012 : Update of BT4, passing PSIC         ( version 4.04 )
!/    13-Jul-2012 : Move GMD (SNL3) and nonlinear filter (SNLS)
!/                  from 3.15 (HLT).                    ( version 4.08 )
!/    28-Aug-2013 : Corrected MLIM application          ( version 4.11 )
!/    10-Sep-2013 : Special treatment for IG band       ( version 4.15 )
!/    14-Nov-2013 : Make orphaned pars in par lst local ( version 4.13 )
!/    17-Nov-2013 : Coupling fraction of ice-free       ( version 4.13 )
!/                  surface to SIN and SDS. (S. Zieger)
!/    01-Avr-2014 : Adding ice thickness and floe size  ( version 4.18 )
!/    23-May-2014 : Adding ice fluxes to W3SRCE         ( version 5.01 )
!/    27-Aug-2015 : Adding inputs to function W3SIS2    ( version 5.10 )
!/    13-Dec-2015 : Exact integration of ice terms      ( version 5.10 )
!/
!/    Copyright 2009-2013 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Calculate and integrate source terms for a single grid point.
!
!  2. Method :
!
!     Physics  : see manual and corresponding subroutines.
!
!     Numerics :
!
!     Dynamic-implicit integration of the source terms based on
!     WW-II (Tolman 1992). The dynamic time step is calculated
!     given a maximum allowed change of spectral densities for
!     frequencies / wavenumbers below the usual cut-off.
!     The maximum change is given by the minimum of a parametric
!     and a relative change. The parametric change relates to a
!     PM type equilibrium range
!
!                                -1  (2pi)**4       1
!       dN(k)     =  Xp alpha  pi   ---------- ------------
!            max                       g**2     k**3 sigma
!
!                              1                                     .
!                 =  FACP ------------                              (1)
!                          k**3 sigma                                .
!
!     where
!           alpha = 0.62e-4                       (set in W3GRID)
!           Xp      fraction of PM shape          (read in W3GRID)
!           FACP    combined factor               (set in W3GRID)
!
!     The maximum relative change is given as
!
!                           /            +-                  -+ \    .
!       dN(k)     =  Xr max | N(k) , max | Nx , Xfilt N(k)    | |   (2)
!            max            \            +-               max-+ /    .
!
!     where
!           Xr      fraction of relative change   (read in W3GRID)
!           Xfilt   filter level                  (read in W3GRID)
!           Nx      Maximum parametric change (1)
!                   for largest wavenumber.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IX,IY   Int.   I   Discrete grid point counters.
!       IMOD    Int.   I   Model number.
!       SPEC    R.A.  I/O  Spectrum (action) in 1-D form.
!       ALPHA   R.A.  I/O  Nondimenional 1-D spectrum corresponding
!                          to above full spectra (Phillip's const.).
!                          Calculated separately for numerical
!                          economy on vector machine (W3SPR2).
!       WN1     R.A.   I   Discrete wavenumbers.
!       CG1     R.A.   I   Id. group velocities.
!       D_INP   Real.  I   Depth. Compared to DMIN to get DEPTH.
!       U10ABS  Real.  I   Wind speed at reference height.
!       U10DIR  Real.  I   Id. wind direction.
!       AS      Real.  I   Air-sea temp. difference.      ( !/ST3 )
!       USTAR   Real. !/O  Friction velocity.
!       USTDIR  Real  !/O  Idem, direction.
!       CX-Y    Real.  I   Current velocity components.   ( !/BS1 )
!       ICE     Real   I   Sea ice concentration
!       ICEH    Real   I   Sea ice thickness
!       ICEF    Real  I/O  Sea ice maximum floe diameter  (updated)
!       ICEDMAX Real  I/O  Sea ice maximum floe diameter
!       BERG    Real   I   Iceberg damping coefficient    ( !/BS1 )
!       REFLEC  R.A.   I   reflection coefficients        ( !/BS1 )
!       REFLED  I.A.   I   reflection direction           ( !/BS1 )
!       TRNX-Y  Real   I   Grid transparency in X and Y   ( !/BS1 )
!       DELX    Real.  I   grid cell size in X direction  ( !/BS1 )
!       DELY    Real.  I   grid cell size in Y direction  ( !/BS1 )
!       DELA    Real.  I   grid cell area                 ( !/BS1 )
!       FPI     Real  I/O  Peak-input frequency.          ( !/ST2 )
!      WHITECAP R.A.   O   Whitecap statisics             ( !/ST4 )
!       DTDYN   Real   O   Average dynamic time step.
!       FCUT    Real   O   Cut-off frequency for tail.
!       DTG     Real   I   Global time step.
!       D50     Real   I   Sand grain size                ( !/BT4 )
!       BEDFORM R.A.  I/O  Bedform parameters             ( !/BT4 )
!       PSIC    Real   I   Critical Shields               ( !/BT4 )
!       PHIBBL  Real   O   Energy flux to BBL             ( !/BTx )
!       TAUBBL  R.A.   O   Momentum flux to BBL           ( !/BTx )
!       TAUICE  R.A.   O   Momentum flux to sea ice       ( !/ICx )
!       PHICE   Real   O   Energy flux to sea ice         ( !/ICx )
!     ----------------------------------------------------------------
!       Note: several pars are set to I/O to avoid compiler warnings.
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3SPRn    Subr. W3SRCnMD Mean wave parameters for use in
!                               source terms.
!      W3FLXn    Subr. W3FLXnMD Flux/stress computation.
!      W3SLNn    Subr. W3SLNnMD Linear input.
!      W3SINn    Subr. W3SRCnMD Input source term.
!      W3SNLn    Subr. W3SNLnMD Nonlinear interactions.
!      W3SNLS    Subr. W3SNLSMD Nonlinear smoother.
!      W3SDSn    Subr. W3SRCnMD Whitecapping source term
!      W3SBTn    Subr. W3SBTnMD Bottom friction source term.
!      W3SDBn    Subr. W3SBTnMD Depth induced breaking source term.
!      W3STRn    Subr. W3STRnMD Triad interaction source term.
!      W3SBSn    Subr. W3SBSnMD Bottom scattering source term.
!      W3REFn    Subr. W3REFnMD Reflexions (shore, icebergs ...).
!      W3SXXn    Subr. W3SXXnMD Unclassified source term.
!      STRACE    Subr. W3SERVMD Subroutine tracing (!/S)
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3WAVE    Subr. W3WAVEMD Actual wave model routine.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!       None.
!
!  7. Remarks :
!
!     - No testing is performed on the status of the grid point.
!
!  8. Structure :
!
!     -----------------------------------------------------------------
!       1.   Preparations
!         a  Set maximum change and wavenumber arrays.
!         b  Prepare dynamic time stepping.
!         c  Compute mean parameters.                       ( W3SPRn )
!         d  Compute stresses (if posible).
!         e  Prepare cut-off
!         f  Test output for !/NNT option.
!     --start-dynamic-integration-loop---------------------------------
!       2.  Calculate source terms
!         a Input.                                  ( W3SLNx, W3SINn )
!         b Nonlinear interactions.                         ( W3SNLn )
!         c Dissipation                                     ( W3SDSn )
!           1 as included in source terms                   ( W3SDSn )
!           2 optional dissipation due to different physics ( W3SWLn )
!         d Bottom friction.                                ( W3SBTn )
!       3.  Calculate cut-off frequencie(s)
!       4.  Summation of source terms and diagonal term and time step.
!       5.  Increment spectrum.
!       6.  Add tail
!         a Mean wave parameters and cut-off                ( W3SPRn )
!         b 'Seeding' of spectrum.                          ( !/SEED )
!         c Add tail
!       7.  Check if integration complete.
!     --end-dynamic-integration-loop-----------------------------------
!       8.  Save integration data.
!     -----------------------------------------------------------------
!
!  9. Switches :
!
!     !/FLX1  Wu (1980) stress computation.              ( Choose one )
!     !/FLX2  T&C (1996) stress computation.
!     !/FLX3  T&C (1996) stress computation with cap.
!     !/FLX4  Hwang (2011) stress computation (2nd order).
!
!     !/LN0   No linear input.                           ( Choose one )
!     !/LNX   User-defined bottom friction.
!
!     !/ST0   No input and dissipation.                  ( Choose one )
!     !/ST1   WAM-3 input and dissipation.
!     !/ST2   Tolman and Chalikov (1996)  input and dissipation.
!     !/ST3   WAM 4+ input and dissipation.
!     !/ST4   Ardhuin et al. (2009, 2010)
!     !/ST6   BYDB source terms after Babanin, Young, Donelan and Banner.
!     !/STX   User-defined input and dissipation.
!
!     !/NL0   No nonlinear interactions.                 ( Choose one )
!     !/NL1   Discrete interaction approximation.
!     !/NL2   Exact nonlinear interactions.
!     !/NL3   Generalized Multiple DIA.
!     !/NL4   Two Scale Approximation
!     !/NLX   User-defined nonlinear interactions.
!     !/NLS   Nonlinear HF smoother.
!
!     !/BT0   No bottom friction.                        ( Choose one )
!     !/BT1   JONSWAP bottom friction.
!     !/BT4   Bottom friction using movable bed roughness
!                  (Tolman 1994, Ardhuin & al. 2003)
!     !/BT8   Muddy bed (Dalrymple & Liu).
!     !/BT9   Muddy bed (Ng).
!     !/BTX   User-defined bottom friction.
!
!     !/IC1   Dissipation via interaction with ice according to simple
!             methods: 1) uniform in frequency or
!     !/IC2            2) Liu et al. model
!     !/IC3   Dissipation via interaction with ice according to a
!             viscoelastic sea ice model (Wang and Shen 2010).
!     !/IC4   Dissipation via interaction with ice as a function of freq.
!             (empirical/parametric methods)
!     !/DB0   No depth-limited breaking.                 ( Choose one )
!     !/DB1   Battjes-Janssen depth-limited breaking.
!     !/DBX   User-defined bottom friction.
!
!     !/TR0   No triad interactions.                     ( Choose one )
!     !/TR1   Lumped Triad Approximation (LTA).
!     !/TRX   User-defined triad interactions.
!
!     !/BS0   No bottom scattering.                      ( Choose one )
!     !/BS1   Scattering term by Ardhuin and Magne (2007).
!     !/BSX   User-defined bottom scattering
!
!     !/XX0   No arbitrary additional source term.       ( Choose one )
!     !/XXX   User-defined bottom friction.
!
!     !/MLIM  Miche style limiter for shallow water and steepness.
!
!     !/SEED  'Seeding' of lowest frequency for suffuciently strong
!             winds.
!
!     !/NNT   Write output to file test_data_NNN.ww3 for NN training.
!
!     !/S     Enable subroutine tracing.
!     !/T     Enable general test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, TH, DMIN, DTMAX,       &
                          DTMIN, FACTI1, FACTI2, FACSD, FACHFA, FACP, &
                          XFC, XFLT, XREL, XFT, FXFM, FXPM, DDEN,     &
                          FTE, FTF, FHMAX, ECOS, ESIN, IICEDISP,      &
                          ICE100WIND
      USE W3WDATMD, ONLY: TIME
      USE W3ODATMD, ONLY: NDSE, NDST, IAPROC
      USE W3IDATMD, ONLY: INFLAGS1, INFLAGS2, ICEP2
      USE W3DISPMD
      USE W3SLN1MD
      USE W3SRC4MD, ONLY : W3SPR4, W3SIN4, W3SDS4
      USE W3GDATMD, ONLY : ZZWND, FFXFM, FFXPM, FFXFA, FFXFI, FFXFD
      USE W3SNL1MD
      USE W3SBT1MD
      USE W3SDB1MD
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IT, IX, IY, IMOD
      REAL, INTENT(IN)        :: D_INP, U10ABS,     &
                                 U10DIR, AS, CX, CY, DTG, D50,PSIC,   &
                                 ICE, ICEH
      INTEGER, INTENT(IN)     :: REFLED(6)
      REAL, INTENT(IN)        :: REFLEC(4), DELX, DELY, DELA,         &
                                 TRNX, TRNY, BERG, ICEDMAX
      REAL, INTENT(INOUT)     :: WN1(NK), CG1(NK), &
                                 SPEC(NSPEC), ALPHA(NK), USTAR,       &
                                 USTDIR, FPI, TAUOX, TAUOY,           &
                                 TAUWX, TAUWY, PHIAW, PHIOC, PHICE,   &
                                 CHARN, BEDFORM(3), PHIBBL,           &
                                 TAUBBL(2), TAUICE(2), WHITECAP(4),   &
                                 TAUWIX, TAUWIY, TAUWNX, TAUWNY,      &
                                 ICEF
      REAL, INTENT(OUT)       :: DTDYN, FCUT
      REAL, INTENT(IN)        :: COEF
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IK, ITH, IS, IS0, NSTEPS,  NKH, NKH1,&
                                 IKS1, IS1, NSPECH, IDT, IERR, NKI, NKD
      REAL                    :: DTTOT, FHIGH, DT, AFILT, DAMAX, AFAC,&
                                 HDT, ZWND, FP, DEPTH, TAUSCX, TAUSCY, FHIGI
      REAL                    :: EMEAN, FMEAN, WNMEAN, AMAX, CD, Z0, SCAT
      REAL                    :: WN_R(NK), CG_ICE(NK),ALPHA_LIU(NK), ICECOEF2
      DOUBLE PRECISION        :: ATT
      REAL                    :: FMEANS, FH1, FH2, FAGE
      REAL                    :: QCERR  = 0.     !/XNL2 and !/NNT
      REAL                    :: HM, EM
      REAL                    :: EBAND, DIFF, EFINISH, HSTOT, PHINL,       &
                                 FMEAN1, FMEANWS, MWXINIT, MWYINIT,        &
                                 FACTOR, FACTOR2, DRAT, TAUWAX, TAUWAY,    &
                                 MWXFINISH, MWYFINISH, A1BAND, B1BAND,     &
                                 COSI(2)
      REAL                    :: SPECINIT(NSPEC), SPEC2(NSPEC)
      REAL                    :: DAM (NSPEC), WN2 (NSPEC),            &
                                 VSLN(NSPEC),                         &
                                 VSIN(NSPEC), VDIN(NSPEC),            &
                                 VSNL(NSPEC), VDNL(NSPEC),            &
                                 VSDS(NSPEC), VDDS(NSPEC),            &
                                 VSBT(NSPEC), VDBT(NSPEC),            &
                                 VSDB(NSPEC), VDDB(NSPEC),            &
                                 VS  (NSPEC), VD  (NSPEC), EB(NK)
      LOGICAL                 :: LLWS(NSPEC)
      REAL                    :: BRLAMBDA(NSPEC)
      REAL                    :: FOUT(NK,NTH), SOUT(NK,NTH), DOUT(NK,NTH)
      LOGICAL, SAVE           :: FLTEST = .FALSE., FLAGNN = .TRUE.
      LOGICAL                 :: SHAVE
      LOGICAL, SAVE           :: FIRST = .TRUE.
 
!/
!/ ------------------------------------------------------------------- /
!/
!
      DEPTH  = MAX ( DMIN , D_INP )
      IKS1 = 1
      IS1=(IKS1-1)*NTH+1
!
      VSBT = 0.
      VDBT = 0.
!
      ZWND   = ZZWND
!
       DRAT  = DAIR / DWAT
!
! 1.  Preparations --------------------------------------------------- *
!
! 1.a Set maximum change and wavenumber arrays.
!
      DO IK=1, NK
        DAM(1+(IK-1)*NTH) = FACP / ( SIG(IK) * WN1(IK)**3 )
        WN2(1+(IK-1)*NTH) = WN1(IK)
        END DO
!
      DO IK=1, NK
        IS0    = (IK-1)*NTH
        DO ITH=2, NTH
          DAM(ITH+IS0) = DAM(1+IS0)
          WN2(ITH+IS0) = WN2(1+IS0)
          END DO
        END DO
!
! 1.b Prepare dynamic time stepping
!
      DTDYN  = 0.
      DTTOT  = 0.
      NSTEPS = 0
      PHIAW  = 0.
      CHARN  = 0.
      PHINL  = 0.
      PHIBBL = 0.
      TAUWIX = 0.
      TAUWIY = 0.
      TAUWNX = 0.
      TAUWNY = 0.
      TAUWAX = 0.
      TAUWAY = 0.
      TAUSCX = 0.
      TAUSCY = 0.
      TAUBBL(:)= 0.
      TAUICE(:)= 0.
      PHICE  = 0.
      BRLAMBDA(:)=0.
!
! 1.c Set mean parameters
!
      TAUWX=0.
      TAUWY=0.
      IF ( FIRST ) THEN
          LLWS(:) = .TRUE.
          USTAR=0.
          USTDIR=0.
      ELSE
        CALL W3SPR4 (SPEC, CG1, WN1, EMEAN, FMEAN, FMEAN1, WNMEAN, &
                   AMAX, U10ABS, U10DIR, USTAR, USTDIR,            &
                   TAUWX, TAUWY, CD, Z0, CHARN, LLWS, FMEANWS)
        CALL W3SIN4 ( SPEC, CG1, WN2, U10ABS, USTAR, DRAT, AS,       &
                 U10DIR, Z0, CD, TAUWX, TAUWY, TAUWAX, TAUWAY,       &
                 VSIN, VDIN, LLWS, IX, IY, BRLAMBDA )
        END IF
      CALL W3SPR4 (SPEC, CG1, WN1, EMEAN, FMEAN, FMEAN1, WNMEAN, &
                   AMAX, U10ABS, U10DIR, USTAR, USTDIR,          &
                   TAUWX, TAUWY, CD, Z0, CHARN, LLWS, FMEANWS)
!
! 1.c2 Stores the initial data
!
      SPECINIT = SPEC
!
! 1.d Stresses
!
! 1.e Prepare cut-off
!
! !/ST4      FAGE   = FFXFA*TANH(0.3*U10ABS*FMEANWS*TPI/GRAV)
      FAGE   = 0.
      FHIGH  = MAX( (FFXFM + FAGE ) * MAX(FMEAN1,FMEANWS),FFXPM / USTAR)
      FHIGI  = FFXFA * FMEAN1
!
! 1.f Prepare output file for !/NNT option
!
! ... Branch point dynamic integration - - - - - - - - - - - - - - - -
!
      DO
!
        NSTEPS = NSTEPS + 1
!
! 2.  Calculate source terms ----------------------------------------- *
!
! 2.a Input.
!
        CALL W3SLN1 (       WN1, FHIGH, USTAR, U10DIR , VSLN       )
!
        CALL W3SIN4 ( SPEC, CG1, WN2, U10ABS, USTAR, DRAT, AS,       &
                 U10DIR, Z0, CD, TAUWX, TAUWY, TAUWAX, TAUWAY,       &
                 VSIN, VDIN, LLWS, IX, IY, BRLAMBDA )
!
! 2.b Nonlinear interactions.
!
        CALL W3SNL1 ( SPEC, CG1, WNMEAN*DEPTH,        VSNL, VDNL )
!
! 2.c Dissipation... except for ST4
! 2.c1   as in source term package
!
        CALL W3SDS4 ( SPEC, WN1, CG1, USTAR, USTDIR, DEPTH, VSDS,    &
                      VDDS, IX, IY, BRLAMBDA, WHITECAP )
!
        CALL W3SDB1 ( SPEC, WN2, DEPTH, EMEAN, FMEAN, WNMEAN,        &
                                                        VSDB, VDDB )
!
! 2.c2   optional dissipation parameterisations
!
! 2.d Bottom interactions.
!
        CALL W3SBT1 ( SPEC, CG1, WN1, DEPTH,            VSBT, VDBT )
!
! 2.e Additional sources.
!
 
!
! 2.f Dump training data if necessary
!
! 3.  Set frequency cut-off ------------------------------------------ *
!
        NKH    = MIN ( NK , INT(FACTI2+FACTI1*LOG(MAX(1.E-7,FHIGH))) )
        NKH1   = MIN ( NK , NKH+1 )
        NSPECH = NKH1*NTH
!
! 4.  Summation of source terms and diagonal term and time step ------ *
!
        DT     = MIN ( DTG-DTTOT , DTMAX )
        AFILT  = MAX ( DAM(NSPEC) , XFLT*AMAX )
!
!     For input and dissipation calculate the fraction of the ice-free
!     surface. In the presence of ice, the effective water surface
!     is reduce to a fraction of the cell size free from ice, and so is
!     input and dissipation:
!             SIN = (1-ICE)*SIN and SDS=(1-ICE)*SDS ------------------ *
!     INFLAGS2(4) is true if ice concentration was ever read during
!             this simulation
        IF ( INFLAGS2(4) ) THEN
          VSLN(1:NSPECH) = (MAX(0.,MIN(1.,1.-ICE)))**2 * VSLN(1:NSPECH)
          VSIN(1:NSPECH) = MAX(0.,MIN(1.,1.-ICE*ICE100WIND)) &
                           * VSIN(1:NSPECH)
          VDIN(1:NSPECH) = MAX(0.,MIN(1.,1.-ICE*ICE100WIND)) &
                           * VDIN(1:NSPECH)
          VSDS(1:NSPECH) = MAX(0.,MIN(1.,1.-ICE)) * VSDS(1:NSPECH)
          VDDS(1:NSPECH) = MAX(0.,MIN(1.,1.-ICE)) * VDDS(1:NSPECH)
          END IF
!
        NKI    = MAX ( 2 , MIN ( NKH1 ,                           &
                 INT ( FACTI2 + FACTI1*LOG(MAX(1.E-7,FFXFI* FMEAN1)) ) ) )
        DO IS=IS1, NSPECH
          VS(IS) = VSLN(IS) + VSIN(IS) + VSNL(IS)  &
                 + VSDS(IS) + VSBT(IS)
          VS(IS) = VS(IS) + VSDB(IS)
          VD(IS) =            VDIN(IS) + VDNL(IS)  &
                 + VDDS(IS) + VDBT(IS)
          VD(IS) = VD(IS) + VDDB(IS)
          DAMAX  = MIN ( DAM(IS) , MAX ( XREL*SPEC(IS) , AFILT ) )
          AFAC   = 1. / MAX( 1.E-10 , ABS(VS(IS)/DAMAX) )
          DT     = MIN ( DT , AFAC / ( MAX ( 1.E-10,                  &
                         1. + OFFSET*AFAC*MIN(0.,VD(IS)) ) ) )
          END DO  ! end of loop on IS
!
        DT     = MAX ( 0.5 , DT )
!
        DTDYN  = DTDYN + DT
        IDT    = 1 + INT ( 0.99*(DTG-DTTOT)/DT )
        DT     = (DTG-DTTOT)/REAL(IDT)
        SHAVE  = DT.LT.DTMIN .AND. DT.LT.DTG-DTTOT
        DT     = MAX ( DT , MIN (DTMIN,DTG-DTTOT) )
        HDT    = OFFSET * DT
        DTTOT  = DTTOT + DT
 
!
! 5.  Increment spectrum --------------------------------------------- *
!
        IF ( SHAVE ) THEN
!
            DO IS=IS1, NSPECH
              VS(IS) = VS(IS) * DT / MAX ( 1. , (1.-HDT*VD(IS)))
              VS(IS) = SIGN ( MIN (DAM(IS),ABS(VS(IS))) , VS(IS) )
              SPEC(IS) = MAX ( 0. , SPEC(IS)+VS(IS) )
              END DO
!
          ELSE
!
            DO IS=IS1, NSPECH
              VS(IS) = VS(IS) * DT / MAX ( 1. , (1.-HDT*VD(IS)))
              SPEC(IS) = MAX ( 0. , SPEC(IS)+VS(IS) )
              END DO
!
          END IF
!
! 5.b  Computes
!              atmos->wave flux PHIAW-------------------------------- *
!              wave ->BBL  flux PHIBBL------------------------------- *
!              wave ->ice  flux PHICE ------------------------------- *
!
       WHITECAP(3)=0.
       HSTOT=0.
       DO IK=IKS1, NK
         FACTOR = DDEN(IK)/CG1(IK)                    !Jacobian to get energy in band
         FACTOR2= FACTOR*GRAV*WN1(IK)/SIG(IK)         ! coefficient to get momentum
 
         ! Wave direction is "direction to"
         ! therefore there is a PLUS sign for the stress
         DO ITH=1, NTH
           IS   = (IK-1)*NTH + ITH
           COSI(1)=ECOS(IS)
           COSI(2)=ESIN(IS)
           PHIAW = PHIAW + (VSIN(IS))* DT * FACTOR                    &
             / MAX ( 1. , (1.-HDT*VDIN(IS))) ! semi-implict integration scheme
 
           PHIBBL= PHIBBL- (VSBT(IS))* DT * FACTOR                    &
             / MAX ( 1. , (1.-HDT*VDBT(IS))) ! semi-implict integration scheme
           PHINL = PHINL + VSNL(IS)* DT * FACTOR                      &
             / MAX ( 1. , (1.-HDT*VDNL(IS))) ! semi-implict integration scheme
           IF (VSIN(IS).GT.0.) THEN
             WHITECAP(3) = WHITECAP(3) + SPEC(IS)  * FACTOR
           ELSE
           ! computes the upward energy flux (counted > 0 upward)
             CHARN = CHARN - (VSIN(IS))* DT * FACTOR                    &
              / MAX ( 1. , (1.-HDT*VDIN(IS))) ! semi-implict integration scheme
             END IF
           HSTOT = HSTOT + SPEC(IS) * FACTOR
           END DO
         END DO
       WHITECAP(3)=4.*SQRT(WHITECAP(3))
       HSTOT=4.*SQRT(HSTOT)
       TAUWIX= TAUWIX+ TAUWX * DRAT *DT
       TAUWIY= TAUWIY+ TAUWY * DRAT *DT
       TAUWNX= TAUWNX+ TAUWAX * DRAT *DT
       TAUWNY= TAUWNY+ TAUWAY * DRAT *DT
       ! MISSING: TAIL TO BE ADDED ?
!
! 6.  Add tail ------------------------------------------------------- *
!   a Mean parameters
!
        CALL W3SPR4 (SPEC, CG1, WN1, EMEAN, FMEAN, FMEAN1, WNMEAN,&
                   AMAX, U10ABS, U10DIR, USTAR, USTDIR,           &
                   TAUWX, TAUWY, CD, Z0, CHARN, LLWS, FMEANWS)
!
        FAGE   = FFXFA*TANH(0.3*U10ABS*FMEANWS*TPI/GRAV)
        FH1    = (FFXFM+FAGE) * FMEAN1
 
        FH2    = FFXPM / USTAR
        FHIGH  = MIN ( SIG(NK) , MAX ( FH1 , FH2 ) )
        NKH    = MAX ( 2 , MIN ( NKH1 ,                           &
                 INT ( FACTI2 + FACTI1*LOG(MAX(1.E-7,FHIGH)) ) ) )
!
! 6.b Limiter for shallow water or Miche style criterion
!     Last time step only !
!     uses true depth (D_INP) instead of limited depth
!
        IF ( DTTOT .GE. 0.9999*DTG ) THEN
            HM     = FHMAX *TANH(WNMEAN*MAX(0.,D_INP)) / MAX(1.E-4,WNMEAN )
            EM     = HM * HM / 16.
            IF ( EMEAN.GT.EM .AND. EMEAN.GT.1.E-30 ) THEN
                SPEC   = SPEC / EMEAN * EM
                EMEAN  = EM
              END IF
          END IF
!
! 6.c Seeding of spectrum
!     alpha = 0.005 , 0.5 in eq., 0.25 for directional distribution
!
! 6.d Add tail
!
        DO IK=NKH+1, NK
          DO ITH=1, NTH
            SPEC(ITH+(IK-1)*NTH) = SPEC(ITH+(IK-2)*NTH) * FACHFA         &
                       + 0.
            END DO
          END DO
!
! 6.e  Update wave-supported stress----------------------------------- *
!
        CALL W3SIN4 ( SPEC, CG1, WN2, U10ABS, USTAR, DRAT, AS,      &
                      U10DIR, Z0, CD, TAUWX, TAUWY, TAUWAX, TAUWAY, &
                      VSIN, VDIN, LLWS, IX, IY, BRLAMBDA )
 
!
! 7.  Check if integration complete ---------------------------------- *
!
        IF ( DTTOT .GE. 0.9999*DTG ) EXIT
        END DO
!
! ... End point dynamic integration - - - - - - - - - - - - - - - - - -
!
! 8.  Save integration data ------------------------------------------ *
!
      DTDYN  = DTDYN / REAL(MAX(1,NSTEPS))
      FCUT   = FHIGH * TPIINV
!
      GOTO 888
!
! Error escape locations
!
  888 CONTINUE
!
! 9.a  Computes PHIOC------------------------------------------ *
!     The wave to ocean flux is the difference between initial energy
!     and final energy, plus wind input plus the SNL flux to high freq.,
!     minus the energy lost to the bottom boundary layer (BBL)
!
      EFINISH  = 0.
      MWXFINISH  = 0.
      MWYFINISH  = 0.
      DO IK=1, NK
        EBAND = 0.
        A1BAND = 0.
        B1BAND = 0.
        DO ITH=1, NTH
          DIFF = SPECINIT(ITH+(IK-1)*NTH)-SPEC(ITH+(IK-1)*NTH)
          EBAND = EBAND + DIFF
          A1BAND = A1BAND + DIFF*ECOS(ITH)
          B1BAND = B1BAND + DIFF*ESIN(ITH)
          END DO
        EFINISH  = EFINISH  + EBAND * DDEN(IK) / CG1(IK)
        MWXFINISH  = MWXFINISH  + A1BAND * DDEN(IK) / CG1(IK)        &
                  * WN1(IK)/SIG(IK)
        MWYFINISH  = MWYFINISH  + B1BAND * DDEN(IK) / CG1(IK)        &
                  * WN1(IK)/SIG(IK)
        END DO
!
! Transformation in momentum flux in m^2 / s^2
!
      TAUOX=(GRAV*MWXFINISH+TAUWIX-TAUBBL(1))/DTG
      TAUOY=(GRAV*MWYFINISH+TAUWIY-TAUBBL(2))/DTG
      TAUWIX=TAUWIX/DTG
      TAUWIY=TAUWIY/DTG
      TAUWNX=TAUWNX/DTG
      TAUWNY=TAUWNY/DTG
      TAUBBL(:)=TAUBBL(:)/DTG
!
! Transformation in wave energy flux in W/m^2=kg / s^3
!
      PHIOC =DWAT*GRAV*(EFINISH+PHIAW-PHIBBL)/DTG
      PHIAW =DWAT*GRAV*PHIAW /DTG
      PHINL =DWAT*GRAV*PHINL /DTG
      PHIBBL=DWAT*GRAV*PHIBBL/DTG
!
! 9.b  Adds ice scattering and dissipation ---------------------------- *
!     INFLAGS2(4) is true if ice concentration was ever read during
!             this simulation
!
 
      IF ( INFLAGS2(4).AND.ICE.GT.0 ) THEN
 
        IF ( INFLAGS2(-6)) THEN
          ICECOEF2 = ICEP2(IX,IY)
        ELSE
          ICECOEF2 = 0.
          END IF
!
        IF (IICEDISP) THEN
          CALL LIU_FORWARD_DISPERSION (ICEH,ICECOEF2,DEPTH, &
                                       SIG,WN_R,CG_ICE,ALPHA_LIU)
         ELSE
            WN_R=WN1
            CG_ICE=CG1
            END IF
!
      SPEC2 = SPEC
!
           TAUICE(:) = 0.
           PHICE = 0.
           DO IK=1,NK
             IS = 1+(IK-1)*NTH
!
! First part of ice term integration: dissipation part ATT.
!
             ATT=1.
             SPEC(1+(IK-1)*NTH:NTH+(IK-1)*NTH) = ATT*SPEC2(1+(IK-1)*NTH:NTH+(IK-1)*NTH)
!
! Second part of ice term integration: scattering (only if there is back-scattering, i.e. IS2PARS(2).NE.0 )
!
 
             FACTOR = DDEN(IK)/CG1(IK)                    !Jacobian to get energy in band
             FACTOR2= FACTOR*GRAV*WN1(IK)/SIG(IK)         ! coefficient to get momentum
             DO ITH = 1,NTH
               IS = ITH+(IK-1)*NTH
               PHICE = PHICE + (SPEC(IS)-SPEC2(IS)) * FACTOR
               COSI(1)=ECOS(IS)
               COSI(2)=ESIN(IS)
               TAUICE(:) = TAUICE(:) - (SPEC(IS)-SPEC2(IS))*FACTOR2*COSI(:)
               END DO
             END DO
           PHICE =-1.*DWAT*GRAV*PHICE /DTG
           TAUICE(:)=TAUICE(:)/DTG
           ELSE
           END IF
!
! - - - - - - - - - - - - - - - - - - - - - -
! 9.c Sea state dependent stress routine calls
! - - - - - - - - - - - - - - - - - - - - - -
!Note the Sea-state dependent stress calculations are primarily for high-wind
!conditions (>10 m/s).  It is not recommended to use these at lower wind
!in their current state.
!
! FLD1/2 requires the calculation of FPI:
!
! 10 includes shoreline reflection --------------------------------------------- *
!
      FIRST  = .FALSE.
      IF (IT.EQ.0) SPEC = SPECINIT
!
      RETURN
!
! Formats
!
 9006 FORMAT (' TEST W3SRCE : FHIGH (3X) : ',3F8.4/                   &
              ' ------------- NEW DYNAMIC INTEGRATION LOOP',          &
              ' ------------- ')
!
 9062 FORMAT (' TEST W3SRCE : FHIGH (3X) : ',3F8.4/                   &
              '               NKH        : ',I3)
!/
!/ End of W3SRCE ----------------------------------------------------- /
!/
      END SUBROUTINE W3SRCE
!/
      SUBROUTINE CALC_FPI( A, CG, FPI, S )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |          Jessica Meixner          |
!/                  |                                   |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         6-July-2016 |
!/                  +-----------------------------------+
!/
!/    06-Jul-2016 : Origination                         ( version 5.12 )
!/
!  1. Purpose :
!
!     Calculate equivalent peak frequency
!
!  2. Method :
!
!     Tolman and Chalikov (1996), equivalent peak frequency from source
 
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D).
!       CG      R.A.  I   Group velocities for k-axis of spectrum.
!       FPI     R.A.  O   Input 'peak' frequency.
!       S       R.A.  I   Source term (1-D version).
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
!      W3SRCE Subr.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!  7. Remarks :
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!       !/S      Enable subroutine tracing.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, XFR, DDEN, SIG,FTE, FTTR
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NSPEC), CG(NK), S(NSPEC)
      REAL, INTENT(OUT)       :: FPI
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS, IK
      REAL                    ::  M0, M1, SIN1A(NK)
!/
!/ ------------------------------------------------------------------- /
!/
!
!     Calculate FPI: equivalent peak frequncy from wind source term
!     input
!
      DO IK=1, NK
        SIN1A(IK) = 0.
        DO IS=(IK-1)*NTH+1, IK*NTH
          SIN1A(IK) = SIN1A(IK) + MAX ( 0. , S(IS) )
        END DO
      END DO
!
      M0     = 0.
      M1     = 0.
      DO IK=1, NK
        SIN1A(IK) = SIN1A(IK) * DDEN(IK) / ( CG(IK) * SIG(IK)**3 )
        M0        = M0 + SIN1A(IK)
        M1        = M1 + SIN1A(IK)/SIG(IK)
      END DO
!
      SIN1A(NK) = SIN1A(NK) / DDEN(NK)
      M0        = M0 + SIN1A(NK) * FTE
      M1        = M1 + SIN1A(NK) * FTTR
      IF ( M1 .LT. 1E-20 ) THEN
          FPI    = XFR * SIG(NK)
      ELSE
          FPI    = M0 / M1
      END IF
 
      END SUBROUTINE CALC_FPI
!
!/ End of module W3SRCEMD -------------------------------------------- /
!/
      END MODULE W3SRCEMD
