!/ ------------------------------------------------------------------- /
      MODULE W3SRCEMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
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
      SUBROUTINE W3SRCE ( IX, IY, IMOD, SPEC, ALPHA, WN1, CG1, D_INP, &
                          U10ABS, U10DIR, AS, USTAR, USTDIR, CX, CY,  &
                          EMEAN, FMEAN, WNMEAN, AMAX, FPI, CD, Z0,    &
                          DTDYN, FCUT, DTG )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
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
!/
!/    Copyright 2009 National Weather Service (NWS),
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
!       DEPTH   Real.  I   Depth.
!       U10ABS  Real.  I   Wind speed at reference height.
!       U10DIR  Real.  I   Id. wind direction.
!       AS      Real.  I   Air-sea temp. difference.      ( !/ST3 )
!       USTAR   Real. !/O  Friction velocity.
!       USTDIR  Real  !/O  Idem, direction.
!       CX-Y    Real.  I   Current velocitu components.   ( !/BS1 )
!       EMEAN   Real. I/O  Mean energy.                   ( !/ST1 )
!       FMEAN   Real. I/O  Mean frequency.                ( !/ST1 )
!       WNMEAN  Real. I/O  Mean wavenumber.       ( !/ST1 , !/ST2 )
!       AMAX    Real.  O   Maximum energy.        ( !/ST1 , !/ST2 )
!       FPI     Real  I/O  Peak-input frequency.          ( !/ST2 )
!       CD      Real  I/O  Drag coefficient.              ( !/ST2 )
!       Z0      Real  I/O  Roughness length.              ( !/ST2 )
!       DTDYN   Real   O   Average dynamic time step.
!       FCUT    Real   O   Cut-off frequency for tail.
!       DTG     Real   I   Global time step.
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
!      W3SDSn    Subr. W3SRCnMD Whitecapping source term
!      W3SBTn    Subr. W3SBTnMD Bottom friction source term.
!      W3SDBn    Subr. W3SBTnMD Depth induced breaking source term.
!      W3STRn    Subr. W3STRnMD Triad interaction source term.
!      W3SBSn    Subr. W3SBSnMD Bottom scattering source term.
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
!         c Dissipation.                                    ( W3SDSn )
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
!
!     !/LN0   No linear input.                           ( Choose one )
!     !/LNX   User-defined bottom friction.
!
!     !/ST0   No input and dissipation.                  ( Choose one )
!     !/ST1   WAM-3 input and sissipation.
!     !/ST2   Tolman and Chalikov (1996)  input and sissipation.
!     !/ST3   WAM 4+ input and dissipation.
!     !/STX   User-defined input and sissipation.
!
!     !/NL0   No nonlinear interactions.                 ( Choose one )
!     !/NL1   Discrete interaction approximation.
!     !/NL2   Exact nonlinear interactions.
!     !/NLX   User-defined nonlinear interactions.
!
!     !/BT0   No bottom friction.                        ( Choose one )
!     !/BT1   JONSWAP bottom friction.
!     !/BTX   User-defined bottom friction.
!
!     !/DB0   No depth-limited breaking.                 ( Choose one )
!     !/DB1   Battjes-Janssen depth-limited breaking.
!     !/DBX   User-defined bottom friction.
!
!     !/TR0   No triad interactions.                     ( Choose one )
!     !/TRX   User-defined bottom friction.
!
!     !/BS0   No bottom scattering.                      ( Choose one )
!     !/BS1   Scattering term by F. Ardhuin.
!     !/BSX   User-defined bottom friction.
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
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, TH, DMIN, DTMAX,       &
                          DTMIN, FACTI1, FACTI2, FACSD, FACHFA, FACP, &
                          XFC, XFLT, XREL, XFT, FXFM, FXPM, DDEN,     &
                          FTE, FTF, FHMAX
      USE W3WDATMD, ONLY: TIME
      USE W3ODATMD, ONLY: NDSE, NDST
      USE W3SLN1MD
      USE W3SRC3MD
      USE W3GDATMD, ONLY : ZZWND
      USE W3SNL1MD
      USE W3SBT1MD
      USE W3SDB1MD
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IX, IY, IMOD
      REAL, INTENT(IN)        :: WN1(NK), CG1(NK), D_INP, U10ABS,     &
                                 U10DIR, AS, CX, CY, DTG
      REAL, INTENT(INOUT)     :: SPEC(NSPEC), ALPHA(NK), EMEAN,       &
                                 FMEAN, WNMEAN, USTAR, USTDIR, FPI,   &
                                 CD, Z0
      REAL, INTENT(OUT)       :: AMAX, DTDYN, FCUT
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IK, ITH, IS, IS0, NSTEPS,  NKH, NKH1,&
                                 NSPECH, IDT, IERR
      REAL                    :: DTTOT, FHIGH, DT, AFILT, DAMAX, AFAC,&
                                 HDT, ZWND, FP, DEPTH
      REAL                    :: FMEANS, FMEANWS, TAUWX, TAUWY, FH1, FH2
      REAL                    :: QCERR  = 0.     !/XNL2 and !/NNT
      REAL                    :: HM, EM
      REAL                    :: DAM (NSPEC), WN2 (NSPEC),            &
                                 VSLN(NSPEC),                         &
                                 VSIN(NSPEC), VDIN(NSPEC),            &
                                 VSNL(NSPEC), VDNL(NSPEC),            &
                                 VSDS(NSPEC), VDDS(NSPEC),            &
                                 VSBT(NSPEC), VDBT(NSPEC),            &
                                 VSDB(NSPEC), VDDB(NSPEC),            &
                                 VS  (NSPEC), VD  (NSPEC), EB(NK)
      LOGICAL                :: LLWS(NSPEC)
      REAL                    :: FOUT(NK,NTH), SOUT(NK,NTH), DOUT(NK,NTH)
      LOGICAL, SAVE           :: FLTEST = .FALSE., FLAGNN = .TRUE.
      LOGICAL                 :: SHAVE
      LOGICAL, SAVE           :: FIRST = .TRUE.
!/
!/ ------------------------------------------------------------------- /
!/
!
      DEPTH  = MAX ( DMIN , D_INP )
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
!
! 1.c Set mean parameters
!
      TAUWX=0.
      TAUWY=0.
      IF ( FIRST ) THEN
          FIRST  = .FALSE.
          LLWS(:) = .TRUE.
      ELSE
        CALL W3SPR3 (SPEC, CG1, WN1, EMEAN, FMEAN, FMEANS, WNMEAN, &
                   AMAX, U10ABS, U10DIR, USTAR, USTDIR,          &
                   TAUWX, TAUWY, CD, Z0, LLWS, FMEANWS)
        CALL W3SIN3 ( SPEC, CG1, WN2, U10ABS, USTAR, DAIR/DWAT, AS,&
                 U10DIR, Z0, CD, TAUWX, TAUWY, VSIN, VDIN, LLWS )
        END IF
      CALL W3SPR3 (SPEC, CG1, WN1, EMEAN, FMEAN, FMEANS, WNMEAN, &
                   AMAX, U10ABS, U10DIR, USTAR, USTDIR,          &
                   TAUWX, TAUWY, CD, Z0, LLWS, FMEANWS)
!
! 1.d Stresses
!
! 1.e Prepare cut-off
!
      FHIGH  = FXFM * MAX(FMEAN,FMEANWS)
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
      CALL W3SIN3 ( SPEC, CG1, WN2, U10ABS, USTAR, DAIR/DWAT, AS,&
                 U10DIR, Z0, CD, TAUWX, TAUWY, VSIN, VDIN, LLWS )
!
! 2.b Nonlinear interactions.
!
        CALL W3SNL1 ( SPEC, CG1, WNMEAN*DEPTH,          VSNL, VDNL )
!
! 2.c Dissipation.
!
        CALL W3SDS3 ( SPEC,  WN1, EMEAN, FMEANS, WNMEAN,              &
                                  USTAR, USTDIR, DEPTH, VSDS, VDDS )
!
        CALL W3SDB1 ( SPEC, WN2, DEPTH, EMEAN, FMEAN, WNMEAN,        &
                                                        VSDB, VDDB )
!
! 2.d Bottom interactions.
!
        CALL W3SBT1 ( SPEC, CG1, WN1, DEPTH,            VSBT, VDBT )
!
! 2.e Additional sources.
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
        DO IS=1, NSPECH
          VS(IS) = VSLN(IS) + VSIN(IS) + VSNL(IS) + VSDS(IS) + VSBT(IS)
          VS(IS) = VS(IS) + VSDB(IS)
          VD(IS) =            VDIN(IS) + VDNL(IS) + VDDS(IS) + VDBT(IS)
          VD(IS) = VD(IS) + VDDB(IS)
          DAMAX  = MIN ( DAM(IS) , MAX ( XREL*SPEC(IS) , AFILT ) )
          AFAC   = 1. / MAX( 1.E-10 , ABS(VS(IS)/DAMAX) )
          DT     = MIN ( DT , AFAC / ( MAX ( 1.E-10,                  &
                         1. + OFFSET*AFAC*MIN(0.,VD(IS)) ) ) )
          END DO
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
            DO IS=1, NSPECH
              VS(IS) = VS(IS) * DT / MAX ( 1. , (1.-HDT*VD(IS)))
              VS(IS) = SIGN ( MIN (DAM(IS),ABS(VS(IS))) , VS(IS) )
              SPEC(IS) = MAX ( 0. , SPEC(IS)+VS(IS) )
              END DO
!
          ELSE
!
            DO IS=1, NSPECH
              VS(IS) = VS(IS) * DT / MAX ( 1. , (1.-HDT*VD(IS)))
              SPEC(IS) = MAX ( 0. , SPEC(IS)+VS(IS) )
              END DO
!
          END IF
!
! 6.  Add tail ------------------------------------------------------- *
!   a Mean parameters
!
        CALL W3SPR3 (SPEC, CG1, WN1, EMEAN, FMEAN, FMEANS,        &
                     WNMEAN, AMAX, U10ABS, U10DIR, USTAR, USTDIR, &
                     TAUWX, TAUWY, CD, Z0, LLWS, FMEANWS)
!
        FH1    = FXFM * FMEAN
        FH2    = FXPM / USTAR
        FHIGH  = MIN ( SIG(NK) , MAX ( FH1 , FH2 ) )
        NKH    = MAX ( 2 , MIN ( NKH1 ,                           &
                 INT ( FACTI2 + FACTI1*LOG(MAX(1.E-7,FHIGH)) ) ) )
!
        IF ( FLTEST ) WRITE (NDST,9062)                           &
                      FH1*TPIINV, FH2*TPIINV, FHIGH*TPIINV, NKH
!
! 6.b Limiter for shallow water or Miche style criterion
!     Last time step only !
!
        IF ( DTTOT .GE. 0.9999*DTG ) THEN
            HM     = FHMAX / WNMEAN * TANH(WNMEAN*MAX(0.,DEPTH))
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
            SPEC(ITH+(IK-1)*NTH) = SPEC(ITH+(IK-2)*NTH) * FACHFA          &
                       + 0.
            END DO
          END DO
!
! 6.e  Update wave-supported stress----------------------------------- *
!
        CALL W3SIN3 ( SPEC, CG1, WN2, U10ABS, USTAR, DAIR/DWAT, AS,&
                 U10DIR, Z0, CD, TAUWX, TAUWY, VSIN, VDIN, LLWS )
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
!/ End of module W3SRCEMD -------------------------------------------- /
!/
      END MODULE W3SRCEMD
