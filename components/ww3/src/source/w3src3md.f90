!/ ------------------------------------------------------------------- /
      MODULE W3SRC3MD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                SHOM |
!/                  !            F. Ardhuin             !
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Oct-2007 |
!/                  +-----------------------------------+
!/
!/    09-Oct-2007 : Origination.                        ( version 3.13 )
!/
!  1. Purpose :
!
!     The 'WAM4+' source terms based on P.A.E.M. Janssen's work, with
!     extensions by him and by J.-R. Bidlot. Converted from the original
!     WAM codes by F. Ardhuin.
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
!      W3SPR3    Subr. Public   Mean parameters from spectrum.
!      W3SIN3    Subr. Public   WAM4+ input source term.
!      INSIN3    Subr. Public   Corresponding initialization routine.
!      TABU_STRESS, TABU_TAUHF, TABU_TAUHF2
!                Subr. Public   Populate various tables.
!      CALC_USTAR
!                Subr. Public   Compute stresses.
!      W3SDS3    Subr. Public   User supplied dissipation.
!      INPTAB, W3BETA           Adapted from !/ST2.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!  6. Switches :
!
!  7. Source code :
!/
!/ ------------------------------------------------------------------- /
!/
      PUBLIC
!/
!/ Public variables
!/
      INTEGER, PARAMETER      :: NHMAX =    25
      INTEGER(KIND=4)         :: NH(3), THO(2,3,NHMAX)
      REAL                    :: HA(NHMAX,3), HD(NHMAX,3), HA2(NHMAX,3)
      REAL,    PARAMETER      :: kappa = 0.40       !Von Karman's constant
      !air kinematic viscosity (used in WAM)
      REAL,    PARAMETER      :: nu_air  = 1.4E-5
      INTEGER, PARAMETER      :: ITAUMAX=200,JUMAX=200
      INTEGER, PARAMETER      :: IUSTAR=100,IALPHA=200, ILEVTAIL=50
      INTEGER, PARAMETER      :: IAB=200
      REAL                    :: TAUT(0:ITAUMAX,0:JUMAX), DELTAUW, DELU
      ! Table for H.F. stress as a function of 2 variables
      REAL                    :: TAUHFT(0:IUSTAR,0:IALPHA), DELUST, DELALP
      ! Table for H.F. stress as a function of 3 variables
      REAL                    :: TAUHFT2(0:IUSTAR,0:IALPHA,0:ILEVTAIL)
      ! Table for swell damping
      REAL                    :: SWELLFT(0:IAB)
      REAL                    :: DELTAIL
      REAL                    :: DELAB
      REAL,    PARAMETER      :: UMAX    = 50.
      REAL,    PARAMETER      :: TAUWMAX = 2.2361 !SQRT(5.)
      REAL,    PARAMETER      :: ABMIN = 0.3
      REAL,    PARAMETER      :: ABMAX = 8.
      REAL, DIMENSION(:)   , ALLOCATABLE :: XSTRESS,YSTRESS
      INTEGER, DIMENSION(:,:)   , ALLOCATABLE :: SATINDICES
!
! variables for negative wind input (beta from ST2)
!
      INTEGER, PARAMETER, PRIVATE :: NRSIGA =  400
      INTEGER, PARAMETER, PRIVATE :: NRDRAG =   20
      REAL, PARAMETER, PRIVATE    :: SIGAMX =   40.
      REAL, PARAMETER, PRIVATE    :: DRAGMX =    1.E-2
!
      REAL, PRIVATE           :: DSIGA, DDRAG,                        &
                                 BETATB(-NRSIGA:NRSIGA+1,NRDRAG+1)
!/
      LOGICAL, SAVE , PRIVATE :: FIRST = .TRUE.
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SPR3 (A, CG, WN, EMEAN, FMEAN, FMEANS, WNMEAN,     &
                    AMAX, U, UDIR, USTAR, USDIR, TAUWX, TAUWY, CD, Z0,&
		              LLWS, FMEANWS)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                SHOM |
!/                  !            F. Ardhuin             !
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Oct-2007 |
!/                  +-----------------------------------+
!/
!/    03-Oct-2007 : Origination.                        ( version 3.13 )
!/
!  1. Purpose :
!
!     Calculate mean wave parameters for the use in the source term
!     routines.
!
!  2. Method :
!
!     See source term routines.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum.
!       CG      R.A.  I   Group velocities.
!       WN      R.A.  I   Wavenumbers.
!       EMEAN   Real  O   Energy
!       FMEAN   Real  O   Mean  frequency for determination of tail
!       FMEANS  Real  O   Mean  frequency for dissipation source term
!       WNMEAN  Real  O   Mean wavenumber.
!       AMAX    Real  O   Maximum of action spectrum.
!       U       Real  I   Wind speed.
!       UDIR    Real  I   Wind direction.
!       USTAR   Real I/O  Friction velocity.
!       USDIR   Real I/O  wind stress direction.
!       TAUWX-Y Real  I   Components of wave-supported stress.
!       CD      Real  O   Drag coefficient at wind level ZWND.
!       Z0      Real  O   Corresponding z0.
!       LLWS    L.A.  I   Wind sea true/false array for each component
!       FMEANWS Real  O   Mean frequency of wind sea, used for tail
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SRCE   Source term integration routine.
!       W3OUTP   Point output program.
!       GXEXPO   GrADS point output program.
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
!       !/T      Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DTH, DDEN, WWNMEANP, &
                          WWNMEANPTAIL, FTE, FTF, SSTXFTF, SSTXFTWN,&
                          SSTXFTFTAIL, SSWELLF ,SSWELLFPAR
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NTH,NK), CG(NK), WN(NK), U, UDIR
      REAL, INTENT(IN)        :: TAUWX, TAUWY
      LOGICAL, INTENT(IN)     :: LLWS(NSPEC)
!!      REAL, INTENT(IN)        :: DI(NSPEC)  TO BE CLEANED UP
      REAL, INTENT(INOUT)     :: USTAR ,USDIR
      REAL, INTENT(OUT)       :: EMEAN, FMEAN, FMEANS, WNMEAN, AMAX,  &
                                 CD, Z0, FMEANWS
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IKP1   ! wind sea peak index
      INTEGER                 :: IS, IK, ITH, I1, ITT
 
      REAL                    :: TAUW, EBAND, EMEANWS, RDCH, FXPMC,    &
                                 WNP, UNZ,    FP,                      &
                                 R1, CP, EB(NK),EB2(NK),ALFA(NK)
!/
!/ ------------------------------------------------------------------- /
!/
!
      UNZ    = MAX ( 0.01 , U )
      USTAR  = MAX ( 0.0001 , USTAR )
!
      EMEAN  = 0.
      EMEANWS= 0.
      FMEANWS= 0.
      FMEAN  = 0.
      FMEANS = 0.
      WNMEAN = 0.
      AMAX   = 0.
!
! 1.  Integral over directions and maximum --------------------------- *
!
      DO IK=1, NK
        EB(IK)  = 0.
        EB2(IK) = 0.
        DO ITH=1, NTH
          IS=ITH+(IK-1)*NTH
          EB(IK) = EB(IK) + A(ITH,IK)
          IF (LLWS(IS)) EB2(IK) = EB2(IK) + A(ITH,IK)
          AMAX   = MAX ( AMAX , A(ITH,IK) )
          END DO
        END DO
!
! 2.  Integrate over directions -------------------------------------- *
!
      DO IK=1, NK
        ALFA(IK) = 2. * DTH * SIG(IK) * EB(IK) * WN(IK)**3
        EB(IK)   = EB(IK) * DDEN(IK) / CG(IK)
        EB2(IK)   = EB2(IK) * DDEN(IK) / CG(IK)
        EMEAN    = EMEAN  + EB(IK)
        FMEAN    = FMEAN  + EB(IK) *(SIG(IK)**(2.*WWNMEANPTAIL))
        FMEANS   = FMEANS + EB(IK) *(SIG(IK)**(2.*WWNMEANP))
        WNMEAN   = WNMEAN + EB(IK) *(WN(IK)**WWNMEANP)
        EMEANWS  = EMEANWS+ EB2(IK)
        FMEANWS  = FMEANWS+ EB2(IK)*(SIG(IK)**(2.*WWNMEANPTAIL))
        END DO
!
! 3.  Add tail beyond discrete spectrum and get mean pars ------------ *
!     ( DTH * SIG absorbed in FTxx )
!
      EBAND  = EB(NK) / DDEN(NK)
      EMEAN  = EMEAN  + EBAND * FTE
      FMEAN  = FMEAN  + EBAND * SSTXFTFTAIL
      FMEANS = FMEANS + EBAND * SSTXFTF
      WNMEAN = WNMEAN + EBAND * SSTXFTWN
      EBAND  = EB2(NK) / DDEN(NK)
      EMEANWS = EMEANWS + EBAND * FTE
      FMEANWS = FMEANWS + EBAND * SSTXFTFTAIL
!
! 4.  Final processing
!
 
      IF (FMEAN.LT.1.E-7) THEN
        FMEAN=TPIINV * SIG(NK)
      ELSE
      FMEAN  = TPIINV *( MAX ( 1.E-7 , FMEAN )                       &
                     / MAX ( 1.E-7 , EMEAN ))**(1/(2.*WWNMEANPTAIL))
        ENDIF
      IF (FMEANS.LT.1.E-7) THEN
        FMEANS=TPIINV * SIG(NK)
      ELSE
        FMEANS = TPIINV *( MAX ( 1.E-7 , FMEANS )                      &
                     / MAX ( 1.E-7 , EMEAN ))**(1/(2.*WWNMEANP))
        ENDIF
      WNMEAN = ( MAX ( 1.E-7 , WNMEAN )                              &
                / MAX ( 1.E-7 , EMEAN ) )**(1/WWNMEANP)
      IF (FMEANWS.LT.1.E-7.OR.EMEANWS.LT.1.E-7) THEN
        FMEANWS=TPIINV * SIG(NK)
      ELSE
        FMEANWS  = TPIINV *( MAX ( 1.E-7 , FMEANWS )                       &
                     / MAX ( 1.E-7 , EMEANWS ))**(1/(2.*WWNMEANPTAIL))
	ENDIF
!
! 5.  Cd and z0 ----------------------------------------------- *
!
      TAUW = SQRT(TAUWX**2+TAUWY**2)
 
      Z0=0.
      CALL CALC_USTAR(U,TAUW,USTAR,Z0)
      UNZ    = MAX ( 0.01 , U )
      CD     = (USTAR/UNZ)**2
      USDIR = UDIR
!
! 6.  Final test output ---------------------------------------------- *
!
      RETURN
!
! Formats
!
!/
!/ End of W3SPR3 ----------------------------------------------------- /
!/
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SIN3 (A, CG, K, U, USTAR, DRAT, AS, USDIR, Z0, CD,    &
                         TAUWX, TAUWY, S, D, LLWS)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III                SHOM |
!/                  !            F. Ardhuin             !
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Oct-2007 |
!/                  +-----------------------------------+
!/
!/    09-Oct-2007 : Origination.                        ( version 3.13 )
!/
!  1. Purpose :
!
!     Calculate diagonal and input source term for WAM4+ approach.
!
!  2. Method :
!
!       WAM-4 : Janssen et al.
!       WAM-"4.5" : gustiness effect (Cavaleri et al. )
!       SAT    : high-frequency input reduction for balance with
!                saturation dissipation (Ardhuin et al., work in progress)
!       SWELL: negative wind input with various parameterizations,
!              including Tolman and Chalikov 1996 (SWELLFPAR=1)
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D).
!       CG      R.A.  I   Group speed                              *)
!       K       R.A.  I   Wavenumber for entire spectrum.          *)
!       U       Real  I   WIND SPEED
!       USTAR   Real  I   Friction velocity.
!       DRAT    Real  I   Air/water density ratio.
!       AS      Real  I   Air-sea temperature difference
!       USDIR   Real  I   wind stress direction
!       Z0      Real  I   Air-side roughness lengh.
!       CD      Real  I   Wind drag coefficient.
!       USDIR   Real  I   Direction of friction velocity
!       TAUWX-Y Real  I   Components of the wave-supported stress.
!       S       R.A.  O   Source term (1-D version).
!       D       R.A.  O   Diagonal term of derivative.             *)
!       LLWS    L.A.  O   Wind sea true/false array for each component
!     ----------------------------------------------------------------
!                         *) Stored as 1-D array with dimension NTH*NK
!
!  4. Subroutines used :
!
!       STRACE    Subroutine tracing.                 ( !/S switch )
!       PRT2DS    Print plot of spectrum.             ( !/T0 switch )
!       OUTMAT    Print out matrix.                   ( !/T1 switch )
!
!  5. Called by :
!
!       W3SRCE   Source term integration.
!       W3EXPO   Point output program.
!       GXEXPO   GrADS point output program.
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
!     !/S   Enable subroutine tracing.
!     !/T   Enable general test output.
!     !/T0  2-D print plot of source term.
!     !/T1  Print arrays.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: NK, NTH, NSPEC, XFR, DDEN, SIG, SIG2, TH,   &
                          ESIN, ECOS, EC2, ZZWND, AALPHA, BBETA, ZZALP,&
                          TTAUWSHELTER, SSWELLF, SSWELLFPAR, &
                          DDEN2, DTH, SSINTHP,ZZ0RAT
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NSPEC), CG(NK), K(NSPEC),Z0,U, CD
      REAL, INTENT(IN)        :: USTAR, USDIR, AS, DRAT
      REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC), TAUWX, TAUWY
      LOGICAL, INTENT(OUT)    :: LLWS(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS,IK,ITH, IOMA, ICL
      REAL                    :: FACLN1, FACLN2, ULAM, CLAM, OMA, &
                                 RD1, RD2, LAMBDA, COSFAC, FACPAR
      REAL                    :: COSU, SINU, TAUX, TAUY, USDIRP, USTP
      REAL                    :: TAUPX, TAUPY, UST2, TAUW, TAUWB
      REAL   , PARAMETER      :: EPS1 = 0.00001, EPS2 = 0.000001
      REAL Usigma       !standard deviation of U due to gustiness
      REAL USTARsigma       !standard deviation of USTAR due to gustiness
      REAL                    :: BETA, mu_janssen, omega_janssen,     &
                                 CM,ZCO,UCO,UCN,ZCN, &
                                 Z0VISC, Z0NOZ,  UORBX, UORBY, EB,  &
				 EBX, EBY, AORB, FW, UORB, M2, TH2, &
                                 UORBT, REYNOLDS, FU, FUD
      REAL XI,DELI1,DELI2
      REAL XJ,DELJ1,DELJ2
      REAL XK,DELK1,DELK2
      REAL CONST,CONST0,CONST2,TAU1
      REAL X,ZARG,ZLOG,UST
      REAL COSWIND,XSTRESS,YSTRESS,TAUHF
      REAL TEMP, TEMP2
      INTEGER IND,J,I,ISTAB
      REAL DSTAB(3,NSPEC)
      REAL STRESSSTAB(3,2),STRESSSTABN(3,2)
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Preparations
!
      Z0VISC = 0.1*nu_air/MAX(USTAR,0.0001)
      Z0NOZ = MAX(Z0VISC,ZZ0RAT*Z0)
      FACLN1 = U / LOG(ZZWND/Z0NOZ)
      FACLN2 = LOG(Z0NOZ)
      IF (SSWELLFPAR.EQ.1) THEN
        FACPAR=1
      ELSE
        FACPAR=0.5
        END IF
! 2.  Diagonal
!
! Here AS is the air-sea temperature difference in degrees. Expression given by
! Abdalla & Cavaleri, JGR 2002 for Usigma. For USTARsigma ... I do not see where
! I got it from, maybe just made up from drag law ...
      Usigma=MAX(0.,-0.025*AS)
      USTARsigma=(1.0+U/(10.+U))*Usigma
      UST=USTAR
      ISTAB=3
      DO ISTAB=1,2
      IF (ISTAB.EQ.1) UST=USTAR*(1.-USTARsigma)
      IF (ISTAB.EQ.2) UST=USTAR*(1.+USTARsigma)
      TAUX = UST**2* COS(USDIR)
      TAUY = UST**2* SIN(USDIR)
 
      IF (SSWELLFPAR.GT.2) THEN
        UORBX=0.
        UORBY=0.
        UORBT=0.
        AORB=0.
        DO IK=1, NK
          EB  = 0.
          EBX = 0.
          EBY = 0.
          DO ITH=1, NTH
             IS=ITH+(IK-1)*NTH
             EB  = EB  + A(IS)
	   !  EBX = EBX + A(IS)*(2*EC2(IS)-1)
           !  EBY = EBY + A(IS)*2*ESIN(IS)*ECOS(IS)
             END DO
	  UORBT=UORBT+EB *SIG(IK)**2 * DDEN(IK) / CG(IK)
	  !UORBX=UORBX+EBX*SIG(IK)**2 * DDEN(IK) / CG(IK)
          !UORBY=UORBY+EBY*SIG(IK)**2 * DDEN(IK) / CG(IK)
          AORB= AORB + EB            * DDEN(IK) / CG(IK)  !deep water only
          END DO
        !M2=SQRT(UORBX**2+UORBY**2)
          UORBT=2*SQRT(UORBT)  ! this is the significant orbital amplitude
	  REYNOLDS = 8*UORBT*AORB / NU_AIR ! this is the Reynolds number
        !M2=2*SQRT(M2)
	!TH2=ATAN2(UORBY,UORBX)/2.
          IF (SSWELLF(2).EQ.0) THEN
	     FW=MAX(ABS(SSWELLF(3)),0.)
	     FU=0.
	     FUD=0.
	  ELSE
	    FU=ABS(SSWELLF(3))
	    FUD=SSWELLF(2)
	    AORB=2*SQRT(AORB)
            XI=(ALOG10(MAX(AORB/Z0NOZ,3.)) &
                                     -ABMIN)/DELAB
            IND  = MIN (IAB-1, INT(XI))
            DELI1= MIN (1. ,XI-FLOAT(IND))
            DELI2= 1. - DELI1
            FW =SWELLFT(IND)*DELI2+SWELLFT(IND+1)*DELI1
            END IF
          UORB=UORBT
       END IF
 !
 ! Loop over the resolved part of the spectrum
 !
      STRESSSTAB(ISTAB,:)=0.
      STRESSSTABN(ISTAB,:)=0.
      DO IK=1, NK
        TAUPX=TAUX-ABS(TTAUWSHELTER)*STRESSSTAB(ISTAB,1)
        TAUPY=TAUY-ABS(TTAUWSHELTER)*STRESSSTAB(ISTAB,2)
        USTP=(TAUPX**2+TAUPY**2)**0.25
        USDIRP=ATAN2(TAUPY,TAUPX)
        COSU   = COS(USDIRP)
        SINU   = SIN(USDIRP)
        IS=1+(IK-1)*NTH
        CM=K(IS)/SIG2(IS) !inverse of phase speed
        UCN=USTP*CM+ZZALP  !this is the inverse wave age
           ! the stress is the real stress (N/m^2) divided by
           ! rho_a, and thus comparable to USTAR**2
           ! it is the integral of rho_w g Sin/C /rho_a
           ! (air-> waves momentum flux)
        CONST2=DDEN2(IS)/CG(IK) &        !Jacobian to get energy in band
              *GRAV/(SIG(IK)/K(IS)*DRAT) ! coefficient to get momentum
        CONST=SIG2(IS)*BBETA*DRAT/(kappa**2)
           ! this CM parameter is 1 / C_phi
           ! this is the "correct" shallow-water expression
           ! here Z0 corresponds to Z0+Z1 of the Janssen eq. 14
        ZCN=ALOG(K(IS)*Z0)
           ! this was the origina WAM version (OK for deep water)  g*z0/C^2
           ! ZCN=ALOG(G*Z0b(I)*CM(I)**2)
        DO ITH=1,NTH
          IS=ITH+(IK-1)*NTH
          COSWIND=(ECOS(IS)*COSU+ESIN(IS)*SINU)
          IF (COSWIND.GT.0.01) THEN
            X=COSWIND*UCN
            ! this ZARG term is the argument of the exponential
            ! in Janssen 1991 eq. 16.
            ZARG=KAPPA/X
            ! ZLOG is ALOG(MU) where MU is defined by Janssen 1991 eq. 15
            ! MU=
            ZLOG=ZCN+ZARG
 
            IF (ZLOG.LT.0.) THEN
              ! The source term Sp is beta * omega * X**2
              ! as given by Janssen 1991 eq. 19
              DSTAB(ISTAB,IS) = CONST*EXP(ZLOG)*ZLOG**4*UCN**2*COSWIND**SSINTHP
	      LLWS(IS)=.TRUE.
            ELSE
              DSTAB(ISTAB,IS) = 0.
	      LLWS(IS)=.FALSE.
              END IF
!
!  Added for consistency with ECWAM implsch.F
!
	    IF (28.*CM*USTAR*COSWIND.GE.1) THEN
	       LLWS(IS)=.TRUE.
	       ENDIF
          ELSE
            DSTAB(ISTAB,IS) = 0.
	    LLWS(IS)=.FALSE.
            END IF
          IF ((SSWELLF(1).NE.0.AND.DSTAB(ISTAB,IS).LT.1E-7*SIG2(IS)) &
              .OR.SSWELLF(3).GT.0) THEN
            IF (SSWELLFPAR.LE.2) THEN
              ! This is done if negative input is activated
              ! First test: adaptation of Tolman and Chalikov
 
              COSFAC = ECOS(IS)*COSU + ESIN(IS)*SINU
              COSFAC = SIGN ( MAX ( 0.0087 , ABS(COSWIND) ) , COSWIND )
              LAMBDA = FACPAR*TPI / ( K(IS) * ABS(COSFAC) )
              ULAM   = FACLN1 * ( LOG(LAMBDA) - FACLN2 )
              CLAM   = CD * ( U*USTAR/MAX(USTAR,0.0001)/ ULAM )**2
              OMA    = K(IS) * ULAM * COSFAC / SIG2(IS)
              IOMA   = INT ( OMA/DSIGA ) +                              &
                    MIN ( 0 , INT ( SIGN ( -1.1 , OMA ) ) )
              ICL    = INT ( CLAM/DDRAG )
              RD1    = OMA/DSIGA - REAL(IOMA)
              RD2    = CLAM/DDRAG - REAL(ICL)
              IOMA   = MAX ( -NRSIGA , MIN ( NRSIGA , IOMA ) )
              ICL    = MAX ( 1 , MIN ( NRDRAG , ICL ) )
              BETA   = (1.-RD1) * (1.-RD2) * BETATB( IOMA , ICL )       &
                   +    RD1   * (1.-RD2) * BETATB(IOMA+1, ICL )       &
               + (1.-RD1) *    RD2   * BETATB( IOMA ,ICL+1)           &
               +    RD1   *    RD2   * BETATB(IOMA+1,ICL+1)
              DSTAB(ISTAB,IS)  = BETA * DRAT * SIG2(IS)
              DSTAB(ISTAB,IS)  = MIN ( 0. , SSWELLF(1)*DSTAB(ISTAB,IS) )   !only for negative values
   	    ELSE
	      IF (REYNOLDS.LT.SSWELLF(4)) THEN
	        DSTAB(ISTAB,IS) =  DSTAB(ISTAB,IS)                    &
		 		-SSWELLF(5)*DRAT*2*K(IS)*SQRT(2*NU_AIR*SIG2(IS))
	      ELSE
	        DSTAB(ISTAB,IS) =  DSTAB(ISTAB,IS) &
                                -DRAT*SSWELLF(1)*                         &
                                (FW*UORB+(FU+FUD*COSWIND)*USTP) &
                                           *16*SIG2(IS)**2/GRAV
                END IF					
              END IF
          ! Wave direction is "direction to"
          ! therefore there is a PLUS sign for the stress
	    ENDIF
          TEMP2=CONST2*DSTAB(ISTAB,IS)*A(IS)
          IF (DSTAB(ISTAB,IS).GE.0) THEN
            STRESSSTAB(ISTAB,1)=STRESSSTAB(ISTAB,1)+TEMP2*ECOS(IS)
            STRESSSTAB(ISTAB,2)=STRESSSTAB(ISTAB,2)+TEMP2*ESIN(IS)
            END IF
          END DO
        END DO
!
        D(:)=DSTAB(3,:)
        XSTRESS=STRESSSTAB(3,1)
        YSTRESS=STRESSSTAB(3,2)
      END DO
      D(:)=0.5*(DSTAB(1,:)+DSTAB(2,:))
      XSTRESS=0.5*(STRESSSTAB(1,1)+STRESSSTAB(2,1))
      YSTRESS=0.5*(STRESSSTAB(1,2)+STRESSSTAB(2,2))
      S = D * A
!
! ... Test output of arrays
!
      ! Computes the high-frequency contribution
      ! the difference in spectal density (kx,ky) to (f,theta)
      ! is integrated in this modified CONST0
      CONST0=DTH*SIG(NK)**5/((GRAV**2)*tpi) &
         *TPI*SIG(NK) / CG(NK)  !conversion WAM (E(f,theta) to WW3 A(k,theta)
      TEMP=0.
      DO ITH=1,NTH
         IS=ITH+(NK-1)*NTH
         COSWIND=(ECOS(IS)*COSU+ESIN(IS)*SINU)
         TEMP=TEMP+A(IS)*(MAX(COSWIND,0.))**3
         END DO
 
      TAUPX=TAUX-ABS(TTAUWSHELTER)*XSTRESS
      TAUPY=TAUY-ABS(TTAUWSHELTER)*YSTRESS
      USTP=(TAUPX**2+TAUPY**2)**0.25
      USDIRP=ATAN2(TAUPY,TAUPX)
 
      UST=USTP
      ! finds the values in the tabulated stress TAUHFT
      XI=UST/DELUST
      IND  = MIN (IUSTAR-1, INT(XI))
      DELI1= MIN (1. ,XI-FLOAT(IND))
      DELI2= 1. - DELI1
      XJ=MAX(0.,(GRAV*Z0/UST**2-AALPHA) / DELALP)
      J    = MIN (IALPHA-1, INT(XJ))
      DELJ1= MIN (1. ,XJ-FLOAT(J))
      DELJ2=1. - DELJ1
      IF (TTAUWSHELTER.GT.0) THEN
        XK = CONST0*TEMP / DELTAIL
         I = MIN (ILEVTAIL-1, INT(XK))
         DELK1= MIN (1. ,XK-FLOAT(I))
         DELK2=1. - DELK1
         TAU1 =((TAUHFT2(IND,J,I)*DELI2+TAUHFT2(IND+1,J,I)*DELI1 )*DELJ2 &
               +(TAUHFT2(IND,J+1,I)*DELI2+TAUHFT2(IND+1,J+1,I)*DELI1)*DELJ1)*DELK2 &
              +((TAUHFT2(IND,J,I+1)*DELI2+TAUHFT2(IND+1,J,I+1)*DELI1 )*DELJ2 &
               +(TAUHFT2(IND,J+1,I+1)*DELI2+TAUHFT2(IND+1,J+1,I+1)*DELI1)*DELJ1)*DELK1
      ELSE
        TAU1 =(TAUHFT(IND,J)*DELI2+TAUHFT(IND+1,J)*DELI1 )*DELJ2 &
         +(TAUHFT(IND,J+1)*DELI2+TAUHFT(IND+1,J+1)*DELI1)*DELJ1
        END IF
      TAUHF = CONST0*TEMP*UST**2*TAU1
      TAUWX = XSTRESS+TAUHF*COS(USDIRP)
      TAUWY = YSTRESS+TAUHF*SIN(USDIRP)
!
! Reduces tail effect to make sure that wave-supported stress
! is less than total stress, this is borrowed from ECWAM Stresso.F
!
      TAUW = SQRT(TAUWX**2+TAUWY**2)
      UST2   = MAX(USTAR,EPS2)**2
      TAUWB = MIN(TAUW,MAX(UST2-EPS1,EPS2**2))
      IF (TAUWB.LT.TAUW) THEN
         TAUWX=TAUWX*TAUWB/TAUW
         TAUWY=TAUWY*TAUWB/TAUW
	 END IF
	
      RETURN
!
! Formats
!
!/
!/ End of W3SIN3 ----------------------------------------------------- /
!/
      END SUBROUTINE
!/ ------------------------------------------------------------------- /
      SUBROUTINE INSIN3
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         21-Aug-2006 |
!/                  +-----------------------------------+
!/
!/    23-Jun-2006 : Origination.                        ( version 3.09 )
!
!  1. Purpose :
!
!     Initialization for source term routine.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
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
!      W3SIN3    Subr. W3SRC3MD Corresponding source term.
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
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS, ONLY: TPIINV, RADE
      USE W3ODATMD,  ONLY: NDSE
      USE W3SERVMD,  ONLY: EXTCDE
      USE W3GDATMD,  ONLY: SIG, NK, NTH, TTAUWSHELTER, SSWELLFPAR, &
                           SSDSDTH, DTH
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
    INTEGER  SDSNTH, ITH, I_INT
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  .... ----------------------------------------------------------- *
!
 
      CALL TABU_STRESS
      CALL TABU_TAUHF(SIG(NK) * TPIINV)   !tabulate high-frequency stress
      IF (SSWELLFPAR.EQ.3) CALL TABU_SWELLFT
      IF (TTAUWSHELTER.GT.0) THEN
        CALL TABU_TAUHF2(SIG(NK) * TPIINV)   !tabulate high-frequency stress
      !ELSE
      !  CALL TABU_TAUHF(SIG(NK) * TPIINV)   !tabulate high-frequency stress
        END IF
!
! Precomputes the indices for integrating the spectrum to get saturation
!
      IF (SSDSDTH.NE.180) THEN
	SDSNTH  = MIN(NINT(SSDSDTH/(DTH*RADE)),NTH/2-1)
        ALLOCATE(SATINDICES(NTH,SDSNTH*2+1))
        DO ITH=1,NTH
          DO I_INT=ITH-SDSNTH, ITH+SDSNTH
            SATINDICES(ITH,I_INT-(ITH-SDSNTH)+1)=I_INT
	    IF (I_INT.LT.1)  &
	    	SATINDICES(ITH,I_INT-(ITH-SDSNTH)+1)=I_INT+NTH
	    IF (I_INT.GT.NTH) &
	    	SATINDICES(ITH,I_INT-(ITH-SDSNTH)+1)=I_INT-NTH
            END DO
	  END DO
	END IF
!/
!/ End of INSIN3 ----------------------------------------------------- /
!/
      END SUBROUTINE INSIN3
! ----------------------------------------------------------------------
      SUBROUTINE TABU_STRESS
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Oct-2007 |
!/                  +-----------------------------------+
!/
!/    23-Jun-2006 : Origination.                        ( version 3.13 )
!/     adapted from WAM, original:P.A.E.M. JANSSEN    KNMI AUGUST 1990
!/     adapted version (subr. STRESS): J. BIDLOT    ECMWF OCTOBER 2004
!/     Table values were checkes against the original f90 result and found to
!/     be identical (at least at 0.001 m/s accuracy)
!/
!  1. Purpose :
!     TO GENERATE friction velocity table TAUT(TAUW,U10)=SQRT(TAU).
!     METHOD.
!       A STEADY STATE WIND PROFILE IS ASSUMED.
!       THE WIND STRESS IS COMPUTED USING THE ROUGHNESSLENGTH
!                  Z1=Z0/SQRT(1-TAUW/TAU)
!       WHERE Z0 IS THE CHARNOCK RELATION , TAUW IS THE WAVE-
!       INDUCED STRESS AND TAU IS THE TOTAL STRESS.
!       WE SEARCH FOR STEADY-STATE SOLUTIONS FOR WHICH TAUW/TAU < 1.
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
!
!     Initialization for source term routine.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
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
!      W3SIN3    Subr. W3SRC3MD Corresponding source term.
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
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: ZZWND, AALPHA, ZZ0MAX
      IMPLICIT NONE
      INTEGER, PARAMETER      :: NITER=10
      REAL   , PARAMETER      :: XM=0.50, EPS1=0.00001
!     VARIABLE.   TYPE.     PURPOSE.
!      *XM*        REAL      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
!      *XNU*       REAL      KINEMATIC VISCOSITY OF AIR.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
!      *EPS1*      REAL      SMALL NUMBER TO MAKE SURE THAT A SOLUTION
!                            IS OBTAINED IN ITERATION WITH TAU>TAUW.
! ----------------------------------------------------------------------
      INTEGER I,J,ITER
      REAL ZTAUW,UTOP,CDRAG,WCD,USTOLD,TAUOLD
      REAL X,UST,ZZ0,ZNU,F,DELF,ZZ00
!
      DELU    = UMAX/FLOAT(JUMAX)
      DELTAUW = TAUWMAX/FLOAT(ITAUMAX)
      DO I=0,ITAUMAX
         ZTAUW   = (REAL(I)*DELTAUW)**2
         DO J=0,JUMAX
            UTOP    = FLOAT(J)*DELU
            CDRAG   = 0.0012875
            WCD     = SQRT(CDRAG)
            USTOLD  = UTOP*WCD
            TAUOLD  = MAX(USTOLD**2, ZTAUW+EPS1)
            DO ITER=1,NITER
               X   = ZTAUW/TAUOLD
               UST = SQRT(TAUOLD)
               ZZ00=AALPHA*TAUOLD/GRAV
	       IF (ZZ0MAX.NE.0) ZZ00=MIN(ZZ00,ZZ0MAX)
               ! Corrects roughness ZZ00 for quasi-linear effect
	       ZZ0 = ZZ00/(1.-X)**XM
	       !ZNU = 0.1*nu_air/UST  ! This was removed by Bidlot in 1996
               !ZZ0 = MAX(ZNU,ZZ0)
               F   = UST-KAPPA*UTOP/(ALOG(ZZWND/ZZ0))
               DELF= 1.-KAPPA*UTOP/(ALOG(ZZWND/ZZ0))**2*2./UST &
                        *(1.-(XM+1)*X)/(1.-X)
               UST = UST-F/DELF
               TAUOLD= MAX(UST**2., ZTAUW+EPS1)
               END DO
            TAUT(I,J)  = SQRT(TAUOLD)
            END DO
         END DO
	 I=ITAUMAX
	 J=JUMAX
!
!  Force zero wind to have zero stress (Bidlot 1996)
!
      DO I=0,ITAUMAX
        TAUT(I,0)=0.0
        END DO
!      DO I=0,ITAUMAX
!         DO J=0,JUMAX
!   WRITE(101,'(A,2I4,G16.8)') 'I,J,TAUT(I,J):',I,J,TAUT(I,J)
!         END DO
!        END DO
      RETURN
      END SUBROUTINE TABU_STRESS
!/ ------------------------------------------------------------------- /
      SUBROUTINE TABU_TAUHF(FRMAX)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update 2006/08/14            |
!/                  +-----------------------------------+
!/
!/    27-Feb-2004 : Origination in WW3                  ( version 2.22.SHOM )
!/     the resulting table was checked to be identical to the original f77 result
!/    14-Aug-2006 : Modified following Bidlot           ( version 2.22.SHOM )
!/    18-Aug-2006 : Ported to version 3.09
!
!  1. Purpose :
!
!     Tabulation of the high-frequency wave-supported stress
!
!  2. Method :
!
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
!     See tech. Memo ECMWF 03 december 2003 by Bidlot & Janssen
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       FRMAX   Real  I   maximum frequency.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SIN3   Wind input Source term routine.
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
!       !/T      Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: AALPHA, BBETA, ZZALP, XFR, FACHFE, ZZ0MAX
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, intent(in) :: FRMAX  !  maximum frequency
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!       USTARM  R.A.  Maximum friction velocity
!       ALPHAM  R.A.  Maximum Charnock Coefficient
!       WLV     R.A.  Water levels.
!       UA      R.A.  Absolute wind speeds.
!       UD      R.A.  Absolute wind direction.
!       U10     R.A.  Wind speed used.
!       U10D    R.A.  Wind direction used.
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      REAL                    :: USTARM, ALPHAM
      REAL                    :: CONST1, OMEGA, OMEGAC
      REAL                    :: UST, ZZ0,OMEGACC, CM
      INTEGER, PARAMETER      :: JTOT=250
      REAL, ALLOCATABLE       :: W(:)
      REAL                    :: ZX,ZARG,ZMU,ZLOG,ZZ00,ZBETA
      REAL                    :: Y,YC,DELY
      INTEGER                 :: I,J,K,L
      REAL                    :: X0
!
      USTARM = 5.
      ALPHAM = 20.*AALPHA
      DELUST = USTARM/REAL(IUSTAR)
      DELALP = ALPHAM/REAL(IALPHA)
      CONST1 = BBETA/KAPPA**2
      OMEGAC = TPI*FRMAX
!
      TAUHFT(0:IUSTAR,0:IALPHA)=0. !table initialization
!
      ALLOCATE(W(JTOT))
      W(2:JTOT-1)=1.
      W(1)=0.5
      W(JTOT)=0.5
      X0 = 0.05
!
      DO L=0,IALPHA
         DO K=0,IUSTAR
            UST      = MAX(REAL(K)*DELUST,0.000001)
            ZZ00       = UST**2*AALPHA/GRAV
            IF (ZZ0MAX.NE.0) ZZ00=MIN(ZZ00,ZZ0MAX)
	    ZZ0       = ZZ00*(1+FLOAT(L)*DELALP/AALPHA)
            OMEGACC  = MAX(OMEGAC,X0*GRAV/UST)
            YC       = OMEGACC*SQRT(ZZ0/GRAV)
            DELY     = MAX((1.-YC)/REAL(JTOT),0.)
            ! For a given value of UST and ALPHA,
            ! the wave-supported stress is integrated all the way
            ! to 0.05*g/UST
            DO J=1,JTOT
               Y        = YC+REAL(J-1)*DELY
               OMEGA    = Y*SQRT(GRAV/ZZ0)
               ! This is the deep water phase speed
               CM       = GRAV/OMEGA
               !this is the inverse wave age, shifted by ZZALP (tuning)
               ZX       = UST/CM +ZZALP
               ZARG     = MIN(KAPPA/ZX,20.)
               ZMU      = MIN(GRAV*ZZ0/CM**2*EXP(ZARG),1.)
               ZLOG     = MIN(ALOG(ZMU),0.)
               ZBETA        = CONST1*ZMU*ZLOG**4
               ! Power of Y in denominator should be FACHFE-4
               TAUHFT(K,L)  = TAUHFT(K,L)+W(J)*ZBETA/Y*DELY
               END DO
               !IF (MOD(K,5).EQ.0.AND.MOD(L,5).EQ.0) &
               !WRITE(102,'(2I4,3G16.8)') L,K,UST,AALPHA+FLOAT(L)*DELALP,TAUHFT(K,L)
         END DO
      END DO
      DEALLOCATE(W)
!      DO L=0,IALPHA
!         DO K=0,IUSTAR
!WRITE(101,'(A,2I4,G16.8)') 'L,K,TAUHFT(K,L):',L,K,TAUHFT(K,L)
!         END DO
!        END DO
!WRITE(101,*) 'TAUHFT:',FRMAX,BBETA,AALPHA,CONST1,OMEGAC,TPI
!WRITE(101,'(20G16.8)') TAUHFT
      RETURN
      END SUBROUTINE TABU_TAUHF
 
!/ ------------------------------------------------------------------- /
      SUBROUTINE TABU_TAUHF2(FRMAX)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update 2006/08/14            |
!/                  +-----------------------------------+
!/
!/    15-May-2007 : Origination in WW3                  ( version 3.10.SHOM )
!
!  1. Purpose :
!
!     Tabulation of the high-frequency wave-supported stress as a function of
!     ustar, alpha (modified Charnock), and tail energy level
!
!  2. Method :
!
!       SEE REFERENCE FOR WAVE STRESS CALCULATION.
!       FOR QUASILINEAR EFFECT SEE PETER A.E.M. JANSSEN,1990.
!     See tech. Memo ECMWF 03 december 2003 by Bidlot & Janssen
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       FRMAX   Real  I   maximum frequency.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SIN3   Wind input Source term routine.
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
!       !/T      Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: AALPHA, BBETA, ZZALP, XFR, FACHFE,  &
                          TTAUWSHELTER, ZZ0MAX
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, intent(in) :: FRMAX  !  maximum frequency
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!       USTARM  R.A.  Maximum friction velocity
!       ALPHAM  R.A.  Maximum Charnock Coefficient
!       WLV     R.A.  Water levels.
!       UA      R.A.  Absolute wind speeds.
!       UD      R.A.  Absolute wind direction.
!       U10     R.A.  Wind speed used.
!       U10D    R.A.  Wind direction used.
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      REAL                    :: USTARM, ALPHAM, LEVTAILM
      REAL                    :: CONST1, OMEGA, OMEGAC, LEVTAIL
      REAL                    :: UST, UST0, ZZ0,OMEGACC, CM
      REAL                    :: TAUW, TAUW0
      INTEGER, PARAMETER      :: JTOT=250
      REAL, ALLOCATABLE       :: W(:)
      REAL                    :: ZX,ZARG,ZMU,ZLOG,ZBETA
      REAL                    :: Y,YC,DELY
      INTEGER                 :: I, J, K, L
      REAL                    :: X0
!
      USTARM = 5.
      ALPHAM = 20.*AALPHA
      LEVTAILM = 0.05
      DELUST  = USTARM/REAL(IUSTAR)
      DELALP  = ALPHAM/REAL(IALPHA)
      DELTAIL = ALPHAM/REAL(ILEVTAIL)
      CONST1  = BBETA/KAPPA**2
      OMEGAC  = TPI*FRMAX
!
      TAUHFT(0:IUSTAR,0:IALPHA)=0.  !table initialization
!
      ALLOCATE(W(JTOT))
      W(2:JTOT-1)=1.
      W(1)=0.5
      W(JTOT)=0.5
      X0 = 0.05
!
      DO L=0,IALPHA
        DO K=0,IUSTAR
          UST0      = MAX(REAL(K)*DELUST,0.000001)
          UST=UST0
          ZZ0       = UST0**2*(AALPHA+FLOAT(L)*DELALP)/GRAV
          OMEGACC  = MAX(OMEGAC,X0*GRAV/UST)
          YC       = OMEGACC*SQRT(ZZ0/GRAV)
          DELY     = MAX((1.-YC)/REAL(JTOT),0.)
          ! For a given value of UST and ALPHA,
          ! the wave-supported stress is integrated all the way
          ! to 0.05*g/UST
          DO I=0,ILEVTAIL
            LEVTAIL=REAL(I)*DELTAIL
            TAUHFT(K,L)=0.
            TAUHFT2(K,L,I)=0.
            TAUW0=UST0**2
            TAUW=TAUW0
            DO J=1,JTOT
               Y        = YC+REAL(J-1)*DELY
               OMEGA    = Y*SQRT(GRAV/ZZ0)
               ! This is the deep water phase speed
               CM       = GRAV/OMEGA
               !this is the inverse wave age, shifted by ZZALP (tuning)
               ZX       = UST0/CM +ZZALP
               ZARG     = MIN(KAPPA/ZX,20.)
               ZMU      = MIN(GRAV*ZZ0/CM**2*EXP(ZARG),1.)
               ZLOG     = MIN(ALOG(ZMU),0.)
               ZBETA        = CONST1*ZMU*ZLOG**4
               ! Power of Y in denominator should be FACHFE-4
               TAUHFT(K,L)  = TAUHFT(K,L)+W(J)*ZBETA/Y*DELY
               ZX       = UST/CM +ZZALP
               ZARG     = MIN(KAPPA/ZX,20.)
               ZMU      = MIN(GRAV*ZZ0/CM**2*EXP(ZARG),1.)
               ZLOG     = MIN(ALOG(ZMU),0.)
               ZBETA        = CONST1*ZMU*ZLOG**4
               ! Power of Y in denominator should be FACHFE-4
               TAUHFT2(K,L,I)  = TAUHFT2(K,L,I)+W(J)*ZBETA*(UST/UST0)**2/Y*DELY
               TAUW=TAUW-W(J)*UST**2*ZBETA*LEVTAIL/Y*DELY
               UST=SQRT(MAX(TAUW,0.))
               END DO
              ! IF (MOD(K,5).EQ.0.AND.MOD(L,5).EQ.0) &
              ! WRITE(101,'(3I4,8G16.8)') L,K,I,UST0,AALPHA+FLOAT(L)*DELALP,TAUW0,UST, &
              ! TAUHFT(K,L),TAUHFT2(K,L,I),TAUHFT(K,L)*LEVTAIL*UST0**2,TAUHFT2(K,L,I)*LEVTAIL*UST0**2
             END DO
           END DO
        END DO
      DEALLOCATE(W)
      RETURN
      END SUBROUTINE TABU_TAUHF2
! ----------------------------------------------------------------------
      SUBROUTINE TABU_SWELLFT
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         17-Oct-2007 |
!/                  +-----------------------------------+
!/
!/    19-Oct-2007 : Origination.                        ( version 3.13 )
!/
!  1. Purpose :
!     TO estimate friction coefficients in oscillatory boundary layers
!     METHOD.
!      tabulation on Kelvin functions
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
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
!      W3SIN3    Subr. W3SRC3MD Corresponding source term.
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
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      IMPLICIT NONE
      INTEGER, PARAMETER      :: NITER=100
      REAL   , PARAMETER      :: XM=0.50, EPS1=0.00001
!     VARIABLE.   TYPE.     PURPOSE.
!      *XM*        REAL      POWER OF TAUW/TAU IN ROUGHNESS LENGTH.
!      *XNU*       REAL      KINEMATIC VISCOSITY OF AIR.
!      *NITER*     INTEGER   NUMBER OF ITERATIONS TO OBTAIN TOTAL STRESS
!      *EPS1*      REAL      SMALL NUMBER TO MAKE SURE THAT A SOLUTION
!                            IS OBTAINED IN ITERATION WITH TAU>TAUW.
! ----------------------------------------------------------------------
      INTEGER I,ITER
      REAL KER, KEI
      REAL ABR,ABRLOG,L10,FACT,FSUBW,FSUBWMEMO,dzeta0,dzeta0memo
!
      DELAB   = (ABMAX-ABMIN)/REAL(IAB)
      L10=ALOG(10.)
      DO I=0,IAB
         ABRLOG=ABMIN+REAL(I)*DELAB
	      ABR=EXP(ABRLOG*L10)
	      FACT=1/ABR/(21.2*KAPPA)
         FSUBW=0.05
         DO ITER=1,NITER
            fsubwmemo=fsubw
            dzeta0memo=dzeta0
            dzeta0=fact*fsubw**(-.5)
            CALL KERKEI(2.*SQRT(dzeta0),ker,kei)
            fsubw=.08/(ker**2+kei**2)
            fsubw=.5*(fsubwmemo+fsubw)
            dzeta0=.5*(dzeta0memo+dzeta0)
            END DO
            SWELLFT(I)  = fsubw
	 !WRITE(994,*) I,ABR,fsubw
         END DO
      RETURN
      END SUBROUTINE TABU_SWELLFT
 
!/ ------------------------------------------------------------------- /
      SUBROUTINE CALC_USTAR(WINDSPEED,TAUW,USTAR,Z0)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update 2006/08/14            |
!/                  +-----------------------------------+
!/
!/    27-Feb-2004 : Origination in WW3                  ( version 2.22-SHOM )
!/     the resulting table was checked to be identical to the original f77 result
!/    14-Aug-2006 : Modified following Bidlot           ( version 2.22-SHOM )
!/    18-Aug-2006 : Ported to version 3.09
!
!  1. Purpose :
!
!     Compute friction velocity based on wind speed U10
!
!  2. Method :
!
!     Computation of u* based on Quasi-linear theory
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       U10,TAUW,USTAR,Z0
!     ----------------------------------------------------------------
!       WINDSPEED Real  I   10-m wind speed ... should be NEUTRAL
!       TAUW      Real  I   Wave-supported stress
!       USTAR     Real  O   Friction velocity.
!       Z0        Real  O   air-side roughness length
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       STRACE   Service routine.
!
!  5. Called by :
!
!       W3SIN3   Wind input Source term routine.
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
!       !/T      Enable test output.
!
! 10. Source code :
!-----------------------------------------------------------------------------!
      USE W3CONSTANTS, ONLY: GRAV
      USE W3GDATMD,  ONLY: ZZWND
      IMPLICIT NONE
      REAL, intent(in) :: WINDSPEED,TAUW
      REAL, intent(out) :: USTAR,Z0
      ! local variables
      real a,b  ! constants of parameterisation
      real cd   ! drag coefficient
      REAL SQRTCDM1
      REAL X,XI,DELI1,DELI2,XJ,delj1,delj2
      REAL UST,DELTOLD,TAUW_LOCAL
      INTEGER IND,J
   !
      TAUW_LOCAL=MAX(MIN(TAUW,TAUWMAX),0.)
      XI      = SQRT(TAUW_LOCAL)/DELTAUW
      IND     = MIN ( ITAUMAX-1, INT(XI)) ! index for stress table
      DELI1   = MIN(1.,XI - REAL(IND))  !interpolation coefficient for stress table
      DELI2   = 1. - DELI1
      XJ      = WINDSPEED/DELU
      J       = MIN ( JUMAX-1, INT(XJ) )
      DELJ1   = MIN(1.,XJ - REAL(J))
      DELJ2   = 1. - DELJ1
      USTAR=(TAUT(IND,J)*DELI2+TAUT(IND+1,J  )*DELI1)*DELJ2 &
       + (TAUT(IND,J+1)*DELI2+TAUT(IND+1,J+1)*DELI1)*DELJ1
      !
      ! Determines roughness length
      !
      SQRTCDM1  = MIN(WINDSPEED/USTAR,100.0)
      Z0  = ZZWND*EXP(-KAPPA*SQRTCDM1)
!
      RETURN
      END SUBROUTINE CALC_USTAR
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SDS3 (A, K, EMEAN, FMEAN, WNMEAN, USTAR, USDIR,    &
                                                          DEPTH, S, D)
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  !            F. Ardhuin             !
!/                  |                        FORTRAN 90 |
!/                  | Last update :         14-Aug-2006 |
!/                  +-----------------------------------+
!/
!/    05-Dec-1996 : Final FORTRAN 77                    ( version 1.18 )
!/    08-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    14-Aug-2006 : Generic WAM4+ dissipation term      ( version 2.22SHOM )
!/
!  1. Purpose :
!
!     Calculate whitecapping source term and diagonal term of derivative.
!
!  2. Method :
!
!       WAM-Cycle 4 and following.
!       The last update (09-May-2005) follows the redefinition of
!       the mean wavenumber as in Bidlot et al. (2005).
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       A       R.A.  I   Action density spectrum (1-D).
!       K       R.A.  I   Wavenumber for entire spectrum.          *)
!       EMEAN   Real  I   Mean wave energy.
!       FMEAN   Real  I   Mean wave frequency.
!       WNMEAN  Real  I   Mean wavenumber.
!       USTAR   Real  I   Friction velocity.
!       USDIR   Real  I   wind stress direction.
!       DEPTH   Real  I   Water depth.
!       S       R.A.  O   Source term (1-D version).
!       D       R.A.  O   Diagonal term of derivative.             *)
!     ----------------------------------------------------------------
!                         *) Stored in 1-D array with dimension NTH*NK
!
!  4. Subroutines used :
!
!       STRACE    Subroutine tracing.                 ( !/S switch )
!       PRT2DS    Print plot of spectrum.             ( !/T0 switch )
!       OUTMAT    Print out matrix.                   ( !/T1 switch )
!
!  5. Called by :
!
!       W3SRCE   Source term integration.
!       W3EXPO   Point output program.
!       GXEXPO   GrADS point output program.
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
!     !/S   Enable subroutine tracing.
!     !/T   Enable general test output.
!     !/T0  2-D print plot of source term.
!     !/T1  Print arrays.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      USE W3GDATMD, ONLY: NSPEC, NTH, NK, SSDSBR, DDELTA1, DDELTA2,   &
                          SSDSC1, SSDSC2, SSDSC3, SSDSC4, SSDSC5, SSDSC6,    &
                          SIG, SSDSP, SSDSLF, ECOS, ESIN, DTH, DSIP,  &
                          SSDSHF, SSDSDTH, SSDSBR2, SSDSBM
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      REAL, INTENT(IN)        :: A(NSPEC), K(NK), DEPTH, USTAR, USDIR, &
                                 EMEAN, FMEAN, WNMEAN
      REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: IS, IS0, IA
      INTEGER                 :: IK, ITH, I_INT, IK2
      REAL                    :: FACTOR, COSWIND
      REAL                    :: ALFAMEAN, WNMEAN2
      REAL                    :: FACTURB, DTURB, DCUMULATIVE, BREAKFRACTION
      REAL                    :: BTH0(NK)  !saturation spectrum
      REAL                    :: BTH(NSPEC)  !saturation spectrum
      REAL                    :: BSIGBAJ, SATURATION, SATURATION2, FACSAT
      REAL                    :: W, P0, MICHE, X
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Common factor
!
      W=26
      ALFAMEAN=WNMEAN**2*EMEAN
      FACTOR = SSDSC1 * TPI*FMEAN * ALFAMEAN**2
!
! 2.  Source term
!
      FACTURB=SSDSC5*USTAR**2/GRAV*DAIR/DWAT
      BREAKFRACTION=0.
      WNMEAN2 = MAX( 1.E-10 , WNMEAN  )
 
      DO  IK=1, NK
        IS0=1+(IK-1)*NTH
	BTH(IS0)=0.
	FACSAT=SIG(IK)*K(IK)**3*DTH
	IS0=(IK-1)*NTH
        BTH0(IK)=SUM(A(IS0+1:IS0+NTH))*FACSAT
        IF (SSDSDTH.GE.180) THEN
	! integrates around full circle
	  BTH(IS0+1:IS0+NTH)=BTH0(IK)
	ELSE
	! partial integration
          DO ITH=1,NTH
            IS=ITH+(IK-1)*NTH
            BTH(IS)=SUM(A(IS0+SATINDICES(ITH,:)))*FACSAT
	    END DO
	  END IF
!
        DO ITH=1,NTH
          IS=ITH+(IK-1)*NTH
          SATURATION2=TANH(10*(((BTH(IS)/SSDSBR)**0.5)-SSDSBR2))
          BREAKFRACTION=0.
          DCUMULATIVE=0.
          !IF (SSDSC6.NE.0) THEN
	  !  DO  IK2=1, IK
          !    IS0=ITH+(IK2-1)*NTH
          !    BREAKFRACTION=BREAKFRACTION + (MAX(BTH(IS0)-SSDSBR,0.))         &
          !                 * DSIP(IK2)/SIG(IK2)*(SIG(IK2)/SIG(IK))**2
          !    ENDDO
          !    DCUMULATIVE=BREAKFRACTION*SSDSC4*(BTH(IS)/SSDSBR)*SIG(IK)
          !  ENDIF
	!DTURB=-2.*SSDSC5*SIG(IK)*K(IK)*FACTURB
!
!    Original WAM4/WAM4+ dissipation term
!
          BSIGBAJ=FACTOR*(DDELTA1*K(IK)/WNMEAN2 + DDELTA2*(K(IK)/WNMEAN2)**2)
          X=TANH(MIN(K(IK)*DEPTH,10.))
!
!   Correction of saturation level for shallow-water kinematics
!
	  IF (SSDSBM(0).EQ.1) THEN
	    MICHE=1.
	  ELSE
	    MICHE=(X/(SSDSBM(1)+X*(SSDSBM(2)+X*(SSDSBM(3)+X*SSDSBM(4)))))**2
	    END IF
!
            COSWIND=(ECOS(IS)*COS(USDIR)+ESIN(IS)*SIN(USDIR))
            DTURB=-2.*SIG(IK)*K(IK)*FACTURB*COSWIND  ! Theory -> stress direction
            P0=SSDSP-0.5*SSDSC3*(1-TANH(W*USTAR*K(IK)/SIG(IK)-0.1))  ! for SDSC3=1 this is vdW et al.
 
            D(IS)=SSDSC2 * SIG(IK)                                           &
                 * ( SSDSC6*(MAX(0.,BTH0(IK)/(SSDSBR*MICHE)-SSDSC4))**P0     &
		    +(1-SSDSC6)*(MAX(0.,BTH(IS)/(SSDSBR*MICHE)-SSDSC4))**P0) &
                 + (DCUMULATIVE+DTURB)
 
           D(IS)=  D(IS) + BSIGBAJ*SSDSLF *0.5*(1-SATURATION2)         &          !low-frequency part
                         + BSIGBAJ*SSDSHF *0.5*(SATURATION2+1)                    !high-frequency part
 
           END DO
        END DO
!
      S = D * A
!
! ... Test output of arrays
!
      RETURN
!
! Formats
!
!/
!/ End of W3SDS3 ----------------------------------------------------- /
!/
      END SUBROUTINE W3SDS3
 
!/ ------------------------------------------------------------------- /
      SUBROUTINE INPTAB
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         21-Feb-2004 |
!/                  +-----------------------------------+
!/
!/    03-Jun-1996 : Final version 1.18 / FORTRAN 77 version.
!/    06-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    21-Feb-2004 : Multiple model version.             ( version 3.06 )
!/
!  1. Purpose :
!
!     Generate an interpolation table for the air-sea interaction
!     parameter of Chalikov and Belevich (1993).
!
!  2. Method :
!
!     The size of the table is set in parameter statements, the range
!     is set by the input parameters of this routine. The first counter
!     of the table corresponds to the nondimensional frequency
!
!                  SIGMA Ul
!        SIGA  =  ----------  COS ( THETA - THETA     )           (1)
!                     g                          wind
!
!     The second counter of the table represents the drag coefficient.
!     The maximum values of both parameters are passed to the routine
!     through the parameter list.
!
!  3. Parameters :
!
!  4. Subroutines used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      W3BETA    Func. Internal Function to calculate the
!                               interaction parameter.
!     ----------------------------------------------------------------
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3IOGR    Subr. W3IOGRMD Model definition IO routine.
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
!     !/S   Enable subroutine tracing.
!     !/T   Enable test output.
!     !/T0  Print table.
!     !/T1  Estimate maximum errors.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3CONSTANTS
      USE W3ODATMD, ONLY: NDST
      USE W3GDATMD, ONLY: SSWELLFPAR!
!/
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: ISIGA, IDRAG
      REAL                    :: SIGA, DRAG
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Determine range and increments of table ------------------------ *
 
!
      DSIGA  = SIGAMX / REAL(NRSIGA)
      DDRAG  = DRAGMX / REAL(NRDRAG)
!
! 2.  Fill table ----------------------------------------------------- *
!
      DO ISIGA=-NRSIGA,NRSIGA+1
        SIGA   = REAL(ISIGA) * DSIGA
        DO IDRAG=1, NRDRAG+1
          DRAG   = REAL(IDRAG) * DDRAG
          IF (SSWELLFPAR.LE.2) BETATB(ISIGA,IDRAG) = W3BETA ( SIGA, DRAG , NDST )
          !WRITE(994,*) 'TEST:',ISIGA,IDRAG,SIGA, DRAG, BETATB(ISIGA,IDRAG)
          END DO
        END DO
!
! 3.  Test output ---------------------------------------------------- *
!
      RETURN
!
! Formats
!
!/
!/    Internal function W3BETA
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      REAL FUNCTION W3BETA ( OMA , CL , NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |            D.Chalikov             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         21-Feb-2004 |
!/                  +-----------------------------------+
!/
!/    06-Dec-1996 : Final version 1.18 / FORTRAN 77 version.
!/    06-Dec-1999 : Upgrade to FORTRAN 90               ( version 2.00 )
!/    21-Feb-2004 : Multiple model version.             ( version 3.06 )
!/
!  1. Purpose :
!
!     Calculate wind-wave interaction parameter beta.
!
!  2. Method :
!
!     Chalikov and Belevich (1992), see also manual.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       W3BETA  Real  O   Wind-wave interaction parameter multiplied
!                         by density ratio.
!       OMA     Real  I   Non-dimensional apparent frequency.
!
!                         OMA = OMEGA | U | cos(theta-theta ) / g
!                                        l                 w
!
!       CL      Real  I   Drag coefficient at height l
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!  5. Called by :
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
!     !/S   Enable subroutine tracing.
!     !/T0  Enable test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: NDST
      REAL, INTENT(IN)        :: OMA, CL
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      REAL                    :: OM1, OM2, A0, A1, A2, A3, A4, A5,    &
                                 A6, A7, A8, A9, A10
!/
!/ ------------------------------------------------------------------- /
!/
!
! calculate Omegas
!
      OM1    =  1.075 +  75.*CL
      OM2    =  1.2   + 300.*CL
!
! calculate factors a
!
      A1     =  0.25  + 395.*CL
      A2     =  0.35  + 150.*CL
      A4     =  0.3   + 300.*CL
      A9     =  0.35  + 240.*CL
      A10    = -0.06  + 470.*CL
!
      A5     =  A4 * OM1
      A0     =  0.25 * A5**2 / A4
      A3     = (A0-A2-A1) / (A0+A4+A5)
      A6     =  A0 * (1.-A3)
      A7     = (A9*(OM2-1)**2+A10) / (OM2-OM1)
      A8     =  A7 * OM1
!
! calculate beta * 1.e4
!
      IF  ( OMA .LT. -1. ) THEN
          W3BETA = -A1 * OMA**2 - A2
        ELSE IF (OMA .LT. OM1/2.) THEN
          W3BETA =  A3 * OMA * ( A4 * OMA - A5 ) - A6
        ELSE IF (OMA .LT. OM1) THEN
          W3BETA =       OMA * ( A4 * OMA - A5 )
        ELSE IF (OMA .LT. OM2) THEN
          W3BETA = A7 * OMA - A8
        ELSE
          W3BETA = A9 * (OMA-1.)**2 + A10
        END IF
        W3BETA=W3BETA*1.E-4*1026/1.225
!
! WARNING: mutiplication by DRAT=DAIR/DWAT is done in W3SIN3
!
      RETURN
!
! Formats
!
!/
!/ End of W3BETA ----------------------------------------------------- /
!/
      END FUNCTION W3BETA
!/
!/ End of INPTAB ----------------------------------------------------- /
!/
!/ ------------------------------------------------------------------- /
!/
      END SUBROUTINE INPTAB
 
!/
!/ End of module W3SRC3MD -------------------------------------------- /
!/
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE KZEONE(X, Y, RE0, IM0, RE1, IM1)
!  June 1999 adaptation to CRESTb, all tests on range of (x,y) have been
!  bypassed, we implicitly expect X to be positive or |x,y| non zero
!
! This subroutine is copyright by ACM
! see http://www.acm.org/pubs/copyright_policy/softwareCRnotice.html
! ACM declines any responsibility of any kind
!
! THE VARIABLES X AND Y ARE THE REAL AND IMAGINARY PARTS OF
! THE ARGUMENT OF THE FIRST TWO MODIFIED BESSEL FUNCTIONS
! OF THE SECOND KIND,K0 AND K1.  RE0,IM0,RE1 AND IM1 GIVE
! THE REAL AND IMAGINARY PARTS OF EXP(X)*K0 AND EXP(X)*K1,
! RESPECTIVELY.  ALTHOUGH THE REAL NOTATION USED IN THIS
! SUBROUTINE MAY SEEM INELEGANT WHEN COMPARED WITH THE
! COMPLEX NOTATION THAT FORTRAN ALLOWS, THIS VERSION RUNS
! ABOUT 30 PERCENT FASTER THAN ONE WRITTEN USING COMPLEX
! VARIABLES.
! ACM Libraries
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   IMPLICIT NONE
   DOUBLE PRECISION X, Y, X2, Y2, RE0, IM0, RE1, IM1, &
      R1, R2, T1, T2, P1, P2, RTERM, ITERM, L
   DOUBLE PRECISION , PARAMETER, DIMENSION(8) :: EXSQ = &
         (/ 0.5641003087264D0,0.4120286874989D0,0.1584889157959D0, &
            0.3078003387255D-1,0.2778068842913D-2,0.1000044412325D-3, &
            0.1059115547711D-5,0.1522475804254D-8 /)
   DOUBLE PRECISION , PARAMETER, DIMENSION(8) :: TSQ = &
         (/ 0.0D0,3.19303633920635D-1,1.29075862295915D0, &
            2.95837445869665D0,5.40903159724444D0,8.80407957805676D0, &
            1.34685357432515D1,2.02499163658709D1 /)
   INTEGER N,M,K
! THE ARRAYS TSQ AND EXSQ CONTAIN THE SQUARE OF THE
! ABSCISSAS AND THE WEIGHT FACTORS USED IN THE GAUSS-
! HERMITE QUADRATURE.
      R2 = X*X + Y*Y
      IF (R2.GE.1.96D2) GO TO 50
      IF (R2.GE.1.849D1) GO TO 30
! THIS SECTION CALCULATES THE FUNCTIONS USING THE SERIES
! EXPANSIONS
      X2 = X/2.0D0
      Y2 = Y/2.0D0
      P1 = X2*X2
      P2 = Y2*Y2
      T1 = -(DLOG(P1+P2)/2.0D0+0.5772156649015329D0)
! THE CONSTANT IN THE PRECEDING STATEMENT IS EULER*S
! CONSTANT
      T2 = -DATAN2(Y,X)
      X2 = P1 - P2
      Y2 = X*Y2
      RTERM = 1.0D0
      ITERM = 0.0D0
      RE0 = T1
      IM0 = T2
      T1 = T1 + 0.5D0
      RE1 = T1
      IM1 = T2
      P2 = DSQRT(R2)
      L = 2.106D0*P2 + 4.4D0
      IF (P2.LT.8.0D-1) L = 2.129D0*P2 + 4.0D0
      DO 20 N=1,L
        P1 = N
        P2 = N*N
        R1 = RTERM
        RTERM = (R1*X2-ITERM*Y2)/P2
        ITERM = (R1*Y2+ITERM*X2)/P2
        T1 = T1 + 0.5D0/P1
        RE0 = RE0 + T1*RTERM - T2*ITERM
        IM0 = IM0 + T1*ITERM + T2*RTERM
        P1 = P1 + 1.0D0
        T1 = T1 + 0.5D0/P1
        RE1 = RE1 + (T1*RTERM-T2*ITERM)/P1
        IM1 = IM1 + (T1*ITERM+T2*RTERM)/P1
   20 CONTINUE
      R1 = X/R2 - 0.5D0*(X*RE1-Y*IM1)
      R2 = -Y/R2 - 0.5D0*(X*IM1+Y*RE1)
      P1 = DEXP(X)
      RE0 = P1*RE0
      IM0 = P1*IM0
      RE1 = P1*R1
      IM1 = P1*R2
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE INTEGRAL
! REPRESENTATION, EQN 3, EVALUATED WITH 15 POINT GAUSS-
! HERMITE QUADRATURE
   30 X2 = 2.0D0*X
      Y2 = 2.0D0*Y
      R1 = Y2*Y2
      P1 = DSQRT(X2*X2+R1)
      P2 = DSQRT(P1+X2)
      T1 = EXSQ(1)/(2.0D0*P1)
      RE0 = T1*P2
      IM0 = T1/P2
      RE1 = 0.0D0
      IM1 = 0.0D0
      DO 40 N=2,8
        T2 = X2 + TSQ(N)
        P1 = DSQRT(T2*T2+R1)
        P2 = DSQRT(P1+T2)
        T1 = EXSQ(N)/P1
        RE0 = RE0 + T1*P2
        IM0 = IM0 + T1/P2
        T1 = EXSQ(N)*TSQ(N)
        RE1 = RE1 + T1*P2
        IM1 = IM1 + T1/P2
   40 CONTINUE
      T2 = -Y2*IM0
      RE1 = RE1/R2
      R2 = Y2*IM1/R2
      RTERM = 1.41421356237309D0*DCOS(Y)
      ITERM = -1.41421356237309D0*DSIN(Y)
! THE CONSTANT IN THE PREVIOUS STATEMENTS IS,OF COURSE,
! SQRT(2.0).
      IM0 = RE0*ITERM + T2*RTERM
      RE0 = RE0*RTERM - T2*ITERM
      T1 = RE1*RTERM - R2*ITERM
      T2 = RE1*ITERM + R2*RTERM
      RE1 = T1*X + T2*Y
      IM1 = -T1*Y + T2*X
      RETURN
! THIS SECTION CALCULATES THE FUNCTIONS USING THE
! ASYMPTOTIC EXPANSIONS
   50 RTERM = 1.0D0
      ITERM = 0.0D0
      RE0 = 1.0D0
      IM0 = 0.0D0
      RE1 = 1.0D0
      IM1 = 0.0D0
      P1 = 8.0D0*R2
      P2 = DSQRT(R2)
      L = 3.91D0+8.12D1/P2
      R1 = 1.0D0
      R2 = 1.0D0
      M = -8
      K = 3
      DO 60 N=1,L
        M = M + 8
        K = K - M
        R1 = FLOAT(K-4)*R1
        R2 = FLOAT(K)*R2
        T1 = FLOAT(N)*P1
        T2 = RTERM
        RTERM = (T2*X+ITERM*Y)/T1
        ITERM = (-T2*Y+ITERM*X)/T1
        RE0 = RE0 + R1*RTERM
        IM0 = IM0 + R1*ITERM
        RE1 = RE1 + R2*RTERM
        IM1 = IM1 + R2*ITERM
   60 CONTINUE
      T1 = DSQRT(P2+X)
      T2 = -Y/T1
      P1 = 8.86226925452758D-1/P2
! THIS CONSTANT IS SQRT(PI)/2.0, WITH PI=3.14159...
      RTERM = P1*DCOS(Y)
      ITERM = -P1*DSIN(Y)
      R1 = RE0*RTERM - IM0*ITERM
      R2 = RE0*ITERM + IM0*RTERM
      RE0 = T1*R1 - T2*R2
      IM0 = T1*R2 + T2*R1
      R1 = RE1*RTERM - IM1*ITERM
      R2 = RE1*ITERM + IM1*RTERM
      RE1 = T1*R1 - T2*R2
      IM1 = T1*R2 + T2*R1
      RETURN
      END SUBROUTINE KZEONE
 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
SUBROUTINE KERKEI(X,KER,KEI)
!**********************************************************************
! Computes the values of the zeroth order Kelvin function Ker and Kei
! These functions are used to determine the friction factor fw as a
! function of the bottom roughness length assuming a linear profile
! of eddy viscosity (See Grant and Madsen, 1979)
!**********************************************************************
   IMPLICIT NONE
 
   DOUBLE PRECISION ZR,ZI,CYR,CYI,CYR1,CYI1
   INTEGER NZ,IERR
   REAL X,KER,KEI
 
   ZR=X*.50D0*SQRT(2.0D0)
   ZI=ZR
   CALL KZEONE(ZR, ZI, CYR, CYI,CYR1,CYI1)
   KER=CYR/EXP(ZR)
   KEI=CYI/EXP(ZR)
END SUBROUTINE KERKEI
 
      END MODULE W3SRC3MD
