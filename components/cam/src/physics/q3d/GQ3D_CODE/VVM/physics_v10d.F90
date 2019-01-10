! This program uses its own variables and parameters.
! Only inputs from the main code are the arguments of the subroutines and
! the parameters mi1, mj1, and nk2.
!
!  Add "LocalP" option: using local pressure to calculate Temp. (and local density) 
!
!  Contains:
!
!  SUBROUTINE INITIALIZE_PHYSICS
!  SUBROUTINE SATURATION
!  SUBROUTINE CONDENSATION
!  SUBROUTINE SUBLIMATION
!  Subroutine MICROPHYSICS
!  SUBROUTINE IMICRO
!  SUBROUTINE CAL_ACR3
!  SUBROUTINE CAL_IDW
!  Subroutine BERGERON
!
!  FUNCTION GGESW
!  FUNCTION GGESI
!  FUNCTION ESATW
!  FUNCTION ESATI
!  FUNCTION DESWDT
!  FUNCTION DESIDT
!-----------------------------------------------------------------------------

!=============================================================================
  SUBROUTINE INITIALIZE_MICROPHYSICS
!=============================================================================
!   INPUT:
!      DT,SFCRHO (module: mphysics_parameters)
!
!   OUTPUT:
!      CRACS,CSACR,CGACR,CGACS     out
!      ACCO(3,4)                   out
!      CSACW,CWACS,CIACR           out
!      CRACI,CSACI,CGACW,CGACI     out
!      CGACW,CGACI,CRACW           out
!      CSSUB(5),CGSUB(5)           out
!      CREVP(5),CGFR(2),CSMLT(5)   out
!      CGMLT(5),CGWET(4),CRAUT(2)  out
!      QI0,QS0,QL0,ES0,CES0        out (constant parameter)
!      C1BRG,C2BRG                 out (constant parameter)
!      RMI50,RMI40                 out (constant parameter)
!      CPLC,CPLS,CPLF              out (thermodynamic constant)
!      TDIF                        out (thermodynamic constant, TICE-TTFRZ = -40 C )
!      ESW00                       out [=GGESW(TICE-TDIF=TTFRZ)]
!      ESI00                       out [=GGESI(TICE-TDIF=TTFRZ)]
!      ESWT(151),ESIT(151)         out (saturation pressure table)
!      DESWT(151),DESIT(151)       out (saturation pressure derivative table)
!      TT                          out (=TICE-101. lower bound  used in saturation table)
!      TTD                         out (=TICE-100.5 lower bound used in
!                                                           saturation derivation table)
!      HLTS,HLTC,HLTF              out (latent heating constants)
!      CH2O                        out (constant)
!      CICE,TICE                   out (constants)
!      CP                          out (=1004.0 spec. heat under const. pres. constant)
!      EPS                         out (=062197 constant)
!      TTFRZ                       out (=233.1K constant)
!      VCONR,VCONS,VCONG           out (terminal velocities of precipitating quantities)
!      TDT                         out (currently =DT)
!      RHOSFC                      out (=SFCRHO)
!
!   INTERNAL:
!      PIE                         internal (=3.14159265 constant)
!      GRAV                        internal (=9.81 gravitational acceleration constant)
!      VDIFU,TCOND                 internal (constant parameter)
!      RVARP                       internal (?)
!      RDRYA                       internal (=287 dry air gas constant)
!      VISK                        internal (constant)
!      SCM3                        internal
!      DTL                         internal (currently =DT)
!      RI50                        internal (constant parameter)
!      PISQ, AA22, GCON            internal
!
!   CALLS:
!      FUNCTION GGESW(T)
!      FUNCTION GGESI(T)
!---------------------------------------------------------------------------------------------------

USE shr_kind_mod, only: r8 => shr_kind_r8

USE constld, only: dt,pi=>piv,gas_const,wvg_const,cp_val=>cp,grav_val=>grav, &
                   hlf,hlm,hls,epsilo,t_melt,cp_h2o,cp_ice,rho_h2o

Use mphysics_parameters, only: ri50,cracs,csacr,cgacr,cgacs,acco, &
                               csacw,cwacs,ciacr,craci,csaci,cgacw,cgaci,cracw, &
                               cssub,cgsub,crevp,cgfr,csmlt,cgmlt,cgwet,craut, &
                               qi0,qs0,ql0,es0,ces0,c1brg,c2brg, &
                               rmi50,rmi40,cplc,cpls,cplf, &
                               tdif,esw00,esi00,eswt,esit,&
                               deswt,desit,tt,ttd, &
                               pie,grav,vdifu,tcond,rvapr,rdrya,visk, &
                               hlts,hltc,hltf,ch2o,cice,tice,cp,eps,ttfrz, &
                               vconr,vcons,vcong, &
                               dtl,tdt,rhosfc,sfcrho

   IMPLICIT NONE
   
   REAL (KIND=r8) ::   &
          ACT(8),ACC(3),RNZR,RNZS,RNZG,RHOR,RHOS,RHOG,ALIN,CLIN, &
          GAM263,GAM275,GAM290,GAM325,GAM350,GAM380,GAM425,      &
          GAM450,GAM480,GAM625,GAM680
   REAL (KIND=r8) :: PISQ,SCM3,AA22,CD,GCON,GGESW, GGESI

   integer :: i, k  ! do loop indices
   integer :: itc

      RNZR = 8.0e6_r8     ! intercept parameter of the raindrop size distribution
      RNZS = 3.0e6_r8     ! intercept parameter of the snowflake size distribution
      RNZG = 4.0e6_r8     ! intercept parameter of the hailstone size distribution

      RHOR = rho_h2o      ! density of liquid water   (kg/m^3)
      RHOS = 0.1e3_r8     ! density of snow    (kg/m^3)
      RHOG = 0.4e3_r8     ! density of graupel (kg/m^3)

      ALIN = 841.99667_r8
      CLIN = 4.836071224_r8

      ACC(1) = 5.0_r8
      ACC(2) = 2.0_r8
      ACC(3) = 0.5_r8

      CRAUT(1) = 1.2e-1_r8
      CRAUT(2) = 1.064e-2_r8

      GAM263 = 1.456943_r8
      GAM275 = 1.608355_r8
      GAM290 = 1.827363_r8
      GAM325 = 2.54925_r8
      GAM350 = 3.323363_r8
      GAM380 = 4.694155_r8
      GAM425 = 8.285063_r8
      GAM450 = 11.631769_r8
      GAM480 = 17.837789_r8
      GAM625 = 184.860962_r8
      GAM680 = 496.604067_r8

      TCOND = 2.36e-2_r8               ! thermal conductivity of air
      RMI50 = 3.84e-6_r8
      RMI40 = 2.46e-7_r8
      RI50  = 1.0e-4_r8
      QS0   = 1.0e-3_r8
      QI0   = 0.6e-3_r8
      QL0   = 0.5e-3_r8

      VDIFU = 2.11e-5_r8               ! diffusivity of water vapor in air (cm^2/s)
      VISK  = 1.259e-5_r8              ! kinematic viscosity of air (cm^2/s)
      SCM3 = (VISK/VDIFU)**(1.0/3.0)   ! Schmidt number

      TTFRZ = 233.1_r8                 ! T = -40.05 C

      PIE   = pi                       ! 3.14159265358979323846_r8
      GRAV  = grav_val                 ! 9.80616_r8
      HLTS  = hls                      ! 2.8336e6_r8
      HLTC  = hlf                      ! 2.500e6_r8
      HLTF  = hlm                      ! 3.336e5_r8
      RDRYA = gas_const                ! 287.042_r8
      RVAPR = wvg_const                ! 461.5_r8
      EPS   = epsilo                   ! 0.62197_r8
      CP    = cp_val                   ! 1.00464e3_r8
      CH2O  = cp_h2o                   ! 4.188e3_r8
      CICE  = cp_ice                   ! 2.11727e3_r8
      TICE  = t_melt                   ! 273.15_r8

      TDT = DT
      DTL = DT
      
      RHOSFC = SFCRHO
      
      CPLC = HLTC/CP
      CPLS = HLTS/CP
      CPLF = HLTF/CP
      PISQ = PIE*PIE

!-------------------------------------------------------------------
!     ACR3:  FOUR LEAD CONSTANTS REQUIRED, THREE FOR EACH SUMMATION
!            FOUR SEPARATE PROCESSES:  RACS,SACR,GACR,GACS
        CRACS = PISQ*RNZR*RNZS*RHOS;  CSACR = PISQ*RNZR*RNZS*RHOR
        CGACR = PISQ*RNZR*RNZG*RHOR;  CGACS = PISQ*RNZG*RNZS*RHOS

!------------------------------------------------
!     ACT:  1-2:RACS(S-R); 3-4:SACR(R-S);
!           5-6:GACR(R-G); 7-8:GACS(S-G)
        ACT(1)=PIE*RNZS*RHOS;  ACT(2)=PIE*RNZR*RHOR; ACT(6)=PIE*RNZG*RHOG
        ACT(3)=ACT(2);  ACT(4)=ACT(1);  ACT(5)=ACT(2);  ACT(7)=ACT(1)
        ACT(8)=ACT(6)

  do I=1,3
    do K=1,4
      ACCO(I,K) = ACC(I)/(ACT(2*K-1)**((7-I)*0.25)*ACT(2*K)**(I*0.25))
    enddo
  enddo

!----- TERMINAL VELOCITY CONSTANTS ----------------
        aa22 = 40.74_r8 * 40.74_r8
        CD = 4.0_r8 * GRAV * RHOG / 3.0_r8 / RHOSFC / aa22
        GCON  = SQRT(4.0_r8*GRAV*RHOG/3.0_r8/CD)

        VCONR = ALIN*GAM480/(6.0_r8*ACT(2)**0.20)
        VCONS = CLIN*GAM425/(6.0_r8*ACT(1)**0.0625)
        VCONG = GAM450*GCON/(6.0_r8*ACT(6)**0.125)

!------------------------------------------------
!     ACR1:  SINGLE CONSTANT REQUIRED
!     FIVE SEPARATE PROCESSES:  SACW,WACS,IACR,RACI,SACI
        CSACW = PIE*RNZS*CLIN*GAM325/(4.0_r8*ACT(1)**0.8125)
        CWACS = PISQ*RHOS*RNZS*CLIN*GAM625/(1.0056D-10*ACT(1)**1.5625)
        CIACR = PISQ*RHOR*RNZR*ALIN*GAM680/(1.0056D-11*ACT(2)**1.7)
        CRACI = PIE*RNZR*ALIN*GAM380/(4.0_r8*ACT(2)**0.95)
        CSACI = CSACW
!------------------------------------------------
!     ACR2:  SINGLE CONSTANT REQUIRED
!     CAGCI IS FOR DRY GROWTH
        CGACW = PIE*RNZG*GAM350*GCON/(4.0_r8*ACT(6)**0.875)
        CGACI = CGACW*0.1_r8
!-----------------------------
!     RACW
        CRACW = CRACI

!-------------------------------------------------------------------
!     SUBL AND REVP:  FIVE CONSTANTS FOR THREE SEPARATE PROCESSES
        CSSUB(1) = 2.0_r8*PIE*VDIFU*TCOND*RVAPR*RNZS
        CGSUB(1) = 2.0_r8*PIE*VDIFU*TCOND*RVAPR*RNZG
        CREVP(1) = 2.0_r8*PIE*VDIFU*TCOND*RVAPR*RNZR
        CSSUB(2) = 0.78_r8/SQRT(ACT(1))
        CGSUB(2) = 0.78_r8/SQRT(ACT(6))
        CREVP(2) = 0.78_r8/SQRT(ACT(2))
        CSSUB(3) = 0.31_r8*SCM3*GAM263*SQRT(CLIN/VISK)/ACT(1)**0.65625
        CGSUB(3) = 0.31_r8*SCM3*GAM275*SQRT(GCON/VISK)/ACT(6)**0.6875
        CREVP(3) = 0.31_r8*SCM3*GAM290*SQRT(ALIN/VISK)/ACT(2)**0.725
        CSSUB(4) = TCOND*RVAPR
        CSSUB(5) = HLTS**2*VDIFU
        CGSUB(4) = CSSUB(4)
        CREVP(4) = CSSUB(4)
        CGSUB(5) = CSSUB(5)
        CREVP(5) = HLTC**2*VDIFU

!------ GFR:  TWO CONSTANTS -------------------------------
      CGFR(1) = 20.0D2*PISQ*RNZR*RHOR/ACT(2)**1.75
      CGFR(2) = 0.66_r8

!----- SMLT:  FIVE CONSTANTS ( Lin et al. 1983 ) ----------
      CSMLT(1) = 2.0_r8*PIE*TCOND*RNZS/HLTF
      CSMLT(2) = 2.0_r8*PIE*VDIFU*RNZS*HLTC/HLTF
      CSMLT(3) = CSSUB(2);  CSMLT(4) = CSSUB(3);  CSMLT(5) = CH2O/HLTF

!----- GMLT:  FIVE CONSTANTS ------------------------------
      CGMLT(1) = 2.0_r8*PIE*TCOND*RNZG/HLTF
      CGMLT(2) = 2.0_r8*PIE*VDIFU*RNZG*HLTC/HLTF
      CGMLT(3) = CGSUB(2);  CGMLT(4) = CGSUB(3);  CGMLT(5) = CH2O/HLTF

!----- GWET:  Two constants +plus -------------------------
!             calculation of saturation vapor pressure at T=0 deg C
!             ventilation coefficient of 10 is included
      CGWET(1) = CGMLT(1);  CGWET(2) = CGMLT(2);     CGWET(3) = CGMLT(3)
      CGWET(4) = CGMLT(4)
      ES0  = 6.107799961_r8
      CES0 = EPS*ES0

!----- BERGRN: TWO CONSTANTS ------------------------------
!              C2BRG HAS CONVERSION FACTOR OF 10**3

      C1BRG = DTL/RMI50;  C2BRG = PIE*RI50**2*1.D3

!------ SATURATION VAPOR PRESSURE VALUES ARE TABULATED EACH WHOLE
!       DEG. (C) FOR -100 < TC < 50 FOR WATER AND -100 < TC < 0 FOR ICE
!       DERIVATIVES ARE CENTERED AT EACH HALF DEG. (C)
      TT  = TICE-101.0_r8
      TTD = TICE-100.5_r8

      ESWT(1) = GGESW(TICE-100.0_r8)
      ESIT(1) = GGESI(TICE-100.0_r8)

    ITC = -99
  do I=2,151
    ESWT(I)    = GGESW(TICE+FLOAT(ITC));  DESWT(I-1) = ESWT(I)-ESWT(I-1)
    ITC        = ITC+1
  enddo
    DESWT(151) = 6.25_r8

    ITC = -99
  do I=2,101
    ESIT(I)    = GGESI(TICE+FLOAT(ITC));  DESIT(I-1) = ESIT(I)-ESIT(I-1)
    ITC = ITC+1
  enddo

    DESIT(101) = DESWT(101)
  do I=102,151
    ESIT(I)  = ESWT(I);  DESIT(I) = DESWT(I)
  enddo

!-------------------------------------------------------------------
!     SATURATION ADJUSTMENT:  RANGE OF MIXED ICE-WATER SATURATION
!     AND SATURATION VAPOR PRESSURE OVER WATER AND ICE AT T = -40 C.
        TDIF  = TICE-TTFRZ;  ESW00 = GGESW(TICE-TDIF);  ESI00 = GGESI(TICE-TDIF)

RETURN
END Subroutine initialize_microphysics

!=============================================================================
  SUBROUTINE INITIALIZE_MICROPHYSICS_LOC (SFC_RHO)
!=============================================================================
! Initialize the local Parameters (Used if LocalP=.T.)

   USE shr_kind_mod, only: r8 => shr_kind_r8
   Use mphysics_parameters, only: cgacw,cgaci,cgsub,cgmlt,cgwet, &
                                  pie,grav,vdifu,visk,vcong,rhosfc

   IMPLICIT NONE

   REAL (kind=r8), INTENT(IN) :: SFC_RHO   
   REAL (KIND=r8) :: ACT_6,RNZG,RHOG,GAM275,GAM350,GAM450
   REAL (KIND=r8) :: SCM3,AA22,CD,GCON

      RHOSFC = SFC_RHO
      
      RNZG = 4.0e6_r8     ! intercept parameter of the hailstone size distribution
      RHOG = 0.4e3_r8     ! density of graupel (kg/m^3)

      GAM275 = 1.608355_r8
      GAM350 = 3.323363_r8
      GAM450 = 11.631769_r8
    
      SCM3 = (VISK/VDIFU)**(1.0/3.0)   ! Schmidt number
      ACT_6 = PIE*RNZG*RHOG

!----- TERMINAL VELOCITY CONSTANTS ----------------
        aa22 = 40.74_r8 * 40.74_r8
        CD = 4.0_r8 * GRAV * RHOG / 3.0_r8 / RHOSFC / aa22
        GCON  = SQRT(4.0_r8*GRAV*RHOG/3.0_r8/CD)

        VCONG = GAM450*GCON/(6.0_r8*ACT_6**0.125)
!------------------------------------------------
        CGACW = PIE*RNZG*GAM350*GCON/(4.0_r8*ACT_6**0.875)
        CGACI = CGACW*0.1_r8
!-------------------------------------------------------------------
!     SUBL AND REVP:  FIVE CONSTANTS FOR THREE SEPARATE PROCESSES
       CGSUB(3) = 0.31_r8*SCM3*GAM275*SQRT(GCON/VISK)/ACT_6**0.6875
        
!----- GMLT:  FIVE CONSTANTS ------------------------------
       CGMLT(4) = CGSUB(3)

!----- GWET:  Two constants +plus -------------------------
       CGWET(4) = CGMLT(4)

RETURN
END Subroutine initialize_microphysics_loc

!=============================================================================
  SUBROUTINE SATURATION (IPOINT,JPOINT,KPOINT,T,P,QVK,QLK,QIK,QFK)
!=============================================================================
!------------------------------------------------------------------
!   INPUT:
!      IPOINT,JPOINT,KPOINT  in  (i,j,k for identification only)
!      T,QVK,QLK,QIK         in  (to be modified by stauration)
!      P                     in  (pressure in mb)
!
!      CPLC,CPLS,CPLF          |
!      TDIF,ESW00,ESI00        |
!      ESWT,ESIT,DESWT,DESIT    >  Through Module mphysics_parameters
!      TT,TTD,HLTC,HLTF        |
!      TICE,EPS,TTFRZ          |
!
!   INTERNAL:
!      IERR, JERR, KERR, PRESS, TO, QVO, QLO, QIO, ISTOP, NOCONV, NLOOP
!      DTMP, DTFR, ESI, QSI, QS, ESW, QSW, DTMLT, QSIW, N, CRIT
!      GAMFAC, SS, GAMMAW, GAMMAI, TEM, EX1, TQIK, DQI, SST, QV1, QL1
!      QI1, SS1, T1, QSIW1, SS2, T2, QV3, DQ, DQL, DQ, DQV, DQLK, QL3
!      QI3, T3, QSIT, QSWT, QSIW3, SS3, QV2, QL2, QI2, T2, TQLK, QVKSV
!      TSV, QSIW2
!
!   OUTPUT:
!      T,QVK,QLK,QIK,QFK         out  (modified by stauration)
!
!   CALLS:
!      Function ESATI (T)
!      Function ESATW (T)
!      Call Condensation (T,P,QVK,QSW,ESW,QLK,ISTOP)
!      Call Sublimation (T,P,QVK,QSI,ESI,QIK,ISTOP)
!------------------------------------------------------------------

USE shr_kind_mod, only: r8 => shr_kind_r8
Use mphysics_parameters, only: cplc,cpls,cplf,tdif,ttfrz, &
                               esw00,esi00,eswt,esit,hltc,hltf,tice,eps

IMPLICIT NONE

    INTEGER, INTENT(IN) :: IPOINT, JPOINT, KPOINT
    REAL(kind=r8), INTENT(IN) :: P  ! pressure in mb
    REAL(kind=r8), INTENT(INOUT) :: &
       T,    & ! temperature to be modified
       QVK,  & ! qv to be modified
       QLK,  & ! qc to be modified
       QIK,  & ! qi to be modified
       QFK     ! amount of freezing of cloud water into cloud ice

    REAL(kind=r8) TO,QVO,QLO,QIO
    REAL(kind=r8) DTMP,DTFR,ESI,ESATI,QSI,QS, &
           ESW,ESATW,QSW,DTMLT,QSIW,CRIT,GAMFAC
    REAL(kind=r8) SS,GAMMAW,DESWAT,GAMMAI,DESIDT,TEM,EX1,TQIK, &
           QVKSV,TSV,DQI,SST,QV1,QL1,QI1,SS1,       &
           T1,QSIW1,DQL,QV2,QL2,QI2,T2,QSIT
    REAL(kind=r8) DESWDT,QSWT,QSIW2,SS2,QV3,DQ,QL3,QI3,T3,  &
           QSIW3,SS3,DQLK, DQV,TQLK, PRESS
    INTEGER :: i,j,k
    INTEGER :: ierr, jerr, kerr,  &  ! location indices
               istop, n, noconv, nloop

!     -------------------------------------------------------------
!     ----                                                     ----
!     ----    PRELIMINARY ADJUSTMENT CONSISTS OF FREEZING      ----
!     ----    SUSPENDED LIQUID WATER (SLW) OR MELTING          ----
!     ----    SUSPENDED ICE CRYSTALS (SIW) AND DEFINING        ----
!     ----    SATURATION MIXING RATIO (QS).                    ----
!     ----      IF    TC > 0    NO SIW --> QS = QSW            ----
!     ----            TC < -40  NO SLW --> QS = QSI            ----
!     ----      -40 < TC < 0    SLW AND SIW EXIST --> QS IS    ----
!     ----                      WEIGHTED AVG. OF QSW AND QSI.  ----
!     ----      -40 < TC < 0    NO SLW AND NO SIW --> QS IS A  ----
!     ----                      WEIGHTED AVG. OF QSW AND QSI   ----
!     ----                      THAT DEPENDS ON TC ONLY.       ----
!     ----                                                     ----
!     -------------------------------------------------------------

      IERR = IPOINT;  JERR = JPOINT;  KERR = KPOINT;  PRESS = P

      TO = T;  QVO = QVK;  QLO = QLK;  QIO = QIK
      ISTOP = 0;  NOCONV = 0;  NLOOP  = 0

      IF (QVK .LT. 1.0D-20) QVK=0.0_r8
      IF (QLK .LT. 1.0D-20) QLK=0.0_r8
      IF (QIK .LT. 1.0D-20) QIK=0.0_r8

      IF(T .GT. TICE)   GOTO 100

      IF(QLK .EQ. 0.0_r8)  GOTO 30
      IF(T   .LT. TTFRZ)  GOTO 10
      IF(QIK .EQ. 0.0_r8)  GOTO 110
      GOTO 130

!----- T < -40 DEG C ---------------------------
!     FREEZE ALL SLW UNLESS HEAT RELEASED WOULD
!----- RAISE TEMPERATURE ABOVE TTFRZ -----------

   10 DTMP = CPLF*QLK

      DTFR = TTFRZ-T
      IF(DTFR .GT. DTMP)  GOTO 20
      QIK = QIK+DTFR/CPLF
      QLK = QLK-DTFR/CPLF

      T = T+DTFR
!     WRITE(6,*) 'FREEZE LIQ TIL -40   DTMP,DTFR',DTMP,DTFR
      GOTO 130

   20 QIK = QIK+QLK
      QLK = 0.0_r8
      T   = T+DTMP
!     WRITE(6,*) 'FREEZE ALL LIQ   DTMP,DTFR ',DTMP,DTFR

   30 IF ( QIK .EQ. 0.0_r8 ) GOTO 40

!--- NO SLW --> QS = QSI ---------------------------
      ESI = ESATI (T)
      QSI = EPS*ESI/(P-ESI)
      QS  = QSI
!     WRITE(6,*) ' ICE BUT NO LIQ WATER '
      GOTO 200

!--- NO SLW AND NO SIW -- QS DEPENDS ONLY ON T ---
   40 ESW = ESATW (T)
      ESI = ESATI (T)
      QSW = EPS*ESW/(P-ESW)
      QSI = EPS*ESI/(P-ESI)
      QS  = ( QSW * MAX( T-TTFRZ ,0.0_r8 ) + QSI * MIN( TICE-T ,TDIF ) ) / TDIF
!     WRITE(6,*) ' NO LIQ WATER AND NO ICE'
      GOTO 200

!--------- T > 0 ------------------------------------
!       MELT ALL SIW UNLESS HEAT ABSORBED WOULD
!       LOWER TEMPERATURE BELOW TICE
  100 IF(QIK .EQ. 0.0_r8)  GO TO 110
      DTMP  =-CPLF*QIK
      DTMLT = TICE-T

      IF(DTMLT .GE. DTMP)  GO TO 120
      QLK = QLK+QIK
      QIK = 0.0_r8
      T = T+DTMP
!     WRITE(6,*) ' MELT ALL ICE   DTMP,DTMLT ',DTMP,DTMLT

!--- NO SIW --> QS = QSW ----------------------------
  110 ESW = ESATW (T)
      QSW = EPS*ESW/(P-ESW)
      QS  = QSW
!     WRITE(6,*) ' LIQ WATER BUT NO ICE '
      GO TO 200

  120 QLK = QLK-DTMLT/CPLF
      QIK = QIK+DTMLT/CPLF
      T = T+DTMLT
!     WRITE(6,*) ' MELT ICE TILL T=0   DTMP,DTMLT ',DTMP,DTMLT

!--- MIXED WATER AND ICE CASE ---------------------
  130 ESW = ESATW (T)
      ESI = ESATI (T)
      QSW = EPS*ESW/(P-ESW)
      QSI = EPS*ESI/(P-ESI)
      QS  = (QLK*QSW+QIK*QSI)/(QLK+QIK)
!     WRITE(6,*) ' ICE AND LIQ WATER '

!     ----------------------------------------------------
!     ----                                            ----
!     ----         ADJUSTMENT TO SATURATION           ----
!     ----                                            ----
!     ----     INITIALLY UNSATURATED                  ----
!     ----          CONDENSATION OF VAPOR TO SLW      ----
!     ----          SUBLIMATION OF SIW                ----
!     ----                                            ----
!     ----     INITIALLY SATURATED                    ----
!     ----          CONDENSATION OF VAPOR TO SLW      ----
!     ----          SUBLIMATION  OF VAPOR TO SIW      ----
!     ----                                            ----
!     ----     T > TICE ......... CONDENSATON ONLY    ----
!     ----     T < TTFRZ ........ SUBLIMATION ONLY    ----
!     ----     TTFRZ < T < TICE . BOTH PROCESSES      ----
!     ----                        OCCUR. PRODUCTION   ----
!     ----                        OF SLW AND SIW      ----
!     ----                        DEPENDS LINEARLY    ----
!     ----                        ON TEMPERATURE      ----
!     ----------------------------------------------------

  200 CONTINUE
!     WRITE(6,201)  QVK,QLK,QIK,QS,T
!  201 FORMAT('0',20X,'200'/10X,'QVK',13X,'QLK',13X,'QIK',13X,'QS',13X,'T',/5X,4E16.8,F12.7)

!-- OBTAIN THE AMOUNT OF FREEZING OF CLOUD WATER INTO CLOUD ICE -------
      QFK= QIK - QIO

      IF(DABS(QVK/QS-1.0_r8) .LE. 1.D-6)  GO TO 1000
      IF(QVK .GT. QS)  GO TO 500

!     ---------------------------------------------
!     ---------     UNSATURATED CASE      ---------
!     ---------------------------------------------

  210 CONTINUE
!     WRITE(6,*) 'UNSATURATED CASE'
      IF(QLK .EQ. 0.0_r8)  GO TO 300

!-- SLW PRESENT --> EVAPORATE ALL SLW -------------
      T = T-CPLC*QLK
      QVK = QVK+QLK
      QLK = 0.0_r8
!     WRITE(6,*) ' EVAPORATE ALL LIQ WATER '
!     WRITE(6,*) 'TEMP',T

!-- TEST FOR SATURATION ---------------------------
      IF(QIK .GT. 0.0_r8)  GO TO 250

!-- WATER SATURATION -----------------------------
!      (NO ICE PRESENT)
      ESW = ESATW (T)
      QSW = EPS*ESW/(P-ESW)
!     WRITE(6,*) ' QSW ',QSW,' QVK ',QVK

!-- ADJUST TO SATURATION (CONDENSATION) -----------
      IF(QVK/QSW-1.0_r8 .LE. 1.D-6)  GO TO 1000

      Call Condensation (T,P,QVK,QSW,ESW,QLK,ISTOP)

      GO TO 280

!-- MIXED WATER AND ICE SATURATION ----------------
  250 ESW = ESATW(T)
      ESI = ESATI (T)
      QSW = EPS*ESW/(P-ESW)
      QSI = EPS*ESI/(P-ESI)
      QSIW = QSI
!     WRITE(6,*) ' QSIW ',QSIW,' QVK ',QVK
      IF(QVK .LT. QSIW)  GO TO 310

!-- ADJUST TO SATURATION (CONDENSATION) -----------
      IF(QVK/QSI-1.0_r8 .LE. 1.D-6)  GO TO 1000
      N = 0
      CRIT = 1.D-5
      GAMFAC = EPS*CPLC*P

  260 N = N+1
!     WRITE(6,261)  QSW,QSI,QS,QSIW,QVK,QLK,QIK
!  261 FORMAT('0','   260'/5X,'QSW',13X,'QSI',13X,'QS',13X,'QSIW',13X,'QVK',13X,'QLK',13X,'QIK'/7E16.8)
      IF(N .GT. 15)  GO TO 290
      SS = QVK-QSIW
      GAMMAW = GAMFAC*DESWDT(T)/(P-ESW)**2
      GAMMAI = GAMFAC*DESIDT(T)/(P-ESI)**2
      TEM = 1.0_r8+(GAMMAW*QLK+GAMMAI*QIK)/(QLK+QIK)+QIK*(QSW-QSI)/(QLK+QIK)**2
      EX1 = SS/TEM
      IF(N .EQ.  1)  GO TO 270
      IF(ABS(EX1/QLK) .LT. CRIT)  GO TO 280
  270 T = T+CPLC*EX1
      ESW = ESATW (T)
      ESI = ESATI (T)
      QSW = EPS*ESW/(P-ESW)
      QSI = EPS*ESI/(P-ESI)
      QVK = QVK-EX1
      QLK = QLK+EX1
      QSIW = (QLK*QSW+QIK*QSI)/(QLK+QIK)
!     WRITE(6,111) N,T,SS,EX1,QSIW,QVK,QLK,TEM
!  111 FORMAT('0',I5,F14.7,6E16.8)
      GO TO 260

!-- CHECK THAT TEMPERATURE DOES NOT FALL BELOW TTFRZ WITH ----------
!      LIQUID WATER PRESENT
  280 IF(T .GE. TTFRZ)  GO TO 1000
!     WRITE(6,*) 'T<-40 WITH LIQ  NLOOP ',NLOOP
      IF(NLOOP .GT. 0)  GO TO 295
      NLOOP = 1
      GO TO 10

  290 ISTOP = 290
      GO TO 280

  295 ISTOP = 295
      GO TO 2000

!-- SUBLIMATION OF ALL SIW ------------------
  300 IF(QIK .EQ. 0.0_r8)  GO TO 1000

!-- SIW PRESENT --> SUBLIMATE ALL SIW -------
  310 T = T-QIK*CPLS
      QVK = QVK+QIK
      QIK = 0.0_r8

!     WRITE(6,*) ' NO WATER---SUBLIME ALL ICE '
!     WRITE(6,*) 'TEMP',T

!---- TEST FOR ICE SATURATION -----------------------
!         (NO SLW PRESENT)
      ESI = ESATI (T)
      QSI = EPS*ESI/(P-ESI)
!      WRITE(*,*) ' QSI ',QSI,' QVK ',QVK

!---- ADJUST TO SATURATION (SUBLIMATION) ------------
      IF(QVK/QSI-1.0_r8 .LE. 1.0D-6)  GO TO 1000

      Call Sublimation (T,P,QVK,QSI,ESI,QIK,ISTOP)

      GO TO 1000

!     ---------------------------------------------
!     ---------      SATURATED CASE       ---------
!     ---------------------------------------------

  500 IF(T .GT. TICE)  GO TO 700
      IF(T .GE. TTFRZ) GO TO 520
!     WRITE(6,*) 'SUPERSAT  T<-40  ICE ONLY'

!-- SUBLIMATE VAPOR UNTIL ICE SATURATION ----------
      IF(QVK/QSI-1.0_r8 .LE. 1.D-6)  GO TO 1000
      TQIK = 0.0_r8
      QVKSV = QVK
      TSV   = T

      Call Sublimation (T,P,QVK,QSI,ESI,TQIK,ISTOP)

!     WRITE(6,*) 'TEMP AFTER VAP-->ICE ',T

      IF(T .GT. TTFRZ)  GO TO 510

      QIK = QIK+TQIK
      GO TO 1000

!--- LATENT HEAT RELEASE PRODUCES INCONSISTENCY IN INITIAL ASSUMPTION -----
!       OF T <= TTFRZ.  THEREFORE SUBLIMATE ENOUGH VAPOR TO BRING
!       T = TTFRZ AND THEN PROCEED WITH CONDENSATION AND SUBLIMATION
!       FOR MIXED WATER AND ICE SATURATION CASE
  510 T = TTFRZ
!     WRITE(6,*) 'RE-DO VAP-->ICE UNTIL T=-40'
      DQI = (T-TSV)/CPLS
      QIK = QIK+DQI
      QVK = QVKSV-DQI
      QSW = EPS*ESW00/(P-ESW00)
      QSI = EPS*ESI00/(P-ESI00)
      QS  = (QLK*QSW+QIK*QSI)/(QLK+QIK)

      IF(ABS(QVK/QS-1.0_r8) .LE. 1.0D-6)  GO TO 1000
      IF(QVK .LT. QS)  GO TO 600

!---- CONDENSE AND SUBLIMATE VAPOR UNTIL MIXED WATER AND ICE SATURATION ----
!     CONDENSATION VERSUS SUBLIMATION DEPENDS ON TEMPERATURE.

!     METHOD IS AS FOLLOWS --

!       INITIAL STEP (TAG = 1) IS SUPERSATURATED (SS1 > 0)
!       SECOND GUESS (TAG = 2) IS MADE BY CONDENSING AND
!            SUBLIMATING SOME ARBITRARY AMOUNT OF VAPOR
!            (SST).  THIS GUESS MAY GIVE SS2 < SS1 (I.E.
!            CONVERGENCE).  IF SO, THIRD GUESS IS CALCU-
!            LATED ITERATIVELY AS EXPLAINED BELOW.  IF
!            SS2 > SS1 (I.E. NON-CONVERGENCE) MORE VAPOR
!            IS CONDENSED UNTIL CONVERGENCE OCCURS.
!       THIRD GUESS (TAG = 3) IS BY LINEAR INTERPOLATION IN
!            MIXING RATIO (QV3).
!       IF SS3 > 0, TAG = 1 IS REPLACED BY TAG = 3.
!       IF SS3 < 0, TAG = 2 IS REPLACED BY TAG = 3.
!---------------------------------------------------------------------------

  520 CRIT = 5.0D-6
!     WRITE(6,*) 'SUPERSAT:  COND VAPOR TO ICE/LIQ SAT '
      N   = 0
      ESI = ESATI (T)
      QSI = EPS*ESI/(P-ESI)
      SST = min(5.0_r8*(QVK-QS),QVK-QSI)
  530 QV1 = QVK
      QL1 = QLK
      QI1 = QIK
      SS1 = QVK-QS
      T1  = T
      QSIW1 = QS

      DQL = SST*MIN(MAX(T-TTFRZ,0.0_r8),TDIF)/TDIF
      DQI = SST*MIN(MAX(TICE-T,0.0_r8),TDIF)/TDIF
      QV2 = QV1-SST
      QL2 = QL1+DQL
      QI2 = QI1+DQI
      T2  = T1+DQL*CPLC+DQI*CPLS
      ESI = ESATI (T2)
      ESW = ESATW(T2)
      QSIT = EPS*ESI/(P-ESI)
      QSWT = EPS*ESW/(P-ESW)
      QSIW2 = (QL2*QSWT+QI2*QSIT)/(QL2+QI2)
      SS2 = QV2-QSIW2
      IF(SS2 .LT. 0.0_r8) GO TO 540

!-- CONDENSATION AND SUBLIMATION PRODUCE LARGER SUPERSATURATION  --
!       (TRY AGAIN)
  531 T   = T2
      QVK = QV2
      QLK = QL2
      QIK = QI2
      QS  = QSIW2
!     WRITE(6,11)  N,T1,T2,SS1,SS2,QV1,QV2,QSIW1,QSIW2,QL1,QL2,
!    1QI1,QI2,DQL,DQI
!   11 FORMAT('0'///'  SS2 > SS1 DIVERGENT CASE'/I6,2F14.7/(2E16.8/))
      IF(SS2/QSIW2 .LE. CRIT)  GO TO 570
      N = N+1
      IF(N .GE. 5)  GO TO 532
      GO TO 530

!-- IF (QVK-QSI) WAS CHOSEN, THERE APPEARS TO BE NO HOPE TO FIND
!     THE SOLUTION.  HOWEVER, IF 5.*(QVK-QS) IS CHOSEN, THERE IS
!     SOME SLIM HOPE THAT THE LARGER VALUE OF (QVK-QSI) WILL PRODUCE
!     A NEGATIVE SS2
  532 IF(SST .EQ. QVK-QSI)  GO TO 533
      SST = QVK-QSI
      GO TO 530

!-- THE SITUATION APPEARS HOPELESS ----------------
  533 ISTOP = 533
      GO TO 570

!-- ITERATION BEGINS ------------------------------
  540 QV3 = (QV2*SS1-QV1*SS2)/(SS1-SS2)
      DQ  = QV1-QV3
      DQL = DQ*MIN(MAX(T1-TTFRZ,0.0_r8),TDIF)/TDIF
      DQI = DQ*MIN(MAX(TICE-T1,0.0_r8),TDIF)/TDIF
      QL3 = QL1+DQL
      QI3 = QI1+DQI
      T3  = T1+DQL*CPLC+DQI*CPLS
      ESI = ESATI (T3)
      ESW = ESATW(T3)
      QSIT = EPS*ESI/(P-ESI)
      QSWT = EPS*ESW/(P-ESW)
      QSIW3 = (QL3*QSWT+QI3*QSIT)/(QL3+QI3)
      SS3 = QV3-QSIW3
!     WRITE(6,112) N,T1,T2,T3,SS1,SS2,SS3,QV1,QV2,QV3,QSIW1,QSIW2
!    1,QSIW3,QL1,QL2,QL3,QI1,QI2,QI3,DQL,DQI
!  112 FORMAT('0',I5,3F14.7/(3E16.8/))

      IF(DABS(SS3)/QSIW3 .LE. CRIT)  GO TO 560

      N = N+1
      IF(N .EQ. 30)  GO TO 565
      IF(SS3 .GT. 0.0_r8)  GO TO 550

      QV2 = QV3
      SS2 = SS3
      GO TO 540

  550 T1  = T3
      QV1 = QV3
      QL1 = QL3
      QI1 = QI3
      SS1 = SS3
      GO TO 540

  560 T = T3
      QVK = QV3
      QLK = QL3
      QIK = QI3
      QSIW = QSIW3
      GO TO 570

!---- NON-CONVERGENT CASE.  SITUATION IS AS FOLLOWS ------------------------
!       SOME SLW AND SIW HAVE BEEN PRODUCED BUT PRESENCE OF SLW
!       RAISES SATURATION VAPOR PRESSURE ABOVE INITIAL (TAG = 1) VALUE.
!       THEREFORE, EVAPORATE ALL SLW AND ITERATE TO SATURATION (I.E. GO TO #210)
  565 QVK = QV3
      QLK = QL3
      QIK = QI3
      T   = T3
!     WRITE(6,22)  T,QVK,QSIW3,QLK,QIK
!   22 FORMAT('0'///10X,'CONDENSATION AND SUBLIMATION TO MIXEDWATER AND ICE SATURATION FAILS'/5X,  &
!             'T',13X,'QVK',13X,'QSIW3',13X,'QLK',13X,'QIK'/F12.7,4E16.8)
      IF(SS3 .LT. 0.0_r8)  GO TO 210
      NOCONV = 1
      ISTOP = 565

  570 IF(T .LE. TICE)  GO TO 1000

!-- LATENT HEAT RELEASE PRODUCED INCONSISTENCY IN INITIAL ----------------
!       ASSUMPTION OF T <= TICE.  THEREFORE MELT ENOUGH ICE TO SET T = TICE
      DQLK = (T-TICE)/CPLF
!     WRITE(6,*) 'T>0  ICE AND WATER  SATD'
!     WRITE(6,*) 'DQLK,QIK',DQLK,QIK
      IF(DQLK .LE. QIK)  GO TO 580

!-- ALL ICE REMELTED BUT T > TICE -----------------
!     WRITE(6,*) ' ALL ICE REMELTED BUT T>0'
      DQLK = QIK
      QIK  = 0.0_r8
      T    = T-DQLK*CPLF
      ESW  = ESATW(T)
      QSW  = EPS*ESW/(P-ESW)
      QS   = QSW
      QLK  = QLK+DQLK
!     WRITE(6,571) T,DQLK,QS,QLK
!  571 FORMAT('0'///'   571',F14.7,3E16.8)

      IF(ABS(QVK/QS-1.0_r8) .LE. 1.0D-6)  GO TO 1000
      IF(QVK .GT. QS)  GO TO 700
      IF(QVK .LT. QS)  GO TO 210

!-- SOME ICE MELTS AND T = TICE -------------------
  580 T = TICE
!     WRITE(6,*) 'MELT ENOUGH ICE TO MAKE T=0'
      QLK = QLK+DQLK
      QIK = QIK-DQLK
      ESW = ESATW (T)
      QSW = EPS*ESW/(P-ESW)
      QS  = QSW
!     WRITE(6,581) T,QLK,QIK,QS
!  581 FORMAT('0'///'   581',F14.7,3E16.8)

      IF(ABS(QVK/QS-1.0_r8) .LE. 1.D-6)  GO TO 1000
      IF(QVK .LT. QS)  GO TO 210

!-- STILL SUPERSATURATED AT T = TICE WITH ICE PRESENT.
!       CONDENSE SUPERSATURATE AND USE HEAT TO MELT ICE
!       (T IS HELD CONSTANT).  IF THERE IS STILL EXCESS HEAT,
!       ADJUST TO SATURATION.
      DQV = QVK-QS
      DQI = DQV*HLTC/HLTF
!     WRITE(6,*) 'SUPERSAT  T=0  ICE AND LIQ '
!     WRITE(6,*) 'DQI,QIK',DQI,QIK
      IF(DQI .GT. QIK)  GO TO 590

!-- SOME ICE MELTS --------------------------------
!     WRITE(6,*) 'MELT ENOUGH ICE TO USE HEAT  SATD'
      QVK = QVK-DQV
      QIK = QIK-DQI
      QLK = QLK+DQV+DQI
!     WRITE(6,589) T,DQV,DQI,QVK,QIK,QLK
!  589 FORMAT('0'///'   589',F14.7,5E16.8)
      GO TO 1000

!-- ALL ICE MELTS ---------------------------------
  590 DQV = QIK*HLTF/HLTC
      QVK = QVK-DQV
      QLK = QLK+DQV+QIK
      QIK = 0.0_r8
!     WRITE(6,*) 'MELT ALL ICE TO USE HEAT  SUPERSAT'
!     WRITE(6,591) T,DQV,QVK,QIK,QLK
!  591 FORMAT('0'///'   591',F14.7,4E16.8)
      GO TO 700

  600 WRITE(6,601)  T,QVK,QS,QSW,QSI,QLK,QIK
  601 FORMAT('1',10X,'WILD AND CRAZY CASE'/F12.7,6E16.8)
      ISTOP = 600
      GO TO 2000

  700 IF(QVK/QS-1.0_r8 .LE. 1.D-6)  GO TO 1000
!     WRITE(6,*)' T>=0  LIQ ONLY  SUPERSAT'
      TQLK = 0.0_r8

      Call Condensation (T,P,QVK,QSW,ESW,TQLK,ISTOP)

      QLK = QLK+TQLK

 1000 CONTINUE
      IF(ISTOP .NE. 0)  NOCONV = 1
      IF (NOCONV.EQ.0) GO TO 1001
!     WRITE(*,*)'NOCONV 1000 ','I,J,K',IERR,JERR,KERR
!     WRITE(*,2001)  ISTOP,TO,P,QVO,QLO,QIO
 1001 CONTINUE
      RETURN

 2000 CONTINUE
      WRITE(6,2001)  ISTOP,TO,P,QVO,QLO,QIO
 2001 FORMAT('0','   PROBLEM',I8,' IN SAT'/15X,'T=',18X,'P=',            &
      17X,'QVO=',17X,'QLO=',17X,'QIO='/5E20.10)
      STOP
      END Subroutine Saturation

!=============================================================================
   SUBROUTINE CONDENSATION (T,P,QVK,QSW,ESW,QLK,ISTOP)
!=============================================================================
!  INPUT:
!     T           in (initial temperature to be modified in subroutine)
!     QVK         in (initial qv to be modified in subroutine)
!     QLK         in (initial qc to be modified in subroutine)
!     P           in (pressure, mb)
!     QSW         in (initial q*)
!     ESW         in (initial water vapor pressure)
!
!     EPS         in (constant =0.622)                  |
!     CPLC        in (constant =Lc/cp)                  |
!     TTD         in (constant)                          > Through Module mphysics_parameters
!     DESWT       in (constant array desat/dT table)    |
!     TT          in (constant)                         |
!     ESWT        in (constant array esat table)        |
!
!  INTERNAL:
!     CRIT, N, GAMFAC, SS         internal
!     GAMMAW, EX,  ESW            internal
!
!  OUTPUT:
!     T           out (adjusted temperature after condensation)
!     QVK         out (adjusted qv after condensation)
!     QLK         out (adjusted qc after condensation)
!     ISTOP       out (Error code)
!
!  CALLS:
!     FUNCTION DESWDT (T)
!     FUNCTION ESATW (T)
!---------------------------------------------------------------------------------------------

USE shr_kind_mod, only: r8 => shr_kind_r8
use mphysics_parameters, only: eswt,deswt,tt,ttd,cplc,eps

IMPLICIT NONE

    REAL(KIND=r8), INTENT(IN) :: P  ! pressure, mb
    REAL(KIND=r8), INTENT(INOUT) :: &
       T,    &  ! temperature to be modified
       QVK,  &  ! qv to be modified
       QLK,  &  ! qc to be modified
       ESW,  &  ! water vapor pressure
       QSW      ! initial q*

    INTEGER, INTENT(OUT) :: ISTOP

! local variables

    REAL (KIND=r8) :: GAMFAC,SS,GAMMAW,EX
    REAL (KIND=r8) :: ESATW,DESWDT
    real (kind=r8), Parameter :: CRIT=1.E-6_r8
    INTEGER :: n

    GAMFAC = EPS * CPLC * P

      GAMMAW = GAMFAC * DESWDT(T) / (P-ESW)**2
      EX     = (QVK - QSW) / (1.0_r8 + GAMMAW)
      T      = T + CPLC * EX
      ESW    = ESATW(T) ;    QSW = EPS * ESW / (P-ESW)
      QVK    = QVK - EX ;    QLK = QLK + EX

  do n=1,9

      GAMMAW = GAMFAC * DESWDT(T) / (P-ESW)**2
      EX     = (QVK - QSW) / (1.0_r8 + GAMMAW)

    IF (DABS(EX/QLK) .LT. CRIT)  RETURN

      T   = T + CPLC * EX
      ESW = ESATW(T) ;   QSW = EPS * ESW / (P-ESW)
      QVK = QVK - EX ;   QLK = QLK + EX

  enddo

      ISTOP = 1111

RETURN
END SUBROUTINE Condensation

!=============================================================================
   SUBROUTINE SUBLIMATION (T,P,QVK,QSI,ESI,QIK,ISTOP)
!=============================================================================
!  INPUT:
!     T           in (initial temperature to be modified in subroutine)
!     QVK         in (initial qv to be modified in subroutine)
!     QIK         in (initial qi to be modified in subroutine)
!     P           in (pressure, mb)
!     QSI         in (initial q*)
!     ESI         in (initial water vapor pressure over ice)
!
!     EPS         in (constant =0.622)               |
!     CPLS        in (constant =Ls/cp)               |
!     TTD         in (constant)                       > Through Module physics_paramaters
!     DESIT       in (constant array desat/dT table) |
!     TT          in (constant)                      |
!     ESIT        in (constant array esat table)     |
!
!  INTERNAL:
!     CRIT, N, GAMFAC, SS        internal
!     GAMMAI, EX, ESI            internal
!
!  OUTPUT:
!     T           out (adjusted temperature after condensation)
!     QVK         out (adjusted qv after sublimation)
!     QIK         out (adjusted qc after sublimation)
!     ISTOP       out (Error code)
!
!  CALLS:
!     FUNCTION DESIDT(T)
!     FUNCTION ESATI (T)
!------------------------------------------------------------------------------------------

USE shr_kind_mod, only: r8 => shr_kind_r8
use mphysics_parameters, only: esit,desit,tt,ttd,cpls,eps

IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: P  ! pressure, mb
    REAL (KIND=r8), INTENT(INOUT) :: &
       T,     &  ! temperature to be modified
       QVK,   &  ! qv to be modified
       QSI,   &  ! q*
       ESI,   &  ! water vapor pressure over ice
       QIK       ! qi to be modified

    INTEGER, INTENT(OUT) :: ISTOP

! local variables
    real(KIND=r8) :: GAMFAC,SS,GAMMAI,EX
    real(KIND=r8) :: ESATI,DESIDT

    REAL (kind=r8), Parameter :: CRIT=0.5E-7_r8
    INTEGER :: n

    GAMFAC = EPS * CPLS * P

      GAMMAI = GAMFAC * DESIDT(T) / (P-ESI)**2
      EX     = (QVK - QSI) / (1.0_r8 + GAMMAI)
      T   = T + CPLS * EX
      ESI = ESATI (T) ;  QSI = EPS * ESI / (P-ESI)
      QVK = QVK - EX  ;  QIK = QIK + EX

  do n=1,9

      GAMMAI = GAMFAC * DESIDT(T) / (P-ESI)**2
      EX     = (QVK - QSI) / (1.0_r8 + GAMMAI)

    IF(DABS(EX/QIK) .LT. CRIT)  RETURN

      T   = T + CPLS * EX
      ESI = ESATI (T) ;  QSI = EPS * ESI / (P-ESI)
      QVK = QVK - EX  ;  QIK = QIK + EX

  enddo

      ISTOP = 2222

RETURN
END SUBROUTINE Sublimation

!=============================================================================
    Subroutine MICROPHYSICS
!=============================================================================
!  INPUT:
!     theta,qv,qc,qi,qr,qs,qg
!     pl0,pil0,rho_int,rhol,dzl,hx_sub
!
!     INTERNAL:
!        T                   temperature entering to Imicro
!        QMIX(1) = qv    |
!        QMIX(2) = qc    |
!        QMIX(3) = qi     >  mixing ratios entering to Imicro
!        QMIX(4) = qs    |
!        QMIX(5) = qg    |
!        QMIX(6) = qr    |
!        TSS(1)       ->  tendency_microphysics_qv     |
!        TSS(2)       ->  tendency_microphysics_qc     |
!        TSS(3)       ->  tendency_microphysics_qi     |
!        TSS(4)       ->  tendency_microphysics_qs      > microphysics tendencies
!        TSS(5)       ->  tendency_microphysics_qg     |     returned from Imicro
!        TSS(6)       ->  tendency_microphysics_qr     |
!        TSS(7)/pil0  ->  tendency_microphysics_theta  |
!        PSS(26)      ->  conversion rates between the species returned from Imicro
!        VTERMS,VTERMG,VTERMR -> terminal velocities returned from Imicro
!        VTR, VTS, VTG           terminal velocities and fluxes
!        LMICRO                  logical
!        THRESH
!        PRESS_mb
!        TLHR         ->  Latent heating rate from all processes except condensational
!                         and sublimational saturation adjustment
!     OUTPUT:
!        tendency_microphysics_theta,tendency_microphysics_qv
!        tendency_microphysics_qc,tendency_microphysics_qi
!        tendency_microphysics_qr,tendency_microphysics_qs
!        tendency_microphysics_qg
!        tendency_rain,tendency_snow, tendency_graupel  (Adams-Bashforth tendencies of
!                                                        precipitation                 )
!        Surface_rain, Surface_snow, Surface_graupel
!        VTR_int,VTS_int,VTG_int
!        latent_heating_rate
!------------------------------------------------------------------------------------------
USE shr_kind_mod,        only: r8 => shr_kind_r8
USE constld,             only: localp
USE mphysics_parameters, only: im,km
USE mphysics_variables,  only: rho_int,rhol,dzl,pl0,pil0,   &
                               theta,qv,qc,qi,qr,qs,qg,     &
                               tendency_microphysics_theta, &
                               tendency_microphysics_qv,    &
                               tendency_microphysics_qc,    &
                               tendency_microphysics_qi,    &
                               tendency_microphysics_qr,    &
                               tendency_microphysics_qs,    &
                               tendency_microphysics_qg,    &
                               tendency_rain,tendency_snow,tendency_graupel, &
                               vtr_int,vts_int,vtg_int,     &
                               surface_rain,surface_snow,surface_graupel, &
                               latent_heating_rate,hx_sub

USE mphysics_variables, only: rho_int_loc,rhol_loc,pl0_loc,pil0_loc       ! JUNG_LocalP

IMPLICIT NONE
      
      ! LOCAL
      REAL (kind=r8) :: T, THRESH, PRESS_mb, VTERMS, VTERMG,  &
                              VTERMR, QMIX(6), TSS(7), PSS(26)
      REAL (kind=r8) :: VTR(km), VTS(km), VTG(km)   
      REAL (kind=r8) :: TLHR
      
      REAL (kind=r8) :: rho_int_v(0:km),rhol_v(km),pl0_v(km),pil0_v(km)   ! JUNG_LocalP

      LOGICAL lmicro
      INTEGER :: I, K, KLOW

      LMICRO = .True.
      THRESH = 1.0D-6

!--------------------------------------
      DO I = 1, IM
!--------------------------------------
      
      klow = hx_sub(i)
      
      IF (LocalP) THEN
       pl0_v(:)     = pl0_loc(I,:)
       pil0_v(:)    = pil0_loc(I,:) 
       rhol_v(:)    = rhol_loc(I,:)
       rho_int_v(:) = rho_int_loc(I,:)
       
      ! Initialize Local Parameters for Microphysics
      CALL initialize_microphysics_loc (rho_int_v(0))

      ELSE
       pl0_v(:)     = pl0(:)
       pil0_v(:)    = pil0(:) 
       rhol_v(:)    = rhol(:)
       rho_int_v(:) = rho_int(:)
      ENDIF

      DO K = klow, km
      T       = pil0_v(K) * theta(I,K)
      QMIX(1) = QV(I,K);  QMIX(2) = QC(I,K);  QMIX(3) = QI(I,K)
      QMIX(4) = QS(I,K);  QMIX(5) = QG(I,K);  QMIX(6) = QR(I,K)

      LMICRO  =  QMIX(2) .GT. THRESH .OR. QMIX(3) .GT. THRESH .OR.   &
                 QMIX(4) .GT. THRESH .OR. QMIX(5) .GT. THRESH .OR.   &
                 QMIX(6) .GT. THRESH

      IF ( LMICRO ) THEN
        PRESS_mb = pl0_v(K) * 0.01_r8
        
        Call Imicro (TSS,PSS,VTERMS,VTERMG,VTERMR,T,PRESS_mb,QMIX,rhol_v(k),TLHR)

        tendency_microphysics_QV(I,K)    = TSS(1)
        tendency_microphysics_QC(I,K)    = TSS(2)
        tendency_microphysics_QI(I,K)    = TSS(3)
        tendency_microphysics_QS(I,K)    = TSS(4)
        tendency_microphysics_QG(I,K)    = TSS(5)
        tendency_microphysics_QR(I,K)    = TSS(6)
        tendency_microphysics_theta(I,K) = TSS(7) / pil0_v(K)

        VTS(K) = VTERMS
        VTG(K) = VTERMG
        VTR(K) = VTERMR

        latent_heating_rate(i,k) = TLHR
      ELSE
        VTS(K) = 0.0_r8
        VTG(K) = 0.0_r8
        VTR(K) = 0.0_r8
      ENDIF

      ENDDO  ! k-loop
      
      ! Rain, snow, graupel and sedimentation fluxes
      VTS_int(i,km) = 0.0_r8
      VTG_int(i,km) = 0.0_r8
      VTR_int(i,km) = 0.0_r8

      VTS_int(i,klow-1) = rhol_v(klow) * VTS(klow) * QS(I,klow)
      VTG_int(i,klow-1) = rhol_v(klow) * VTG(klow) * QG(I,klow)
      VTR_int(i,klow-1) = rhol_v(klow) * VTR(klow) * QR(I,klow)

      do K = km-1,klow,-1
        VTS_int(i,K) = rho_int_v(k) * 0.5_r8 * (VTS(K+1)+VTS(K)) * QS(I,K+1)
        VTG_int(i,K) = rho_int_v(k) * 0.5_r8 * (VTG(K+1)+VTG(K)) * QG(I,K+1)
        VTR_int(i,K) = rho_int_v(k) * 0.5_r8 * (VTR(K+1)+VTR(K)) * QR(I,K+1)
      enddo

      do K = klow,km
        tendency_rain(I,K)    = tendency_rain(I,K)                                    &
                                      +(VTR_int(i,K)-VTR_int(i,K-1))/(rhol_v(k)*dzl(K))
        tendency_snow(I,K)    = tendency_snow(I,K)                                    &
                                      +(VTS_int(i,K)-VTS_int(i,K-1))/(rhol_v(k)*dzl(K))
        tendency_graupel(I,K) = tendency_graupel(I,K)                                 &
                                      +(VTG_int(i,K)-VTG_int(i,K-1))/(rhol_v(k)*dzl(K))
      enddo
      do K = 1,klow-1
        tendency_rain(I,K)    = 0.0_r8
        tendency_snow(I,K)    = 0.0_r8
        tendency_graupel(I,K) = 0.0_r8
      enddo

      ! Surface rain
      Surface_rain(I)    = VTR_int(I,klow-1)
      Surface_snow(I)    = VTS_int(I,klow-1)
      Surface_graupel(I) = VTG_int(I,klow-1)

!--------------------------------------
      ENDDO   ! i-loop
!--------------------------------------

RETURN
END Subroutine microphysics

!=============================================================================
    SUBROUTINE IMICRO (TSS,PSS,VTS,VTG,VTR,T,PRESR,QMIX,RHO,TLHR)
!=============================================================================
!---------------------------------------------------------------------------
! MODIFICATIONS BELOW ARE MADE TO THE ORIGINAL CODE (Q. Fu & S. Krueger)
!---------------------------------------------------------------------------
!	Updates (Jun. 8, 1993, cfu comments):
!	  1. RMINUC = QI * RHO / NC;  this  will give the depositionalTSS
!            growth of cloud ice at expense of cloud water  properly (
!            PIDW).
!         2. Replacing R50 by R100; then DT1 is the growth time needed
!            for ice crystal to grow from 40 to 100um (PSFI and PSFW).
!         3. Setting Psfi = 0.0 when QL is less than 5e-6 (PSFI).
!         4. A1 and A2 are obtained by using interpolation; A1T(14) is
!            0.1725e-9 for function IDW ( PSFI, PSFW and PIDW).
!         5. Setting QI0 = .6e-3 (threshold for cloud ice aggregation)
!            and QS0 = 1.0e-3 (threshold for snow aggregation ) (PSAUT
!            and PGAUT) (private communication with S. Chin).
!         6. GAM625 = 184.860962 (PWACS)
!         7. Esi and Egs ( T < T0 ) are set to be  0.1 following  RH84
!            ( PSACI and PGACS ).
!---------------------------------------------------------------------------
!       Updates (Jun. 9, 1993, crh comments):
!         1. Modifying intercept parameters and densities.
!         2. Replacing hail fall speed relation (by Lin) with an average
!            graupel relation (Ug=40.74*Dg**0.5) based on Locatelli &
!            Hobbs (1974).
!---------------------------------------------------------------------------
!       Updates (Jun. 10, 1993, clin comments):
!         1. Correcting C2BRG.
!         2. Egs = 1.0 for T > T0 and wet growth. Note:  modifications
!            here should be consistent to (7) under the comment 'cfu'.
!---------------------------------------------------------------------------
!       Updates (Jun. 10, 1993, csk comments):
!         1. Lin et al. (1983) version is  installed for the following
!            processes: PSMLT, PGWET, PGMLT (constants only).
!            Change DIMENSION CSMLT(3),CGWET(2) => CSMLT(5), CGWET(4).
!         2. PWACS = 0.0 following Lin et al. (1983).
!         3. And more.
!---------------------------------------------------------------------------
!    INPUT:
!       T                      in (Temperature, K)
!       PRESR                  in (Pressure, mb)
!       QMIX(1) = qv    |
!       QMIX(2) = qc    |
!       QMIX(3) = qi     >  mixing ratios entering to Imicro
!       QMIX(4) = qs    |
!       QMIX(5) = qg    |
!       QMIX(6) = qr    |
!       RHO                    in (density)
!       CSACW, CWACS,CIACR,CRACS              |
!       CSACR,CGACR,CGACS,ACCO(3,4)            >        in (constant from initial_physics)
!       CRACI,CSACI,CGACW,CGACI,CGSUB,CREVP   |
!       CGFR,CSMLT,CGMLT,CGWET,CRAUT          |
!       QI0, QS0, QL0, ES0, CES0, C1BRG                 in
!       C2BRG, RMI50, RMI40, ESWT(151)                  in
!       DESWT(151),ESIT(151),DESIT(151)                 in
!       TT,CH2O,HLTS,HLTC,HLTF,CICE                     in
!       TICE,CP,EPS,VCONR,VCONS,VCONG                   in
!       TDT, RHOS                                       in
!
!   INTERNAL:
!       THRSH(6)          internal  (constant, threshold of mixing ratio of species)
!       QLOGIC(6)         internal (logical, defining existence of species exceeding a threshold)
!       QVK,QLK,QIK,QSK
!       QGK,QRK,QVR,QLR
!       QIR,QSR,QGR,QRR
!       VAPR,LIQ,ICE,SNO
!       GRAUP,RAIN
!       TSSV,TSSL,TSSI
!       TSSS,TSSG,TSSR,TTMP
!
!       SVMAX=SMAX(1),SLMAX=SMAX(2)
!       SIMAX=SMAX(3),SSMAX=SMAX(4)
!       SGMAX=SMAX(5),SRMAX=SMAX(6)
!
!       C3RACS(1)=ACCO(1,1),C3RACS(2)=ACCO(2,1),C3RACS(3) = ACCO(3,1)
!       C3SACR(1)=ACCO(1,2),C3SACR(2)=ACCO(2,2),C3SACR(3) = ACCO(3,2)
!       C3GACR(1)=ACCO(1,3),C3GACR(2)=ACCO(2,3),C3GACR(3) = ACCO(3,3)
!       C3GACS(1)=ACCO(1,4),C3GACS(2)=ACCO(2,4),C3GACS(3) = ACCO(3,4)
!       ARC1,ARC2,ARC3,QEXRHO
!       QL0RHO,RAUT,RHOFAC
!       CTGACS,DQS,DQS0,QSW,SW
!       REVP,TSQ,DQS0,SMLT,GMLT
!       QIKT,EXPT1,EXPT2,C,AUT
!       GFR,ESI,SI,SSUB,GSUB
!       PGDRY
!
!   OUTPUTS:
!     Primary outputs are tendencies and terminal velocity of rain, snow and graupel
!
!       TSS(1)    = TSSV    : tendency_microphysics_qv
!       TSS(2)    = TSSL    : tendency_microphysics_qc
!       TSS(3)    = TSSI    : tendency_microphysics_qi
!       TSS(4)    = TSSS    : tendency_microphysics_qs
!       TSS(5)    = TSSG    : tendency_microphysics_qg
!       TSS(6)    = TSSR    : tendency_microphysics_qr
!       TSS(7)/pi = TTMP    : tendency_microphysics_theta
!
!
!       VTR       : terminal velocity of rain
!       VTS       : terminal velocity of snow
!       VTG       : terminal velocity of graupel
!
!     Secondary output is source/sinks of microphysical species (for diagnostic purposes)
!
!      PSS(1)  = PRAUT; PSS(2)  = PRACW;  PSS(3)  = PRACS; PSS(4)  = PSACW; PSS(5)  = PWACS
!      PSS(6)  = PGACW; PSS(7)  = PGMLT;  PSS(8)  = PSMLT; PSS(9)  = PREVP; PSS(10) = PIACR
!      PSS(11) = PSACR; PSS(12) = PGACR;  PSS(13) = PGFR;  PSS(14) = PGACS; PSS(15) = PGAUT
!      PSS(16) = PGACI; PSS(17) = PGWORD; PSS(18) = PRACI; PSS(19) = PSAUT; PSS(20) = PSACI
!      PSS(21) = PSSUB; PSS(22) = PGSUB;  PSS(23) = PSFW;  PSS(24) = PSFI;  PSS(25) = PIDW
!      PSS(26) = PGWET
!
!  CALLS:
!      FUNCTION ESATI(T)
!      FUNCTION ESATW(T)
!      SUBROUTINE Cal_ACR3 (ACR3,V1,V2,QRHO1,QRHO2,C,CAC,rho)
!      SUBROUTINE Cal_IDW (IDW,TC,QIK,RHO)
!      Subroutine Bergeron (PSFW,PSFI,TC,QL,QI,QLRHO)
!---------------------------------------------------------------------------------------------------

      USE shr_kind_mod, only: r8 => shr_kind_r8
      USE mphysics_parameters, RHOS=>RHOSFC

      IMPLICIT NONE

      REAL (KIND=r8), INTENT(IN) ::     &
         T,          &
         PRESR,      &
         QMIX(6),    &
         RHO
      REAL (KIND=r8), INTENT(OUT) ::     &
         VTS,        &
         VTG,        &
         VTR,        &
         TSS(7),     &
         PSS(26),    &
         TLHR

! local variables
      REAL (KIND=r8) :: RHOFAC, QL0RHO, TC, TSQ, EXPT1, EXPT2, ESW, ESATW
      REAL (KIND=r8) :: QSW, DQS, DQS0
      REAL (KIND=r8) :: SVMAX,SLMAX,SIMAX,SSMAX,SGMAX,SRMAX
      REAL (KIND=r8) :: TSSV,TSSL,TSSI,TSSS,TSSG,TSSR,TTMP
      REAL (KIND=r8) :: QVK,QLK,QIK,QSK,QGK,QRK,QVR,QLR,QIR,QSR,QGR,QRR

      REAL (KIND=r8) :: C3RACS(3),C3SACR(3),C3GACR(3),C3GACS(3),QRHO(6),SMAX(6)
      REAL (KIND=r8) :: PRAUT,PRACW,PRACS,PSACW,PWACS,PGACW,PGMLT,PSMLT,PREVP,                       &
             PIACR,PSACR,PGACR,PGFR,PGACS,PGAUT,PGACI,PGWORD,PRACI,PSAUT,                 &
             PSACI,PSSUB,PGSUB,PSFW,PSFI,PIDW,PGWET
      REAL (KIND=r8) :: ACR1,ACR2,QEXRHO,RAUT,ACR3,CTGACS,SW,REVP,SMLT,GMLT,TEM,                     &
             QIKT,C,AUT,GFR,ESI,ESATI,QSI,SI,SSUB,GSUB,PGDRY,GWET,IDW,X,Y

      LOGICAL QLOGIC(6),VAPR,LIQ,ICE,SNO,GRAUP,RAIN
      INTEGER :: k   ! do loop index

      REAL (KIND=r8), DIMENSION(6) :: thrsh

!     LOCAL VARIABLES: TEMPERATURE, DENSITY, TEMP. DEPENDENT EFFICIENCIES
!     AND SATURATION VARIABLES

      thrsh = (/0.0_r8,1.0D-6,1.0D-6,1.0D-6,1.0D-6,1.0D-6/)

      QVK = QMIX(1); QLK = QMIX(2); QIK = QMIX(3); QSK = QMIX(4)
      QGK = QMIX(5); QRK = QMIX(6)

      C3RACS(1) = ACCO(1,1);  C3RACS(2) = ACCO(2,1);  C3RACS(3) = ACCO(3,1)
      C3SACR(1) = ACCO(1,2);  C3SACR(2) = ACCO(2,2);  C3SACR(3) = ACCO(3,2)
      C3GACR(1) = ACCO(1,3);  C3GACR(2) = ACCO(2,3);  C3GACR(3) = ACCO(3,3)
      C3GACS(1) = ACCO(1,4);  C3GACS(2) = ACCO(2,4);  C3GACS(3) = ACCO(3,4)

      RHOFAC = SQRT(RHOS/RHO)
      QL0RHO = QL0 * RHO
      TC     = T-TICE
      TSQ    = T**2

      EXPT1  = min(EXP(0.090_r8*TC),1.0_r8)
      EXPT2  = min(EXP(0.025_r8*TC),1.0_r8)

      ESW    = ESATW(T)
      QSW    = EPS*ESW/(PRESR-ESW)
      DQS    = QVK-QSW
      DQS0   = CES0/(PRESR-ES0)-QVK

!----- ZERO SOURCE AND SINK TERMS--------------------------

      PSS = 0.0_r8

      PRAUT = PSS(1);  PRACW  = PSS(2);  PRACS = PSS(3);  PSACW = PSS(4);  PWACS = PSS(5)
      PGACW = PSS(6);  PGMLT  = PSS(7);  PSMLT = PSS(8);  PREVP = PSS(9);  PIACR = PSS(10)
      PSACR = PSS(11); PGACR  = PSS(12); PGFR  = PSS(13); PGACS = PSS(14); PGAUT = PSS(15)
      PGACI = PSS(16); PGWORD = PSS(17); PRACI = PSS(18); PSAUT = PSS(19); PSACI = PSS(20)
      PSSUB = PSS(21); PGSUB  = PSS(22); PSFW  = PSS(23); PSFI  = PSS(24); PIDW  = PSS(25)
      PGWET = PSS(26)

!----- ZERO OUTPUT TENDENCIES -----------------------------

       TSS = 0.0_r8

       TSSV  = TSS(1);    TSSL  = TSS(2);    TSSI  = TSS(3)
       TSSS  = TSS(4);    TSSG  = TSS(5);    TSSR  = TSS(6);   TTMP = TSS(7)

!----------------------------------------------------------
!     DEFINE MIXING RATIOS GREATER THAN THRESHOLD FOR
!     ACCRETION PROCESSES AND MASS DENSITIES

!           1:  WATER VAPOR                     (=qv)
!           2:  CLOUD (SUSPENDED LIQUID) WATER  (=qc)
!           3:  CLOUD ICE CRYSTALS              (=qi)
!           4:  SNOW                            (=qs)
!           5:  GRAUPEL                         (=qg)
!           6:  RAIN                            (=qr)

do K=1,6
  QLOGIC(K) = QMIX(K) .GT. THRSH(K)
  SMAX(K)   = QMIX(K) / TDT
  QRHO(K)   = QMIX(K) * RHO
enddo

!-------------------------------------------------------------------------------------

      VAPR  = QLOGIC(1); LIQ   = QLOGIC(2); ICE   = QLOGIC(3)
      SNO   = QLOGIC(4); GRAUP = QLOGIC(5); RAIN  = QLOGIC(6)

      SVMAX = SMAX(1);   SLMAX = SMAX(2);   SIMAX = SMAX(3)
      SSMAX = SMAX(4);   SGMAX = SMAX(5);   SRMAX = SMAX(6)

      QVR   = QRHO(1);   QLR   = QRHO(2);   QIR   = QRHO(3)
      QSR   = QRHO(4);   QGR   = QRHO(5);   QRR   = QRHO(6)

      VTS   = 0.0_r8
      VTG   = 0.0_r8
      VTR   = 0.0_r8

!----- TERMINAL VELOCITIES ------------------------------------------------

      IF(SNO)    VTS =       VCONS * QSR**0.0625 * RHOFAC
      IF(GRAUP)  VTG = min ( VCONG * QGR**0.125 / DSQRT(RHO) , 20.0_r8 )
      IF(RAIN)   VTR = min ( VCONR * QRR**0.2 * RHOFAC , 10.0_r8 )

!     ---------------------------------------------------------------------
!     ----   PROCESSES CALLED INDEPENDENT OF   ----
!     ----   TEMPERATURE ARE:                  ----
!     ----                                     ----
!     ----        GACW    SACW    RAUT         ----
!     ----        RACW    RACS    GACS         ----
!     ----                REVP                 ----
!     ---------------------------------------------

      IF(.NOT. LIQ)    GOTO 150
      IF(.NOT. GRAUP)  GOTO 110

!     PGACW: Accretion of cloud water by graupel.

      ACR2  = CGACW * QLK * QGR**0.875 / DSQRT(RHO)
      PGACW = min(ACR2,SLMAX)

  110 IF(.NOT. SNO)  GOTO 120

!     PSACW: Accretion of cloud water by snow.  Produces snow if T<TICE and rain if T>=TICE.
!            Enhanced snow melting if T>=TICE.
        ACR1  = CSACW * QLK * QSR**0.8125 * RHOFAC
        PSACW = min(ACR1,SLMAX)

  120 IF(QLK .LE. QL0)  GOTO 130

!     PRAUT: Autoconversion of cloud water to form rain.
!     SEE ORVILLE AND KOPP (1977)
        QEXRHO = QLR - QL0RHO
        RAUT   = QEXRHO**3 / ((CRAUT(1)*QEXRHO+CRAUT(2))*RHO)

        PRAUT  = min(RAUT,SLMAX)

  130 IF(.NOT. RAIN)  GOTO 200

!     PRACW: Accretion of cloud water by rain.
        ACR1  = CRACW * QLK * QRR**0.95 * RHOFAC
        PRACW = min(ACR1,SLMAX)

  150 IF (.NOT. RAIN) GOTO 200
  170 IF (.NOT. SNO)  GOTO 290

!     PRACS (Accretion of snow by rain. Produces graupel if rain and snow exceeds threshold
!            and T<TICE)

        Call Cal_ACR3(ACR3,VTR,VTS,QSR,QRR,CRACS,C3RACS,RHO)
        PRACS = min(ACR3,SSMAX)

      GOTO 210

  200 IF (.NOT. SNO)   GOTO 290
  210 IF (.NOT. GRAUP) GOTO 290

!     PGACS: Accretion of snow by graupel.
        CTGACS = CGACS*0.1_r8

!     PGACS: Accretion of snow by graupel.
        Call cal_ACR3(ACR3,VTG,VTS,QSR,QGR,CTGACS,C3GACS,rho)
        PGACS = min(ACR3,SSMAX)

  290 IF(QRK .EQ. 0.0_r8 .OR. ICE .OR. DQS .GE. 0.0_r8)  GO TO 300

!     PREVP: Evaporation of rain.
        SW = QVK/QSW
        REVP = CREVP(1)*TSQ*QSW*(1.0_r8-SW)*(CREVP(2)*DSQRT(QRR)+CREVP(3)*QRR**0.725) &
             /(CREVP(4)*TSQ+CREVP(5)*QSW*RHO)
        PREVP = min(REVP,SRMAX)

  300 IF (TC .LT. 0.0_r8)  GOTO 400


!     -----------------------------------
!     ----     TC >= 0 PROCESSES     ----
!     ----                           ----
!     ----    SMLT   WACS   GMLT     ----
!     -----------------------------------

      IF(QSK .EQ. 0.0_r8) GO TO 305

      IF ( SNO .AND. RAIN ) THEN

!     PSACR CALLED FOR SMLT HERE (T < 0 PROCESS) (Accretion of rain by snow.)
        Call cal_ACR3(ACR3,VTS,VTR,QRR,QSR,CSACR,C3SACR,rho)
        PSACR = min(ACR3,SRMAX)
      END IF

!   PSMLT: Melting of snow to form rain if T>=TICE. Follow Lin et al. (1983)
        SMLT = (CSMLT(1)*TC/RHO-CSMLT(2)*DQS0) * (CSMLT(3)*DSQRT(QSR)+                         &
                CSMLT(4)*QSR**0.65625*SQRT(RHOFAC)) + CSMLT(5)*TC*(PSACW+PSACR)

        PSMLT = MIN(SMLT,SSMAX)
        PSMLT = MAX(PSMLT,0.0_r8)
        PSACR = 0.0_r8                         ! Not used except in PSMLT

  305 IF(.NOT. LIQ .OR. .NOT. SNO)  GO TO 310

!     PWACS:
        PWACS = 0.0_r8

  310 IF(.NOT. GRAUP .OR. .NOT. RAIN)  GO TO 320

!     PGACR CALLED FOR PGMLT HERE (Accretion of rain by graupel.)
        Call cal_ACR3(ACR3,VTG,VTR,QRR,QGR,CGACR,C3GACR,RHO)
        PGACR = min(ACR3,SRMAX)

!     GMLT:  GACW HAS ALREADY BEEN CALLED IF APPROPRIATE
!            GUARD AGAINST NEGATIVE VALUES AT TEMP CLOSE TO 0 DEG C

  320 IF(QGK .EQ. 0.0_r8)  GOTO 330

!     PGMLT: Melting of graupel to form rain.
        GMLT = (CGMLT(1)*TC/RHO-CGMLT(2)*DQS0)*(CGMLT(3)*SQRT(QGR)  &
               +CGMLT(4)*QGR**0.6875/RHO**0.25)+CGMLT(5)*TC*(PGACW+PGACR)
        PGMLT = MIN(GMLT,SGMAX)
        PGMLT = MAX(PGMLT,0.0_r8)

!     PGACR   (Accretion of rain by graupel.)
        PGACR = 0.0_r8                        ! Not used except in PGMLT

!     PGACS: Accretion of snow by graupel.
        PGACS = min(10.0_r8*PGACS,SSMAX)

!     --------------------------------------
!     ----     ADD SOURCES AND SINKS    ----
!     ----             T>=0             ----
!     --------------------------------------

  330 CONTINUE

        TSSV = PREVP
        TSSL = -(PRAUT+PRACW+PSACW+PGACW)
        TSSS = -(PGACS+PRACS+PWACS+PSMLT)
        TSSG = PGACS-PGMLT
        TSSR = PRAUT+PRACW+PRACS+PSACW+PWACS+PGACW+PSMLT+PGMLT-PREVP
        TTMP = (-HLTF*(PRACS+PWACS+PGMLT+PSMLT)-HLTC*PREVP)/CP
        
      TSS(1) = TSSV;  TSS(2) = TSSL;  TSS(3) = TSSI
      TSS(4) = TSSS;  TSS(5) = TSSG;  TSS(6) = TSSR; TSS(7) = TTMP

     TEM = 0.0_r8
   do K=1,6
     TEM = TEM + TSS(K)
   enddo
      GOTO 1000

!     ------------------------------------
!     ----      TC < 0 PROCESSES      ----
!     ----                            ----
!     ----   IDW     BERGRN   RACI    ----
!     ----   GACI    SACI     SAUT    ----
!     ----   GFR     SACR     GACR    ----
!     ----   SSUB    GAUT     GSUB    ----
!     ----           WORD             ----
!     ----   NOTE: IACR IS NOT USED   ----
!     ----                            ----
!     ------------------------------------

  400 IF(.NOT. LIQ)  GO TO 410

!     PIDW:  Deposition on cloud ice at expense of cloud water.
        Call Cal_IDW (IDW,TC,QIK,RHO)
        PIDW = min(IDW,SLMAX)

  410 IF(.NOT. ICE)  GO TO 450

!     BERGERON PROCESSES -- PSFW (Transfer of cloud water to form snow) and
!                           PSFI (Transfer of cloud ice to form snow)
      Call BERGERON (PSFW,PSFI,TC,QLK,QIK,QLR)

      PSFW = min(PSFW,SLMAX)
      PSFI = min(PSFI,SIMAX)

      IF(.NOT. RAIN)  GO TO 420

!     PRACI: Accretion of cloud ice by rain.
        ACR1  = CRACI * QIK * QRR**0.95 * RHOFAC
        PRACI = min(ACR1,SIMAX)

  420 IF(.NOT. GRAUP)  GO TO 430

!     PGACI: Accretion of cloud ice by graupel.
        ACR2  = CGACI * QIK * QGR**0.875 / SQRT(RHO)
        PGACI = min(ACR2,SIMAX)

  430 IF(.NOT. SNO)  GO TO 440

!     PSACI: Accretion of cloud ice by snow.
        QIKT  = QIK*0.1_r8
        ACR1  = CSACI * QIKT * QSR**0.8125 * RHOFAC
        PSACI = min(ACR1,SIMAX)

  440 IF(QIK .LE. QI0)  GO TO 450

!     PSAUT: Autoconversion of cloud ice to form snow.
        C     = 0.1_r8*EXPT2
        AUT   = MAX(C*(QIK-QI0),0.0_r8)
        PSAUT = MIN(AUT,SIMAX)

  450 IF(.NOT. RAIN)  GO TO 470

!     PGFR: Probabilistic freezing of rain to form graupel.
        GFR  = CGFR(1)*(EXP(-CGFR(2)*TC)-1.0_r8)*QRR**1.75/RHO
        PGFR = min(GFR,SRMAX)

      IF(.NOT. SNO)  GO TO 460

!     PSACR: Acceretion of rain by snow.
        Call cal_ACR3(ACR3,VTS,VTR,QRR,QSR,CSACR,C3SACR,rho)
        PSACR = min(ACR3,SRMAX)

  460 IF(.NOT. GRAUP)  GO TO 470

!     PGACR: Acceretion of rain by graupel.
        Call cal_ACR3(ACR3,VTG,VTR,QRR,QGR,CGACR,C3GACR,rho)
        PGACR = min(ACR3,SRMAX)

  470 ESI = ESATI(T)
      QSI = EPS*ESI/(PRESR-ESI)
      SI = QVK/QSI

      IF(QSK .EQ. 0.0_r8)  GO TO 480

!     PSSUB: CAN BE EITHER POSITIVE OR NEGATIVE (Sublimation of snow.)
        SSUB  = CSSUB(1)*TSQ*QSI*(SI-1.0_r8) * (CSSUB(2)*SQRT(QSR)    &
              + CSSUB(3)*QSR**0.65625*SQRT(RHOFAC)) / (CSSUB(4)*TSQ+CSSUB(5)*QSI*RHO)
        PSSUB = MIN(SSUB,SVMAX)
        PSSUB = MAX(PSSUB,-SSMAX)

      IF(QSK .LE. QS0)  GO TO 480

!     PGAUT: Autoconversion of snow to form graupel.
        C     = 1.D-3*EXPT1
        AUT   = MAX(C*(QSK-QS0),0.0_r8)
        PGAUT = AUT

  480 IF(QGK .EQ. 0.0_r8)  GO TO 520
      IF(QIK .NE. 0.0_r8 .OR. QLK .NE. 0.0_r8 .OR. QVK .GE. QSI)  GO TO 500

!     PGSUB: NEGATIVE VALUES ONLY (Sublimation of graupel.)
        GSUB  = CGSUB(1)*TSQ*QSI*(SI-1.0_r8)*(CGSUB(2)*SQRT(QGR)          &
              + CGSUB(3)*QGR**0.6875/RHO**0.25)/(CGSUB(4)*TSQ+CGSUB(5)*QSI*RHO)
        PGSUB = MAX(GSUB,-SGMAX)

!     PGDRY (OR PGWET): Dry and wet growth of graupel.
  500   PGDRY = PGACW + PGACI + PGACR + PGACS

      IF(.NOT.LIQ .AND. .NOT.RAIN)  GO TO 510

        x    = min(PGACI*10.0_r8,SIMAX)
        y    = min(10.0_r8*PGACS,SSMAX)
        GWET = (CGWET(2)*DQS0-CGWET(1)*TC/RHO)*HLTF/(HLTF+CH2O*TC)   &
              *(CGWET(3)*SQRT(QGR)+CGWET(4)*QGR**0.6875/RHO**0.25)   &
              +(x+y)*(1.0_r8-CICE*TC/(HLTF+CH2O*TC))
        PGWET = MAX(GWET,0.0_r8)

      IF(PGDRY .LE. PGWET)  GO TO 510

!  PGWET INVOLVES REDEFINITION OF PGACI, PGACS (INCREASE OF
!     COLLECTION EFFICIENCY BY FACTOR OF 10), AND RECALCULATION
!     OF PGACR (WHICH CAN NOW BECOME NEGATIVE DUE TO SHEDDING
!     OF CLOUD LIQUID WATER WHICH IS CONVERTED INTO RAINWATER)

        PGACI  = min(PGACI * 10.0_r8,SIMAX)
        PGACS  = min(10.0_r8*PGACS,SSMAX)
        PGACR  = min(PGWET - PGACW - PGACI - PGACS , SRMAX)
        PGWORD = PGWET
      GO TO 520

  510 PGWORD = PGDRY

!     --------------------------------------
!     ----     ADD SOURCES AND SINKS    ----
!     ----              T<0             ----
!     --------------------------------------

  520 CONTINUE
        TSSV = -(PSSUB+PGSUB)+PREVP
        TSSL = -(PSACW+PSFW+PRAUT+PRACW+PGACW+PIDW)
        TSSI = -(PSAUT+PSACI+PSFI+PRACI+PGACI)+PIDW
        TSSS = PSAUT+PSACI+PSACW+PSFW+PSFI+PSSUB-(PGACS+PGAUT)

      IF(QRK .LT. 1.0D-4 .AND. QSK .LT. 1.D-4)  GO TO 530
        TSSS = TSSS-PRACS
        TSSG = PSACR+PRACS
      GO TO 540
  530   TSSS = TSSS+PSACR
  540 IF(QRK .GE. 1.0D-4)  GO TO 550
        TSSS = TSSS+PRACI
      GO TO 560
  550   TSSG = PRACI+TSSG
  560   TSSG = TSSG+PGSUB+PGAUT+PGFR+PGWORD
        TSSR = PRAUT+PRACW-(PSACR+PGACR+PREVP+PGFR)
        TTMP = (HLTF*(PSACW+PGACW+PSACR+PGACR+PGFR+PSFW+PIDW)             &
               -HLTC*PREVP+HLTS*(PSSUB+PGSUB))/CP

      TSS(1) = TSSV; TSS(2) = TSSL; TSS(3) = TSSI
      TSS(4) = TSSS; TSS(5) = TSSG; TSS(6) = TSSR; TSS(7) = TTMP
            
      TEM = 0.0_r8
do K=1,6
  TEM = TEM + TSS(K)
enddo

 1000 CONTINUE

      PSS(1)  = PRAUT; PSS(2)  = PRACW;  PSS(3)  = PRACS; PSS(4)  = PSACW; PSS(5)  = PWACS
      PSS(6)  = PGACW; PSS(7)  = PGMLT;  PSS(8)  = PSMLT; PSS(9)  = PREVP; PSS(10) = PIACR
      PSS(11) = PSACR; PSS(12) = PGACR;  PSS(13) = PGFR;  PSS(14) = PGACS; PSS(15) = PGAUT
      PSS(16) = PGACI; PSS(17) = PGWORD; PSS(18) = PRACI; PSS(19) = PSAUT; PSS(20) = PSACI
      PSS(21) = PSSUB; PSS(22) = PGSUB;  PSS(23) = PSFW;  PSS(24) = PSFI;  PSS(25) = PIDW
      PSS(26) = PGWET

!   LATENT HEATING RATE FROM ALL PROCESSES EXCEPT CONDENSATIONAL
!   AND SUBLIMATIONAL SATURATION ADJUSTMENT.  THE CONTRIBUTIONS
!   FROM THOSE PROCESSES ARE FOUND IN Q_CHK_3D.
      TLHR = - PSS(9)*(2.25D6/1004.0_r8) &
             + (PSS(25)+PSS(4)+PSS(6)+PSS(11)+PSS(12)+PSS(13)+PSS(23) &
             - (PSS(9)+PSS(3)+PSS(7)+PSS(8)))*(3.34D5/1004.0_r8) &
             - (PSS(21)+PSS(22))*(2.584D6/1004.0_r8)

RETURN
END SUBROUTINE imicro

!=============================================================================
   SUBROUTINE CAL_ACR3 (ACR3,V1,V2,QRHO1,QRHO2,C,CAC,rho)
!=============================================================================
   USE shr_kind_mod, only: r8 => shr_kind_r8

   IMPLICIT NONE

   REAL (KIND=r8), INTENT(IN) :: &
      CAC(3),    &
      V1,        &
      V2,        &
      QRHO1,     &
      QRHO2,     &
      C,         &
      RHO

   REAL (KIND=r8), INTENT(OUT) :: ACR3

!  local variables
   INTEGER :: K
   REAL (KIND=r8) :: A

!  SOURCE TERMS:  PRACS,PSACR,PGACR,PGACS

      A=0.0_r8
    do K=1,3
      A = A + CAC(K) * (QRHO1**((7-K)*0.25) * QRHO2**(K*0.25))
    enddo
      ACR3 = C * DABS(V1-V2) * A / RHO

END SUBROUTINE cal_acr3

!=============================================================================
    SUBROUTINE CAL_IDW (IDW,TC,QIK,RHO)
!=============================================================================
    USE shr_kind_mod, only: r8 => shr_kind_r8

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: &
       TC,    &
       QIK,   &
       RHO
    REAL (KIND=r8), INTENT(OUT) :: &
       IDW

! local variables

    REAL (KIND=r8) ::    &
       TC1,A1,A2,FNC,FMI,A1T(0:31),A2T(0:31)

    DATA A1T/ 0.0_r8, &
              0.7939D-12, 0.7841D-11, 0.3369D-10, 0.4336D-10, 0.5285D-10,     &
              0.3728D-10, 0.1852D-10, 0.2991D-11, 0.4248D-11, 0.7434D-11,     &
              0.1812D-10, 0.4394D-10, 0.9145D-10, 0.1725D-09, 0.3348D-09,     &
              0.1725D-09, 0.9175D-10, 0.4412D-10, 0.2252D-10, 0.9115D-11,     &
              0.4876D-11, 0.3473D-11, 0.4758D-11, 0.6306D-11, 0.8573D-11,     &
              0.7868D-11, 0.7192D-11, 0.6513D-11, 0.5956D-11, 0.5333D-11,     &
              0.4834D-11 /

    DATA A2T/ 0.0_r8,    0.4006_r8, 0.4831_r8, 0.5320_r8, &
              0.5307_r8, 0.5319_r8, 0.5249_r8, 0.4888_r8, &
              0.3894_r8, 0.4047_r8, 0.4318_r8, 0.4771_r8, &
              0.5183_r8, 0.5463_r8, 0.5651_r8, 0.5813_r8, &
              0.5655_r8, 0.5478_r8, 0.5203_r8, 0.4906_r8, &
              0.4447_r8, 0.4126_r8, 0.3960_r8, 0.4149_r8, &
              0.4320_r8, 0.4506_r8, 0.4483_r8, 0.4460_r8, &
              0.4433_r8, 0.4413_r8, 0.4382_r8, 0.4361_r8 /

    TC1 = MAX(TC,-30.0_r8)

    A1  = (A1T(-INT(TC1))-A1T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0)+A1T(-INT(TC1)+1)
    A2  = (A2T(-INT(TC1))-A2T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0)+A2T(-INT(TC1)+1)

    fnc = 0.01_r8 * exp ( - 0.6_r8 * TC )
    fmi = rho * qik / fnc * 1000.0_r8
    IDW = EXP(-0.6_r8*TC)*A1*fmi**A2/RHO

RETURN
END SUBROUTINE cal_idw

!=============================================================================
    Subroutine BERGERON (PSFW,PSFI,TC,QL,QI,QLRHO)
!=============================================================================
    USE shr_kind_mod, only: r8 => shr_kind_r8
    Use mphysics_parameters, only: c1brg,c2brg,rmi50,rmi40

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: &
       TC,      &
       QL,      &
       QI,      &
       QLRHO
    REAL (KIND=r8), INTENT(OUT) :: &
       PSFW,    &
       PSFI

! local variables
    real (KIND=r8) :: TC1,A1,A2,A21,DT1
    real (KIND=r8) :: A1T(0:31),A2T(0:31)

    DATA A1T/ 0.0_r8,  &
              0.7939D-7, 0.7841D-6, 0.3369D-5, 0.4336D-5, 0.5285D-5, &
              0.3728D-5, 0.1852D-5, 0.2991D-6, 0.4248D-6, 0.7434D-6, &
              0.1812D-5, 0.4394D-5, 0.9145D-5, 0.1725D-4, 0.3348D-4, &
              0.1725D-4, 0.9175D-5, 0.4412D-5, 0.2252D-5, 0.9115D-6, &
              0.4876D-6, 0.3473D-6, 0.4758D-6, 0.6306D-6, 0.8573D-6, &
              0.7868D-6, 0.7192D-6, 0.6513D-6, 0.5956D-6, 0.5333D-6, &
              0.4834D-6 /
    DATA A2T/ 0.0_r8,    0.4006_r8, 0.4831_r8, 0.5320_r8, &
              0.5307_r8, 0.5319_r8, 0.5249_r8, 0.4888_r8, &
              0.3894_r8, 0.4047_r8, 0.4318_r8, 0.4771_r8, &
              0.5183_r8, 0.5463_r8, 0.5651_r8, 0.5813_r8, &
              0.5655_r8, 0.5478_r8, 0.5203_r8, 0.4906_r8, &
              0.4447_r8, 0.4126_r8, 0.3960_r8, 0.4149_r8, &
              0.4320_r8, 0.4506_r8, 0.4483_r8, 0.4460_r8, &
              0.4433_r8, 0.4413_r8, 0.4382_r8, 0.4361_r8 /

    TC1 = MAX(TC,-30.0_r8)

    if ( tc1 .gt. -1.0_r8 ) then
      a1 = a1t(1) ;   a2 = a2t(1)
    else
      A1  = (A1T(-INT(TC1))-A1T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0_r8)+A1T(-INT(TC1)+1)
      A2  = (A2T(-INT(TC1))-A2T(-INT(TC1)+1))*(TC1-INT(TC1)+1.0_r8)+A2T(-INT(TC1)+1)
    endif

    A21 = 1.0_r8 - A2
    DT1 = (RMI50**A21 - RMI40**A21) / (A1*A21)

!---- NOTE:  Units are MKS

    PSFW = C1BRG*QI/DT1*(A1*RMI50**A2+QLRHO*C2BRG)
    PSFI = QI/DT1

END Subroutine bergeron

!=============================================================================
   FUNCTION GGESW(T)  RESULT(ggesw_result)
!=============================================================================
!----- Saturation vapor pressure over water in mb
! (Goff and Gratch, 1946; and Smithsonian Tables, 1984)
!----- T  is temperature in deg K

    USE shr_kind_mod, only: r8 => shr_kind_r8

    IMPLICIT NONE

    REAL (KIND=r8), INTENT(IN) :: T

    REAL (KIND=r8) :: ggesw_result
    REAL (KIND=r8) :: TS,F5,E1,E2,F1,F2,F3,F4,F

      TS = 373.16_r8
      F5 = 3.0057149_r8
      E1 = 11.344_r8*(1.0_r8 - T/TS)
      E2 =-3.49149_r8*(TS/T - 1.0_r8)
      F1 =-7.90298_r8*(TS/T - 1.0_r8)
      F2 = 5.02808_r8*LOG10(TS/T)
      F3 =-1.3816D-7*(10.0_r8**E1 - 1.0_r8)
      F4 = 8.1328D-3*(10.0_r8**E2 - 1.0_r8)
      F  = F1 + F2 + F3 + F4 + F5
      ggesw_result = 10.0_r8**F

END Function ggesw

!=============================================================================
   FUNCTION GGESI(T) RESULT(ggesi_result)
!=============================================================================
!----- Saturation vapor pressure over water in mb
! (Goff and Gratch, 1946; and Smithsonian Tables, 1984)
!----- T  is temperature in deg K

   USE shr_kind_mod, only: r8 => shr_kind_r8

   IMPLICIT NONE

   REAL (KIND=r8), INTENT(IN) :: T

   REAL (KIND=r8) :: GGESI_RESULT
   REAL (KIND=r8) :: C1,C2,C3,C4,A,B

      C1 = -9.09718_r8
      C2 = -3.56654_r8
      C3 = 0.876793_r8
      C4 = 0.78583503_r8
      A = 273.16_r8/T
      B = C1*(A-1.0_r8)+C2*LOG10(A)+C3*(1.0_r8-1.0_r8/A)+C4
      GGESI_result = 10.0_r8**B

END Function ggesi

!=============================================================================
   FUNCTION ESATW (T) RESULT(esatw_result)
!=============================================================================
!---- TT=172.16  and 173 <= T <= 323 --------

   USE shr_kind_mod, only: r8 => shr_kind_r8
   USE mphysics_parameters, only: eswt,deswt,tt

   IMPLICIT NONE

   REAL (KIND=r8), INTENT(IN) :: t
   REAL (KIND=r8) :: ESATW_RESULT
   REAL(KIND=r8) :: TTT
   INTEGER :: I

      TTT = T-TT;  I = INT(TTT)
      IF(I.LT.1 .OR. I.GT.151) THEN
       stop 131
      ENDIF
      ESATW_result = ESWT(I)+(TTT-AINT(TTT))*DESWT(I)

END Function esatw

!=============================================================================
   FUNCTION ESATI (T) RESULT(esati_result)
!=============================================================================
!---- TT=172.16  and 173 <= T <= 323 --------

   USE shr_kind_mod, only: r8 => shr_kind_r8
   Use mphysics_parameters, only: esit,desit,tt

   IMPLICIT NONE

   REAL (KIND=r8), INTENT(IN) :: T

   REAL (KIND=r8) :: ESATI_RESULT
   REAL(KIND=r8) :: TTT
   INTEGER :: I

      TTT = T-TT;  I = INT(TTT)
      IF(I.LT.1 .OR. I.GT.151) THEN
       stop 132
      ENDIF
      ESATI_result = ESIT(I)+(TTT-AINT(TTT))*DESIT(I)

END FUNCTION esati

!=============================================================================
   FUNCTION DESWDT (T) RESULT(deswdt_result)
!=============================================================================
!---- TTD = 172.66  and 173 <= T <= 323 --------

   USE shr_kind_mod, only: r8 => shr_kind_r8
   Use mphysics_parameters, only: ttd, deswt

   IMPLICIT NONE

   REAL (KIND=r8), INTENT(IN) :: T

   REAL (KIND=r8) :: DESWDT_RESULT
   REAL(KIND=r8) :: TTT
   INTEGER :: I

      TTT = T-TTD;  I = INT(TTT)
      IF(I.LT.1 .OR. I.GT.150) THEN
       stop 133
      ENDIF
      deswdt_result = DESWT(I)+(TTT-AINT(TTT))*(DESWT(I+1)-DESWT(I))

END FUNCTION deswdt

!=============================================================================
   FUNCTION DESIDT (T)  RESULT(desidt_result)
!=============================================================================
!---- TTD=172.66  and 173 <= T <= 323 --------

   USE shr_kind_mod, only: r8 => shr_kind_r8
   Use mphysics_parameters, only: ttd,desit

   IMPLICIT NONE

   REAL (KIND=r8), INTENT(IN) :: T

   REAL (KIND=r8) :: DESIDT_RESULT
   REAL(KIND=r8) :: TTT
   INTEGER :: I

      TTT = T-TTD;  I = INT(TTT)
      IF(I.LT.1 .OR. I.GT.150) THEN
       stop 132
      ENDIF
      desidt_result = DESIT(I)+(TTT-AINT(TTT))*(DESIT(I+1)-DESIT(I))

END FUNCTION desidt
