#define CLDFRAC

!WRF:MODEL_LAYER:PHYSICS
!HM: This is version 2 of Hugh Morrison's two moment, five class scheme.
!

! THIS MODULE CONTAINS THE TWO-MOMENT MICROPHYSICS CODE DESCRIBED BY
!     MORRISON ET AL. (2009, MWR)
! recent changes with respect to V1.4

! V1.5
! 1) more pathways to allow hail to form (only affects IHAIL=1 option), from collisions of snow/cloud water
! 2) bug fix to PGAM calculation (multiplication instead of division by air density)

! V1.6
! 1) added parameter TMELT for all calculations involving melting point
! 2) replaced hard-wired gas constant for air with parameter value 'R'

! V1.7
! 1) modification to minimum mixing ratio in dry conditions, change from 10^-6 to 10^-8 kg/kg
!   to improve reflectivity at low mixing ratio amounts
! 2) bug fix to prevent possible division by zero error involving LAMI
! 3) change for liquid saturation vapor pressure, replace old formula with Flatau et al. 1992

! V2
! 1) bug fix to maximum-allowed particle fallspeeds (air density correction factor considered)
! 2) change to comments

! *** Changes incorporated from WRF: ***
! CHANGES FOR V3.2, RELATIVE TO MOST RECENT (BUG-FIX) CODE FOR V3.1

! 1) ADDED ACCELERATED MELTING OF GRAUPEL/SNOW DUE TO COLLISION WITH RAIN, FOLLOWING LIN ET AL. (1983)
! 2) INCREASED MINIMUM LAMBDA FOR RAIN, AND ADDED RAIN DROP BREAKUP FOLLOWING MODIFIED VERSION
!     OF VERLINDE AND COTTON (1993)
! 3) CHANGE MINIMUM ALLOWED MIXING RATIOS IN DRY CONDITIONS (RH < 90%), THIS IMPROVES RADAR REFLECTIIVITY
!     IN LOW REFLECTIVITY REGIONS
! 4) BUG FIX TO MAXIMUM ALLOWED PARTICLE FALLSPEEDS AS A FUNCTION OF AIR DENSITY
! 5) BUG FIX TO CALCULATION OF LIQUID WATER SATURATION VAPOR PRESSURE (CHANGE IS VERY MINOR)

! bug fix, 5/12/10
! 6) bug fix for saturation vapor pressure in low pressure, to avoid division by zero

! CHANGES FOR V3.3
! 1) MODIFY FALLSPEED BELOW THE LOWEST LEVEL OF PRECIPITATION, WHICH PREVENTS
!      POTENTIAL FOR SPURIOUS ACCUMULATION OF PRECIPITATION DURING SUB-STEPPING FOR SEDIMENTATION
! 2) BUG FIX TO LATENT HEAT RELEASE DUE TO COLLISIONS OF CLOUD ICE WITH RAIN
! 3) CLEAN UP OF COMMENTS IN THE CODE
! additional minor bug fixes and small changes, 5/30/2011 (CLUBB/SAM-CLUBB as of 5 Oct 2011)
! minor revisions by A. Ackerman April 2011:
! 1) replaced kinematic with dynamic viscosity 
! 2) replaced scaling by air density for cloud droplet sedimentation
!    with viscosity-dependent Stokes expression
! 3) use Ikawa and Saito (1991) air-density scaling for cloud ice
! 4) corrected typo in 2nd digit of ventilation constant F2R

! Additional fixes
! 5) TEMPERATURE FOR ACCELERATED MELTING DUE TO COLLIIONS OF SNOW AND GRAUPEL
!    WITH RAIN SHOULD USE CELSIUS, NOT KELVIN (BUG REPORTED BY K. VAN WEVERBERG)
! 6) NPRACS IS NO SUBTRACTED SUBTRACTED FROM SNOW NUMBER CONCENTRATION, SINCE
!    DECREASE IN SNOW NUMBER IS ALREADY ACCOUNTED FOR BY NSMLTS 
! 7) MODIFY FALLSPEED BELOW THE LOWEST LEVEL OF PRECIPITATION, WHICH PREVENTS
!      POTENTIAL FOR SPURIOUS ACCUMULATION OF PRECIPITATION DURING SUB-STEPPING FOR SEDIMENTATION
! 8) BUG FIX TO LATENT HEAT RELEASE DUE TO COLLISIONS OF CLOUD ICE WITH RAIN
! 9) BUG FIX TO IGRAUP SWITCH FOR NO GRAUPEL/HAIL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL.

MODULE module_mp_GRAUPEL
!bloss   USE     module_wrf_error
!bloss      USE module_utility, ONLY: WRFU_Clock, WRFU_Alarm  ! GT
!bloss      USE module_domain, ONLY : HISTORY_ALARM, Is_alarm_tstep  ! GT

!  USE module_state_description
#ifdef CLUBB_CRM
   use constants_clubb, only: Lv, Ls, Cp, Rv, Rd, T_freeze_K, rho_lw, grav, EP_2 => ep
#else
  ! parameters from SAM and options from wrapper routine.
   use params, only: lcond, lsub, cp, rgas, rv
#endif /*CLUBB_CRM*/

#if (defined CRM && defined MODAL_AERO)
   use drop_activation, only: drop_activation_ghan
   use abortutils, only: endrun
#endif

   IMPLICIT NONE

! Adding coefficient term for clex9_oct14 case. This will reduce NNUCCD and NNUCCC
! by some factor to allow cloud to persist at realistic time intervals.

#ifdef CLUBB_CRM
!   REAL, SAVE :: NNUCCD_REDUCE_COEF = 1.0, NNUCCC_REDUCE_COEF = 1.0
   REAL, SAVE :: NNUCCD_REDUCE_COEF = 1.0, NNUCCC_REDUCE_COEF = 1.0, NNUCHOM_REDUCE_COEF=1.0e-2
#endif

! Change by Marc Pilon on 11/16/11


   REAL, PARAMETER :: PI = 3.1415926535897932384626434
   REAL, PARAMETER :: SQRTPI = 0.9189385332046727417803297

   PUBLIC  ::  MP_GRAUPEL
   PUBLIC  ::  POLYSVP

   PRIVATE :: GAMMA, DERF1
   PRIVATE :: PI, SQRTPI
   PUBLIC :: M2005MICRO_GRAUPEL !bloss

   !bloss: added options that may be set in prm file namelist 
   !         -- initialized in micrphysics.f90
   logical, public :: &
        doicemicro, &         ! use ice species (snow/cloud ice/graupel)
        dograupel, &          ! use graupel
        dohail, &             ! make graupel species have properties of hail
        dosb_warm_rain, &     ! use Seifert & Beheng (2001) warm rain parameterization
        dopredictNc, &        ! prediction of cloud droplet number
        dosubgridw, &         ! input estimate of subgrid w to microphysics
        doarcticicenucl, &    ! use arctic parameter values for ice nucleation
        docloudedgeactivation,& ! activate cloud droplets throughout the cloud
        dofix_pgam            ! option to fix value of pgam (exponent in cloud water gamma distn)

   integer, public :: &
        aerosol_mode          ! determines aerosol mode used
                              ! 0 = no aerosol mode
                              ! 1 = power-law
                              ! 2 = lognormal
#if (defined CRM && defined MODAL_AERO)
   logical, public :: domodal_aero    ! use modal aerosol from the CAM
#endif

   real, public :: &
        Nc0, &                ! specified cloud droplet number conc (#/cm3)
        ccnconst, ccnexpnt, & ! dospecifyaerosol=.false. params (powerlaw CCN)
        aer_rm1, aer_rm2, &   ! two modes of aerosol for dospecifyaer...=.true.
        aer_n1, aer_n2, &     ! rm=geom mean radius (um), n=aer conc. (#/cm3)
        aer_sig1, aer_sig2, & ! sig=geom standard deviation of aer size distn.
        pgam_fixed            ! fixed value of pgam used if dofix_pgam=.true.

! SWITCHES FOR MICROPHYSICS SCHEME
! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA
! There's no IACT = 3 in SAM / SAM-CLUBB as per WRF
#if (defined CRM && defined MODAL_AERO)
! IACT = 3, USE MULTIMODE AEROSOL SIZE DIST to DERIVER CCN SPECTRA
#endif

     INTEGER, PRIVATE ::  IACT

! INUM = 0, PREDICT DROPLET CONCENTRATION
! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION   

     INTEGER, PRIVATE ::  INUM

! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (CM-3)
     REAL, PRIVATE ::      NDCNST

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

     INTEGER, PRIVATE ::  ILIQ

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS

     INTEGER, PRIVATE ::  INUC

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING 
!             NON-EQULIBRIUM SUPERSATURATION, 
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION, BASED ON THE 
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY 
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

     INTEGER, PRIVATE ::  IBASE

! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

     INTEGER, PRIVATE ::  ISUB      

! SWITCH FOR GRAUPEL/NO GRAUPEL
! IGRAUP = 0, INCLUDE GRAUPEL
! IGRAUP = 1, NO GRAUPEL

     INTEGER, PRIVATE ::  IGRAUP

! HM ADDED NEW OPTION FOR HAIL V1.3
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING GICE IS HAIL

     INTEGER, PRIVATE ::  IHAIL

! HM ADDED 8/1/08, v1.4
! SWITCH FOR WARM RAIN SCHEME
! IRAIN = 0, WARM RAIN (AUTO, ACC, SELF-COLL) FROM KHAIROUTIDNOV AND KOGAN (2000)
! IRAIN = 1, WARM RAIN (AUTO, ACC, SELF-COLL) FROM SEIFERT AND BEHENG (2001)

     INTEGER, PRIVATE ::  IRAIN      

! PB ADDED 4/13/09
! SWITCH TO TURN ON/OFF CLOUD LIQUID WATER SATURATION ADJUSTMENT
! WHEN USING TOTAL WATER FORMULATION IN SAM, THE SATURATION 
! ADJUSTMENT IS PERFORMED BEFORE CALLING M2005MICRO_GRAUPEL.
! THIS OPTION ALLOWS US TO AVOID PERFORMING IT IN M2005MICRO_GRAUPEL
! UNDER THE THEORY THAT THE OTHER MICROPHYSICAL PROCESSES WILL NOT
! DRIVE IT FAR FROM SATURATION.
! ISATADJ = 0, SATURATION ADJUSTMENT PEROFORMED IN M2005MICRO_GRAUPEL
! ISATADJ = 1, SATURATION ADJUSTMENT _NOT_ PEROFORMED IN M2005MICRO_GRAUPEL

     INTEGER, PRIVATE :: ISATADJ

! CLOUD MICROPHYSICS CONSTANTS

     REAL, PRIVATE ::      AI,AC,AS,AR,AG ! 'A' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
     REAL, PRIVATE ::      BI,BC,BS,BR,BG ! 'B' PARAMETER IN FALLSPEED-DIAM RELATIONSHIP
     REAL, PRIVATE ::      R           ! GAS CONSTANT FOR AIR
!bloss     REAL, PRIVATE ::      RV          ! GAS CONSTANT FOR WATER VAPOR
!bloss     REAL, PRIVATE ::      CP          ! SPECIFIC HEAT AT CONSTANT PRESSURE FOR DRY AIR
     REAL, PRIVATE ::      RHOSU       ! STANDARD AIR DENSITY AT 850 MB
     REAL, PRIVATE ::      RHOW        ! DENSITY OF LIQUID WATER
     REAL, PRIVATE ::      RHOI        ! BULK DENSITY OF CLOUD ICE
     REAL, PRIVATE ::      RHOSN       ! BULK DENSITY OF SNOW
     REAL, PRIVATE ::      RHOG        ! BULK DENSITY OF GRAUPEL
     REAL, PRIVATE ::      AIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
     REAL, PRIVATE ::      BIMM        ! PARAMETER IN BIGG IMMERSION FREEZING
     REAL, PRIVATE ::      ECR         ! COLLECTION EFFICIENCY BETWEEN DROPLETS/RAIN AND SNOW/RAIN
     REAL, PRIVATE ::      DCS         ! THRESHOLD SIZE FOR CLOUD ICE AUTOCONVERSION
     REAL, PRIVATE ::      MI0         ! INITIAL SIZE OF NUCLEATED CRYSTAL
     REAL, PRIVATE ::      MG0         ! MASS OF EMBRYO GRAUPEL
     REAL, PRIVATE ::      F1S         ! VENTILATION PARAMETER FOR SNOW
     REAL, PRIVATE ::      F2S         ! VENTILATION PARAMETER FOR SNOW
     REAL, PRIVATE ::      F1R         ! VENTILATION PARAMETER FOR RAIN
     REAL, PRIVATE ::      F2R         ! VENTILATION PARAMETER FOR RAIN
     REAL, PRIVATE ::      G           ! GRAVITATIONAL ACCELERATION
     REAL, PRIVATE ::      QSMALL      ! SMALLEST ALLOWED HYDROMETEOR MIXING RATIO
     REAL, PRIVATE ::      CI,DI,CS,DS,CG,DG ! SIZE DISTRIBUTION PARAMETERS FOR CLOUD ICE, SNOW, GRAUPEL
     REAL, PRIVATE ::      EII         ! COLLECTION EFFICIENCY, ICE-ICE COLLISIONS
     REAL, PRIVATE ::      ECI         ! COLLECTION EFFICIENCY, ICE-DROPLET COLLISIONS
     REAL, PRIVATE ::      RIN     ! RADIUS OF CONTACT NUCLEI (M)
! V1.6
     REAL, PRIVATE ::      TMELT     ! melting temp (K)
! hm, add for V2.1
     REAL, PRIVATE ::      CPW     ! SPECIFIC HEAT OF LIQUID WATER

! CCN SPECTRA FOR IACT = 1

     REAL, PRIVATE ::      C1     ! 'C' IN NCCN = CS^K (CM-3)
     REAL, PRIVATE ::      K1     ! 'K' IN NCCN = CS^K

! AEROSOL PARAMETERS FOR IACT = 2

     REAL, PRIVATE ::      MW      ! MOLECULAR WEIGHT WATER (KG/MOL)
     REAL, PRIVATE ::      OSM     ! OSMOTIC COEFFICIENT
     REAL, PRIVATE ::      VI      ! NUMBER OF ION DISSOCIATED IN SOLUTION
     REAL, PRIVATE ::      EPSM    ! AEROSOL SOLUBLE FRACTION
     REAL, PRIVATE ::      RHOA    ! AEROSOL BULK DENSITY (KG/M3)
     REAL, PRIVATE ::      MAP     ! MOLECULAR WEIGHT AEROSOL (KG/MOL)
     REAL, PRIVATE ::      MA      ! MOLECULAR WEIGHT OF 'AIR' (KG/MOL)
     REAL, PRIVATE ::      RR      ! UNIVERSAL GAS CONSTANT
     REAL, PRIVATE ::      BACT    ! ACTIVATION PARAMETER
     REAL, PRIVATE ::      RM1     ! GEOMETRIC MEAN RADIUS, MODE 1 (M)
     REAL, PRIVATE ::      RM2     ! GEOMETRIC MEAN RADIUS, MODE 2 (M)
     REAL, PRIVATE ::      NANEW1  ! TOTAL AEROSOL CONCENTRATION, MODE 1 (M^-3)
     REAL, PRIVATE ::      NANEW2  ! TOTAL AEROSOL CONCENTRATION, MODE 2 (M^-3)
     REAL, PRIVATE ::      SIG1    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 1
     REAL, PRIVATE ::      SIG2    ! STANDARD DEVIATION OF AEROSOL S.D., MODE 2
     REAL, PRIVATE ::      F11     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     REAL, PRIVATE ::      F12     ! CORRECTION FACTOR FOR ACTIVATION, MODE 1
     REAL, PRIVATE ::      F21     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2
     REAL, PRIVATE ::      F22     ! CORRECTION FACTOR FOR ACTIVATION, MODE 2     
     REAL, PRIVATE ::      MMULT   ! MASS OF SPLINTERED ICE PARTICLE
     REAL, PRIVATE ::      LAMMAXI,LAMMINI,LAMMAXR,LAMMINR,LAMMAXS,LAMMINS,LAMMAXG,LAMMING

! CONSTANTS TO IMPROVE EFFICIENCY

     REAL, PRIVATE :: CONS1,CONS2,CONS3,CONS4,CONS5,CONS6,CONS7,CONS8,CONS9,CONS10
     REAL, PRIVATE :: CONS11,CONS12,CONS13,CONS14,CONS15,CONS16,CONS17,CONS18,CONS19,CONS20
     REAL, PRIVATE :: CONS21,CONS22,CONS23,CONS24,CONS25,CONS26,CONS27,CONS28,CONS29,CONS30
     REAL, PRIVATE :: CONS31,CONS32,CONS33,CONS34,CONS35,CONS36,CONS37,CONS38,CONS39,CONS40
     REAL, PRIVATE :: CONS41

! v1.4
     REAL, PRIVATE :: dnu(16)

!..Various radar related variables, from GT

!..Lookup table dimensions
      INTEGER, PARAMETER, PRIVATE:: nbins = 100
      INTEGER, PARAMETER, PRIVATE:: nbr = nbins
      INTEGER, PARAMETER, PRIVATE:: nbs = nbins
      INTEGER, PARAMETER, PRIVATE:: nbg = nbins
      REAL(8), DIMENSION(nbins+1):: ddx
      REAL(8), DIMENSION(nbr):: Dr, dtr
      REAL(8), DIMENSION(nbs):: Dds, dts
      REAL(8), DIMENSION(nbg):: Ddg, dtg
      REAL(8), PARAMETER, PRIVATE:: lamda_radar = 0.10         ! in meters
      REAL(8), PRIVATE:: K_w, PI5, lamda4
      COMPLEX*16, PRIVATE:: m_w_0, m_i_0
      REAL(8), DIMENSION(nbins+1), PRIVATE:: simpson
      REAL(8), DIMENSION(3), PARAMETER, PRIVATE:: basis =      &
                           (/1.d0/3.d0, 4.d0/3.d0, 1.d0/3.d0/)

      INTEGER, PARAMETER, PRIVATE:: slen = 20
      CHARACTER(len=slen), PRIVATE::                                    &
              mixingrulestring_s, matrixstring_s, inclusionstring_s,    &
              hoststring_s, hostmatrixstring_s, hostinclusionstring_s,  &
              mixingrulestring_g, matrixstring_g, inclusionstring_g,    &
              hoststring_g, hostmatrixstring_g, hostinclusionstring_g

      REAL, PARAMETER, PRIVATE:: D0r = 50.E-6
      REAL, PARAMETER, PRIVATE:: D0s = 100.E-6
      REAL, PARAMETER, PRIVATE:: D0g = 100.E-6
      CHARACTER*256:: mp_debug
#ifdef CLUBB_CRM
      REAL, PARAMETER, PUBLIC :: cloud_frac_thresh = 0.005
#endif /* CLUBB_CRM */

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE GRAUPEL_INIT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE INITIALIZES ALL PHYSICAL CONSTANTS AMND PARAMETERS 
! NEEDED BY THE MICROPHYSICS SCHEME.
! NEEDS TO BE CALLED AT FIRST TIME STEP, PRIOR TO CALL TO MAIN MICROPHYSICS INTERFACE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IMPLICIT NONE

      integer n,i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THE FOLLOWING PARAMETERS ARE USER-DEFINED SWITCHES AND NEED TO BE
! SET PRIOR TO CODE COMPILATION

! INUM = 0, PREDICT DROPLET CONCENTRATION
! INUM = 1, ASSUME CONSTANT DROPLET CONCENTRATION   

      INUM = 1 !bloss: use flag in prm file
      if(dopredictNc) then
         INUM = 0
      end if

! FOR INUM = 1, SET CONSTANT DROPLET CONCENTRATION (UNITS OF CM-3)

      NDCNST = Nc0 !bloss: use value from prm file (default=100.)

! IACT = 1, USE POWER-LAW CCN SPECTRA, NCCN = CS^K
! IACT = 2, USE LOGNORMAL AEROSOL SIZE DIST TO DERIVE CCN SPECTRA
! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)
#if (defined CRM && defined MODAL_AERO)
! IACT = 3, USE MULTIMODE AEROSOL SIZE DIST to DERIVER CCN SPECTRA
#endif

      if( aerosol_mode == 2 ) then !bloss: specify using flag from prm file
#if (defined CRM && defined MODAL_AERO)
        if(domodal_aero) then
          IACT = 3
        else 
#endif
         IACT = 2
#if (defined CRM && defined MODAL_AERO)
        endif
#endif
      else if( aerosol_mode == 1 ) then
         IACT = 1
      else
         IACT = 0
      end if

! IBASE = 1, NEGLECT DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING, ACTIVATE
!             AT CLOUD BASE OR IN REGION WITH LITTLE CLOUD WATER USING 
!             NON-EQULIBRIUM SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, 
!             IN CLOUD INTERIOR ACTIVATE USING EQUILIBRIUM SUPERSATURATION
! IBASE = 2, ASSUME DROPLET ACTIVATION AT LATERAL CLOUD EDGES DUE TO 
!             UNRESOLVED ENTRAINMENT AND MIXING DOMINATES,
!             ACTIVATE DROPLETS EVERYWHERE IN THE CLOUD USING NON-EQUILIBRIUM
!             SUPERSATURATION ASSUMING NO INITIAL CLOUD WATER, BASED ON THE 
!             LOCAL SUB-GRID AND/OR GRID-SCALE VERTICAL VELOCITY 
!             AT THE GRID POINT

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      if(docloudedgeactivation) then
         IBASE = 2
      else
         IBASE = 1
      end if

! INCLUDE SUB-GRID VERTICAL VELOCITY IN DROPLET ACTIVATION
! ISUB = 0, INCLUDE SUB-GRID W (RECOMMENDED FOR LOWER RESOLUTION)
! ISUB = 1, EXCLUDE SUB-GRID W, ONLY USE GRID-SCALE W

! NOTE: ONLY USED FOR PREDICTED DROPLET CONCENTRATION (INUM = 0)

      if(dosubgridw) then
         ISUB = 0
      else
         ISUB = 1      
      end if

! SWITCH FOR LIQUID-ONLY RUN
! ILIQ = 0, INCLUDE ICE
! ILIQ = 1, LIQUID ONLY, NO ICE

      if(doicemicro) then !bloss: specify using flag from prm file
         ILIQ = 0
      else
         ILIQ = 1
      end if

! SWITCH FOR ICE NUCLEATION
! INUC = 0, USE FORMULA FROM RASMUSSEN ET AL. 2002 (MID-LATITUDE)
!      = 1, USE MPACE OBSERVATIONS (ARCTIC ONLY)

      if(doarcticicenucl) then !bloss: specify using flag from prm file
         INUC = 1
      else
         INUC = 0
      end if

! SWITCH FOR GRAUPEL/NO GRAUPEL
! IGRAUP = 0, INCLUDE GRAUPEL
! IGRAUP = 1, NO GRAUPEL

      if(dograupel) then
         IGRAUP = 0
      else
         IGRAUP = 1
      end if

! HM ADDED 11/7/07, V1.3
! SWITCH FOR HAIL/GRAUPEL
! IHAIL = 0, DENSE PRECIPITATING ICE IS GRAUPEL
! IHAIL = 1, DENSE PRECIPITATING ICE IS HAIL

      if(dohail) then
         IHAIL = 1
      else
         IHAIL = 0
      end if
 
! HM ADDED 8/1/08, v1.4
! SWITCH FOR WARM RAIN SCHEME
! IRAIN = 0, WARM RAIN (AUTO, ACC, SELF-COLL) FROM KHAIROUTIDNOV AND KOGAN (2000)
! IRAIN = 1, WARM RAIN (AUTO, ACC, SELF-COLL) FROM SEIFERT AND BEHENG (2001)

      if(dosb_warm_rain) then
        IRAIN = 1
      else
        IRAIN = 0
      end if

! PB ADDED 4/13/09.  TURN OFF SATURATION ADJUSTMENT WITHIN M2005MICRO_GRAUPEL
! IN TOTAL WATER VERSION.  IT NOW TAKES PLACE BEFORE M2005MICRO_GRAUPEL IS CALLED.

#ifdef CLUBB_CRM
!      ISATADJ = 0 ! Enable for CLUBB
      ISATADJ = 1  ! When CLUBB is called, saturation adjustment is done in CLUBB, 
                   ! so should we set ISATADJ=1 here? test by Minghuai Wang +++mhwang
#else
      ISATADJ = 1
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SET PHYSICAL CONSTANTS

! FALLSPEED PARAMETERS (V=AD^B)
         AI = 700.
         AC = 3.E7
         AS = 11.72
         AR = 841.99667
         BI = 1.
         BC = 2.
         BS = 0.41
         BR = 0.8
! V1.3
         IF (IHAIL.EQ.0) THEN
	 AG = 19.3
	 BG = 0.37
         ELSE ! (MATSUN AND HUGGINS 1980)
         AG = 114.5 
         BG = 0.5
         END IF

#ifdef CLUBB_CRM
         ! Use CLUBB values for constants
         R = Rd
         RHOW  = rho_lw
         TMELT = T_freeze_K
         RHOSU = 85000./(R*TMELT)
#else
! CONSTANTS AND PARAMETERS
         !bloss: use values from params module
         R = rgas
!bloss         R = 287.15
!bloss         RV = 465.5
!bloss         CP = 1005.
! V1.6
         TMELT = 273.15
#endif
! V1.6
         RHOSU = 85000./(R*TMELT)
         RHOW = 997.
         RHOI = 500.
         RHOSN = 100.
! V1.3
         IF (IHAIL.EQ.0) THEN
	 RHOG = 400.
         ELSE
         RHOG = 900.
         END IF
         AIMM = 0.66
         BIMM = 100.
         ECR = 1.
         DCS = 125.E-6
         MI0 = 4./3.*PI*RHOI*(10.E-6)**3
	 MG0 = 1.6E-10
         F1S = 0.86
         F2S = 0.28
         F1R = 0.78
! V3 5/27/11
!        F2R = 0.32
! AA revision 4/1/11
         F2R = 0.308

#ifdef CLUBB_CRM
         G = grav
         ! Should this be set to SAM's ggr if CLUBB is not defined?
#else
         G = 9.806
#endif
         QSMALL = 1.E-14
         EII = 0.1
         ECI = 0.7
! HM, ADD FOR V3.2
         CPW = 4218.

! SIZE DISTRIBUTION PARAMETERS

         CI = RHOI*PI/6.
         DI = 3.
         CS = RHOSN*PI/6.
         DS = 3.
         CG = RHOG*PI/6.
         DG = 3.

! RADIUS OF CONTACT NUCLEI
         RIN = 0.1E-6

         MMULT = 4./3.*PI*RHOI*(5.E-6)**3

! SIZE LIMITS FOR LAMBDA

         LAMMAXI = 1./1.E-6
         LAMMINI = 1./(2.*DCS+100.E-6)
         LAMMAXR = 1./20.E-6
!         LAMMINR = 1./500.E-6
         LAMMINR = 1./2800.E-6
         LAMMAXS = 1./10.E-6
         LAMMINS = 1./2000.E-6
         LAMMAXG = 1./20.E-6
         LAMMING = 1./2000.E-6

! CCN SPECTRA FOR IACT = 1

! MARITIME
! MODIFIED FROM RASMUSSEN ET AL. 2002
! NCCN = C*S^K, NCCN IS IN CM-3, S IS SUPERSATURATION RATIO IN %

              K1 = ccnexpnt !bloss: specify using values from prm file
              C1 = ccnconst !bloss

!bloss              K1 = 0.4
!bloss              C1 = 120. 

! CONTINENTAL

!              K1 = 0.5
!              C1 = 1000. 

! AEROSOL ACTIVATION PARAMETERS FOR IACT = 2
! PARAMETERS CURRENTLY SET FOR AMMONIUM SULFATE

         MW = 0.018
         OSM = 1.
         VI = 3.
         EPSM = 0.7
         RHOA = 1777.
         MAP = 0.132
         MA = 0.0284
         RR = 8.3187
         BACT = VI*OSM*EPSM*MW*RHOA/(MAP*RHOW)

! AEROSOL SIZE DISTRIBUTION PARAMETERS CURRENTLY SET FOR MPACE 
! (see morrison et al. 2007, JGR)
! MODE 1

         RM1 = aer_rm1 !bloss: specify using values from prm file
         SIG1 = aer_sig1
         NANEW1 = aer_n1
!bloss         RM1 = 0.052E-6
!bloss         SIG1 = 2.04
!bloss         NANEW1 = 100.0E6
         F11 = 0.5*EXP(2.5*(LOG(SIG1))**2)
         F21 = 1.+0.25*LOG(SIG1)

! MODE 2

         RM2 = aer_rm2 !bloss: specify using values from prm file
         SIG2 = aer_sig2
         NANEW2 = aer_n2
!bloss         RM2 = 1.3E-6
!bloss         SIG2 = 2.5
!bloss         NANEW2 = 1.E6
         F12 = 0.5*EXP(2.5*(LOG(SIG2))**2)
         F22 = 1.+0.25*LOG(SIG2)

! CONSTANTS FOR EFFICIENCY

         CONS1=GAMMA(1.+DS)*CS
         CONS2=GAMMA(1.+DG)*CG
         CONS3=GAMMA(4.+BS)/6.
         CONS4=GAMMA(4.+BR)/6.
         CONS5=GAMMA(1.+BS)
         CONS6=GAMMA(1.+BR)
         CONS7=GAMMA(4.+BG)/6.
         CONS8=GAMMA(1.+BG)
         CONS9=GAMMA(5./2.+BR/2.)
         CONS10=GAMMA(5./2.+BS/2.)
         CONS11=GAMMA(5./2.+BG/2.)
         CONS12=GAMMA(1.+DI)*CI
         CONS13=GAMMA(BS+3.)*PI/4.*ECI
         CONS14=GAMMA(BG+3.)*PI/4.*ECI
         CONS15=-1108.*EII*PI**((1.-BS)/3.)*RHOSN**((-2.-BS)/3.)/(4.*720.)
         CONS16=GAMMA(BI+3.)*PI/4.*ECI
         CONS17=4.*2.*3.*RHOSU*PI*ECI*ECI*GAMMA(2.*BS+2.)/(8.*(RHOG-RHOSN))
         CONS18=RHOSN*RHOSN
         CONS19=RHOW*RHOW
         CONS20=20.*PI*PI*RHOW*BIMM
         CONS21=4./(DCS*RHOI)
         CONS22=PI*RHOI*DCS**3/6.
         CONS23=PI/4.*EII*GAMMA(BS+3.)
         CONS24=PI/4.*ECR*GAMMA(BR+3.)
         CONS25=PI*PI/24.*RHOW*ECR*GAMMA(BR+6.)
         CONS26=PI/6.*RHOW
         CONS27=GAMMA(1.+BI)
         CONS28=GAMMA(4.+BI)/6.
         CONS29=4./3.*PI*RHOW*(25.E-6)**3
         CONS30=4./3.*PI*RHOW
         CONS31=PI*PI*ECR*RHOSN
         CONS32=PI/2.*ECR
         CONS33=PI*PI*ECR*RHOG
         CONS34=5./2.+BR/2.
         CONS35=5./2.+BS/2.
         CONS36=5./2.+BG/2.
         CONS37=4.*PI*1.38E-23/(6.*PI*RIN)
         CONS38=PI*PI/3.*RHOW
         CONS39=PI*PI/36.*RHOW*BIMM
         CONS40=PI/6.*BIMM
         CONS41=PI*PI*ECR*RHOW

! v1.4
         dnu(1) = -0.557
         dnu(2) = -0.557
         dnu(3) = -0.430
         dnu(4) = -0.307
         dnu(5) = -0.186
         dnu(6) = -0.067
         dnu(7) = 0.050
         dnu(8) = 0.167
         dnu(9) = 0.282
         dnu(10) = 0.397
         dnu(11) = 0.512
         dnu(12) = 0.626
         dnu(13) = 0.739
         dnu(14) = 0.853
         dnu(15) = 0.966
         dnu(16) = 0.966

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables for radar reflecitivity calculations
!..Create bins of rain (from min diameter up to 5 mm).
      ddx(1) = D0r*1.0d0
      ddx(nbr+1) = 0.005d0
      do n = 2, nbr
         ddx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbr) &
                  *DLOG(ddx(nbr+1)/ddx(1)) +DLOG(ddx(1)))
      enddo
      do n = 1, nbr
         Dr(n) = DSQRT(ddx(n)*ddx(n+1))
         dtr(n) = ddx(n+1) - ddx(n)
      enddo

!..Create bins of snow (from min diameter up to 2 cm).
      Ddx(1) = D0s*1.0d0
      Ddx(nbs+1) = 0.02d0
      do n = 2, nbs
         Ddx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbs) &
                  *DLOG(Ddx(nbs+1)/Ddx(1)) +DLOG(Ddx(1)))
      enddo
      do n = 1, nbs
         Dds(n) = DSQRT(Ddx(n)*Ddx(n+1))
         dts(n) = Ddx(n+1) - Ddx(n)
      enddo

!..Create bins of graupel (from min diameter up to 5 cm).
      Ddx(1) = D0g*1.0d0
      Ddx(nbg+1) = 0.05d0
      do n = 2, nbg
         Ddx(n) = DEXP(DFLOAT(n-1)/DFLOAT(nbg) &
                  *DLOG(Ddx(nbg+1)/Ddx(1)) +DLOG(Ddx(1)))
      enddo   
      do n = 1, nbg
         Ddg(n) = DSQRT(Ddx(n)*Ddx(n+1))
         dtg(n) = Ddx(n+1) - Ddx(n)
      enddo

      do i = 1, 256
         mp_debug(i:i) = char(0)
      enddo

      call radar_init
#ifndef CLUBB_CRM
!      WRITE(0,*) "WARNING: This version of the Morrison microphysics ", &
!        "incorporates changes from WRF V3.3 not found in standard SAM."
!      STOP "Comment out this stop if you want to run this code anyway."
#endif /* not CLUBB_CRM */

END SUBROUTINE GRAUPEL_INIT

!interface copied from new thompson interface
!and added NC, NS, NR, and NG variables.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS SUBROUTINE IS MAIN INTERFACE WITH THE TWO-MOMENT MICROPHYSICS SCHEME
! THIS INTERFACE TAKES IN 3D VARIABLES FROM DRIVER MODEL, CONVERTS TO 1D FOR
! CALL TO THE MAIN MICROPHYSICS SUBROUTINE (SUBROUTINE M2005MICRO_GRAUPEL) 
! WHICH OPERATES ON 1D VERTICAL COLUMNS.
! 1D VARIABLES FROM THE MAIN MICROPHYSICS SUBROUTINE ARE THEN REASSIGNED BACK TO 3D FOR OUTPUT
! BACK TO DRIVER MODEL USING THIS INTERFACE

! ******IMPORTANT******
! THIS CODE ASSUMES THE DRIVER MODEL USES PROCESS-SPLITTING FOR SOLVING THE TIME-DEPENDENT EQS.
! THUS, MODEL VARIABLES ARE UPDATED WITH MICROPHYSICS TENDENCIES INSIDE OF THE MICROPHYSICS
! SCHEME. THESE UPDATED VARIABLES ARE PASSED BACK TO DRIVER MODEL. THIS IS WHY THERE
! ARE NO TENDENCIES PASSED BACK AND FORTH BETWEEN DRIVER AND THE INTERFACE SUBROUTINE

! AN EXCEPTION IS THE TURBULENT MIXING TENDENCIES FOR DROPLET AND CLOUD ICE NUMBER CONCENTRATIONS
! (NCTEND, NITEND BELOW). FOR APPLICATION IN MODELS OTHER THAN WRF, TURBULENT MIXING TENDENCIES
! CAN BE ADDED TO THE VARIABLES ELSEWHERE (IN DRIVER OR PBL ROUTINE), AND THEN DON'T
! NEED TO BE PASSED INTO THE SUBROUTINE HERE.....

! FOR QUESTIONS, CONTACT: HUGH MORRISON, E-MAIL: MORRISON@UCAR.EDU, PHONE:303-497-8916

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE MP_GRAUPEL(ITIMESTEP,                       &
                TH, QV, QC, QR, QI, QS, QG, NI, NC, NS, NR, NG, TKE, NCTEND, &
                NITEND,KZH,  &
                RHO, PII, P, DT_IN, DZ, HT, W,          &
                RAINNC, RAINNCV, SR                     &
               ,EFFCS,EFFIS                             & ! HM ADD 4/13/07 
               ,refl_10cm                           & ! GT
!bloss               ,grid_clock                          & ! GT
!bloss               ,grid_alarms                         & ! GT
               ,IDS,IDE, JDS,JDE, KDS,KDE               & ! domain dims
               ,IMS,IME, JMS,JME, KMS,KME               & ! memory dims
               ,ITS,ITE, JTS,JTE, KTS,KTE               & ! tile   dims            )
                                            )
 
! QV - water vapor mixing ratio (kg/kg)
! QC - cloud water mixing ratio (kg/kg)
! QR - rain water mixing ratio (kg/kg)
! QI - cloud ice mixing ratio (kg/kg)
! QS - snow mixing ratio (kg/kg)
! QG - graupel mixing ratio (KG/KG)
! NI - cloud ice number concentration (1/kg)
! NC - Droplet Number concentration (1/kg)
! NS - Snow Number concentration (1/kg)
! NR - Rain Number concentration (1/kg)
! NG - Graupel number concentration (1/kg)
! NOTE: RHO AND HT NOT USED BY THIS SCHEME AND DO NOT NEED TO BE PASSED INTO SCHEME!!!!
! P - AIR PRESSURE (PA)
! W - VERTICAL AIR VELOCITY (M/S)
! TH - POTENTIAL TEMPERATURE (K)
! PII - exner function - used to convert potential temp to temp
! DZ - difference in height over interface (m)
! DT_IN - model time step (sec)
! ITIMESTEP - time step counter
! RAINNC - accumulated grid-scale precipitation (mm)
! RAINNCV - one time step grid scale precipitation (mm/time step)
! SR - one time step mass ratio of snow to total precip
! TKE - turbulence kinetic energy (m^2 s-2), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! NCTEND - droplet concentration tendency from pbl (kg-1 s-1)
! NCTEND - CLOUD ICE concentration tendency from pbl (kg-1 s-1)
! KZH - heat eddy diffusion coefficient from YSU scheme (M^2 S-1), NEEDED FOR DROPLET ACTIVATION (SEE CODE BELOW)
! EFFCS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! EFFIS - CLOUD DROPLET EFFECTIVE RADIUS OUTPUT TO RADIATION CODE (micron)
! REFL_10CM - CALCULATED RADAR REFLECTIVITY AT 10 CM (DBZ)
!................................
! GRID_CLOCK, GRID_ALARMS - parameters to limit radar reflectivity calculation only when needed
! otherwise radar reflectivity calculation every time step is too slow
! only needed for coupling with WRF, see code below for details

! EFFC - DROPLET EFFECTIVE RADIUS (MICRON)
! EFFR - RAIN EFFECTIVE RADIUS (MICRON)
! EFFS - SNOW EFFECTIVE RADIUS (MICRON)
! EFFI - CLOUD ICE EFFECTIVE RADIUS (MICRON)

! ADDITIONAL OUTPUT FROM MICRO - SEDIMENTATION TENDENCIES, NEEDED FOR LIQUID-ICE STATIC ENERGY

! QGSTEN - GRAUPEL SEDIMENTATION TEND (KG/KG/S)
! QRSTEN - RAIN SEDIMENTATION TEND (KG/KG/S)
! QISTEN - CLOUD ICE SEDIMENTATION TEND (KG/KG/S)
! QNISTEN - SNOW SEDIMENTATION TEND (KG/KG/S)
! QCSTEN - CLOUD WATER SEDIMENTATION TEND (KG/KG/S)

! ADDITIONAL INPUT NEEDED BY MICRO
! ********NOTE: WVAR IS SHOULD BE USED IN DROPLET ACTIVATION
! FOR CASES WHEN UPDRAFT IS NOT RESOLVED, EITHER BECAUSE OF
! LOW MODEL RESOLUTION OR CLOUD TYPE

! WVAR - STANDARD DEVIATION OF SUB-GRID VERTICAL VELOCITY (M/S)

   IMPLICIT NONE

   INTEGER,      INTENT(IN   )    ::   ids, ide, jds, jde, kds, kde , &
                                       ims, ime, jms, jme, kms, kme , &
                                       its, ite, jts, jte, kts, kte
! Temporary changed from INOUT to IN

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT):: &
                          qv, qc, qr, qi, qs, qg, ni, nc, ns, nr, TH, NG, effcs, effis

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(IN):: &
                          pii, p, dz, rho, w, tke, nctend, nitend,kzh
   REAL, INTENT(IN):: dt_in
   INTEGER, INTENT(IN):: ITIMESTEP

   REAL, DIMENSION(ims:ime, jms:jme), INTENT(INOUT):: &
                          RAINNC, RAINNCV, SR
   REAL, DIMENSION(ims:ime, kms:kme, jms:jme), INTENT(INOUT)::       &  ! GT
                          refl_10cm

   REAL , DIMENSION( ims:ime , jms:jme ) , INTENT(IN) ::       ht

!bloss      TYPE (WRFU_Clock):: grid_clock                  ! GT
!bloss      TYPE (WRFU_Alarm), POINTER:: grid_alarms(:)     ! GT

   ! LOCAL VARIABLES

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)::                     &
                      effi, effs, effr, EFFG

   REAL, DIMENSION(ims:ime, kms:kme, jms:jme)::                     &
                      T, WVAR, EFFC

   REAL, DIMENSION(kts:kte) ::                                                                & 
                            QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D, NC_TEND1D,           &
                            NI_TEND1D, NS_TEND1D, NR_TEND1D,                                  &
                            QC1D, QI1D, QR1D, NC1D,NI1D, NS1D, NR1D, QS1D,                    &
                            T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, RHO1D, W1D, WVAR1D,         &
                            EFFC1D, EFFI1D, EFFS1D, EFFR1D,DZ1D,   &
   ! HM ADD GRAUPEL
                            QG_TEND1D, NG_TEND1D, QG1D, NG1D, EFFG1D, &

! ADD SEDIMENTATION TENDENCIES (UNITS OF KG/KG/S)
                            QGSTEN,QRSTEN, QISTEN, QNISTEN, QCSTEN, &

! HM add reflectivity      
                            dbz
                          
   REAL PRECPRT1D, SNOWRT1D

   INTEGER I,K,J
   
   REAL DT
   LOGICAL:: dBZ_tstep ! GT

! set dbz logical based on grid_clock
!+---+
! only calculate reflectivity when it is needed for output
! in this instance, logical dbz_tstep is set to .true.
! *******NOTE: FOR COUPLING WITH DRIVER MODEL OTHER THAN WRF,
! THIS BLOCK OF CODE WILL NEED TO BE MODIFIED TO CORRECTLY
! SET WHEN REFLECTIVIITY CALCULATION IS MADE

      dBZ_tstep = .false.
!bloss      if ( Is_alarm_tstep(grid_clock, grid_alarms(HISTORY_ALARM)) ) then
!bloss         dBZ_tstep = .true.
!bloss      endif

   ! Initialize tendencies (all set to 0) and transfer
   ! array to local variables
   DT = DT_IN   
   do I=ITS,ITE
   do J=JTS,JTE
   DO K=KTS,KTE
       T(I,K,J)        = TH(i,k,j)*PII(i,k,j)

! wvar is the ST. DEV. OF sub-grid vertical velocity, used for calculating droplet 
! activation rates.
! WVAR BE DERIVED EITHER FROM PREDICTED TKE (AS IN MYJ PBL SCHEME),
! OR FROM EDDY DIFFUSION COEFFICIENT KZH (AS IN YSU PBL SCHEME),
! DEPENDING ON THE PARTICULAR pbl SCHEME DRIVER MODEL IS COUPLED WITH
! NOTE: IF MODEL HAS HIGH ENOUGH RESOLUTION TO RESOLVE UPDRAFTS, WVAR IS 
! PROBABLY NOT NEEDED 

! for MYJ pbl scheme:
!       WVAR(I,K,J)     = (0.667*tke(i,k,j))**0.5
! for YSU pbl scheme:
       WVAR(I,K,J) = KZH(I,K,J)/20.
       WVAR(I,K,J) = MAX(0.1,WVAR(I,K,J))
       WVAR(I,K,J) = MIN(4.,WVAR(I,K,J))

! add tendency from pbl to droplet and cloud ice concentration
! NEEDED FOR WRF TEMPORARILY!!!!
! OTHER DRIVER MODELS MAY ADD TURBULENT DIFFUSION TENDENCY FOR
! SCALARS SOMEWHERE ELSE IN THE MODEL (I.E, NOT IN THE MICROPHYSICS)
! IN THIS CASE THESE 2 LINES BELOW MAY BE REMOVED
       nc(i,k,j) = nc(i,k,j)+nctend(i,k,j)*dt
       ni(i,k,j) = ni(i,k,j)+nitend(i,k,j)*dt
   END DO
   END DO
   END DO

   do i=its,ite      ! i loop (east-west)
   do j=jts,jte      ! j loop (north-south)
   !
   ! Transfer 3D arrays into 1D for microphysical calculations
   !

! hm , initialize 1d tendency arrays to zero

      do k=kts,kte   ! k loop (vertical)

          QC_TEND1D(k)  = 0.
          QI_TEND1D(k)  = 0.
          QNI_TEND1D(k) = 0.
          QR_TEND1D(k)  = 0.
          NC_TEND1D(k)  = 0.
          NI_TEND1D(k)  = 0.
          NS_TEND1D(k)  = 0.
          NR_TEND1D(k)  = 0.
          T_TEND1D(k)   = 0.
          QV_TEND1D(k)  = 0.

          QC1D(k)       = QC(i,k,j)
          QI1D(k)       = QI(i,k,j)
          QS1D(k)       = QS(i,k,j)
          QR1D(k)       = QR(i,k,j)

          NC1D(k)       = NC(i,k,j)
          NI1D(k)       = NI(i,k,j)

          NS1D(k)       = NS(i,k,j)
          NR1D(k)       = NR(i,k,j)
! HM ADD GRAUPEL
          QG1D(K)       = QG(I,K,j)
          NG1D(K)       = NG(I,K,j)
          QG_TEND1D(K)  = 0.
          NG_TEND1D(K)  = 0.

          T1D(k)        = T(i,k,j)
          QV1D(k)       = QV(i,k,j)
          P1D(k)        = P(i,k,j)
          RHO1D(k)      = P1D(K)/(R*T1D(K))
          DZ1D(k)       = DZ(i,k,j)
          W1D(k)        = W(i,k,j)
          WVAR1D(k)     = WVAR(i,k,j)
      end do

      !bloss: add extra argument for rho for consistency with below subroutine.
      !         done by repeating p1z.
      !         diable routine to make sure it is not used.
      STOP 'in mp_graupel wrapper routine.  Only use m2005micro_graupel()'

#ifndef CLUBB_CRM
!      call m2005micro_graupel(QC_TEND1D, QI_TEND1D, QNI_TEND1D, QR_TEND1D, NC_TEND1D,            &
!       NI_TEND1D, NS_TEND1D, NR_TEND1D,                                                  &
!       QC1D, QI1D, QS1D, QR1D, NC1D,NI1D, NS1D, NR1D,                                    &
!       T_TEND1D,QV_TEND1D, T1D, QV1D, P1D, RHO1D, DZ1D, W1D, WVAR1D,                   &
!       PRECPRT1D,SNOWRT1D,                                                               &
!       EFFC1D,EFFI1D,EFFS1D,EFFR1D,DT,                                                   &
!                                            IMS,IME, JMS,JME, KMS,KME,                   &
!                                            ITS,ITE, JTS,JTE, KTS,KTE,                   & ! HM ADD GRAUPEL
!                                    QG_TEND1D,NG_TEND1D,QG1D,NG1D,EFFG1D, &
! ADD SEDIMENTATION TENDENCIES
!                                  QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN)
#endif /*CLUBB_CRM*/
   !
   ! Transfer 1D arrays back into 3D arrays
   !
      do k=kts,kte

! hm, add tendencies to update global variables 
! HM, TENDENCIES FOR Q AND N NOW ADDED IN M2005MICRO, SO WE
! ONLY NEED TO TRANSFER 1D VARIABLES BACK TO 3D

          QC(i,k,j)        = QC1D(k)
          QI(i,k,j)        = QI1D(k)
          QS(i,k,j)        = QS1D(k)
          QR(i,k,j)        = QR1D(k)
          NC(i,k,j)        = NC1D(k)
          NI(i,k,j)        = NI1D(k)
          NS(i,k,j)        = NS1D(k)          
          NR(i,k,j)        = NR1D(k)
	  QG(I,K,j)        = QG1D(K)
          NG(I,K,j)        = NG1D(K)

          T(i,k,j)         = T1D(k)
          TH(I,K,J)        = T(i,k,j)/PII(i,k,j) ! CONVERT TEMP BACK TO POTENTIAL TEMP
          QV(i,k,j)        = QV1D(k)

          EFFC(i,k,j)      = EFFC1D(k)
          EFFI(i,k,j)      = EFFI1D(k)
          EFFS(i,k,j)      = EFFS1D(k)
          EFFR(i,k,j)      = EFFR1D(k)
	  EFFG(I,K,j)      = EFFG1D(K)

! EFFECTIVE RADIUS FOR RADIATION CODE
! HM, ADD LIMIT TO PREVENT BLOWING UP OPTICAL PROPERTIES, 8/18/07
! LIMITS ARE FROM THE CAM MODEL APPLIED BY ANDREW GETTELMAN
          EFFCS(I,K,J)     = MIN(EFFC(I,K,J),16.)
          EFFCS(I,K,J)     = MAX(EFFCS(I,K,J),4.)
          EFFIS(I,K,J)     = MIN(EFFI(I,K,J),130.)
          EFFIS(I,K,J)     = MAX(EFFIS(I,K,J),13.)

      end do

! hm modified so that m2005 precip variables correctly match wrf precip variables
      RAINNC(i,j) = RAINNC(I,J)+PRECPRT1D
      RAINNCV(i,j) = PRECPRT1D
      SR(i,j) = SNOWRT1D/(PRECPRT1D+1.E-12)

! add reflectivity calculations
! only calculate if logical parameter dbz_tstep = .true.

         if (dBZ_tstep) then
          call calc_refl10cm (qv1d, qr1d, qs1d, qg1d, t1d, p1d, dBZ,    &
                      kts, kte, i, j, nr1d, ns1d, ng1d)
          do k = kts, kte
             refl_10cm(i,k,j) = dBZ(k)
          enddo
         endif

   end do
   end do   

END SUBROUTINE MP_GRAUPEL

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#ifdef CLUBB_CRM
      SUBROUTINE M2005MICRO_GRAUPEL(QC3DTEN,QI3DTEN,QNI3DTEN,QR3DTEN,NC3DTEN,    &
       NI3DTEN,NS3DTEN,NR3DTEN,QC3D,QI3D,QNI3D,QR3D,NC3D,NI3D,NS3D,NR3D,         &
       T3DTEN,QV3DTEN,T3D,QV3D,PRES,RHO,DZQ,W3D,WVAR, &
! hm 7/26/11, new output
       acc1d,aut1d,evpc1d,evpr1d,mlt1d,sub1d,dep1d,con1d, &
       PRECRT,SNOWRT,            &
       EFFC,EFFI,EFFS,EFFR,DT,                                                   &
                                            IMS,IME, JMS,JME, KMS,KME,           &
                                            ITS,ITE, JTS,JTE, KTS,KTE,           & ! ADD GRAUPEL
                        QG3DTEN,NG3DTEN,QG3D,NG3D,EFFG,QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN, &
                        CF3D & ! Cloud fraction from clubb
#ifdef CLDFRAC           ! enable fractional cloudiness in Morrison microphysics 
                        ,CFL3D, CFI3D, CLDMAX &  ! Cloud fraction for liquid and ice condensate
#endif
#ifdef ECPP
                        ,C2PREC,QSINK,CSED,ISED,SSED,GSED,RSED, RH3D   & ! mhwang added, for ECPP
#endif /*ECPP*/
                       )
#else
      SUBROUTINE M2005MICRO_GRAUPEL(QC3DTEN,QI3DTEN,QNI3DTEN,QR3DTEN,NC3DTEN,    &
       NI3DTEN,NS3DTEN,NR3DTEN,QC3D,QI3D,QNI3D,QR3D,NC3D,NI3D,NS3D,NR3D,         &
       T3DTEN,QV3DTEN,T3D,QV3D,PRES,RHO,DZQ,W3D,WVAR, &
! hm 7/26/11, new output
       acc1d,aut1d,evpc1d,evpr1d,mlt1d,sub1d,dep1d,con1d, &
       PRECRT,SNOWRT,            &
       EFFC,EFFI,EFFS,EFFR,DT,                                                   &
                                            IMS,IME, JMS,JME, KMS,KME,           &
                                            ITS,ITE, JTS,JTE, KTS,KTE,           & ! ADD GRAUPEL
                        QG3DTEN,NG3DTEN,QG3D,NG3D,EFFG,QGSTEN,QRSTEN,QISTEN,QNISTEN,QCSTEN   &
#ifdef ECPP
                        ,C2PREC,QSINK,CSED,ISED,SSED,GSED,RSED, RH3D   & ! mhwang added, for ECPP
#endif /*ECPP*/
                                )
#endif /*CLUBB_CRM*/
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THIS PROGRAM IS THE MAIN TWO-MOMENT MICROPHYSICS SUBROUTINE DESCRIBED BY
! MORRISON ET AL. 2005 JAS; MORRISON AND PINTO 2005 JAS.
! ADDITIONAL CHANGE IS ADDITION OF GRAUPEL MICROPHYSICS.
! SCHEME IS DESCRIBED IN DETAIL BY MORRISON ET AL. (MONTHLY WEATHER REVIEW, IN PREP.)

! THIS SCHEME IS A BULK DOUBLE-MOMENT SCHEME THAT PREDICTS MIXING
! RATIOS AND NUMBER CONCENTRATIONS OF FIVE HYDROMETEOR SPECIES:
! CLOUD DROPLETS, CLOUD (SMALL) ICE, RAIN, SNOW, AND GRAUPEL.

! CODE STRUCTURE: MAIN SUBROUTINE IS 'M2005MICRO_GRAUPEL'. ALSO INCLUDED IN THIS FILE IS
! 'FUNCTION POLYSVP', 'FUNCTION DERF1', AND
! 'FUNCTION GAMMA'.

! NOTE: THIS SUBROUTINE USES 1D ARRAY IN VERTICAL (COLUMN), EVEN THOUGH VARIABLES ARE CALLED '3D'......

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! DECLARATIONS

      IMPLICIT NONE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! THESE VARIABLES BELOW MUST BE LINKED WITH THE MAIN MODEL.
! DEFINE ARRAY SIZES

! INPUT NUMBER OF GRID CELLS

! INPUT/OUTPUT PARAMETERS                                 ! DESCRIPTION (UNITS)
      INTEGER, INTENT( IN)  :: IMS,IME, JMS,JME, KMS,KME,          &
                               ITS,ITE, JTS,JTE, KTS,KTE

      REAL, DIMENSION(KMS:KME) ::  QC3DTEN            ! CLOUD WATER MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QI3DTEN            ! CLOUD ICE MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QNI3DTEN           ! SNOW MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QR3DTEN            ! RAIN MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NC3DTEN            ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NI3DTEN            ! CLOUD ICE NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NS3DTEN            ! SNOW NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NR3DTEN            ! RAIN NUMBER CONCENTRATION (1/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QC3D               ! CLOUD WATER MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QI3D               ! CLOUD ICE MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QNI3D              ! SNOW MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QR3D               ! RAIN MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  NC3D               ! CLOUD DROPLET NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  NI3D               ! CLOUD ICE NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  NS3D               ! SNOW NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  NR3D               ! RAIN NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  T3DTEN             ! TEMPERATURE TENDENCY (K/S)
      REAL, DIMENSION(KMS:KME) ::  QV3DTEN            ! WATER VAPOR MIXING RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  T3D                ! TEMPERATURE (K)
      REAL, DIMENSION(KMS:KME) ::  QV3D               ! WATER VAPOR MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  PRES               ! ATMOSPHERIC PRESSURE (PA)
!bloss: make rho an input argument
      REAL, DIMENSION(KMS:KME), INTENT(IN) ::  RHO   ! AIR DENSITY 
      REAL, DIMENSION(KMS:KME) ::  DZQ                ! DIFFERENCE IN HEIGHT ACROSS LEVEL (m)
      REAL, DIMENSION(KMS:KME) ::  W3D                ! GRID-SCALE VERTICAL VELOCITY (M/S)
      REAL, DIMENSION(KMS:KME) ::  WVAR               ! SUB-GRID VERTICAL VELOCITY (M/S)

! hm 7/26/11, new output
      REAL, DIMENSION(KMS:KME) ::  aut1d               ! 
      REAL, DIMENSION(KMS:KME) ::  acc1d               ! 
      REAL, DIMENSION(KMS:KME) ::  evpc1d               ! 
      REAL, DIMENSION(KMS:KME) ::  evpr1d               ! 
      REAL, DIMENSION(KMS:KME) ::  mlt1d               ! 
      REAL, DIMENSION(KMS:KME) ::  sub1d               ! 
      REAL, DIMENSION(KMS:KME) ::  dep1d               ! 
      REAL, DIMENSION(KMS:KME) ::  con1d               ! 

! HM ADDED GRAUPEL VARIABLES
      REAL, DIMENSION(KMS:KME) ::  QG3DTEN            ! GRAUPEL MIX RATIO TENDENCY (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NG3DTEN            ! GRAUPEL NUMB CONC TENDENCY (1/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QG3D            ! GRAUPEL MIX RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  NG3D            ! GRAUPEL NUMBER CONC (1/KG)

! HM, ADD 1/16/07, SEDIMENTATION TENDENCIES FOR MIXING RATIO

      REAL, DIMENSION(KMS:KME) ::  QGSTEN            ! GRAUPEL SED TEND (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QRSTEN            ! RAIN SED TEND (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QISTEN            ! CLOUD ICE SED TEND (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QNISTEN           ! SNOW SED TEND (KG/KG/S)
      REAL, DIMENSION(KMS:KME) ::  QCSTEN            ! CLOUD WAT SED TEND (KG/KG/S)      

      REAL, DIMENSION(KMS:KME) ::  NGSTEN            ! GRAUPEL SED TEND (#KG/S)
      REAL, DIMENSION(KMS:KME) ::  NRSTEN            ! RAIN SED TEND (#/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NISTEN            ! CLOUD ICE SED TEND (#/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NSSTEN           ! SNOW SED TEND (#/KG/S)
      REAL, DIMENSION(KMS:KME) ::  NCSTEN            ! CLOUD WAT SED TEND (#/KG/S)      

#ifdef CLUBB_CRM
! ADDED BY UWM JAN 7 2008
      REAL, INTENT(IN), DIMENSION(KMS:KME) ::  CF3D  ! SUBGRID SCALE CLOUD FRACTION (%)
#ifdef CLDFRAC 
      REAL, INTENT(IN), DIMENSION(KMS:KME) :: CFL3D  ! Cloud fraction for liquid condensate
      REAL, INTENT(IN), DIMENSION(KMS:KME) :: CFI3D  ! Cloud fraction for liquid condensate
#endif 
#endif
! OUTPUT VARIABLES

        REAL PRECRT                ! TOTAL PRECIP PER TIME STEP (mm)
        REAL SNOWRT                ! SNOW PER TIME STEP (mm)

        REAL, DIMENSION(KMS:KME) ::   EFFC            ! DROPLET EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KMS:KME) ::   EFFI            ! CLOUD ICE EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KMS:KME) ::   EFFS            ! SNOW EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KMS:KME) ::   EFFR            ! RAIN EFFECTIVE RADIUS (MICRON)
        REAL, DIMENSION(KMS:KME) ::   EFFG            ! GRAUPEL EFFECTIVE RADIUS (MICRON)

! MODEL INPUT PARAMETERS (FORMERLY IN COMMON BLOCKS)

        REAL DT         ! MODEL TIME STEP (SEC)

#ifdef ECPP
        REAL, DIMENSION(KMS:KME) ::   C2PREC          ! CLOUD WATER SINK rate FROM PRECIPITATION (kg/kg/s)
        REAL, DIMENSION(KMS:KME) ::   QSINK           ! CLOUD WATER SINK rate FROM PRECIPITATION (/s)
        REAL, DIMENSION(KMS:KME) ::   CSED            ! sedimentation flux of cloud water (kg/m2/s)
        REAL, DIMENSION(KMS:KME) ::   ISED            ! sedimentation flux of cloud ice (kg/m2/s)
        REAL, DIMENSION(KMS:KME) ::   SSED            ! sedimentation flux of snow (kg/m2/s)
        REAL, DIMENSION(KMS:KME) ::   GSED            ! sedimentation flux of graupel (kg/m2/s)
        REAL, DIMENSION(KMS:KME) ::   RSED            ! sedimentation flux of rain (kg/m2/s)
        REAL, DIMENSION(KMS:KME) ::   RH3D            ! relative humidity w.r.t water. 
#endif /*ECPP*/

!.....................................................................................................
! LOCAL VARIABLES: ALL PARAMETERS BELOW ARE LOCAL TO SCHEME AND DON'T NEED TO COMMUNICATE WITH THE
! REST OF THE MODEL.

! SIZE PARAMETER VARIABLES

     REAL, DIMENSION(KMS:KME) :: LAMC          ! SLOPE PARAMETER FOR DROPLETS (M-1)
     REAL, DIMENSION(KMS:KME) :: LAMI          ! SLOPE PARAMETER FOR CLOUD ICE (M-1)
     REAL, DIMENSION(KMS:KME) :: LAMS          ! SLOPE PARAMETER FOR SNOW (M-1)
     REAL, DIMENSION(KMS:KME) :: LAMR          ! SLOPE PARAMETER FOR RAIN (M-1)
     REAL, DIMENSION(KMS:KME) :: LAMG          ! SLOPE PARAMETER FOR GRAUPEL (M-1)
     REAL, DIMENSION(KMS:KME) :: CDIST1        ! PSD PARAMETER FOR DROPLETS
     REAL, DIMENSION(KMS:KME) :: N0I           ! INTERCEPT PARAMETER FOR CLOUD ICE (KG-1 M-1)
     REAL, DIMENSION(KMS:KME) :: N0S           ! INTERCEPT PARAMETER FOR SNOW (KG-1 M-1)
     REAL, DIMENSION(KMS:KME) :: N0RR          ! INTERCEPT PARAMETER FOR RAIN (KG-1 M-1)
     REAL, DIMENSION(KMS:KME) :: N0G           ! INTERCEPT PARAMETER FOR GRAUPEL (KG-1 M-1)
     REAL, DIMENSION(KMS:KME) :: PGAM          ! SPECTRAL SHAPE PARAMETER FOR DROPLETS

! MICROPHYSICAL PROCESSES

     REAL, DIMENSION(KMS:KME) ::  NCACT     ! NC TENDENCY DUE TO DROPLET ACTIVATION
     REAL, DIMENSION(KMS:KME) ::  NSUBC     ! LOSS OF NC DURING EVAP
     REAL, DIMENSION(KMS:KME) ::  NSUBI     ! LOSS OF NI DURING SUB.
     REAL, DIMENSION(KMS:KME) ::  NSUBS     ! LOSS OF NS DURING SUB.
     REAL, DIMENSION(KMS:KME) ::  NSUBR     ! LOSS OF NR DURING EVAP
     REAL, DIMENSION(KMS:KME) ::  PRD       ! DEP CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  PRE       ! EVAP OF RAIN
     REAL, DIMENSION(KMS:KME) ::  PRDS      ! DEP SNOW

! Conversion rate for the squeential splitting (only for droplet activation, 
! ice nulceation, and converion of liquid to ice, snow and graupel)
     REAL, DIMENSION(KMS:KME) ::  NSUBC0     ! LOSS OF NC DURING EVAP
     REAL, DIMENSION(KMS:KME) ::  NSUBI0     ! LOSS OF NI DURING SUB.
     REAL, DIMENSION(KMS:KME) ::  NSUBS0     ! LOSS OF NS DURING SUB.
     REAL, DIMENSION(KMS:KME) ::  NSUBG0     ! LOSS OF NG DURING EVAP
     REAL, DIMENSION(KMS:KME) ::  PRD0       ! DEP CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  PRDS0      ! DEP SNOW
     REAL, DIMENSION(KMS:KME) ::  PRDG0    ! DEP OF GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  EPRD0       ! EVAP CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  EPRDS0      ! EVAP SNOW
     REAL, DIMENSION(KMS:KME) ::  EPRDG0    ! EVAP OF GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  BERFD0       ! CONVERSION FROM LIQUID TO CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  BERFDS0       ! CONVERSION FROM LIQUID TO SNOW
     REAL, DIMENSION(KMS:KME) ::  BERFDG0      ! CONVERSION FROM LIQUID TO GRAUPEL 
     REAL, DIMENSION(KMS:KME) ::  NNUCCD0    ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL, DIMENSION(KMS:KME) ::  MNUCCD0    ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)

     REAL, DIMENSION(KMS:KME) ::  NNUCCC    ! CHANGE N DUE TO CONTACT FREEZ DROPLETS
     REAL, DIMENSION(KMS:KME) ::  MNUCCC    ! CHANGE Q DUE TO CONTACT FREEZ DROPLETS
     REAL, DIMENSION(KMS:KME) ::  PRA       ! ACCRETION DROPLETS BY RAIN
     REAL, DIMENSION(KMS:KME) ::  PRC       ! AUTOCONVERSION DROPLETS
     REAL, DIMENSION(KMS:KME) ::  PCC       ! COND/EVAP DROPLETS
     REAL, DIMENSION(KMS:KME) ::  NNUCCD    ! CHANGE N FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL, DIMENSION(KMS:KME) ::  MNUCCD    ! CHANGE Q FREEZING AEROSOL (PRIM ICE NUCLEATION)
     REAL, DIMENSION(KMS:KME) ::  MNUCCR    ! CHANGE Q DUE TO CONTACT FREEZ RAIN
     REAL, DIMENSION(KMS:KME) ::  NNUCCR    ! CHANGE N DUE TO CONTACT FREEZ RAIN
     REAL, DIMENSION(KMS:KME) ::  NPRA      ! CHANGE IN N DUE TO DROPLET ACC BY RAIN
     REAL, DIMENSION(KMS:KME) ::  NRAGG     ! SELF-COLLECTION OF RAIN
     REAL, DIMENSION(KMS:KME) ::  NSAGG     ! SELF-COLLECTION OF SNOW
     REAL, DIMENSION(KMS:KME) ::  NPRC      ! CHANGE NC AUTOCONVERSION DROPLETS
     REAL, DIMENSION(KMS:KME) ::  NPRC1      ! CHANGE NR AUTOCONVERSION DROPLETS
     REAL, DIMENSION(KMS:KME) ::  PRAI      ! CHANGE Q ACCRETION CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  PRCI      ! CHANGE Q AUTOCONVERSION CLOUD ICE BY SNOW
     REAL, DIMENSION(KMS:KME) ::  PSACWS    ! CHANGE Q DROPLET ACCRETION BY SNOW
     REAL, DIMENSION(KMS:KME) ::  NPSACWS   ! CHANGE N DROPLET ACCRETION BY SNOW
     REAL, DIMENSION(KMS:KME) ::  PSACWI    ! CHANGE Q DROPLET ACCRETION BY CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  NPSACWI   ! CHANGE N DROPLET ACCRETION BY CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  NPRCI     ! CHANGE N AUTOCONVERSION CLOUD ICE BY SNOW
     REAL, DIMENSION(KMS:KME) ::  NPRAI     ! CHANGE N ACCRETION CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  NMULTS    ! ICE MULT DUE TO RIMING DROPLETS BY SNOW
     REAL, DIMENSION(KMS:KME) ::  NMULTR    ! ICE MULT DUE TO RIMING RAIN BY SNOW
     REAL, DIMENSION(KMS:KME) ::  QMULTS    ! CHANGE Q DUE TO ICE MULT DROPLETS/SNOW
     REAL, DIMENSION(KMS:KME) ::  QMULTR    ! CHANGE Q DUE TO ICE RAIN/SNOW
     REAL, DIMENSION(KMS:KME) ::  PRACS     ! CHANGE Q RAIN-SNOW COLLECTION
     REAL, DIMENSION(KMS:KME) ::  NPRACS    ! CHANGE N RAIN-SNOW COLLECTION
     REAL, DIMENSION(KMS:KME) ::  PCCN      ! CHANGE Q DROPLET ACTIVATION
     REAL, DIMENSION(KMS:KME) ::  PSMLT     ! CHANGE Q MELTING SNOW TO RAIN
     REAL, DIMENSION(KMS:KME) ::  EVPMS     ! CHNAGE Q MELTING SNOW EVAPORATING
     REAL, DIMENSION(KMS:KME) ::  NSMLTS    ! CHANGE N MELTING SNOW
     REAL, DIMENSION(KMS:KME) ::  NSMLTR    ! CHANGE N MELTING SNOW TO RAIN
! HM ADDED 12/13/06
     REAL, DIMENSION(KMS:KME) ::  PIACR     ! CHANGE QR, ICE-RAIN COLLECTION
     REAL, DIMENSION(KMS:KME) ::  NIACR     ! CHANGE N, ICE-RAIN COLLECTION
     REAL, DIMENSION(KMS:KME) ::  PRACI     ! CHANGE QI, ICE-RAIN COLLECTION
     REAL, DIMENSION(KMS:KME) ::  PIACRS     ! CHANGE QR, ICE RAIN COLLISION, ADDED TO SNOW
     REAL, DIMENSION(KMS:KME) ::  NIACRS     ! CHANGE N, ICE RAIN COLLISION, ADDED TO SNOW
     REAL, DIMENSION(KMS:KME) ::  PRACIS     ! CHANGE QI, ICE RAIN COLLISION, ADDED TO SNOW
     REAL, DIMENSION(KMS:KME) ::  EPRD      ! SUBLIMATION CLOUD ICE
     REAL, DIMENSION(KMS:KME) ::  EPRDS     ! SUBLIMATION SNOW
! HM ADDED GRAUPEL PROCESSES
     REAL, DIMENSION(KMS:KME) ::  PRACG    ! CHANGE IN Q COLLECTION RAIN BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  PSACWG    ! CHANGE IN Q COLLECTION DROPLETS BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  PGSACW    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     REAL, DIMENSION(KMS:KME) ::  PGRACS    ! CONVERSION Q TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     REAL, DIMENSION(KMS:KME) ::  PRDG    ! DEP OF GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  EPRDG    ! SUB OF GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  EVPMG    ! CHANGE Q MELTING OF GRAUPEL AND EVAPORATION
     REAL, DIMENSION(KMS:KME) ::  PGMLT    ! CHANGE Q MELTING OF GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  NPRACG    ! CHANGE N COLLECTION RAIN BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  NPSACWG    ! CHANGE N COLLECTION DROPLETS BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  NSCNG    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION DROPLETS BY SNOW
     REAL, DIMENSION(KMS:KME) ::  NGRACS    ! CHANGE N CONVERSION TO GRAUPEL DUE TO COLLECTION RAIN BY SNOW
     REAL, DIMENSION(KMS:KME) ::  NGMLTG    ! CHANGE N MELTING GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  NGMLTR    ! CHANGE N MELTING GRAUPEL TO RAIN
     REAL, DIMENSION(KMS:KME) ::  NSUBG    ! CHANGE N SUB/DEP OF GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  PSACR    ! CONVERSION DUE TO COLL OF SNOW BY RAIN
     REAL, DIMENSION(KMS:KME) ::  NMULTG    ! ICE MULT DUE TO ACC DROPLETS BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  NMULTRG    ! ICE MULT DUE TO ACC RAIN BY GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  QMULTG    ! CHANGE Q DUE TO ICE MULT DROPLETS/GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  QMULTRG    ! CHANGE Q DUE TO ICE MULT RAIN/GRAUPEL
     REAL, DIMENSION(KMS:KME) ::  MELTQI2QR  ! Instananouc melting of QI to QR

! TIME-VARYING ATMOSPHERIC PARAMETERS

     REAL, DIMENSION(KMS:KME) ::   KAP   ! THERMAL CONDUCTIVITY OF AIR
     REAL, DIMENSION(KMS:KME) ::   EVS   ! SATURATION VAPOR PRESSURE
     REAL, DIMENSION(KMS:KME) ::   EIS   ! ICE SATURATION VAPOR PRESSURE
     REAL, DIMENSION(KMS:KME) ::   QVS   ! SATURATION MIXING RATIO
     REAL, DIMENSION(KMS:KME) ::   QVI   ! ICE SATURATION MIXING RATIO
     REAL, DIMENSION(KMS:KME) ::   QVQVS ! SAUTRATION RATIO
     REAL, DIMENSION(KMS:KME) ::   QVQVSI! ICE SATURAION RATIO
     REAL, DIMENSION(KMS:KME) ::   DV    ! DIFFUSIVITY OF WATER VAPOR IN AIR
     REAL, DIMENSION(KMS:KME) ::   XXLS  ! LATENT HEAT OF SUBLIMATION
     REAL, DIMENSION(KMS:KME) ::   XXLV  ! LATENT HEAT OF VAPORIZATION
     REAL, DIMENSION(KMS:KME) ::   CPM   ! SPECIFIC HEAT AT CONST PRESSURE FOR MOIST AIR
     REAL, DIMENSION(KMS:KME) ::   MU    ! VISCOCITY OF AIR
     REAL, DIMENSION(KMS:KME) ::   SC    ! SCHMIDT NUMBER
     REAL, DIMENSION(KMS:KME) ::   XLF   ! LATENT HEAT OF FREEZING
!bloss     REAL, DIMENSION(KMS:KME) ::   RHO   ! AIR DENSITY
     REAL, DIMENSION(KMS:KME) ::   AB    ! CORRECTION TO CONDENSATION RATE DUE TO LATENT HEATING
     REAL, DIMENSION(KMS:KME) ::   ABI    ! CORRECTION TO DEPOSITION RATE DUE TO LATENT HEATING

! TIME-VARYING MICROPHYSICS PARAMETERS

     REAL, DIMENSION(KMS:KME) ::   DAP    ! DIFFUSIVITY OF AEROSOL
     REAL    NACNT                    ! NUMBER OF CONTACT IN
     REAL    FMULT                    ! TEMP.-DEP. PARAMETER FOR RIME-SPLINTERING
     REAL    COFFI                    ! ICE AUTOCONVERSION PARAMETER

! FALL SPEED WORKING VARIABLES (DEFINED IN CODE)

      REAL, DIMENSION(KMS:KME) ::    DUMI,DUMR,DUMFNI,DUMG,DUMFNG
      REAL UNI, UMI,UMR
      REAL, DIMENSION(KMS:KME) ::    FR, FI, FNI,FG,FNG
      REAL RGVM
      REAL, DIMENSION(KMS:KME) ::   FALOUTR,FALOUTI,FALOUTNI
      REAL FALTNDR,FALTNDI,FALTNDNI,RHO2
      REAL, DIMENSION(KMS:KME) ::   DUMQS,DUMFNS
      REAL UMS,UNS
      REAL, DIMENSION(KMS:KME) ::   FS,FNS, FALOUTS,FALOUTNS,FALOUTG,FALOUTNG
      REAL FALTNDS,FALTNDNS,UNR,FALTNDG,FALTNDNG
      REAL, DIMENSION(KMS:KME) ::    DUMC,DUMFNC
      REAL UNC,UMC,UNG,UMG
      REAL, DIMENSION(KMS:KME) ::   FC,FALOUTC,FALOUTNC
      REAL FALTNDC,FALTNDNC
      REAL, DIMENSION(KMS:KME) ::   FNC,DUMFNR,FALOUTNR
      REAL FALTNDNR
      REAL, DIMENSION(KMS:KME) ::   FNR

! FALL-SPEED PARAMETER 'A' WITH AIR DENSITY CORRECTION

      REAL, DIMENSION(KMS:KME) ::    AIN,ARN,ASN,ACN,AGN

! EXTERNAL FUNCTION CALL RETURN VARIABLES

!      REAL GAMMA,      ! EULER GAMMA FUNCTION
!      REAL POLYSVP,    ! SAT. PRESSURE FUNCTION
!      REAL DERF1        ! ERROR FUNCTION

! DUMMY VARIABLES

     REAL DUM,DUM1,DUM2,DUMT,DUMQV,DUMQSS,DUMQSI,DUMS

! PROGNOSTIC SUPERSATURATION

     REAL DQSDT    ! CHANGE OF SAT. MIX. RAT. WITH TEMPERATURE
     REAL DQSIDT   ! CHANGE IN ICE SAT. MIXING RAT. WITH T
     REAL EPSI     ! 1/PHASE REL. TIME (SEE M2005), ICE
     REAL EPSS     ! 1/PHASE REL. TIME (SEE M2005), SNOW
     REAL EPSR     ! 1/PHASE REL. TIME (SEE M2005), RAIN
     REAL EPSG     ! 1/PHASE REL. TIME (SEE M2005), GRAUPEL

! NEW DROPLET ACTIVATION VARIABLES
     REAL TAUC     ! PHASE REL. TIME (SEE M2005), DROPLETS
     REAL TAUR     ! PHASE REL. TIME (SEE M2005), RAIN
     REAL TAUI     ! PHASE REL. TIME (SEE M2005), CLOUD ICE
     REAL TAUS     ! PHASE REL. TIME (SEE M2005), SNOW
     REAL TAUG     ! PHASE REL. TIME (SEE M2005), GRAUPEL
     REAL DUMACT,DUM3

! COUNTING/INDEX VARIABLES

     INTEGER K,NSTEP,N ! ,I

! LTRUE IS ONLY USED TO SPEED UP THE CODE !!
! LTRUE, SWITCH = 0, NO HYDROMETEORS IN COLUMN, 
!               = 1, HYDROMETEORS IN COLUMN

      INTEGER LTRUE

! DROPLET ACTIVATION/FREEZING AEROSOL


     REAL    CT      ! DROPLET ACTIVATION PARAMETER
     REAL    TEMP1   ! DUMMY TEMPERATURE
     REAL    SAT1    ! DUMMY SATURATION
     REAL    SIGVL   ! SURFACE TENSION LIQ/VAPOR
     REAL    KEL     ! KELVIN PARAMETER
     REAL    KC2     ! TOTAL ICE NUCLEATION RATE

       REAL CRY,KRY   ! AEROSOL ACTIVATION PARAMETERS

! MORE WORKING/DUMMY VARIABLES

     REAL DUMQI,DUMNI,DC0,DS0,DG0
     REAL DUMQC,DUMQR,RATIO,SUM_DEP,FUDGEF

! EFFECTIVE VERTICAL VELOCITY  (M/S)
     REAL WEF

! WORKING PARAMETERS FOR ICE NUCLEATION

      REAL ANUC,BNUC

! WORKING PARAMETERS FOR AEROSOL ACTIVATION

        REAL AACT,GAMM,GG,PSI,ETA1,ETA2,SM1,SM2,SMAX,UU1,UU2,ALPHA

! DUMMY SIZE DISTRIBUTION PARAMETERS

        REAL DLAMS,DLAMR,DLAMI,DLAMC,DLAMG,LAMMAX,LAMMIN

        INTEGER IDROP

#if (defined CRM && defined MODAL_AERO)
        INTEGER INES
#endif

! v1.4
! new variables for seifert and beheng warm rain scheme
      REAL, DIMENSION(KMS:KME) :: nu
      integer dumii


      REAL, DIMENSION(KMS:KME) ::  QC3D_INIT               ! Tempary variable for CLOUD WATER MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QI3D_INIT               ! Tempary variable for CLOUD ICE MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QNI3D_INIT              ! Tempary variable for SNOW MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QR3D_INIT               ! Tempary variable for RAIN MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  NC3D_INIT               ! Tempary variable for CLOUD DROPLET NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  NI3D_INIT               ! Tempary variable for CLOUD ICE NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  NS3D_INIT               ! Tempary variable for SNOW NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  NR3D_INIT               ! Tempary variable for RAIN NUMBER CONCENTRATION (1/KG)
      REAL, DIMENSION(KMS:KME) ::  T3D_TEMP                ! Tempary variable for TEMPERATURE (K)
      REAL, DIMENSION(KMS:KME) ::  QV3D_INIT               ! Tempary variable for WATER VAPOR MIXING RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  QG3D_INIT               ! Tempary variable for GRAUPEL MIX RATIO (KG/KG)
      REAL, DIMENSION(KMS:KME) ::  NG3D_INIT               ! Tempary variable for GRAUPEL NUMBER CONC (1/KG)

      REAL, DIMENSION(KMS:KME) :: QC3DTEN0, QI3DTEN0, QNI3DTEN0, QG3DTEN0, QR3DTEN0, T3DTEN0, QV3DTEN0

#ifdef CLUBB_CRM
      REAL :: QV_INIT ! Temporary variable for vapor
      REAL :: QSAT_INIT ! Temporary variable for saturation
      REAL :: TMPQSMALL ! Temporary variable for QSMALL (a lower bound in kg/kg)
      REAL :: T3D_INIT ! Temporary variable for T3D (absolute temperature in [K] )

#ifdef CLDFRAC
      REAL, DIMENSION(KMS:KME) :: CLDMAX  ! cloud fraction for precipitating hydrometers by assuming 
                                          ! maximum cloud overlap
      REAL, DIMENSION(KMS:KME) ::  CF3D_TEMP               ! Tempary variables for cloud fraction
      REAL, DIMENSION(KMS:KME) ::  CFL3D_TEMP               ! Tempary variables for cloud fraction
      REAL, DIMENSION(KMS:KME) ::  CFI3D_TEMP               ! Tempary variables for cloud fraction
      REAL, DIMENSION(KMS:KME) ::  QVCLR3D                 ! Water vapor mixing ratio in the clear-sky part
      REAL, DIMENSION(KMS:KME) ::  QV_ICLD                 ! Water vapor mixing ratio in ice clouds for mixed-phase clouds 
      REAL, DIMENSION(KMS:KME) ::  QISEVAP                 ! Evaroration of cloud ice that falls into the clear-sky grid from sedimenation 
      REAL, DIMENSION(KMS:KME) ::  QCSEVAP                 ! Evaporation of cloud water that falls into the clear-sky grid from sedimenation 
      REAL, DIMENSION(KMS:KME) ::  BERFD                   ! B-F process from liquid to ice
      REAL, DIMENSION(KMS:KME) ::  BERFDS                   ! B-F process from liquid to snow 
      REAL, DIMENSION(KMS:KME) ::  BERFDG                   ! B-F process from liquid to graupel 

      REAL :: DUMCFL, DUMCFI
      REAL :: TOTAL_WATER_BEFORE, TOTAL_WATER_AFTER
      REAL :: TOTAL_ENERGY_BEFORE, TOTAL_ENERGY_AFTER
      REAL :: QTOT_AFTER(6), QTOT_BEFORE(6)
      REAL :: QTOTE_AFTER(10), QTOTE_BEFORE(10)
      REAL :: KC2_INLIQ                                    ! ice nucleation rate inside of liquid clouds
      REAL, DIMENSION(KMS:KME) :: QTOT_TEND, QTOTE_TEND, QTOTE_TEND2, NITEN_HOMC
      REAL, DIMENSION(KMS:KME) :: NI3D_T1, NI3D_T2, NI3D_T3, NI3D_T4, QC3D_T2, QC3D_T3, NC3D_T3
      REAL, DIMENSION(KMS:KME) :: NR3D_T1, NR3D_T2, NR3D_T3, NR3D_T4
      REAL, DIMENSION(KMS:KME) :: PRD1, PRD2, PRD3
#endif 

#else 
      REAL ::EP_2 ! Dry air gas constant over water vapor gas constant [-]
      EP_2 = rgas / rv
#endif
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! SET LTRUE INITIALLY TO 0

         LTRUE = 0

! V.13 initialize effective radii to default values (from P. Blossey)
         effc(kts:kte) = 25.         
         effi(kts:kte) = 25.         
         effs(kts:kte) = 25.         
         effr(kts:kte) = 25.         
         effg(kts:kte) = 25.         

! 09/19/2011 mhwang Initialize the micropysics process rate for output 
         acc1d(kms:kme) = 0.0
         aut1d(kms:kme) = 0.0
         evpc1d(kms:kme) = 0.0 
         evpr1d(kms:kme) = 0.0 
         mlt1d(kms:kme) = 0.0 
         sub1d(kms:kme) = 0.0  
         dep1d(kms:kme) = 0.0 
         con1d(kms:kme) = 0.0

         PRC(KMS:KME) = 0.
         NPRC(KMS:KME) = 0.
         NPRC1(KMS:KME) = 0.
         PRA(KMS:KME) = 0.
         NPRA(KMS:KME) = 0.
         NRAGG(KMS:KME) = 0.
         PSMLT(KMS:KME) = 0.
         NSMLTS(KMS:KME) = 0.
         NSMLTR(KMS:KME) = 0.
         EVPMS(KMS:KME) = 0.
         PCC(KMS:KME) = 0.
         PRE(KMS:KME) = 0.
         NCACT(KMS:KME) = 0.
         NSUBC(KMS:KME) = 0.
         NSUBR(KMS:KME) = 0.
         PRACG(KMS:KME) = 0.
         NPRACG(KMS:KME) = 0.
         PSMLT(KMS:KME) = 0.
         EVPMS(KMS:KME) = 0.
         PGMLT(KMS:KME) = 0.
         EVPMG(KMS:KME) = 0.
         PRACS(KMS:KME) = 0.
         NPRACS(KMS:KME) = 0.
         NGMLTG(KMS:KME) = 0.
         NGMLTR(KMS:KME) = 0.
         MELTQI2QR(KMS:KME) = 0.

         MNUCCD0=0.
         NNUCCD0=0.
         PRD0 = 0.0
         EPRD0 = 0.0
         PRDS0 = 0.0
         EPRDS0 = 0.0
         NSUBI0 =0.0
         NSUBS0 = 0.0
         BERFD0 = 0.0
         BERFDS0 = 0.0
         NSUBC0 = 0.0
         PRDG0 = 0.
         EPRDG0 = 0.
         NSUBG0 = 0.
         BERFDG0 = 0.

         NPRC1 = 0.0
         NRAGG = 0.0
         NPRACG =  0.0
         NSUBR =  0.0
         NSMLTR =  0.
         NGMLTR = 0.
         NPRACS =  0.
         NNUCCR =  0.
         NIACR = 0.
         NIACRS = 0.
         NGRACS = 0.

         QC3DTEN0 = 0.
         QI3DTEN0 = 0.
         QNI3DTEN0 = 0.
         QG3DTEN0 = 0.
         QR3DTEN0 = 0.
         QV3DTEN0 = 0.
         T3DTEN0 = 0.

! save the grid-mean hydrometer mixing ratios 
         QC3D_INIT = QC3D
         QI3D_INIT = QI3D
         QNI3D_INIT = QNI3D
         QR3D_INIT = QR3D
         NC3D_INIT = NC3D
         NI3D_INIT = NI3D
         NS3D_INIT = NS3D
         NR3D_INIT = NR3D
         QG3D_INIT = QG3D
         NG3D_INIT = NG3D
         T3D_TEMP = T3D
         QV3D_INIT = QV3D

#ifdef CLUBB_CRM 
#ifdef CLDFRAC
       QVCLR3D = QV3D
        
!    calculate precipitation fraction based on the maximum cloud overlap 
!    This follows Morrison and Gettelman scheme in CAM5
        DO K=KTE, KTS, -1
         CF3D_TEMP(K) = min(0.999, max(CF3D(K), cloud_frac_thresh))
         CFL3D_TEMP(K) = min(0.999, max(CFL3D(K), cloud_frac_thresh))
         CFI3D_TEMP(K) = min(0.999, max(CFI3D(K), cloud_frac_thresh))
        END DO
        CLDMAX(KTE)=CF3D_TEMP(KTE)
        DO K=KTE-1,KTS,-1
         ! if rain, snow, or grauple mixing ratio is smaller than threshold, set cldmax
         ! to cloud fraction at current level
!         if(QR3D(K+1).ge.QSMALL.OR.QNI3D(K+1).ge.QSMALL.OR.QG3D(K+1).ge.QSMALL) then 
         if(QR3D(K+1).ge.1.0e-8.OR.QNI3D(K+1).ge.1.0e-8.OR.QG3D(K+1).ge.1.0e-8) then
           CLDMAX(K) = max(CLDMAX(K+1), CF3D_TEMP(K))
         else 
           CLDMAX(K) = CF3D_TEMP(K)
         end if
        END DO

        DO K=KTS, KTE
!   calculate in-cloud values for hydrometer mixing ratios used in the microphsical schemes 
         QC3D(K) = QC3D(K)/CFL3D_TEMP(K)
         NC3D(K) = NC3D(K)/CFL3D_TEMP(K)
         QI3D(K) = QI3D(K)/CFI3D_TEMP(K)
         NI3D(K) = NI3D(K)/CFI3D_TEMP(K)
         QR3D(K) = QR3D(K)/CLDMAX(K)
         NR3D(K) = NR3D(K)/CLDMAX(K)
         QNI3D(K) = QNI3D(K)/CLDMAX(K)
         NS3D(K) = NS3D(K)/CLDMAX(K)
         QG3D(K) = QG3D(K)/CLDMAX(K)
         NG3D(K) = NG3D(K)/CLDMAX(K)
        END DO
#endif
#endif


! ATMOSPHERIC PARAMETERS THAT VARY IN TIME AND HEIGHT
         DO K = KTS,KTE

#ifdef ECPP
! INITIALIZE VARIABLES FOR ECPP OUTPUT TO ZERO
                C2PREC(K)=0.
                QSINK(K)=0.
                CSED(K)=0.
                ISED(K)=0.
                SSED(K)=0.
                GSED(K)=0.
                RSED(K)=0.
                RH3D(K)=0.
#endif /*ECPP*/

#ifdef CLUBB_CRM
            XXLV    = Lv
            XXLS(K) = Ls
            CPM(K)  = Cp
#else
! LATENT HEAT OF VAPORATION

            XXLV(K) = lcond !bloss 3.1484E6-2370.*T3D(K)

! LATENT HEAT OF SUBLIMATION

            XXLS(K) = lsub !bloss 3.15E6-2370.*T3D(K)+0.3337E6

            CPM(K) = cp !bloss CP*(1.+0.887*QV3D(K))

#endif
! SATURATION VAPOR PRESSURE AND MIXING RATIO

! hm, add fix for low pressure, 5/12/10
            EVS(K) = min(0.99*pres(k),POLYSVP(T3D(K),0))   ! PA
            EIS(K) = min(0.99*pres(k),POLYSVP(T3D(K),1))   ! PA

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

            QVS(K) = EP_2*EVS(K)/(PRES(K)-EVS(K))
            QVI(K) = EP_2*EIS(K)/(PRES(K)-EIS(K))

#ifdef CLUBB_CRM
! ADDITION BY UWM TO WEIGHT BY SGS CLOUD FRACTION
! We assume that Morrison microphysics only acts within cloud
!             IF ( CF3D(K) > cloud_frac_thresh ) THEN
!               T3D_INIT = T3D(K) ! SAVE TEMPERATURE
!               QV_INIT = QV3D(K) ! SAVE VAPOR

               ! We now set QV3D to be saturated w.r.t liquid at all
               ! temperatures -dschanen 15 May 2009
!              IF ( T3D(K) < 273.15 ) THEN
!                QV3D(K) = QVI(K) ! SET VAPOR TO ICE SATURATION WITHIN CLOUD
!                TMPQSAT = QVI(K) ! Save value
!              ELSE
!                 QV3D(K) = QVS(K) ! SET VAPOR TO LIQUID SATURATION WITHIN CLOUD
!                 QSAT_INIT = QVS(K) ! Save value
!              END IF 
!               QV3D(k) = QV_INIT 
!               QSAT_INIT = QV_INIT 

!               QC3D(K) = QC3D(K) / CF3D(K) ! Within cloud cloud water mix ratio

!               IF ( INUM == 0 ) THEN
!                 NC3D(K) = NC3D(K) / CF3D(K) ! Cloud drop num conc
!               END IF

!               QR3D(K) = QR3D(K) / CF3D(K) ! Rain mix ratio
!               NR3D(K) = NR3D(K) / CF3D(K) ! Rain num conc

!               IF ( ILIQ == 0 ) THEN
!                 QI3D(K) = QI3D(K) / CF3D(K) ! Ice mix ratio
!                 NI3D(K) = NI3D(K) / CF3D(K) ! Ice num conc
!                 QNI3D(K) = QNI3D(K) / CF3D(K) ! Snow mix ratio
!                 NS3D(K) = NS3D(K) / CF3D(K) ! Snow num conc
!               END IF
!               IF ( IGRAUP == 0 ) THEN
!                 QG3D(K) = QG3D(K) / CF3D(K) ! Graupel mix ratio
!                 NG3D(K) = NG3D(K) / CF3D(K) ! Graupel num conc
!               END IF
!             END IF
#endif
            QVQVS(K) = QV3D(K)/QVS(K)
            QVQVSI(K) = QV3D(K)/QVI(K)

! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER
! V1.3, change limit from 10^-7 to 10^-6
! V1.7 7/9/09 change limit from 10^-6 to 10^-8
! this improves reflectivity at low mixing ratios

#ifndef CLUBB_CRM ! 
#ifndef CLDFRAC  ! For fractional cloudiness, QVQVS and/or QVQVSI can be less than 0.9 
             IF (QVQVS(K).LT.0.9) THEN
               IF (QR3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QR3D(K)
                  T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                  QR3D(K)=0.
               END IF
               IF (QC3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QC3D(K)
                  T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                  QC3D(K)=0.
               END IF
             END IF

             IF (QVQVSI(K).LT.0.9) THEN
               IF (QI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QI3D(K)
                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  QI3D(K)=0.
               END IF
               IF (QNI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QNI3D(K)
                  T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                  QNI3D(K)=0.
               END IF
               IF (QG3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QG3D(K)
                  T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                  QG3D(K)=0.
               END IF
             END IF

#endif
#endif

! AIR DENSITY

!bloss: now an input argument            RHO(K) = PRES(K)/(R*T3D(K))

! HEAT OF FUSION

            XLF(K) = XXLS(K)-XXLV(K)

!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       IF (QC3D(K).LT.QSMALL) THEN
!+++mhwang
#ifndef  CLDFRAC
         QV3D(K)=QV3D(K)+QC3D(K)
         T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
#else
!   QV3D and T3D are grid-mean values while QC3D is in-cloud value, 
!   so cloud fraction is needed here
         QV3D(K)=QV3D(K)+QC3D(K)*CFL3D_TEMP(K)
         T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)*CFL3D_TEMP(K)
#endif
!---mhwang
         QC3D(K) = 0.
         NC3D(K) = 0.
         EFFC(K) = 0.
       END IF
       IF (QR3D(K).LT.QSMALL) THEN
!+++mhwang
#ifndef CLDFRAC
         QV3D(K)=QV3D(K)+QR3D(K)
         T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
#else
         QV3D(K)=QV3D(K)+QR3D(K)*CLDMAX(K)
         T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)*CLDMAX(K)
#endif
!---mhwang
         QR3D(K) = 0.
         NR3D(K) = 0.
         EFFR(K) = 0.
       END IF
       IF (QI3D(K).LT.QSMALL) THEN
!+++mhwang
#ifndef CLDFRAC
         QV3D(K)=QV3D(K)+QI3D(K)
         T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
#else
         QV3D(K)=QV3D(K)+QI3D(K)*CFI3D_TEMP(K) 
         T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)*CFI3D_TEMP(K)
#endif 
!+++mhwang
         QI3D(K) = 0.
         NI3D(K) = 0.
         EFFI(K) = 0.
       END IF
       IF (QNI3D(K).LT.QSMALL) THEN
!+++mhwang
#ifndef CLDFRAC
         QV3D(K)=QV3D(K)+QNI3D(K)
         T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
#else
         QV3D(K)=QV3D(K)+QNI3D(K)*CLDMAX(K)
         T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)*CLDMAX(K)
#endif
!+++mhwang
         QNI3D(K) = 0.
         NS3D(K) = 0.
         EFFS(K) = 0.
       END IF
       IF (QG3D(K).LT.QSMALL) THEN
!+++mhwang
#ifndef CLDFRAC
         QV3D(K)=QV3D(K)+QG3D(K)
         T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
#else
         QV3D(K)=QV3D(K)+QG3D(K)*CLDMAX(K)
         T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)*CLDMAX(K)
#endif
!+++mhwang
         QG3D(K) = 0.
         NG3D(K) = 0.
         EFFG(K) = 0.
       END IF

! INITIALIZE SEDIMENTATION TENDENCIES FOR MIXING RATIO

      QRSTEN(K) = 0.
      QISTEN(K) = 0.
      QNISTEN(K) = 0.
      QCSTEN(K) = 0.
      QGSTEN(K) = 0.

      NRSTEN(K) = 0.
      NISTEN(K) = 0.
      NSSTEN(K) = 0.
      NCSTEN(K) = 0.
      NGSTEN(K) = 0.

#ifdef CLDFRAC
      QISEVAP(K) = 0.
      QCSEVAP(K) = 0.
#endif

!..................................................................
! MICROPHYSICS PARAMETERS VARYING IN TIME/HEIGHT

! DYNAMIC VISCOSITY OF AIR
! fix 053011
            MU(K) = 1.496E-6*T3D(K)**1.5/(T3D(K)+120.)

! FALL SPEED WITH DENSITY CORRECTION (HEYMSFIELD AND BENSSEMER 2006)

            DUM = (RHOSU/RHO(K))**0.54

! fix 053011
!            AIN(K) = DUM*AI
! AA revision 4/1/11: Ikawa and Saito 1991 air-density correction 
!           AIN(K) = (RHOSU/RHO(K))**0.35
! HM bug fix 10/32/2011 
            AIN(K) = (RHOSU/RHO(K))**0.35*AI
            ARN(K) = DUM*AR
            ASN(K) = DUM*AS
!            ACN(K) = DUM*AC
! AA revision 4/1/11: temperature-dependent Stokes fall speed
            ACN(K) = G*RHOW/(18.*MU(K))
! HM ADD GRAUPEL 8/28/06
            AGN(K) = DUM*AG

! V1.7
! bug fix 7/10/09 
!hm 4/15/09 bug fix, initialize lami to prevent later division by zero
            LAMI(K)=0.

!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, AND IF SUBSATURATED, THEN SKIP MICROPHYSICS
! FOR THIS LEVEL

            IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
                 .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) THEN
                 IF (T3D(K).LT.TMELT.AND.QVQVSI(K).LT.0.999) GOTO 200
                 IF (T3D(K).GE.TMELT.AND.QVQVS(K).LT.0.999) GOTO 200
            END IF

! THERMAL CONDUCTIVITY FOR AIR

! fix 053011
            KAP(K) = 1.414E3*MU(K)

! DIFFUSIVITY OF WATER VAPOR

            DV(K) = 8.794E-5*T3D(K)**1.81/PRES(K)

! SCHMIT NUMBER

! fix 053011
            SC(K) = MU(K)/(RHO(K)*DV(K))

! PSYCHOMETIC CORRECTIONS

! RATE OF CHANGE SAT. MIX. RATIO WITH TEMPERATURE

            DUM = (RV*T3D(K)**2)

            DQSDT = XXLV(K)*QVS(K)/DUM
            DQSIDT =  XXLS(K)*QVI(K)/DUM

            ABI(K) = 1.+DQSIDT*XXLS(K)/CPM(K)
            AB(K) = 1.+DQSDT*XXLV(K)/CPM(K)

! 
!.....................................................................
!.....................................................................
! CASE FOR TEMPERATURE ABOVE FREEZING

            IF (T3D(K).GE.TMELT) THEN

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

         IF (INUM.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
            NC3D(K)=NDCNST*1.E6/RHO(K)
         END IF

   NC3D(K) = MAX(0.,NC3D(K))
!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         !bloss: option for fixing pgam
         if(dofix_pgam) then
            pgam(k) = pgam_fixed
         else

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714

         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

         end if
! v1.4
! interpolate
         dumii=int(pgam(k))
         nu(k)=dnu(dumii)+(dnu(dumii+1)-dnu(dumii))* &
               (pgam(k)-real(dumii))

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

      END IF

!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF
      END IF

! CALCULATE EPSR
      IF (QR3D(K).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
      ELSE
      EPSR = 0.
      END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PCC(K) = 0.0
#ifndef CLDFRAC ! for fractional cloudiness, no saturation adjustment is allowed. 
      IF(ISATADJ.EQ.0) THEN !PB 4/13/09

! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(K)
      DUMQV = QV3D(K)
! hm, add fix for low pressure, 5/12/10
      dum=min(0.99*pres(k),POLYSVP(DUMT,0))
      DUMQSS = EP_2*dum/(PRES(K)-dum)
      DUMQC = QC3D(K)
      DUMQC = MAX(DUMQC,0.)

! SATURATION ADJUSTMENT FOR LIQUID

      DUMS = DUMQV-DUMQSS
      PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
!      IF (PCC(K)*DT+DUMQC.LT.0.) THEN
!           PCC(K) = -DUMQC/DT
!      END IF
!+++mhwang
      IF (PCC(K)*DT+QC3D(K).LT.0.) THEN
           PCC(K) = -(QC3D(K))/DT
      END IF
!---mhwang

!      QV3DTEN(K) = QV3DTEN(K)-PCC(K)
!      T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
!      QC3DTEN(K) = QC3DTEN(K)+PCC(K)
       QV3D(K) = QV3D(K)-PCC(K)*DT
       QC3D(K) = QC3D(K)+PCC(K)*DT
       T3D(K) = T3D(K)+PCC(K)*XXLV(K)/CPM(K)*DT
      END IF
#endif

      QC3DTEN(K) = 0.0
      NC3DTEN(K) = 0.0
      T3DTEN(K) = 0.0
      NCACT(K) = 0.0
!.......................................................................
! ACTIVATION OF CLOUD DROPLETS

!bloss: only do activation if droplet number is predicted
!bloss      IF (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL) THEN
      IF (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL.AND.INUM.EQ.0) THEN

! EFFECTIVE VERTICAL VELOCITY (M/S)

      IF (ISUB.EQ.0) THEN
! ADD SUB-GRID VERTICAL VELOCITY
         DUM = W3D(K)+WVAR(K)

! ASSUME MINIMUM EFF. SUB-GRID VELOCITY 0.10 M/S
          DUM = MAX(DUM,0.10)

      ELSE IF (ISUB.EQ.1) THEN
          DUM=W3D(K)
      END IF

! ONLY ACTIVATE IN REGIONS OF UPWARD MOTION
      IF (DUM.GE.0.001) THEN

      IF (IBASE.EQ.1) THEN

! ACTIVATE ONLY IF THERE IS LITTLE CLOUD WATER
! OR IF AT CLOUD BASE, OR AT LOWEST MODEL LEVEL (K=1)

         IDROP=0

! V1.3 USE CURRENT VALUE OF QC FOR IDROP
         IF (QC3D(K).LE.0.05E-3/RHO(K)) THEN
            IDROP=1
         END IF
         IF (K.EQ.1) THEN
            IDROP=1
         ELSE IF (K.GE.2) THEN
            IF (QC3D(K).GT.0.05E-3/RHO(K).AND. &
             QC3D(K-1).LE.0.05E-3/RHO(K-1)) THEN
            IDROP=1
            END IF
         END IF

         IF (IDROP.EQ.1) THEN
! ACTIVATE AT CLOUD BASE OR REGIONS WITH VERY LITTLE LIQ WATER

           IF (IACT.EQ.1) THEN
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

            DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
            DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else 
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif

           ELSE IF (IACT.EQ.2) THEN
! DROPLET ACTIVATION FROM ABDUL-RAZZAK AND GHAN (2000)

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#if (defined CRM && defined MODAL_AERO)
            ELSE if (IACT.EQ.3) then
              INES = 0 
              CALL DROP_ACTIVATION_GHAN(DUM, T3D(k), RHO(k),  &
                   DUM2, INES, SMAX, K)  
              DUM2 = (DUM2-NC3D(K))/DT 
              DUM2 = MAX(0., DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#endif
           END IF  ! IACT

!.............................................................................
        ELSE IF (IDROP.EQ.0) THEN
! ACTIVATE IN CLOUD INTERIOR
! FIND EQUILIBRIUM SUPERSATURATION

           TAUC=1./(2.*PI*RHO(k)*DV(K)*NC3D(K)*(PGAM(K)+1.)/LAMC(K))
           IF (EPSR.GT.1.E-8) THEN
             TAUR=1./EPSR
           ELSE
             TAUR=1.E8
           END IF

           DUM3=(QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
           DUM3=DUM3*TAUC*TAUR/(TAUC+TAUR)

           IF (DUM3/QVS(K).GE.1.E-6) THEN
           IF (IACT.EQ.1) THEN

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

            DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
            DUMACT = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))

! USE POWER LAW CCN SPECTRA

! CONVERT FROM ABSOLUTE SUPERSATURATION TO SUPERSATURATION RATIO IN %
            DUM3=DUM3/QVS(K)*100.

            DUM2=C1*DUM3**K1
! MAKE SURE VALUE DOESN'T EXCEED THAT FOR NON-EQUILIBRIUM SS
            DUM2=MIN(DUM2,DUMACT)
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif

           ELSE IF (IACT.EQ.2) THEN

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

           DUMACT = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

! USE LOGNORMAL AEROSOL
           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

! GET SUPERSATURATION RATIO FROM ABSOLUTE SUPERSATURATION
           SMAX = DUM3/QVS(K)

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

! MAKE SURE ISN'T GREATER THAN NON-EQUIL. SS
            DUM2=MIN(DUM2,DUMACT)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#if (defined CRM && defined MODAL_AERO)
            ELSE if (IACT.EQ.3) then
              INES =1 
! GET SUPERSATURATION RATIO FROM ABSOLUTE SUPERSATURATION
              SMAX = DUM3/QVS(K)
              CALL DROP_ACTIVATION_GHAN(DUM, T3D(k), RHO(k),  &
                   DUM2, INES, SMAX, K)
              DUM2 = (DUM2-NC3D(K))/DT
              DUM2 = MAX(0., DUM2)
#ifndef CLDFRAC
              NC3DTEN(K) = NC3DTEN(K)+DUM2
#else
              NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#endif

           END IF ! IACT
           END IF ! DUM3/QVS > 1.E-6
        END IF  ! IDROP = 1

!.......................................................................
      ELSE IF (IBASE.EQ.2) THEN

           IF (IACT.EQ.1) THEN
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

            DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
            DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif

           ELSE IF (IACT.EQ.2) THEN

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#if (defined CRM && defined MODAL_AERO)
            ELSE if (IACT.EQ.3) then
              INES = 0
              CALL DROP_ACTIVATION_GHAN(DUM, T3D(k), RHO(k),  &
                   DUM2, INES, SMAX, K)
              DUM2 = (DUM2-NC3D(K))/DT
              DUM2 = MAX(0., DUM2)
#ifndef CLDFRAC
              NC3DTEN(K) = NC3DTEN(K)+DUM2
#else
              NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
           END IF  ! IACT
        END IF  ! IBASE
        END IF  ! W > 0.001
#ifndef CLDFRAC
          NC3D(K) = NC3D(K) + NC3DTEN(K)*DT
#else
          NC3D(K) = (NC3D(K)* CFL3D_TEMP(K)+ NC3DTEN(K) * DT)/CFL3D_TEMP(K)
#endif
          NCACT(K) = NC3DTEN(K)
          NC3DTEN(K) = 0.0
        END IF  ! QC3D > QSMALL


! GET SIZE DISTRIBUTION PARAMETERS

! MELT VERY SMALL SNOW AND GRAUPEL MIXING RATIOS, ADD TO RAIN
       IF (QNI3D(K).LT.1.E-6) THEN
          QR3D(K)=QR3D(K)+QNI3D(K)
          NR3D(K)=NR3D(K)+NS3D(K)
#ifndef CLDFRAC
          T3D(K)=T3D(K)-QNI3D(K)*XLF(K)/CPM(K)
#else
!   QR3D and QNI3D are in-cloud value, and assumed to have the same fraction.
!   But T3D is grid-mean values  
          T3D(K)=T3D(K)-QNI3D(K)*XLF(K)/CPM(K)*CLDMAX(K)
#endif 
          QNI3D(K) = 0.
          NS3D(K) = 0.
       END IF
       IF (QG3D(K).LT.1.E-6) THEN
          QR3D(K)=QR3D(K)+QG3D(K)
          NR3D(K)=NR3D(K)+NG3D(K)
#ifndef CLDFRAC
          T3D(K)=T3D(K)-QG3D(K)*XLF(K)/CPM(K)
#else
          T3D(K)=T3D(K)-QG3D(K)*XLF(K)/CPM(K)*CLDMAX(K)
#endif 
          QG3D(K) = 0.
          NG3D(K) = 0.
       END IF

       IF (QC3D(K).LT.QSMALL.AND.QNI3D(K).LT.1.E-8.AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.1.E-8) GOTO 300

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NS3D(K) = MAX(0.,NS3D(K))
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))

      NR3D_T1(K) = NR3D(K)*CLDMAX(K)  !+++mhwang test

!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF
      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         !bloss: option for fixing pgam
         if(dofix_pgam) then
            pgam(k) = pgam_fixed
         else

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714

         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

         end if
! v1.4
! interpolate
         dumii=int(pgam(k))
         nu(k)=dnu(dumii)+(dnu(dumii+1)-dnu(dumii))* &
               (pgam(k)-real(dumii))

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

      END IF

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
      N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(K).LT.LAMMINS) THEN
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      ELSE IF (LAMS(K).GT.LAMMAXS) THEN

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)
      END IF
      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
      N0G(K) = NG3D(K)*LAMG(K)

! ADJUST VARS

      IF (LAMG(K).LT.LAMMING) THEN
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      ELSE IF (LAMG(K).GT.LAMMAXG) THEN

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      END IF
      END IF

!.....................................................................
! ZERO OUT PROCESS RATES

            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            PSMLT(K) = 0.
            NSMLTS(K) = 0.
            NSMLTR(K) = 0.
            EVPMS(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            NSUBC(K) = 0.
            NSUBR(K) = 0.
            PRACG(K) = 0.
            NPRACG(K) = 0.
            PSMLT(K) = 0.
            EVPMS(K) = 0.
            PGMLT(K) = 0.
            EVPMG(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NGMLTG(K) = 0.
            NGMLTR(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES, T > 273.15 K

!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

         IF (QC3D(K).GE.1.E-6) THEN

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

            IF (IRAIN.EQ.0) THEN

                PRC(K)=1350.*QC3D(K)**2.47*  &
           (NC3D(K)/1.e6*RHO(K))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

        NPRC1(K) = PRC(K)/CONS29
        NPRC(K) = PRC(K)/(QC3D(k)/NC3D(K))

                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)

            ELSE IF (IRAIN.EQ.1) THEN

! v1.4
! replace with seifert and beheng

        dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
        dum1 = 600.*dum**0.68*(1.-dum**0.68)**3

        prc(k) = 9.44e9/(20.*2.6e-7)* &
        (nu(k)+2.)*(nu(k)+4.)/(nu(k)+1.)**2* &
        (rho(k)*qc3d(k)/1000.)**4/(rho(k)*nc3d(k)/1.e6)**2* &
        (1.+dum1/(1.-dum)**2)*1000./rho(k)

        nprc(k) = prc(k)*2./2.6e-7*1000.
        nprc1(k) = 0.5*nprc(k)

        END IF
         END IF

!.......................................................................
! HM ADD 12/13/06, COLLECTION OF SNOW BY RAIN ABOVE FREEZING
! FORMULA FROM IKAWA AND SAITO (1991)

         IF (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) THEN

            UMS = ASN(K)*CONS3/(LAMS(K)**BS)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNS = ASN(K)*CONS5/LAMS(K)**BS
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS

! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACS(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
                 N0RR(K)*N0S(K)/LAMS(K)**3*                    &
                  (5./(LAMS(K)**3*LAMR(K))+                    &
                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
                  0.5/(LAMS(K)*LAMR(K)**3)))

! fix 053011, npracs no longer subtracted from snow
!            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
!                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
!                (1./(LAMR(K)**3*LAMS(K))+                      &
!                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
!                 1./(LAMR(K)*LAMS(K)**3))

         END IF

! ADD COLLECTION OF GRAUPEL BY RAIN ABOVE FREEZING
! ASSUME ALL RAIN COLLECTION BY GRAUPEL ABOVE FREEZING IS SHED
! ASSUME SHED DROPS ARE 1 MM IN SIZE

         IF (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) THEN

            UMG = AGN(K)*CONS7/(LAMG(K)**BG)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNG = AGN(K)*CONS8/LAMG(K)**BG
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

! PRACG IS MIXING RATIO OF RAIN PER SEC COLLECTED BY GRAUPEL/HAIL
            PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMG(K))+                    &
                  2./(LAMR(K)**2*LAMG(K)**2)+				   &
				  0.5/(LAMR(k)*LAMG(k)**3)))

! ASSUME 1 MM DROPS ARE SHED, GET NUMBER CONC (KG-1) SHED PER SEC

            DUM = PRACG(K)/5.2E-7

! GET NUMBER CONC OF RAIN DROPS COLLECTED

            NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                (1./(LAMR(K)**3*LAMG(K))+                      &
                 1./(LAMR(K)**2*LAMG(K)**2)+                   &
                 1./(LAMR(K)*LAMG(K)**3))

            NPRACG(K)=MAX(NPRACG(K)-DUM,0.)

	    END IF

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

            IF (IRAIN.EQ.0) THEN

           DUM=(QC3D(K)*QR3D(K))
           PRA(K) = 67.*(DUM)**1.15
           NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

           ELSE IF (IRAIN.EQ.1) THEN

! v1.4
! seifert and beheng (2001) formulation

           dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
           dum1 = (dum/(dum+5.e-4))**4
           pra(k) = 5.78e3*rho(k)/1000.*qc3d(k)*qr3d(k)*dum1
           npra(k) = pra(k)*rho(k)/1000.*(nc3d(k)*rho(k)/1.e6)/ &
           (qc3d(k)*rho(k)/1000.)*1.e6/rho(k)

         END IF
         END IF
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

! v1.4, replace with seifert and beheng (2001)

         IF (QR3D(K).GE.1.E-8) THEN
! include breakup add 10/09/09
            dum1=300.e-6
            if (1./lamr(k).lt.dum1) then
            dum=1.
            else if (1./lamr(k).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(k)-dum1))
            end if
!            NRAGG(K) = -8.*NR3D(K)*QR3D(K)*RHO(K)
            NRAGG(K) = -5.78*dum*NR3D(K)*QR3D(K)*RHO(K)
!            NRAGG(K) = 0.0
         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP OF RAIN (RUTLEDGE AND HOBBS 1983)

      IF (QR3D(K).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
      ELSE
      EPSR = 0.
      END IF

! NO CONDENSATION ONTO RAIN, ONLY EVAP ALLOWED

#ifndef CLDFRAC
           IF (QV3D(K).LT.QVS(K)) THEN
              PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
              PRE(K) = MIN(PRE(K),0.)
           ELSE
              PRE(K) = 0.
           END IF
#else 
!         calcluate QV3D for out of cloud region, assuming water saturaiton in-cloud
           QVCLR3D(K) = (QV3D(K)-QVS(K)*CFL3D_TEMP(K))/(1.-CFL3D_TEMP(K))
!         only calculate if there is some precipitation fraction > cloud fraction
           IF (QVCLR3D(K).LT.QVS(K)) THEN
              PRE(K) = EPSR*(QVCLR3D(K)-QVS(K))/AB(K)
!         only evaporate in out-of-cloud region and distribute across cldmax
              PRE(K) = MIN(PRE(K)*(CLDMAX(K)-CFL3D_TEMP(K)),0.)
              PRE(K) = PRE(K)/CLDMAX(K)
           ELSE
              PRE(K) = 0.
           END IF

! test +++mhwang
!           IF(QC3D(K).GT.QSMALL) PRE(K) = 0.0   ! When QC is present, no rain evaporation
    
#endif

!.......................................................................
! MELTING OF SNOW

! SNOW MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, SNOW MELTS TO FORM RAIN

          IF (QNI3D(K).GE.1.E-8) THEN

! fix 053011
! HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
!             DUM = -CPW/XLF(K)*T3D(K)*PRACS(K)
!             DUM = -CPW/XLF(K)*(T3D(K)-TMELT)*PRACS(K)
             DUM = -CPW/XLF(K)*max(T3D(K)-TMELT, 0.0)*PRACS(K)     !+++mhwang 09/20/2011

!             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*(TMELT-T3D(K))/       &
             PSMLT(K)=2.*PI*N0S(K)*KAP(K)*min(TMELT-T3D(K), 0.0)/       &    !+++mhwang 09/20/2011
                    XLF(K)*RHO(K)*(F1S/(LAMS(K)*LAMS(K))+        &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
                   (LAMS(K)**CONS35))+DUM

! IN WATER SUBSATURATION, SNOW MELTS AND EVAPORATES

#ifndef CLDFRAC
      IF (QVQVS(K).LT.1.) THEN
        EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                   (F1S/(LAMS(K)*LAMS(K))+                       &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
               (LAMS(K)**CONS35))
! bug fix V1.4
        EVPMS(K) = (QV3D(K)-QVS(K))*EPSS/AB(K)    
        EVPMS(K) = MAX(EVPMS(K),PSMLT(K))
        PSMLT(K) = PSMLT(K)-EVPMS(K)
      END IF
#else
      IF (QVCLR3D(K).LT.QVS(K)) THEN
        EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                   (F1S/(LAMS(K)*LAMS(K))+                       &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
               (LAMS(K)**CONS35))
        EVPMS(K) = (QVCLR3D(K)-QVS(K))*EPSS/AB(K)
        EVPMS(K) = min(EVPMS(K)*(CLDMAX(K)-CFL3D_TEMP(K)), 0.0)
        EVPMS(K) = EVPMS(K)/CLDMAX(K)
        EVPMS(K) = MAX(EVPMS(K),PSMLT(K))
        PSMLT(K) = PSMLT(K)-EVPMS(K)
      END IF
#endif
      END IF

!.......................................................................
! MELTING OF GRAUPEL

! GRAUPEL MAY PERSITS ABOVE FREEZING, FORMULA FROM RUTLEDGE AND HOBBS, 1984
! IF WATER SUPERSATURATION, GRAUPEL MELTS TO FORM RAIN

          IF (QG3D(K).GE.1.E-8) THEN

! fix 053011
! HM, MODIFY FOR V3.2, ADD ACCELERATED MELTING DUE TO COLLISION WITH RAIN
!             DUM = -CPW/XLF(K)*T3D(K)*PRACG(K)
!             DUM = -CPW/XLF(K)*(T3D(K)-273.15)*PRACG(K)
             DUM = -CPW/XLF(K)*max(T3D(K)-TMELT, 0.0)*PRACG(K)     !+++mhwang 10/17/2011

             PGMLT(K)=2.*PI*N0G(K)*KAP(K)*(TMELT-T3D(K))/ 		 &
                    XLF(K)*RHO(K)*(F1S/(LAMG(K)*LAMG(K))+                &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
                   (LAMG(K)**CONS36))+DUM

! IN WATER SUBSATURATION, GRAUPEL MELTS AND EVAPORATES

#ifndef CLDFRAC
      IF (QVQVS(K).LT.1.) THEN
        EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                   (F1S/(LAMG(K)*LAMG(K))+                               &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
               (LAMG(K)**CONS36))
! bug fix V1.4
        EVPMG(K) = (QV3D(K)-QVS(K))*EPSG/AB(K)
        EVPMG(K) = MAX(EVPMG(K),PGMLT(K))
        PGMLT(K) = PGMLT(K)-EVPMG(K)
      END IF
#else
      IF (QVCLR3D(K).LT.QVS(K)) THEN
        EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                   (F1S/(LAMG(K)*LAMG(K))+                               &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
               (LAMG(K)**CONS36))
! bug fix V1.4
        EVPMG(K) = (QVCLR3D(K)-QVS(K))*EPSG/AB(K)
        EVPMG(K) = min(EVPMG(K)*(CLDMAX(K)-CFL3D_TEMP(K)), 0.0)
        EVPMG(K) = EVPMG(K)/CLDMAX(K)
        EVPMG(K) = MAX(EVPMG(K),PGMLT(K))
        PGMLT(K) = PGMLT(K)-EVPMG(K)
      END IF
#endif
      END IF

! HM, V3.2
! RESET PRACG AND PRACS TO ZERO, THIS IS DONE BECAUSE THERE IS NO
! TRANSFER OF MASS FROM SNOW AND GRAUPEL TO RAIN DIRECTLY FROM COLLECTION
! ABOVE FREEZING, IT IS ONLY USED FOR ENHANCEMENT OF MELTING AND SHEDDING

      PRACG(K) = 0.
      PRACS(K) = 0.

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! FOR CLOUD ICE, ONLY PROCESSES OPERATING AT T > 273.15 IS
! MELTING, WHICH IS ALREADY CONSERVED DURING PROCESS
! CALCULATION

! CONSERVATION OF QC

      DUM = (PRC(K)+PRA(K))*DT

      IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN

        RATIO = QC3D(K)/DUM

        PRC(K) = PRC(K)*RATIO
        PRA(K) = PRA(K)*RATIO

        END IF

! CONSERVATION OF SNOW

        DUM = (-PSMLT(K)-EVPMS(K)+PRACS(K))*DT

        IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

! NO SOURCE TERMS FOR SNOW AT T > FREEZING
        RATIO = QNI3D(K)/DUM

        PSMLT(K) = PSMLT(K)*RATIO
        EVPMS(K) = EVPMS(K)*RATIO
        PRACS(K) = PRACS(K)*RATIO

        END IF

! CONSERVATION OF GRAUPEL

        DUM = (-PGMLT(K)-EVPMG(K)+PRACG(K))*DT

        IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN

! NO SOURCE TERM FOR GRAUPEL ABOVE FREEZING
        RATIO = QG3D(K)/DUM

        PGMLT(K) = PGMLT(K)*RATIO
        EVPMG(K) = EVPMG(K)*RATIO
        PRACG(K) = PRACG(K)*RATIO

        END IF

! CONSERVATION OF QR
! HM 12/13/06, ADDED CONSERVATION OF RAIN SINCE PRE IS NEGATIVE

#ifndef CLDFRAC
        DUM = (-PRACS(K)-PRACG(K)-PRE(K)-PRA(K)-PRC(K)+PSMLT(K)+PGMLT(K))*DT
#else 
! the fractional area is CLDMAX for PRACS, PRACG, PRE, PSMLT, PGMLT, but is CFL3D for PRA and PRC
! so distribute PRA and PRC across CLDMAX
        DUM = (-PRACS(K)-PRACG(K)-PRE(K)-(PRA(K)+PRC(K))*CFL3D_TEMP(K)/CLDMAX(K)+PSMLT(K)+PGMLT(K))*DT
#endif

        IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN

#ifndef CLDFRAC
        RATIO = (QR3D(K)/DT+PRACS(K)+PRACG(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K))/ &
                        (-PRE(K))
#else
! the fractional area is CLDMAX for PRACS, PRACG, PRE, PSMLT, PGMLT, but is CFL3D for PRA and PRC
! so distribute PRA and PRC across CLDMAX
        RATIO = (QR3D(K)/DT+PRACS(K)+PRACG(K)+(PRA(K)+PRC(K))*CFL3D_TEMP(K)/CLDMAX(K)-PSMLT(K)-PGMLT(K))/ &
                        (-PRE(K))
#endif 
        PRE(K) = PRE(K)*RATIO
        
        END IF

!....................................

#ifndef CLDFRAC
      QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-EVPMS(K)-EVPMG(K))

      T3DTEN(K) = T3DTEN(K)+(PRE(K)*XXLV(K)+(EVPMS(K)+EVPMG(K))*XXLS(K)+&
                    (PSMLT(K)+PGMLT(K)-PRACS(K)-PRACG(K))*XLF(K))/CPM(K)
#else
!     QV3DTEN and T3DTEN are grid-mean tendency, but
!     for hydrometers, they are in-cloud tendency 
      QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-EVPMS(K)-EVPMG(K))*CLDMAX(K)

      T3DTEN(K) = T3DTEN(K)+(PRE(K)*XXLV(K)+(EVPMS(K)+EVPMG(K))*XXLS(K)+&
                    (PSMLT(K)+PGMLT(K)-PRACS(K)-PRACG(K))*XLF(K))/CPM(K)*CLDMAX(K)
#endif 

      QC3DTEN(K) = QC3DTEN(K)+(-PRA(K)-PRC(K))
#ifndef CLDFRAC
      QR3DTEN(K) = QR3DTEN(K)+(PRE(K)+PRA(K)+PRC(K)-PSMLT(K)-PGMLT(K)+PRACS(K)+PRACG(K))
#else
! the fractional area is CLDMAX for PRACS, PRACG, PRE, PSMLT, PGMLT, but is CFL3D for PRA and PRC
! so distribute PRA and PRC across CLDMAX
      QR3DTEN(K) = QR3DTEN(K)+(PRE(K)+(PRA(K)+PRC(K))*CFL3D_TEMP(K)/CLDMAX(K)-PSMLT(K)-PGMLT(K)+PRACS(K)+PRACG(K))
#endif
      QNI3DTEN(K) = QNI3DTEN(K)+(PSMLT(K)+EVPMS(K)-PRACS(K))
      QG3DTEN(K) = QG3DTEN(K)+(PGMLT(K)+EVPMG(K)-PRACG(K))
! fix 053011
!      NS3DTEN(K) = NS3DTEN(K)-NPRACS(K)
! HM, bug fix 5/12/08, npracg is subtracted from nr not ng
!      NG3DTEN(K) = NG3DTEN(K)
      NC3DTEN(K) = NC3DTEN(K)+ (-NPRA(K)-NPRC(K))
#ifndef CLDFRAC
      NR3DTEN(K) = NR3DTEN(K)+ (NPRC1(K)+NRAGG(K)-NPRACG(K))
#else
! The fractional area is CLDMAX for NRAGG, NPRACG, but is CFL3D for NPRC1, 
! so distribute NPRC1 across CLDMAX
      NR3DTEN(K) = NR3DTEN(K)+ (NPRC1(K)*CFL3D_TEMP(K)/CLDMAX(K)+NRAGG(K)-NPRACG(K))
#endif

      IF (PRE(K).LT.0.) THEN
         DUM = PRE(K)*DT/QR3D(K)
           DUM = MAX(-1.,DUM)
         NSUBR(K) = DUM*NR3D(K)/DT
      END IF

! V1.3 move code below to before saturation adjustment
        IF (EVPMS(K)+PSMLT(K).LT.0.) THEN
         DUM = (EVPMS(K)+PSMLT(K))*DT/QNI3D(K)
           DUM = MAX(-1.,DUM)
         NSMLTS(K) = DUM*NS3D(K)/DT
        END IF
        IF (PSMLT(K).LT.0.) THEN
          DUM = PSMLT(K)*DT/QNI3D(K)
          DUM = MAX(-1.0,DUM)
          NSMLTR(K) = DUM*NS3D(K)/DT
        END IF
        IF (EVPMG(K)+PGMLT(K).LT.0.) THEN
         DUM = (EVPMG(K)+PGMLT(K))*DT/QG3D(K)
           DUM = MAX(-1.,DUM)
         NGMLTG(K) = DUM*NG3D(K)/DT
        END IF
        IF (PGMLT(K).LT.0.) THEN
          DUM = PGMLT(K)*DT/QG3D(K)
          DUM = MAX(-1.0,DUM)
          NGMLTR(K) = DUM*NG3D(K)/DT
        END IF

!        nsubr(k)=0.
!        nsubs(k)=0.
!        nsubg(k)=0.

         NS3DTEN(K) = NS3DTEN(K)+(NSMLTS(K))
         NG3DTEN(K) = NG3DTEN(K)+(NGMLTG(K))
         NR3DTEN(K) = NR3DTEN(K)+(NSUBR(K)-NSMLTR(K)-NGMLTR(K))

#ifdef CLDFRAC
! the tendency of hydrometer mixing ratios are grid-mean values
        QC3DTEN(K) = QC3DTEN(K) * CFL3D_TEMP(K)
        QR3DTEN(K) = QR3DTEN(K) * CLDMAX(K)
        QNI3DTEN(K) = QNI3DTEN(K) * CLDMAX(K)
        QG3DTEN(K) = QG3DTEN(K) * CLDMAX(K)
        NC3DTEN(K) = NC3DTEN(K) * CFL3D_TEMP(K)
        NR3DTEN(K) = NR3DTEN(K) * CLDMAX(K)
        NS3DTEN(K) = NS3DTEN(K) * CLDMAX(K)
        NG3DTEN(K) = NG3DTEN(K) * CLDMAX(K)
#endif

 300  CONTINUE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! SUBLIMATE, MELT, OR EVAPORATE NUMBER CONCENTRATION
! THIS FORMULATION ASSUMES 1:1 RATIO BETWEEN MASS LOSS AND
! LOSS OF NUMBER CONCENTRATION

!     IF (PCC(K).LT.0.) THEN
!        DUM = PCC(K)*DT/QC3D(K)
!           DUM = MAX(-1.,DUM)
!        NSUBC(K) = DUM*NC3D(K)/DT
!     END IF

! UPDATE TENDENCIES

!        NC3DTEN(K) = NC3DTEN(K)+NSUBC(K)

!.....................................................................
!.....................................................................
         ELSE  ! TEMPERATURE < 273.15

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

         IF (INUM.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
            NC3D(K)=NDCNST*1.E6/RHO(K)
         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      
! Calculate size distribution of ice, snow, and graupel
!......................................................................
! CLOUD ICE

      IF (QI3D(K).GE.QSMALL) THEN
         LAMI(K) = (CONS12*                 &
              NI3D(K)/QI3D(K))**(1./DI)
         N0I(K) = NI3D(K)*LAMI(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMI(K).LT.LAMMINI) THEN

      LAMI(K) = LAMMINI

      N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      ELSE IF (LAMI(K).GT.LAMMAXI) THEN
      LAMI(K) = LAMMAXI
      N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      END IF
      END IF

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
      N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(K).LT.LAMMINS) THEN
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      ELSE IF (LAMS(K).GT.LAMMAXS) THEN

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)
      END IF
      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
      N0G(K) = NG3D(K)*LAMG(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMG(K).LT.LAMMING) THEN
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      ELSE IF (LAMG(K).GT.LAMMAXG) THEN

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      END IF
      END IF
!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF
      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         !bloss: option for fixing pgam
         if(dofix_pgam) then
            pgam(k) = pgam_fixed
         else

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

         end if
! v1.4
! interpolate
         dumii=int(pgam(k))
         nu(k)=dnu(dumii)+(dnu(dumii+1)-dnu(dumii))* &
               (pgam(k)-real(dumii))

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

! TO CALCULATE DROPLET FREEZING

        CDIST1(K) = NC3D(K)/GAMMA(PGAM(K)+1.)

      END IF

! NO VENTILATION FOR CLOUD ICE

      IF (QI3D(K).GE.QSMALL) THEN

         EPSI = 2.*PI*N0I(K)*RHO(K)*DV(K)/(LAMI(K)*LAMI(K))

      ELSE
         EPSI = 0.
      END IF

      IF (QNI3D(K).GE.QSMALL) THEN
        EPSS = 2.*PI*N0S(K)*RHO(K)*DV(K)*                            &
                   (F1S/(LAMS(K)*LAMS(K))+                       &
                    F2S*(ASN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS10/                   &
               (LAMS(K)**CONS35))
      ELSE
      EPSS = 0.
      END IF

      IF (QG3D(K).GE.QSMALL) THEN
        EPSG = 2.*PI*N0G(K)*RHO(K)*DV(K)*                                &
                   (F1S/(LAMG(K)*LAMG(K))+                               &
                    F2S*(AGN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS11/                   &
               (LAMG(K)**CONS36))
      ELSE
      EPSG = 0.
      END IF

      IF (QR3D(K).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
      ELSE
      EPSR = 0.
      END IF

      PCC(K) = 0.0
#ifndef CLDFRAC
      IF(ISATADJ.EQ.0) THEN !PB 4/13/09

! NOW CALCULATE SATURATION ADJUSTMENT TO CONDENSE EXTRA VAPOR ABOVE
! WATER SATURATION

      DUMT = T3D(K)
      DUMQV = QV3D(K)
! hm, add fix for low pressure, 5/12/10
      dum=min(0.99*pres(k),POLYSVP(DUMT,0))
      DUMQSS = EP_2*dum/(PRES(K)-dum)
      DUMQC = QC3D(K)
      DUMQC = MAX(DUMQC,0.)

! SATURATION ADJUSTMENT FOR LIQUID

      DUMS = DUMQV-DUMQSS
      PCC(K) = DUMS/(1.+XXLV(K)**2*DUMQSS/(CPM(K)*RV*DUMT**2))/DT
!      IF (PCC(K)*DT+DUMQC.LT.0.) THEN
!           PCC(K) = -DUMQC/DT
!      END IF
!+++mhwang
      IF (PCC(K)*DT+QC3D(K).LT.0.) THEN
           PCC(K) = -(QC3D(K))/DT
      END IF
!---mhwang

!      QV3DTEN(K) = QV3DTEN(K)-PCC(K)
!      T3DTEN(K) = T3DTEN(K)+PCC(K)*XXLV(K)/CPM(K)
!      QC3DTEN(K) = QC3DTEN(K)+PCC(K)
       QV3D(K) = QV3D(K)-PCC(K)*DT
       QC3D(K) = QC3D(K)+PCC(K)*DT
       T3D(K) = T3D(K)+PCC(K)*XXLV(K)/CPM(K)*DT
      END IF
#endif

      QC3DTEN(K) = 0.0
      NC3DTEN(K) = 0.0
      T3DTEN(K) = 0.0
      NCACT(K) = 0.0

!.......................................................................
! ACTIVATION OF CLOUD DROPLETS

!bloss: only do activation if droplet number is predicted
!bloss      IF (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL) THEN
      IF (QC3D(K)+QC3DTEN(K)*DT.GE.QSMALL.AND.INUM.EQ.0) THEN

! EFFECTIVE VERTICAL VELOCITY (M/S)

      IF (ISUB.EQ.0) THEN
! ADD SUB-GRID VERTICAL VELOCITY
         DUM = W3D(K)+WVAR(K)

! ASSUME MINIMUM EFF. SUB-GRID VELOCITY 0.10 M/S
         DUM = MAX(DUM,0.10)

      ELSE IF (ISUB.EQ.1) THEN
         DUM=W3D(K)
      END IF

! ONLY ACTIVATE IN REGIONS OF UPWARD MOTION
      IF (DUM.GE.0.001) THEN

      IF (IBASE.EQ.1) THEN

! ACTIVATE ONLY IF THERE IS LITTLE CLOUD WATER
! OR IF AT CLOUD BASE, OR AT LOWEST MODEL LEVEL (K=1)

         IDROP=0

! V1.3 USE CURRENT VALUE OF QC FOR IDROP
         IF (QC3D(K).LE.0.05E-3/RHO(K)) THEN
            IDROP=1
         END IF
         IF (K.EQ.1) THEN
            IDROP=1
         ELSE IF (K.GE.2) THEN
            IF (QC3D(K).GT.0.05E-3/RHO(K).AND. &
             QC3D(K-1).LE.0.05E-3/RHO(K-1)) THEN
            IDROP=1
            END IF
         END IF

         IF (IDROP.EQ.1) THEN
! ACTIVATE AT CLOUD BASE OR REGIONS WITH VERY LITTLE LIQ WATER

           IF (IACT.EQ.1) THEN
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

            DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
            DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else 
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif

           ELSE IF (IACT.EQ.2) THEN
! DROPLET ACTIVATION FROM ABDUL-RAZZAK AND GHAN (2000)

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else 
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#if (defined CRM && defined MODAL_AERO)
            ELSE if (IACT.EQ.3) then
              INES = 0
              CALL DROP_ACTIVATION_GHAN(DUM, T3D(k), RHO(k),  &
                   DUM2, INES, SMAX, K)
              DUM2 = (DUM2-NC3D(K))/DT
              DUM2 = MAX(0., DUM2)
#ifndef CLDFRAC
              NC3DTEN(K) = NC3DTEN(K)+DUM2
#else 
              NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#endif
           END IF  ! IACT

!.............................................................................
        ELSE IF (IDROP.EQ.0) THEN
! ACTIVATE IN CLOUD INTERIOR
! FIND EQUILIBRIUM SUPERSATURATION

           TAUC=1./(2.*PI*RHO(k)*DV(K)*NC3D(K)*(PGAM(K)+1.)/LAMC(K))
           IF (EPSR.GT.1.E-8) THEN
             TAUR=1./EPSR
           ELSE
             TAUR=1.E8
           END IF
           IF (EPSI.GT.1.E-8) THEN
             TAUI=1./EPSI
           ELSE
             TAUI=1.E8
           END IF
           IF (EPSS.GT.1.E-8) THEN
             TAUS=1./EPSS
           ELSE
             TAUS=1.E8
           END IF
           IF (EPSG.GT.1.E-8) THEN
             TAUG=1./EPSG
           ELSE
             TAUG=1.E8
           END IF

! EQUILIBRIUM SS INCLUDING BERGERON EFFECT

           DUM3=(QVS(K)*RHO(K)/(PRES(K)-EVS(K))+DQSDT/CP)*G*DUM
           DUM3=(DUM3*TAUC*TAUR*TAUI*TAUS*TAUG- &
           (QVS(K)-QVI(K))*(TAUC*TAUR*TAUI*TAUG+TAUC*TAUR*TAUS*TAUG+TAUC*TAUR*TAUI*TAUS))/ &
           (TAUC*TAUR*TAUI*TAUG+TAUC*TAUR*TAUS*TAUG+TAUC*TAUR*TAUI*TAUS+ &
            TAUR*TAUI*TAUS*TAUG+TAUC*TAUI*TAUS*TAUG)

           IF (DUM3/QVS(K).GE.1.E-6) THEN
           IF (IACT.EQ.1) THEN

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

            DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
            DUMACT = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))

! USE POWER LAW CCN SPECTRA

! CONVERT FROM ABSOLUTE SUPERSATURATION TO SUPERSATURATION RATIO IN %
            DUM3=DUM3/QVS(K)*100.

            DUM2=C1*DUM3**K1
! MAKE SURE VALUE DOESN'T EXCEED THAT FOR NON-EQUILIBRIUM SS
            DUM2=MIN(DUM2,DUMACT)
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else 
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif

           ELSE IF (IACT.EQ.2) THEN

! FIND MAXIMUM ALLOWED ACTIVATION WITH NON-EQULIBRIUM SS

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

           DUMACT = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

! USE LOGNORMAL AEROSOL
           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

! GET SUPERSATURATION RATIO FROM ABSOLUTE SUPERSATURATION
           SMAX = DUM3/QVS(K)

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

! MAKE SURE ISN'T GREATER THAN NON-EQUIL. SS
            DUM2=MIN(DUM2,DUMACT)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else 
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#if (defined CRM && defined MODAL_AERO)
            ELSE if (IACT.EQ.3) then
! GET SUPERSATURATION RATIO FROM ABSOLUTE SUPERSATURATION
              SMAX = DUM3/QVS(K)
              INES = 1 
              CALL DROP_ACTIVATION_GHAN(DUM, T3D(k), RHO(k),  &
                   DUM2, INES, SMAX, K)
              DUM2 = (DUM2-NC3D(K))/DT
              DUM2 = MAX(0., DUM2)
#ifndef CLDFRAC
              NC3DTEN(K) = NC3DTEN(K)+DUM2
#else 
              NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#endif

           END IF ! IACT
           END IF ! DUM3/QVS > 1.E-6
        END IF  ! IDROP = 1

!.......................................................................
      ELSE IF (IBASE.EQ.2) THEN

           IF (IACT.EQ.1) THEN
! USE ROGERS AND YAU (1989) TO RELATE NUMBER ACTIVATED TO W
! BASED ON TWOMEY 1959

            DUM=DUM*100.  ! CONVERT FROM M/S TO CM/S
            DUM2 = 0.88*C1**(2./(K1+2.))*(7.E-2*DUM**1.5)**(K1/(K1+2.))
            DUM2=DUM2*1.E6 ! CONVERT FROM CM-3 TO M-3
            DUM2=DUM2/RHO(K)  ! CONVERT FROM M-3 TO KG-1
            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else 
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif

           ELSE IF (IACT.EQ.2) THEN

           SIGVL = 0.0761-1.55E-4*(T3D(K)-TMELT)
           AACT = 2.*MW/(RHOW*RR)*SIGVL/T3D(K)
           ALPHA = G*MW*XXLV(K)/(CPM(K)*RR*T3D(K)**2)-G*MA/(RR*T3D(K))
           GAMM = RR*T3D(K)/(EVS(K)*MW)+MW*XXLV(K)**2/(CPM(K)*PRES(K)*MA*T3D(K))

           GG = 1./(RHOW*RR*T3D(K)/(EVS(K)*DV(K)*MW)+ XXLV(K)*RHOW/(KAP(K)*T3D(K))*(XXLV(K)*MW/ &
              (T3D(K)*RR)-1.))

           PSI = 2./3.*(ALPHA*DUM/GG)**0.5*AACT

           ETA1 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW1)
           ETA2 = (ALPHA*DUM/GG)**1.5/(2.*PI*RHOW*GAMM*NANEW2)

           SM1 = 2./BACT**0.5*(AACT/(3.*RM1))**1.5
           SM2 = 2./BACT**0.5*(AACT/(3.*RM2))**1.5

           DUM1 = 1./SM1**2*(F11*(PSI/ETA1)**1.5+F21*(SM1**2/(ETA1+3.*PSI))**0.75)
           DUM2 = 1./SM2**2*(F12*(PSI/ETA2)**1.5+F22*(SM2**2/(ETA2+3.*PSI))**0.75)

           SMAX = 1./(DUM1+DUM2)**0.5

           UU1 = 2.*LOG(SM1/SMAX)/(4.242*LOG(SIG1))
           UU2 = 2.*LOG(SM2/SMAX)/(4.242*LOG(SIG2))
           DUM1 = NANEW1/2.*(1.-DERF1(UU1))
           DUM2 = NANEW2/2.*(1.-DERF1(UU2))

           DUM2 = (DUM1+DUM2)/RHO(K)  !CONVERT TO KG-1

! MAKE SURE THIS VALUE ISN'T GREATER THAN TOTAL NUMBER OF AEROSOL

            DUM2 = MIN((NANEW1+NANEW2)/RHO(K),DUM2)

            DUM2 = (DUM2-NC3D(K))/DT
            DUM2 = MAX(0.,DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else 
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#if (defined CRM && defined MODAL_AERO)
            ELSE if (IACT.EQ.3) then
              INES = 0
              CALL DROP_ACTIVATION_GHAN(DUM, T3D(k), RHO(k),  &
                   DUM2, INES, SMAX, K)
              DUM2 = (DUM2-NC3D(K))/DT
              DUM2 = MAX(0., DUM2)
#ifndef CLDFRAC
            NC3DTEN(K) = NC3DTEN(K)+DUM2
#else 
            NC3DTEN(K) = NC3DTEN(K)+DUM2 * CFL3D_TEMP(K)
#endif
#endif
           END IF  ! IACT
        END IF  ! IBASE
        END IF  ! W > 0.001
#ifndef CLDFRAC
          NC3D(K) = NC3D(K) + NC3DTEN(K)*DT
#else
          NC3D(K) = (NC3D(K)* CFL3D_TEMP(K)+ NC3DTEN(K) * DT)/CFL3D_TEMP(K)
#endif
          NCACT(K) = NC3DTEN(K)
          NC3DTEN(K) = 0.0
        END IF  ! QC3D > QSMALL



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! NUCLEATION OF CLOUD ICE FROM HOMOGENEOUS AND HETEROGENEOUS FREEZING ON AEROSOL

        NNUCCD0(K) = 0.0
        MNUCCD0(K) = 0.0
#ifdef  CLDFRAC  ! the new scheme with fractional cloudiness
         IF (INUC.EQ.0) THEN

! FREEZING OF AEROSOL ONLY ALLOWED BELOW -5 C
! AND ABOVE DELIQUESCENCE THRESHOLD OF 80%
! AND ABOVE ICE SATURATION

! add threshold according to Greg Thomspon
          kc2_inliq = 0.0
          kc2 = 0.0
! For ice nucleation inside liquid clouds
          IF(CFL3D_TEMP(K).gt.cloud_frac_thresh) then
            IF(T3D(K).le.265.15 .or. QVS(K)/QVI(K).ge.1.08) then 
! hm, modify dec. 5, 2006, replace with cooper curve
             kc2_inliq = 0.005*exp(0.304*(TMELT-T3D(K)))*1000. ! convert from L-1 to m-3
! limit to 500 L-1
             kc2_inliq = min(kc2_inliq,500.e3)
             kc2_inliq=MAX(kc2_inliq/rho(k),0.)  ! convert to kg-1
            ENDIF
          ENDIF
          
! for ice nucleation in ice clouds outside of liquid clouds
          if ((QVQVS(K).GE.0.999.and.T3D(K).le.265.15).or. &
              QVQVSI(K).ge.1.08) then
! hm, modify dec. 5, 2006, replace with cooper curve
           kc2 = 0.005*exp(0.304*(TMELT-T3D(K)))*1000. ! convert from L-1 to m-3
! limit to 500 L-1
           kc2 = min(kc2,500.e3)
           kc2=MAX(kc2/rho(k),0.)  ! convert to kg-1
          end if

          kc2 = (kc2_inliq * CFL3D_TEMP(K) + kc2*(CFI3D_TEMP(K)-CFL3D_TEMP(K)))/CFI3D_TEMP(K)

          IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
             NNUCCD0(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
             MNUCCD0(K) = NNUCCD0(K)*MI0
          END IF

         ELSE IF (INUC.EQ.1) THEN

          KC2_INLIQ = 0.0
          KC2 = 0.0
! For ice nucleation inside liquid clouds
          IF (T3D(K).LT.TMELT .and. CFL3D_TEMP(K).GT.cloud_frac_thresh) THEN 
              KC2_INLIQ = 0.16*1000./RHO(K)  ! CONVERT FROM L-1 TO KG-1
          END IF
! For ice nucleation outside liquid clouds
          IF (T3D(K).LT.TMELT.AND.QVQVSI(K).GT.1.) THEN
             KC2 = 0.16*1000./RHO(K)  ! CONVERT FROM L-1 TO KG-1
          ENDIF 
          KC2 = (KC2_INLIQ * CFL3D_TEMP(K) + KC2*(CFI3D_TEMP(K)-CFL3D_TEMP(K)))/CFI3D_TEMP(K)
          IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
             NNUCCD0(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
             MNUCCD0(K) = NNUCCD0(K)*MI0
          END IF

         END IF
#else  ! the original scheme 
         IF (INUC.EQ.0) THEN

! FREEZING OF AEROSOL ONLY ALLOWED BELOW -5 C
! AND ABOVE DELIQUESCENCE THRESHOLD OF 80%
! AND ABOVE ICE SATURATION

! add threshold according to Greg Thomspon

         if ((QVQVS(K).GE.0.999.and.T3D(K).le.265.15).or. &
              QVQVSI(K).ge.1.08) then

! hm, modify dec. 5, 2006, replace with cooper curve
      kc2 = 0.005*exp(0.304*(TMELT-T3D(K)))*1000. ! convert from L-1 to m-3
! limit to 500 L-1
      kc2 = min(kc2,500.e3)
      kc2=MAX(kc2/rho(k),0.)  ! convert to kg-1

          IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
             NNUCCD0(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
             MNUCCD0(K) = NNUCCD0(K)*MI0
          END IF

          END IF

          ELSE IF (INUC.EQ.1) THEN

          IF (T3D(K).LT.TMELT.AND.QVQVSI(K).GT.1.) THEN

             KC2 = 0.16*1000./RHO(K)  ! CONVERT FROM L-1 TO KG-1
          IF (KC2.GT.NI3D(K)+NS3D(K)+NG3D(K)) THEN
             NNUCCD0(K) = (KC2-NI3D(K)-NS3D(K)-NG3D(K))/DT
             MNUCCD0(K) = NNUCCD0(K)*MI0
          END IF
          END IF

         END IF
#endif  

#ifdef CLUBB_CRM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the case of clex9_oct14, we need to decrease the ice    !
! nucleation in order for the cloud to persist for realistic  !
! lengths. It is suggested to reduce by a factor of 100       !
! This coefficent can be changed in subroutine init_microphys !
! in the microphys_driver subroutine.                         !
!                                                             !
        NNUCCD0(K)=NNUCCD0(K)*NNUCCD_REDUCE_COEF
!
! Change made by Marc Pilon on 11/16/11                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif /* CLUBB */

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP/SUB/DEP TERMS FOR QI,QNI,QR



#ifndef CLDFRAC
! ONLY INCLUDE REGION OF ICE SIZE DIST < DCS
! DUM IS FRACTION OF D*N(D) < DCS

! LOGIC BELOW FOLLOWS THAT OF HARRINGTON ET AL. 1995 (JAS)
              IF (QI3D(K).GE.QSMALL) THEN              
              DUM=(1.-EXP(-LAMI(K)*DCS)*(1.+LAMI(K)*DCS))
              PRD0(K) = EPSI*(QV3D(K)-QVI(K))/ABI(K)*DUM
              ELSE
              DUM=0.
              END IF
! ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
              IF (QNI3D(K).GE.QSMALL) THEN
              PRDS0(K) = EPSS*(QV3D(K)-QVI(K))/ABI(K)+ &
                EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
! OTHERWISE ADD TO CLOUD ICE
              ELSE
              PRD0(K) = PRD0(K)+EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)
              END IF
! VAPOR DPEOSITION ON GRAUPEL
              PRDG0(K) = EPSG*(QV3D(K)-QVI(K))/ABI(K)
#else
! ONLY INCLUDE REGION OF ICE SIZE DIST < DCS
! DUM IS FRACTION OF D*N(D) < DCS

! LOGIC BELOW FOLLOWS THAT OF HARRINGTON ET AL. 1995 (JAS)
              IF (QI3D(K).GE.QSMALL) THEN

! ADD DEPOSITION IN TAIL OF ICE SIZE DIST TO SNOW IF SNOW IS PRESENT
              IF (QNI3D(K).GE.QSMALL) THEN
              DUM=(1.-EXP(-LAMI(K)*DCS)*(1.+LAMI(K)*DCS))
!  Fractional area for snow is CLDMAX, while it is CFI3D for ice
! inside liquid clouds, use saturaiton vapor pressure over liquid 
! First, calculate B-F processes (conversion from liquid to ice)
              BERFDS0(K) = min(EPSS*(QVS(K)-QVI(K))/ABI(K)+EPSI*(QVS(K)-QVI(K))/ABI(K)*(1.-DUM), QC3D(K)/DT)
! test set BERFDS to be zero
!              BERFDS(K) = 0.0  
! Second, calculate the vapor deposition after liquid water is depleted inside liquid clouds
              PRDS0(K) = (EPSS*(QVS(K)-QVI(K))/ABI(K)+EPSI*(QVS(K)-QVI(K))/ABI(K)*(1.-DUM) - BERFDS0(K))*CFL3D(K)
!              IF(CFL3D_TEMP(K).LT.0.01 .OR. QV3D(K).GT.QVI(K)) THEN  ! IF CFL3D > 0.01, no evaporation
              IF(CFL3D_TEMP(K).LT.0.9 .OR. QV3D(K).GT.QVI(K)) THEN  ! IF CFL3D > 0.90, no evaporation
              PRDS0(K) = (PRDS0(K)+EPSS*(QV3D(K)-QVI(K))/ABI(K)*(CLDMAX(K)-CFL3D_TEMP(K))+ &
                        EPSI*(QV3D(K)-QVI(K))/ABI(K)*(1.-DUM)*(CFI3D_TEMP(K)-CFL3D_TEMP(K)))/CLDMAX(K)
! OTHERWISE ADD TO CLOUD ICE
              ELSE 
              PRDS0(K) = PRDS0(K)/CLDMAX(K)
              ENDIF 
              ELSE
              DUM=1.0
              ENDIF 

!              PRD(K) = EPSI*(QV3D(K)-QVI(K))/ABI(K)*DUM
! Inside liquid clouds, use saturation varpor pressure over liquid 
!
! First, calculate B-F process (conversion from liquid to ice)
              BERFD0(K) = min(EPSI*(QVS(K)-QVI(K))/ABI(K)*DUM, QC3D(K)/DT)
! test set BERFDS to be zero
!              BERFD(K) = 0.0
! Second, the vapor postion to ice after liquid water is depleted: 
              PRD0(K) = (EPSI*(QVS(K)-QVI(K))/ABI(K)*DUM-BERFD0(K)) * CFL3D(K)           

              PRD1(K) = PRD0(K) ! +++mhwang test

! Outside liquid clouds, use in-situ saturation vapor pressure
!              IF(CFL3D_TEMP(K).LT.0.01.OR. QV3D(K).GT.QVI(K)) THEN  ! IF CFL3D > 0.01, no evaporation
              IF(CFL3D_TEMP(K).LT.0.90.OR. QV3D(K).GT.QVI(K)) THEN  ! IF CFL3D > 0.90, no evaporation
              PRD0(K) = (PRD0(K)+(EPSI*(QV3D(K)-QVI(K))/ABI(K)*DUM) * (CFI3D_TEMP(K)-CFL3D_TEMP(K)))/CFI3D_TEMP(K)   
              ELSE
              PRD0(K) = PRD0(K)/CFI3D_TEMP(K)
              ENDIF
              PRD2(K) = PRD0(K) ! +++mhwang test
              END IF

! VAPOR DPEOSITION ON GRAUPEL
!              PRDG(K) = EPSG*(QV3D(K)-QVI(K))/ABI(K)
! Inside liquid clouds use saturaiton varpor pressure over liquid 
! First, calculate B-F process (conversion from liquid to ice)
              BERFDG0(K) = min(EPSG*(QVS(K)-QVI(K))/ABI(K), QC3D(K)/DT)
! test set BERFDS to be zero
!              BERFDG(K) = 0.0
! Second, the vapor postion to ice after liquid water is depleted:
              PRDG0(K) = (EPSG*(QVS(K)-QVI(K))/ABI(K)-BERFDG0(K)) * CFL3D(K)
!              IF(CFL3D_TEMP(K).LT.0.01.OR. QV3D(K).GT.QVI(K)) THEN  ! IF CFL3D > 0.01, no evaporation
              IF(CFL3D_TEMP(K).LT.0.90.OR. QV3D(K).GT.QVI(K)) THEN  ! IF CFL3D > 0.90, no evaporation
              PRDG0(K) = (PRDG0(K)+(EPSG*(QV3D(K)-QVI(K))/ABI(K)*(CLDMAX(K)-CFL3D_TEMP(K))))/CLDMAX(K)
              ELSE
              PRDG0(K) = PRDG0(K)/CLDMAX(K)
              ENDIF

! Test, if QC is present, no evaporation of ice, snow, and graupel
!              IF(QC3D(K).GT.QSMALL) THEN 
!                 IF(PRD(K).LT.0.0) PRD(K) = 0.0
!                 IF(PRDS(K).LT.0.0) PRDS(K) = 0.0
!                 IF(PRDG(K).LT.0.0) PRDG(K) = 0.0
!              ENDIF  

! Make sure B-F process does not exceed liquid water content
              DUM = (BERFD0(K)+BERFDS0(K)+BERFDG0(K)) * DT
              IF(DUM.gt.QC3D(K) .and. QC3D(K).GT.QSMALL) THEN
                BERFD0(K) = BERFD0(K) * QC3D(K)/DUM
                BERFDS0(K) = BERFDS0(K) * QC3D(K)/DUM
                BERFDG0(K) = BERFDG0(K) * QC3D(K)/DUM
              END IF 
#endif

! MAKE SURE NOT PUSHED INTO ICE SUPERSAT/SUBSAT
! FORMULA FROM REISNER 2 SCHEME

#ifndef CLDFRAC
           DUM = (QV3D(K)-QVI(K))/DT

           FUDGEF = 0.9999
           SUM_DEP = PRD0(K)+PRDS0(K)+MNUCCD0(K)+PRDG0(K)
           IF( (DUM.GT.0. .AND. SUM_DEP.GT.DUM*FUDGEF) .OR.                      &
               (DUM.LT.0. .AND. SUM_DEP.LT.DUM*FUDGEF) ) THEN
               MNUCCD0(K) = FUDGEF*MNUCCD0(K)*DUM/SUM_DEP
               PRD0(K) = FUDGEF*PRD0(K)*DUM/SUM_DEP
               PRDS0(K) = FUDGEF*PRDS0(K)*DUM/SUM_DEP
	       PRDG0(K) = FUDGEF*PRDG0(K)*DUM/SUM_DEP
           ENDIF
#else 

           FUDGEF = 0.9999
! with fractional cloudiness, QV3D is the grid-mean value, so do sum_dep
           SUM_DEP = ((PRD0(K)+MNUCCD0(K))*CFI3D_TEMP(K)+(PRDS0(K)+PRDG0(K))*CLDMAX(K))*DT

! For net evaporation
           IF (SUM_DEP.LT. (-1.*QSMALL)) THEN
               DUM = min((QV3D(K)-QVI(K)), 0.0)*CLDMAX(K)
               IF(SUM_DEP.LT.DUM) THEN 
                 MNUCCD0(K) = FUDGEF*MNUCCD0(K)*DUM/SUM_DEP
                 PRD0(K) = FUDGEF*PRD0(K)*DUM/SUM_DEP
                 PRDS0(K) = FUDGEF*PRDS0(K)*DUM/SUM_DEP
                 PRDG0(K) = FUDGEF*PRDG0(K)*DUM/SUM_DEP
               ENDIF 
           ENDIF
! For net depositon
           IF (SUM_DEP.GT.QSMALL) THEN
               IF(CFL3D(K).gt.1.0e-5) then
                 DUM = ((QVS(K)-QVI(K))*CFL3D(K)+(max((QV3D(K)-QVI(K)), 0.0)*(CLDMAX(K)-CFL3D_TEMP(K))))
                 DUM = min(max((QV3D(K)-0.5*QVI(K)), 0.0)*CLDMAX(K), DUM) ! when CFL3D > 0, vapor deposition does not 
                                                                          ! allow to reduce the grid-mean ice saturaiton below 0.5
               ELSE
                 DUM = max((QV3D(K)-QVI(K)), 0.0)*CLDMAX(K) 
               ENDIF
               IF(SUM_DEP.GT.DUM) THEN
                 MNUCCD0(K) = FUDGEF*MNUCCD0(K)*DUM/SUM_DEP
                 PRD0(K) = FUDGEF*PRD0(K)*DUM/SUM_DEP
                 PRDS0(K) = FUDGEF*PRDS0(K)*DUM/SUM_DEP
                 PRDG0(K) = FUDGEF*PRDG0(K)*DUM/SUM_DEP
               ENDIF
           ENDIF

           PRD3(K) = PRD(K) ! +++mhwang test
#endif

! IF CLOUD ICE/SNOW/GRAUPEL VAP DEPOSITION IS NEG, THEN ASSIGN TO SUBLIMATION PROCESSES

           IF (PRD0(K).LT.0.) THEN
              EPRD0(K)=PRD0(K)
              PRD0(K)=0.
           END IF
           IF (PRDS0(K).LT.0.) THEN
              EPRDS0(K)=PRDS0(K)
              PRDS0(K)=0.
           END IF
           IF (PRDG0(K).LT.0.) THEN
              EPRDG0(K)=PRDG0(K)
              PRDG0(K)=0.
           END IF

      IF (EPRD0(K).LT.0.) THEN
         DUM = EPRD0(K)*DT/QI3D(K)
            DUM = MAX(-1.,DUM)
         NSUBI0(K) = DUM*NI3D(K)/DT
      END IF
      IF (EPRDS0(K).LT.0.) THEN
         DUM = EPRDS0(K)*DT/QNI3D(K)
           DUM = MAX(-1.,DUM)
         NSUBS0(K) = DUM*NS3D(K)/DT
      END IF
      IF (EPRDG0(K).LT.0.) THEN
         DUM = EPRDG0(K)*DT/QG3D(K)
           DUM = MAX(-1.,DUM)
         NSUBG0(K) = DUM*NG3D(K)/DT
      END IF

#ifdef CLDFRAC
      IF (BERFD0(K)+BERFDS0(K)+BERFDG0(K).GT.0) THEN
          DUM = -1. * (BERFD0(K)+BERFDS0(K)+BERFDG0(K))*DT/QC3D(K)
          DUM = MAX(-1., DUM)
          NSUBC0(K) = DUM * NC3D(K)/DT
      ENDIF
#endif

! ****SENSITIVITY - NO ICE
      IF (ILIQ.EQ.1) THEN
        MNUCCD0(K)=0.
        NNUCCD0(K)=0.
        PRD0(K) = 0.0
        EPRD0(K) = 0.0
        PRDS0(K) = 0.0
        EPRDS0(K) = 0.0
        NSUBI0(K) =0.0
        NSUBS0(K) = 0.0
#ifdef CLDFRAC
        BERFD0(K) = 0.0
        BERFDS0(K) = 0.0
        NSUBC0(K) = 0.0
#endif
      END IF
! ****SENSITIVITY - NO GRAUPEL
      IF (IGRAUP.EQ.1) THEN
            PRDG0(K) = 0.
            EPRDG0(K) = 0.
            NSUBG0(K) = 0.
#ifdef CLDFRAC
            BERFDG0(K) = 0.
#endif
       END IF


!  Update QV, T, QC, NC, QI, NI, QS, NS, QG, NG
#ifdef CLDFRAC
           QG3D(K) = (QG3D(K)*CLDMAX(K)+(EPRDG0(K)+PRDG0(K))*DT*CLDMAX(K)  &
                     +BERFDG0(K)*CFL3D_TEMP(K)*DT)/CLDMAX(K)
           NG3D(K) = (NG3D(K)*CLDMAX(K)+NSUBG0(K)*DT*CLDMAX(K))/CLDMAX(K)
           QNI3D(K) = (QNI3D(K)*CLDMAX(K)+(EPRDS0(K)+PRDS0(K))*DT*CLDMAX(K)  &
                     +BERFDS0(K)*CFL3D_TEMP(K)*DT)/CLDMAX(K)
           NS3D(K) = (NS3D(K)*CLDMAX(K)+NSUBS0(K)*DT*CLDMAX(K))/CLDMAX(K)
           QI3D(K) = (QI3D(K)*CFI3D_TEMP(K)+(EPRD0(K)+PRD0(K))*DT*CFI3D_TEMP(K)  &
                     +BERFD0(K)*CFL3D_TEMP(K)*DT+MNUCCD0(K)*CFI3D_TEMP(K)*DT)/CFI3D_TEMP(K)
           NI3D(K) = (NI3D(K)*CFI3D_TEMP(K)+NSUBI0(K)*DT*CFI3D_TEMP(K)  &
                     +NNUCCD0(K)*CFI3D_TEMP(K)*DT)/CFI3D_TEMP(K)
           QC3D(K) = (QC3D(K)*CFL3D_TEMP(K)-(BERFD0(K)+BERFDS0(K)+BERFDG0(K))*CFL3D_TEMP(K)*DT)/CFL3D_TEMP(K)
           QV3D(K) = QV3D(K)-(MNUCCD0(K)+EPRD0(K)+PRD0(K))*DT*CFI3D_TEMP(K)    &
                     -(EPRDG0(K)+PRDG0(K)+EPRDS0(K)+PRDS0(K))*DT*CLDMAX(K)
           NC3D(K) = NC3D(K)+NSUBC0(K)*DT
           T3D(K) = T3D(K)+((MNUCCD0(K)*CFI3D_TEMP(K)+(EPRDG0(K)+PRDG0(K)+EPRDS0(K)+PRDS0(K))*CLDMAX(K)  &
                    +(EPRD0(K)+PRD0(K))*CFI3D_TEMP(K))*DT*XXLS(K)  &
                    +(BERFD0(K)+BERFDS0(K)+BERFDG0(K))*CFL3D_TEMP(K)*DT*XLF(K))/CPM(K) 
          
           QC3DTEN0(K) = -(BERFD0(K)+BERFDS0(K)+BERFDG0(K))*CFL3D_TEMP(K)
           QI3DTEN0(K) = (EPRD0(K)+PRD0(K))*CFI3D_TEMP(K)  &
                     +BERFD0(K)*CFL3D_TEMP(K)+MNUCCD0(K)*CFI3D_TEMP(K)
           QNI3DTEN0(K) = (EPRDS0(K)+PRDS0(K))*CLDMAX(K)  &
                     +BERFDS0(K)*CFL3D_TEMP(K)
           QG3DTEN0(K) = (EPRDG0(K)+PRDG0(K))*CLDMAX(K)  &
                     +BERFDG0(K)*CFL3D_TEMP(K)
           QV3DTEN0(K) = -(MNUCCD0(K)+EPRD0(K)+PRD0(K))*CFI3D_TEMP(K)    &
                     -(EPRDG0(K)+PRDG0(K)+EPRDS0(K)+PRDS0(K))*CLDMAX(K)
           T3DTEN0(K) = ((MNUCCD0(K)*CFI3D_TEMP(K)+(EPRDG0(K)+PRDG0(K)+EPRDS0(K)+PRDS0(K))*CLDMAX(K)  &
                    +(EPRD0(K)+PRD0(K))*CFI3D_TEMP(K))*XXLS(K)  &
                    +(BERFD0(K)+BERFDS0(K)+BERFDG0(K))*CFL3D_TEMP(K)*XLF(K))/CPM(K)
#else
           QG3D(K) = (QG3D(K)+(EPRDG0(K)+PRDG0(K))*DT)
           NG3D(K) = (NG3D(K)+NSUBG0(K)*DT)
           QNI3D(K) = (QNI3D(K)+(EPRDS0(K)+PRDS0(K))*DT)
           NS3D(K) = (NS3D(K)+NSUBS0(K)*DT)
           QI3D(K) = (QI3D(K)+(EPRD0(K)+PRD0(K)+MNUCCD0(K))*DT)
           NI3D(K) = (NI3D(K)+(NSUBI0(K)+NNUCCD0(K))*DT)
           QV3D(K) = QV3D(K)-(MNUCCD0(K)+EPRD0(K)+PRD0(K)+EPRDS0(K)+PRDS0(K)  &
                     +EPRDG0(K)+PRDG0(K))*DT
           T3D(K) = T3D(K)+(MNUCCD0(K)+   &
                    +EPRDG0(K)+PRDG0(K)+EPRDS0(K)+PRDS0(K)+EPRD0(K)+PRD0(K))*DT*XXLS(K)/CPM(K)
#endif

!......................................................................
!HM ADD, ALLOW FOR CONSTANT DROPLET NUMBER
! INUM = 0, PREDICT DROPLET NUMBER
! INUM = 1, SET CONSTANT DROPLET NUMBER

      IF (INUM.EQ.1) THEN
! CONVERT NDCNST FROM CM-3 TO KG-1
         NC3D(K)=NDCNST*1.E6/RHO(K)
      END IF

! CALCULATE SIZE DISTRIBUTION PARAMETERS
! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(K) = MAX(0.,NI3D(K))
      NS3D(K) = MAX(0.,NS3D(K))
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))

!......................................................................
! CLOUD ICE

      IF (QI3D(K).GE.QSMALL) THEN
         LAMI(K) = (CONS12*                 &
              NI3D(K)/QI3D(K))**(1./DI)
         N0I(K) = NI3D(K)*LAMI(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMI(K).LT.LAMMINI) THEN

      LAMI(K) = LAMMINI

      N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      ELSE IF (LAMI(K).GT.LAMMAXI) THEN
      LAMI(K) = LAMMAXI
      N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      END IF
      END IF

!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)
      N0RR(K) = NR3D(K)*LAMR(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF
      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         !bloss: option for fixing pgam
         if(dofix_pgam) then
            pgam(k) = pgam_fixed
         else

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

         end if
! v1.4
! interpolate
         dumii=int(pgam(k))
         nu(k)=dnu(dumii)+(dnu(dumii+1)-dnu(dumii))* &
               (pgam(k)-real(dumii))

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN

      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26
      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

! TO CALCULATE DROPLET FREEZING

        CDIST1(K) = NC3D(K)/GAMMA(PGAM(K)+1.)

      END IF

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)
      N0S(K) = NS3D(K)*LAMS(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(K).LT.LAMMINS) THEN
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      ELSE IF (LAMS(K).GT.LAMMAXS) THEN

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)
      END IF
      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)
      N0G(K) = NG3D(K)*LAMG(K)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMG(K).LT.LAMMING) THEN
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      ELSE IF (LAMG(K).GT.LAMMAXG) THEN

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      END IF
      END IF

!.....................................................................
! ZERO OUT PROCESS RATES

            MNUCCC(K) = 0.
            NNUCCC(K) = 0.
            PRC(K) = 0.
            NPRC(K) = 0.
            NPRC1(K) = 0.
            NSAGG(K) = 0.
            PSACWS(K) = 0.
            NPSACWS(K) = 0.
            PSACWI(K) = 0.
            NPSACWI(K) = 0.
            PRACS(K) = 0.
            NPRACS(K) = 0.
            NMULTS(K) = 0.
            QMULTS(K) = 0.
            NMULTR(K) = 0.
            QMULTR(K) = 0.
            NMULTG(K) = 0.
            QMULTG(K) = 0.
            NMULTRG(K) = 0.
            QMULTRG(K) = 0.
            MNUCCR(K) = 0.
            NNUCCR(K) = 0.
            PRA(K) = 0.
            NPRA(K) = 0.
            NRAGG(K) = 0.
            PRCI(K) = 0.
            NPRCI(K) = 0.
            PRAI(K) = 0.
            NPRAI(K) = 0.
            NNUCCD(K) = 0.
            MNUCCD(K) = 0.
            PCC(K) = 0.
            PRE(K) = 0.
            PRD(K) = 0.
            PRDS(K) = 0.
            EPRD(K) = 0.
            EPRDS(K) = 0.
            NSUBC(K) = 0.
            NSUBI(K) = 0.
            NSUBS(K) = 0.
            NSUBR(K) = 0.
            PIACR(K) = 0.
            NIACR(K) = 0.
            PRACI(K) = 0.
            PIACRS(K) = 0.
            NIACRS(K) = 0.
            PRACIS(K) = 0.
! HM: ADD GRAUPEL PROCESSES
            PRACG(K) = 0.
            PSACR(K) = 0.
	    PSACWG(K) = 0.
	    PGSACW(K) = 0.
            PGRACS(K) = 0.
	    PRDG(K) = 0.
	    EPRDG(K) = 0.
	    NPRACG(K) = 0.
	    NPSACWG(K) = 0.
	    NSCNG(K) = 0.
 	    NGRACS(K) = 0.
	    NSUBG(K) = 0.

            BERFD(K) = 0.0
            BERFDS(K) = 0.0
            BERFDG(K) = 0.0

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATION OF MICROPHYSICAL PROCESS RATES
! ACCRETION/AUTOCONVERSION/FREEZING/MELTING/COAG.
!.......................................................................
! FREEZING OF CLOUD DROPLETS
! ONLY ALLOWED BELOW -4 C
        IF (QC3D(K).GE.QSMALL .AND. T3D(K).LT.269.15) THEN

! NUMBER OF CONTACT NUCLEI (M^-3) FROM MEYERS ET AL., 1992
! FACTOR OF 1000 IS TO CONVERT FROM L^-1 TO M^-3

! MEYERS CURVE

           NACNT = EXP(-2.80+0.262*(TMELT-T3D(K)))*1000.

! COOPER CURVE
!        NACNT =  5.*EXP(0.304*(TMELT-T3D(K)))

! FLECTHER
!     NACNT = 0.01*EXP(0.6*(TMELT-T3D(K)))

! CONTACT FREEZING

! MEAN FREE PATH

            DUM = 7.37*T3D(K)/(288.*10.*PRES(K))/100.

! EFFECTIVE DIFFUSIVITY OF CONTACT NUCLEI
! BASED ON BROWNIAN DIFFUSION

            DAP(K) = CONS37*T3D(K)*(1.+DUM/RIN)/MU(K)
 
           MNUCCC(K) = CONS38*DAP(K)*NACNT*EXP(LOG(CDIST1(K))+   &
                   LOG(GAMMA(PGAM(K)+5.))-4.*LOG(LAMC(K)))
           NNUCCC(K) = 2.*PI*DAP(K)*NACNT*CDIST1(K)*           &
                    GAMMA(PGAM(K)+2.)/                         &
                    LAMC(K)

! IMMERSION FREEZING (BIGG 1953)

           MNUCCC(K) = MNUCCC(K)+CONS39*                   &
                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
                   EXP(AIMM*(TMELT-T3D(K)))

           NNUCCC(K) = NNUCCC(K)+                                  &
            CONS40*EXP(LOG(CDIST1(K))+LOG(GAMMA(PGAM(K)+4.))-3.*LOG(LAMC(K)))              &
                *EXP(AIMM*(TMELT-T3D(K)))

! PUT IN A CATCH HERE TO PREVENT DIVERGENCE BETWEEN NUMBER CONC. AND
! MIXING RATIO, SINCE STRICT CONSERVATION NOT CHECKED FOR NUMBER CONC

           NNUCCC(K) = MIN(NNUCCC(K),NC3D(K)/DT)

! PUT IN A CATCH HERE FOR MNUCC too +++mhwang
           MNUCCC(K) = MIN(MNUCCC(K),QC3D(K)/DT)

        END IF

! test the ratio of MNUCC and QC +++mhwang
!        IF (QC3D(K).gt.1.0e-5) then
!          IF((MNUCCC(K)*3600/QC3D(K)).gt.50) then 

!           write(955, *) 'MNUCC', K, T3D(K), MNUCCC(K)*3600/QC3D(K),  &
!                  MNUCCC(K)*3600*CFL3D_TEMP(K)/max(1.0e-12, (QI3D(K)*CFI3D_TEMP(K))),            &
!                  CONS39*EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
!                   EXP(AIMM*(TMELT-T3D(K)))*3600/QC3D(K)
!           write(955, *) 'MNUCC2', DAP(K), PGAM(K), NACNT, DUM, RIN, MU(K), LAMC(K), CDIST1(K)
!           write(955, *) 'Immersion freezing', QC3D(K), QI3D(K), CONS39*                   &
!                  EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K)))*             &
!                   EXP(AIMM*(TMELT-T3D(K))), & 
!                          EXP(LOG(CDIST1(K))+LOG(GAMMA(7.+PGAM(K)))-6.*LOG(LAMC(K))), &
!                          EXP(AIMM*(TMELT-T3D(K))) 
!           write(955, *) 'immersion terms', CONS39, LOG(CDIST1(K)), LOG(GAMMA(7.+PGAM(K))), -6.*LOG(LAMC(K))
!           write(955, *) 'immersion terms2', CDIST1(K), PGAM(K), LAMC(K), GAMMA(PGAM(K)+4.)/                        &
!             GAMMA(PGAM(K)+3.)/LAMC(K)/2.*1.E6, CFL3D_TEMP(K), CFI3D_TEMP(K), CLDMAX(K), NNUCCC(K)*DT*RHO(K)*1.0e-6
!          END IF
!        ENDIF

#ifdef CLUBB_CRM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For the case of clex9_oct14, we need to decrease the ice    !
! nucleation in order for the cloud to persist for realistic  !
! lengths. It is suggested to reduce by a factor of 100       !
! This coefficient can be changed in the subroutine           !
! init_microphys of the microphys_driver subroutine           !
!                                                             !
        NNUCCC(K)=NNUCCC(K)*NNUCCC_REDUCE_COEF
!                                                             ! 
! Change made by Marc Pilon on 11/16/11                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif /* CLUBB */




!.................................................................
!.......................................................................
! AUTOCONVERSION OF CLOUD LIQUID WATER TO RAIN
! FORMULA FROM BEHENG (1994)
! USING NUMERICAL SIMULATION OF STOCHASTIC COLLECTION EQUATION
! AND INITIAL CLOUD DROPLET SIZE DISTRIBUTION SPECIFIED
! AS A GAMMA DISTRIBUTION

! USE MINIMUM VALUE OF 1.E-6 TO PREVENT FLOATING POINT ERROR

         IF (QC3D(K).GE.1.E-6) THEN

! HM ADD 12/13/06, REPLACE WITH NEWER FORMULA
! FROM KHAIROUTDINOV AND KOGAN 2000, MWR

            IF (IRAIN.EQ.0) THEN
                PRC(K)=1350.*QC3D(K)**2.47*  &
           (NC3D(K)/1.e6*RHO(K))**(-1.79)

! note: nprc1 is change in Nr,
! nprc is change in Nc

        NPRC1(K) = PRC(K)/CONS29
        NPRC(K) = PRC(K)/(QC3D(K)/NC3D(K))

                NPRC(K) = MIN(NPRC(K),NC3D(K)/DT)

             ELSE IF (IRAIN.EQ.1) THEN

! v1.4
! replace with seifert and beheng

        dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
        dum1 = 600.*dum**0.68*(1.-dum**0.68)**3

        prc(k) = 9.44e9/(20.*2.6e-7)* &
        (nu(k)+2.)*(nu(k)+4.)/(nu(k)+1.)**2* &
        (rho(k)*qc3d(k)/1000.)**4/(rho(k)*nc3d(k)/1.e6)**2* &
        (1.+dum1/(1.-dum)**2)*1000./rho(k)

        nprc(k) = prc(k)*2./2.6e-7*1000.
        nprc1(k) = 0.5*nprc(k)

         END IF
         END IF

!.......................................................................
! SELF-COLLECTION OF DROPLET NOT INCLUDED IN KK2000 SCHEME

! SNOW AGGREGATION FROM PASSARELLI, 1978, USED BY REISNER, 1998
! THIS IS HARD-WIRED FOR BS = 0.4 FOR NOW

         IF (QNI3D(K).GE.1.E-8) THEN
             NSAGG(K) = CONS15*ASN(K)*RHO(K)**            &
            ((2.+BS)/3.)*QNI3D(K)**((2.+BS)/3.)*                  &
            (NS3D(K)*RHO(K))**((4.-BS)/3.)/                       &
            (RHO(K))
         END IF

!.......................................................................
! ACCRETION OF CLOUD DROPLETS ONTO SNOW/GRAUPEL
! HERE USE CONTINUOUS COLLECTION EQUATION WITH
! SIMPLE GRAVITATIONAL COLLECTION KERNEL IGNORING

! SNOW

         IF (QNI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

           PSACWS(K) = CONS13*ASN(K)*QC3D(K)*RHO(K)*               &
                  N0S(K)/                        &
                  LAMS(K)**(BS+3.)
           NPSACWS(K) = CONS13*ASN(K)*NC3D(K)*RHO(K)*              &
                  N0S(K)/                        &
                  LAMS(K)**(BS+3.)

         END IF

!............................................................................
! COLLECTION OF CLOUD WATER BY GRAUPEL

         IF (QG3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

           PSACWG(K) = CONS14*AGN(K)*QC3D(K)*RHO(K)*               &
                  N0G(K)/                        &
                  LAMG(K)**(BG+3.)
           NPSACWG(K) = CONS14*AGN(K)*NC3D(K)*RHO(K)*              &
                  N0G(K)/                        &
                  LAMG(K)**(BG+3.)
	    END IF

!.......................................................................
! HM, ADD 12/13/06
! CLOUD ICE COLLECTING DROPLETS, ASSUME THAT CLOUD ICE MEAN DIAM > 100 MICRON
! BEFORE RIMING CAN OCCUR
! ASSUME THAT RIME COLLECTED ON CLOUD ICE DOES NOT LEAD
! TO HALLET-MOSSOP SPLINTERING

         IF (QI3D(K).GE.1.E-8 .AND. QC3D(K).GE.QSMALL) THEN

! PUT IN SIZE DEPENDENT COLLECTION EFFICIENCY BASED ON STOKES LAW
! FROM THOMPSON ET AL. 2004, MWR

            IF (1./LAMI(K).GE.100.E-6) THEN

           PSACWI(K) = CONS16*AIN(K)*QC3D(K)*RHO(K)*               &
                  N0I(K)/                        &
                  LAMI(K)**(BI+3.)
           NPSACWI(K) = CONS16*AIN(K)*NC3D(K)*RHO(K)*              &
                  N0I(K)/                        &
                  LAMI(K)**(BI+3.)
           END IF
         END IF

!.......................................................................
! ACCRETION OF RAIN WATER BY SNOW
! FORMULA FROM IKAWA AND SAITO, 1991, USED BY REISNER ET AL, 1998

         IF (QR3D(K).GE.1.E-8.AND.QNI3D(K).GE.1.E-8) THEN

            UMS = ASN(K)*CONS3/(LAMS(K)**BS)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNS = ASN(K)*CONS5/LAMS(K)**BS
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMS=MIN(UMS,1.2*dum)
            UNS=MIN(UNS,1.2*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACS(K) = CONS41*(((1.2*UMR-0.95*UMS)**2+                   &
                  0.08*UMS*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0S(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMS(K))+                    &
                  2./(LAMR(K)**2*LAMS(K)**2)+                  &				 
                  0.5/(LAMR(k)*LAMS(k)**3)))

            NPRACS(K) = CONS32*RHO(K)*(1.7*(UNR-UNS)**2+            &
                0.3*UNR*UNS)**0.5*N0RR(K)*N0S(K)*              &
                (1./(LAMR(K)**3*LAMS(K))+                      &
                 1./(LAMR(K)**2*LAMS(K)**2)+                   &
                 1./(LAMR(K)*LAMS(K)**3))

! MAKE SURE PRACS DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACS(K) = MIN(PRACS(K),QR3D(K)/DT)

! COLLECTION OF SNOW BY RAIN - NEEDED FOR GRAUPEL CONVERSION CALCULATIONS
! ONLY CALCULATE IF SNOW AND RAIN MIXING RATIOS EXCEED 0.1 G/KG

! V1.3
! ASSUME COLLECTION OF SNOW BY RAIN PRODUCES GRAUPEL NOT HAIL

! V1.5
!            IF (IHAIL.EQ.0) THEN
            IF (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) THEN
            PSACR(K) = CONS31*(((1.2*UMR-0.95*UMS)**2+              &
                  0.08*UMS*UMR)**0.5*RHO(K)*                     &
                 N0RR(K)*N0S(K)/LAMS(K)**3*                               &
                  (5./(LAMS(K)**3*LAMR(K))+                    &
                  2./(LAMS(K)**2*LAMR(K)**2)+                  &
                  0.5/(LAMS(K)*LAMR(K)**3)))            
            END IF
!            END IF

         END IF

!.......................................................................

! COLLECTION OF RAINWATER BY GRAUPEL, FROM IKAWA AND SAITO 1990, 
! USED BY REISNER ET AL 1998
         IF (QR3D(K).GE.1.E-8.AND.QG3D(K).GE.1.E-8) THEN

            UMG = AGN(K)*CONS7/(LAMG(K)**BG)
            UMR = ARN(K)*CONS4/(LAMR(K)**BR)
            UNG = AGN(K)*CONS8/LAMG(K)**BG
            UNR = ARN(K)*CONS6/LAMR(K)**BR

! SET REASLISTIC LIMITS ON FALLSPEEDS
! bug fix, 10/08/09
            dum=(rhosu/rho(k))**0.54
            UMG=MIN(UMG,20.*dum)
            UNG=MIN(UNG,20.*dum)
            UMR=MIN(UMR,9.1*dum)
            UNR=MIN(UNR,9.1*dum)

            PRACG(K) = CONS41*(((1.2*UMR-0.95*UMG)**2+                   &
                  0.08*UMG*UMR)**0.5*RHO(K)*                      &
                  N0RR(K)*N0G(K)/LAMR(K)**3*                              &
                  (5./(LAMR(K)**3*LAMG(K))+                    &
                  2./(LAMR(K)**2*LAMG(K)**2)+				   &
				  0.5/(LAMR(k)*LAMG(k)**3)))

            NPRACG(K) = CONS32*RHO(K)*(1.7*(UNR-UNG)**2+            &
                0.3*UNR*UNG)**0.5*N0RR(K)*N0G(K)*              &
                (1./(LAMR(K)**3*LAMG(K))+                      &
                 1./(LAMR(K)**2*LAMG(K)**2)+                   &
                 1./(LAMR(K)*LAMG(K)**3))

! MAKE SURE PRACG DOESN'T EXCEED TOTAL RAIN MIXING RATIO
! AS THIS MAY OTHERWISE RESULT IN TOO MUCH TRANSFER OF WATER DURING
! RIME-SPLINTERING

            PRACG(K) = MIN(PRACG(K),QR3D(K)/DT)

	    END IF

!.......................................................................
! RIME-SPLINTERING - SNOW
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW AND DROPLET MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS
! THESE THRESHOLDS CORRESPOND WITH GRAUPEL THRESHOLDS IN RH 1984

!v1.4
         IF (QNI3D(K).GE.0.1E-3) THEN
         IF (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) THEN
         IF (PSACWS(K).GT.0..OR.PRACS(K).GT.0.) THEN
            IF (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) THEN

               IF (T3D(K).GT.270.16) THEN
                  FMULT = 0.
               ELSE IF (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  THEN
                  FMULT = (270.16-T3D(K))/2.
               ELSE IF (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   THEN
                  FMULT = (T3D(K)-265.16)/3.
               ELSE IF (T3D(K).LT.265.16) THEN
                  FMULT = 0.
               END IF

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO SNOW

               IF (PSACWS(K).GT.0.) THEN
                  NMULTS(K) = 35.E4*PSACWS(K)*FMULT*1000.
                  QMULTS(K) = NMULTS(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                  QMULTS(K) = MIN(QMULTS(K),PSACWS(K))
                  PSACWS(K) = PSACWS(K)-QMULTS(K)

               END IF

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               IF (PRACS(K).GT.0.) THEN
                   NMULTR(K) = 35.E4*PRACS(K)*FMULT*1000.
                   QMULTR(K) = NMULTR(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM SNOW TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO SNOW

                   QMULTR(K) = MIN(QMULTR(K),PRACS(K))

                   PRACS(K) = PRACS(K)-QMULTR(K)

               END IF

            END IF
            END IF
         END IF
         END IF

!.......................................................................
! RIME-SPLINTERING - GRAUPEL 
! HALLET-MOSSOP (1974)
! NUMBER OF SPLINTERS FORMED IS BASED ON MASS OF RIMED WATER

! DUM1 = MASS OF INDIVIDUAL SPLINTERS

! HM ADD THRESHOLD SNOW MIXING RATIO FOR RIME-SPLINTERING
! TO LIMIT RIME-SPLINTERING IN STRATIFORM CLOUDS

! V1.3
! ONLY CALCULATE FOR GRAUPEL NOT HAIL
! V1.5
!         IF (IHAIL.EQ.0) THEN
! v1.4
         IF (QG3D(K).GE.0.1E-3) THEN
         IF (QC3D(K).GE.0.5E-3.OR.QR3D(K).GE.0.1E-3) THEN
         IF (PSACWG(K).GT.0..OR.PRACG(K).GT.0.) THEN
            IF (T3D(K).LT.270.16 .AND. T3D(K).GT.265.16) THEN

               IF (T3D(K).GT.270.16) THEN
                  FMULT = 0.
               ELSE IF (T3D(K).LE.270.16.AND.T3D(K).GT.268.16)  THEN
                  FMULT = (270.16-T3D(K))/2.
               ELSE IF (T3D(K).GE.265.16.AND.T3D(K).LE.268.16)   THEN
                  FMULT = (T3D(K)-265.16)/3.
               ELSE IF (T3D(K).LT.265.16) THEN
                  FMULT = 0.
               END IF

! 1000 IS TO CONVERT FROM KG TO G

! SPLINTERING FROM DROPLETS ACCRETED ONTO GRAUPEL

               IF (PSACWG(K).GT.0.) THEN
                  NMULTG(K) = 35.E4*PSACWG(K)*FMULT*1000.
                  QMULTG(K) = NMULTG(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                  QMULTG(K) = MIN(QMULTG(K),PSACWG(K))
                  PSACWG(K) = PSACWG(K)-QMULTG(K)

               END IF

! RIMING AND SPLINTERING FROM ACCRETED RAINDROPS

               IF (PRACG(K).GT.0.) THEN
                   NMULTRG(K) = 35.E4*PRACG(K)*FMULT*1000.
                   QMULTRG(K) = NMULTRG(K)*MMULT

! CONSTRAIN SO THAT TRANSFER OF MASS FROM GRAUPEL TO ICE CANNOT BE MORE MASS
! THAN WAS RIMED ONTO GRAUPEL

                   QMULTRG(K) = MIN(QMULTRG(K),PRACG(K))
                   PRACG(K) = PRACG(K)-QMULTRG(K)

               END IF

            END IF
            END IF
            END IF
            END IF
!         END IF

!........................................................................
! CONVERSION OF RIMED CLOUD WATER ONTO SNOW TO GRAUPEL
! ASSUME CONVERTED SNOW FORMS GRAUPEL NOT HAIL
! HAIL ASSUMED TO ONLY FORM BY FREEZING OF RAIN
! OR COLLISIONS OF RAIN WITH CLOUD ICE

! V1.3
! V1.5
!           IF (IHAIL.EQ.0) THEN
	   IF (PSACWS(K).GT.0.) THEN
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QC > 0.5 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              IF (QNI3D(K).GE.0.1E-3.AND.QC3D(K).GE.0.5E-3) THEN

! PORTION OF RIMING CONVERTED TO GRAUPEL (REISNER ET AL. 1998, ORIGINALLY IS1991)
	     PGSACW(K) = MIN(PSACWS(K),CONS17*DT*N0S(K)*QC3D(K)*QC3D(K)* &
                          ASN(K)*ASN(K)/ &
                           (RHO(K)*LAMS(K)**(2.*BS+2.))) 

! MIX RAT CONVERTED INTO GRAUPEL AS EMBRYO (REISNER ET AL. 1998, ORIG M1990)
	     DUM = MAX(RHOSN/(RHOG-RHOSN)*PGSACW(K),0.) 

! NUMBER CONCENTRAITON OF EMBRYO GRAUPEL FROM RIMING OF SNOW
	     NSCNG(K) = DUM/MG0*RHO(K)
! LIMIT MAX NUMBER CONVERTED TO SNOW NUMBER
             NSCNG(K) = MIN(NSCNG(K),NS3D(K)/DT)

! PORTION OF RIMING LEFT FOR SNOW
             PSACWS(K) = PSACWS(K) - PGSACW(K)
             END IF
	   END IF

! CONVERSION OF RIMED RAINWATER ONTO SNOW CONVERTED TO GRAUPEL

	   IF (PRACS(K).GT.0.) THEN
! ONLY ALLOW CONVERSION IF QNI > 0.1 AND QR > 0.1 G/KG FOLLOWING RUTLEDGE AND HOBBS (1984)
              IF (QNI3D(K).GE.0.1E-3.AND.QR3D(K).GE.0.1E-3) THEN
! PORTION OF COLLECTED RAINWATER CONVERTED TO GRAUPEL (REISNER ET AL. 1998)
	      DUM = CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3 &    
                   /(CONS18*(4./LAMS(K))**3*(4./LAMS(K))**3+ &  
                   CONS19*(4./LAMR(K))**3*(4./LAMR(K))**3)
              DUM=MIN(DUM,1.)
              DUM=MAX(DUM,0.)
	      PGRACS(K) = (1.-DUM)*PRACS(K)
            NGRACS(K) = (1.-DUM)*NPRACS(K)
! LIMIT MAX NUMBER CONVERTED TO MIN OF EITHER RAIN OR SNOW NUMBER CONCENTRATION
            NGRACS(K) = MIN(NGRACS(K),NR3D(K)/DT)
            NGRACS(K) = MIN(NGRACS(K),NS3D(K)/DT)

! AMOUNT LEFT FOR SNOW PRODUCTION
            PRACS(K) = PRACS(K) - PGRACS(K)
            NPRACS(K) = NPRACS(K) - NGRACS(K)
! CONVERSION TO GRAUPEL DUE TO COLLECTION OF SNOW BY RAIN
            PSACR(K)=PSACR(K)*(1.-DUM)
            END IF
	   END IF
!           END IF

!.......................................................................
! FREEZING OF RAIN DROPS
! FREEZING ALLOWED BELOW -4 C

         IF (T3D(K).LT.269.15.AND.QR3D(K).GE.QSMALL) THEN

! IMMERSION FREEZING (BIGG 1953)
            MNUCCR(K) = CONS20*NR3D(K)*EXP(AIMM*(TMELT-T3D(K)))/LAMR(K)**3 &
                 /LAMR(K)**3

            NNUCCR(K) = PI*NR3D(K)*BIMM*EXP(AIMM*(TMELT-T3D(K)))/LAMR(K)**3

! PREVENT DIVERGENCE BETWEEN MIXING RATIO AND NUMBER CONC
            NNUCCR(K) = MIN(NNUCCR(K),NR3D(K)/DT)

! set limit on MNUCCR +++mhwang
            MNUCCR(K) = MIN(MNUCCR(K), QR3D(K)/DT)

         END IF

!.......................................................................
! ACCRETION OF CLOUD LIQUID WATER BY RAIN
! CONTINUOUS COLLECTION EQUATION WITH
! GRAVITATIONAL COLLECTION KERNEL, DROPLET FALL SPEED NEGLECTED

         IF (QR3D(K).GE.1.E-8 .AND. QC3D(K).GE.1.E-8) THEN

! 12/13/06 HM ADD, REPLACE WITH NEWER FORMULA FROM
! KHAIROUTDINOV AND KOGAN 2000, MWR

            IF (IRAIN.EQ.0) THEN

           DUM=(QC3D(K)*QR3D(K))
           PRA(K) = 67.*(DUM)**1.15
           NPRA(K) = PRA(K)/(QC3D(K)/NC3D(K))

           ELSE IF (IRAIN.EQ.1) THEN

! v1.4
! seifert and beheng (2001) formulation

           dum = 1.-qc3d(k)/(qc3d(k)+qr3d(k))
           dum1 = (dum/(dum+5.e-4))**4
           pra(k) = 5.78e3*rho(k)/1000.*qc3d(k)*qr3d(k)*dum1
           npra(k) = pra(k)*rho(k)/1000.*(nc3d(k)*rho(k)/1.e6)/ &
           (qc3d(k)*rho(k)/1000.)*1.e6/rho(k)

         END IF
         END IF
!.......................................................................
! SELF-COLLECTION OF RAIN DROPS
! FROM BEHENG(1994)
! FROM NUMERICAL SIMULATION OF THE STOCHASTIC COLLECTION EQUATION
! AS DESCRINED ABOVE FOR AUTOCONVERSION

! v1.4 replace with seifert and beheng (2001)

         IF (QR3D(K).GE.1.E-8) THEN
! include breakup add 10/09/09
            dum1=300.e-6
            if (1./lamr(k).lt.dum1) then
            dum=1.
            else if (1./lamr(k).ge.dum1) then
            dum=2.-exp(2300.*(1./lamr(k)-dum1))
            end if
!            NRAGG(K) = -8.*NR3D(K)*QR3D(K)*RHO(K)
            NRAGG(K) = -5.78*dum*NR3D(K)*QR3D(K)*RHO(K)
!            NRAGG(K) = 0.0             
         END IF

!.......................................................................
! AUTOCONVERSION OF CLOUD ICE TO SNOW
! FOLLOWING HARRINGTON ET AL. (1995) WITH MODIFICATION
! HERE IT IS ASSUMED THAT AUTOCONVERSION CAN ONLY OCCUR WHEN THE
! ICE IS GROWING, I.E. IN CONDITIONS OF ICE SUPERSATURATION

#ifndef CLDFRAC
         IF (QI3D(K).GE.1.E-8 .AND.QVQVSI(K).GE.1.) THEN

!           COFFI = 2./LAMI(K)
!           IF (COFFI.GE.DCS) THEN
              NPRCI(K) = CONS21*(QV3D(K)-QVI(K))*RHO(K)                         &
                *N0I(K)*EXP(-LAMI(K)*DCS)*DV(K)/ABI(K)
              PRCI(K) = CONS22*NPRCI(K)
              NPRCI(K) = MIN(NPRCI(K),NI3D(K)/DT)

!           END IF
         END IF
#else 
         IF (QI3D(K).GE.1.E-8) THEN
!           IF(QC3D(K).GE.QSMALL) then
           IF(CFL3D_TEMP(K).GT.(1.01*cloud_frac_thresh)) THEN 
! For mixed-phase clouds, QV_ICLD is set to the hybrid 
! saturaiton of liquid water and ice water. 
!
!              QV_ICLD(K)=(QC3D(K)*QVS(K)+QI3D(K)*QVI(K))/(QC3D(K)+QI3D(K))
              QV_ICLD(K)=(QC3D(K)*QVS(K)*CFL3D_TEMP(K)+QI3D(K)*QVI(K)*CFI3D_TEMP(K))  &
                         /(QC3D(K)*CFL3D_TEMP(K)+QI3D(K)*CFI3D_TEMP(K))
              QV_ICLD(K)= max(QV3D(K), QV_ICLD(K))
           ELSE
!              QV_ICLD(K) = QV3D(K) 
              QV_ICLD(K) = QV3D(K) * 1.0/0.9 ! ice clouds are allowed when RHI=0.9
           ENDIF
!           QV_ICLD(K)= QVS(K)  ! test +++mhwang

           IF(QV_ICLD(K).GT.QVI(K)) THEN
!           COFFI = 2./LAMI(K)
!           IF (COFFI.GE.DCS) THEN
              NPRCI(K) = CONS21*(QV_ICLD(K)-QVI(K))*RHO(K)                         &
                *N0I(K)*EXP(-LAMI(K)*DCS)*DV(K)/ABI(K)

! New formula from CAM5.2: 
! note: assumes autoconversion timescale of 180 sec
!              NPRCI(K) = N0I(K)/(LAMI(K)*180.)*exp(-LAMI(K)*DCS)
! test mhwang  reducing the autoconversoin rate of ice by a factor of 10.
!              NPRCI(K) = 0.1 * NPRCI(K)


              PRCI(K) = CONS22*NPRCI(K)
              NPRCI(K) = MIN(NPRCI(K),NI3D(K)/DT)

!           END IF
           END IF
         END IF
#endif

!.......................................................................
! ACCRETION OF CLOUD ICE BY SNOW
! FOR THIS CALCULATION, IT IS ASSUMED THAT THE VS >> VI
! AND DS >> DI FOR CONTINUOUS COLLECTION

         IF (QNI3D(K).GE.1.E-8 .AND. QI3D(K).GE.QSMALL) THEN
            PRAI(K) = CONS23*ASN(K)*QI3D(K)*RHO(K)*N0S(K)/     &
                     LAMS(K)**(BS+3.)
            NPRAI(K) = CONS23*ASN(K)*NI3D(K)*                                       &
                  RHO(K)*N0S(K)/                                 &
                  LAMS(K)**(BS+3.)
            NPRAI(K)=MIN(NPRAI(K),NI3D(K)/DT)
         END IF

!.......................................................................
! HM, ADD 12/13/06, COLLISION OF RAIN AND ICE TO PRODUCE SNOW OR GRAUPEL
! FOLLOWS REISNER ET AL. 1998
! ASSUMED FALLSPEED AND SIZE OF ICE CRYSTAL << THAN FOR RAIN

         IF (QR3D(K).GE.1.E-8.AND.QI3D(K).GE.1.E-8.AND.T3D(K).LE.TMELT) THEN

! ALLOW GRAUPEL FORMATION FROM RAIN-ICE COLLISIONS ONLY IF RAIN MIXING RATIO > 0.1 G/KG,
! OTHERWISE ADD TO SNOW

            IF (QR3D(K).GE.0.1E-3) THEN
            NIACR(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)*RHO(K)
            PIACR(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
            PRACI(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                LAMR(K)**(BR+3.)*RHO(K)
            NIACR(K)=MIN(NIACR(K),NR3D(K)/DT)
            NIACR(K)=MIN(NIACR(K),NI3D(K)/DT)
            ELSE 
            NIACRS(K)=CONS24*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)*RHO(K)
            PIACRS(K)=CONS25*NI3D(K)*N0RR(K)*ARN(K) &
                /LAMR(K)**(BR+3.)/LAMR(K)**3*RHO(K)
            PRACIS(K)=CONS24*QI3D(K)*N0RR(K)*ARN(K)/ &
                LAMR(K)**(BR+3.)*RHO(K)
            NIACRS(K)=MIN(NIACRS(K),NR3D(K)/DT)
            NIACRS(K)=MIN(NIACRS(K),NI3D(K)/DT)
            END IF
         END IF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 101      CONTINUE

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE EVAP OF QR

      IF (QR3D(K).GE.QSMALL) THEN
        EPSR = 2.*PI*N0RR(K)*RHO(K)*DV(K)*                           &
                   (F1R/(LAMR(K)*LAMR(K))+                       &
                    F2R*(ARN(K)*RHO(K)/MU(K))**0.5*                      &
                    SC(K)**(1./3.)*CONS9/                   &
                (LAMR(K)**CONS34))
      ELSE
      EPSR = 0.
      END IF

! NO CONDENSATION ONTO RAIN, ONLY EVAP

#ifndef CLDFRAC
           IF (QV3D(K).LT.QVS(K)) THEN
              PRE(K) = EPSR*(QV3D(K)-QVS(K))/AB(K)
              PRE(K) = MIN(PRE(K),0.)
           ELSE
              PRE(K) = 0.
           END IF
#else
!         calcluate QV3D for out of cloud region, assuming water saturaiton in-cloud
           QVCLR3D(K) = (QV3D(K)-QVS(K)*CFL3D_TEMP(K))/(1.-CFL3D_TEMP(K))
!         only calculate if there is some precipitation fraction > cloud fraction
           IF (QVCLR3D(K).LT.QVS(K)) THEN
              PRE(K) = EPSR*(QVCLR3D(K)-QVS(K))/AB(K)
              PRE(K) = MIN(PRE(K)*(CLDMAX(K)-CFL3D_TEMP(K)),0.)
              PRE(K) = PRE(K)/CLDMAX(K)
           ELSE
              PRE(K) = 0.
           END IF
! test +++mhwang
!           IF(QC3D(K).GT.QSMALL) PRE(K) = 0.0   ! When QC is present, no rain evaporation

#endif

!.......................................................................
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

! CONSERVATION OF WATER
! THIS IS ADOPTED LOOSELY FROM MM5 RESINER CODE. HOWEVER, HERE WE
! ONLY ADJUST PROCESSES THAT ARE NEGATIVE, RATHER THAN ALL PROCESSES.
! THIS SECTION IS SEPARATED INTO TWO PARTS, IF T < 0 C, T > 0 C
! DUE TO DIFFERENT PROCESSES THAT ACT DEPENDING ON FREEZING/ABOVE FREEZING

! IF MIXING RATIOS LESS THAN QSMALL, THEN NO DEPLETION OF WATER
! THROUGH MICROPHYSICAL PROCESSES, SKIP CONSERVATION

! NOTE: CONSERVATION CHECK NOT APPLIED TO NUMBER CONCENTRATION SPECIES. ADDITIONAL CATCH
! BELOW WILL PREVENT NEGATIVE NUMBER CONCENTRATION
! FOR EACH MICROPHYSICAL PROCESS WHICH PROVIDES A SOURCE FOR NUMBER, THERE IS A CHECK
! TO MAKE SURE THAT CAN'T EXCEED TOTAL NUMBER OF DEPLETED SPECIES WITH THE TIME
! STEP

!****SENSITIVITY - NO ICE

      IF (ILIQ.EQ.1) THEN
      MNUCCC(K)=0.
      NNUCCC(K)=0.
      MNUCCR(K)=0.
      NNUCCR(K)=0.
      MNUCCD(K)=0.
      NNUCCD(K)=0.
      END IF

! ****SENSITIVITY - NO GRAUPEL
      IF (IGRAUP.EQ.1) THEN
            PRACG(K) = 0.
            PSACR(K) = 0.
	    PSACWG(K) = 0.
	    PGSACW(K) = 0.
            PGRACS(K) = 0.
	    PRDG(K) = 0.
	    EPRDG(K) = 0.
            EVPMG(K) = 0.
            PGMLT(K) = 0.
	    NPRACG(K) = 0.
	    NPSACWG(K) = 0.
	    NSCNG(K) = 0.
 	    NGRACS(K) = 0.
	    NSUBG(K) = 0.
	    NGMLTG(K) = 0.
            NGMLTR(K) = 0.
! fix 053011
            PIACRS(K)=PIACRS(K)+PIACR(K)
            PIACR(K) = 0.
#ifdef CLDFRAC
            BERFDG(K) = 0.
#endif
       END IF

! CONSERVATION OF QC

#ifndef CLDFRAC
      DUM = (PRC(K)+PRA(K)+MNUCCC(K)+PSACWS(K)+PSACWI(K)+QMULTS(K)+PSACWG(K)+PGSACW(K)+QMULTG(K))*DT
#else
      DUM = (PRC(K)+PRA(K)+MNUCCC(K)+PSACWS(K)+PSACWI(K)+QMULTS(K)+PSACWG(K)+PGSACW(K)+QMULTG(K)+BERFD(K)+BERFDS(K)+BERFDG(K))*DT
#endif

      IF (DUM.GT.QC3D(K).AND.QC3D(K).GE.QSMALL) THEN
        RATIO = QC3D(K)/DUM

        PRC(K) = PRC(K)*RATIO
        PRA(K) = PRA(K)*RATIO
        MNUCCC(K) = MNUCCC(K)*RATIO
        PSACWS(K) = PSACWS(K)*RATIO
        PSACWI(K) = PSACWI(K)*RATIO
        QMULTS(K) = QMULTS(K)*RATIO
        QMULTG(K) = QMULTG(K)*RATIO
        PSACWG(K) = PSACWG(K)*RATIO
	PGSACW(K) = PGSACW(K)*RATIO
#ifdef CLDFRAC
        BERFD(K) = BERFD(K)*RATIO 
        BERFDS(K) = BERFDS(K)*RATIO
        BERFDG(K) = BERFDG(K)*RATIO
#endif
        END IF
 
! CONSERVATION OF QI

#ifndef CLDFRAC
      DUM = (-PRD(K)-MNUCCC(K)+PRCI(K)+PRAI(K)-QMULTS(K)-QMULTG(K)-QMULTR(K)-QMULTRG(K) &
                -MNUCCD(K)+PRACI(K)+PRACIS(K)-EPRD(K)-PSACWI(K))*DT
#else
      DUM = (-BERFD(K)*CFL3D_TEMP(K)-PRD(K)*CFI3D_TEMP(K)-MNUCCC(K)*CFL3D_TEMP(K)+(PRCI(K)+PRAI(K))*CFI3D_TEMP(K)-(QMULTS(K)+QMULTG(K))*CFL3D_TEMP(K) & 
             -(QMULTR(K)+QMULTRG(K))*CLDMAX(K) &
               +(-MNUCCD(K)+PRACI(K)+PRACIS(K)-EPRD(K))*CFI3D_TEMP(K)-PSACWI(K)*CFL3D_TEMP(K))*DT
#endif

#ifndef CLDFRAC
      IF (DUM.GT.QI3D(K).AND.QI3D(K).GE.QSMALL) THEN

        RATIO = (QI3D(K)/DT+PRD(K)+MNUCCC(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+ &
                     MNUCCD(K)+PSACWI(K))/ &
                      (PRCI(K)+PRAI(K)+PRACI(K)+PRACIS(K)-EPRD(K))
#else
      IF (DUM.GT.(QI3D(K)*CFI3D_TEMP(K)).AND.QI3D(K).GE.QSMALL) THEN

        RATIO = (QI3D(K)*CFI3D_TEMP(K)/DT+PRD(K)*CFI3D_TEMP(K)+(BERFD(K)+MNUCCC(K))*CFL3D_TEMP(K)+(QMULTS(K)+QMULTG(K))*CFL3D_TEMP(K) &
                     +(QMULTR(K)+QMULTRG(K))*CLDMAX(K)+ &
                     MNUCCD(K)*CFI3D_TEMP(K)+PSACWI(K)*CFL3D_TEMP(K))/ &
                      ((PRCI(K)+PRAI(K)+PRACI(K)+PRACIS(K)-EPRD(K))*CFI3D_TEMP(K))
#endif

        PRCI(K) = PRCI(K)*RATIO
        PRAI(K) = PRAI(K)*RATIO
        PRACI(K) = PRACI(K)*RATIO
        PRACIS(K) = PRACIS(K)*RATIO
        EPRD(K) = EPRD(K)*RATIO

        END IF

! CONSERVATION OF QR

#ifndef CLDFRAC
      DUM=((PRACS(K)-PRE(K))+(QMULTR(K)+QMULTRG(K)-PRC(K))+(MNUCCR(K)-PRA(K))+ &
             PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))*DT

      IF (DUM.GT.QR3D(K).AND.QR3D(K).GE.QSMALL) THEN

        RATIO = (QR3D(K)/DT+PRC(K)+PRA(K))/ &
             (-PRE(K)+QMULTR(K)+QMULTRG(K)+PRACS(K)+MNUCCR(K)+PIACR(K)+PIACRS(K)+PGRACS(K)+PRACG(K))
#else
      DUM=((PRACS(K)-PRE(K))*CLDMAX(K)+((QMULTR(K)+QMULTRG(K))*CLDMAX(K)-PRC(K)*CFL3D_TEMP(K))+(MNUCCR(K)*CLDMAX(K)-PRA(K)*CFL3D_TEMP(K))+ &
             (PIACR(K)+PIACRS(K))*CFI3D_TEMP(K)+(PGRACS(K)+PRACG(K))*CLDMAX(K))*DT

      IF (DUM.GT.(QR3D(K)*CLDMAX(K)).AND.QR3D(K).GE.QSMALL) THEN

        RATIO = (QR3D(K)*CLDMAX(K)/DT+(PRC(K)+PRA(K))*CFL3D_TEMP(K))/ &
             ((-PRE(K)+QMULTR(K)+QMULTRG(K)+PRACS(K)+MNUCCR(K))*CLDMAX(K)+(PIACR(K)+PIACRS(K))*CFI3D_TEMP(K)+(PGRACS(K)+PRACG(K))*CLDMAX(K))
#endif

        PRE(K) = PRE(K)*RATIO
        PRACS(K) = PRACS(K)*RATIO
        QMULTR(K) = QMULTR(K)*RATIO
        QMULTRG(K) = QMULTRG(K)*RATIO
        MNUCCR(K) = MNUCCR(K)*RATIO
        PIACR(K) = PIACR(K)*RATIO
        PIACRS(K) = PIACRS(K)*RATIO
        PGRACS(K) = PGRACS(K)*RATIO
        PRACG(K) = PRACG(K)*RATIO

        END IF

! CONSERVATION OF QNI
! CONSERVATION FOR GRAUPEL SCHEME

        IF (IGRAUP.EQ.0) THEN

#ifndef CLDFRAC
      DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K))*DT

      IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

        RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K))/(-EPRDS(K)+PSACR(K))
#else
      DUM = (-BERFDS(K)*CFL3D_TEMP(K)-PRDS(K)*CLDMAX(K)+(-PSACWS(K)-PRAI(K)-PRCI(K))*CFI3D_TEMP(K)+  &
            (-PRACS(K)-EPRDS(K)+PSACR(K))*CLDMAX(K)-(PIACRS(K)+PRACIS(K))*CFI3D_TEMP(K))*DT

      IF (DUM.GT.(QNI3D(K)*CLDMAX(K)).AND.QNI3D(K).GE.QSMALL) THEN

        RATIO = (QNI3D(K)*CLDMAX(K)/DT+BERFDS(K)*CF3D_TEMP(K)+PRDS(K)*CLDMAX(K)+(PSACWS(K)+PRAI(K)+PRCI(K))*CFI3D_TEMP(K)+  &
         PRACS(K)*CLDMAX(K)+(PIACRS(K)+PRACIS(K))*CFI3D_TEMP(K))/((-EPRDS(K)+PSACR(K))*CLDMAX(K))
#endif

       EPRDS(K) = EPRDS(K)*RATIO
       PSACR(K) = PSACR(K)*RATIO

       END IF

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
       ELSE IF (IGRAUP.EQ.1) THEN

#ifndef CLDFRAC
      DUM = (-PRDS(K)-PSACWS(K)-PRAI(K)-PRCI(K)-PRACS(K)-EPRDS(K)+PSACR(K)-PIACRS(K)-PRACIS(K)-MNUCCR(K))*DT

      IF (DUM.GT.QNI3D(K).AND.QNI3D(K).GE.QSMALL) THEN

       RATIO = (QNI3D(K)/DT+PRDS(K)+PSACWS(K)+PRAI(K)+PRCI(K)+PRACS(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))/(-EPRDS(K)+PSACR(K))
#else
      DUM = (-BERFDS(K)*CFL3D_TEMP(K)-PRDS(K)*CLDMAX(K)+(-PSACWS(K)-PRAI(K)-PRCI(K))*CFI3D_TEMP(K)+  &
            (-PRACS(K)-EPRDS(K)+PSACR(K))*CLDMAX(K)-(PIACRS(K)+PRACIS(K))*CFI3D_TEMP(K)-MNUCCR(K)*CLDMAX(K))*DT

      IF (DUM.GT.(QNI3D(K)*CLDMAX(K)).AND.QNI3D(K).GE.QSMALL) THEN

        RATIO = (QNI3D(K)*CLDMAX(K)/DT+BERFDS(K)*CF3D_TEMP(K)+PRDS(K)*CLDMAX(K)+(PSACWS(K)+PRAI(K)+PRCI(K))*CFI3D_TEMP(K)+  &         
         PRACS(K)*CLDMAX(K)+(PIACRS(K)+PRACIS(K))*CFI3D_TEMP(K)+MNUCCR(K)*CLDMAX(K))/((-EPRDS(K)+PSACR(K))*CLDMAX(K))
#endif

       EPRDS(K) = EPRDS(K)*RATIO
       PSACR(K) = PSACR(K)*RATIO

       END IF

       END IF

! CONSERVATION OF QG

#ifndef CLDFRAC
      DUM = (-PSACWG(K)-PRACG(K)-PGSACW(K)-PGRACS(K)-PRDG(K)-MNUCCR(K)-EPRDG(K)-PIACR(K)-PRACI(K)-PSACR(K))*DT

      IF (DUM.GT.QG3D(K).AND.QG3D(K).GE.QSMALL) THEN

        RATIO = (QG3D(K)/DT+PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PRDG(K)+MNUCCR(K)+PSACR(K)+&
                  PIACR(K)+PRACI(K))/(-EPRDG(K))
#else
      DUM = (-BERFDG(K)*CFL3D_TEMP(K)-PSACWG(K)*CFL3D_TEMP(K)-PRACG(K)*CLDMAX(K)-PGSACW(K)*CFL3D_TEMP(K)-(PGRACS(K)+PRDG(K)+MNUCCR(K)+EPRDG(K))*CLDMAX(K) &
            -(PIACR(K)+PRACI(K))*CFI3D_TEMP(K)-PSACR(K)*CLDMAX(K))*DT

      IF (DUM.GT.QG3D(K)*CLDMAX(K).AND.QG3D(K).GE.QSMALL) THEN

        RATIO = (QG3D(K)*CLDMAX(K)/DT+BERFDG(K)*CFL3D_TEMP(K)+PSACWG(K)*CFL3D_TEMP(K)+PRACG(K)*CLDMAX(K)+PGSACW(K)*CFL3D_TEMP(K) &
                 +(PGRACS(K)+PRDG(K)+MNUCCR(K)+PSACR(K))*CLDMAX(K)+&
                  (PIACR(K)+PRACI(K))*CFI3D_TEMP(K))/(-EPRDG(K)*CLDMAX(K))
#endif

       EPRDG(K) = EPRDG(K)*RATIO

      END IF

#ifndef CLDFRAC
! TENDENCIES

      QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)-PRD(K)-PRDS(K)-MNUCCD(K)-EPRD(K)-EPRDS(K)-PRDG(K)-EPRDG(K))

! BUG FIX HM, 3/1/11, INCLUDE PIACR AND PIACRS
      T3DTEN(K) = T3DTEN(K)+(PRE(K)                                 &
               *XXLV(K)+(PRD(K)+PRDS(K)+                            &
                MNUCCD(K)+EPRD(K)+EPRDS(K)+PRDG(K)+EPRDG(K))*XXLS(K)+         &
               (PSACWS(K)+PSACWI(K)+MNUCCC(K)+MNUCCR(K)+                      &
                QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+PRACS(K) &
                +PSACWG(K)+PRACG(K)+PGSACW(K)+PGRACS(K)+PIACR(K)+PIACRS(K))*XLF(K))/CPM(K)

      QC3DTEN(K) = QC3DTEN(K)+                                      &
                 (-PRA(K)-PRC(K)-MNUCCC(K)+PCC(K)-                  &
                  PSACWS(K)-PSACWI(K)-QMULTS(K)-QMULTG(K)-PSACWG(K)-PGSACW(K))
      QI3DTEN(K) = QI3DTEN(K)+                                      &
         (PRD(K)+EPRD(K)+PSACWI(K)+MNUCCC(K)-PRCI(K)-                                 &
                  PRAI(K)+QMULTS(K)+QMULTG(K)+QMULTR(K)+QMULTRG(K)+MNUCCD(K)-PRACI(K)-PRACIS(K))
      QR3DTEN(K) = QR3DTEN(K)+                                      &
                 (PRE(K)+PRA(K)+PRC(K)-PRACS(K)-MNUCCR(K)-QMULTR(K)-QMULTRG(K) &
             -PIACR(K)-PIACRS(K)-PRACG(K)-PGRACS(K))

      IF (IGRAUP.EQ.0) THEN

      QNI3DTEN(K) = QNI3DTEN(K)+                                    &
           (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K))
      NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K))
      QG3DTEN(K) = QG3DTEN(K)+(PRACG(K)+PSACWG(K)+PGSACW(K)+PGRACS(K)+ &
                    PRDG(K)+EPRDG(K)+MNUCCR(K)+PIACR(K)+PRACI(K)+PSACR(K))
      NG3DTEN(K) = NG3DTEN(K)+(NSCNG(K)+NGRACS(K)+NNUCCR(K)+NIACR(K))

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
      ELSE IF (IGRAUP.EQ.1) THEN

      QNI3DTEN(K) = QNI3DTEN(K)+                                    &
           (PRAI(K)+PSACWS(K)+PRDS(K)+PRACS(K)+PRCI(K)+EPRDS(K)-PSACR(K)+PIACRS(K)+PRACIS(K)+MNUCCR(K))
      NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)+NPRCI(K)-NSCNG(K)-NGRACS(K)+NIACRS(K)+NNUCCR(K))

      END IF

      NC3DTEN(K) = NC3DTEN(K)+(-NNUCCC(K)-NPSACWS(K)                &
            -NPRA(K)-NPRC(K)-NPSACWI(K)-NPSACWG(K))

      NI3DTEN(K) = NI3DTEN(K)+                                      &
       (NNUCCC(K)-NPRCI(K)-NPRAI(K)+NMULTS(K)+NMULTG(K)+NMULTR(K)+NMULTRG(K)+ &
               NNUCCD(K)-NIACR(K)-NIACRS(K))

      NR3DTEN(K) = NR3DTEN(K)+(NPRC1(K)-NPRACS(K)-NNUCCR(K)      &
                   +NRAGG(K)-NIACR(K)-NIACRS(K)-NPRACG(K)-NGRACS(K))
#else  ! frational cloudiness
! with fractional cloudiness, all tendencies are grid-mean
! TENDENCIES

      QV3DTEN(K) = QV3DTEN(K)+(-PRE(K)*CLDMAX(K)-PRD(K)*CFI3D_TEMP(K)-PRDS(K)*CLDMAX(K) &
                   -(MNUCCD(K)+EPRD(K))*CFI3D_TEMP(K)-(EPRDS(K)+PRDG(K)+EPRDG(K))*CLDMAX(K))

! BUG FIX HM, 3/1/11, INCLUDE PIACR AND PIACRS
      T3DTEN(K) = T3DTEN(K)+(PRE(K)*CLDMAX(K)                                &
               *XXLV(K)+(PRD(K)*CFI3D_TEMP(K)+PRDS(K)*CLDMAX(K)+                            &
                (MNUCCD(K)+EPRD(K))*CFI3D_TEMP(K)+(EPRDS(K)+PRDG(K)+EPRDG(K))*CLDMAX(K))*XXLS(K)+         &
               ((PSACWS(K)+PSACWI(K)+MNUCCC(K)+BERFD(K)+BERFDS(K)+BERFDG(K))*CFL3D_TEMP(K)+MNUCCR(K)*CLDMAX(K)+                      &
                (QMULTS(K)+QMULTG(K))*CFL3D_TEMP(K)+(QMULTR(K)+QMULTRG(K)+PRACS(K))*CLDMAX(K) &
                +PSACWG(K)*CFL3D_TEMP(K)+PRACG(K)*CLDMAX(K)+PGSACW(K)*CFL3D_TEMP(K)  &
                +PGRACS(K)*CLDMAX(K)+(PIACR(K)+PIACRS(K))*CFI3D_TEMP(K))*XLF(K))/CPM(K)

      QC3DTEN(K) = QC3DTEN(K)+                                      &
                 (-PRA(K)-PRC(K)-MNUCCC(K)+PCC(K)-                  &
                  PSACWS(K)-PSACWI(K)-QMULTS(K)-QMULTG(K)-PSACWG(K)-PGSACW(K)-BERFD(K)-BERFDS(K)-BERFDG(K)) * CFL3D_TEMP(K)
      QI3DTEN(K) = QI3DTEN(K)+                                      &
         ((PRD(K)+EPRD(K))*CFI3D_TEMP(K)+(PSACWI(K)+MNUCCC(K)+BERFD(K))*CFL3D_TEMP(K)-(PRCI(K)+                                 &
                  PRAI(K))*CFI3D_TEMP(K)+(QMULTS(K)+QMULTG(K))*CFL3D_TEMP(K)   & 
                  +(QMULTR(K)+QMULTRG(K))*CLDMAX(K)+(MNUCCD(K)-PRACI(K)-PRACIS(K))*CFI3D_TEMP(K))
      QR3DTEN(K) = QR3DTEN(K)+                                      &
                 (PRE(K)*CLDMAX(K)+(PRA(K)+PRC(K))*CFL3D_TEMP(K)-(PRACS(K)+MNUCCR(K)+QMULTR(K)+QMULTRG(K))*CLDMAX(K) &
             -(PIACR(K)+PIACRS(K))*CFI3D_TEMP(K)-(PRACG(K)+PGRACS(K))*CLDMAX(K))

      IF (IGRAUP.EQ.0) THEN

      QNI3DTEN(K) = QNI3DTEN(K)+                                    &
           (BERFDS(K)*CFL3D_TEMP(K)+PRAI(K)*CFI3D_TEMP(K)+PSACWS(K)*CFL3D_TEMP(K)+(PRDS(K)+PRACS(K))*CLDMAX(K)+PRCI(K)*CFI3D_TEMP(K)     &
           +(EPRDS(K)-PSACR(K))*CLDMAX(K)+(PIACRS(K)+PRACIS(K))*CFI3D_TEMP(K))
      NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)*CLDMAX(K)+NPRCI(K)*CFI3D_TEMP(K)-NSCNG(K)*CFL3D_TEMP(K) &
                   -NGRACS(K)*CLDMAX(K)+NIACRS(K)*CFI3D_TEMP(K))
      QG3DTEN(K) = QG3DTEN(K)+(BERFDG(K)*CFL3D_TEMP(K)+PRACG(K)*CLDMAX(K)+(PSACWG(K)+PGSACW(K))*CFL3D_TEMP(K)+(PGRACS(K)+ &
                    PRDG(K)+EPRDG(K)+MNUCCR(K))*CLDMAX(K)+(PIACR(K)+PRACI(K))*CFI3D_TEMP(K)+PSACR(K)*CLDMAX(K))
      NG3DTEN(K) = NG3DTEN(K)+(NSCNG(K)*CFL3D_TEMP(K)+(NGRACS(K)+NNUCCR(K))*CLDMAX(K)+NIACR(K)*CFI3D_TEMP(K))

! FOR NO GRAUPEL, NEED TO INCLUDE FREEZING OF RAIN FOR SNOW
      ELSE IF (IGRAUP.EQ.1) THEN

      QNI3DTEN(K) = QNI3DTEN(K)+                                    &
           (BERFDS(K)*CFL3D_TEMP(K)+PRAI(K)*CFI3D_TEMP(K)+PSACWS(K)*CFL3D_TEMP(K)+(PRDS(K)+PRACS(K))*CLDMAX(K)+PRCI(K)*CFI3D_TEMP(K)     &
           +(EPRDS(K)-PSACR(K))*CLDMAX(K)+(PIACRS(K)+PRACIS(K))*CFI3D_TEMP(K)+MNUCCR(K)*CLDMAX(K))
      NS3DTEN(K) = NS3DTEN(K)+(NSAGG(K)*CLDMAX(K)+NPRCI(K)*CFI3D_TEMP(K)-NSCNG(K)*CFL3D_TEMP(K)  &
                   -NGRACS(K)*CLDMAX(K)+NIACRS(K)*CFI3D_TEMP(K)+NNUCCR(K)*CLDMAX(K))

      END IF

      NC3DTEN(K) = NC3DTEN(K)+(-NNUCCC(K)-NPSACWS(K)                &
            -NPRA(K)-NPRC(K)-NPSACWI(K)-NPSACWG(K))*CFL3D_TEMP(K)

      NI3DTEN(K) = NI3DTEN(K)+                                      &
       (NNUCCC(K)*CFL3D_TEMP(K)-(NPRCI(K)+NPRAI(K))*CFI3D_TEMP(K)+(NMULTS(K)+NMULTG(K))*CFL3D_TEMP(K)  & 
       +(NMULTR(K)+NMULTRG(K))*CLDMAX(K)+ &
               (NNUCCD(K)-NIACR(K)-NIACRS(K))*CFI3D_TEMP(K))

      NR3DTEN(K) = NR3DTEN(K)+(NPRC1(K)*CFL3D_TEMP(K)+(-NPRACS(K)-NNUCCR(K)      &
                   +NRAGG(K))*CLDMAX(K)-(NIACR(K)+NIACRS(K))*CFI3D_TEMP(K)-(NPRACG(K)+NGRACS(K))*CLDMAX(K))
#endif

! V1.3 move code below to before saturation adjustment
      IF (PRE(K).LT.0.) THEN
         DUM = PRE(K)*DT/QR3D(K)
           DUM = MAX(-1.,DUM)
         NSUBR(K) = DUM*NR3D(K)/DT
      END IF

!        nsubr(k)=0.
!        nsubs(k)=0.
!        nsubg(k)=0.

#ifndef CLDFRAC
         NR3DTEN(K) = NR3DTEN(K)+NSUBR(K)
#else
         NR3DTEN(K) = NR3DTEN(K)+NSUBR(K)*CLDMAX(K)
#endif

         END IF !!!!!! TEMPERATURE

! SWITCH LTRUE TO 1, SINCE HYDROMETEORS ARE PRESENT
         LTRUE = 1

 200     CONTINUE
#ifdef CLUBB_CRM
! ADDITION BY UWM TO WEIGHT BY SGS CLOUD FRACTION
!         IF ( CF3D(K) > cloud_frac_thresh ) THEN

!           T3D(K)  = T3D_INIT + ( T3D(K) - T3D_INIT ) * CF3D(K) ! Absolute temp.
!           T3DTEN(K) = T3DTEN(K) * CF3D(K) ! Absolute temperature tendency

!           QV3D(K) = QV_INIT + ( QV3D(K) - QSAT_INIT ) * CF3D(K) ! Vapor
!           QV3DTEN(K) = QV3DTEN(K) * CF3D(K) ! Vapor mix ratio time tendency

!           QC3D(K) = QC3D(K) * CF3D(K) ! Cloud mix ratio
!           QC3DTEN(K) = QC3DTEN(K) * CF3D(K) ! Cloud mix ratio time tendency

!           IF ( INUM == 0 ) THEN
!             NC3D(K) = NC3D(K) * CF3D(K) ! Cloud drop num conc
!             NC3DTEN(K) = NC3DTEN(K) * CF3D(K) ! Cloud drop num conc time tendency
!           END IF

!           QR3D(K) = QR3D(K) * CF3D(K) ! Rain mix ratio
!           QR3DTEN(K) = QR3DTEN(K) * CF3D(K) ! Rain mix ratio time tendency

!           NR3D(K) = NR3D(K) * CF3D(K) ! Rain num conc
!           NR3DTEN(K) = NR3DTEN(K) * CF3D(K) ! Rain num conc time tendency

!           IF ( ILIQ == 0 ) THEN
!             QI3D(K) = QI3D(K) * CF3D(K) ! Ice mix ratio
!             QI3DTEN(K) = QI3DTEN(K) * CF3D(K) ! Ice mix ratio time tendency

!             NI3D(K) = NI3D(K) * CF3D(K) ! Ice num conc
!             NI3DTEN(K) = NI3DTEN(K) * CF3D(K) ! Ice num conc time tendency

!             QNI3D(K) = QNI3D(K) * CF3D(K) ! Snow mix ratio
!             QNI3DTEN(K) = QNI3DTEN(K) * CF3D(K) ! Snow mix ratio time tendency

!             NS3D(K) = NS3D(K) * CF3D(K) ! Snow num conc
!             NS3DTEN(K) = NS3DTEN(K) * CF3D(K) ! Snow num conc time tendency
!           END IF
!           IF ( IGRAUP == 0 ) THEN
!             QG3D(K)    = QG3D(K) * CF3D(K) ! Graupel mix ratio
!             QG3DTEN(K) = QG3DTEN(K) * CF3D(K) ! Graupel mix ratio time tendency

!             NG3D(K) = NG3D(K) * CF3D(K) ! Graupel num conc
!             NG3DTEN(K) = NG3DTEN(K) * CF3D(K) ! Graupel num conc time tendency
!           END IF
!         END IF ! CF3D(K) > 0.01
#endif /*CLUBB_CRM*/

        END DO

! V1.3 move precip initialization to here
! INITIALIZE PRECIP AND SNOW RATES

      PRECRT = 0.
      SNOWRT = 0.

! IF THERE ARE NO HYDROMETEORS, THEN SKIP TO END OF SUBROUTINE

        IF (LTRUE.EQ.0) GOTO 400

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!.......................................................................
! CALCULATE SEDIMENATION
! THE NUMERICS HERE FOLLOW FROM REISNER ET AL. (1998)
! FALLOUT TERMS ARE CALCULATED ON SPLIT TIME STEPS TO ENSURE NUMERICAL
! STABILITY, I.E. COURANT# < 1

!.......................................................................

      NSTEP = 1

! v3 5/27/11
      DO K = KTE,KTS,-1

#ifndef CLDFRAC
        DUMI(K) = QI3D(K)+QI3DTEN(K)*DT
        DUMQS(K) = QNI3D(K)+QNI3DTEN(K)*DT
        DUMR(K) = QR3D(K)+QR3DTEN(K)*DT
        DUMFNI(K) = NI3D(K)+NI3DTEN(K)*DT
        DUMFNS(K) = NS3D(K)+NS3DTEN(K)*DT
        DUMFNR(K) = NR3D(K)+NR3DTEN(K)*DT
        DUMC(K) = QC3D(K)+QC3DTEN(K)*DT
        DUMFNC(K) = NC3D(K)+NC3DTEN(K)*DT
	DUMG(K) = QG3D(K)+QG3DTEN(K)*DT
	DUMFNG(K) = NG3D(K)+NG3DTEN(K)*DT
#else  ! with fractional cloudiness, in-cloud values are calculated here
        DUMI(K) = (QI3D(K)*CFI3D_TEMP(K)+QI3DTEN(K)*DT)/CFI3D_TEMP(K)
        DUMQS(K) = (QNI3D(K)*CLDMAX(K)+QNI3DTEN(K)*DT)/CLDMAX(K)
        DUMR(K) = (QR3D(K)*CLDMAX(K)+QR3DTEN(K)*DT)/CLDMAX(K)
        DUMFNI(K) = (NI3D(K)*CFI3D_TEMP(K)+NI3DTEN(K)*DT)/CFI3D_TEMP(K)
        DUMFNS(K) = (NS3D(K)*CLDMAX(K)+NS3DTEN(K)*DT)/CLDMAX(K)
        DUMFNR(K) = (NR3D(K)*CLDMAX(K)+NR3DTEN(K)*DT)/CLDMAX(K)
        DUMC(K) = (QC3D(K)*CFL3D_TEMP(K)+QC3DTEN(K)*DT)/CFL3D_TEMP(K)
        DUMFNC(K) = (NC3D(K)*CFL3D_TEMP(K)+NC3DTEN(K)*DT)/CFL3D_TEMP(K)
        DUMG(K) = (QG3D(K)*CLDMAX(K)+QG3DTEN(K)*DT)/CLDMAX(K)
        DUMFNG(K) = (NG3D(K)*CLDMAX(K)+NG3DTEN(K)*DT)/CLDMAX(K)
#endif

! SWITCH FOR CONSTANT DROPLET NUMBER
        IF (INUM.EQ.1) THEN
        DUMFNC(K) = NC3D(K)
        END IF

! GET DUMMY LAMDA FOR SEDIMENTATION CALCULATIONS

! MAKE SURE NUMBER CONCENTRATIONS ARE POSITIVE
      DUMFNI(K) = MAX(0.,DUMFNI(K))
      DUMFNS(K) = MAX(0.,DUMFNS(K))
      DUMFNC(K) = MAX(0.,DUMFNC(K))
      DUMFNR(K) = MAX(0.,DUMFNR(K))
      DUMFNG(K) = MAX(0.,DUMFNG(K))

!......................................................................
! CLOUD ICE

      IF (DUMI(K).GE.QSMALL) THEN
        DLAMI = (CONS12*DUMFNI(K)/DUMI(K))**(1./DI)
        DLAMI=MAX(DLAMI,LAMMINI)
        DLAMI=MIN(DLAMI,LAMMAXI)
      END IF
!......................................................................
! RAIN

      IF (DUMR(K).GE.QSMALL) THEN
        DLAMR = (PI*RHOW*DUMFNR(K)/DUMR(K))**(1./3.)
        DLAMR=MAX(DLAMR,LAMMINR)
        DLAMR=MIN(DLAMR,LAMMAXR)
      END IF
!......................................................................
! CLOUD DROPLETS

      IF (DUMC(K).GE.QSMALL) THEN
         !bloss: option for fixing pgam
         if(dofix_pgam) then
            pgam(k) = pgam_fixed
         else

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

         end if

        DLAMC = (CONS26*DUMFNC(K)*GAMMA(PGAM(K)+4.)/(DUMC(K)*GAMMA(PGAM(K)+1.)))**(1./3.)
        LAMMIN = (PGAM(K)+1.)/60.E-6
        LAMMAX = (PGAM(K)+1.)/1.E-6
        DLAMC=MAX(DLAMC,LAMMIN)
        DLAMC=MIN(DLAMC,LAMMAX)
      END IF
!......................................................................
! SNOW

      IF (DUMQS(K).GE.QSMALL) THEN
        DLAMS = (CONS1*DUMFNS(K)/ DUMQS(K))**(1./DS)
        DLAMS=MAX(DLAMS,LAMMINS)
        DLAMS=MIN(DLAMS,LAMMAXS)
      END IF
!......................................................................
! GRAUPEL

      IF (DUMG(K).GE.QSMALL) THEN
        DLAMG = (CONS2*DUMFNG(K)/ DUMG(K))**(1./DG)
        DLAMG=MAX(DLAMG,LAMMING)
        DLAMG=MIN(DLAMG,LAMMAXG)
      END IF

!......................................................................
! CALCULATE NUMBER-WEIGHTED AND MASS-WEIGHTED TERMINAL FALL SPEEDS

! CLOUD WATER

      IF (DUMC(K).GE.QSMALL) THEN
      UNC =  ACN(K)*GAMMA(1.+BC+PGAM(K))/ (DLAMC**BC*GAMMA(PGAM(K)+1.))
      UMC = ACN(K)*GAMMA(4.+BC+PGAM(K))/  (DLAMC**BC*GAMMA(PGAM(K)+4.))
      ELSE
      UMC = 0.
      UNC = 0.
      END IF

      IF (DUMI(K).GE.QSMALL) THEN
      UNI =  AIN(K)*CONS27/DLAMI**BI
      UMI = AIN(K)*CONS28/(DLAMI**BI)
      ELSE
      UMI = 0.
      UNI = 0.
      END IF

      IF (DUMR(K).GE.QSMALL) THEN
      UNR = ARN(K)*CONS6/DLAMR**BR
      UMR = ARN(K)*CONS4/(DLAMR**BR)
      ELSE
      UMR = 0.
      UNR = 0.
      END IF

      IF (DUMQS(K).GE.QSMALL) THEN
      UMS = ASN(K)*CONS3/(DLAMS**BS)
      UNS = ASN(K)*CONS5/DLAMS**BS
      ELSE
      UMS = 0.
      UNS = 0.
      END IF

      IF (DUMG(K).GE.QSMALL) THEN
      UMG = AGN(K)*CONS7/(DLAMG**BG)
      UNG = AGN(K)*CONS8/DLAMG**BG
      ELSE
      UMG = 0.
      UNG = 0.
      END IF

! SET REALISTIC LIMITS ON FALLSPEED

! bug fix, 10/08/09
        dum=(rhosu/rho(k))**0.54
        UMS=MIN(UMS,1.2*dum)
        UNS=MIN(UNS,1.2*dum)
! v3 5/27/11
! fix for correction by AA 4/6/11
        UMI=MIN(UMI,1.2*(rhosu/rho(k))**0.35)
        UNI=MIN(UNI,1.2*(rhosu/rho(k))**0.35)
        UMR=MIN(UMR,9.1*dum)
        UNR=MIN(UNR,9.1*dum)
        UMG=MIN(UMG,20.*dum)
        UNG=MIN(UNG,20.*dum)

      FR(K) = UMR
      FI(K) = UMI
      FNI(K) = UNI
      FS(K) = UMS
      FNS(K) = UNS
      FNR(K) = UNR
      FC(K) = UMC
      FNC(K) = UNC
      FG(K) = UMG
      FNG(K) = UNG

! V3.3 MODIFY FALLSPEED BELOW LEVEL OF PRECIP

	IF (K.LE.KTE-1) THEN
        IF (FR(K).LT.1.E-10) THEN
	FR(K)=FR(K+1)
	END IF
        IF (FI(K).LT.1.E-10) THEN
	FI(K)=FI(K+1)
	END IF
        IF (FNI(K).LT.1.E-10) THEN
	FNI(K)=FNI(K+1)
	END IF
        IF (FS(K).LT.1.E-10) THEN
	FS(K)=FS(K+1)
	END IF
        IF (FNS(K).LT.1.E-10) THEN
	FNS(K)=FNS(K+1)
	END IF
        IF (FNR(K).LT.1.E-10) THEN
	FNR(K)=FNR(K+1)
	END IF
        IF (FC(K).LT.1.E-10) THEN
	FC(K)=FC(K+1)
	END IF
        IF (FNC(K).LT.1.E-10) THEN
	FNC(K)=FNC(K+1)
	END IF
        IF (FG(K).LT.1.E-10) THEN
	FG(K)=FG(K+1)
	END IF
        IF (FNG(K).LT.1.E-10) THEN
	FNG(K)=FNG(K+1)
	END IF
	END IF ! K LE KTE-1

! CALCULATE NUMBER OF SPLIT TIME STEPS

      RGVM = MAX(FR(K),FI(K),FS(K),FC(K),FNI(K),FNR(K),FNS(K),FNC(K),FG(K),FNG(K))
! VVT CHANGED IFIX -> INT (GENERIC FUNCTION)
      NSTEP = MAX(INT(RGVM*DT/DZQ(K)+1.),NSTEP)

! MULTIPLY VARIABLES BY RHO
#ifndef CLDFRAC
      DUMR(k) = DUMR(k)*RHO(K)
      DUMI(k) = DUMI(k)*RHO(K)
      DUMFNI(k) = DUMFNI(K)*RHO(K)
      DUMQS(k) = DUMQS(K)*RHO(K)
      DUMFNS(k) = DUMFNS(K)*RHO(K)
      DUMFNR(k) = DUMFNR(K)*RHO(K)
      DUMC(k) = DUMC(K)*RHO(K)
      DUMFNC(k) = DUMFNC(K)*RHO(K)
      DUMG(k) = DUMG(K)*RHO(K)
      DUMFNG(k) = DUMFNG(K)*RHO(K)
#else  ! with fractional cloudiness, sedimentatin calculation is applied to grid-mean mixing ratios
      DUMR(k) = DUMR(k)*RHO(K)*CLDMAX(K)
      DUMI(k) = DUMI(k)*RHO(K)*CFI3D_TEMP(K)
      DUMFNI(k) = DUMFNI(K)*RHO(K)*CFI3D_TEMP(K)
      DUMQS(k) = DUMQS(K)*RHO(K)*CLDMAX(K)
      DUMFNS(k) = DUMFNS(K)*RHO(K)*CLDMAX(K)
      DUMFNR(k) = DUMFNR(K)*RHO(K)*CLDMAX(K)
      DUMC(k) = DUMC(K)*RHO(K)*CFL3D_TEMP(K)
      DUMFNC(k) = DUMFNC(K)*RHO(K)*CFL3D_TEMP(K)
      DUMG(k) = DUMG(K)*RHO(K)*CLDMAX(K)
      DUMFNG(k) = DUMFNG(K)*RHO(K)*CLDMAX(K)
#endif

      END DO

      DO N = 1,NSTEP

      DO K = KTS,KTE
      FALOUTR(K) = FR(K)*DUMR(K)
      FALOUTI(K) = FI(K)*DUMI(K)
      FALOUTNI(K) = FNI(K)*DUMFNI(K)
      FALOUTS(K) = FS(K)*DUMQS(K)
      FALOUTNS(K) = FNS(K)*DUMFNS(K)
      FALOUTNR(K) = FNR(K)*DUMFNR(K)
      FALOUTC(K) = FC(K)*DUMC(K)
      FALOUTNC(K) = FNC(K)*DUMFNC(K)
      FALOUTG(K) = FG(K)*DUMG(K)
      FALOUTNG(K) = FNG(K)*DUMFNG(K)
      END DO

! TOP OF MODEL

      K = KTE
      FALTNDR = FALOUTR(K)/DZQ(k)
      FALTNDI = FALOUTI(K)/DZQ(k)
      FALTNDNI = FALOUTNI(K)/DZQ(k)
      FALTNDS = FALOUTS(K)/DZQ(k)
      FALTNDNS = FALOUTNS(K)/DZQ(k)
      FALTNDNR = FALOUTNR(K)/DZQ(k)
      FALTNDC = FALOUTC(K)/DZQ(k)
      FALTNDNC = FALOUTNC(K)/DZQ(k)
      FALTNDG = FALOUTG(K)/DZQ(k)
      FALTNDNG = FALOUTNG(K)/DZQ(k)
! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

      QRSTEN(K) = QRSTEN(K)-FALTNDR/NSTEP/RHO(k)
      QISTEN(K) = QISTEN(K)-FALTNDI/NSTEP/RHO(k)
      NI3DTEN(K) = NI3DTEN(K)-FALTNDNI/NSTEP/RHO(k)
      QNISTEN(K) = QNISTEN(K)-FALTNDS/NSTEP/RHO(k)
      NS3DTEN(K) = NS3DTEN(K)-FALTNDNS/NSTEP/RHO(k)
      NR3DTEN(K) = NR3DTEN(K)-FALTNDNR/NSTEP/RHO(k)
      QCSTEN(K) = QCSTEN(K)-FALTNDC/NSTEP/RHO(k)
      NC3DTEN(K) = NC3DTEN(K)-FALTNDNC/NSTEP/RHO(k)
      QGSTEN(K) = QGSTEN(K)-FALTNDG/NSTEP/RHO(k)
      NG3DTEN(K) = NG3DTEN(K)-FALTNDNG/NSTEP/RHO(k)

      NISTEN(K) = NISTEN(K)-FALTNDNI/NSTEP/RHO(k)
      NSSTEN(K) = NSSTEN(K)-FALTNDNS/NSTEP/RHO(k)
      NRSTEN(K) = NRSTEN(K)-FALTNDNR/NSTEP/RHO(k)
      NCSTEN(K) = NCSTEN(K)-FALTNDNC/NSTEP/RHO(k)
      NGSTEN(K) = NGSTEN(K)-FALTNDNG/NSTEP/RHO(k)

      DUMR(K) = DUMR(K)-FALTNDR*DT/NSTEP
      DUMI(K) = DUMI(K)-FALTNDI*DT/NSTEP
      DUMFNI(K) = DUMFNI(K)-FALTNDNI*DT/NSTEP
      DUMQS(K) = DUMQS(K)-FALTNDS*DT/NSTEP
      DUMFNS(K) = DUMFNS(K)-FALTNDNS*DT/NSTEP
      DUMFNR(K) = DUMFNR(K)-FALTNDNR*DT/NSTEP
      DUMC(K) = DUMC(K)-FALTNDC*DT/NSTEP
      DUMFNC(K) = DUMFNC(K)-FALTNDNC*DT/NSTEP
      DUMG(K) = DUMG(K)-FALTNDG*DT/NSTEP
      DUMFNG(K) = DUMFNG(K)-FALTNDNG*DT/NSTEP

      DO K = KTE-1,KTS,-1
#ifndef CLDFRAC
      FALTNDR = (FALOUTR(K+1)-FALOUTR(K))/DZQ(K)
      FALTNDI = (FALOUTI(K+1)-FALOUTI(K))/DZQ(K)
      FALTNDNI = (FALOUTNI(K+1)-FALOUTNI(K))/DZQ(K)
      FALTNDS = (FALOUTS(K+1)-FALOUTS(K))/DZQ(K)
      FALTNDNS = (FALOUTNS(K+1)-FALOUTNS(K))/DZQ(K)
      FALTNDNR = (FALOUTNR(K+1)-FALOUTNR(K))/DZQ(K)
      FALTNDC = (FALOUTC(K+1)-FALOUTC(K))/DZQ(K)
      FALTNDNC = (FALOUTNC(K+1)-FALOUTNC(K))/DZQ(K)
      FALTNDG = (FALOUTG(K+1)-FALOUTG(K))/DZQ(K)
      FALTNDNG = (FALOUTNG(K+1)-FALOUTNG(K))/DZQ(K)
#else
! For fractional cloudniess, cloud ice and cloud water falls into clear-sky is 
! assumed to evaporate instantly
      DUMCFL =  CFL3D_TEMP(K)/CFL3D_TEMP(K+1)
      DUMCFL = min(DUMCFL, 1.)
      DUMCFI =  CFI3D_TEMP(K)/CFI3D_TEMP(K+1)
      DUMCFI = min(DUMCFI, 1.)

      FALTNDR = (FALOUTR(K+1)-FALOUTR(K))/DZQ(K)

      FALTNDI = (FALOUTI(K+1)*DUMCFI-FALOUTI(K))/DZQ(K)
      FALTNDNI = (FALOUTNI(K+1)*DUMCFI-FALOUTNI(K))/DZQ(K)

      FALTNDS = (FALOUTS(K+1)-FALOUTS(K))/DZQ(K)
      FALTNDNS = (FALOUTNS(K+1)-FALOUTNS(K))/DZQ(K)
      FALTNDNR = (FALOUTNR(K+1)-FALOUTNR(K))/DZQ(K)

      FALTNDC = (FALOUTC(K+1)*DUMCFL-FALOUTC(K))/DZQ(K)
      FALTNDNC = (FALOUTNC(K+1)*DUMCFL-FALOUTNC(K))/DZQ(K)

      FALTNDG = (FALOUTG(K+1)-FALOUTG(K))/DZQ(K)
      FALTNDNG = (FALOUTNG(K+1)-FALOUTNG(K))/DZQ(K)
      
! Evaporation of cloud ice and cloud liquid  
      QISEVAP(K) = QISEVAP(K)+FALOUTC(K+1)*(1-DUMCFL)/DZQ(K)/NSTEP/RHO(k)
      QCSEVAP(K) = QCSEVAP(K)+FALOUTI(K+1)*(1-DUMCFI)/DZQ(K)/NSTEP/RHO(k)
#endif

! ADD FALLOUT TERMS TO EULERIAN TENDENCIES

      QRSTEN(K) = QRSTEN(K)+FALTNDR/NSTEP/RHO(k)
      QISTEN(K) = QISTEN(K)+FALTNDI/NSTEP/RHO(k)
      NI3DTEN(K) = NI3DTEN(K)+FALTNDNI/NSTEP/RHO(k)
      QNISTEN(K) = QNISTEN(K)+FALTNDS/NSTEP/RHO(k)
      NS3DTEN(K) = NS3DTEN(K)+FALTNDNS/NSTEP/RHO(k)
      NR3DTEN(K) = NR3DTEN(K)+FALTNDNR/NSTEP/RHO(k)
      QCSTEN(K) = QCSTEN(K)+FALTNDC/NSTEP/RHO(k)
      NC3DTEN(K) = NC3DTEN(K)+FALTNDNC/NSTEP/RHO(k)
      QGSTEN(K) = QGSTEN(K)+FALTNDG/NSTEP/RHO(k)
      NG3DTEN(K) = NG3DTEN(K)+FALTNDNG/NSTEP/RHO(k)

      NISTEN(K) = NISTEN(K)+FALTNDNI/NSTEP/RHO(k)
      NSSTEN(K) = NSSTEN(K)+FALTNDNS/NSTEP/RHO(k)
      NRSTEN(K) = NRSTEN(K)+FALTNDNR/NSTEP/RHO(k)
      NCSTEN(K) = NCSTEN(K)+FALTNDNC/NSTEP/RHO(k)
      NGSTEN(K) = NGSTEN(K)+FALTNDNG/NSTEP/RHO(k)

      DUMR(K) = DUMR(K)+FALTNDR*DT/NSTEP
      DUMI(K) = DUMI(K)+FALTNDI*DT/NSTEP
      DUMFNI(K) = DUMFNI(K)+FALTNDNI*DT/NSTEP
      DUMQS(K) = DUMQS(K)+FALTNDS*DT/NSTEP
      DUMFNS(K) = DUMFNS(K)+FALTNDNS*DT/NSTEP
      DUMFNR(K) = DUMFNR(K)+FALTNDNR*DT/NSTEP
      DUMC(K) = DUMC(K)+FALTNDC*DT/NSTEP
      DUMFNC(K) = DUMFNC(K)+FALTNDNC*DT/NSTEP
      DUMG(K) = DUMG(K)+FALTNDG*DT/NSTEP
      DUMFNG(K) = DUMFNG(K)+FALTNDNG*DT/NSTEP

#ifdef ECPP
      RSED(K)=RSED(K)+FALOUTR(K)/NSTEP
      ISED(K)=ISED(K)+FALOUTI(K)/NSTEP
      CSED(K)=CSED(K)+FALOUTC(K)/NSTEP
      SSED(K)=SSED(K)+FALOUTS(K)/NSTEP
      GSED(K)=GSED(K)+FALOUTG(K)/NSTEP
#endif

      END DO

! GET PRECIPITATION AND SNOWFALL ACCUMULATION DURING THE TIME STEP
! FACTOR OF 1000 CONVERTS FROM M TO MM, BUT DIVISION BY DENSITY
! OF LIQUID WATER CANCELS THIS FACTOR OF 1000

        PRECRT = PRECRT+(FALOUTR(KTS)+FALOUTC(KTS)+FALOUTS(KTS)+FALOUTI(KTS)+FALOUTG(KTS))  &
                     *DT/NSTEP
        SNOWRT = SNOWRT+(FALOUTS(KTS)+FALOUTI(KTS)+FALOUTG(KTS))*DT/NSTEP

      END DO

      DO K=KTS,KTE

! ADD ON SEDIMENTATION TENDENCIES FOR MIXING RATIO TO REST OF TENDENCIES

        QR3DTEN(K)=QR3DTEN(K)+QRSTEN(K)
        QI3DTEN(K)=QI3DTEN(K)+QISTEN(K)
        QC3DTEN(K)=QC3DTEN(K)+QCSTEN(K)
        QG3DTEN(K)=QG3DTEN(K)+QGSTEN(K)
        QNI3DTEN(K)=QNI3DTEN(K)+QNISTEN(K)
#ifdef CLDFRAC
        QV3DTEN(K)=QV3DTEN(K)+QISEVAP(K)+QCSEVAP(K)
        T3DTEN(K)=T3DTEN(K)-(QISEVAP(K)*XXLS(K)+QCSEVAP(K)*XXLV(K))/CPM(K)
#endif

! PUT ALL CLOUD ICE IN SNOW CATEGORY IF MEAN DIAMETER EXCEEDS 2 * dcs

! V1.7
!hm 7/9/09 bug fix
#ifndef CLDFRAC
!        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.273.15) THEN
        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.TMELT.AND.LAMI(K).GE.1.E-10) THEN

        IF (1./LAMI(K).GE.2.*DCS) THEN
           QNI3DTEN(K) = QNI3DTEN(K)+QI3D(K)/DT+ QI3DTEN(K)
           NS3DTEN(K) = NS3DTEN(K)+NI3D(K)/DT+   NI3DTEN(K)
           QI3DTEN(K) = -QI3D(K)/DT
           NI3DTEN(K) = -NI3D(K)/DT
        END IF
        END IF
#else
!        IF (QI3D(K)*CFI3D_TEMP(K).GE.QSMALL.AND.T3D(K).LT.273.15) THEN
        IF (QI3D(K).GE.QSMALL.AND.T3D(K).LT.TMELT.AND.LAMI(K).GE.1.E-10) THEN

        IF (1./LAMI(K).GE.2.*DCS) THEN
           QNI3DTEN(K) = QNI3DTEN(K)+QI3D(K)*CFI3D_TEMP(K)/DT+ QI3DTEN(K)
           NS3DTEN(K) = NS3DTEN(K)+NI3D(K)*CFI3D_TEMP(K)/DT+   NI3DTEN(K)
           QI3DTEN(K) = -QI3D(K)*CFI3D_TEMP(K)/DT
           NI3DTEN(K) = -NI3D(K)*CFI3D_TEMP(K)/DT
        END IF
        END IF
#endif

! hm add tendencies here, then call sizeparameter
! to ensure consisitency between mixing ratio and number concentration

#ifndef CLDFRAC
          QC3D(k)        = QC3D(k)+QC3DTEN(k)*DT
          QI3D(k)        = QI3D(k)+QI3DTEN(k)*DT
          QNI3D(k)        = QNI3D(k)+QNI3DTEN(k)*DT
          QR3D(k)        = QR3D(k)+QR3DTEN(k)*DT
          NC3D(k)        = NC3D(k)+NC3DTEN(k)*DT
          NI3D(k)        = NI3D(k)+NI3DTEN(k)*DT
          NS3D(k)        = NS3D(k)+NS3DTEN(k)*DT
          NR3D(k)        = NR3D(k)+NR3DTEN(k)*DT

          IF (IGRAUP.EQ.0) THEN
          QG3D(k)        = QG3D(k)+QG3DTEN(k)*DT
          NG3D(k)        = NG3D(k)+NG3DTEN(k)*DT
          END IF
#else
! with fractional cloudiness, hydrometer mixing ratios are still in-cloud values
          QC3D(k)        = (QC3D(k)*CFL3D_TEMP(K)+QC3DTEN(k)*DT)/CFL3D_TEMP(K)
          QI3D(k)        = (QI3D(k)*CFI3D_TEMP(K)+QI3DTEN(k)*DT)/CFI3D_TEMP(K)
          QNI3D(k)        = (QNI3D(k)*CLDMAX(K)+QNI3DTEN(k)*DT)/CLDMAX(K)
          QR3D(k)        = (QR3D(k)*CLDMAX(K)+QR3DTEN(k)*DT)/CLDMAX(K)
          NC3D(k)        = (NC3D(k)*CFL3D_TEMP(K)+NC3DTEN(k)*DT)/CFL3D_TEMP(K)
          NI3D(k)        = (NI3D(k)*CFI3D_TEMP(K)+NI3DTEN(k)*DT)/CFI3D_TEMP(K)
          NS3D(k)        = (NS3D(k)*CLDMAX(K)+NS3DTEN(k)*DT)/CLDMAX(K)
          NR3D(k)        = (NR3D(k)*CLDMAX(K)+NR3DTEN(k)*DT)/CLDMAX(K)

          IF (IGRAUP.EQ.0) THEN
          QG3D(k)        = (QG3D(k)*CLDMAX(K)+QG3DTEN(k)*DT)/CLDMAX(K)
          NG3D(k)        = (NG3D(k)*CLDMAX(K)+NG3DTEN(k)*DT)/CLDMAX(K)
          END IF
! for testing purpose +++mhwang
          NI3D_T1(K) = NI3D(k)
          NI3D_T2(K) = NI3DTEN(k)*DT
          NR3D_T2(K) = NR3D(K) * CLDMAX(K)
#endif

! ADD TEMPERATURE AND WATER VAPOR TENDENCIES FROM MICROPHYSICS
          T3D(K)         = T3D(K)+T3DTEN(k)*DT
          QV3D(K)        = QV3D(K)+QV3DTEN(k)*DT

! SATURATION VAPOR PRESSURE AND MIXING RATIO

! hm, add fix for low pressure, 5/12/10
            EVS(K) = min(0.99*pres(k),POLYSVP(T3D(K),0))   ! PA
            EIS(K) = min(0.99*pres(k),POLYSVP(T3D(K),1))   ! PA

! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING

            IF (EIS(K).GT.EVS(K)) EIS(K) = EVS(K)

            QVS(K) = EP_2*EVS(K)/(PRES(K)-EVS(K))
            QVI(K) = EP_2*EIS(K)/(PRES(K)-EIS(K))

            QVQVS(K) = QV3D(K)/QVS(K)
            QVQVSI(K) = QV3D(K)/QVI(K)

! AT SUBSATURATION, REMOVE SMALL AMOUNTS OF CLOUD/PRECIP WATER

! V1.3, change limit from 10^-7 to 10^-6
! V1.7 7/9/09 change limit from 10^-6 to 10^-8

#ifndef CLDFRAC    !For fractional cloudiness, QVQVS can be less than 0.9
             IF (QVQVS(K).LT.0.9) THEN
               IF (QR3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QR3D(K)
                  T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
                  QR3D(K)=0.
               END IF
               IF (QC3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QC3D(K)
                  T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
                  QC3D(K)=0.
               END IF
             END IF

             IF (QVQVSI(K).LT.0.9) THEN
               IF (QI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QI3D(K)
                  T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
                  QI3D(K)=0.
               END IF
               IF (QNI3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QNI3D(K)
                  T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
                  QNI3D(K)=0.
               END IF
               IF (QG3D(K).LT.1.E-8) THEN
                  QV3D(K)=QV3D(K)+QG3D(K)
                  T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
                  QG3D(K)=0.
               END IF
             END IF
#endif

!..................................................................
! IF MIXING RATIO < QSMALL SET MIXING RATIO AND NUMBER CONC TO ZERO

       IF (QC3D(K).LT.QSMALL) THEN
!+++mhwang
#ifndef CLDFRAC
         QV3D(K)=QV3D(K)+QC3D(K)
         T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)
#else
!   QV3D and T3D are grid-mean values while QC3D is in-cloud value, 
!   so cloud fraction is needed here
         QV3D(K)=QV3D(K)+QC3D(K)*CFL3D_TEMP(K)
         T3D(K)=T3D(K)-QC3D(K)*XXLV(K)/CPM(K)*CFL3D_TEMP(K)
#endif
!---mhwang
         QC3D(K) = 0.
         NC3D(K) = 0.
         EFFC(K) = 0.
       END IF
       IF (QR3D(K).LT.QSMALL) THEN
!+++mhwang
#ifndef CLDFRAC
         QV3D(K)=QV3D(K)+QR3D(K)
         T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)
#else
         QV3D(K)=QV3D(K)+QR3D(K)*CLDMAX(K)
         T3D(K)=T3D(K)-QR3D(K)*XXLV(K)/CPM(K)*CLDMAX(K)
#endif
!---mhwang
         QR3D(K) = 0.
         NR3D(K) = 0.
         EFFR(K) = 0.
       END IF
       IF (QI3D(K).LT.QSMALL) THEN
!+++mhwang
#ifndef CLDFRAC
         QV3D(K)=QV3D(K)+QI3D(K)
         T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)
#else
         QV3D(K)=QV3D(K)+QI3D(K)*CFI3D_TEMP(K)
         T3D(K)=T3D(K)-QI3D(K)*XXLS(K)/CPM(K)*CFI3D_TEMP(K)
#endif 
!+++mhwang
         QI3D(K) = 0.
         NI3D(K) = 0.
         EFFI(K) = 0.
       END IF
       IF (QNI3D(K).LT.QSMALL) THEN
!+++mhwang
#ifndef CLDFRAC
         QV3D(K)=QV3D(K)+QNI3D(K)
         T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)
#else
         QV3D(K)=QV3D(K)+QNI3D(K)*CLDMAX(K)
         T3D(K)=T3D(K)-QNI3D(K)*XXLS(K)/CPM(K)*CLDMAX(K)
#endif
!+++mhwang
         QNI3D(K) = 0.
         NS3D(K) = 0.
         EFFS(K) = 0.
       END IF
       IF (QG3D(K).LT.QSMALL) THEN
!+++mhwang
#ifndef CLDFRAC
         QV3D(K)=QV3D(K)+QG3D(K)
         T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)
#else
         QV3D(K)=QV3D(K)+QG3D(K)*CLDMAX(K)
         T3D(K)=T3D(K)-QG3D(K)*XXLS(K)/CPM(K)*CLDMAX(K)
#endif
!+++mhwang
         QG3D(K) = 0.
         NG3D(K) = 0.
         EFFG(K) = 0.
       END IF

!..................................
! IF THERE IS NO CLOUD/PRECIP WATER, THEN SKIP CALCULATIONS

            IF (QC3D(K).LT.QSMALL.AND.QI3D(K).LT.QSMALL.AND.QNI3D(K).LT.QSMALL &
                 .AND.QR3D(K).LT.QSMALL.AND.QG3D(K).LT.QSMALL) GOTO 500

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CALCULATE INSTANTANEOUS PROCESSES

! ADD MELTING OF CLOUD ICE TO FORM RAIN

        IF (QI3D(K).GE.QSMALL.AND.T3D(K).GE.TMELT) THEN
#ifndef CLDFRAC
           QR3D(K) = QR3D(K)+QI3D(K)
           T3D(K) = T3D(K)-QI3D(K)*XLF(K)/CPM(K)
           QI3D(K) = 0.
           NR3D(K) = NR3D(K)+NI3D(K)
           NI3D(K) = 0.
#else
           QR3D(K) = QR3D(K)+QI3D(K)*CFI3D_TEMP(K)/CLDMAX(K)
           T3D(K) = T3D(K)-QI3D(K)*XLF(K)/CPM(K)*CFI3D_TEMP(K)
           QI3D(K) = 0.
           NR3D(K) = NR3D(K)+NI3D(K)*CFI3D_TEMP(K)/CLDMAX(K)
           NI3D(K) = 0.

           NR3D_T3(K) = NR3D(K) * CLDMAX(K)  !+++mhwang test
#endif
           MELTQI2QR(K) = QI3D(K)/DT
        END IF

! ****SENSITIVITY - NO ICE
        IF (ILIQ.EQ.1) GOTO 778
!
#ifdef CLDFRAC
! This additional check for Rain and droplet number is added by Minghuai Wang
! With fractional cloudiness, droplet number concentration may not match with mass
      NC3D(K) = MAX(0., NC3D(K))
      NR3D(K) = MAX(0., NR3D(K))
!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF

      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         !bloss: option for fixing pgam
         if(dofix_pgam) then
            pgam(k) = pgam_fixed
         else

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

         end if

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

      END IF
      QC3D_T3(K) = 0.0
      NITEN_HOMC(K) =0.0
      NC3D_T3(K) =0.0
      NI3D_T3(K) = 0.0
      QC3D_T2(K) = QC3D(K)
#endif


! HOMOGENEOUS FREEZING OF CLOUD WATER

        IF (T3D(K).LE.233.15.AND.QC3D(K).GE.QSMALL) THEN
#ifndef CLDFRAC
           QI3D(K)=QI3D(K)+QC3D(K)
           T3D(K)=T3D(K)+QC3D(K)*XLF(K)/CPM(K)
           QC3D(K)=0.
           NI3D(K)=NI3D(K)+NC3D(K)
           NC3D(K)=0.
#else
!+++mhwang
! for testing purpose
           QC3D_T3(K) = QC3D(K)
!---mhwang

           QI3D(K)=QI3D(K)+QC3D(K)*CFL3D_TEMP(K)/CFI3D_TEMP(K)
           T3D(K)=T3D(K)+QC3D(K)*XLF(K)/CPM(K)*CFL3D_TEMP(K)
           QC3D(K)=0.
!+++mhwang
! for testing purpose
#ifdef CLUBB_CRM
           NITEN_HOMC(K) = NNUCHOM_REDUCE_COEF*NC3D(K)*CFL3D_TEMP(K)/CFI3D_TEMP(K)/DT
#else
           NITEN_HOMC(K) = NC3D(K)*CFL3D_TEMP(K)/CFI3D_TEMP(K)/DT
#endif

           NC3D_T3(K) = NC3D(K)
!---mhwang

#ifdef CLUBB_CRM
           NI3D(K)=NI3D(K)+NNUCHOM_REDUCE_COEF*NC3D(K)*CFL3D_TEMP(K)/CFI3D_TEMP(K)
#else
           NI3D(K)=NI3D(K)+NC3D(K)*CFL3D_TEMP(K)/CFI3D_TEMP(K)
#endif
           NC3D(K)=0.

! for testing purpose
           NI3D_T3(K) = NI3D(K)
#endif
        END IF

! HOMOGENEOUS FREEZING OF RAIN

        IF (IGRAUP.EQ.0) THEN

        IF (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) THEN
           QG3D(K) = QG3D(K)+QR3D(K)
#ifndef CLDFRAC
           T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
#else
           T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K) * CLDMAX(K)
#endif
           QR3D(K) = 0.
           NG3D(K) = NG3D(K)+ NR3D(K)
           NR3D(K) = 0.
        END IF

        ELSE IF (IGRAUP.EQ.1) THEN

        IF (T3D(K).LE.233.15.AND.QR3D(K).GE.QSMALL) THEN
           QNI3D(K) = QNI3D(K)+QR3D(K)
#ifndef CLDFRAC
           T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)
#else
           T3D(K) = T3D(K)+QR3D(K)*XLF(K)/CPM(K)*CLDMAX(K)
#endif
           QR3D(K) = 0.
           NS3D(K) = NS3D(K)+NR3D(K)
           NR3D(K) = 0.
        END IF

        END IF

 778    CONTINUE

! MAKE SURE NUMBER CONCENTRATIONS AREN'T NEGATIVE

      NI3D(K) = MAX(0.,NI3D(K))
      NS3D(K) = MAX(0.,NS3D(K))
      NC3D(K) = MAX(0.,NC3D(K))
      NR3D(K) = MAX(0.,NR3D(K))
      NG3D(K) = MAX(0.,NG3D(K))

!......................................................................
! CLOUD ICE

      IF (QI3D(K).GE.QSMALL) THEN
         LAMI(K) = (CONS12*                 &
              NI3D(K)/QI3D(K))**(1./DI)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMI(K).LT.LAMMINI) THEN

      LAMI(K) = LAMMINI

      N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      ELSE IF (LAMI(K).GT.LAMMAXI) THEN
      LAMI(K) = LAMMAXI
      N0I(K) = LAMI(K)**(DI+1.)*QI3D(K)/CONS12

      NI3D(K) = N0I(K)/LAMI(K)
      END IF
      END IF

!+++mhwang
      NI3D_T4(K) = NI3D(K)
!---mhwang 

!......................................................................
! RAIN

      IF (QR3D(K).GE.QSMALL) THEN
      LAMR(K) = (PI*RHOW*NR3D(K)/QR3D(K))**(1./3.)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMR(K).LT.LAMMINR) THEN

      LAMR(K) = LAMMINR

      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      ELSE IF (LAMR(K).GT.LAMMAXR) THEN
      LAMR(K) = LAMMAXR
      N0RR(K) = LAMR(K)**4*QR3D(K)/(PI*RHOW)

      NR3D(K) = N0RR(K)/LAMR(K)
      END IF

      END IF

!......................................................................
! CLOUD DROPLETS

! MARTIN ET AL. (1994) FORMULA FOR PGAM

      IF (QC3D(K).GE.QSMALL) THEN

         !bloss: option for fixing pgam
         if(dofix_pgam) then
            pgam(k) = pgam_fixed
         else

!         DUM = PRES(K)/(R*T3D(K))
! V1.5
         PGAM(K)=0.0005714*(NC3D(K)/1.E6*RHO(K))+0.2714
         PGAM(K)=1./(PGAM(K)**2)-1.
         PGAM(K)=MAX(PGAM(K),2.)
         PGAM(K)=MIN(PGAM(K),10.)

         end if

! CALCULATE LAMC

      LAMC(K) = (CONS26*NC3D(K)*GAMMA(PGAM(K)+4.)/   &
                 (QC3D(K)*GAMMA(PGAM(K)+1.)))**(1./3.)

! LAMMIN, 60 MICRON DIAMETER
! LAMMAX, 1 MICRON

      LAMMIN = (PGAM(K)+1.)/60.E-6
      LAMMAX = (PGAM(K)+1.)/1.E-6

      IF (LAMC(K).LT.LAMMIN) THEN
      LAMC(K) = LAMMIN
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      ELSE IF (LAMC(K).GT.LAMMAX) THEN
      LAMC(K) = LAMMAX
      NC3D(K) = EXP(3.*LOG(LAMC(K))+LOG(QC3D(K))+              &
                LOG(GAMMA(PGAM(K)+1.))-LOG(GAMMA(PGAM(K)+4.)))/CONS26

      END IF

      END IF

!......................................................................
! SNOW

      IF (QNI3D(K).GE.QSMALL) THEN
      LAMS(K) = (CONS1*NS3D(K)/QNI3D(K))**(1./DS)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMS(K).LT.LAMMINS) THEN
      LAMS(K) = LAMMINS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1

      NS3D(K) = N0S(K)/LAMS(K)

      ELSE IF (LAMS(K).GT.LAMMAXS) THEN

      LAMS(K) = LAMMAXS
      N0S(K) = LAMS(K)**(DS+1.)*QNI3D(K)/CONS1
      NS3D(K) = N0S(K)/LAMS(K)
      END IF

      END IF

!......................................................................
! GRAUPEL

      IF (QG3D(K).GE.QSMALL) THEN
      LAMG(K) = (CONS2*NG3D(K)/QG3D(K))**(1./DG)

! CHECK FOR SLOPE

! ADJUST VARS

      IF (LAMG(K).LT.LAMMING) THEN
      LAMG(K) = LAMMING
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)

      ELSE IF (LAMG(K).GT.LAMMAXG) THEN

      LAMG(K) = LAMMAXG
      N0G(K) = LAMG(K)**(DG+1.)*QG3D(K)/CONS2

      NG3D(K) = N0G(K)/LAMG(K)
      END IF

      END IF

 500  CONTINUE

! CALCULATE EFFECTIVE RADIUS

!#ifdef CLUBB_CRM
      ! Account for subgrid scale effective droplet radii 
!      IF ( CF3D(K) > cloud_frac_thresh ) THEN 
!        TMPQSMALL = QSMALL / CF3D(K)
!      ELSE
!        TMPQSMALL = QSMALL
!      END IF 

!      IF (QI3D(K).GE.TMPQSMALL) THEN
!         EFFI(K) = 3./LAMI(K)/2.*1.E6
!      ELSE
!         EFFI(K) = 25.
!      END IF

!      IF (QNI3D(K).GE.TMPQSMALL) THEN
!         EFFS(K) = 3./LAMS(K)/2.*1.E6
!      ELSE
!         EFFS(K) = 25.
!      END IF

!      IF (QR3D(K).GE.TMPQSMALL) THEN
!         EFFR(K) = 3./LAMR(K)/2.*1.E6
!      ELSE
!         EFFR(K) = 25.
!      END IF

!      IF (QC3D(K).GE.TMPQSMALL) THEN
!      EFFC(K) = GAMMA(PGAM(K)+4.)/                        &
!             GAMMA(PGAM(K)+3.)/LAMC(K)/2.*1.E6
!      ELSE
!      EFFC(K) = 25.
!      END IF

!      IF (QG3D(K).GE.TMPQSMALL) THEN
!         EFFG(K) = 3./LAMG(K)/2.*1.E6
!      ELSE
!         EFFG(K) = 25.
!      END IF
!#else
      IF (QI3D(K).GE.QSMALL) THEN
         EFFI(K) = 3./LAMI(K)/2.*1.E6
      ELSE
         EFFI(K) = 25.
      END IF

      IF (QNI3D(K).GE.QSMALL) THEN
         EFFS(K) = 3./LAMS(K)/2.*1.E6
      ELSE
         EFFS(K) = 25.
      END IF

      IF (QR3D(K).GE.QSMALL) THEN
         EFFR(K) = 3./LAMR(K)/2.*1.E6
      ELSE
         EFFR(K) = 25.
      END IF

      IF (QC3D(K).GE.QSMALL) THEN
      EFFC(K) = GAMMA(PGAM(K)+4.)/                        &
             GAMMA(PGAM(K)+3.)/LAMC(K)/2.*1.E6
      ELSE
      EFFC(K) = 25.
      END IF

      IF (QG3D(K).GE.QSMALL) THEN
         EFFG(K) = 3./LAMG(K)/2.*1.E6
      ELSE
         EFFG(K) = 25.
      END IF
!#endif /*CLUBB*/

! HM ADD 1/10/06, ADD UPPER BOUND ON ICE NUMBER, THIS IS NEEDED
! TO PREVENT VERY LARGE ICE NUMBER DUE TO HOMOGENEOUS FREEZING
! OF DROPLETS, ESPECIALLY WHEN INUM = 1, SET MAX AT 10 CM-3
          NI3D(K) = MIN(NI3D(K),10.E6/RHO(K))
! ADD BOUND ON DROPLET NUMBER - CANNOT EXCEED AEROSOL CONCENTRATION
          IF (INUM.EQ.0.AND.IACT.EQ.2) THEN
          NC3D(K) = MIN(NC3D(K),(NANEW1+NANEW2)/RHO(K))
          END IF
! SWITCH FOR CONSTANT DROPLET NUMBER
          IF (INUM.EQ.1) THEN
! CHANGE NDCNST FROM CM-3 TO KG-1
             NC3D(K) = NDCNST*1.E6/RHO(K)
          END IF

      END DO !!! K LOOP

 400         CONTINUE

      DO K=KTS, KTE
#ifdef CLDFRAC
          QC3D(K) = QC3D(K)*CFL3D_TEMP(K)
          QI3D(K) = QI3D(K)*CFI3D_TEMP(K)
          QNI3D(K) = QNI3D(K)*CLDMAX(K)
          QR3D(K) = QR3D(K)*CLDMAX(K)
          QG3D(K) = QG3D(K)*CLDMAX(K)
          NC3D(K) = NC3D(K)*CFL3D_TEMP(K)
          NI3D(K) = NI3D(K)*CFI3D_TEMP(K)
          NS3D(K) = NS3D(K)*CLDMAX(K)
          NR3D(K) = NR3D(K)*CLDMAX(K)
          NG3D(K) = NG3D(K)*CLDMAX(K)

! individual process rates
          PRC(K) = PRC(K) * CFL3D_TEMP(K)
          PRA(K) = PRA(K) * CFL3D_TEMP(K)
          PSMLT(K) = PSMLT(K) * CLDMAX(K)
          EVPMS(K) = EVPMS(K) * CLDMAX(K)
          PRACS(K) = PRACS(K) * CLDMAX(K)
          EVPMG(K) = EVPMG(K) * CLDMAX(K)
          PRACG(K) = PRACG(K) * CLDMAX(K)
          PRE(K) = PRE(K) * CLDMAX(K)
          PGMLT(K) = PGMLT(K) * CLDMAX(K) 
  
          MNUCCC(K) = MNUCCC(K) * CFL3D_TEMP(K)
          PSACWS(K) = PSACWS(K) * CFL3D_TEMP(K)
          PSACWI(K) = PSACWI(k) * CFL3D_TEMP(K)
          QMULTS(K) = QMULTS(K) * CFL3D_TEMP(K)
          QMULTG(K) = QMULTG(K) * CFL3D_TEMP(K)
          PSACWG(K) = PSACWG(K) * CFL3D_TEMP(K)
          PGSACW(K) = PGSACW(K) * CFL3D_TEMP(K)

          PRD(K) = PRD(K) * CFI3D_TEMP(K)
          PRCI(K) = PRCI(K) * CFI3D_TEMP(K)
          PRAI(K) = PRAI(K) * CFI3D_TEMP(K) 
          QMULTR(K) = QMULTR(K) * CLDMAX(K)
          QMULTRG(K) = QMULTRG(K) * CLDMAX(K)
          MNUCCD(K) = MNUCCD(K) * CFI3D_TEMP(K)
          PRACI(K) = PRACI(K) * CFI3D_TEMP(K)
          PRACIS(K) = PRACIS(K) * CFI3D_TEMP(K)
          EPRD(K) = EPRD(K) * CFI3D_TEMP(K)

          MNUCCR(K) = MNUCCR(K) * CLDMAX(K)
          PIACR(K) = PIACR(K) * CFI3D_TEMP(K)
          PIACRS(K) = PIACRS(K) * CFI3D_TEMP(K) 
          PGRACS(K) = PGRACS(K) * CLDMAX(K)

          PSACR(K) = PSACR(K) * CLDMAX(K)
          MELTQI2QR(K) = MELTQI2QR(K) * CFI3D_TEMP(K)

!    Rain drop number process rates
          NPRC1(K) = NPRC1(K)* CFL3D_TEMP(K)
          NRAGG(K) = NRAGG(K) * CLDMAX(K)
          NPRACG(K) = NPRACG(K) * CLDMAX(K)
          NSUBR(K) = NSUBR(K) * CLDMAX(K)
          NSMLTR(K) = NSMLTR(K) * CLDMAX(K)
          NGMLTR(K) = NGMLTR(K) * CLDMAX(K)
          NPRACS(K) = NPRACS(K) * CLDMAX(K)
          NNUCCR(K) = NNUCCR(K) * CLDMAX(K)
          NIACR(K) = NIACR(K) * CFI3D_TEMP(K)
          NIACRS(K) = NIACRS(K) * CFI3D_TEMP(K) 
          NGRACS(K) = NGRACS(K) * CLDMAX(K)
#endif

#ifdef CLDFRAC
          MNUCCD(K) = MNUCCD0(K) * CFI3D_TEMP(K)
          PRD(K) = PRD0(K) * CFI3D_TEMP(K)
          EPRD(K) = EPRD0(K) * CFI3D_TEMP(K)
          PRDS(K) = PRDS0(K) * CLDMAX(K)
          EPRDS(K) = EPRDS0(K) * CLDMAX(K)
          PRDG(K) = PRDG0(K) * CLDMAX(K)  
          EPRDG(K) = EPRDG0(K) * CLDMAX(K)
          BERFD(K) = BERFD0(K) * CFL3D_TEMP(K)
          BERFDS(K) = BERFDS0(K) * CFL3D_TEMP(K)
          BERFDG(K) = BERFDG0(K) * CFL3D_TEMP(K)
#else
          MNUCCD(K) = MNUCCD0(K)
          PRD(K) = PRD0(K)
          EPRD(K) = EPRD0(K)
          PRDS(K) = PRDS0(K)
          EPRDS(K) = EPRDS0(K)
          PRDG(K) = PRDG0(K)
          EPRDG(K) = EPRDG0(K)
#endif

#ifdef CLUBB_CRM
! ADDITION BY UWM TO ENSURE THE POSITIVE DEFINITENESS OF VAPOR WATER MIXING RATIO
        CALL POSITIVE_QV_ADJ( QV3D(K), QC3D(K), QR3D(K), QI3D(K), &
                              QNI3D(K), QG3D(K), T3D(K) )
#endif /*CLUBB_CRM*/

#ifdef ECPP
! calculate relative humidity
!
         ! SATURATION VAPOR PRESSURE AND MIXING RATIO

            EVS(K) = POLYSVP(T3D(K),0)   ! PA
! MAKE SURE ICE SATURATION DOESN'T EXCEED WATER SAT. NEAR FREEZING
            QVS(K) = .622*EVS(K)/(PRES(K)-EVS(K))
            QVQVS(K) = QV3D(K)/QVS(K)
            RH3D(K)= min(1.0, QVQVS(K))

! WRF-CHEM, ADD TENDENCIES FOR C2PREC
!            C2PREC(K) = PRA(K)+PRC(K)
            C2PREC(K) = PRA(K)+PRC(K)+PSACWS(K)+QMULTS(K)+QMULTG(K)+PSACWG(K)+ &
                    PGSACW(K)+MNUCCC(K)+PSACWI(K)
            if(QC3D_INIT(K).gt.1.0e-10) then
              QSINK(K) = min(1.0, C2PREC(K)/QC3D_INIT(K))
            else 
              QSINK(K) = 0.0
            end if
#endif  /*ECPP*/

! hm 7/26/11, new output
           aut1d(k)=PRC(k)
           acc1d(k)=PRA(k)
           mlt1d(k)=-PSMLT(K)-PGMLT(K)+PRACS(K)+PRACG(K)+MELTQI2QR(K)
           evpr1d(k)=-PRE(K)-EVPMS(K)-EVPMG(K)
           if (PCC(k).lt.0.) then
             evpc1d(k)=-PCC(k)
           else if (PCC(k).gt.0.) then
            con1d(k)=PCC(k)
           end if
           sub1d(k)=-EPRD(K)-EPRDS(K)-EPRDG(K)
           dep1d(k)=PRD(K)+PRDS(K)+MNUCCD(K)+PRDG(K)

#endif /*ECPP*/

        QC3DTEN(K) = (QC3D(K)-QC3D_INIT(K))/DT
        NC3DTEN(K) = (NC3D(K)-NC3D_INIT(K))/DT
        QI3DTEN(K) = (QI3D(K)-QI3D_INIT(K))/DT
        NI3DTEN(K) = (NI3D(K)-NI3D_INIT(K))/DT
        QR3DTEN(K) = (QR3D(K)-QR3D_INIT(K))/DT
        NR3DTEN(K) = (NR3D(K)-NR3D_INIT(K))/DT
        QNI3DTEN(K) = (QNI3D(K)-QNI3D_INIT(K))/DT
        NS3DTEN(K) = (NS3D(K)-NS3D_INIT(K))/DT
        QG3DTEN(K) = (QG3D(K)-QG3D_INIT(K))/DT
        NG3DTEN(K) = (NG3D(K)-NG3D_INIT(K))/DT
        QV3DTEN(K) = (QV3D(K)-QV3D_INIT(K))/DT
        T3DTEN(K) = (T3D(K)-T3D_TEMP(K))/DT
       END DO  ! END K LOOP
 
#ifdef CLDFRAC
! check to ensure water conservations
      TOTAL_WATER_BEFORE = 0.0
      TOTAL_WATER_AFTER =0.0
      QTOT_BEFORE = 0.0
      QTOT_AFTER = 0.0
      QTOT_TEND = 0.0
      DO K=KMS, KME
          TOTAL_WATER_BEFORE = TOTAL_WATER_BEFORE+(QV3D_INIT(K)+QI3D_INIT(K)+QC3D_INIT(K)+QR3D_INIT(K)+QNI3D_INIT(K)+QG3D_INIT(K))*RHO(K)*DZQ(K)
          TOTAL_WATER_AFTER = TOTAL_WATER_AFTER+(QV3D(K)+QI3D(K)+QC3D(K)+QR3D(K)+QNI3D(K)+QG3D(K))*RHO(K)*DZQ(K)
          QTOT_AFTER(1) = QTOT_AFTER(1) + QV3D(K)*RHO(K)*DZQ(K)
          QTOT_AFTER(2) = QTOT_AFTER(2) + QI3D(K)*RHO(K)*DZQ(K)
          QTOT_AFTER(3) = QTOT_AFTER(3) + QC3D(K)*RHO(K)*DZQ(K)
          QTOT_AFTER(4) = QTOT_AFTER(4) + QR3D(K)*RHO(K)*DZQ(K)
          QTOT_AFTER(5) = QTOT_AFTER(5) + QNI3D(K)*RHO(K)*DZQ(K)
          QTOT_AFTER(6) = QTOT_AFTER(6) + QG3D(K)*RHO(K)*DZQ(K)
          QTOT_TEND(K) = QV3DTEN0(K)*RHO(K)*DZQ(K)*DT
          QTOT_BEFORE(1) = QTOT_BEFORE(1) + QV3D_INIT(K)*RHO(K)*DZQ(K)
          QTOT_BEFORE(2) = QTOT_BEFORE(2) + QI3D_INIT(K)*RHO(K)*DZQ(K)
          QTOT_BEFORE(3) = QTOT_BEFORE(3) + QC3D_INIT(K)*RHO(K)*DZQ(K)
          QTOT_BEFORE(4) = QTOT_BEFORE(4) + QR3D_INIT(K)*RHO(K)*DZQ(K)
          QTOT_BEFORE(5) = QTOT_BEFORE(5) + QNI3D_INIT(K)*RHO(K)*DZQ(K)
          QTOT_BEFORE(6) = QTOT_BEFORE(6) + QG3D_INIT(K)*RHO(K)*DZQ(K)
      END DO
      TOTAL_WATER_AFTER =TOTAL_WATER_AFTER+PRECRT
      if(abs(TOTAL_WATER_AFTER-TOTAL_WATER_BEFORE)/TOTAL_WATER_BEFORE.gt.1.0e-6) then
        write(901, *) 'water conservation issue in the morrison microphysics', TOTAL_WATER_BEFORE, TOTAL_WATER_AFTER, PRECRT, &
           abs(TOTAL_WATER_AFTER-TOTAL_WATER_BEFORE)/TOTAL_WATER_BEFORE
        write(901, *) 'water speices after', QTOT_AFTER(1:6)
        write(901, *) 'water species before', QTOT_BEFORE(1:6)
        write(901, *) 'Qtot_TEnd', QTOT_TEND(:)
!        stop
      end if


! Check to ensure energy conservation
      TOTAL_ENERGY_BEFORE = 0.0
      TOTAL_ENERGY_AFTER =0.0
      QTOTE_AFTER=0.
      QTOTE_BEFORE=0.
      DO K=KMS, KME 
        TOTAL_ENERGY_BEFORE = TOTAL_ENERGY_BEFORE + (T3D_TEMP(K)*CPM(K)-(QC3D_INIT(K)+QR3D_INIT(K))*XXLV(K)  &
                              -(QI3D_INIT(K)+QNI3D_INIT(K)+QG3D_INIT(K))*XXLS(K)) * RHO(K)*DZQ(K)
        TOTAL_ENERGY_AFTER = TOTAL_ENERGY_AFTER + (T3D(K)*CPM(K)-(QC3D(K)+QR3D(K))*XXLV(K)  &
                              -(QI3D(K)+QNI3D(K)+QG3D(K))*XXLS(K)) * RHO(K)*DZQ(K)
        QTOTE_BEFORE(1) = QTOTE_BEFORE(1) + T3D_TEMP(K)*CPM(K)*RHO(K)*DZQ(K)
        QTOTE_BEFORE(2) = QTOTE_BEFORE(2) + (-(QC3D_INIT(K)+QR3D_INIT(K))*XXLV(K)  &
                              -(QI3D_INIT(K)+QNI3D_INIT(K)+QG3D_INIT(K))*XXLS(K)) * RHO(K)*DZQ(K) 
        QTOTE_AFTER(1) = QTOTE_AFTER(1) + T3D(K)*CPM(K)*RHO(K)*DZQ(K)
        QTOTE_AFTER(2) = QTOTE_AFTER(2) + (-(QC3D(K)+QR3D(K))*XXLV(K)  &
                             -(QI3D(K)+QNI3D(K)+QG3D(K))*XXLS(K)) * RHO(K)*DZQ(K)
        QTOTE_AFTER(3) = QTOTE_AFTER(3) + T3DTEN(K)*DT*CPM(K)*RHO(K)*DZQ(K)
        QTOTE_TEND(K) = T3DTEN(K)*DT*CPM(K)*RHO(K)*DZQ(K)
        QTOTE_TEND2(K) = (-(QC3D(K)+QR3D(K))*XXLV(K)  &
                             -(QI3D(K)+QNI3D(K)+QG3D(K))*XXLS(K)) * RHO(K)*DZQ(K) - &
                         (-(QC3D_INIT(K)+QR3D_INIT(K))*XXLV(K)  &
                              -(QI3D_INIT(K)+QNI3D_INIT(K)+QG3D_INIT(K))*XXLS(K)) * RHO(K)*DZQ(K)
      END DO
      TOTAL_ENERGY_AFTER = TOTAL_ENERGY_AFTER-(PRECRT-SNOWRT)*XXLV(KMS)-SNOWRT*XXLS(KMS)
      if(abs(TOTAL_ENERGY_AFTER-TOTAL_ENERGY_BEFORE)/TOTAL_ENERGY_BEFORE.gt.1.0e-6) then
        write(0, *) 'energy conservation issue in the morrison microphysics', TOTAL_ENERGY_BEFORE, TOTAL_ENERGY_AFTER, &
               -(PRECRT-SNOWRT)*XXLV(KMS), -SNOWRT*XXLS(KMS),  &
           abs(TOTAL_ENERGY_AFTER-TOTAL_ENERGY_BEFORE)/TOTAL_ENERGY_BEFORE
        write(0, *) 'energy after', QTOTE_AFTER(1:2), QTOTE_AFTER(1)+QTOTE_AFTER(2), QTOTE_AFTER(3)
        write(0, *) 'energy before', QTOTE_BEFORE(1:2), QTOTE_BEFORE(1)+QTOTE_BEFORE(2)
      end if
#endif

! ALL DONE !!!!!!!!!!!
      RETURN
      END  SUBROUTINE M2005MICRO_GRAUPEL

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      REAL FUNCTION POLYSVP (T,TYPE)

!-------------------------------------------

!  COMPUTE SATURATION VAPOR PRESSURE

!  POLYSVP RETURNED IN UNITS OF PA.
!  T IS INPUT IN UNITS OF K.
!  TYPE REFERS TO SATURATION WITH RESPECT TO LIQUID (0) OR ICE (1)

      IMPLICIT NONE

      REAL DUM
      REAL T
      INTEGER TYPE

! REPLACE GOFF-GRATCH WITH FASTER FORMULATION FROM FLATAU ET AL. 1992, TABLE 4 (RIGHT-HAND COLUMN)

! ice
      real a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i 
      data a0i,a1i,a2i,a3i,a4i,a5i,a6i,a7i,a8i /&
	6.11147274, 0.503160820, 0.188439774e-1, &
        0.420895665e-3, 0.615021634e-5,0.602588177e-7, &
        0.385852041e-9, 0.146898966e-11, 0.252751365e-14/	

! liquid
      real a0,a1,a2,a3,a4,a5,a6,a7,a8 

! V1.7
      data a0,a1,a2,a3,a4,a5,a6,a7,a8 /&
	6.11239921, 0.443987641, 0.142986287e-1, &
        0.264847430e-3, 0.302950461e-5, 0.206739458e-7, &
        0.640689451e-10,-0.952447341e-13,-0.976195544e-15/
      real dt

! ICE

      IF (TYPE.EQ.1) THEN

!         POLYSVP = 10.**(-9.09718*(273.16/T-1.)-3.56654*                &
!          LOG10(273.16/T)+0.876793*(1.-T/273.16)+						&
!          LOG10(6.1071))*100.


      dt = max(-80.,t-273.16)
      polysvp = a0i + dt*(a1i+dt*(a2i+dt*(a3i+dt*(a4i+dt*(a5i+dt*(a6i+dt*(a7i+a8i*dt))))))) 
      polysvp = polysvp*100.

      END IF

! LIQUID

      IF (TYPE.EQ.0) THEN

       dt = max(-80.,t-273.16)
       polysvp = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
       polysvp = polysvp*100.

!         POLYSVP = 10.**(-7.90298*(373.16/T-1.)+                        &
!             5.02808*LOG10(373.16/T)-									&
!             1.3816E-7*(10**(11.344*(1.-T/373.16))-1.)+				&
!             8.1328E-3*(10**(-3.49149*(373.16/T-1.))-1.)+				&
!             LOG10(1013.246))*100.

         END IF


      END FUNCTION POLYSVP

!------------------------------------------------------------------------------

      REAL FUNCTION GAMMA(X)
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
!*******************************************************************
!*******************************************************************
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
!----------------------------------------------------------------------
      implicit none
      INTEGER I,N
      LOGICAL PARITY
      REAL                                                          &
          CONV,EPS,FACT,HALF,ONE,RES,SUM,TWELVE,                    &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      REAL, DIMENSION(7) :: C
      REAL, DIMENSION(8) :: P
      REAL, DIMENSION(8) :: Q
!----------------------------------------------------------------------
!  MATHEMATICAL CONSTANTS
!----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/


!----------------------------------------------------------------------
!  MACHINE DEPENDENT PARAMETERS
!----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,XINF/3.4E38/
!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,  &
             -3.79804256470945635097577E+2,6.29331155312818442661052E+2,  &
             8.66966202790413211295064E+2,-3.14512729688483675254357E+4,  &
             -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,  &
             -1.01515636749021914166146E+3,-3.10777167157231109440444E+3, &
              2.25381184209801510330112E+4,4.75584627752788110767815E+3,  &
            -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,                      &
           -5.952379913043012E-04,7.93650793500350248E-04,				   &
           -2.777777777777681622553E-03,8.333333333333333331554247E-02,	   &
            5.7083835261E-03/
!----------------------------------------------------------------------
!  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
!----------------------------------------------------------------------
      CONV(I) = REAL(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
        END DO
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
          DO I=1,N
            RES=RES*Y
            Y=Y+ONE
          END DO
        ENDIF
      ELSE
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO I=1,6
            SUM=SUM/YSQ+C(I)
          END DO
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
      RETURN
! ---------- LAST LINE OF GAMMA ----------
      END FUNCTION GAMMA


      REAL FUNCTION DERF1(X)
      IMPLICIT NONE
      REAL X
      REAL, DIMENSION(0 : 64) :: A, B
      REAL W,T,Y
      INTEGER K,I
      DATA A/                                                 &
         0.00000000005958930743E0, -0.00000000113739022964E0, &
         0.00000001466005199839E0, -0.00000016350354461960E0, &
         0.00000164610044809620E0, -0.00001492559551950604E0, &
         0.00012055331122299265E0, -0.00085483269811296660E0, &
         0.00522397762482322257E0, -0.02686617064507733420E0, &
         0.11283791670954881569E0, -0.37612638903183748117E0, &
         1.12837916709551257377E0,	                          &
         0.00000000002372510631E0, -0.00000000045493253732E0, &
         0.00000000590362766598E0, -0.00000006642090827576E0, &
         0.00000067595634268133E0, -0.00000621188515924000E0, &
         0.00005103883009709690E0, -0.00037015410692956173E0, &
         0.00233307631218880978E0, -0.01254988477182192210E0, &
         0.05657061146827041994E0, -0.21379664776456006580E0, &
         0.84270079294971486929E0,							  &
         0.00000000000949905026E0, -0.00000000018310229805E0, &
         0.00000000239463074000E0, -0.00000002721444369609E0, &
         0.00000028045522331686E0, -0.00000261830022482897E0, &
         0.00002195455056768781E0, -0.00016358986921372656E0, &
         0.00107052153564110318E0, -0.00608284718113590151E0, &
         0.02986978465246258244E0, -0.13055593046562267625E0, &
         0.67493323603965504676E0, 							  &
         0.00000000000382722073E0, -0.00000000007421598602E0, &
         0.00000000097930574080E0, -0.00000001126008898854E0, &
         0.00000011775134830784E0, -0.00000111992758382650E0, &
         0.00000962023443095201E0, -0.00007404402135070773E0, &
         0.00050689993654144881E0, -0.00307553051439272889E0, &
         0.01668977892553165586E0, -0.08548534594781312114E0, &
         0.56909076642393639985E0,							  &
         0.00000000000155296588E0, -0.00000000003032205868E0, &
         0.00000000040424830707E0, -0.00000000471135111493E0, &
         0.00000005011915876293E0, -0.00000048722516178974E0, &
         0.00000430683284629395E0, -0.00003445026145385764E0, &
         0.00024879276133931664E0, -0.00162940941748079288E0, &
         0.00988786373932350462E0, -0.05962426839442303805E0, &
         0.49766113250947636708E0 /
      DATA (B(I), I = 0, 12) /                                  &
         -0.00000000029734388465E0,  0.00000000269776334046E0, 	&
         -0.00000000640788827665E0, -0.00000001667820132100E0,  &
         -0.00000021854388148686E0,  0.00000266246030457984E0, 	&
          0.00001612722157047886E0, -0.00025616361025506629E0, 	&
          0.00015380842432375365E0,  0.00815533022524927908E0, 	&
         -0.01402283663896319337E0, -0.19746892495383021487E0,  &
          0.71511720328842845913E0 /
      DATA (B(I), I = 13, 25) /                                 &
         -0.00000000001951073787E0, -0.00000000032302692214E0,  &
          0.00000000522461866919E0,  0.00000000342940918551E0, 	&
         -0.00000035772874310272E0,  0.00000019999935792654E0, 	&
          0.00002687044575042908E0, -0.00011843240273775776E0, 	&
         -0.00080991728956032271E0,  0.00661062970502241174E0, 	&
          0.00909530922354827295E0, -0.20160072778491013140E0, 	&
          0.51169696718727644908E0 /
      DATA (B(I), I = 26, 38) /                                 &
         0.00000000003147682272E0, -0.00000000048465972408E0,   &
         0.00000000063675740242E0,  0.00000003377623323271E0, 	&
        -0.00000015451139637086E0, -0.00000203340624738438E0, 	&
         0.00001947204525295057E0,  0.00002854147231653228E0, 	&
        -0.00101565063152200272E0,  0.00271187003520095655E0, 	&
         0.02328095035422810727E0, -0.16725021123116877197E0, 	&
         0.32490054966649436974E0 /
      DATA (B(I), I = 39, 51) /                                 &
         0.00000000002319363370E0, -0.00000000006303206648E0,   &
        -0.00000000264888267434E0,  0.00000002050708040581E0, 	&
         0.00000011371857327578E0, -0.00000211211337219663E0, 	&
         0.00000368797328322935E0,  0.00009823686253424796E0, 	&
        -0.00065860243990455368E0, -0.00075285814895230877E0, 	&
         0.02585434424202960464E0, -0.11637092784486193258E0, 	&
         0.18267336775296612024E0 /
      DATA (B(I), I = 52, 64) /                                 &
        -0.00000000000367789363E0,  0.00000000020876046746E0, 	&
        -0.00000000193319027226E0, -0.00000000435953392472E0, 	&
         0.00000018006992266137E0, -0.00000078441223763969E0, 	&
        -0.00000675407647949153E0,  0.00008428418334440096E0, 	&
        -0.00017604388937031815E0, -0.00239729611435071610E0, 	&
         0.02064129023876022970E0, -0.06905562880005864105E0,   &
         0.09084526782065478489E0 /
      W = ABS(X)
      IF (W .LT. 2.2D0) THEN
          T = W * W
          K = INT(T)
          T = T - K
          K = K * 13
          Y = ((((((((((((A(K) * T + A(K + 1)) * T +              &
              A(K + 2)) * T + A(K + 3)) * T + A(K + 4)) * T +     &
              A(K + 5)) * T + A(K + 6)) * T + A(K + 7)) * T +     &
              A(K + 8)) * T + A(K + 9)) * T + A(K + 10)) * T + 	  &
              A(K + 11)) * T + A(K + 12)) * W
      ELSE IF (W .LT. 6.9D0) THEN
          K = INT(W)
          T = W - K
          K = 13 * (K - 2)
          Y = (((((((((((B(K) * T + B(K + 1)) * T +               &
              B(K + 2)) * T + B(K + 3)) * T + B(K + 4)) * T + 	  &
              B(K + 5)) * T + B(K + 6)) * T + B(K + 7)) * T + 	  &
              B(K + 8)) * T + B(K + 9)) * T + B(K + 10)) * T + 	  &
              B(K + 11)) * T + B(K + 12)
          Y = Y * Y
          Y = Y * Y
          Y = Y * Y
          Y = 1 - Y * Y
      ELSE
          Y = 1
      END IF
      IF (X .LT. 0) Y = -Y
      DERF1 = Y
      END FUNCTION DERF1

!+---+-----------------------------------------------------------------+
!
      subroutine radar_init

      IMPLICIT NONE
      INTEGER:: n
      PI5 = PI*PI*PI*PI*PI
      lamda4 = lamda_radar*lamda_radar*lamda_radar*lamda_radar
      m_w_0 = m_complex_water_ray (lamda_radar, 0.0d0)
      m_i_0 = m_complex_ice_maetzler (lamda_radar, 0.0d0)
      K_w = (ABS( (m_w_0*m_w_0 - 1.0) /(m_w_0*m_w_0 + 2.0) ))**2

      do n = 1, nbins+1
         simpson(n) = 0.0d0
      enddo
      do n = 1, nbins-1, 2
         simpson(n) = simpson(n) + basis(1)
         simpson(n+1) = simpson(n+1) + basis(2)
         simpson(n+2) = simpson(n+2) + basis(3)
      enddo

      do n = 1, slen
         mixingrulestring_s(n:n) = char(0)
         matrixstring_s(n:n) = char(0)
         inclusionstring_s(n:n) = char(0)
         hoststring_s(n:n) = char(0)
         hostmatrixstring_s(n:n) = char(0)
         hostinclusionstring_s(n:n) = char(0)
         mixingrulestring_g(n:n) = char(0)
         matrixstring_g(n:n) = char(0)
         inclusionstring_g(n:n) = char(0)
         hoststring_g(n:n) = char(0)
         hostmatrixstring_g(n:n) = char(0)
         hostinclusionstring_g(n:n) = char(0)
      enddo

      mixingrulestring_s = 'maxwellgarnett'
      hoststring_s = 'air'
      matrixstring_s = 'water'
      inclusionstring_s = 'spheroidal'
      hostmatrixstring_s = 'icewater'
      hostinclusionstring_s = 'spheroidal'

      mixingrulestring_g = 'maxwellgarnett'
      hoststring_g = 'air'
      matrixstring_g = 'water'
      inclusionstring_g = 'spheroidal'
      hostmatrixstring_g = 'icewater'
      hostinclusionstring_g = 'spheroidal'

      end subroutine radar_init
!+---+-----------------------------------------------------------------+

      COMPLEX*16 FUNCTION m_complex_water_ray(lambda,T)

!      Complex refractive Index of Water as function of Temperature T
!      [deg C] and radar wavelength lambda [m]; valid for
!      lambda in [0.001,1.0] m; T in [-10.0,30.0] deg C
!      after Ray (1972)

      IMPLICIT NONE
      REAL(8), INTENT(IN):: T,lambda
      REAL(8):: epsinf,epss,epsr,epsi
      REAL(8):: alpha,lambdas,sigma,nenner
      COMPLEX*16, PARAMETER:: i = (0d0,1d0)

      epsinf  = 5.27137d0 + 0.02164740d0 * T - 0.00131198d0 * T*T
      epss    = 78.54d+0 * (1.0 - 4.579d-3 * (T - 25.0)                 &
              + 1.190d-5 * (T - 25.0)*(T - 25.0)                        &
              - 2.800d-8 * (T - 25.0)*(T - 25.0)*(T - 25.0))
      alpha   = -16.8129d0/(T+273.16) + 0.0609265d0
      lambdas = 0.00033836d0 * exp(2513.98d0/(T+273.16)) * 1e-2

      nenner = 1.d0+2.d0*(lambdas/lambda)**(1d0-alpha)*sin(alpha*PI*0.5) &
             + (lambdas/lambda)**(2d0-2d0*alpha)
      epsr = epsinf + ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)   &
           * sin(alpha*PI*0.5)+1d0)) / nenner
      epsi = ((epss-epsinf) * ((lambdas/lambda)**(1d0-alpha)            &
           * cos(alpha*PI*0.5)+0d0)) / nenner                           &
           + lambda*1.25664/1.88496
      
      m_complex_water_ray = SQRT(CMPLX(epsr,-epsi))
      
      END FUNCTION m_complex_water_ray

!+---+-----------------------------------------------------------------+
      
      COMPLEX*16 FUNCTION m_complex_ice_maetzler(lambda,T)
      
!      complex refractive index of ice as function of Temperature T
!      [deg C] and radar wavelength lambda [m]; valid for
!      lambda in [0.0001,30] m; T in [-250.0,0.0] C
!      Original comment from the Matlab-routine of Prof. Maetzler:
!      Function for calculating the relative permittivity of pure ice in
!      the microwave region, according to C. Maetzler, "Microwave
!      properties of ice and snow", in B. Schmitt et al. (eds.) Solar
!      System Ices, Astrophys. and Space Sci. Library, Vol. 227, Kluwer
!      Academic Publishers, Dordrecht, pp. 241-257 (1998). Input:
!      TK = temperature (K), range 20 to 273.15
!      f = frequency in GHz, range 0.01 to 3000
         
      IMPLICIT NONE
      REAL(8), INTENT(IN):: T,lambda
      REAL(8):: f,c,TK,B1,B2,b,deltabeta,betam,beta,theta,alfa

      c = 2.99d8
      TK = T + 273.16
      f = c / lambda * 1d-9

      B1 = 0.0207
      B2 = 1.16d-11
      b = 335.0d0
      deltabeta = EXP(-10.02 + 0.0364*(TK-273.16))
      betam = (B1/TK) * ( EXP(b/TK) / ((EXP(b/TK)-1)**2) ) + B2*f*f
      beta = betam + deltabeta
      theta = 300. / TK - 1.
      alfa = (0.00504d0 + 0.0062d0*theta) * EXP(-22.1d0*theta)
      m_complex_ice_maetzler = 3.1884 + 9.1e-4*(TK-273.16)
      m_complex_ice_maetzler = m_complex_ice_maetzler                   &
                             + CMPLX(0.0d0, (alfa/f + beta*f)) 
      m_complex_ice_maetzler = SQRT(CONJG(m_complex_ice_maetzler))
      
      END FUNCTION m_complex_ice_maetzler
!+---+-----------------------------------------------------------------+

      subroutine rayleigh_soak_wetgraupel (x_g, a_geo, b_geo, fmelt,    &
                     meltratio_outside, m_w, m_i, lambda, C_back,       &
                     mixingrule,matrix,inclusion,                       &
                     host,hostmatrix,hostinclusion)

      IMPLICIT NONE

      REAL(8), INTENT(in):: x_g, a_geo, b_geo, fmelt, lambda,  &
                                     meltratio_outside
      REAL(8), INTENT(out):: C_back
      COMPLEX*16, INTENT(in):: m_w, m_i
      CHARACTER(len=*), INTENT(in):: mixingrule, matrix, inclusion,     &
                                     host, hostmatrix, hostinclusion

      COMPLEX*16:: m_core, m_air
      REAL(8):: D_large, D_g, rhog, x_w, xw_a, fm, fmgrenz,    &
                         volg, vg, volair, volice, volwater,            &
                         meltratio_outside_grenz, mra
      INTEGER:: error
      real :: rho_i, rho_w

      rho_i = 900.
      rho_w = 1000.


!     refractive index of air:
      m_air = (1.0d0,0.0d0)

!     Limiting the degree of melting --- for safety: 
      fm = DMAX1(DMIN1(fmelt, 1.0d0), 0.0d0)
!     Limiting the ratio of (melting on outside)/(melting on inside):
      mra = DMAX1(DMIN1(meltratio_outside, 1.0d0), 0.0d0)

!    ! The relative portion of meltwater melting at outside should increase
!    ! from the given input value (between 0 and 1)
!    ! to 1 as the degree of melting approaches 1,
!    ! so that the melting particle "converges" to a water drop.
!    ! Simplest assumption is linear:
      mra = mra + (1.0d0-mra)*fm

      x_w = x_g * fm

      D_g = a_geo * x_g**b_geo

      if (D_g .ge. 1d-12) then

       vg = PI/6. * D_g**3
       rhog = DMAX1(DMIN1(x_g / vg, DBLE(rho_i)), 10.0d0)
       vg = x_g / rhog
      
       meltratio_outside_grenz = 1.0d0 - rhog / rho_w

       if (mra .le. meltratio_outside_grenz) then
        !..In this case, it cannot happen that, during melting, all the
        !.. air inclusions within the ice particle get filled with
        !.. meltwater. This only happens at the end of all melting.
        volg = vg * (1.0d0 - mra * fm)
 
       else
        !..In this case, at some melting degree fm, all the air
        !.. inclusions get filled with meltwater.
        fmgrenz=(rho_i-rhog)/(mra*rho_i-rhog+rho_i*rhog/rho_w)

        if (fm .le. fmgrenz) then
         !.. not all air pockets are filled:
         volg = (1.0 - mra * fm) * vg
        else
         !..all air pockets are filled with meltwater, now the
         !.. entire ice sceleton melts homogeneously:
         volg = (x_g - x_w) / rho_i + x_w / rho_w
        endif

       endif

       D_large  = (6.0 / PI * volg) ** (1./3.)
       volice = (x_g - x_w) / (volg * rho_i)
       volwater = x_w / (rho_w * volg)
       volair = 1.0 - volice - volwater
      
       !..complex index of refraction for the ice-air-water mixture
       !.. of the particle:
       m_core = get_m_mix_nested (m_air, m_i, m_w, volair, volice,      &
                         volwater, mixingrule, host, matrix, inclusion, &
                         hostmatrix, hostinclusion, error)
       if (error .ne. 0) then
        C_back = 0.0d0
        return
       endif

       !..Rayleigh-backscattering coefficient of melting particle: 
       C_back = (ABS((m_core**2-1.0d0)/(m_core**2+2.0d0)))**2           &
                * PI5 * D_large**6 / lamda4

      else
       C_back = 0.0d0
      endif

      end subroutine rayleigh_soak_wetgraupel
!+---+-----------------------------------------------------------------+

      complex*16 function get_m_mix_nested (m_a, m_i, m_w, volair,      &
                     volice, volwater, mixingrule, host, matrix,        &
                     inclusion, hostmatrix, hostinclusion, cumulerror)

      IMPLICIT NONE

      REAL(8), INTENT(in):: volice, volair, volwater
      COMPLEX*16, INTENT(in):: m_a, m_i, m_w
      CHARACTER(len=*), INTENT(in):: mixingrule, host, matrix,          &
                     inclusion, hostmatrix, hostinclusion
      INTEGER, INTENT(out):: cumulerror

      REAL(8):: vol1, vol2
      COMPLEX*16:: mtmp
      INTEGER:: error

      !..Folded: ( (m1 + m2) + m3), where m1,m2,m3 could each be
      !.. air, ice, or water

      cumulerror = 0
      get_m_mix_nested = CMPLX(1.0d0,0.0d0)

      if (host .eq. 'air') then

       if (matrix .eq. 'air') then
        write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
        !bloss CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volice / MAX(volice+volwater,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, 0.0d0, vol1, vol2,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error
          
        if (hostmatrix .eq. 'air') then
         get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,              &
                         volair, (1.0d0-volair), 0.0d0, mixingrule,     &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'icewater') then
         get_m_mix_nested = get_m_mix (m_a, mtmp, 2.0*m_a,              &
                         volair, (1.0d0-volair), 0.0d0, mixingrule,     &
                         'ice', hostinclusion, error)
         cumulerror = cumulerror + error
        else
         write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                           hostmatrix
         !bloss CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'ice') then

       if (matrix .eq. 'ice') then
        write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
        !bloss CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volair / MAX(volair+volwater,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, vol1, 0.0d0, vol2,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'ice') then
         get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,              &
                         (1.0d0-volice), volice, 0.0d0, mixingrule,     &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'airwater') then
         get_m_mix_nested = get_m_mix (mtmp, m_i, 2.0*m_a,              &
                         (1.0d0-volice), volice, 0.0d0, mixingrule,     &
                         'air', hostinclusion, error)
         cumulerror = cumulerror + error          
        else
         write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',        &
                           hostmatrix
         !bloss CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'water') then

       if (matrix .eq. 'water') then
        write(mp_debug,*) 'GET_M_MIX_NESTED: bad matrix: ', matrix
        !bloss CALL wrf_debug(150, mp_debug)
        cumulerror = cumulerror + 1
       else
        vol1 = volair / MAX(volice+volair,1d-10)
        vol2 = 1.0d0 - vol1
        mtmp = get_m_mix (m_a, m_i, m_w, vol1, vol2, 0.0d0,             &
                         mixingrule, matrix, inclusion, error)
        cumulerror = cumulerror + error

        if (hostmatrix .eq. 'water') then
         get_m_mix_nested = get_m_mix (2.0d0*m_a, mtmp, m_w,            &
                         0.0d0, (1.0d0-volwater), volwater, mixingrule, &
                         hostmatrix, hostinclusion, error)
         cumulerror = cumulerror + error
        elseif (hostmatrix .eq. 'airice') then
         get_m_mix_nested = get_m_mix (2.0d0*m_a, mtmp, m_w,            &
                         0.0d0, (1.0d0-volwater), volwater, mixingrule, &
                         'ice', hostinclusion, error)
         cumulerror = cumulerror + error          
        else
         write(mp_debug,*) 'GET_M_MIX_NESTED: bad hostmatrix: ',         &
                           hostmatrix
         !bloss CALL wrf_debug(150, mp_debug)
         cumulerror = cumulerror + 1
        endif
       endif

      elseif (host .eq. 'none') then

       get_m_mix_nested = get_m_mix (m_a, m_i, m_w,                     &
                       volair, volice, volwater, mixingrule,            &
                       matrix, inclusion, error)
       cumulerror = cumulerror + error
        
      else
       write(mp_debug,*) 'GET_M_MIX_NESTED: unknown matrix: ', host
       !bloss CALL wrf_debug(150, mp_debug)
       cumulerror = cumulerror + 1
      endif

      IF (cumulerror .ne. 0) THEN
       write(mp_debug,*) 'GET_M_MIX_NESTED: error encountered'
       !bloss CALL wrf_debug(150, mp_debug)
       get_m_mix_nested = CMPLX(1.0d0,0.0d0)    
      endif

      end function get_m_mix_nested

!+---+-----------------------------------------------------------------+

      COMPLEX*16 FUNCTION get_m_mix (m_a, m_i, m_w, volair, volice,     &
                     volwater, mixingrule, matrix, inclusion, error)

      IMPLICIT NONE

      REAL(8), INTENT(in):: volice, volair, volwater
      COMPLEX*16, INTENT(in):: m_a, m_i, m_w
      CHARACTER(len=*), INTENT(in):: mixingrule, matrix, inclusion
      INTEGER, INTENT(out):: error

      error = 0
      get_m_mix = CMPLX(1.0d0,0.0d0)

      if (mixingrule .eq. 'maxwellgarnett') then
       if (matrix .eq. 'ice') then
        get_m_mix = m_complex_maxwellgarnett(volice, volair, volwater,  &
                           m_i, m_a, m_w, inclusion, error)
       elseif (matrix .eq. 'water') then
        get_m_mix = m_complex_maxwellgarnett(volwater, volair, volice,  &
                           m_w, m_a, m_i, inclusion, error)
       elseif (matrix .eq. 'air') then
        get_m_mix = m_complex_maxwellgarnett(volair, volwater, volice,  &
                           m_a, m_w, m_i, inclusion, error)
       else
        write(mp_debug,*) 'GET_M_MIX: unknown matrix: ', matrix
        !bloss CALL wrf_debug(150, mp_debug)
        error = 1
       endif

      else
       write(mp_debug,*) 'GET_M_MIX: unknown mixingrule: ', mixingrule
       !bloss CALL wrf_debug(150, mp_debug)
       error = 2
      endif

      if (error .ne. 0) then
       write(mp_debug,*) 'GET_M_MIX: error encountered'
       !bloss CALL wrf_debug(150, mp_debug)
      endif

      END FUNCTION get_m_mix

!+---+-----------------------------------------------------------------+

      COMPLEX*16 FUNCTION m_complex_maxwellgarnett(vol1, vol2, vol3,    &
                     m1, m2, m3, inclusion, error)

      IMPLICIT NONE

      COMPLEX*16 :: m1, m2, m3
      REAL(8) :: vol1, vol2, vol3
      CHARACTER(len=*) :: inclusion

      COMPLEX*16 :: beta2, beta3, m1t, m2t, m3t
      INTEGER, INTENT(out) :: error

      error = 0

      if (DABS(vol1+vol2+vol3-1.0d0) .gt. 1d-6) then
       write(mp_debug,*) 'M_COMPLEX_MAXWELLGARNETT: sum of the ',       &
              'partial volume fractions is not 1...ERROR'
       !bloss CALL wrf_debug(150, mp_debug)
       m_complex_maxwellgarnett=CMPLX(-999.99d0,-999.99d0)
       error = 1
       return
      endif

      m1t = m1**2
      m2t = m2**2
      m3t = m3**2

      if (inclusion .eq. 'spherical') then
       beta2 = 3.0d0*m1t/(m2t+2.0d0*m1t)
       beta3 = 3.0d0*m1t/(m3t+2.0d0*m1t)
      elseif (inclusion .eq. 'spheroidal') then
       beta2 = 2.0d0*m1t/(m2t-m1t) * (m2t/(m2t-m1t)*LOG(m2t/m1t)-1.0d0)
       beta3 = 2.0d0*m1t/(m3t-m1t) * (m3t/(m3t-m1t)*LOG(m3t/m1t)-1.0d0)
      else
       write(mp_debug,*) 'M_COMPLEX_MAXWELLGARNETT: ',                  &
                         'unknown inclusion: ', inclusion
       !bloss CALL wrf_debug(150, mp_debug)
       m_complex_maxwellgarnett=DCMPLX(-999.99d0,-999.99d0)
       error = 1
       return
      endif

      m_complex_maxwellgarnett = &
       SQRT(((1.0d0-vol2-vol3)*m1t + vol2*beta2*m2t + vol3*beta3*m3t) / &
       (1.0d0-vol2-vol3+vol2*beta2+vol3*beta3))

      END FUNCTION m_complex_maxwellgarnett

!+---+-----------------------------------------------------------------+
!..Compute radar reflectivity assuming 10 cm wavelength radar and using
!.. Rayleigh approximation.  Only complication is melted snow/graupel
!.. which we treat as water-coated ice spheres and use Uli Blahak's
!.. library of routines.  The meltwater fraction is simply the amount
!.. of frozen species remaining from what initially existed at the
!.. melting level interface.
!+---+-----------------------------------------------------------------+
      subroutine calc_refl10cm (qv1d, qr1d, qs1d, qg1d, t1d, p1d, dBZ,  &
                          kts, kte, ii, jj, nr1d, ns1d, ng1d)

      IMPLICIT NONE

!..Sub arguments
      INTEGER, INTENT(IN):: kts, kte, ii, jj
      REAL, DIMENSION(kts:kte), INTENT(IN)::                            &
                qv1d, qr1d, qs1d, qg1d, t1d, p1d, nr1d, ns1d, ng1d
      REAL, DIMENSION(kts:kte), INTENT(INOUT):: dBZ

!..Local variables
      REAL, DIMENSION(kts:kte):: temp, pres, qv, rho
      REAL, DIMENSION(kts:kte):: rr, rs, rg,rnr,rns,rng

      REAL(8), DIMENSION(kts:kte):: ilamr, ilamg, N0_r, N0_g,ilams,n0_s

      REAL, DIMENSION(kts:kte):: ze_rain, ze_snow, ze_graupel

      REAL(8):: lamg
      REAL(8):: fmelt_s, fmelt_g

      INTEGER:: i, k, k_0
      LOGICAL:: melti
      LOGICAL, DIMENSION(kts:kte):: L_qr, L_qs, L_qg

!..Single melting snow/graupel particle 70% meltwater on external sfc
      REAL(8), PARAMETER:: melt_outside_s = 0.7d0
      REAL(8), PARAMETER:: melt_outside_g = 0.7d0

      REAL(8):: cback, x, eta, f_d

! hm added parameter
      REAL R1,t_0,dumlams,dumlamr,dumlamg,dumn0s,dumn0r,dumn0g,ocms,obms,ocmg,obmg

      integer n

      R1 = 1.E-12
      t_0 = 273.15

!+---+

      do k = kts, kte
         dBZ(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kts, kte
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))
         if (qr1d(k) .gt. R1) then
            rr(k) = qr1d(k)*rho(k)
            L_qr(k) = .true.
         else
            rr(k) = R1
            L_qr(k) = .false.
         endif
         if (qs1d(k) .gt. R1) then
            rs(k) = qs1d(k)*rho(k)
            L_qs(k) = .true.
         else
            rs(k) = R1
            L_qs(k) = .false.
         endif
         if (qg1d(k) .gt. R1) then
            rg(k) = qg1d(k)*rho(k)
            L_qg(k) = .true.
         else
            rg(k) = R1
            L_qg(k) = .false.
         endif

! hm add number concentration
         if (nr1d(k) .gt. R1) then
            rnr(k) = nr1d(k)*rho(k)
         else
            rnr(k) = R1
         endif
         if (ns1d(k) .gt. R1) then
            rns(k) = ns1d(k)*rho(k)
         else
            rns(k) = R1
         endif
         if (ng1d(k) .gt. R1) then
            rng(k) = ng1d(k)*rho(k)
         else
            rng(k) = R1
         endif

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope, and useful moments for snow.
!+---+-----------------------------------------------------------------+
      do k = kts, kte

! compute moments for snow

! calculate slope and intercept parameter

      dumLAMS = (CONS1*rns(K)/rs(K))**(1./DS)
      dumN0S = rns(K)*dumLAMS/rho(k)

! CHECK FOR SLOPE to make sure min/max bounds are not exceeded

! ADJUST VARS

      IF (dumLAMS.LT.LAMMINS) THEN
      dumLAMS = LAMMINS
      dumN0S = dumLAMS**4*rs(K)/CONS1
      ELSE IF (dumLAMS.GT.LAMMAXS) THEN
      dumLAMS = LAMMAXS
      dumN0S = dumLAMS**4*rs(k)/CONS1
      end if

      ilams(k)=1./dumlams
      n0_s(k)=dumn0s

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept, slope values for graupel.
!+---+-----------------------------------------------------------------+

      do k = kte, kts, -1


! calculate slope and intercept parameter

      dumLAMg = (CONS2*rng(K)/rg(K))**(1./Dg)
      dumN0g = rng(K)*dumLAMg/rho(k)

! CHECK FOR SLOPE to make sure min/max bounds are not exceeded

! ADJUST VARS

      IF (dumLAMg.LT.LAMMINg) THEN
      dumLAMg = LAMMINg
      dumN0g = dumLAMg**4*rg(K)/CONS2
      ELSE IF (dumLAMg.GT.LAMMAXg) THEN
      dumLAMg = LAMMAXg
      dumN0g = dumLAMg**4*rg(k)/CONS2
      end if

      ilamg(k)=1./dumlamg
      n0_g(k)=dumn0g

      enddo

!+---+-----------------------------------------------------------------+
!..Calculate y-intercept & slope values for rain.
!+---+-----------------------------------------------------------------+

      do k = kte, kts, -1

! calculate slope and intercept parameter

      dumLAMr = (PI*RHOW*rnr(K)/rr(K))**(1./3.)
      dumN0r = rnr(K)*dumLAMr/rho(k)

! CHECK FOR SLOPE to make sure min/max bounds are not exceeded

! ADJUST VARS

      IF (dumLAMr.LT.LAMMINr) THEN
      dumLAMr = LAMMINr
      dumN0r = dumLAMr**4*rr(K)/(PI*RHOW)
      ELSE IF (dumLAMr.GT.LAMMAXr) THEN
      dumLAMr = LAMMAXr
      dumN0r = dumLAMr**4*rr(k)/(PI*RHOW)
      end if

      ilamr(k)=1./dumlamr
      n0_r(k)=dumn0r

      enddo

      melti = .false.
      k_0 = kts
      do k = kte-1, kts, -1
         if ( (temp(k).gt. T_0) .and. (rr(k).gt. 0.001e-3) &
                   .and. ((rs(k+1)+rg(k+1)).gt. 0.01e-3) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue

!+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

      do k = kts, kte
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*720.*ilamr(k)**7

         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                                 * (pi*rhosn/6./900.)*(pi*rhosn/6./900.) &
                                    * N0_s(k)*720.*ilams(k)**7
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
                                    * (pi*rhog/6./900.)* (pi*rhog/6./900.)        &
                                    * N0_g(k)*720.*ilamg(k)**7
      enddo

!+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (melti .and. k_0.ge.2) then
       do k = k_0-1, 1, -1

!..Reflectivity contributed by melting snow
          fmelt_s = DMIN1(1.0d0-rs(k)/rs(k_0), 1.0d0)
          if (fmelt_s.gt.0.01d0 .and. fmelt_s.lt.0.99d0 .and.           &
                         rs(k).gt.R1) then
           eta = 0.d0
           obms = 1./ds
           ocms = (1./(pi*rhosn/6.))**obms
           do n = 1, nbs
              x = pi*rhosn/6. * Dds(n)**3
              call rayleigh_soak_wetgraupel (x, DBLE(ocms), DBLE(obms), &
                    fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_s, matrixstring_s,          &
                    inclusionstring_s, hoststring_s,                    &
                    hostmatrixstring_s, hostinclusionstring_s)
              f_d = N0_s(k)* DEXP(-Dds(n)/ilams(k))
              eta = eta + f_d * CBACK * simpson(n) * dts(n)

           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif


!..Reflectivity contributed by melting graupel

          fmelt_g = DMIN1(1.0d0-rg(k)/rg(k_0), 1.0d0)
          if (fmelt_g.gt.0.01d0 .and. fmelt_g.lt.0.99d0 .and.           &
                         rg(k).gt.R1) then
           eta = 0.d0
           lamg = 1./ilamg(k)
           obmg = 1./dg
           ocmg = (1./(pi*rhog/6.))**obmg
           do n = 1, nbg
              x = pi*rhog/6. * Ddg(n)**3
              call rayleigh_soak_wetgraupel (x, DBLE(ocmg), DBLE(obmg), &
                    fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_g, matrixstring_g,          &
                    inclusionstring_g, hoststring_g,                    &
                    hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)* DEXP(-lamg*Ddg(n))
              eta = eta + f_d * CBACK * simpson(n) * dtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

       enddo
      endif

      do k = kte, kts, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo


      end subroutine calc_refl10cm
#ifdef CLUBB_CRM
!-------------------------------------------------------------------------------
  SUBROUTINE POSITIVE_QV_ADJ( QV, QC, QR, QI, &
                              QS, QG, T_IN_K )
! Description:
!   The following was produced by UW-Milwaukee to prevent vapor water mixing
!   ratio from becoming negative.  This is necessary in the event that a
!   process, e.g. depositional growth of ice, causes negative vapor. This
!   appears to happen in some circumstances due to the code that will set 
!   vapor to saturation w.r.t to liquid when we have subgrid scale cloud 
!   fraction greater than our 1% threshold.

! References:
!   None
!-------------------------------------------------------------------------------
    use constants_clubb, only: Lv, Ls, Cp ! Constant(s)

    IMPLICIT NONE

    ! Constant Parameters
    ! The value of epsilon was picked based on how small a 4 bytes float we can
    ! add to vapor without it being lost to catastophic round-off.  For an 8
    ! byte float a smaller value might be used -dschanen 5 Oct 2009.
    REAL, PARAMETER :: &
      EPS = 1.E-12 ! Small value of vapor [kg/kg]

    ! Input/Output Variables
    REAL, INTENT(INOUT) :: &
      QV,  & ! Vapor water mixing ratio       [kg/kg]
      QC,  & ! Cloud water mixing ratio       [kg/kg]
      QR,  & ! Rain water mixing ratio        [kg/kg]
      QI,  & ! Ice water mixing ratio         [kg/kg]
      QS,  & ! Snow water mixing ratio        [kg/kg]
      QG     ! Graupel water mixing ratio     [kg/kg]

    REAL, INTENT(INOUT) :: &
      T_IN_K ! Absolute Temperature     [K]

    ! Local Variables
    REAL :: &
      QT_COND_LIQ, & ! Total water in liquid phase      [kg/kg]
      QT_COND_ICE, & ! Total water in ice phase         [kg/kg]
      QT_TOTAL       ! Total water ice + liquid         [kg/kg]

    REAL :: &
      DELTA_QV, DELTA_QT_COND_LIQ, DELTA_QT_COND_ICE, REDUCE_COEF

    ! ---- Begin Code ----

    ! If vapor is greater than or equal to epsilon, then exit.
    IF ( QV >= EPS ) RETURN

!   PRINT *, "BEFORE", QV, QC, QR, QI, QS, QG, T_IN_K

    ! Determine total water
    QT_COND_LIQ = QC + QR

    QT_COND_ICE = 0.0
    ! Add ice if it is enabled
    IF ( ILIQ == 0 ) THEN
      QT_COND_ICE = QT_COND_ICE + QS + QI
    END IF

    ! Add graupel if it is enabled
    IF ( IGRAUP == 0 ) THEN
      QT_COND_ICE = QT_COND_ICE + QG
    END IF

    ! Total water mixing ratio = vapor + liquid + ice
    QT_TOTAL = QV + QT_COND_LIQ + QT_COND_ICE

    ! If the total water available at this altitude is too small,
    ! then we need to apply hole-filling globally instead.
    IF ( QT_TOTAL < 2 * EPS ) RETURN 

    ! Determine delta qv, the amount to change vapor water mixing ratio by.
    DELTA_QV = EPS - QV

    ! Set QV to the minimum value
    QV = EPS

    ! Reduce other variables according to the amount we've increased vapor by,
    ! in order to conserve total water.
    REDUCE_COEF = 1. - ( DELTA_QV / (QT_COND_LIQ + QT_COND_ICE) )

    ! Compute total change in warm-phase variables
    QC = QC * REDUCE_COEF
    QR = QR * REDUCE_COEF 

    DELTA_QT_COND_LIQ = QT_COND_LIQ - ( QC + QR )

    ! Compute total change in ice-phase variables

    DELTA_QT_COND_ICE = 0.0
    IF ( ILIQ == 0 ) THEN
      QI = QI * REDUCE_COEF
      QS = QS * REDUCE_COEF

      IF ( IGRAUP /= 0 ) THEN
        DELTA_QT_COND_ICE = QT_COND_ICE - ( QI + QS )
      END IF
    END IF

    IF ( IGRAUP == 0 ) THEN
      QG = QG * REDUCE_COEF

      DELTA_QT_COND_ICE = QT_COND_ICE - ( QI + QS + QG )
    END IF

    ! Adjust absolute temperature
    T_IN_K = T_IN_K - ( Lv / Cp * ( DELTA_QT_COND_LIQ ) ) &
                    - ( Ls / Cp * ( DELTA_QT_COND_ICE ) )

!   PRINT *, "AFTER", QV, QC, QR, QI, QS, QG, T_IN_K
    RETURN
  END SUBROUTINE POSITIVE_QV_ADJ
#endif /*CLUBB_CRM*/

END MODULE module_mp_GRAUPEL
