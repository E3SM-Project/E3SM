!------------------------------------------------------------------------------
!    'fjx_cmn_mod.f90'  for Solar/Cloud/Fast-J code v 7.6c (prather 06/2019)
!------------------------------------------------------------------------------

      MODULE FJX_CMN_MOD

      implicit none
      public

!------------------------------------------------------------------------------
      logical  LRRTMG, LCLIRAD, LGGLLNL  ! set in fjx_init_mod based on W_rrtmg
      integer  NWBIN, NSBIN     ! dimensions of readin spec - NWBIN can zero out strat-wavels
!-------------basic vertical grid----------these can be changed as needed, needed for dimensions
      integer, parameter :: &  !  these must e set to exact atmospheric dimensions
            LPAR = PLEV,      &  !  # layers in model (one layer is added to go to ZTOP+ZZHT)
            LWEPAR = PLEV !47    !  # layers that have clouds (LWEtPAR < LPAR)
      integer, parameter :: L_ = LPAR, L1_=L_+1, L2_ = L_+2
!  L_ = # of CTM layers, L1_=L_+1 = # of CTM layer edges (radii)
!  L2_ = L_+2 = total # of layer edges counting top (TAU=0)

!----------these below should be set by the calling program
      integer, parameter :: JVL_ = LPAR  ! vertical(levels) dim for J-values sent to CTM
      integer, parameter :: JVN_ = 101   ! max # of J-values
      integer, parameter :: AN_=25       ! max # FJX aerosols in layer (needs NDX for each)

!--------------------------------------------------------------------------
      ! JXL_: vertical(levels) dim for J-values computed within fast-JX
      integer, parameter ::  JXL_=100, JXL1_=JXL_+1
      ! JXL2_: 2*JXL_ + 2 = mx no.levels in basic FJX grid (mid-level)
      integer, parameter ::  JXL2_=2*JXL_+2

!  case 1 RRTMG super bins S_ , SX_(for tables) =27
      integer, parameter ::  WX_= 18    ! used for table dimensions
      integer, parameter ::  SX_= 27
      ! W_   = dim = no. of Fast-J Wavelength bins:  currenly only 18, TROP-ONLY is done by zeroing FL fluxes
      integer, parameter ::  W_=18
      ! S_   = dim = number of wavelength bins INLCUDING the Solar-J extensions (RRTMG value = 27)
      integer, parameter ::  S_=27

      integer, parameter ::  W_r = S_-W_  ! # of bins that is added on top of W_
!SJ! this is defined in fjx_init_mod      integer :: W_r
      integer, parameter ::  W_RRTMG = 0  !  = 82 for std RRTMG
      integer, parameter ::  W_CLIRAD = 0
      integer, parameter ::  W_LLNL = 0
      integer, parameter, dimension(27) :: NGC = &
          (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &      ! these are Cloud-J, no sub-bins
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
            1, 1, 1, 1, 1, 1, 1/)
! SJ          (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &   !these are use in Solar-J for RRTMG
! SJ           1, 1, 1, 1, 1, 1, 1, 5,10, 2, &
! SJ          10,10, 8, 8,12, 6,12/)
!
      ! X_   = dim = max no. of X-section data sets (input data)
      integer, parameter ::  X_=72
      ! A_   = dim = max no. of Aerosol Mie sets (input data) not including clouds and SSA
      integer, parameter ::  A_=40
      ! SSA_ & GGA_ = dim = no. of strat sulfate aerosol types (input data)
      integer, parameter ::  SSA_=18, GGA_=15
      ! C_   = dim = no. of cld-data sets (input data):  liquid-water,irregular-ice, hexagonal ice
      integer, parameter ::  C_=3
      ! CR_   = dim = no. of effective radii in each cld-data sets
      integer, parameter ::  CR_=6
      ! N_  = no. of levels in Mie scattering arrays = 2*(L_+1) + 1`+ 2*added-cld-layers
      integer, parameter ::  N_=601
      ! M_  = no. of Gauss points used, must = 4 in fast_JX (no option)
      integer, parameter ::  M_=4
      ! M2_ = 2*M_ = 8, replaces MFIT
      integer, parameter ::  M2_=2*M_

!-----------------------------------------------------------------------
      ! 4 Gauss pts = 8-stream
      real*8, DIMENSION(M_), parameter  ::  &
                         EMU = [.06943184420297d0, .33000947820757d0, &
                                .66999052179243d0, .93056815579703d0]
      real*8, DIMENSION(M_), parameter  :: &
                         WT  = [.17392742256873d0, .32607257743127d0, &
                                .32607257743127d0, .17392742256873d0]
!-----------------------------------------------------------------------
!     input to osa_sub_mod module   - need single precision
      real, DIMENSION(5) :: ANGLES
!-----------------------------------------------------------------------

! Radiation Field, Cloud Cover & Other fixed parameters
      ! MASFAC: Conversion factor for pressure to column density: [PJC: P(mbar) * MASFAC = molecules/cm2]
      real*8, parameter:: MASFAC = 100.d0*6.022d+23/(28.97d0*9.8d0*10.d0)
      ! HeatFac_: convert watt/m2 to K/day
      real*8, parameter:: HeatFac_ = 86400.d0*9.80616d0/1.00464d5

! these key parameters are read in at initialization and not locked in at compile
      ! ZZHT: scale height (cm) used above top of CTM ZHL(LPAR+1)
!      real*8, parameter   :: ZZHT = 5.d5
      ! RAD: Radius of Earth (cm)
!      real*8, parameter   :: RAD = 6375.d5
      ! ATAU: Factor increase in cloud optical depth (OD) from layer to next below
!      real*8, parameter   :: ATAU = 1.050d0
      ! ATAU0: Minimum cloud OD in uppermost inserted layer
!      real*8, parameter   :: ATAU0 = 0.005d0
      ! ATM0: Option for spherical corrections: 0=flat 1=sphr 2=refr 3=geom
!      integer, parameter  :: ATM0 = 3
      real*8 ZZHT, RAD, ATAU, ATAU0, CLDCOR
      integer ATM0, NRANDO, LNRG, CLDFLAG
      character*25, dimension(8), parameter :: TITCLD =  &
      ['clear sky - no clouds    ', &
       'avg cloud cover          ', &
       'avg cloud cover^3/2      ', &
       'ICAs - avg direct beam   ', &
       'ICAs - random N ICAs     ', &
       'QCAs - midpt of bins     ', &
       'QCAs - avg clouds in bins', &
       'ICAs - use all ICAs***   ']

!---- Variables in file 'FJX_spec.dat' (RD_XXX)
      ! WL: Centres of wavelength bins - 'effective wavelength'  (nm)
      real*8  WL(S_)
      ! WBIN: Boundaries of wavelength bins                  (microns)
      real*8  WBIN(S_+1)
      ! FL: Solar flux incident on top of atmosphere (cm-2.s-1)
      ! FW: Solar flux in W/m2
      ! FP: PAR quantum action spectrum
      real*8  FL(S_),FW(S_),FP(S_)
      ! QRAYL: Rayleigh parameters (effective cross-section) (cm2)
      real*8  QRAYL(S_)       ! ????? SJ had this set ant S_+1, it should not needed
      ! SJSUB:  intended for breakdown of the super-bins (1:27) into smaller sub-bins.
      real*8  SJSUB(S_,16)
      integer KDOKR(100), LDOKR(100) ! RRTMG =18+82 max set at 100
      integer NSJSUB(SX_)

      real*8  QO2(WX_,3)   ! QO2: O2 cross-sections
      real*8  QO3(WX_,3)   ! QO3: O3 cross-sections
      real*8  Q1D(WX_,3)   ! Q1D: O3 => O(1D) quantum yield

      ! QQQ: Supplied cross sections in each wavelength bin (cm2)
      real*8  QQQ(WX_,3,X_)
      ! TQQ: Temperature for supplied cross sections
      real*8  TQQ(3,X_)
      ! LQQ = 1, 2, or 3 to determine interpolation with T or P
      integer LQQ(X_)

      ! TITLEJX: Title (short & long) for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*6  TITLEJX(X_)
      CHARACTER*16 TITLEJL(X_)
      ! SQQ: Flag for supplied cross sections, from 'FJX_spec.dat'
      CHARACTER*1  SQQ(X_)

!---- Variables in file 'FJX_scat-aer.dat' (RD_MIE)
      ! TITLAA: Aerosol Mie Titles
      character*12  TITLAA(A_)
      ! QAA: Aerosol scattering phase functions
      real*8  QAA(5,A_)
      ! WAA: 5 Wavelengths for the supplied phase functions
      real*8  WAA(5,A_)
      ! PAA: Phase function: first 8 terms of expansion
      real*8  PAA(8,5,A_)
      ! RAA: Effective radius associated with aerosol type
      real*8  RAA(A_)
      ! SAA: Single scattering albedo
      real*8  SAA(5,A_)
      ! DAA: density (g/cm^3)
      real*8  DAA(A_)
      ! NAA: Number of categories for scattering phase functions
      integer NAA

!---- Variables in file 'FJX_scat-cld.dat' (RD_CLD)  ***major update for v75
      ! NCC: Number of categories for cloud scattering phase functions, MCC: no. of R_eff
      integer NCC,MCC
      ! TITLCC: Cloud type titles
      character*12  TITLCC(C_)
      ! RCC: Effective radius associated with cloud type:  should be indep of S-bin
      real*8  RCC(CR_,C_)
      ! GCC: Effective geometric cross section:  should be indep of S-bin
      real*8  GCC(CR_,C_)
      ! DCC: density (g/cm^3):  should be indep of S-bin, and eff radius
      real*8  DCC(C_)
      ! QCC: Cloud Q-ext
      real*8  QCC(SX_,CR_,C_)
      ! WCC: Wavelengths for supplied phase functions, should be std S-bins
      real*8  WCC(SX_,C_)
      ! SCC: Single scattering albedo
      real*8  SCC(SX_,CR_,C_)
      ! PCC: Phase function: first 8 terms of expansion
      real*8  PCC(8,SX_,CR_,C_)

!---- Variables in file 'FJX_scat-ssa.dat' (RD_SSA)
      ! NSS: Number of categories for Strat Sulf Aerosol scattering phase functions
      integer NSS
      ! TITLSS: Cloud type titles
      character*12  TITLSS(SSA_)
      ! RSS: Effective radius associated with cloud type
      real*8  RSS(SSA_)
      ! GSS: Effective geometric cross section
      real*8  GSS(SSA_)
      ! DSS: density (g/cm^3)
      real*8  DSS(SSA_)
      ! TSS: temperature (K)
      real*8  TSS(SSA_)
      ! WSS: weight percent sulfuric acid (%)
      real*8  WSS(SSA_)
      ! QSS: Q-ext      ----begin wavelength dependent quantities
      real*8  QSS(SX_,SSA_)
      ! SSS: Single scattering albedo
      real*8  SSS(SX_,SSA_)
      ! PSS: Phase function: first 8 terms of expansion
      real*8  PSS(8,SX_,SSA_)

!---- Variables in file 'FJX_scat-geo.dat' (RD_GEO)
      ! NGG: Number of categories for Strat Sulf Aerosol scattering phase functions
      integer NGG
      ! RGG: Effective radius associated with cloud type
      real*8  RGG(GGA_)
      ! DGG: density (g/cm^3)
      real*8  DGG(GGA_)
      ! QGG: Q-ext      ----begin wavelength dependent quantities
      real*8  QGG(SX_,GGA_)
      ! SGG: Single scattering albedo
      real*8  SGG(SX_,GGA_)
      ! PGG: Phase function: first 8 terms of expansion
      real*8  PGG(8,SX_,GGA_)

!---- Variables in file 'FJX_scat-UMa.dat' (RD_CLD)
      ! WMM: U Michigan aerosol wavelengths
      real*8  WMM(6)
      ! UMAER: U Michigan aerosol data sets
      real*8  UMAER(3,6,21,33)

!---- Variables in file 'atmos_std.dat' (RD_PROF) and 'atmos_h2och4.dat' (RD_TRPROF)
      integer, parameter ::  LREF=51   ! layer dim. in reference profiles
      integer, parameter ::  JREF=18   ! latitude dim. in reference profiles
!----- T and O3, H2O, CH4, reference profiles added underscore _ because TREF used in RRTMG_SW
      real*8, DIMENSION(LREF,JREF,12) :: T_REF, O_REF, H2O_REF, CH4_REF
      integer NJX,NW1,NW2,NS1,NS2

!----- Reference monthly zonal mean profiles for GEOMIP SSA 'atmos_geomip.dat'
      integer, parameter ::  LGREF=19   ! layer dim. in reference profiles
      real*8, dimension(64,19,12) :: R_GREF,X_GREF,A_GREF   ! R = Reff (microns), X = micro-g-H2SO4/kg-air
      real*8, dimension(64)       :: Y_GREF                 !  A = 4 pi R^2 = microns^2/cm^3
      real*8, dimension(19)       :: P_GREF

      real*8  JFACTA(JVN_)  ! multiplication factor for fast-JX calculated J
      integer JIND(JVN_)    ! index arrays that map Jvalue(j) onto rates
      integer NRATJ         ! number of Photolysis reactions in CTM chemistry, NRATJ <= JVN_
      character*6 JVMAP(JVN_) !label of J-value used to match w/FJX J's
      character*50 JLABEL(JVN_) ! label of J-value used in the chem model

! Cloud overlap parameters
!-----------------------------------------------------------------------
      integer, parameter :: CBIN_ = 10     ! # of quantized cloud fration bins
      integer, parameter :: ICA_ = 20000   ! Max # of indep colm atmospheres
      integer, parameter :: NQD_ = 4       ! # of cloud fraction bins (4)

      real*8,  parameter ::  CPI    = 3.141592653589793d0
      real*8,  parameter ::  C2PI   = 2.d0*CPI
      real*8,  parameter ::  CPI180 = CPI/180.d0
      real*8,  parameter ::  G0     = 9.80665d0
      real*8,  parameter ::  G100 = 100.d0/G0
!-------data to set up the random number sequence for use in cloud-JX
      integer, parameter :: NRAN_ = 10007  ! dimension for random number
      real*4   RAN4(NRAN_)      ! Random number set

!----extras for Cloud-J
      integer, parameter::  mxlay = L1_, ngptsw = S_+1

      END MODULE FJX_CMN_MOD
