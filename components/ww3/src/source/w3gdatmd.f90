!/ ------------------------------------------------------------------- /
      MODULE W3GDATMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  !           J. H. Alves             !
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         29-May-2009 |
!/                  +-----------------------------------+
!/
!/    24-Jun-2005 : Origination.                        ( version 3.07 )
!/    09-Nov-2005 : Remove soft boundary options.       ( version 3.08 )
!/    23-Jun-2006 : Add data for W3SLN1.                ( version 3.09 )
!/    18-Jul-2006 : Add input grids.                    ( version 3.10 )
!/    05-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    02-Feb-2007 : Add FLAGST.                         ( version 3.10 )
!/    14-Apr-2007 : Add Miche style limiter.            ( version 3.11 )
!/                  ( J. H. Alves )
!/    25-Apr-2007 : Adding Battjes-Janssen Sdb.         ( version 3.11 )
!/                  ( J. H. Alves )
!/    06-Aug-2007 : Fixing SLNP !/SEED bug.             ( version 3.13 )
!/    18-Sep-2007 : Adding WAM4 source terms.           ( version 3.13 )
!/                  ( F. Ardhuin )
!/    15-Apr-2008 : Clean up for distribution.          ( version 3.14 )
!/    27-Jun-2008 : Expand WAM4 variants namelist       ( version 3.14 )
!/                  ( F. Ardhuin )
!/    29-May-2009 : Preparing distribution version.     ( version 3.14 )
!/
!/    Copyright 2009 National Weather Service (NWS),
!/       National Oceanic and Atmospheric Administration.  All rights
!/       reserved.  WAVEWATCH III is a trademark of the NWS.
!/       No unauthorized use without permission.
!/
!  1. Purpose :
!
!     Define data structures to set up wave model grids and aliases
!     to use individual grids transparently. Also includes subroutines
!     to manage data structure and pointing to individual models.
!     Definition of grids and model set up.
!
!  2. Variables and types :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NGRIDS    Int.  Public   Number of grids, initialized at -1
!                               to check proper model initialization.
!      NAUXGR    Int.  Public   Auxiliary grids.
!      IGRID     Int.  Public   Selected spatial grid, init. at -1.
!      ISGRD     Int.  Public   Selected spectral grid, init. at -1.
!      IPARS     Int.  Public   Selected num. and ph. pars, init. at -1.
!      FLAGLL    L.P.  Public   Flag for LL grid (otherwise XY).
!      GRID      TYPE  Public   Data structure defining grid.
!      GRIDS     GRID  Public   Array of grids.
!      SGRD      TYPE  Public   Data structure defining spectral grid.
!      SGRDS     GRID  Public   Array of spectral grids.
!      MPAR      TYPE  Public   Data structure with all other model
!                               parameters.
!      MPARS     GRID  Public   Array of MPAR.
!     ----------------------------------------------------------------
!
!     All elements of GRID are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NX, NY    Int.  Public   Discrete dimensions of spatial grid.
!      NSEA(L)   Int.  Public   Number of sea points (local for MPP).
!      TRFLAG    Int.  Public   Flag for use of transparencies
!                                0: No sub-grid obstacles.
!                                1: Obstructions at cell boundaries.
!                                2: Obstructions at cell centers.
!                                3: Like 1 with continuous ice.
!                                4: Like 2 with continuous ice.
!      MAPSTA    I.A.  Public   Grid status map.
!      MAPST2    I.A.  Public   Second grid status map.
!      MAPxx     I.A.  Public   Storage grid maps.
!      SX,SY     Real  Public   Spatial grid increments.
!      X0,Y0     Real  Public   Lower left corner of spatial grid.
!      DTCFL     Real  Public   Maximum CFL time step X-Y propagation.
!      DTCFLI    Real  Public   Id. intra-spectral.
!      DTMAX     Real  Public   Maximum overall time step.
!      DTMIN     Real  Public   Minimum dynamic time step for source
!                               term integration.
!      DMIN      Real  Public   Minimum water depth.
!      CTMAX     Real  Public   Maximum CFL number for depth refr.
!      FICE0/N   Real  Public   Cut-off ice conc. for ice coverage.
!      PFMOVE    Real  Public   Tunable parameter in GSE correction
!                               for moving grids.
!      ZB        R.A.  Public   Bottom levels on storage grid.
!      CLAT(I)   R.A.  Public   (Inverse) cosine of latitude.
!      CLATS     R.A.  Public   Id.
!      CTHG0     R.A.  Public   Constant in great-circle refr. term.
!      TRNX/Y    R.A.  Public   Transparencies in X/Y for sub-grid
!      GINIT     Log.  Public   Flag identifying grid initialization.
!      GLOBAL    Log.  Public   Flag for global grid (closed parallels).
!      FLDRY     Log.  Public   Flag for 'dry' run (IO and data
!                               processing only).
!      FLCx      Log.  Public   Flags for prop. is different spaces.
!      FLSOU     Log.  Public   Flag for source term calcualtion.
!      FLAGST    L.A.  Public   Flag for source term computations
!                               for individual grid points.
!      GNAME     C*30  Public   Grid name.
!      FILEXT    C*10  Public   Extension of WAVEWATCH III file names
!                               default in 'ww3'.
!     ----------------------------------------------------------------
!
!     All elements of SGRD are aliased to pointers with the same
!     name. These pointers are defined as :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      NK        Int.  Public   Number of discrete wavenumbers.
!      NK2       Int.  Public   Extended wavenumber range.
!      NTH       Int.  Public   Number of discrete directions.
!      NSPEC     Int.  Public   Number of discrete spectral bins.
!      MAPxx     I.A.  Public   Spectral maps.
!      DTH       Real  Public   Directional increments (radians).
!      XFR       Real  Public   Frequency multiplication factor.
!      FR1       Real  Public   Lowest frequency                 (Hz)
!      FTE       Real  Public   Factor in tail integration energy.
!      FTF       Real  Public   Id. frequency.
!      FTWN      Real  Public   Id. wavenumber.
!      FTTR      Real  Public   Id. wave period.
!      FTWL      Real  Public   Id. wave length.
!      FACTIn    Real  Public   Factors for obtaining integer cut-off
!                               frequency.
!      FACHFx    Real  Public   Factor for tail.
!      TH        R.A   Public   Directions (radians).
!      ESIN      R.A   Public   Sine of discrete directions.
!      ECOS      R.A   Public   Cosine of discrete directions.
!      ES2, ESC, EC2
!                R.A   Public   Sine and cosine products
!      SIG       R.A   Public   Relative frequencies (invariant
!                                                     in grid). (rad)
!      SIG2      R.A   Public   Id. for full 2-D spectrum.
!      DSIP      R.A   Public   Frequency bandwidths (prop.)    (rad)
!      DSII      R.A   Public   Frequency bandwidths (int.)     (rad)
!      DDEN      R.A   Public   DSII * DTH * SIG (for integration
!                               based on energy)
!      DDEN2     R.A   Public   Idem, full spectrum.
!      SINIT     Log.  Public   Flag identifying grid initialization.
!     ----------------------------------------------------------------
!
!     The structure MPAR contains all other model parameters for
!     numerical methods and physical parameterizations. It contains
!     itself several structures as outlined below.
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      PINIT     Log.  Public   Flag identifying initialization.
!      NPARS     NPAR  Public   Numerical parameters,
!      PROPS     PROP  Public   Parameters propagatrion schemes.
!      SFLPS     SFLP  Public   Parameters for flux computation.
!      SLNPS     SLNP  Public   Parameters Sln.
!      SRCPS     SRCP  Public   Parameters Sin and Sds.
!      SNLPS     SNLP  Public   Parameters Snl.
!      SBTPS     SBTP  Public   Parameters Sbt.
!      SDBPS     SDBP  Public   Parameters Sdb.
!      STRPS     STRP  Public   Parameters Str.
!      SBSPS     SBSP  Public   Parameters Sbs.
!      SXXPS     SXXP  Public   Parameters Sxx.
!     ----------------------------------------------------------------
!
!     The structure NPAR contains numerical parameters and is aliased
!     as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      FACP      Real  Public   Constant in maximum par. change in
!                               dynamic integration scheme (depends
!                               upon Xp).
!      XREL      Real  Public   Id. relative change.
!      XFLT      Real  Public   Id. filter level.
!      FXFM      Real  Public   Constant for mean frequency in
!                               cut-off.                       (!/ST1)
!      FXPM      Real  Public   Id. PM.
!      XFT       Real  Public   Constant for cut-off freq.     (!/ST2)
!      XFC       Real  Public   Id.
!      FACSD     Real  Public   Constant in seeding algorithm.
!      FHMAX     Real  Public   Hs/depth ratio in limmiter    (!/MLIM)
!     ----------------------------------------------------------------
!
!     The structure PROP contains parameters for the propagation
!     schemes and is aliased as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      DTME      Real  Public   Swell age in disp. corr.      (!/PR2)
!      CLATMN    Real  Public   Id. minimum cosine of lat.    (!/PR2)
!
!      WDCG      Real  Public   Factors in width of av. Cg.   (!/PR3)
!      WDTH      Real  Public   Factors in width of av. Th.   (!/PR3)
!     ----------------------------------------------------------------
!
!     The structure SFLP contains parameters for the fluxes
!     and is aliased as above:
!     ----------------------------------------------------------------
!                                                            (!/FLX2)
!      NITTIN    Int.  Public   Number of itterations for drag calc.
!      CINXSI    Real  Public   Constant in parametric description
!                                                            (!/FLX3)
!      NITTIN    Int.  Public   Number of itterations for drag calc.
!      CAP_ID    Int   Public   Type of cap used.
!      CINXSI    Real  Public   Constant in parametric description
!      CD_MAX    Real  Public   Cap on Cd.
!     ----------------------------------------------------------------
!
!     The structure SLNP contains parameters for the linear input
!     source terms and is aliased as above:
!
!     ----------------------------------------------------------------
!                                                             (!/LN1)
!      SLNC1     Real  Public   Proportionality and other constants in
!                               input source term.
!      FSPM      Real  Public   Factor for fPM in filter.
!      FSHF      Real  Public   Factor for fh in filter.
!     ----------------------------------------------------------------
!
!     The structure SRCP contains parameters for the input and dis,
!     source terms and is aliased as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      WWNMEANPTAIL R  Public   Power of tail for WNMEAN calculation
!      SSTXFTFTAIL  R  Public   Tail factor for  WNMEAN calculation
!                                                             (!/ST1)
!      SINC1     Real  Public   Proportionality and other constants in
!                               input source term.
!      SDSC1     Real  Public   Combined constant in dissipation
!                               source term.
!                                                             (!/ST2)
!      ZWIND     Real  Public   Height at which the wind is defined
!                               of drag.
!      FSWELL    Real  Public   Reduction factor of negative input
!                               for swell.
!      SHSTAB, OFSTAB, CCNG, CCPS, FFNG, FFPS
!                Real  Public   Factors in effective wind speed.
!      CDSAn     Real  Public   Constants in high-freq. dis.
!      SDSALN    Real  Public   Factor for nondimensional 1-D spectrum.
!      CDSBn     Real  Public   Constants in parameterization of PHI.
!      XFH       Real  Public   Constant for turbulent length scale.
!      XFn       Real  Public   Constants in combining low and high
!                               frequency dissipation.
!                                                             (!/ST3)
!      ZZWND     Real  Public   Height at which the wind is defined
!      AALPHA    Real  Public   Minimum value of charnock parameter
!      BBETA     Real  Public   Wind-wave coupling coefficient
!      ZZALP     Real  Public   Wave age tuning coefficient in Sin
!      TTAUWSHELTER Real  Public Sheltering coefficient for short waves
!      ZZ0MAX    Real  Public   Maximum value of air-side roughness
!      ZZ0RAT    Real  Public   ratio of roughness for mean and
!                               oscillatory flows
!      SSINTHP   Real  Public   Power in cosine of wind input
!      SSWELLF   R.A.  Public   Swell damping coefficients
!      SSDSCn    Real  Public   Dissipation parameters
!      SSDSBR    Real  Public   Threshold in saturation spectrum for Sds
!      SSDSP     Real  Public   Power of B(k) in Sds
!      WWNMEANP  Real  Public   Power that defines the mean wavenumber
!                               in Sds
!      SSTXFTF, SSTXFTWN Real  Public   Tail constants
!      SSDSC4,   Real  Public   Threshold shift in saturation diss.
!      SSDSC5,   Real  Public   Wave-turbulence dissipation factor
!      SSDSC6,   Real  Public   dissipation parameter
!      DDELTA1   Real  Public   Low-frequency dissipation coefficient
!                               in WAM4
!      DDELTA2   Real  Public   High-frequency dissipation coefficient
!                               in WAM4
!      SSDSLF    Real  Public   Activator of WAM4 for non-saturated
!                               part of spectrum
!      SSDSHF    Real  Public   Activator of WAM4 for saturated part
!                               of spectrum
!      SSDSDTH   Real  Public   Maximum angular sector for saturation
!                               spectrum
!      SSDSCOS   Real  Public   Power of cosine in saturation integral
!      SSDSISO   Int.  Public   Choice of definition of the isotropic
!                               saturation
!     ----------------------------------------------------------------
!
!     The structure SNLP contains parameters for the nonl. inter.
!     source term and is aliased as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!                                                             (!/NL1)
!      SNLC1     Real  Public   Scaled proportionality constant.
!      LAM       Real  Public   Factor defining quadruplet.
!      KDCON     Real  Public   Conversion factor for relative depth.
!      KDMN      Real  Public   Minimum relative depth.
!      SNLSn     Real  Public   Constants in shallow water factor.
!                                                             (!/NL2)
!      IQTPE     Int.  Public   Type of depth treatment
!                                1 : Deep water
!                                2 : Deep water / WAM scaling
!                                3 : Finite water depth
!      NDPTHS    Int.  Public   Number of depth for which integration
!                               space needs to be computed.
!      NLTAIL    Real  Public   Tail factor for parametric tail.
!      DPTHNL    R.A.  Public   Depths corresponding to NDPTHS.
!                               *** NOTE: This array is not allocated
!                                         in the W3DIMP routine ***
!     ----------------------------------------------------------------
!
!     The structure SBTP contains parameters for the bottom friction
!     source term and is aliased as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      SBTC1     Real  Public   Proportionality constant.    (!/BT1)
!     ----------------------------------------------------------------
!
!     The structure SDBP contains parameters for the depth incduced
!     breaking source term and is aliased as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      SDBC1     Real  Public   Proportionality constant.    (!/DB1)
!      SDBC2     Real  Public   Hmax/d ratio.                (!/DB1)
!      FDONLY    Log.  Public   Flag for checking depth only (!/DB1)
!                               otherwise Miche criterion.
!     ----------------------------------------------------------------
!
!     The structure STRP contains parameters for the triad interaction
!     source term and is aliased as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!     The structure SBSP contains parameters for the bottom scattering
!     source term and is aliased as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!     The structure SXXP contains parameters for arbitrary source
!     term and is aliased as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!     ----------------------------------------------------------------
!
!  3. Subroutines and functions :
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      W3NMOD    Subr. Public   Set number of grids.
!      W3DIMX    Subr. Public   Set dimensions of spatial grid.
!      W3DIMS    Subr. Public   Set dimensions of spectral grid.
!      W3SETG    Subr. Public   Point to selected grid / model.
!     ----------------------------------------------------------------
!
!  4. Subroutines and functions used :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      STRACE    Subr. W3SERVMD Subroutine tracing.
!      EXTCDE    Subr. W3SERVMD Abort program with exit code.
!     ----------------------------------------------------------------
!
!  5. Remarks :
!
!     - In model versions before 3.06 the parameters in the grid
!       structure were stored in the module W3IOGR.
!     - No subroutine DIMP is provided, instead, arrays are set
!       one-by-one in W3IOGR.
!
!  6. Switches :
!
!     !/LLG  Grid type
!     !/XYG
!
!     !/PRn  Select propagation scheme
!
!     !/LNn  Select source terms
!     !/STn
!     !/NLn
!     !/BTn
!     !/DBn
!     !/TRn
!     !/BSn
!     !/XXn
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
!  7. Source code :
!
!/ ------------------------------------------------------------------- /
      PUBLIC
!/
!/ Conventional declarations
!/
      INTEGER                 :: NGRIDS = -1, IGRID = -1, ISGRD = -1, &
                                 IPARS = -1, NAUXGR
!
      LOGICAL, PARAMETER      :: FLAGLL = .TRUE.
!/
!/ Data structures
!/
      TYPE GRID
        INTEGER               :: NX, NY, NSEA, NSEAL, TRFLAG
        INTEGER, POINTER      :: MAPSTA(:,:), MAPST2(:,:),            &
                                 MAPFS(:,:), MAPSF(:,:)
        REAL                  :: SX, SY, X0, Y0, DTCFL, DTCFLI,       &
                                 DTMAX, DTMIN, DMIN, CTMAX,           &
                                 FICE0, FICEN, PFMOVE
        REAL, POINTER         :: ZB(:), CLAT(:), CLATI(:), CLATS(:),  &
                                 CTHG0(:), TRNX(:,:), TRNY(:,:)
        LOGICAL               :: GINIT, GLOBAL, FLDRY, FLCX, FLCY,    &
                                 FLCTH, FLCK, FLSOU
        LOGICAL, POINTER      :: FLAGST(:)
        CHARACTER(LEN=30)     :: GNAME
        CHARACTER(LEN=10)     :: FILEXT
      END TYPE GRID
!
      TYPE SGRD
        INTEGER               :: NK, NK2, NTH, NSPEC
        INTEGER, POINTER      :: MAPWN(:), MAPTH(:)
        REAL                  :: DTH, XFR, FR1, FTE, FTF, FTWN, FTTR, &
                                 FTWL, FACTI1, FACTI2, FACHFA, FACHFE
        REAL, POINTER         :: TH(:), ESIN(:), ECOS(:), ES2(:),     &
                                 ESC(:), EC2(:), SIG(:), SIG2(:),     &
                                 DSIP(:), DSII(:), DDEN(:), DDEN2(:)
        LOGICAL               :: SINIT
      END TYPE SGRD
!
      TYPE NPAR
        REAL                  :: FACP, XREL, XFLT, FXFM, FXPM,        &
                                 XFT, XFC, FACSD, FHMAX
      END TYPE NPAR
!
      TYPE PROP
        REAL                  :: WDCG, WDTH
      END TYPE PROP
!
      TYPE SFLP
        REAL                  :: DUMMY
      END TYPE SFLP
!
      TYPE SLNP
        REAL                  :: SLNC1, FSPM, FSHF
      END TYPE SLNP
!
      TYPE SRCP
        REAL                       :: WWNMEANPTAIL, SSTXFTFTAIL
        INTEGER               :: SSWELLFPAR, SSDSISO
        REAL                  :: AALPHA, BBETA, ZZ0MAX, ZZ0RAT, ZZALP,&
                                 SSINTHP, TTAUWSHELTER, SSWELLF(1:5), &
                                 SSDSC1, SSDSC2, SSDSC3, SSDSBR,      &
                                 SSDSP, WWNMEANP, SSTXFTF, SSTXFTWN,  &
                                 SSDSC4, SSDSC5, SSDSC6, DDELTA1
        REAL                  :: DDELTA2, SSDSLF, SSDSHF, ZZWND
        REAL                  :: SSDSCOS, SSDSDTH, SSDSBR2, SSDSBM(0:4)
      END TYPE SRCP
!
      TYPE SNLP
        REAL                  :: SNLC1, LAM, KDCON, KDMN,             &
                                 SNLS1, SNLS2, SNLS3
      END TYPE SNLP
!
      TYPE SBTP
        REAL                  :: SBTC1
      END TYPE SBTP
!
      TYPE SDBP
        REAL                  :: SDBC1, SDBC2
        LOGICAL               :: FDONLY
      END TYPE SDBP
!
      TYPE STRP
        REAL                  :: DUMMY
      END TYPE STRP
!
      TYPE SBSP
        REAL                  :: DUMMY
      END TYPE SBSP
!
      TYPE SXXP
        REAL                  :: DUMMY
      END TYPE SXXP
!
      TYPE MPAR
        LOGICAL               :: PINIT
        TYPE(NPAR)            :: NPARS
        TYPE(PROP)            :: PROPS
        TYPE(SFLP)            :: SFLPS
        TYPE(SLNP)            :: SLNPS
        TYPE(SRCP)            :: SRCPS
        TYPE(SNLP)            :: SNLPS
        TYPE(SBTP)            :: SBTPS
        TYPE(SDBP)            :: SDBPS
        TYPE(STRP)            :: STRPS
        TYPE(SBSP)            :: SBSPS
        TYPE(SXXP)            :: SXXPS
      END TYPE MPAR
!/
!/ Data storage
!/
      TYPE(GRID), TARGET, ALLOCATABLE :: GRIDS(:)
      TYPE(SGRD), TARGET, ALLOCATABLE :: SGRDS(:)
      TYPE(MPAR), TARGET, ALLOCATABLE :: MPARS(:)
!/
!/ Data aliasses for structure GRID(S)
!/
      INTEGER, POINTER        :: NX, NY, NSEA, NSEAL, TRFLAG
      INTEGER, POINTER        :: MAPSTA(:,:), MAPST2(:,:),            &
                                 MAPFS(:,:), MAPSF(:,:)
      REAL, POINTER           :: SX, SY, X0, Y0, DTCFL, DTCFLI,       &
                                 DTMAX, DTMIN, DMIN, CTMAX,           &
                                 FICE0, FICEN, PFMOVE
      REAL, POINTER           :: ZB(:), CLAT(:), CLATI(:), CLATS(:),  &
                                 CTHG0(:), TRNX(:,:), TRNY(:,:)
      LOGICAL, POINTER        :: GINIT, GLOBAL, FLDRY, FLCX, FLCY,    &
                                 FLCTH, FLCK, FLSOU, FLAGST(:)
      CHARACTER(LEN=30), POINTER :: GNAME
      CHARACTER(LEN=10), POINTER :: FILEXT
!/
!/ Data aliasses for structure SGRD(S)
!/
      INTEGER, POINTER        :: NK, NK2, NTH, NSPEC
      INTEGER, POINTER        :: MAPWN(:), MAPTH(:)
      REAL, POINTER           :: DTH, XFR, FR1, FTE, FTF, FTWN, FTTR, &
                                 FTWL, FACTI1, FACTI2, FACHFA, FACHFE
      REAL, POINTER           :: TH(:), ESIN(:), ECOS(:), ES2(:),     &
                                 ESC(:), EC2(:), SIG(:), SIG2(:),     &
                                 DSIP(:), DSII(:), DDEN(:), DDEN2(:)
      LOGICAL, POINTER        :: SINIT
!/
!/ Data aliasses for structure MPAR(S)
!/
      LOGICAL, POINTER        :: PINIT
!/
!/ Data aliasses for structure NPAR(S)
!/
      REAL, POINTER           :: FACP, XREL, XFLT, FXFM, FXPM,        &
                                 XFT, XFC, FACSD, FHMAX
!/
!/ Data aliasses for structure PROP(S)
!/
      REAL, POINTER           :: WDCG, WDTH
!/
!/ Data aliasses for structure SFLP(S)
!/
!/
!/ Data aliasses for structure SLNP(S)
!/
      REAL, POINTER           :: SLNC1, FSPM, FSHF
!/
!/ Data aliasses for structure SRCP(S)
!/
      INTEGER, POINTER        :: SSWELLFPAR, SSDSISO
      REAL, POINTER           :: ZZWND, AALPHA, BBETA, ZZ0MAX, ZZ0RAT, ZZALP, &
                                 SSINTHP, TTAUWSHELTER, SSWELLF(:),   &
                                 SSDSC1, SSDSC2, SSDSC3, SSDSBR,      &
                                 SSDSP, WWNMEANP, SSTXFTF, SSTXFTWN,  &
                                 SSDSC4, SSDSC5, SSDSC6, SSDSBR2,     &
                                 DDELTA1, DDELTA2, SSDSLF, SSDSHF,    &
                                 SSDSCOS, SSDSDTH, SSDSBM(:)
      REAL, POINTER           :: WWNMEANPTAIL, SSTXFTFTAIL
!/
!/ Data aliasses for structure SNLP(S)
!/
      REAL, POINTER           :: SNLC1, LAM, KDCON, KDMN,             &
                                 SNLS1, SNLS2, SNLS3
!/
!/ Data aliasses for structure SBTP(S)
!/
      REAL, POINTER           :: SBTC1
!/
!/ Data aliasses for structure SDBP(S)
!/
      REAL, POINTER           :: SDBC1, SDBC2
      LOGICAL, POINTER        :: FDONLY
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3NMOD ( NUMBER, NDSE, NDST, NAUX )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         18-Jul-2006 !
!/                  +-----------------------------------+
!/
!/    24-Feb-2004 : Origination.                        ( version 3.06 )
!/    18-Jul-2006 : Add input grids.                    ( version 3.10 )
!/
!  1. Purpose :
!
!     Set up the number of grids to be used.
!
!  2. Method :
!
!     Store in NGRIDS and allocate GRIDS.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       NUMBER  Int.   I   Number of grids to be used.
!       NDSE    Int.   I   Error output unit number.
!       NDST    Int.   I   Test output unit number.
!       NAUX    Int.   I   Number of auxiliary grids to be used.
!                          Grids -NAUX:NUBMER are defined, optional
!                          parameters.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation.
!
!  5. Called by :
!
!     Any program that uses this grid structure.
!
!  6. Error messages :
!
!     - Error checks on previous setting of variable.
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)           :: NUMBER, NDSE, NDST
      INTEGER, INTENT(IN), OPTIONAL :: NAUX
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER                 :: I, NLOW
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( NGRIDS .NE. -1 ) THEN
          WRITE (NDSE,1001) NGRIDS
          CALL EXTCDE (1)
        END IF
!
      IF ( NUMBER .LT. 1 ) THEN
          WRITE (NDSE,1002) NUMBER
          CALL EXTCDE (2)
        END IF
!
      IF ( PRESENT(NAUX) ) THEN
          NLOW   = -NAUX
        ELSE
          NLOW   = 1
        END IF
!
      IF ( NLOW .GT. 1 ) THEN
          WRITE (NDSE,1003) -NLOW
          CALL EXTCDE (3)
        END IF
!
! -------------------------------------------------------------------- /
! 1.  Set variable and allocate arrays
!
      NGRIDS = NUMBER
      NAUXGR = - NLOW
      ALLOCATE ( GRIDS(NLOW:NUMBER) )
      ALLOCATE ( SGRDS(NLOW:NUMBER) )
      ALLOCATE ( MPARS(NLOW:NUMBER) )
!
! -------------------------------------------------------------------- /
! 2.  Initialize GINIT adn SINIT
!
      DO I=NLOW, NUMBER
        GRIDS(I)%GINIT  = .FALSE.
        SGRDS(I)%SINIT  = .FALSE.
        MPARS(I)%PINIT  = .FALSE.
        END DO
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3NMOD : GRIDS ALREADY INITIALIZED *** '/  &
               '                    NGRIDS = ',I10/)
 1002 FORMAT (/' *** ERROR W3NMOD : ILLEGAL NUMBER OF GRIDS *** '/    &
               '                    NUMBER = ',I10/)
 1003 FORMAT (/' *** ERROR W3NMOD : ILLEGAL NUMBER OF AUX GRIDS *** '/&
               '                    NUMBER = ',I10/)
!
!/
!/ End of W3NMOD ----------------------------------------------------- /
!/
      END SUBROUTINE W3NMOD
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DIMX  ( IMOD, MX, MY, MSEA, NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         02-Feb-2007 !
!/                  +-----------------------------------+
!/
!/    24-Jun-2005 : Origination.                        ( version 3.07 )
!/    18-Jul-2006 : Add input grids.                    ( version 3.10 )
!/    05-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    02-Feb-2007 : Add FLAGST.                         ( version 3.10 )
!/
!  1. Purpose :
!
!     Initialize an individual spatial grid at the proper dimensions.
!
!  2. Method :
!
!     Allocate directly into the structure array GRIDS. Note that
!     this cannot be done through the pointer alias!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number to point to.
!       NDSE    Int.   I   Error output unit number.
!       NDST    Int.   I   Test output unit number.
!       MX, MY, MSEA       Like NX, NY, NSEA in data structure.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!       See module documentation.
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3IOGR    Subr. W3IOGRMD Model definition file IO program.
!      WW3_GRID  Prog.   N/A    Model set up program.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - Check on input parameters.
!     - Check on previous allocation.
!
!  7. Remarks :
!
!     - Grid dimensions apre passed through parameter list and then
!       locally stored to assure consistency between allocation and
!       data in structure.
!     - W3SETG needs to be called after allocation to point to
!       proper allocated arrays.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMOD, MX, MY, MSEA, NDSE, NDST
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( NGRIDS .EQ. -1 ) THEN
          WRITE (NDSE,1001)
          CALL EXTCDE (1)
        END IF
!
      IF ( IMOD.LT.-NAUXGR .OR. IMOD.GT.NGRIDS ) THEN
          WRITE (NDSE,1002) IMOD, -NAUXGR, NGRIDS
          CALL EXTCDE (2)
        END IF
!
      IF ( MX.LT.3 .OR. MY.LT.3 .OR. MSEA.LT.1 ) THEN
          WRITE (NDSE,1003) MX, MY, NSEA, NSEAL
          CALL EXTCDE (3)
        END IF
!
      IF ( GRIDS(IMOD)%GINIT ) THEN
          WRITE (NDSE,1004)
          CALL EXTCDE (4)
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays
!
      ALLOCATE ( GRIDS(IMOD)%MAPSTA(MY,MX),                           &
                 GRIDS(IMOD)%MAPST2(MY,MX),                           &
                 GRIDS(IMOD)%MAPFS(MY,MX),                            &
                 GRIDS(IMOD)%MAPSF(MSEA,3),                           &
                 GRIDS(IMOD)%FLAGST(MSEA),                            &
                 GRIDS(IMOD)%ZB(MSEA),                                &
                 GRIDS(IMOD)%CLAT(MY),                                &
                 GRIDS(IMOD)%CLATI(MY),                               &
                 GRIDS(IMOD)%CLATS(MSEA),                             &
                 GRIDS(IMOD)%CTHG0(MY),                               &
                 GRIDS(IMOD)%TRNX(MY,MX),                             &
                 GRIDS(IMOD)%TRNY(MY,MX) )
!
      GRIDS(IMOD)%FLAGST = .TRUE.
      GRIDS(IMOD)%GINIT  = .TRUE.
!
! -------------------------------------------------------------------- /
! 3.  Point to allocated arrays
!
      CALL W3SETG ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 4.  Update counters in grid
!
      NX     = MX
      NY     = MY
      NSEA   = MSEA
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3DIMX : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3DIMX : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NAUXGR = ',I10/                   &
               '                    NGRIDS = ',I10/)
 1003 FORMAT (/' *** ERROR W3DIMX : ILLEGAL GRID DIMENSION(S) *** '/  &
               '                    INPUT = ',I10/)
 1004 FORMAT (/' *** ERROR W3DIMX : ARRAY(S) ALREADY ALLOCATED *** ')
!
!/
!/ End of W3DIMX ----------------------------------------------------- /
!/
      END SUBROUTINE W3DIMX
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DIMS  ( IMOD, MK, MTH, NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         05-Oct-2006 !
!/                  +-----------------------------------+
!/
!/    19-Feb-2004 : Origination.                        ( version 3.06 )
!/    18-Jul-2006 : Add input grids.                    ( version 3.10 )
!/    05-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/
!  1. Purpose :
!
!     Initialize an individual spatial grid at the proper dimensions.
!
!  2. Method :
!
!     Allocate directly into the structure array GRIDS. Note that
!     this cannot be done through the pointer alias!
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number to point to.
!       NDSE    Int.   I   Error output unit number.
!       MK,MTH  Int.   I   Spectral dimensions.
!       NDST    Int.   I   Test output unit number.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation.
!
!  5. Called by :
!
!      Name      Type  Module   Description
!     ----------------------------------------------------------------
!      W3IOGR    Subr. W3IOGRMD Model definition file IO program.
!      WW3_GRID  Prog.   N/A    Model set up program.
!     ----------------------------------------------------------------
!
!  6. Error messages :
!
!     - Check on input parameters.
!     - Check on previous allocation.
!
!  7. Remarks :
!
!     - Grid dimensions apre passed through parameter list and then
!       locally stored to assure consistency between allocation and
!       data in structure.
!     - W3SETG needs to be called after allocation to point to
!       proper allocated arrays.
!
!  8. Structure :
!
!     See source code.
!
!  9. Switches :
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMOD, MK, MTH, NDSE, NDST
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      INTEGER, SAVE           :: MK2, MSPEC
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( NGRIDS .EQ. -1 ) THEN
          WRITE (NDSE,1001)
          CALL EXTCDE (1)
        END IF
!
      IF ( IMOD.LT.-NAUXGR .OR. IMOD.GT.NGRIDS ) THEN
          WRITE (NDSE,1002) IMOD, -NAUXGR, NGRIDS
          CALL EXTCDE (2)
        END IF
!
      IF ( MK.LT.3 .OR. MTH.LT.4 ) THEN
          WRITE (NDSE,1003) MK, MTH
          CALL EXTCDE (3)
        END IF
!
      IF ( SGRDS(IMOD)%SINIT ) THEN
          WRITE (NDSE,1004)
          CALL EXTCDE (4)
        END IF
!
      MK2    = MK + 2
      MSPEC  = MK * MTH
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays
!
      ALLOCATE ( SGRDS(IMOD)%MAPWN(MSPEC+MTH),                        &
                 SGRDS(IMOD)%MAPTH(MSPEC+MTH),                        &
                 SGRDS(IMOD)%TH(MTH),                                 &
                 SGRDS(IMOD)%ESIN(MSPEC+MTH),                         &
                 SGRDS(IMOD)%ECOS(MSPEC+MTH),                         &
                 SGRDS(IMOD)%ES2(MSPEC+MTH),                          &
                 SGRDS(IMOD)%ESC(MSPEC+MTH),                          &
                 SGRDS(IMOD)%EC2(MSPEC+MTH),                          &
                 SGRDS(IMOD)%SIG(0:MK+1),                             &
                 SGRDS(IMOD)%SIG2(MSPEC),                             &
                 SGRDS(IMOD)%DSIP(0:MK+1),                            &
                 SGRDS(IMOD)%DSII(MK),                                &
                 SGRDS(IMOD)%DDEN(MK),                                &
                 SGRDS(IMOD)%DDEN2(MSPEC)  )
!
      SGRDS(IMOD)%SINIT  = .TRUE.
!
! -------------------------------------------------------------------- /
! 3.  Point to allocated arrays
!
      CALL W3SETG ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 4.  Update counters in grid
!
      NK     = MK
      NK2    = MK + 2
      NTH    = MTH
      NSPEC  = MK * MTH
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3DIMS : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3DIMS : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NAUXGR = ',I10/                   &
               '                    NGRIDS = ',I10/)
 1003 FORMAT (/' *** ERROR W3DIMS : ILLEGAL GRID DIMENSION(S) *** '/  &
               '                    INPUT = ',4I10/)
 1004 FORMAT (/' *** ERROR W3DIMS : ARRAY(S) ALREADY ALLOCATED *** ')
!
!/
!/ End of W3DIMS ----------------------------------------------------- /
!/
      END SUBROUTINE W3DIMS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SETG ( IMOD, NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  !           J. H. Alves             !
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         27-Jun-2008 |
!/                  +-----------------------------------+
!/
!/    24-Jun-2005 : Origination.                        ( version 3.07 )
!/    09-Nov-2005 : Remove soft boundary options.       ( version 3.08 )
!/    23-Jun-2006 : Add data for W3SLN1.                ( version 3.09 )
!/    18-Jul-2006 : Add input grids.                    ( version 3.10 )
!/    05-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    02-Feb-2007 : Add FLAGST.                         ( version 3.10 )
!/    14-Apr-2007 : Add Miche style limiter.            ( version 3.11 )
!/                  ( J. H. Alves )
!/    25-Apr-2007 : Adding Battjes-Janssen Sdb.         ( version 3.11 )
!/                  ( J. H. Alves )
!/    18-Sep-2007 : Adding WAM4 source terms.           ( version 3.13 )
!/                  ( F. Ardhuin )
!/    27-Jun-2008 : Expand WAM4 variants namelist       ( version 3.14 )
!/                  ( F. Ardhuin )
!/
!  1. Purpose :
!
!     Select one of the WAVEWATCH III grids / models.
!
!  2. Method :
!
!     Point pointers to the proper variables in the proper element of
!     the GRIDS array.
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number to point to.
!       NDSE    Int.   I   Error output unit number.
!       NDST    Int.   I   Test output unit number.
!     ----------------------------------------------------------------
!
!  4. Subroutines used :
!
!     See module documentation.
!
!  5. Called by :
!
!     Many subroutines in eth WAVEWATCH system.
!
!  6. Error messages :
!
!     Checks on parameter list IMOD.
!
!  7. Remarks :
!
!  8. Structure :
!
!  9. Switches :
!
!     !/PRn  Select propagation scheme
!
!     !/STn  Select source terms
!     !/NLn
!     !/BTn
!
!     !/S    Enable subroutine tracing.
!     !/T    Enable test output
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE W3SERVMD, ONLY: EXTCDE
!
      IMPLICIT NONE
!
!/
!/ ------------------------------------------------------------------- /
!/ Parameter list
!/
      INTEGER, INTENT(IN)     :: IMOD, NDSE, NDST
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
!/
!
! -------------------------------------------------------------------- /
! 1.  Test input and module status
!
      IF ( NGRIDS .EQ. -1 ) THEN
          WRITE (NDSE,1001)
          CALL EXTCDE (1)
        END IF
!
      IF ( IMOD.LT.-NAUXGR .OR. IMOD.GT.NGRIDS ) THEN
          WRITE (NDSE,1002) IMOD, -NAUXGR, NGRIDS
          CALL EXTCDE (2)
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Set model numbers
!
      IGRID  = IMOD
      ISGRD  = IMOD
      IPARS  = IMOD
!
! -------------------------------------------------------------------- /
! 3.  Set pointers in structure GRID
!
      NX     => GRIDS(IMOD)%NX
      NY     => GRIDS(IMOD)%NY
      NSEA   => GRIDS(IMOD)%NSEA
      NSEAL  => GRIDS(IMOD)%NSEAL
      TRFLAG => GRIDS(IMOD)%TRFLAG
!
      SX     => GRIDS(IMOD)%SX
      SY     => GRIDS(IMOD)%SY
      X0     => GRIDS(IMOD)%X0
      Y0     => GRIDS(IMOD)%Y0
      DTCFL  => GRIDS(IMOD)%DTCFL
      DTCFLI => GRIDS(IMOD)%DTCFLI
      DTMAX  => GRIDS(IMOD)%DTMAX
      DTMIN  => GRIDS(IMOD)%DTMIN
      DMIN   => GRIDS(IMOD)%DMIN
      CTMAX  => GRIDS(IMOD)%CTMAX
      FICE0  => GRIDS(IMOD)%FICE0
      FICEN  => GRIDS(IMOD)%FICEN
      PFMOVE => GRIDS(IMOD)%PFMOVE
!
      GINIT  => GRIDS(IMOD)%GINIT
      GLOBAL => GRIDS(IMOD)%GLOBAL
      FLDRY  => GRIDS(IMOD)%FLDRY
      FLCX   => GRIDS(IMOD)%FLCX
      FLCY   => GRIDS(IMOD)%FLCY
      FLCTH  => GRIDS(IMOD)%FLCTH
      FLCK   => GRIDS(IMOD)%FLCK
      FLSOU  => GRIDS(IMOD)%FLSOU
!
      GNAME  => GRIDS(IMOD)%GNAME
      FILEXT => GRIDS(IMOD)%FILEXT
!
      IF ( GINIT ) THEN
!
          MAPSTA => GRIDS(IMOD)%MAPSTA
          MAPST2 => GRIDS(IMOD)%MAPST2
          MAPFS  => GRIDS(IMOD)%MAPFS
          MAPSF  => GRIDS(IMOD)%MAPSF
          FLAGST => GRIDS(IMOD)%FLAGST
!
          ZB     => GRIDS(IMOD)%ZB
          CLAT   => GRIDS(IMOD)%CLAT
          CLATI  => GRIDS(IMOD)%CLATI
          CLATS  => GRIDS(IMOD)%CLATS
          CTHG0  => GRIDS(IMOD)%CTHG0
          TRNX   => GRIDS(IMOD)%TRNX
          TRNY   => GRIDS(IMOD)%TRNY
!
        END IF
!
! -------------------------------------------------------------------- /
! 4.  Set pointers in structure SGRD
!
      NK     => SGRDS(IMOD)%NK
      NK2    => SGRDS(IMOD)%NK2
      NTH    => SGRDS(IMOD)%NTH
      NSPEC  => SGRDS(IMOD)%NSPEC
!
      DTH    => SGRDS(IMOD)%DTH
      XFR    => SGRDS(IMOD)%XFR
      FR1    => SGRDS(IMOD)%FR1
      FTE    => SGRDS(IMOD)%FTE
      FTF    => SGRDS(IMOD)%FTF
      FTWN   => SGRDS(IMOD)%FTWN
      FTTR   => SGRDS(IMOD)%FTTR
      FTWL   => SGRDS(IMOD)%FTWL
      FACTI1 => SGRDS(IMOD)%FACTI1
      FACTI2 => SGRDS(IMOD)%FACTI2
      FACHFA => SGRDS(IMOD)%FACHFA
      FACHFE => SGRDS(IMOD)%FACHFE
!
      SINIT  => SGRDS(IMOD)%SINIT
!
      IF ( SINIT ) THEN
!
          MAPWN  => SGRDS(IMOD)%MAPWN
          MAPTH  => SGRDS(IMOD)%MAPTH
!
          TH     => SGRDS(IMOD)%TH
          ESIN   => SGRDS(IMOD)%ESIN
          ECOS   => SGRDS(IMOD)%ECOS
          ES2    => SGRDS(IMOD)%ES2
          ESC    => SGRDS(IMOD)%ESC
          EC2    => SGRDS(IMOD)%EC2
          SIG    => SGRDS(IMOD)%SIG
          SIG2   => SGRDS(IMOD)%SIG2
          DSIP   => SGRDS(IMOD)%DSIP
          DSII   => SGRDS(IMOD)%DSII
          DDEN   => SGRDS(IMOD)%DDEN
          DDEN2  => SGRDS(IMOD)%DDEN2
!
        END IF
!
! -------------------------------------------------------------------- /
! 5.  Set pointers in structure MPAR
!
      PINIT  => MPARS(IMOD)%PINIT
!
!     Structure NPARS
!
      FACP   => MPARS(IMOD)%NPARS%FACP
      XREL   => MPARS(IMOD)%NPARS%XREL
      XFLT   => MPARS(IMOD)%NPARS%XFLT
      FXFM   => MPARS(IMOD)%NPARS%FXFM
      FXPM   => MPARS(IMOD)%NPARS%FXPM
      XFT    => MPARS(IMOD)%NPARS%XFT
      XFC    => MPARS(IMOD)%NPARS%XFC
      FACSD  => MPARS(IMOD)%NPARS%FACSD
      FHMAX  => MPARS(IMOD)%NPARS%FHMAX
!
!     Structure PROPS
!
      WDCG   => MPARS(IMOD)%PROPS%WDCG
      WDTH   => MPARS(IMOD)%PROPS%WDTH
!
!     Structure SFLPS
!
!     Structure SLNPS
!
      SLNC1  => MPARS(IMOD)%SLNPS%SLNC1
      FSPM   => MPARS(IMOD)%SLNPS%FSPM
      FSHF   => MPARS(IMOD)%SLNPS%FSHF
!
!     Structure SRCPS
!
      WWNMEANPTAIL=> MPARS(IMOD)%SRCPS%WWNMEANPTAIL
      SSTXFTFTAIL => MPARS(IMOD)%SRCPS%SSTXFTFTAIL
      ZZWND  => MPARS(IMOD)%SRCPS%ZZWND
      AALPHA => MPARS(IMOD)%SRCPS%AALPHA
      BBETA  => MPARS(IMOD)%SRCPS%BBETA
      SSINTHP  => MPARS(IMOD)%SRCPS%SSINTHP
      ZZ0MAX  => MPARS(IMOD)%SRCPS%ZZ0MAX
      ZZ0RAT  => MPARS(IMOD)%SRCPS%ZZ0RAT
      ZZALP  => MPARS(IMOD)%SRCPS%ZZALP
      TTAUWSHELTER  => MPARS(IMOD)%SRCPS%TTAUWSHELTER
      SSWELLFPAR  => MPARS(IMOD)%SRCPS%SSWELLFPAR
      SSWELLF  => MPARS(IMOD)%SRCPS%SSWELLF
      SSDSC1 => MPARS(IMOD)%SRCPS%SSDSC1
      SSDSC2 => MPARS(IMOD)%SRCPS%SSDSC2
      SSDSC3 => MPARS(IMOD)%SRCPS%SSDSC3
      SSDSC4 => MPARS(IMOD)%SRCPS%SSDSC4
      SSDSC5 => MPARS(IMOD)%SRCPS%SSDSC5
      SSDSC6 => MPARS(IMOD)%SRCPS%SSDSC6
      SSDSBR => MPARS(IMOD)%SRCPS%SSDSBR
      SSDSBR2 => MPARS(IMOD)%SRCPS%SSDSBR2
      SSDSBM => MPARS(IMOD)%SRCPS%SSDSBM
      SSDSP  => MPARS(IMOD)%SRCPS%SSDSP
      WWNMEANP => MPARS(IMOD)%SRCPS%WWNMEANP
      DDELTA1 => MPARS(IMOD)%SRCPS%DDELTA1
      DDELTA2 => MPARS(IMOD)%SRCPS%DDELTA2
      SSDSLF => MPARS(IMOD)%SRCPS%SSDSLF
      SSDSHF => MPARS(IMOD)%SRCPS%SSDSHF
      SSDSDTH => MPARS(IMOD)%SRCPS%SSDSDTH
      SSTXFTF => MPARS(IMOD)%SRCPS%SSTXFTF
      SSTXFTWN => MPARS(IMOD)%SRCPS%SSTXFTWN
      SSDSCOS => MPARS(IMOD)%SRCPS%SSDSCOS
      SSDSISO => MPARS(IMOD)%SRCPS%SSDSISO
!
!     Structure SRNLS
!
      SNLC1  => MPARS(IMOD)%SNLPS%SNLC1
      LAM    => MPARS(IMOD)%SNLPS%LAM
      KDCON  => MPARS(IMOD)%SNLPS%KDCON
      KDMN   => MPARS(IMOD)%SNLPS%KDMN
      SNLS1  => MPARS(IMOD)%SNLPS%SNLS1
      SNLS2  => MPARS(IMOD)%SNLPS%SNLS2
      SNLS3  => MPARS(IMOD)%SNLPS%SNLS3
!
!     Structure SBTPS
!
      SBTC1  => MPARS(IMOD)%SBTPS%SBTC1
!
!     Structure SDBPS
!
      SDBC1  => MPARS(IMOD)%SDBPS%SDBC1
      SDBC2  => MPARS(IMOD)%SDBPS%SDBC2
      FDONLY => MPARS(IMOD)%SDBPS%FDONLY
!
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3SETG : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3SETG : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NAUXGR = ',I10/                   &
               '                    NGRIDS = ',I10/)
!
!/
!/ End of W3SETG ----------------------------------------------------- /
!/
      END SUBROUTINE W3SETG
!/
!/ End of module W3GDATMD -------------------------------------------- /
!/
      END MODULE W3GDATMD
