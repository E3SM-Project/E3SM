#include "w3macros.h"
!/ ------------------------------------------------------------------- /
      MODULE W3GDATMD
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  !           J. H. Alves             !
!/                  |            F. Ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         04-May-2016 |
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
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    30-Oct-2009 : Implement curvilinear grid type.    ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    29-Oct-2010 : Implement unstructured grids        ( version 3.14.1 )
!/                  (A. Roland and F. Ardhuin)
!/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
!/                  specify index closure for a grid.   ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    23-Dec-2010 : Fix HPFAC and HQFAC by including the COS(YGRD)
!/                  factor with DXDP and DXDQ terms.    ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    05-Apr-2011 : Implement interations for DTMAX < 1s( version 3.14.1 )
!/                  (F. Ardhuin)
!/    01-Jul-2011 : Movable bed bottom friction BT4     ( version 4.01 )
!/    03-Nov-2011 : Bug fix: GUGINIT initialization     ( version 4.04 )
!/    29-Nov-2011 : Adding ST6 source term option.      ( version 4.04 )
!/                  (S. Zieger)
!/    14-Mar-2012 : Add PSIC for BT4                    ( version 4.04 )
!/    12-Jun-2012 : Add /RTD option or rotated grid variables.
!/                  (Jian-Guo Li)                       ( version 4.06 )
!/    13-Jul-2012 : Move data structures GMD (SNL3) and nonlinear
!/                  filter (SNLS) from 3.15 (HLT).      ( version 4.08 )
!/    03-Sep-2012 : Clean up of UG grids                ( version 4.08 )
!/    12-Dec-2012 : Adding SMC grid.  JG_Li             ( version 4.09 )
!/    16-Sep-2013 : Add Arctic part SMC grid.           ( version 4.11 )
!/    11-Nov-2013 : SMC and rotated grid incorporated in the main
!/                  trunk                               ( version 4.13 )
!/    16-Nov-2013 : Allows reflection on curvi grids    ( version 4.14 )
!/    26-Jul-2013 : Adding IG waves                     ( version 4.16 )
!/    18-Dec-2013 : Moving FLAGLL into GRID TYPE        ( version 4.16 )
!/    11-Jun-2014 : Changed reflection for subgrid      ( version 5.01 )
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
!/    21-Aug-2015 : Add SMC FUNO3, FVERG options. JGLi  ( version 5.09 )
!/    04-May-2016 : Add IICEDISP                  GB&FA ( version 5.10 )
!/
!/
!/    Copyright 2009-2013 National Weather Service (NWS),
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
!      RLGTYPE   I.P.  Public   Named constant for rectilinear grid type
!      CLGTYPE   I.P.  Public   Named constant for curvilinear grid type
!      FLAGLL    Log.  Public   Flag to indicate coordinate system for all grids
!                               .TRUE.: Spherical (lon/lat in degrees)
!                               .FALSE.: Cartesian (meters)
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
!      GTYPE     Int.  Public   Flag for type of grid
!                               RLGTYPE: Rectilinear grid
!                               CLGTYPE: Curvilinear grid
!                               UNGTYPE: Unstructured triangular grid
!      ICLOSE    Int.  Public   Parameter indicating type of index closure of grid.
!                               ICLOSE_NONE: No grid closure
!                               ICLOSE_SMPL: Simple grid closure
!                                 Grid is periodic in the i-index and wraps at
!                                 I=NX+1. In other words, (NX+1,J) => (1,J).
!                               ICLOSE_TRPL: Tripole grid closure
!                                 Grid is periodic in the i-index and and wraps at
!                                 I=NX+1 and has closure at J=NY+1. In other words,
!                                 (NX+1,J<=NY) => (1,J) and
!                                 (I,NY+1) => (MOD(NX-I+1,NX)+1,NY). The tripole
!                                 closure requires that NX be even.
!      NX, NY    Int.  Public   Discrete dimensions of spatial grid.
!      NSEA(L)   Int.  Public   Number of sea points (local for MPP).
!      NU/VFc    Int.  Public   Number of U/V faces for SMC grid.
!      NRLv      Int.  Public   Number of refined levels for SMC grid.
!      NGLO      Int.  Public   Number of cells in global part for SMC grid.
!      NARC      Int.  Public   Number of cells in Arctic part for SMC grid.
!      NBAC      Int.  Public   Number of boundary cells in Arctic part.
!      NBGL      Int.  Public   Number of boundary cells in global part.
!      TRFLAG    Int.  Public   Flag for use of transparencies
!                                0: No sub-grid obstacles.
!                                1: Obstructions at cell boundaries.
!                                2: Obstructions at cell centers.
!                                3: Like 1 with continuous ice.
!                                4: Like 2 with continuous ice.
!      MAPSTA    I.A.  Public   Grid status map.
!      MAPST2    I.A.  Public   Second grid status map.
!      MAPxx     I.A.  Public   Storage grid maps.
!      IJKCel    I.A.  Public   Cell info array for SMC grid.
!      IJKU/VFc  I.A.  Public   U/V-Face arrays for SMC grid.
!      NLv*      I.A.  Public   Cell, U/V-Face numbers of refine levels.
!      ICLBAC    I.A.  Public   Mapping index for Arctic boundary cells.
!      SX,SY     Real  Public   Spatial (rectilinear) grid increments.
!      X0,Y0     Real  Public   Lower left corner of spatial (rectilinear) grid.
!      DTCFL     Real  Public   Maximum CFL time step X-Y propagation.
!      DTCFLI    Real  Public   Id. intra-spectral.
!      DTMAX     Real  Public   Maximum overall time step.
!      DTMIN     Real  Public   Minimum dynamic time step for source
!      NITERSEC1 Real  Public   Number of interations when DTMAX < 1s
!      DMIN      Real  Public   Minimum water depth.
!      CTMAX     Real  Public   Maximum CFL number for depth refr.
!      FICE0/N   Real  Public   Cut-off ice conc. for ice coverage.
!      FICEL     Real  Public   Length scale for sea ice damping
!      IICEHMIN  Real  Public   Minimum thickness of sea ice
!      IICEHINIT Real  Public   Initial value of ice thickness
!      ICE100WIND Real  Publ.   Reduction of wind input if ice fraction=1.0
!                               Default is 1.0, meaning that 100% ice
!                               concentration result in zero wind input.
!                               Sin_in_ice=Sin_in_open_water * (1-ICE*ICEWIND)
!                               If set to 0.0, then ice has no direct impact on
!                               wind input source function.
!      IC3PARS   R.A.  Public   various parameters for use in IC4, handled as
!                               an array for simplicity
!      IC4_KI    R.A.  Public   KI (dissipation rate) values for use in IC4
!      IC4_FC    R.A.  Public   FC (frequency bin separators) for use in IC4
!      PFMOVE    Real  Public   Tunable parameter in GSE correction
!                               for moving grids.
!      GRIDSHIFT Real  Public   Grid offset for multi-grid w/SCRIP
!      PoLat/Lon R.A.  Public   Rotated N-Pole standard latitude/longitude.
!      AnglD     R.A.  Public   Rotation angle in degree to turn rotated grid
!                               back to standard grid.  JGLi12Jun2012
!      STEXU     Real  Public   Length-scale (X) for space-time extreme averaging
!      STEYU     Real  Public   Length-scale (Y) for space-time extreme averaging
!      STEDU     Real  Public   Time-scale for space-time extreme averaging
!      ZB        R.A.  Public   Bottom levels on storage grid.
!      CLATS(I)  R.A.  Public   (Inverse) cosine of latitude at sea points.
!      CTHG0S    R.A.  Public   Constant in great-circle refr. term at sea points.
!      TRNX/Y    R.A.  Public   Transparencies in X/Y for sub-grid
!      CTRNX/Y   R.A.  Public   Sub-grid transparencies for SMC grid.
!      ANGARC    R.A.  Public   Rotation angle in degree for Arctic cells.
!      SPCBAC    R.A.  Public   Full 2-D spectra for Arctic boundary cells.
!      X/YGRD    R.A.  Public   Spatial grid coordinate arrays.
!      SX/SYGRD  R.A.  Public   Spatial grid increment arrays.
!      GINIT     Log.  Public   Flag identifying grid initialization.
!      FLDRY     Log.  Public   Flag for 'dry' run (IO and data
!                               processing only).
!      FLCx      Log.  Public   Flags for prop. is different spaces.
!      FLSOU     Log.  Public   Flag for source term calcualtion.
!      FUNO3     Log.  Public   Flag for 3rd order UNO3 scheme on SMC grid.
!      FVERG     Log.  Public   Flag for 1-2-1 averaging smoothing on SMC grid.
!      FLAGST    L.A.  Public   Flag for source term computations
!                               for individual grid points.
!      IICEDISP   Log.  Public   Flag for use of the ice covered dispertion relation.
!
!      GNAME     C*30  Public   Grid name.
!      FILEXT    C*13  Public   Extension of WAVEWATCH III file names
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
!      FHMAX     Real  Public   Hs/depth ratio in limiter     (!/MLIM)
!      RWINDC    Real  Public   Coefficient for current in relative
!                               wind                          (!/RWND)
!      WWCOR     R.A.  Public   Wind correction factors       (!/WCOR)
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
!                                                            (!/FLX4)
!      FLX4A0    Real  Public   Scaling value in parametric description
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
!                                                             (!/NL3)
!      NFR       Int.  Public   Number of frequencies or wavenumbers
!                               in discrete spectral space (NFR=>NK).
!      NFRMIN    Int.  Public   Minimum discrete frequency in the
!                               expanded frequency space.
!      NFRMAX    Int.  Public   Idem maximum for first part.
!      NFRCUT    Int.  Public   Idem maximum for second part.
!      NTHMAX    Int.  Public   Extension of directional space.
!      NTHEXP    Int   Public   Number of bins in extended dir. space.
!      NSPMIN, NSPMAX, NSPMX2
!                Int.  Public   1D spectral space range.
!      FRQ       R.A.  Public   Expanded frequency range (Hz).
!      XSI       R.A.  Public   Expanded frequency range (rad/s).
!      NQA       Int.  Public   Number of actual quadruplets.
!      QST1      I.A.  Public   Spectral offsets for compuation of
!                               quadruplet spectral desnities.
!      QST2      R.A.  Public   Idem weights.
!      QST3      R.A.  Public   Proportionality constants and k factors
!                               in diagonal strength.
!      QST4      I.A.  Public   Spectral offsets for combining of
!                               interactions and diagonal.
!      QST5      R.A.  Public   Idem weights for interactions.
!      QST6      R.A.  Public   Idem weights for diagonal.
!      SNLNQ     Int.  Public   Number of quadruplet definitions.
!      SNLMSC    Real  Public   Tuning power 'deep' scaling.
!      SNLNSC    Real  Public   Tuning power 'shallow' scaling.
!      SNLSFD    Real  Public   'Deep' nondimensional filer freq.
!      SNLSFS    Real  Public   'Shallow' nondimensional filer freq.
!      SNLL      R.A.  Public   Array with lambda for quadruplet.
!      SNLM      R.A.  Public   Array with mu for quadruplet.
!      SNLT      R.A.  Public   Array with Dtheta for quadruplet.
!      SNLCD     R.A.  Public   Array with Cd for quadruplet.
!      SNLCS     R.A.  Public   Array with Cs for quadruplet.
!                                                             (!/NL4)
!      ITSA      Int.  Public   Integer indicating TSA (1) or FBI (0)
!      IALT      Int.  Public   Integer determining alternating looping
!                                                             (!/NLS)
!      NTHX      Int.  Public   Expanded discrete direction range.
!      NFRX      Int.  Public   Expanded discrete frequency range.
!      NSPL-H    Int.  Public   Range of 1D spectrum.
!      SNSST     R.A.  Public   Array with interpolation weights.
!      CNLSA     Real  Public   a34 in quadruplet definition.
!      CNLSC     Real  Public   C in Snl definition.
!      CNLSFM    Real  Public   Maximum relative spectral change.
!      CNLSC1/3  Real  Public   Constant in frequency filter.
!     ----------------------------------------------------------------
!
!     The structure SBTP contains parameters for the bottom friction
!     source term and is aliased as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!      SBTC1     Real  Public   Proportionality constant.    (!/BT1)
!      SBTCX     R.A.  Public   Parameters for bottom fric.  (!/BT4)
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
!     The structure SICP contains parameters for arbitrary source
!     term and is aliased as above:
!
!      Name      Type  Scope    Description
!     ----------------------------------------------------------------
!     IS1C1      Real  Public   Scale factor for icecon.     (!/ISx)
!     IS1C2      Real  Public   Offset for ice concentration (!/ISx)
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
!      W3GNTX    Subr. Public   Construct grid arrays
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
!     !/PRn  Select propagation scheme
!     !/SMC  UNO2 propagation on SMC grid.
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
!/
!/ Required modules
!/
      USE W3GSRUMD
!/
!/ Specify default accessibility
!/
      PUBLIC
!/
!/ Module private variable for checking error returns
!/
      INTEGER, PRIVATE        :: ISTAT
!/
!/ Conventional declarations
!/
      INTEGER                 :: NGRIDS = -1, IGRID = -1, ISGRD = -1, &
                                 IPARS = -1, NAUXGR
!
      INTEGER, PARAMETER      :: RLGTYPE = 1
      INTEGER, PARAMETER      :: CLGTYPE = 2
      INTEGER, PARAMETER      :: UNGTYPE = 3
 
      INTEGER, PARAMETER      :: ICLOSE_NONE = ICLO_NONE
      INTEGER, PARAMETER      :: ICLOSE_SMPL = ICLO_SMPL
      INTEGER, PARAMETER      :: ICLOSE_TRPL = ICLO_TRPL
!
! Dimensions of tables for pre-computing of dissipation
!
!/
!/ Data structures
!/
!/ Grid type
      TYPE GRID          ! this is the geographical grid with all associated parameters
        INTEGER          :: GTYPE
        INTEGER          :: ICLOSE
        INTEGER          :: NX, NY, NSEA, NSEAL, TRFLAG
        INTEGER, POINTER :: MAPSTA(:,:), MAPST2(:,:),            &
                            MAPFS(:,:), MAPSF(:,:)
!
        REAL             :: SX, SY, X0, Y0, DTCFL, DTCFLI, DTMAX,      &
                            DTMIN, DMIN, CTMAX, FICE0, FICEN, FICEL,   &
                            PFMOVE, STEXU, STEYU, STEDU, IICEHMIN,     &
                            IICEHINIT, ICE100WIND
        REAL(8)          :: GRIDSHIFT ! see notes in WMGHGH
 
 
        REAL   , POINTER :: ZB(:)     ! BOTTOM GRID, DEFINED ON ISEA
        REAL   , POINTER :: CLATS(:)  ! COS(LAT), DEFINED ON SEA POINTS
        REAL   , POINTER :: CLATIS(:) ! INVERSE OF COS(LAT) DEFINED ON ISEA
        REAL   , POINTER :: CTHG0S(:) ! TAN(Y)/R, DEFINED ON ISEA
 
        REAL   , POINTER :: TRNX(:,:), TRNY(:,:) ! TRANSPARENCY INFORMATION ON IX,IY
        REAL, POINTER         :: SPCBAC(:,:), ANGARC(:)
        REAL   , POINTER :: XGRD(:,:), YGRD(:,:) ! X AND Y DEFINED ON IX,IY
        REAL   , POINTER :: DXDP(:,:), DXDQ(:,:) ! DX/DP & DX/DQ DEFINED ON IX,IY
        REAL   , POINTER :: DYDP(:,:), DYDQ(:,:) ! DY/DP & DY/DQ DEFINED ON IX,IY
        REAL   , POINTER :: DPDX(:,:), DPDY(:,:) ! DP/DX & DP/DY DEFINED ON IX,IY
        REAL   , POINTER :: DQDX(:,:), DQDY(:,:) ! DQ/DX & DQ/DY DEFINED ON IX,IY
        REAL   , POINTER :: GSQRT(:,:) ! SQRT(G) DEFINED ON IX,IY
        REAL   , POINTER :: HPFAC(:,:) ! H_P = SQRT(G_PP) DEFINED ON IX,IY
        REAL   , POINTER :: HQFAC(:,:) ! H_Q = SQRT(G_QQ) DEFINED ON IX,IY
 
        LOGICAL          :: GINIT, FLDRY, FLCX, FLCY, FLCTH, FLCK, FLSOU, IICEDISP
        LOGICAL          :: FLAGLL
        LOGICAL, POINTER :: FLAGST(:)
        CHARACTER(LEN=30):: GNAME
        CHARACTER(LEN=13):: FILEXT
        LOGICAL          :: GUGINIT
        REAL, POINTER    :: REFLC(:,:)  ! reflection coefficient
        INTEGER, POINTER :: REFLD(:,:)  ! reflection direction
        INTEGER          :: E3DF(3,5), P2MSF(3), US3DF(3) ! freq. indices for 3D output
!
        TYPE(T_GSU) :: GSU ! Grid search utility object
!
        REAL                  :: FFACBERG    ! mutiplicative factor for iceberg mask
!
! unstructured data
!
        INTEGER               :: NTRI
        DOUBLE PRECISION, POINTER         :: XYB(:,:)
        INTEGER, POINTER      :: TRIGP(:,:)
        DOUBLE PRECISION, POINTER         :: LEN(:,:),SI(:), IEN(:,:)
 
        REAL                  :: MAXX, MAXY, DXYMAX
        REAL, POINTER         :: ANGLE(:,:),ANGLE0(:,:)
        INTEGER               :: COUNTRI,COUNTOT,NNZ
        INTEGER, POINTER      :: CCON(:), COUNTCON(:), IE_CELL(:), VNEIGH(:,:),  &
                                 POS_CELL(:), IOBP(:), IOBPD(:,:), IOBPA(:),     &
                                 IAA(:), JAA(:), POSI(:,:), INDEX_CELL(:)
        DOUBLE PRECISION, POINTER         :: TRIA(:)
        REAL, POINTER         :: CROSSDIFF(:,:)
      END TYPE GRID
!
      TYPE SGRD   ! this is the spectral grid with all parameters that vary with freq. and direction
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
      TYPE FLDP
         REAL :: DUMMY
      END TYPE FLDP
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
        INTEGER               :: SSDSISO, SSDSBRFDF
        REAL                  :: AALPHA, BBETA, ZZ0MAX, ZZ0RAT, ZZALP,&
                                 SSINTHP, TTAUWSHELTER, SSWELLF(1:6), &
                                 SSDSC1, SSDSC2, SSDSC3, SSDSBR,      &
                                 SSDSP, WWNMEANP, SSTXFTF, SSTXFTWN,  &
                                 FFXPM, FFXFM,                        &
                                 SSDSC4, SSDSC5, SSDSC6, DDELTA1,     &
                                 DDELTA2, ZZWND
!
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
      TYPE SICP
        REAL                  :: DUMMY
      END TYPE SICP
!
      TYPE SXXP
        REAL                  :: DUMMY
      END TYPE SXXP
 
! specific type for unstructured scheme
      TYPE SCHM
         LOGICAL              :: FSN,FSPSI,FSFCT,FSNIMP
      END TYPE SCHM
!
      TYPE MPAR
        LOGICAL               :: PINIT
        TYPE(NPAR)            :: NPARS
        TYPE(PROP)            :: PROPS
        TYPE(FLDP)            :: FLDPS
        TYPE(SFLP)            :: SFLPS
        TYPE(SLNP)            :: SLNPS
        TYPE(SRCP)            :: SRCPS
        TYPE(SNLP)            :: SNLPS
        TYPE(SBTP)            :: SBTPS
        TYPE(SDBP)            :: SDBPS
        TYPE(STRP)            :: STRPS
        TYPE(SBSP)            :: SBSPS
        TYPE(SICP)            :: SICPS
        TYPE(SXXP)            :: SXXPS
        TYPE(SCHM)            :: SCHMS
      END TYPE MPAR
!/
!/ Data storage
!/
      TYPE(GRID), TARGET, ALLOCATABLE :: GRIDS(:)
      TYPE(SGRD), TARGET, ALLOCATABLE :: SGRDS(:)
      TYPE(MPAR), TARGET, ALLOCATABLE :: MPARS(:)
!/
!/ Data aliases for structure GRID(S)
!/
      INTEGER, POINTER :: GTYPE
      INTEGER, POINTER :: ICLOSE
      INTEGER, POINTER        :: NX, NY, NSEA, NSEAL, TRFLAG
      INTEGER, POINTER        :: E3DF(:,:), P2MSF(:), US3DF(:)
      REAL,    POINTER        :: REFLC(:,:)
      INTEGER, POINTER        :: REFLD(:,:)
!
! Variables for unstructured grids
!
      INTEGER, POINTER        :: NTRI,COUNTRI,COUNTOT,NNZ
!  XYB may not be necessary now that we have XGRD and YGRD
!  but these XGRD and YGRD should probably be double precision
      DOUBLE PRECISION, POINTER  ::     XYB(:,:)
      INTEGER, POINTER        :: TRIGP(:,:)
      DOUBLE PRECISION, POINTER  :: IEN(:,:), LEN(:,:), SI(:)
      REAL, POINTER           :: ANGLE(:,:),ANGLE0(:,:)
      INTEGER,  POINTER       :: CCON(:),  COUNTCON(:), IE_CELL(:),    &
                                 POS_CELL(:), IOBP(:), IOBPD(:,:),     &
                                 IAA(:), JAA(:), POSI(:,:),            &
                                 VNEIGH(:,:), INDEX_CELL(:), IOBPA(:)
      DOUBLE PRECISION, POINTER  :: TRIA(:)
      REAL, POINTER           :: CROSSDIFF(:,:)
      REAL,POINTER            :: MAXX, MAXY, DXYMAX
      LOGICAL, POINTER        :: GUGINIT
!
      REAL,    POINTER        :: FFACBERG
      REAL,    POINTER        :: REFPARS(:)
      INTEGER, POINTER        :: MAPSTA(:,:), MAPST2(:,:),            &
                                 MAPFS(:,:), MAPSF(:,:)
!
      REAL, POINTER           :: SX, SY, X0, Y0, DTCFL, DTCFLI, DTMAX, &
                                 DTMIN, DMIN, CTMAX, FICE0, FICEN,     &
                                 FICEL, PFMOVE, STEXU, STEYU, STEDU,   &
                                 IICEHMIN, IICEHINIT, ICE100WIND
      REAL(8),POINTER         :: GRIDSHIFT ! see notes in WMGHGH
      REAL, POINTER           :: ZB(:), CLATS(:)
      REAL   , POINTER :: CLATIS(:) ! INVERSE OF COS(LAT) DEFINED ON ISEA
      REAL   , POINTER :: CTHG0S(:) ! TAN(Y)/R, DEFINED ON ISEA
 
      REAL   , POINTER :: TRNX(:,:), TRNY(:,:) ! TRANSPARENCY INFORMATION ON IX,IY
      REAL, POINTER         :: SPCBAC(:,:), ANGARC(:)
      REAL   , POINTER :: XGRD(:,:), YGRD(:,:) ! X AND Y DEFINED ON IX,IY
      REAL   , POINTER :: DXDP(:,:), DXDQ(:,:) ! DX/DP & DX/DQ DEFINED ON IX,IY
      REAL   , POINTER :: DYDP(:,:), DYDQ(:,:) ! DY/DP & DY/DQ DEFINED ON IX,IY
      REAL   , POINTER :: DPDX(:,:), DPDY(:,:) ! DP/DX & DP/DY DEFINED ON IX,IY
      REAL   , POINTER :: DQDX(:,:), DQDY(:,:) ! DQ/DX & DQ/DY DEFINED ON IX,IY
      REAL   , POINTER :: GSQRT(:,:) ! SQRT(G) DEFINED ON IX,IY
      REAL   , POINTER :: HPFAC(:,:) ! H_P = SQRT(G_PP) DEFINED ON IX,IY
      REAL   , POINTER :: HQFAC(:,:) ! H_Q = SQRT(G_QQ) DEFINED ON IX,IY
 
      LOGICAL, POINTER :: GINIT, FLDRY, FLCX, FLCY, FLCTH, FLCK, FLSOU, IICEDISP
      LOGICAL, POINTER :: FLAGLL
      LOGICAL, POINTER :: FLAGST(:)
 
      CHARACTER(LEN=30), POINTER :: GNAME
      CHARACTER(LEN=13), POINTER :: FILEXT
 
      TYPE(T_GSU), POINTER :: GSU ! Grid search utility object
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
!/ Data aliasses for structure FLDP(S)
!/
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
      REAL, POINTER           :: ZZWND, AALPHA, BBETA, ZZ0MAX, ZZ0RAT,&
                                 ZZALP, FFXFM, FFXPM,                 &
                                 SSINTHP, TTAUWSHELTER, SSWELLF(:),   &
                                 SSDSC1, SSDSC2, SSDSC3, SSDSBR,      &
                                 SSDSP, WWNMEANP, SSTXFTF, SSTXFTWN,  &
                                 SSDSC4, SSDSC5, SSDSC6, SSDSBR2,     &
                                 DDELTA1, DDELTA2,                    &
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
!/ Data aliasing for structure SCHM(S)
      LOGICAL, POINTER        :: FSN,FSPSI,FSFCT,FSNIMP
!/
!/ Data aliasing for structure SICP(S)
!/
      CONTAINS
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3NMOD ( NUMBER, NDSE, NDST, NAUX )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Dec-2014 !
!/                  +-----------------------------------+
!/
!/    24-Feb-2004 : Origination.                        ( version 3.06 )
!/    18-Jul-2006 : Add input grids.                    ( version 3.10 )
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
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
      ALLOCATE ( GRIDS(NLOW:NUMBER), &
                 SGRDS(NLOW:NUMBER), &
                 MPARS(NLOW:NUMBER), &
                 STAT=ISTAT )
      CHECK_ALLOC_STATUS ( ISTAT )
!
! -------------------------------------------------------------------- /
! 2.  Initialize GINIT and SINIT
!
      DO I=NLOW, NUMBER
        GRIDS(I)%GINIT  = .FALSE.
        GRIDS(I)%GUGINIT  = .FALSE.
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
      SUBROUTINE W3DIMX  ( IMOD, MX, MY, MSEA, NDSE, NDST   &
                         )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           H. L. Tolman            |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         10-Dec-2014 |
!/                  +-----------------------------------+
!/
!/    24-Jun-2005 : Origination.                        ( version 3.07 )
!/    18-Jul-2006 : Add input grids.                    ( version 3.10 )
!/    05-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    02-Feb-2007 : Add FLAGST.                         ( version 3.10 )
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    30-Oct-2009 : Implement curvilinear grid type.    ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    30-Oct-2009 : Implement unstructured grids        ( version 3.14.1)
!/    03-Sep-2012 : Clean up of UG grids                ( version 4.08 )
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
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
      IF ( MX.LT.3 .OR. (MY.LT.3.AND.GTYPE.NE.UNGTYPE) .OR. MSEA.LT.1 ) THEN
        WRITE (NDSE,1003) MX, MY, MSEA, GTYPE
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
! NB: Some array start at 0 because MAPFS(IY,IX)=0 for missing points
!
      ALLOCATE ( GRIDS(IMOD)%MAPSTA(MY,MX),  &
                 GRIDS(IMOD)%MAPST2(MY,MX),  &
                 GRIDS(IMOD)%MAPFS(MY,MX),   &
                 GRIDS(IMOD)%MAPSF(MSEA,3),  &
                 GRIDS(IMOD)%FLAGST(MSEA),   &
                 GRIDS(IMOD)%ZB(MSEA),       &
                 GRIDS(IMOD)%CLATS(0:MSEA),  &
                 GRIDS(IMOD)%CLATIS(0:MSEA), &
                 GRIDS(IMOD)%CTHG0S(0:MSEA), &
                 GRIDS(IMOD)%TRNX(MY,MX),    &
                 GRIDS(IMOD)%TRNY(MY,MX),    &
                 GRIDS(IMOD)%XGRD(MY,MX),    &
                 GRIDS(IMOD)%YGRD(MY,MX),    &
                 GRIDS(IMOD)%DXDP(MY,MX),    &
                 GRIDS(IMOD)%DXDQ(MY,MX),    &
                 GRIDS(IMOD)%DYDP(MY,MX),    &
                 GRIDS(IMOD)%DYDQ(MY,MX),    &
                 GRIDS(IMOD)%DPDX(MY,MX),    &
                 GRIDS(IMOD)%DPDY(MY,MX),    &
                 GRIDS(IMOD)%DQDX(MY,MX),    &
                 GRIDS(IMOD)%DQDY(MY,MX),    &
                 GRIDS(IMOD)%GSQRT(MY,MX),   &
                 GRIDS(IMOD)%HPFAC(MY,MX),   &
                 GRIDS(IMOD)%HQFAC(MY,MX),   &
                 STAT=ISTAT                  )
      CHECK_ALLOC_STATUS ( ISTAT )
!
      GRIDS(IMOD)%FLAGST = .TRUE.
      GRIDS(IMOD)%GINIT  = .TRUE.
      GRIDS(IMOD)%MAPSF(:,3)=0.
      GRIDS(IMOD)%CLATS(0)=1.
      GRIDS(IMOD)%CLATIS(0)=1.
      GRIDS(IMOD)%CTHG0S(0)=1.
!
 
!
! -------------------------------------------------------------------- /
! 2.  Update counters in grid
!
      GRIDS(IMOD)%NX     = MX
      GRIDS(IMOD)%NY     = MY
      GRIDS(IMOD)%NSEA   = MSEA
!
! -------------------------------------------------------------------- /
! 3.  Point to allocated arrays
!
      CALL W3SETG ( IMOD, NDSE, NDST )
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
               '                    INPUT = ',4I10 /)
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
!/                  | Last update :         10-Dec-2014 !
!/                  +-----------------------------------+
!/
!/    19-Feb-2004 : Origination.                        ( version 3.06 )
!/    18-Jul-2006 : Add input grids.                    ( version 3.10 )
!/    05-Oct-2006 : Add filter to array pointers.       ( version 3.10 )
!/    10-Dec-2014 : Add checks for allocate status      ( version 5.04 )
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
                 SGRDS(IMOD)%DDEN2(MSPEC),                            &
                 STAT=ISTAT                                           )
      CHECK_ALLOC_STATUS ( ISTAT )
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
!/                  |                        FORTRAN 90 |
!/                  | Last update :         03-Sep-2012 |
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
!/    30-Oct-2009 : Implement run-time grid selection.  ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    30-Oct-2009 : Implement curvilinear grid type.    ( version 3.14 )
!/                  (W. E. Rogers & T. J. Campbell, NRL)
!/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
!/                  specify index closure for a grid.   ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    13-Jul-2012 : Move data structures GMD (SNL3) and nonlinear
!/                  filter (SNLS) from 3.15 (HLT).      ( version 4.08 )
!/    03-Sep-2012 : Clean up of UG grids                ( version 4.08 )
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
      GTYPE  => GRIDS(IMOD)%GTYPE
      ICLOSE => GRIDS(IMOD)%ICLOSE
!
      NX     => GRIDS(IMOD)%NX
      NY     => GRIDS(IMOD)%NY
      NSEA   => GRIDS(IMOD)%NSEA
      NSEAL  => GRIDS(IMOD)%NSEAL
      TRFLAG => GRIDS(IMOD)%TRFLAG
      FLAGLL => GRIDS(IMOD)%FLAGLL
!
      E3DF   => GRIDS(IMOD)%E3DF
      P2MSF  => GRIDS(IMOD)%P2MSF
      US3DF  => GRIDS(IMOD)%US3DF
      REFLC  => GRIDS(IMOD)%REFLC
      REFLD  => GRIDS(IMOD)%REFLD
      FFACBERG => GRIDS(IMOD)%FFACBERG
      SX     => GRIDS(IMOD)%SX
      SY     => GRIDS(IMOD)%SY
      X0     => GRIDS(IMOD)%X0
      Y0     => GRIDS(IMOD)%Y0
!
      DTCFL  => GRIDS(IMOD)%DTCFL
      DTCFLI => GRIDS(IMOD)%DTCFLI
      DTMAX  => GRIDS(IMOD)%DTMAX
      DTMIN  => GRIDS(IMOD)%DTMIN
      DMIN   => GRIDS(IMOD)%DMIN
      CTMAX  => GRIDS(IMOD)%CTMAX
      FICE0  => GRIDS(IMOD)%FICE0
      GRIDSHIFT  => GRIDS(IMOD)%GRIDSHIFT
      FICEN  => GRIDS(IMOD)%FICEN
      FICEL  => GRIDS(IMOD)%FICEL
      IICEHMIN  => GRIDS(IMOD)%IICEHMIN
      IICEHINIT  => GRIDS(IMOD)%IICEHINIT
      ICE100WIND  => GRIDS(IMOD)%ICE100WIND
      PFMOVE => GRIDS(IMOD)%PFMOVE
      STEXU  => GRIDS(IMOD)%STEXU
      STEYU  => GRIDS(IMOD)%STEYU
      STEDU  => GRIDS(IMOD)%STEDU
!
      GINIT  => GRIDS(IMOD)%GINIT
      GUGINIT  => GRIDS(IMOD)%GUGINIT
      FLDRY  => GRIDS(IMOD)%FLDRY
      FLCX   => GRIDS(IMOD)%FLCX
      FLCY   => GRIDS(IMOD)%FLCY
      FLCTH  => GRIDS(IMOD)%FLCTH
      FLCK   => GRIDS(IMOD)%FLCK
      FLSOU  => GRIDS(IMOD)%FLSOU
      IICEDISP => GRIDS(IMOD)%IICEDISP
!
      GNAME  => GRIDS(IMOD)%GNAME
      FILEXT => GRIDS(IMOD)%FILEXT
      XYB    => GRIDS(IMOD)%XYB
      TRIGP  => GRIDS(IMOD)%TRIGP
      NTRI     => GRIDS(IMOD)%NTRI
      COUNTRI     => GRIDS(IMOD)%COUNTRI
      SI     => GRIDS(IMOD)%SI
      COUNTOT    => GRIDS(IMOD)%COUNTOT
      IEN     => GRIDS(IMOD)%IEN
      LEN     => GRIDS(IMOD)%LEN
      ANGLE     => GRIDS(IMOD)%ANGLE
      ANGLE0     => GRIDS(IMOD)%ANGLE0
      CCON     => GRIDS(IMOD)%CCON
      COUNTCON     => GRIDS(IMOD)%COUNTCON
      VNEIGH       => GRIDS(IMOD)%VNEIGH
      INDEX_CELL  => GRIDS(IMOD)%INDEX_CELL
      IE_CELL     => GRIDS(IMOD)%IE_CELL
      POS_CELL     => GRIDS(IMOD)%POS_CELL
      IOBP     => GRIDS(IMOD)%IOBP
      IAA      => GRIDS(IMOD)%IAA
      JAA      => GRIDS(IMOD)%JAA
      POSI     => GRIDS(IMOD)%POSI
      NNZ      => GRIDS(IMOD)%NNZ
      IOBPD     => GRIDS(IMOD)%IOBPD
      IOBPA    => GRIDS(IMOD)%IOBPA
      TRIA     => GRIDS(IMOD)%TRIA
      CROSSDIFF => GRIDS(IMOD)%CROSSDIFF
      MAXX     => GRIDS(IMOD)%MAXX
      MAXY     => GRIDS(IMOD)%MAXY
      DXYMAX   => GRIDS(IMOD)%DXYMAX
 
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
          CLATS  => GRIDS(IMOD)%CLATS
          CLATIS => GRIDS(IMOD)%CLATIS
          CTHG0S => GRIDS(IMOD)%CTHG0S
          TRNX   => GRIDS(IMOD)%TRNX
          TRNY   => GRIDS(IMOD)%TRNY
!
          XGRD   => GRIDS(IMOD)%XGRD
          YGRD   => GRIDS(IMOD)%YGRD
          DXDP   => GRIDS(IMOD)%DXDP
          DXDQ   => GRIDS(IMOD)%DXDQ
          DYDP   => GRIDS(IMOD)%DYDP
          DYDQ   => GRIDS(IMOD)%DYDQ
          DPDX   => GRIDS(IMOD)%DPDX
          DPDY   => GRIDS(IMOD)%DPDY
          DQDX   => GRIDS(IMOD)%DQDX
          DQDY   => GRIDS(IMOD)%DQDY
          GSQRT  => GRIDS(IMOD)%GSQRT
          HPFAC  => GRIDS(IMOD)%HPFAC
          HQFAC  => GRIDS(IMOD)%HQFAC
!
          GSU  => GRIDS(IMOD)%GSU
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
!     Structure FLDP
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
!
      ZZWND  => MPARS(IMOD)%SRCPS%ZZWND
      AALPHA => MPARS(IMOD)%SRCPS%AALPHA
      BBETA  => MPARS(IMOD)%SRCPS%BBETA
      SSINTHP  => MPARS(IMOD)%SRCPS%SSINTHP
      ZZ0MAX  => MPARS(IMOD)%SRCPS%ZZ0MAX
      ZZ0RAT  => MPARS(IMOD)%SRCPS%ZZ0RAT
      ZZALP  => MPARS(IMOD)%SRCPS%ZZALP
      TTAUWSHELTER  => MPARS(IMOD)%SRCPS%TTAUWSHELTER
      SSWELLF  => MPARS(IMOD)%SRCPS%SSWELLF
      SSDSC1 => MPARS(IMOD)%SRCPS%SSDSC1
      WWNMEANP => MPARS(IMOD)%SRCPS%WWNMEANP
      FFXFM => MPARS(IMOD)%SRCPS%FFXFM
      FFXPM => MPARS(IMOD)%SRCPS%FFXPM
      DDELTA1 => MPARS(IMOD)%SRCPS%DDELTA1
      DDELTA2 => MPARS(IMOD)%SRCPS%DDELTA2
      SSTXFTF => MPARS(IMOD)%SRCPS%SSTXFTF
      SSTXFTWN => MPARS(IMOD)%SRCPS%SSTXFTWN
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
!     Structure SICPS
!
!    Structure SCHM
       FSN => MPARS(IMOD)%SCHMS%FSN
       FSPSI => MPARS(IMOD)%SCHMS%FSPSI
       FSFCT => MPARS(IMOD)%SCHMS%FSFCT
       FSNIMP => MPARS(IMOD)%SCHMS%FSNIMP
 
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
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3GNTX ( IMOD, NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |           T. J. Campbell          |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         20-Jul-2011 |
!/                  +-----------------------------------+
!/
!/    30-Oct-2009 : Origination.                        ( version 3.13 )
!/    06-Dec-2010 : Change from GLOBAL (logical) to ICLOSE (integer) to
!/                  specify index closure for a grid.   ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    23-Dec-2010 : Fix HPFAC and HQFAC by including the COS(YGRD)
!/                  factor with DXDP and DXDQ terms.    ( version 3.14 )
!/                  (T. J. Campbell, NRL)
!/    20-Jul-2011 : HPFAC and HQFAC are now calculated using W3DIST.
!/                  Result should be very similar except near pole.
!/                  Due to precision issues, HPFAC and HQFAC revert
!/                  to SX and SY in case of regular grids.
!/                  (W. E. Rogers, NRL)                 ( version 3.14 )
!/
!  1. Purpose :
!
!     Construct required spatial grid quantities for curvilinear grids.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       IMOD    Int.   I   Model number to point to.
!       NDSE    Int.   I   Error output unit number.
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
!     - Check on previous initialization of grids.
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
      INTEGER, INTENT(IN)     :: IMOD, NDSE, NDST
!/
!/ ------------------------------------------------------------------- /
!/ Local parameters
!/
      REAL(8), PARAMETER :: PI = 3.14159265358979323846D0
      REAL(8), PARAMETER :: D360 = 360D0
      REAL(8), PARAMETER :: D180 = 180D0
      REAL(8), PARAMETER :: D2R = PI/D180
      REAL(8), PARAMETER :: R2D = 1D0/D2R
      REAL(8), PARAMETER :: HALF = 0.5D0
      REAL(8)            :: DXCIRC(3)
      INTEGER            :: IX, IY
!     REAL               :: UVECPX(NY,NX), UVECPY(NY,NX)
!     REAL               :: UVECQX(NY,NX), UVECQY(NY,NX)
 
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
! 2.  Create grid search utility object
!
      GRIDS(IMOD)%GSU = W3GSUC( .FALSE., FLAGLL, GRIDS(IMOD)%ICLOSE, &
                                GRIDS(IMOD)%NX, GRIDS(IMOD)%NY, &
                                GRIDS(IMOD)%XGRD, GRIDS(IMOD)%YGRD )
 
!
! -------------------------------------------------------------------- /
! 3.  Reset grid pointers
!
      CALL W3SETG ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 4.  Construct curvilinear grid derivatives and metric
!
!  x = x(p,q);  y = y(p,q);  dp = dq = 1;
!
!  When using spherical coordinates x = longitude and y = latitude.
!  Also, the derivatives with respect to p are adjusted by the spherical
!  coordinate metric.  In other words (with radius = 1),
!
!      dx/dp <= cos(y)*(dx/dp)  and  dy/dp <= cos(y)*(dy/dp).
!
!  The 2D curvilinear coordinates metric is
!
!      g_pq = (dx/dp)*(dx/dq) + (dy/dp)*(dy/dq),
!      g_pp = (dx/dp)*(dx/dp) + (dy/dp)*(dy/dp),
!      g_qq = (dx/dq)*(dx/dq) + (dy/dq)*(dy/dq).
!
!  Derivatives are computed with central differences for all non-boundary
!  cells.  One-sided differences are used at boundary cells.  If the grid
!  is a lat/lon grid, then the i-differences at the i=1 and i=nx boundary
!  cells are computed with central differences using the wrapped cell
!  x-coordinates.
!
!  Admittedly, in case of branch cut, we could use centered difference
!  to get a more accurate derivative at the boundaries (IX=1 or NX), but
!  we don't do this. Note that it would get more complicated for a
!  tripole grid (IY=NY) ....
 
!
!  Compute dx/dp & dy/dp & dm/dp
 
      IX = 1
      IF(FLAGLL)THEN
         DO IY=1,NY
            DXCIRC(1)=       ( XGRD(IY ,IX+1) - XGRD(IY ,IX  ) )-D360
            DXCIRC(2)=       ( XGRD(IY ,IX+1) - XGRD(IY ,IX  ) )
            DXCIRC(3)=       ( XGRD(IY ,IX+1) - XGRD(IY ,IX  ) )+D360
            DXDP(IY,IX)=DXCIRC(MINLOC(ABS(DXCIRC),1))
         END DO
      ELSE
         DXDP(: ,IX) =        ( XGRD(: ,IX+1) - XGRD(: ,IX  ) )
      ENDIF
      DYDP(: ,IX) =        ( YGRD(: ,IX+1) - YGRD(: ,IX  ) )
      DO IY=1,NY
         HPFAC(IY,IX)=W3DIST(FLAGLL,XGRD(IY,IX+1),YGRD(IY,IX+1), &
              XGRD(IY,IX),YGRD(IY,IX))
      END DO
 
      IX = NX
      IF(FLAGLL)THEN
         DO IY=1,NY
            DXCIRC(1)=       ( XGRD(IY ,IX  ) - XGRD(IY ,IX-1) )-D360
            DXCIRC(2)=       ( XGRD(IY ,IX  ) - XGRD(IY ,IX-1) )
            DXCIRC(3)=       ( XGRD(IY ,IX  ) - XGRD(IY ,IX-1) )+D360
            DXDP(IY,IX)=DXCIRC(MINLOC(ABS(DXCIRC),1))
         END DO
      ELSE
         DXDP(: ,IX) =        ( XGRD(: ,IX  ) - XGRD(: ,IX-1) )
      ENDIF
      DYDP(: ,IX) =        ( YGRD(: ,IX  ) - YGRD(: ,IX-1) )
      DO IY=1,NY
         HPFAC(IY,IX)=W3DIST(FLAGLL,XGRD(IY,IX),YGRD(IY,IX), &
              XGRD(IY,IX-1),YGRD(IY,IX-1))
      END DO
 
      DO IX=2, NX-1
 
         IF(FLAGLL)THEN
            DO IY=1,NY
               DXCIRC(1)= ( XGRD(IY ,IX+1) - XGRD(IY ,IX-1) )-D360
               DXCIRC(2)= ( XGRD(IY ,IX+1) - XGRD(IY ,IX-1) )
               DXCIRC(3)= ( XGRD(IY ,IX+1) - XGRD(IY ,IX-1) )+D360
               DXDP(IY,IX)=HALF * DXCIRC(MINLOC(ABS(DXCIRC),1))
            END DO
         ELSE
            DXDP(: ,IX)   = HALF * ( XGRD(: ,IX+1) - XGRD(: ,IX-1) )
         ENDIF
 
         DYDP(: ,IX)   = HALF * ( YGRD(: ,IX+1) - YGRD(: ,IX-1) )
         DO IY=1,NY
            HPFAC(IY,IX)=W3DIST(FLAGLL,XGRD(IY,IX+1),YGRD(IY,IX+1), &
                 XGRD(IY,IX-1),YGRD(IY,IX-1))*HALF
         END DO
 
      END DO
 
!
!  Compute dx/dq & dy/dq & dm/dq
 
      IY = 1
      IF(FLAGLL)THEN
         DO IX=1,NX
            DXCIRC(1)=       ( XGRD(IY+1,IX) - XGRD(IY  ,IX) )-D360
            DXCIRC(2)=       ( XGRD(IY+1,IX) - XGRD(IY  ,IX) )
            DXCIRC(3)=       ( XGRD(IY+1,IX) - XGRD(IY  ,IX) )+D360
            DXDQ(IY,IX)=DXCIRC(MINLOC(ABS(DXCIRC),1))
         END DO
      ELSE
         DXDQ(IY,: )   =        ( XGRD(IY+1,: ) - XGRD(IY  ,: ) )
      ENDIF
      DYDQ(IY,: )   =        ( YGRD(IY+1,: ) - YGRD(IY  ,: ) )
      DO IX=1,NX
         HQFAC(IY,IX)=W3DIST(FLAGLL,XGRD(IY+1,IX),YGRD(IY+1,IX), &
              XGRD(IY,IX),YGRD(IY,IX))
      END DO
 
      IY = NY
      IF(FLAGLL)THEN
         DO IX=1,NX
            DXCIRC(1)=       ( XGRD(IY  ,IX) - XGRD(IY-1,IX) )-D360
            DXCIRC(2)=       ( XGRD(IY  ,IX) - XGRD(IY-1,IX) )
            DXCIRC(3)=       ( XGRD(IY  ,IX) - XGRD(IY-1,IX) )+D360
            DXDQ(IY,IX)=DXCIRC(MINLOC(ABS(DXCIRC),1))
         END DO
      ELSE
         DXDQ(IY,: )   =        ( XGRD(IY  ,: ) - XGRD(IY-1,: ) )
      ENDIF
      DYDQ(IY,: )   =        ( YGRD(IY  ,: ) - YGRD(IY-1,: ) )
      DO IX=1,NX
         HQFAC(IY,IX)=W3DIST(FLAGLL,XGRD(IY,IX),YGRD(IY,IX), &
              XGRD(IY-1,IX),YGRD(IY-1,IX))
      END DO
 
      DO IY=2, NY-1
         IF(FLAGLL)THEN
            DO IX=1,NX
               DXCIRC(1)= ( XGRD(IY+1,IX) - XGRD(IY-1,IX) )-D360
               DXCIRC(2)= ( XGRD(IY+1,IX) - XGRD(IY-1,IX) )
               DXCIRC(3)= ( XGRD(IY+1,IX) - XGRD(IY-1,IX) )+D360
               DXDQ(IY,IX)=HALF * DXCIRC(MINLOC(ABS(DXCIRC),1))
            END DO
         ELSE
            DXDQ(IY,: ) = HALF * ( XGRD(IY+1,: ) - XGRD(IY-1,: ) )
         ENDIF
         DYDQ(IY,: ) = HALF * ( YGRD(IY+1,: ) - YGRD(IY-1,: ) )
         DO IX=1,NX
            HQFAC(IY,IX)=W3DIST(FLAGLL,XGRD(IY+1,IX),YGRD(IY+1,IX), &
                 XGRD(IY-1,IX),YGRD(IY-1,IX))*HALF
         END DO
      END DO
 
!......notes: Though we *should* recover SX and SY using calculations on XGRD and YGRD,
!.............it turns out that there are small differences which cause problems for
!.............SCRIP. Therefore, we add the if-then here to use SX and SY in case of
!.............regular grids, to overwrite DXDP etc. calculated above.
!.............Note 1: that this may be unnecessary if we go to double precision.
!.............Note 2: HQFAC calculated above is still used. Thus, we don't have the entire
!.....................block of code above inside an "if not regular grid"
 
      IF(GRIDS(IMOD)%GTYPE.EQ.RLGTYPE)THEN
         DXDP=SX
         DXDQ=0.0D0
         DYDP=0.0D0
         DYDQ=SY
      ENDIF   !    IF(GRIDS(IMOD)%GTYPE.EQ.RLGTYPE)THEN
 
!
!  Account for spherical metric if lat/lon grid (only if the metric
!  is to be w.r.t. the Cartesian space)
!     IF ( FLAGLL ) THEN
!         DXDP = DXDP * COS(YGRD*D2R)
!         DXDQ = DXDQ * COS(YGRD*D2R)
!       END IF
!  If we do this, then the COS(YGRD*D2R) factor must be removed from
!  the HPFAC & HQFAC.
!
!  GSQRT: sqrt(g) = (dx/dp)(dy/dq) - (dx/dq)(dy/dp)
      GSQRT = DXDP*DYDQ - DXDQ*DYDP
!
!  Compute curvilinear derivatives
!  DPDX: dp/dx =  (1/sqrt(g))*(dy/dq)
!  DPDY: dp/dy = -(1/sqrt(g))*(dx/dq)
!  DQDX: dq/dx = -(1/sqrt(g))*(dy/dp)
!  DQDY: dq/dy =  (1/sqrt(g))*(dx/dp)
      DPDX =  DYDQ / GSQRT
      DPDY = -DXDQ / GSQRT
      DQDX = -DYDP / GSQRT
      DQDY =  DXDP / GSQRT
!
!  Check orthogonality of grid by computing angle between the
!  curvilinear coordinate unit vectors.
!    cos(a) = uvec_p \dot uvec_q = g_pq/sqrt(g_pp*g_qq)
!
!  Compute cell unit vectors
!  uvec_p = (1/h_p)*(dx/dp)*uvec_x + (1/h_p)*(dy/dp)*uvec_y
!  uvec_q = (1/h_q)*(dx/dq)*uvec_x + (1/h_q)*(dy/dq)*uvec_y
!     UVECPX = DXDP / HPFAC
!     UVECPY = DYDP / HPFAC
!     UVECQX = DXDQ / HQFAC
!     UVECQY = DYDQ / HQFAC
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3GNTX : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3GNTX : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NAUXGR = ',I10/                   &
               '                    NGRIDS = ',I10/)
!
!/
!/ End of W3GNTX ----------------------------------------------------- /
!/
      END SUBROUTINE W3GNTX
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3DIMUG  ( IMOD, MTRI, MX, COUNTOTA, NNZ, NDSE, NDST )
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH-III           NOAA/NCEP |
!/                  |             F.ardhuin             |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         15-Mar-2007 !
!/                  +-----------------------------------+
!/
!/    15-Mar-2007 : Origination.                        ( version 3.14 )
!/    11-May-2015 : Updates to 2-ways nestings for UG   ( version 5.08 )
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
!       MX, MTRI, MSEA       Like NX, NTRI, NSEA in data structure.
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
      INTEGER, INTENT(IN)     :: IMOD, MTRI, MX, COUNTOTA, NNZ, NDSE, NDST
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
          WRITE (NDSE,1002) IMOD, NGRIDS
          CALL EXTCDE (2)
        END IF
      IF ( GRIDS(IMOD)%GUGINIT ) THEN
        WRITE (NDSE,1004)
        CALL EXTCDE (4)
        END IF
!
! -------------------------------------------------------------------- /
! 2.  Allocate arrays
!
      ALLOCATE ( GRIDS(IMOD)%TRIGP(MTRI,3),                         &
                 GRIDS(IMOD)%XYB(MX,3),                             &
                 GRIDS(IMOD)%SI(MX),                                &
                 GRIDS(IMOD)%TRIA(MTRI),                            &
                 GRIDS(IMOD)%CROSSDIFF(6,MTRI),                     &
                 GRIDS(IMOD)%IEN(MTRI,6),                           &
                 GRIDS(IMOD)%LEN(MTRI,3),                           &
                 GRIDS(IMOD)%ANGLE(MTRI,3),                         &
                 GRIDS(IMOD)%ANGLE0(MTRI,3),                        &
                 GRIDS(IMOD)%CCON(MX),                              &
                 GRIDS(IMOD)%COUNTCON(MX),                          &
                 GRIDS(IMOD)%VNEIGH(MX,2*COUNTRI),                  &
                 GRIDS(IMOD)%INDEX_CELL(MX+1),                      &
                 GRIDS(IMOD)%IE_CELL(COUNTOTA),                     &
                 GRIDS(IMOD)%POS_CELL(COUNTOTA),                    &
                 GRIDS(IMOD)%IAA(NX+1),                             &
                 GRIDS(IMOD)%JAA(NNZ),                              &
                 GRIDS(IMOD)%POSI(3,COUNTOTA),                      &
                 GRIDS(IMOD)%IOBP(MX),                              &
                 GRIDS(IMOD)%IOBPD(NTH,MX),                         &
                 GRIDS(IMOD)%IOBPA(MX),                             &
                 STAT=ISTAT                                         )
      CHECK_ALLOC_STATUS ( ISTAT )
!
                 GRIDS(IMOD)%IOBP(:)=1
!
!some segmentation troubles can appear, they are related with the allocation of
!normal(1st dimension) and the nesting of the triangulated grid.
! -------------------------------------------------------------------- /
! 3.  Point to allocated arrays
!
      CALL W3SETG ( IMOD, NDSE, NDST )
!
! -------------------------------------------------------------------- /
! 4.  Update counters in grid
!
      NTRI   = MTRI
      COUNTOT=COUNTOTA
      GUGINIT  = .FALSE.
      RETURN
!
! Formats
!
 1001 FORMAT (/' *** ERROR W3DIMUG : GRIDS NOT INITIALIZED *** '/      &
               '                    RUN W3NMOD FIRST '/)
 1002 FORMAT (/' *** ERROR W3DIMUG : ILLEGAL MODEL NUMBER *** '/       &
               '                    IMOD   = ',I10/                   &
               '                    NGRIDS = ',I10/)
 1003 FORMAT (/' *** ERROR W3DIMUG : ILLEGAL GRID DIMENSION(S) *** '/  &
               '                    INPUT = ',I10/)
 1004 FORMAT (/' *** ERROR W3DIMUG : ARRAY(S) ALREADY ALLOCATED *** ')
!
!/
!/ End of W3DIMUG ----------------------------------------------------- /
!/
      END SUBROUTINE W3DIMUG
!/ ------------------------------------------------------------------- /
      SUBROUTINE W3SETREF
!/
!/                  +-----------------------------------+
!/                  | WAVEWATCH III           NOAA/NCEP |
!/                  |           F. Ardhuin              |
!/                  |                        FORTRAN 90 |
!/                  | Last update :         13-Nov-2013 |
!/                  +-----------------------------------+
!/
!/    13-Nov-2013 : Origination.                        ( version 4.13 )
!/
!  1. Purpose :
!
!     Update reflection directions at shoreline.
!
!  2. Method :
!
!  3. Parameters :
!
!     Parameter list
!     ----------------------------------------------------------------
!       None
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
!      WW3_GRID  Prog. WW3_GRID Grid preprocessor
!      W3ULEV    Subr. W3UPDTMD Water level update
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
!       !/S      Enable subroutine tracing.
!       !/T      Basic test output.
!
! 10. Source code :
!
!/ ------------------------------------------------------------------- /
      USE CONSTANTS
!
      IMPLICIT NONE
!/
!/ ------------------------------------------------------------------- /
!/
      INTEGER                 :: ISEA, IX, IY, IXY, IXN, IXP, IYN, IYP
      INTEGER                 :: J, K, NEIGH1(0:7)
 
      REAL                    :: TRIX(NY*NX), TRIY(NY*NX), DX, DY,    &
                                 COSAVG, SINAVG, THAVG, ANGLES(0:7), CLAT
!/
!/ ------------------------------------------------------------------- /
!/
!
! 1.  Preparations --------------------------------------------------- *
!
      RETURN
!
! Formats
!
!/
!/ End of W3SETREF ----------------------------------------------------- /
!/
      END SUBROUTINE W3SETREF
 
!/
!/ End of module W3GDATMD -------------------------------------------- /
!/
      END MODULE W3GDATMD
