module micro_mg1_5

!---------------------------------------------------------------------------------
! Purpose:
!   MG microphysics version 1.5 - Update of MG microphysics and jumping-off
!       point for the development of MG2
!
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Peter Caldwell, Xiaohong Liu and Steve Ghan
! Version 2: Development begun: September 2011
! invoked in CAM by specifying -microphys=mg1.5
!
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!
! NOTE: Modified to allow other microphysics packages (e.g. CARMA) to do ice
! microphysics in cooperation with the MG liquid microphysics. This is
! controlled by the do_cldice variable.
!
! NOTE: If do_cldice is false, then MG microphysics should not update CLDICE
! or NUMICE; however, it is assumed that the other microphysics scheme will have
! updated CLDICE and NUMICE. The other microphysics should handle the following
! processes that would have been done by MG:
!   - Detrainment (liquid and ice)
!   - Homogeneous ice nucleation
!   - Heterogeneous ice nucleation
!   - Bergeron process
!   - Melting of ice
!   - Freezing of cloud drops
!   - Autoconversion (ice -> snow)
!   - Growth/Sublimation of ice
!   - Sedimentation of ice
!---------------------------------------------------------------------------------
! Based on micro_mg (restructuring of former cldwat2m_micro)
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------
! Code comments added by HM, 093011
! General code structure:
!
! Code is divided into two main subroutines:
!   subroutine micro_mg_init --> initializes microphysics routine, should be called
!                                  once at start of simulation
!   subroutine micro_mg_tend --> main microphysics routine to be called each time step
!                                this also calls several smaller subroutines to calculate
!                                microphysical processes and other utilities
!
! List of external functions:
!   qsat_water --> for calculating saturation vapor pressure with respect to liquid water
!   qsat_ice --> for calculating saturation vapor pressure with respect to ice
!   gamma   --> standard mathematical gamma function
! .........................................................................
! List of inputs through use statement in fortran90:
! Variable Name                      Description                Units
! .........................................................................
! gravit          acceleration due to gravity                    m s-2
! rair            dry air gas constant for air                  J kg-1 K-1
! tmelt           temperature of melting point for water          K
! cpair           specific heat at constant pressure for dry air J kg-1 K-1
! rh2o            gas constant for water vapor                  J kg-1 K-1
! latvap          latent heat of vaporization                   J kg-1
! latice          latent heat of fusion                         J kg-1
! qsat_water      external function for calculating liquid water
!                 saturation vapor pressure/humidity              -
! qsat_ice        external function for calculating ice
!                 saturation vapor pressure/humidity              pa
! rhmini          relative humidity threshold parameter for
!                 nucleating ice                                  -
! .........................................................................
! NOTE: List of all inputs/outputs passed through the call/subroutine statement
!       for micro_mg_tend is given below at the start of subroutine micro_mg_tend.
!---------------------------------------------------------------------------------

! Procedures required:
! 1) An implementation of the gamma function (if not intrinsic).
! 2) saturation vapor pressure and specific humidity over water
! 3) svp over ice

#ifndef HAVE_GAMMA_INTRINSICS
use shr_spfn_mod, only: gamma => shr_spfn_gamma
#endif

use wv_sat_methods, only: &
     qsat_water => wv_sat_qsat_water, &
     qsat_ice => wv_sat_qsat_ice

implicit none
private
save

public :: &
     micro_mg_init, &
     micro_mg_get_cols, &
     micro_mg_tend

! switch for specification rather than prediction of droplet and crystal number
! note: number will be adjusted as needed to keep mean size within bounds,
! even when specified droplet or ice number is used

! If constant cloud ice number is set (nicons = .true.),
! then all microphysical processes except mass transfer due to ice nucleation
! (mnuccd) are based on the fixed cloud ice number. Calculation of
! mnuccd follows from the prognosed ice crystal number ni.

! nccons = .true. to specify constant cloud droplet number
! nicons = .true. to specify constant cloud ice number

logical, parameter, public :: nccons = .false.
logical, parameter, public :: nicons = .false.

!=========================================================
! Private module parameters
!=========================================================

integer, parameter :: r8 = selected_real_kind(12)   ! 8 byte real

real(r8), parameter :: pi = 3.14159265358979323846_r8 ! To 20 digits; more than enough
                                                          ! to reach the limit of double precision

real(r8), parameter :: omsm   = 0.99999_r8    ! number near unity for round-off issues

! parameters for specified ice and droplet number concentration
! note: these are local in-cloud values, not grid-mean
real(r8), parameter :: ncnst = 100.e6_r8    ! droplet num concentration when nccons=.true. (m-3)
real(r8), parameter :: ninst = 0.1e6_r8     ! ice num concentration when nicons=.true. (m-3)

! rhow used to be passed in in init, but then overwritten.
! For now, just setting it as a parameter here
real(r8), parameter :: rhosn = 250._r8  ! bulk density snow
real(r8), parameter :: rhoi = 500._r8   ! bulk density ice
real(r8), parameter :: rhow = 1000._r8  ! bulk density liquid

! fall speed parameters, V = aD^b (V is in m/s)
! droplets
real(r8), parameter :: ac = 3.e7_r8
real(r8), parameter :: bc = 2._r8
! snow
real(r8), parameter :: as = 11.72_r8
real(r8), parameter :: bs = 0.41_r8
! cloud ice
real(r8), parameter :: ai = 700._r8
real(r8), parameter :: bi = 1._r8
! rain
real(r8), parameter :: ar = 841.99667_r8
real(r8), parameter :: br = 0.8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d
! exponent
real(r8), parameter :: dsph = 3._r8

! ventilation parameters
! for snow
real(r8), parameter :: f1s = 0.86_r8
real(r8), parameter :: f2s = 0.28_r8
! for rain
real(r8), parameter :: f1r = 0.78_r8
real(r8), parameter :: f2r = 0.308_r8

! collection efficiencies
! aggregation of cloud ice and snow
real(r8), parameter :: eii = 0.1_r8

! autoconversion size threshold for cloud ice to snow (m)
real(r8) :: dcs 

!!== KZ_DCS 
logical :: dcs_tdep 
!!== KZ_DCS 

! smallest mixing ratio considered in microphysics
real(r8), parameter :: qsmall = 1.e-18_r8

! alternate threshold used for some in-cloud mmr
real(r8), parameter :: icsmall = 1.e-8_r8

! immersion freezing parameters, bigg 1953
real(r8), parameter :: bimm = 100._r8
real(r8), parameter :: aimm = 0.66_r8

! mass of new crystal due to aerosol freezing and growth (kg)
real(r8), parameter :: mi0 = 4._r8/3._r8*pi*rhoi*(10.e-6_r8)**3

! minimum mass of new crystal due to freezing of cloud droplets done
! externally (kg)
real(r8), parameter :: mi0l_min = 4._r8/3._r8*pi*rhow*(4.e-6_r8)**3

!Range of cloudsat reflectivities (dBz) for analytic simulator
real(r8), parameter :: csmin = -30._r8
real(r8), parameter :: csmax = 26._r8
real(r8), parameter :: mindbz = -99._r8
real(r8), parameter :: minrefl = 1.26e-10_r8    ! minrefl = 10._r8**(mindbz/10._r8)

!=========================================================
! Constants set in initialization
!=========================================================

! Set using arguments to micro_mg_init
real(r8) :: g           ! gravity
real(r8) :: r           ! dry air gas constant
real(r8) :: rv          ! water vapor gas constant
real(r8) :: cpp         ! specific heat of dry air
real(r8) :: tmelt       ! freezing point of water (K)

! latent heats of:
real(r8) :: xxlv        ! vaporization
real(r8) :: xlf         ! freezing
real(r8) :: xxls        ! sublimation

real(r8) :: rhmini      ! Minimum rh for ice cloud fraction > 0.

! flags
logical :: microp_uniform
logical :: do_cldice
logical :: use_hetfrz_classnuc

real(r8) :: rhosu       ! typical 850mn air density

real(r8) :: icenuct     ! ice nucleation temperature: currently -5 degrees C

real(r8) :: snowmelt    ! what temp to melt all snow: currently 2 degrees C
real(r8) :: rainfrze    ! what temp to freeze all rain: currently -5 degrees C

! additional constants to help speed up code
real(r8) :: cons1
real(r8) :: cons4
real(r8) :: cons5
real(r8) :: cons7
real(r8) :: cons8
real(r8) :: cons11
real(r8) :: cons13
real(r8) :: cons14
real(r8) :: cons16
real(r8) :: cons17
real(r8) :: cons22
real(r8) :: cons23
real(r8) :: cons24
real(r8) :: cons25
real(r8) :: cons27
real(r8) :: cons28

character(len=16)  :: micro_mg_precip_frac_method  ! type of precipitation fraction method
real(r8)           :: micro_mg_berg_eff_factor     ! berg efficiency factor


!===============================================================================
contains
!===============================================================================

subroutine micro_mg_init( &
     kind, gravit, rair, rh2o, cpair,    &
     tmelt_in, latvap, latice,           &
!!== KZ_DCS 
     rhmini_in, micro_mg_dcs, micro_mg_dcs_tdep, microp_uniform_in, do_cldice_in, use_hetfrz_classnuc_in, &
!!== KZ_DCS 
     micro_mg_precip_frac_method_in, micro_mg_berg_eff_factor_in, errstring)

  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! initialize constants for MG microphysics
  !
  ! Author: Andrew Gettelman Dec 2005
  !
  !-----------------------------------------------------------------------

  integer,  intent(in)  :: kind         ! Kind used for reals
  real(r8), intent(in)  :: gravit
  real(r8), intent(in)  :: rair
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: tmelt_in     ! Freezing point of water (K)
  real(r8), intent(in)  :: latvap
  real(r8), intent(in)  :: latice
  real(r8), intent(in)  :: rhmini_in    ! Minimum rh for ice cloud fraction > 0.
  real(r8), intent(in)  :: micro_mg_dcs
!!== KZ_DCS 
  logical,  intent(in)  :: micro_mg_dcs_tdep 
!!== KZ_DCS 

  logical,  intent(in)  :: microp_uniform_in    ! .true. = configure for sub-columns
                                            ! .false. = use w/o sub-columns (standard)
  logical,  intent(in)  :: do_cldice_in     ! .true. = do all processes (standard)
                                            ! .false. = skip all processes affecting
                                            !           cloud ice
  logical,  intent(in)  :: use_hetfrz_classnuc_in ! use heterogeneous freezing

  character(len=16),intent(in)  :: micro_mg_precip_frac_method_in  ! type of precipitation fraction method
  real(r8),         intent(in)  :: micro_mg_berg_eff_factor_in     ! berg efficiency factor


  character(128), intent(out) :: errstring    ! Output status (non-blank for error return)


  !-----------------------------------------------------------------------

  dcs = micro_mg_dcs

  dcs_tdep = micro_mg_dcs_tdep

  errstring = ' '

  if( kind .ne. r8 ) then
     errstring = 'micro_mg_init: KIND of reals does not match'
     return
  endif

  ! declarations for MG code (transforms variable names)

  g= gravit                 ! gravity
  r= rair                   ! dry air gas constant: note units(phys_constants are in J/K/kmol)
  rv= rh2o                  ! water vapor gas constant
  cpp = cpair               ! specific heat of dry air
  tmelt = tmelt_in
  rhmini = rhmini_in
  micro_mg_precip_frac_method = micro_mg_precip_frac_method_in
  micro_mg_berg_eff_factor    = micro_mg_berg_eff_factor_in


  ! latent heats

  xxlv = latvap         ! latent heat vaporization
  xlf  = latice         ! latent heat freezing
  xxls = xxlv + xlf     ! latent heat of sublimation

  ! flags
  microp_uniform = microp_uniform_in
  do_cldice  = do_cldice_in
  use_hetfrz_classnuc = use_hetfrz_classnuc_in

  ! typical air density at 850 mb

  rhosu = 85000._r8/(rair * tmelt)

  ! Maximum temperature at which snow is allowed to exist
  snowmelt = tmelt + 2._r8
  ! Minimum temperature at which rain is allowed to exist
  rainfrze = tmelt - 5._r8

  ! Ice nucleation temperature
  icenuct  = tmelt - 5._r8

  ! Define constants to help speed up code (this limits calls to gamma function)
  ! Unused names: cons6, cons15, cons21, cons26
  cons1=gamma(1._r8+dsph)
  cons4=gamma(1._r8+br)
  cons5=gamma(4._r8+br)
  cons7=gamma(1._r8+bs)
  cons8=gamma(4._r8+bs)
  cons11=gamma(3._r8+bs)
  cons13=gamma(5._r8/2._r8+br/2._r8)
  cons14=gamma(5._r8/2._r8+bs/2._r8)
  cons16=gamma(1._r8+bi)
  cons17=gamma(4._r8+bi)
  cons22=(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
  cons23=dcs**3
  cons24=dcs**2
  cons25=dcs**bs
  cons27=xxlv**2
  cons28=xxls**2

end subroutine micro_mg_init

!===============================================================================
!microphysics routine for each timestep goes here...

subroutine micro_mg_tend ( &
     mgncol,   nlev,     deltatin,           &
     t,                 q,                                     &
     qcn,                          qin,                          &
     ncn,                          nin,                          &
     relvar,            accre_enhan,                           &
     p,                 pdel,              pint,               &
     cldn,               liqcldf,            icecldf,            &
     qcsinksum_rate1ord,  naai,    npccn,  rndst,   nacon,  &
     tlat,    qvlat,   qctend,  qitend,  nctend,  nitend,  &
     effc,    effc_fn, effi,              prect,   preci,   &
     nevapr, evapsnow, prain,  prodsnow, cmeout,  deffi,   &
     pgamrad, lamcrad, qsout,   dsout,   rflx,    sflx,    &
     qrout,             reff_rain,         reff_snow,         &
     qcsevap, qisevap, qvres,   cmeitot,  vtrmc,   vtrmi,   &
     qcsedten,qisedten,pratot,     prctot,     mnuccctot,  mnuccttot,  &
     msacwitot,  psacwstot,  bergstot,   bergtot,    melttot,    homotot,    &
     qcrestot,             prcitot,    praitot,    qirestot,             &
     mnuccrtot,  pracstot,   meltsdttot, frzrdttot,  mnuccdtot,            &
     nrout,   nsout,   refl,    arefl,   areflz,  frefl,   &
     csrfl,   acsrfl,  fcsrfl,            rercld,            &
     ncai,    ncal,    qrout2,  qsout2,  nrout2,  nsout2,  &
     drout2,  dsout2,  freqs,   freqr,   nfice,   qcrat,   &
     errstring, &
     tnd_qsnow,         tnd_nsnow,         re_ice,            &
     prer_evap, &
     frzimm,   frzcnt,   frzdep)

  !Authors: Hugh Morrison, Andrew Gettelman, NCAR, Peter Caldwell, LLNL
  ! e-mail: morrison@ucar.edu, andrew@ucar.edu

  ! input arguments
  integer,  intent(in) :: mgncol                ! number of microphysics columns
  integer,  intent(in) :: nlev                  ! number of layers
  real(r8), intent(in) :: deltatin              ! time step (s)
  real(r8), intent(in) :: t(:,:)               ! input temperature (K)
  real(r8), intent(in) :: q(:,:)               ! input h20 vapor mixing ratio (kg/kg)
  real(r8), intent(in) :: relvar(:,:)          ! relative variance of cloud water (-)
  real(r8), intent(in) :: accre_enhan(:,:)     ! optional accretion enhancement factor (-)

  ! note: all input cloud variables are grid-averaged
  real(r8), intent(in) :: qcn(:,:)       ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(:,:)       ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: ncn(:,:)       ! cloud water number conc (1/kg)
  real(r8), intent(in) :: nin(:,:)       ! cloud ice number conc (1/kg)
  real(r8), intent(in) :: p(:,:)         ! air pressure (pa)
  real(r8), intent(in) :: pdel(:,:)      ! pressure difference across level (pa)
  ! hm add 11-16-11, interface pressure
  real(r8), intent(in) :: pint(:,:)    ! level interface pressure (pa)
  real(r8), intent(in) :: cldn(:,:)      ! cloud fraction (no units)
  real(r8), intent(in) :: liqcldf(:,:)   ! liquid cloud fraction (no units)
  real(r8), intent(in) :: icecldf(:,:)   ! ice cloud fraction (no units)
  ! used for scavenging
  ! Inputs for aerosol activation
  real(r8), intent(in) :: naai(:,:)     ! ice nucleation number (from microp_aero_ts) (1/kg)
  real(r8), intent(in) :: npccn(:,:)   ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)

  ! Note that for these variables, the dust bin is assumed to be the last index.
  ! (For example, in CAM, the last dimension is always size 4.)
  real(r8), intent(in) :: rndst(:,:,:)  ! radius of each dust bin, for contact freezing (from microp_aero_ts) (m)
  real(r8), intent(in) :: nacon(:,:,:) ! number in each dust bin, for contact freezing  (from microp_aero_ts) (1/m^3)

  ! output arguments

  real(r8), intent(out) :: qcsinksum_rate1ord(:,:)    ! 1st order rate for
  ! direct cw to precip conversion
  real(r8), intent(out) :: tlat(:,:)         ! latent heating rate       (W/kg)
  real(r8), intent(out) :: qvlat(:,:)        ! microphysical tendency qv (1/s)
  real(r8), intent(out) :: qctend(:,:)       ! microphysical tendency qc (1/s)
  real(r8), intent(out) :: qitend(:,:)       ! microphysical tendency qi (1/s)
  real(r8), intent(out) :: nctend(:,:)       ! microphysical tendency nc (1/(kg*s))
  real(r8), intent(out) :: nitend(:,:)       ! microphysical tendency ni (1/(kg*s))
  real(r8), intent(out) :: effc(:,:)         ! droplet effective radius (micron)
  real(r8), intent(out) :: effc_fn(:,:)      ! droplet effective radius, assuming nc = 1.e8 kg-1
  real(r8), intent(out) :: effi(:,:)         ! cloud ice effective radius (micron)
  real(r8), intent(out) :: prect(:)          ! surface precip rate (m/s)
  real(r8), intent(out) :: preci(:)          ! cloud ice/snow precip rate (m/s)
  real(r8), intent(out) :: nevapr(:,:)       ! evaporation rate of rain + snow (1/s)
  real(r8), intent(out) :: evapsnow(:,:)     ! sublimation rate of snow (1/s)
  real(r8), intent(out) :: prain(:,:)        ! production of rain + snow (1/s)
  real(r8), intent(out) :: prodsnow(:,:)     ! production of snow (1/s)
  real(r8), intent(out) :: cmeout(:,:)       ! evap/sub of cloud (1/s)
  real(r8), intent(out) :: deffi(:,:)        ! ice effective diameter for optics (radiation) (micron)
  real(r8), intent(out) :: pgamrad(:,:)      ! ice gamma parameter for optics (radiation) (no units)
  real(r8), intent(out) :: lamcrad(:,:)      ! slope of droplet distribution for optics (radiation) (1/m)
  real(r8), intent(out) :: qsout(:,:)        ! snow mixing ratio (kg/kg)
  real(r8), intent(out) :: dsout(:,:)        ! snow diameter (m)
  real(r8), intent(out) :: rflx(:,:)         ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8), intent(out) :: sflx(:,:)         ! grid-box average snow flux (kg m^-2 s^-1)
  real(r8), intent(out) :: qrout(:,:)        ! grid-box average rain mixing ratio (kg/kg)
  real(r8), intent(out) :: reff_rain(:,:)    ! rain effective radius (micron)
  real(r8), intent(out) :: reff_snow(:,:)    ! snow effective radius (micron)
  real(r8), intent(out) :: qcsevap(:,:)      ! cloud water evaporation due to sedimentation (1/s)
  real(r8), intent(out) :: qisevap(:,:)      ! cloud ice sublimation due to sublimation (1/s)
  real(r8), intent(out) :: qvres(:,:)        ! residual condensation term to ensure RH < 100% (1/s)
  real(r8), intent(out) :: cmeitot(:,:)       ! grid-mean cloud ice sub/dep (1/s)
  real(r8), intent(out) :: vtrmc(:,:)        ! mass-weighted cloud water fallspeed (m/s)
  real(r8), intent(out) :: vtrmi(:,:)        ! mass-weighted cloud ice fallspeed (m/s)
  real(r8), intent(out) :: qcsedten(:,:)     ! qc sedimentation tendency (1/s)
  real(r8), intent(out) :: qisedten(:,:)     ! qi sedimentation tendency (1/s)

  ! microphysical process rates for output (mixing ratio tendencies) (all have units of 1/s)
  real(r8), intent(out) :: pratot(:,:)         ! accretion of cloud by rain
  real(r8), intent(out) :: prctot(:,:)         ! autoconversion of cloud to rain
  real(r8), intent(out) :: mnuccctot(:,:)      ! mixing ratio tend due to immersion freezing
  real(r8), intent(out) :: mnuccttot(:,:)      ! mixing ratio tend due to contact freezing
  real(r8), intent(out) :: msacwitot(:,:)      ! mixing ratio tend due to H-M splintering
  real(r8), intent(out) :: psacwstot(:,:)      ! collection of cloud water by snow
  real(r8), intent(out) :: bergstot(:,:)       ! bergeron process on snow
  real(r8), intent(out) :: bergtot(:,:)        ! bergeron process on cloud ice
  real(r8), intent(out) :: melttot(:,:)        ! melting of cloud ice
  real(r8), intent(out) :: homotot(:,:)        ! homogeneous freezing cloud water
  real(r8), intent(out) :: qcrestot(:,:)       ! residual cloud condensation due to removal of excess supersat
  real(r8), intent(out) :: prcitot(:,:)        ! autoconversion of cloud ice to snow
  real(r8), intent(out) :: praitot(:,:)        ! accretion of cloud ice by snow
  real(r8), intent(out) :: qirestot(:,:)       ! residual ice deposition due to removal of excess supersat
  real(r8), intent(out) :: mnuccrtot(:,:)      ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
  real(r8), intent(out) :: pracstot(:,:)       ! mixing ratio tendency due to accretion of rain by snow (1/s)
  real(r8), intent(out) :: meltsdttot(:,:)     ! latent heating rate due to melting of snow  (W/kg)
  real(r8), intent(out) :: frzrdttot(:,:)      ! latent heating rate due to homogeneous freezing of rain (W/kg)
  real(r8), intent(out) :: mnuccdtot(:,:)      ! mass tendency from ice nucleation
  real(r8), intent(out) :: nrout(:,:)        ! rain number concentration (1/m3)
  real(r8), intent(out) :: nsout(:,:)        ! snow number concentration (1/m3)
  real(r8), intent(out) :: refl(:,:)         ! analytic radar reflectivity
  real(r8), intent(out) :: arefl(:,:)        ! average reflectivity will zero points outside valid range
  real(r8), intent(out) :: areflz(:,:)       ! average reflectivity in z.
  real(r8), intent(out) :: frefl(:,:)        ! fractional occurrence of radar reflectivity
  real(r8), intent(out) :: csrfl(:,:)        ! cloudsat reflectivity
  real(r8), intent(out) :: acsrfl(:,:)       ! cloudsat average
  real(r8), intent(out) :: fcsrfl(:,:)       ! cloudsat fractional occurrence of radar reflectivity
  real(r8), intent(out) :: rercld(:,:)       ! effective radius calculation for rain + cloud
  real(r8), intent(out) :: ncai(:,:)        ! output number conc of ice nuclei available (1/m3)
  real(r8), intent(out) :: ncal(:,:)        ! output number conc of CCN (1/m3)
  real(r8), intent(out) :: qrout2(:,:)       ! copy of qrout as used to compute drout2
  real(r8), intent(out) :: qsout2(:,:)       ! copy of qsout as used to compute dsout2
  real(r8), intent(out) :: nrout2(:,:)       ! copy of nrout as used to compute drout2
  real(r8), intent(out) :: nsout2(:,:)       ! copy of nsout as used to compute dsout2
  real(r8), intent(out) :: drout2(:,:)       ! mean rain particle diameter (m)
  real(r8), intent(out) :: dsout2(:,:)       ! mean snow particle diameter (m)
  real(r8), intent(out) :: freqs(:,:)        ! fractional occurrence of snow
  real(r8), intent(out) :: freqr(:,:)        ! fractional occurrence of rain
  real(r8), intent(out) :: nfice(:,:)        ! fractional occurrence of ice
  real(r8), intent(out) :: qcrat(:,:)        ! limiter for qc process rates (1=no limit --> 0. no qc)

  real(r8), intent(out) :: prer_evap(:,:)

  character(128),   intent(out) :: errstring  ! output status (non-blank for error return)

  ! Tendencies calculated by external schemes that can replace MG's native
  ! process tendencies.

  ! Used with CARMA cirrus microphysics
  ! (or similar external microphysics model)
  real(r8), intent(in), pointer :: tnd_qsnow(:,:) ! snow mass tendency (kg/kg/s)
  real(r8), intent(in), pointer :: tnd_nsnow(:,:) ! snow number tendency (#/kg/s)
  real(r8), intent(in), pointer :: re_ice(:,:)    ! ice effective radius (m)

  ! From external ice nucleation.
  real(r8), intent(in), pointer :: frzimm(:,:) ! Number tendency due to immersion freezing (1/cm3)
  real(r8), intent(in), pointer :: frzcnt(:,:) ! Number tendency due to contact freezing (1/cm3)
  real(r8), intent(in), pointer :: frzdep(:,:) ! Number tendency due to deposition nucleation (1/cm3)

  ! local workspace
  ! all units mks unless otherwise stated

  ! parameters
  real(r8), parameter :: mincld = 0.0001_r8     ! minimum allowed cloud fraction
  real(r8), parameter :: cdnl   = 0.e6_r8       ! cloud droplet number limiter

  ! local copies of input variables
  real(r8) :: qc(mgncol,nlev)          ! cloud liquid mixing ratio (kg/kg)
  real(r8) :: qi(mgncol,nlev)          ! cloud ice mixing ratio (kg/kg)
  real(r8) :: nc(mgncol,nlev)          ! cloud liquid number concentration (1/kg)
  real(r8) :: ni(mgncol,nlev)          ! cloud liquid number concentration (1/kg)

  real(r8) :: dcst(mgncol,nlev)        ! t-dependent dcs

  real(r8) :: nevapr2(mgncol,nlev) 

  ! general purpose variables
  real(r8) :: deltat            ! sub-time step (s)
  real(r8) :: mtime             ! the assumed ice nucleation timescale

  real(r8) :: dz(mgncol,nlev)         ! height difference across model vertical level
  real(r8) :: int_to_mid(mgncol, nlev) ! Coefficients for linear interpolation from
                                       ! interface to mid-level

  ! physical properties of the air at a given point
  real(r8) :: rho(mgncol,nlev)    ! density (kg m-3)
  real(r8) :: dv(mgncol,nlev)     ! diffusivity of water vapor
  real(r8) :: mu(mgncol,nlev)     ! viscosity
  real(r8) :: sc(mgncol,nlev)     ! schmidt number
  real(r8) :: rhof(mgncol,nlev)   ! density correction factor for fallspeed

  ! cloud fractions
  real(r8) :: cldmax(mgncol,nlev) ! precip fraction assuming maximum overlap
  real(r8) :: cldm(mgncol,nlev)   ! cloud fraction
  real(r8) :: icldm(mgncol,nlev)  ! ice cloud fraction
  real(r8) :: lcldm(mgncol,nlev)  ! liq cloud fraction

  ! mass mixing ratios
  real(r8) :: qcic(mgncol,nlev)   ! in-cloud cloud liquid
  real(r8) :: qiic(mgncol,nlev)   ! in-cloud cloud ice
  real(r8) :: qsic(mgncol,nlev)   ! in-precip snow
  real(r8) :: qric(mgncol,nlev)   ! in-precip rain

  ! number concentrations
  real(r8) :: ncic(mgncol,nlev)   ! in-cloud droplet
  real(r8) :: niic(mgncol,nlev)   ! in-cloud cloud ice
  real(r8) :: nsic(mgncol,nlev)   ! in-precip snow
  real(r8) :: nric(mgncol,nlev)   ! in-precip rain
  ! maximum allowed ni value
  real(r8) :: nimax(mgncol,nlev)

  ! Size distribution parameters for:
  ! cloud ice
  real(r8) :: lami(mgncol,nlev)   ! slope
  real(r8) :: n0i(mgncol,nlev)    ! intercept
  ! cloud liquid
  real(r8) :: lamc(mgncol,nlev)   ! slope
  real(r8) :: pgam(mgncol,nlev)   ! spectral width parameter
  real(r8) :: cdist1(mgncol,nlev) ! droplet freezing calculation
  ! snow
  real(r8) :: lams(mgncol,nlev)   ! slope
  real(r8) :: n0s(mgncol,nlev)    ! intercept
  ! rain
  real(r8) :: lamr(mgncol,nlev)   ! slope
  real(r8) :: n0r(mgncol,nlev)    ! intercept

  ! Rates/tendencies due to:
  ! deposition of cloud ice
  real(r8) :: vap_dep(mgncol,nlev)    ! deposition from vapor to ice PMC 12/3/12
  ! sublimation of cloud ice
  real(r8) :: ice_sublim(mgncol,nlev) ! sublimation from ice to vapor PMC 12/3/12
  ! ice nucleation
  real(r8) :: nnuccd(mgncol,nlev) ! number rate from deposition/cond.-freezing
  real(r8) :: mnuccd(mgncol,nlev) ! mass mixing ratio
  ! freezing of cloud water
  real(r8) :: mnuccc(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnuccc(mgncol,nlev) ! number concentration
  ! contact freezing of cloud water
  real(r8) :: mnucct(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnucct(mgncol,nlev) ! number concentration
  ! deposition nucleation in mixed-phase clouds (from external scheme)
  real(r8) :: mnudep(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnudep(mgncol,nlev) ! number concentration
  ! HM ice multiplication
  real(r8) :: msacwi(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nsacwi(mgncol,nlev) ! number concentration
  ! autoconversion of cloud droplets
  real(r8) :: prc(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: nprc(mgncol,nlev)   ! number concentration (rain)
  real(r8) :: nprc1(mgncol,nlev)  ! number concentration (cloud droplets)
  ! self-aggregation of snow
  real(r8) :: nsagg(mgncol,nlev)  ! number concentration
  ! self-collection of rain
  real(r8) :: nragg(mgncol,nlev)  ! number concentration
  ! collection of droplets by snow
  real(r8) :: psacws(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: npsacws(mgncol,nlev)    ! number concentration
  ! collection of rain by snow
  real(r8) :: pracs(mgncol,nlev)  ! mass mixing ratio
  real(r8) :: npracs(mgncol,nlev) ! number concentration
  ! freezing of rain
  real(r8) :: mnuccr(mgncol,nlev) ! mass mixing ratio
  real(r8) :: nnuccr(mgncol,nlev) ! number concentration
  ! accretion of droplets by rain
  real(r8) :: pra(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: npra(mgncol,nlev)   ! number concentration
  ! autoconversion of cloud ice to snow
  real(r8) :: prci(mgncol,nlev)   ! mass mixing ratio
  real(r8) :: nprci(mgncol,nlev)  ! number concentration
  ! accretion of cloud ice by snow
  real(r8) :: prai(mgncol,nlev)   ! mass mixing ratio
  real(r8) :: nprai(mgncol,nlev)  ! number concentration
  ! evaporation of rain
  real(r8) :: pre(mgncol,nlev)    ! mass mixing ratio
  ! sublimation of snow
  real(r8) :: prds(mgncol,nlev)   ! mass mixing ratio
  ! number evaporation
  real(r8) :: nsubi(mgncol,nlev)  ! cloud ice
  real(r8) :: nsubc(mgncol,nlev)  ! droplet
  real(r8) :: nsubs(mgncol,nlev)  ! snow
  real(r8) :: nsubr(mgncol,nlev)  ! rain
  ! bergeron process
  real(r8) :: berg(mgncol,nlev)   ! mass mixing ratio (cloud ice)
  real(r8) :: bergs(mgncol,nlev)  ! mass mixing ratio (snow)

  ! fallspeeds
  ! number-weighted
  real(r8) :: uns(mgncol,nlev)    ! snow
  real(r8) :: unr(mgncol,nlev)    ! rain
  ! mass-weighted
  real(r8) :: ums(mgncol,nlev)    ! snow
  real(r8) :: umr(mgncol,nlev)    ! rain
  ! air density corrected fallspeed parameters
  real(r8) :: arn(mgncol,nlev)    ! rain
  real(r8) :: asn(mgncol,nlev)    ! snow
  real(r8) :: acn(mgncol,nlev)    ! cloud droplet
  real(r8) :: ain(mgncol,nlev)    ! cloud ice

  ! Mass of liquid droplets used with external heterogeneous freezing.
  real(r8) :: mi0l(mgncol)

  ! saturation vapor pressures
  real(r8) :: esl(mgncol,nlev)    ! liquid
  real(r8) :: esi(mgncol,nlev)    ! ice
  real(r8) :: esn               ! checking for RH after rain evap

  ! saturation vapor mixing ratios
  real(r8) :: qvl(mgncol,nlev)        ! liquid
  real(r8) :: qvi(mgncol,nlev)        ! ice
  real(r8) :: qvn                   ! checking for RH after rain evap

  ! relative humidity
  real(r8) :: relhum(mgncol,nlev)

  ! parameters for cloud water and cloud ice sedimentation calculations
  real(r8) :: fc(nlev)
  real(r8) :: fnc(nlev)
  real(r8) :: fi(nlev)
  real(r8) :: fni(nlev)

  real(r8) :: faloutc(nlev)
  real(r8) :: faloutnc(nlev)
  real(r8) :: falouti(nlev)
  real(r8) :: faloutni(nlev)

  real(r8) :: faltndc
  real(r8) :: faltndnc
  real(r8) :: faltndi
  real(r8) :: faltndni
  real(r8) :: faltndqie
  real(r8) :: faltndqce

  ! sum of source/sink terms for diagnostic precip
  real(r8) :: qstend(mgncol,nlev)    ! snow mixing ratio
  real(r8) :: nstend(mgncol,nlev)     ! snow number concentration
  real(r8) :: qrtend(mgncol,nlev)     ! rain mixing ratio
  real(r8) :: nrtend(mgncol,nlev)     ! rain number concentration
  ! vertically integrated source/sink terms
  real(r8) :: qrtot(mgncol)           ! rain mixing ratio
  real(r8) :: nrtot(mgncol)           ! rain number concentration
  real(r8) :: qstot(mgncol)           ! snow mixing ratio
  real(r8) :: nstot(mgncol)           ! snow number concentration

  ! for calculation of rate1ord
  real(r8) :: qcsum_rate1ord(mgncol,nlev)     ! sum over iterations of cloud water

  real(r8) :: rainrt(mgncol,nlev)    ! rain rate for reflectivity calculation

  ! dummy variables
  real(r8) :: dum
  real(r8) :: dum1
  real(r8) :: dum2
  ! dummies for checking RH
  real(r8) :: qtmp
  real(r8) :: ttmp
  real(r8) :: qtmp1
  real(r8) :: ttmp1
  ! dummies for conservation check
  real(r8) :: ratio
  real(r8) :: tmpfrz
  ! dummies for in-cloud variables
  real(r8) :: dumc(mgncol,nlev)   ! qc
  real(r8) :: dumnc(mgncol,nlev)  ! nc
  real(r8) :: dumi(mgncol,nlev)   ! qi
  real(r8) :: dumni(mgncol,nlev)  ! ni
  real(r8) :: dumr(mgncol,nlev)   ! rain mixing ratio
  real(r8) :: dumnr(mgncol,nlev)  ! rain number concentration
  ! Array dummy variable
  real(r8) :: dum_2D(mgncol,nlev)

  ! loop array variables
  ! "i" and "k" are column/level iterators for internal (MG) variables
  ! "n" is used for other looping (currently just sedimentation)
  integer i, k, n

  ! number of sub-steps for loops over "n" (for sedimentation)
  integer nstep

  ! Whether or not to limit evaporation/sublimation of precip at each grid
  ! point.
  logical :: limit_precip_evap_sublim

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! default return error message
  errstring = ' '

  if (.not. (do_cldice .or. &
       (associated(tnd_qsnow) .and. associated(tnd_nsnow) .and. associated(re_ice)))) then
     errstring = "MG's native cloud ice processes are disabled, but &
          &no replacement values were passed in."
  end if

  if (use_hetfrz_classnuc .and. (.not. &
       (associated(frzimm) .and. associated(frzcnt) .and. associated(frzdep)))) then
     errstring = "External heterogeneous freezing is enabled, but the &
          &required tendencies were not all passed in."
  end if

  ! Process inputs

  ! assign variable deltat to deltatin
  deltat = deltatin

  ! Copies of input concentrations that may be changed internally.
  qc = qcn
  nc = ncn
  qi = qin
  ni = nin

  ! pint: used to set int_to_mid
  ! interface to mid-level linear interpolation
  do k = 1,nlev
     int_to_mid(:,k) = (p(:,k) - pint(:,k)) / (pint(:,k+1) - pint(:,k))
  end do

  ! cldn: used to set cldm, unused for subcolumns
  ! liqcldf: used to set lcldm, unused for subcolumns
  ! icecldf: used to set icldm, unused for subcolumns

  if (microp_uniform) then
     ! subcolumns, set cloud fraction variables to one
     ! if cloud water or ice is present, if not present
     ! set to mincld (mincld used instead of zero, to prevent
     ! possible division by zero errors).

     where (qc >= qsmall)
        lcldm = 1._r8
     elsewhere
        lcldm = mincld
     end where

     where (qi >= qsmall)
        icldm = 1._r8
     elsewhere
        icldm = mincld
     end where

     cldm = max(icldm, lcldm)

  else
     ! get cloud fraction, check for minimum
     cldm = max(cldn,mincld)
     lcldm = max(liqcldf,mincld)
     icldm = max(icecldf,mincld)
  end if

  ! Initialize local variables

  ! local physical properties
  rho = p/(r*t)
  dv = 8.794E-5_r8 * t**1.81_r8 / p
  mu = 1.496E-6_r8 * t**1.5_r8 / (t + 120._r8)
  sc = mu/(rho*dv)

  ! get dz from dp and hydrostatic approx
  ! keep dz positive (define as layer k-1 - layer k)
  dz = pdel/(rho*g)

  ! air density adjustment for fallspeed parameters
  ! includes air density correction factor to the
  ! power of 0.54 following Heymsfield and Bansemer 2007

  rhof=(rhosu/rho)**0.54_r8

  arn=ar*rhof
  asn=as*rhof
  acn=g*rhow/(18._r8*mu)
  ain=ai*(rhosu/rho)**0.35_r8

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Get humidity and saturation vapor pressures

  do k=1,nlev
     do i=1,mgncol

        call qsat_water(t(i,k), p(i,k), esl(i,k), qvl(i,k))

        ! hm fix, make sure when above freezing that esi=esl, not active yet
        if (t(i,k) >= tmelt) then
           esi(i,k)=esl(i,k)
           qvi(i,k)=qvl(i,k)
        else
           call qsat_ice(t(i,k), p(i,k), esi(i,k), qvi(i,k))
        end if

     end do
  end do

  relhum = q / max(qvl, qsmall)

  !===============================================

  ! hm, set mtime here to avoid answer-changing
  mtime=deltat

  ! initialize microphysics output
  qcsevap=0._r8
  qisevap=0._r8
  qvres  =0._r8
  cmeitot =0._r8
  vtrmc =0._r8
  vtrmi =0._r8
  qcsedten =0._r8
  qisedten =0._r8

  pratot=0._r8
  prctot=0._r8
  mnuccctot=0._r8
  mnuccttot=0._r8
  msacwitot=0._r8
  psacwstot=0._r8
  bergstot=0._r8
  bergtot=0._r8
  melttot=0._r8
  homotot=0._r8
  qcrestot=0._r8
  prcitot=0._r8
  praitot=0._r8
  qirestot=0._r8
  mnuccrtot=0._r8
  pracstot=0._r8
  meltsdttot=0._r8
  frzrdttot=0._r8
  mnuccdtot=0._r8

  rflx=0._r8
  sflx=0._r8

  ! initialize precip output

  qrout=0._r8
  qsout=0._r8
  nrout=0._r8
  nsout=0._r8

  ! for refl calc
  rainrt = 0._r8

  ! initialize rain size
  rercld=0._r8

  qcsinksum_rate1ord = 0._r8
  qcsum_rate1ord     = 0._r8

  ! initialize variables for trop_mozart
  nevapr = 0._r8
  nevapr2 = 0._r8
  evapsnow = 0._r8
  prain = 0._r8
  prodsnow = 0._r8
  cmeout = 0._r8

  cldmax = mincld

  lamc=0._r8

  ! initialize microphysical tendencies

  tlat=0._r8
  qvlat=0._r8
  qctend=0._r8
  qitend=0._r8
  qstend = 0._r8
  qrtend = 0._r8
  nctend=0._r8
  nitend=0._r8
  nrtend = 0._r8
  nstend = 0._r8

  ! initialize diagnostic precipitation to zero
  qcic  = 0._r8
  qiic  = 0._r8
  qsic  = 0._r8
  qric  = 0._r8

  ncic  = 0._r8
  niic  = 0._r8
  nsic  = 0._r8
  nric  = 0._r8

  ! initialize precip at surface

  prect = 0._r8
  preci = 0._r8

  ! initialize vertically-integrated rain and snow tendencies

  qrtot = 0._r8
  nrtot = 0._r8
  qstot = 0._r8
  nstot = 0._r8

  ! initialize precip fallspeeds to zero
  ums = 0._r8
  uns = 0._r8
  umr = 0._r8
  unr = 0._r8

  ! initialize limiter for output
  qcrat = 1._r8


!!== KZ_DCS
  call get_dcst(mgncol,nlev,t,dcst)
!!== KZ_DCS

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! droplet activation
  ! hm, modify 5/12/11
  ! get provisional droplet number after activation. This is used for
  ! all microphysical process calculations, for consistency with update of
  ! droplet mass before microphysics

  ! calculate potential for droplet activation if cloud water is present
  ! tendency from activation (npccn) is read in from companion routine

  ! output activated liquid and ice (convert from #/kg -> #/m3)
  !--------------------------------------------------
  where (qc >= qsmall)
     nc = max(nc + npccn*deltat, cdnl*lcldm/rho)
     ncal = nc*rho/lcldm ! sghan minimum in #/cm3
  elsewhere
     ncal = 0._r8
  end where

  where (t < icenuct)
     ncai = naai*rho
  elsewhere
     ncai = 0._r8
  end where

  !===============================================

  ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
  !-------------------------------------------------------

  if (do_cldice) then
     where (naai > 0._r8 .and. t < icenuct .and. &
       relhum*esl/esi > rhmini+0.05_r8)

        !if NAAI > 0. then set numice = naai (as before)
        !note: this is gridbox averaged
        ! hm, modify to use mtime
        nnuccd = (naai-ni/icldm)/mtime*icldm
        nnuccd = max(nnuccd,0._r8)
        nimax = naai*icldm

        !Calc mass of new particles using new crystal mass...
        !also this will be multiplied by mtime as nnuccd is...

        mnuccd = nnuccd * mi0

     elsewhere
        nnuccd = 0._r8
        nimax = 0._r8
        mnuccd = 0._r8
     end where

  end if

  !=============================================================================
     pre_vert_loop: do k=1,nlev

        pre_col_loop: do i=1,mgncol

           ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
           !-------------------------------------------------------
           ! for microphysical process calculations
           ! units are kg/kg for mixing ratio, 1/kg for number conc

           if (qc(i,k).ge.qsmall) then
              ! limit in-cloud values to 0.005 kg/kg
              qcic(i,k)=min(qc(i,k)/lcldm(i,k),5.e-3_r8)
              ncic(i,k)=max(nc(i,k)/lcldm(i,k),0._r8)

              ! hm add 6/2/11 specify droplet concentration
              if (nccons) then
                 ncic(i,k)=ncnst/rho(i,k)
              end if
           else
              qcic(i,k)=0._r8
              ncic(i,k)=0._r8
           end if

           if (qi(i,k).ge.qsmall) then
              ! limit in-cloud values to 0.005 kg/kg
              qiic(i,k)=min(qi(i,k)/icldm(i,k),5.e-3_r8)
              niic(i,k)=max(ni(i,k)/icldm(i,k),0._r8)

              ! hm add 6/2/11 switch for specification of cloud ice number
              if (nicons) then
                 niic(i,k)=ninst/rho(i,k)
              end if
           else
              qiic(i,k)=0._r8
              niic(i,k)=0._r8
           end if

        end do pre_col_loop
     end do pre_vert_loop

  !========================================================================

     ! for sub-columns cldm has already been set to 1 if cloud
     ! water or ice is present, so cldmax will be correctly set below
     ! and nothing extra needs to be done here

     cldmax = cldm

     micro_vert_loop: do k=1,nlev

        if (trim(micro_mg_precip_frac_method) == 'in_cloud') then 

           do i=1, mgncol

              if (k .eq. 1) then
                 cldmax(i,k)=cldm(i,k)
              else
      
                 if (qc(i,k).ge.qsmall.or.qi(i,k).ge.qsmall) then
                    cldmax(i,k)=cldm(i,k)
                 else
                    cldmax(i,k)=cldmax(i,k-1)
                 end if

              endif

           enddo
        else if(trim(micro_mg_precip_frac_method) == 'max_overlap') then

           ! calculate precip fraction based on maximum overlap assumption
   
           ! if rain or snow mix ratios are smaller than threshold,
           ! then leave cldmax as cloud fraction at current level
           if (k /= 1) then
              where (qric(:,k-1).ge.qsmall .or. qsic(:,k-1).ge.qsmall)
                 cldmax(:,k)=max(cldmax(:,k-1),cldmax(:,k))
              end where
           end if

        endif

        do i = 1, mgncol

           !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           ! get size distribution parameters based on in-cloud cloud water
           ! these calculations also ensure consistency between number and mixing ratio
           !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           ! cloud liquid
           !-------------------------------------------

           ! ".true." below turns on adjustment of ncic for consistency.
           call size_dist_param_liq(qcic(i,k), ncic(i,k), cdnl, rho(i,k), .true., &
                pgam(i,k), lamc(i,k))

           if (lamc(i,k) > 0._r8) then

              ! parameter to calculate droplet freezing
              cdist1(i,k) = ncic(i,k)/gamma(pgam(i,k)+1._r8)

           else
              cdist1(i,k) = 0._r8
           end if

        end do

        !========================================================================
        ! autoconversion of cloud liquid water to rain
        ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
        ! minimum qc of 1 x 10^-8 prevents floating point error

        call kk2000_liq_autoconversion(qcic(:,k), ncic(:,k), rho(:,k), &
             relvar(:,k), prc(:,k), nprc(:,k), nprc1(:,k))

        ! add autoconversion to precip from above to get provisional rain mixing ratio
        ! and number concentration (qric and nric)

        ! hm 11-16-11, modify, divide dz by 2 to use modified mid-point method
        ! This estimates rain and snow mass and number mixing ratios at
        ! mid-point to calculate process rates at mid-point, with final
        ! values of rain and snow mass and number mixing ratios calculated
        ! on interfaces

        if (k .eq. 1) then
           dum=0.45_r8

           qric(:,k)= prc(:,k)*lcldm(:,k)*dz(:,k)/2._r8/cldmax(:,k)/dum
           nric(:,k)=nprc(:,k)*lcldm(:,k)*dz(:,k)/2._r8/cldmax(:,k)/dum
        else

           ! no autoconversion of rain number if rain/snow falling from above
           ! this assumes that new drizzle drops formed by autoconversion are rapidly collected
           ! by the existing rain/snow particles from above

           where (qric(:,k-1).ge.1.e-9_r8 .or. qsic(:,k-1).ge.1.e-9_r8)
              nprc(:,k) = 0._r8
           end where

           do i = 1,mgncol
              if (qric(i,k-1).ge.qsmall) then
                 dum=umr(i,k-1)
                 dum1=unr(i,k-1)
              else
                 ! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)

                 dum=0.45_r8
                 dum1=0.45_r8
              end if

              qric(i,k) = (rho(i,k-1)*umr(i,k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
                   (rho(i,k)*dz(i,k)/2._r8*((pra(i,k-1)+prc(i,k))*lcldm(i,k)+ &
                   (pre(i,k-1)-pracs(i,k-1)-mnuccr(i,k-1))*cldmax(i,k)))) &
                   /(dum*rho(i,k)*cldmax(i,k))
              nric(i,k) = (rho(i,k-1)*unr(i,k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
                   (rho(i,k)*dz(i,k)/2._r8*(nprc(i,k)*lcldm(i,k)+ &
                   (nsubr(i,k-1)-npracs(i,k-1)-nnuccr(i,k-1)+nragg(i,k-1))*cldmax(i,k))))&
                   /(dum1*rho(i,k)*cldmax(i,k))

           end do

        end if

        ! if precip mix ratio is zero so should number concentration

        where (qric(:,k).lt.qsmall)
           qric(:,k)=0._r8
           nric(:,k)=0._r8
        end where

        ! make sure number concentration is a positive number to avoid
        ! taking root of negative later

        nric(:,k)=max(nric(:,k),0._r8)

        ! Get size distribution parameters for cloud ice

!!== KZ_DCS
        call size_dist_param_ice(qiic(:,k), dcst(:,k), niic(:,k), lami(:,k), n0i(:,k))
!!== KZ_DCS

        !.......................................................................
        ! Autoconversion of cloud ice to snow
        ! similar to Ferrier (1994)

        if (do_cldice) then
!!== KZ_DCS
           call ice_autoconversion(t(:,k), qiic(:,k), lami(:,k), n0i(:,k), dcst(:,k),  &
!!== KZ_DCS
                prci(:,k), nprci(:,k))
        else
           ! Add in the particles that we have already converted to snow, and
           ! don't do any further autoconversion of ice.
           prci(:,k)  = tnd_qsnow(:,k) / cldm(:,k)
           nprci(:,k) = tnd_nsnow(:,k) / cldm(:,k)
        end if

        do i=1,mgncol

           ! add autoconversion to flux from level above to get provisional snow mixing ratio
           ! and number concentration (qsic and nsic)

           ! hm 11-16-11 modify for mid-point-type method, see comments above

           if (k == 1) then
              if(dcs_tdep) then 
                 dum=(asn(i,k)*dcst(i,k)**bs)
              else
                 dum=(asn(i,k)*cons25)
              end if 
              qsic(i,k)=prci(i,k)*icldm(i,k)*dz(i,k)/2._r8/cldmax(i,k)/dum
              nsic(i,k)=nprci(i,k)*icldm(i,k)*dz(i,k)/2._r8/cldmax(i,k)/dum
           else
              if (qsic(i,k-1) >= qsmall) then
                 dum=ums(i,k-1)
                 dum1=uns(i,k-1)
              else
              if(dcs_tdep) then 
                 dum = asn(i,k)*dcst(i,k)**bs
              else
                 dum = asn(i,k)*cons25
              end if 
                 dum1 = dum
              end if

              qsic(i,k) = (rho(i,k-1)*ums(i,k-1)*qsic(i,k-1)*cldmax(i,k-1)+ &
                   (rho(i,k)*dz(i,k)/2._r8*((prci(i,k)+prai(i,k-1)+psacws(i,k-1)+bergs(i,k-1))*icldm(i,k)+ &
                   (prds(i,k-1)+pracs(i,k-1)+mnuccr(i,k-1))*cldmax(i,k))))&
                   /(dum*rho(i,k)*cldmax(i,k))

              nsic(i,k) = (rho(i,k-1)*uns(i,k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
                   (rho(i,k)*dz(i,k)/2._r8*(nprci(i,k)*icldm(i,k)+ &
                   (nsubs(i,k-1)+nsagg(i,k-1)+nnuccr(i,k-1))*cldmax(i,k)))) &
                   /(dum1*rho(i,k)*cldmax(i,k))

           end if

        end do

        ! if precip mix ratio is zero so should number concentration

        where (qsic(:,k) < qsmall)
           qsic(:,k)=0._r8
           nsic(:,k)=0._r8
        end where

        ! make sure number concentration is a positive number to avoid
        ! taking root of negative later

        nsic(:,k)=max(nsic(:,k),0._r8)

        !.......................................................................
        ! get size distribution parameters for precip
        !......................................................................
        ! rain

        call size_dist_param_rain(qric(:,k), nric(:,k), lamr(:,k), n0r(:,k))

        where (lamr(:,k) >= qsmall)

           ! provisional rain number and mass weighted mean fallspeed (m/s)

           unr(:,k) = min(arn(:,k)*cons4/lamr(:,k)**br,9.1_r8*rhof(:,k))
           umr(:,k) = min(arn(:,k)*cons5/(6._r8*lamr(:,k)**br),9.1_r8*rhof(:,k))

        elsewhere
           umr(:,k) = 0._r8
           unr(:,k) = 0._r8
        end where

        !......................................................................
        ! snow

        call size_dist_param_snow(qsic(:,k), nsic(:,k), lams(:,k), n0s(:,k))

        where (lams(:,k) > 0._r8)

           ! provisional snow number and mass weighted mean fallspeed (m/s)

           ums(:,k) = min(asn(:,k)*cons8/(6._r8*lams(:,k)**bs),1.2_r8*rhof(:,k))
           uns(:,k) = min(asn(:,k)*cons7/lams(:,k)**bs,1.2_r8*rhof(:,k))

        elsewhere
           ums(:,k) = 0._r8
           uns(:,k) = 0._r8
        end where

        if (do_cldice) then
           if (.not. use_hetfrz_classnuc) then

              ! heterogeneous freezing of cloud water
              !----------------------------------------------

              call immersion_freezing(t(:,k), pgam(:,k), lamc(:,k), cdist1(:,k), qcic(:,k), &
                   relvar(:,k), mnuccc(:,k), nnuccc(:,k))

              ! make sure number of droplets frozen does not exceed available ice nuclei concentration
              ! this prevents 'runaway' droplet freezing


              where (qcic(:,k).ge.qsmall .and. t(:,k).lt.269.15_r8)
                 where (nnuccc(:,k)*lcldm(:,k).gt.nnuccd(:,k))
                    ! scale mixing ratio of droplet freezing with limit
                    mnuccc(:,k)=mnuccc(:,k)*(nnuccd(:,k)/(nnuccc(:,k)*lcldm(:,k)))
                    nnuccc(:,k)=nnuccd(:,k)/lcldm(:,k)
                 end where
              end where

              call contact_freezing(t(:,k), p(:,k), rndst(:,k,:), nacon(:,k,:), &
                   pgam(:,k), lamc(:,k), cdist1(:,k), qcic(:,k), &
                   relvar(:,k), mnucct(:,k), nnucct(:,k))

              mnudep(:,k)=0._r8
              nnudep(:,k)=0._r8

           else

              mi0l = qcic(:,k)/max(ncic(:,k), 1.0e6_r8/rho(:,k))
              mi0l = max(mi0l_min, mi0l)

              where (qcic(:,k) >= qsmall)
                 nnuccc(:,k) = frzimm(:,k)*1.0e6_r8/rho(:,k)
                 mnuccc(:,k) = nnuccc(:,k)*mi0l

                 nnucct(:,k) = frzcnt(:,k)*1.0e6_r8/rho(:,k)
                 mnucct(:,k) = nnucct(:,k)*mi0l

                 nnudep(:,k) = frzdep(:,k)*1.0e6_r8/rho(:,k)
                 mnudep(:,k) = nnudep(:,k)*mi0
              elsewhere
                 nnuccc(:,k) = 0._r8
                 mnuccc(:,k) = 0._r8

                 nnucct(:,k) = 0._r8
                 mnucct(:,k) = 0._r8

                 nnudep(:,k) = 0._r8
                 mnudep(:,k) = 0._r8
              end where

           end if
        else
              mnuccc(:,k)=0._r8
              nnuccc(:,k)=0._r8
              mnucct(:,k)=0._r8
              nnucct(:,k)=0._r8
              mnudep(:,k)=0._r8
              nnudep(:,k)=0._r8
        end if

        call snow_self_aggregation(t(:,k), rho(:,k), asn(:,k), qsic(:,k), nsic(:,k), &
             nsagg(:,k))

        call accrete_cloud_water_snow(t(:,k), rho(:,k), asn(:,k), uns(:,k), mu(:,k), &
             qcic(:,k), ncic(:,k), qsic(:,k), pgam(:,k), lamc(:,k), lams(:,k), n0s(:,k), &
             psacws(:,k), npsacws(:,k))

        if (do_cldice) then
           call secondary_ice_production(t(:,k), psacws(:,k), msacwi(:,k), nsacwi(:,k))
        else
           nsacwi(:,k) = 0.0_r8
           msacwi(:,k) = 0.0_r8
        end if

        call accrete_rain_snow(t(:,k), rho(:,k), umr(:,k), ums(:,k), unr(:,k), uns(:,k), &
             qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), &
             pracs(:,k), npracs(:,k))

        call heterogeneous_rain_freezing(t(:,k), qric(:,k), nric(:,k), lamr(:,k), &
             mnuccr(:,k), nnuccr(:,k))

        call accrete_cloud_water_rain(qric(:,k), qcic(:,k), ncic(:,k), &
             relvar(:,k), accre_enhan(:,k), pra(:,k), npra(:,k))

        call self_collection_rain(rho(:,k), qric(:,k), nric(:,k), nragg(:,k))

        if (do_cldice) then
           call accrete_cloud_ice_snow(t(:,k), rho(:,k), asn(:,k), qiic(:,k), niic(:,k), &
                qsic(:,k), lams(:,k), n0s(:,k), prai(:,k), nprai(:,k))
        else
           prai(:,k) = 0._r8
           nprai(:,k) = 0._r8
        end if

        call evaporate_sublimate_precip(deltat, t(:,k), p(:,k), rho(:,k), &
             dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k), &
             lcldm(:,k), cldmax(:,k), arn(:,k), asn(:,k), qcic(:,k), qiic(:,k), &
             qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), &
             pre(:,k), prds(:,k))

        call bergeron_process_snow(t(:,k), rho(:,k), dv(:,k), mu(:,k), sc(:,k), &
             qvl(:,k), qvi(:,k), asn(:,k), qcic(:,k), qsic(:,k), lams(:,k), n0s(:,k), &
             bergs(:,k))

         !+++PMC 12/3/12 - NEW VAPOR DEP/SUBLIMATION GOES HERE!!!
        if (do_cldice) then

           call ice_deposition_sublimation(deltat, t(:,k), q(:,k), qc(:,k), qi(:,k), ni(:,k), &
!!== KZ_DCS
                lcldm(:,k),icldm(:,k), naai(:,k), rho(:,k), dv(:,k), qvl(:,k), qvi(:,k), dcst(:,k), &
!!== KZ_DCS
                berg(:,k), vap_dep(:,k), ice_sublim(:,k))

           where (vap_dep(:,k) < 0._r8 .and. qi(:,k) > qsmall .and. icldm(:,k) > mincld)
              nsubi(:,k) = vap_dep(:,k) / qi(:,k) * ni(:,k) / icldm(:,k)
           elsewhere
              nsubi(:,k) = 0._r8
           end where

           ! bergeron process should not reduce nc unless
           ! all ql is removed (which is handled elsewhere)
           !in fact, nothing in this entire file makes nsubc nonzero.
           nsubc(:,k) = 0._r8

        end if !do_cldice
         !---PMC 12/3/12

        ! Big "administration" loop enforces conservation, updates variables
        ! that accumulate over substeps, and sets output variables.

        do i=1,mgncol

           ! conservation to ensure no negative values of cloud water/precipitation
           ! in case microphysical process rates are large
           !===================================================================

           ! note: for check on conservation, processes are multiplied by omsm
           ! to prevent problems due to round off error

           ! conservation of qc
           !-------------------------------------------------------------------

           dum = ((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+ &
                psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k))*deltat

           if (dum.gt.qc(i,k)) then
              ratio = qc(i,k)/deltat/((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+ &
               msacwi(i,k)+psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k))*omsm
              prc(i,k) = prc(i,k)*ratio
              pra(i,k) = pra(i,k)*ratio
              mnuccc(i,k) = mnuccc(i,k)*ratio
              mnucct(i,k) = mnucct(i,k)*ratio
              msacwi(i,k) = msacwi(i,k)*ratio
              psacws(i,k) = psacws(i,k)*ratio
              bergs(i,k) = bergs(i,k)*ratio
              berg(i,k) = berg(i,k)*ratio
              qcrat(i,k) = ratio
           else
              qcrat(i,k) = 1._r8
           end if

           !PMC 12/3/12: ratio is also frac of step w/ liquid.
           !thus we apply berg for "ratio" of timestep and vapor
           !deposition for the remaining frac of the timestep.
           if (qc(i,k) >= qsmall) then
              vap_dep(i,k) = vap_dep(i,k)*(1._r8-qcrat(i,k))
           end if

           !=================================================================
           ! apply limiter to ensure that ice/snow sublimation and rain evap
           ! don't push conditions into supersaturation, and ice deposition/nucleation don't
           ! push conditions into sub-saturation
           ! note this is done after qc conservation since we don't know how large
           ! vap_dep is before then
           ! estimates are only approximate since other process terms haven't been limited
           ! for conservation yet

           ! first limit ice deposition/nucleation vap_dep + mnuccd
           dum1 = vap_dep(i,k) + mnuccd(i,k)
           if (dum1 > 1.e-20_r8) then
              dum = (q(i,k)-qvi(i,k))/(1._r8 + cons28*qvi(i,k)/(cpp*rv*t(i,k)**2))/deltat
              dum = max(dum,0._r8)
              if (dum1 > dum) then
                 dum1=mnuccd(i,k)/(vap_dep(i,k)+mnuccd(i,k))
                 ! don't divide by cloud fraction since grid-mean rate
                 mnuccd(i,k)=dum*dum1/deltat

                 ! don't divide by cloud fraction since grid-mean rate
                 vap_dep(i,k)=dum*(1._r8-dum1)/deltat
              end if
           end if

           ! next limit ice and snow sublimation and rain evaporation
           ! get estimate of q and t at end of time step
           ! don't include other microphysical processes since they haven't
           ! been limited via conservation checks yet

           if ((pre(i,k)+prds(i,k))*cldmax(i,k)+ice_sublim(i,k) < -1.e-20_r8) then

              qtmp=q(i,k)-(ice_sublim(i,k)+vap_dep(i,k)+mnuccd(i,k)+ &
                   (pre(i,k)+prds(i,k))*cldmax(i,k))*deltat
              ttmp=t(i,k)+((pre(i,k)*cldmax(i,k))*xxlv+ &
                   (prds(i,k)*cldmax(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k))*xxls)*deltat/cpp

              ! If the unlimited tendencies are so large that ttmp is
              ! extremely low, qsat_water may hit a floating point
              ! exception, so just automatically limit temperatures below
              ! 50 K to prevent this.
              if (ttmp <= 50._r8) then
                 limit_precip_evap_sublim = .true.
              else
                 ! Else, check to see if we are pushing temperature down
                 ! and q up enough to become super-saturated.
                 call qsat_water(ttmp, p(i,k), esn, qvn)
                 limit_precip_evap_sublim = (qtmp > qvn)
              end if

              ! modify ice/precip evaporation rate if q > qsat
              if (limit_precip_evap_sublim) then

                 dum1=pre(i,k)*cldmax(i,k)/((pre(i,k)+prds(i,k))*cldmax(i,k)+ice_sublim(i,k))
                 dum2=prds(i,k)*cldmax(i,k)/((pre(i,k)+prds(i,k))*cldmax(i,k)+ice_sublim(i,k))
                 ! recalculate q and t after vap_dep and mnuccd but without evap or sublim
                 qtmp=q(i,k)-(vap_dep(i,k)+mnuccd(i,k))*deltat
                 ttmp=t(i,k)+((vap_dep(i,k)+mnuccd(i,k))*xxls)*deltat/cpp

                 ! use rhw to allow ice supersaturation
                 call qsat_water(ttmp, p(i,k), esn, qvn)

                 dum=(qtmp-qvn)/(1._r8 + cons27*qvn/(cpp*rv*ttmp**2))
                 dum=min(dum,0._r8)

                 ! modify rates if needed, divide by cldmax to get local (in-precip) value
                 pre(i,k)=dum*dum1/deltat/cldmax(i,k)

                 ! do separately using RHI for prds and ice_sublim
                 call qsat_ice(ttmp, p(i,k), esn, qvn)

                 dum=(qtmp-qvn)/(1._r8 + cons28*qvn/(cpp*rv*ttmp**2))
                 dum=min(dum,0._r8)

                 ! modify rates if needed, divide by cldmax to get local (in-precip) value
                 prds(i,k) = dum*dum2/deltat/cldmax(i,k)

                 ! don't divide ice_sublim by cloud fraction since it is grid-averaged
                 dum1 = (1._r8-dum1-dum2)
                 ice_sublim(i,k) = dum*dum1/deltat
              end if
           end if

           !===================================================================
           ! conservation of nc
           !-------------------------------------------------------------------
           dum = (nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+ &
                npsacws(i,k)-nsubc(i,k))*lcldm(i,k)*deltat

           if (dum.gt.nc(i,k)) then
              ratio = nc(i,k)/deltat/((nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+&
                   npsacws(i,k)-nsubc(i,k))*lcldm(i,k))*omsm

              nprc1(i,k) = nprc1(i,k)*ratio
              npra(i,k) = npra(i,k)*ratio
              nnuccc(i,k) = nnuccc(i,k)*ratio
              nnucct(i,k) = nnucct(i,k)*ratio
              npsacws(i,k) = npsacws(i,k)*ratio
              nsubc(i,k)=nsubc(i,k)*ratio
           end if

           if (do_cldice) then

              ! conservation of qi
              !-------------------------------------------------------------------
              dum = ((-mnuccc(i,k)-mnucct(i,k)-mnudep(i,k)-msacwi(i,k))*lcldm(i,k)+(prci(i,k)+ &
                   prai(i,k))*icldm(i,k)-ice_sublim(i,k)-vap_dep(i,k)-berg(i,k)-mnuccd(i,k))*deltat

              if (dum.gt.qi(i,k)) then
                 ratio = (qi(i,k)/deltat+vap_dep(i,k)+berg(i,k)+mnuccd(i,k)+ &
                      (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k))*lcldm(i,k))/ &
                      ((prci(i,k)+prai(i,k))*icldm(i,k)-ice_sublim(i,k))*omsm
                 prci(i,k) = prci(i,k)*ratio
                 prai(i,k) = prai(i,k)*ratio
                 ice_sublim(i,k) = ice_sublim(i,k)*ratio
              end if

              ! conservation of ni
              !-------------------------------------------------------------------
              if (use_hetfrz_classnuc) then
                 tmpfrz = nnuccc(i,k)
              else
                 tmpfrz = 0._r8
              end if
              dum = ((-nnucct(i,k)-tmpfrz-nnudep(i,k)-nsacwi(i,k))*lcldm(i,k)+(nprci(i,k)+ &
                   nprai(i,k)-nsubi(i,k))*icldm(i,k)-nnuccd(i,k))*deltat

              if (dum.gt.ni(i,k)) then
                 ratio = (ni(i,k)/deltat+nnuccd(i,k)+ &
                      (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k))*lcldm(i,k))/ &
                      ((nprci(i,k)+nprai(i,k)-nsubi(i,k))*icldm(i,k))*omsm
                 nprci(i,k) = nprci(i,k)*ratio
                 nprai(i,k) = nprai(i,k)*ratio
                 nsubi(i,k) = nsubi(i,k)*ratio
              end if

           end if

           ! for precipitation conservation, use logic that vertical integral
           ! of tendency from current level to top of model (i.e., qrtot) cannot be negative

           ! conservation of rain mixing rat
           !-------------------------------------------------------------------
           if (((prc(i,k)+pra(i,k))*lcldm(i,k)+(-mnuccr(i,k)+pre(i,k)-pracs(i,k))*&
                cldmax(i,k))*dz(i,k)*rho(i,k)+qrtot(i).lt.0._r8) then

              if (-pre(i,k)+pracs(i,k)+mnuccr(i,k).ge.qsmall) then

                 ratio = (qrtot(i)/(dz(i,k)*rho(i,k))+(prc(i,k)+pra(i,k))*lcldm(i,k))/&
                      ((-pre(i,k)+pracs(i,k)+mnuccr(i,k))*cldmax(i,k))*omsm

                 pre(i,k) = pre(i,k)*ratio
                 pracs(i,k) = pracs(i,k)*ratio
                 mnuccr(i,k) = mnuccr(i,k)*ratio
              end if
           end if

           ! conservation of nr
           !-------------------------------------------------------------------
           ! for now neglect evaporation of nr
           nsubr(i,k)=0._r8

           if ((nprc(i,k)*lcldm(i,k)+(-nnuccr(i,k)+nsubr(i,k)-npracs(i,k)&
                +nragg(i,k))*cldmax(i,k))*dz(i,k)*rho(i,k)+nrtot(i).lt.0._r8) then

              if (-nsubr(i,k)-nragg(i,k)+npracs(i,k)+nnuccr(i,k).ge.qsmall) then

                 ratio = (nrtot(i)/(dz(i,k)*rho(i,k))+nprc(i,k)*lcldm(i,k))/&
                      ((-nsubr(i,k)-nragg(i,k)+npracs(i,k)+nnuccr(i,k))*cldmax(i,k))*omsm

                 nsubr(i,k) = nsubr(i,k)*ratio
                 npracs(i,k) = npracs(i,k)*ratio
                 nnuccr(i,k) = nnuccr(i,k)*ratio
                 nragg(i,k) = nragg(i,k)*ratio
              end if
           end if

           ! conservation of snow mix ratio
           !-------------------------------------------------------------------
           if (((bergs(i,k)+psacws(i,k))*lcldm(i,k)+(prai(i,k)+prci(i,k))*icldm(i,k)+(pracs(i,k)+&
                mnuccr(i,k)+prds(i,k))*cldmax(i,k))*dz(i,k)*rho(i,k)+qstot(i).lt.0._r8) then

              if (-prds(i,k).ge.qsmall) then

                 ratio = (qstot(i)/(dz(i,k)*rho(i,k))+(bergs(i,k)+psacws(i,k))*lcldm(i,k)+(prai(i,k)+prci(i,k))*icldm(i,k)+&
                      (pracs(i,k)+mnuccr(i,k))*cldmax(i,k))/(-prds(i,k)*cldmax(i,k))*omsm

                 prds(i,k) = prds(i,k)*ratio
              else
                 prds(i,k) = 0._r8
              end if
           end if

           ! conservation of ns
           !-------------------------------------------------------------------
           ! calculate loss of number due to sublimation
           ! for now neglect sublimation of ns
           nsubs(i,k)=0._r8

           if ((nprci(i,k)*icldm(i,k)+(nnuccr(i,k)+nsubs(i,k)+nsagg(i,k))*cldmax(i,k))*&
                dz(i,k)*rho(i,k)+nstot(i).lt.0._r8) then

              if (-nsubs(i,k)-nsagg(i,k).ge.qsmall) then

                 ratio = (nstot(i)/(dz(i,k)*rho(i,k))+nprci(i,k)*icldm(i,k)+&
                      nnuccr(i,k)*cldmax(i,k))/((-nsubs(i,k)-nsagg(i,k))*cldmax(i,k))*omsm

                 nsubs(i,k) = nsubs(i,k)*ratio
                 nsagg(i,k) = nsagg(i,k)*ratio
              end if
           end if

           ! get tendencies due to microphysical conversion processes
           !==========================================================
           ! note: tendencies are multiplied by appropriate cloud/precip
           ! fraction to get grid-scale values
           ! note: vap_dep is already grid-average values

           qvlat(i,k) = qvlat(i,k)-(pre(i,k)+prds(i,k))*cldmax(i,k)-&
                vap_dep(i,k)-ice_sublim(i,k)-mnuccd(i,k)-mnudep(i,k)*lcldm(i,k)

           tlat(i,k) = tlat(i,k)+((pre(i,k)*cldmax(i,k)) &
                *xxlv+(prds(i,k)*cldmax(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k)+mnudep(i,k)*lcldm(i,k))*xxls+ &
                ((bergs(i,k)+psacws(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k))*lcldm(i,k)+(mnuccr(i,k)+ &
                pracs(i,k))*cldmax(i,k)+berg(i,k))*xlf)

           qctend(i,k) = qctend(i,k)+ &
                (-pra(i,k)-prc(i,k)-mnuccc(i,k)-mnucct(i,k)-msacwi(i,k)- &
                psacws(i,k)-bergs(i,k))*lcldm(i,k)-berg(i,k)

           if (do_cldice) then
              qitend(i,k) = qitend(i,k)+ &
                   (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k))*lcldm(i,k)+(-prci(i,k)- &
                   prai(i,k))*icldm(i,k)+vap_dep(i,k)+berg(i,k)+ice_sublim(i,k)+mnuccd(i,k)
           end if

           qrtend(i,k) = qrtend(i,k)+ &
                (pra(i,k)+prc(i,k))*lcldm(i,k)+(pre(i,k)-pracs(i,k)- &
                mnuccr(i,k))*cldmax(i,k)

           qstend(i,k) = qstend(i,k)+ &
                (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(prds(i,k)+ &
                pracs(i,k)+mnuccr(i,k))*cldmax(i,k)

           cmeout(i,k) = cmeout(i,k) + vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)

           ! add output for cmei (accumulate)
           cmeitot(i,k) = cmeitot(i,k) + vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)

           ! assign variables for trop_mozart, these are grid-average
           !-------------------------------------------------------------------
           ! evaporation/sublimation is stored here as positive term

           evapsnow(i,k) = evapsnow(i,k)-prds(i,k)*cldmax(i,k)
           nevapr(i,k) = nevapr(i,k)-pre(i,k)*cldmax(i,k)
           nevapr2(i,k) = nevapr2(i,k)-pre(i,k)*cldmax(i,k)

           ! change to make sure prain is positive: do not remove snow from
           ! prain used for wet deposition
           prain(i,k) = prain(i,k)+(pra(i,k)+prc(i,k))*lcldm(i,k)+(-pracs(i,k)- &
                mnuccr(i,k))*cldmax(i,k)
           prodsnow(i,k) = prodsnow(i,k)+(prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(&
                pracs(i,k)+mnuccr(i,k))*cldmax(i,k)

           ! following are used to calculate 1st order conversion rate of cloud water
           !    to rain and snow (1/s), for later use in aerosol wet removal routine
           ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
           !    used to calculate pra, prc, ... in this routine
           ! qcsinksum_rate1ord = { rate of direct transfer of cloud water to rain & snow }
           !                      (no cloud ice or bergeron terms)
           ! qcsum_rate1ord     = { qc used in calculation of the transfer terms }

           qcsinksum_rate1ord(i,k) = qcsinksum_rate1ord(i,k) + (pra(i,k)+prc(i,k)+psacws(i,k))*lcldm(i,k)
           qcsum_rate1ord(i,k) = qcsum_rate1ord(i,k) + qc(i,k)

           ! microphysics output, note this is grid-averaged
           pratot(i,k)=pratot(i,k)+pra(i,k)*lcldm(i,k)
           prctot(i,k)=prctot(i,k)+prc(i,k)*lcldm(i,k)
           mnuccctot(i,k)=mnuccctot(i,k)+mnuccc(i,k)*lcldm(i,k)
           mnuccttot(i,k)=mnuccttot(i,k)+mnucct(i,k)*lcldm(i,k)
           msacwitot(i,k)=msacwitot(i,k)+msacwi(i,k)*lcldm(i,k)
           psacwstot(i,k)=psacwstot(i,k)+psacws(i,k)*lcldm(i,k)
           bergstot(i,k)=bergstot(i,k)+bergs(i,k)*lcldm(i,k)

           bergtot(i,k)=bergtot(i,k)+berg(i,k)

           prcitot(i,k)=prcitot(i,k)+prci(i,k)*icldm(i,k)
           praitot(i,k)=praitot(i,k)+prai(i,k)*icldm(i,k)
           mnuccdtot(i,k)=mnuccdtot(i,k)+mnuccd(i,k)*icldm(i,k)

           pracstot(i,k)=pracstot(i,k)+pracs(i,k)*cldmax(i,k)
           mnuccrtot(i,k)=mnuccrtot(i,k)+mnuccr(i,k)*cldmax(i,k)

           nctend(i,k) = nctend(i,k)+&
                (-nnuccc(i,k)-nnucct(i,k)-npsacws(i,k)+nsubc(i,k) &
                -npra(i,k)-nprc1(i,k))*lcldm(i,k)

           if (do_cldice) then
              if (use_hetfrz_classnuc) then
                 tmpfrz = nnuccc(i,k)
              else
                 tmpfrz = 0._r8
              end if
              nitend(i,k) = nitend(i,k)+ nnuccd(i,k)+ &
                   (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k))*lcldm(i,k)+(nsubi(i,k)-nprci(i,k)- &
                   nprai(i,k))*icldm(i,k)
           end if

           nstend(i,k) = nstend(i,k)+(nsubs(i,k)+ &
                nsagg(i,k)+nnuccr(i,k))*cldmax(i,k)+nprci(i,k)*icldm(i,k)

           nrtend(i,k) = nrtend(i,k)+ &
                nprc(i,k)*lcldm(i,k)+(nsubr(i,k)-npracs(i,k)-nnuccr(i,k) &
                +nragg(i,k))*cldmax(i,k)

           ! make sure that ni at advanced time step does not exceed
           ! maximum (existing N + source terms*dt), which is possible if mtime < deltat
           ! note that currently mtime = deltat
           !================================================================

           if (do_cldice .and. nitend(i,k).gt.0._r8.and.ni(i,k)+nitend(i,k)*deltat.gt.nimax(i,k)) then
              nitend(i,k)=max(0._r8,(nimax(i,k)-ni(i,k))/deltat)
           end if

        end do

        ! End of "administration" loop

        ! get final values for precipitation q and N, based on
        ! flux of precip from above, source/sink term, and terminal fallspeed
        ! see eq. 15-16 in MG2008

        do i = 1, mgncol

           ! rain

           if (qric(i,k).ge.qsmall) then
              if (k .eq. 1) then
                 qric(i,k)=qrtend(i,k)*dz(i,k)/cldmax(i,k)/umr(i,k)
                 nric(i,k)=nrtend(i,k)*dz(i,k)/cldmax(i,k)/unr(i,k)
              else
                 qric(i,k) = (rho(i,k-1)*umr(i,k-1)*qric(i,k-1)*cldmax(i,k-1)+ &
                      (rho(i,k)*dz(i,k)*qrtend(i,k)))/(umr(i,k)*rho(i,k)*cldmax(i,k))
                 nric(i,k) = (rho(i,k-1)*unr(i,k-1)*nric(i,k-1)*cldmax(i,k-1)+ &
                      (rho(i,k)*dz(i,k)*nrtend(i,k)))/(unr(i,k)*rho(i,k)*cldmax(i,k))

              end if
           else
              qric(i,k)=0._r8
              nric(i,k)=0._r8
           end if

           ! snow

           if (qsic(i,k).ge.qsmall) then
              if (k .eq. 1) then
                 qsic(i,k)=qstend(i,k)*dz(i,k)/cldmax(i,k)/ums(i,k)
                 nsic(i,k)=nstend(i,k)*dz(i,k)/cldmax(i,k)/uns(i,k)
              else
                 qsic(i,k) = (rho(i,k-1)*ums(i,k-1)*qsic(i,k-1)*cldmax(i,k-1)+ &
                      (rho(i,k)*dz(i,k)*qstend(i,k)))/(ums(i,k)*rho(i,k)*cldmax(i,k))
                 nsic(i,k) = (rho(i,k-1)*uns(i,k-1)*nsic(i,k-1)*cldmax(i,k-1)+ &
                      (rho(i,k)*dz(i,k)*nstend(i,k)))/(uns(i,k)*rho(i,k)*cldmax(i,k))
              end if
           else
              qsic(i,k)=0._r8
              nsic(i,k)=0._r8
           end if

           ! calculate precipitation flux at surface
           !=========================================================
           ! divide by density of water to get units of m/s

           prect(i) = prect(i)+(qrtend(i,k)*dz(i,k)*rho(i,k)+&
                qstend(i,k)*dz(i,k)*rho(i,k))/rhow
           preci(i) = preci(i)+qstend(i,k)*dz(i,k)*rho(i,k)/rhow

           ! convert rain rate from m/s to mm/hr

           rainrt(i,k)=rainrt(i,k) + (qric(i,k)*rho(i,k)*umr(i,k)/rhow*3600._r8*1000._r8)

           ! vertically-integrated precip source/sink terms (note: grid-averaged)

           qrtot(i) = max(qrtot(i)+qrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
           qstot(i) = max(qstot(i)+qstend(i,k)*dz(i,k)*rho(i,k),0._r8)
           nrtot(i) = max(nrtot(i)+nrtend(i,k)*dz(i,k)*rho(i,k),0._r8)
           nstot(i) = max(nstot(i)+nstend(i,k)*dz(i,k)*rho(i,k),0._r8)

           ! calculate melting and freezing of precip
           !=========================================================

           ! melt snow at +2 C

           ! Difference in amount of heat between temperature at which
           ! all snow melts, and the current state.
           dum = (snowmelt - t(i,k) - tlat(i,k)/cpp*deltat) * cpp

           ! Test if temperature is above threshold and snow is present.
           if (dum < 0._r8 .and. qstot(i) > 0._r8) then

              ! (negative) heat produced if all snow is melted
              dum1 = -xlf * qstot(i)/(dz(i,k)*rho(i,k))*deltat

              ! ratio of heating needed to get to the threshold, to
              ! the total heating from melting everything, is equal
              ! to the proportion of snow that actually melts
              dum = min(1._r8, dum/dum1)
              dum = max(0._r8, dum)

              ! Melt snow
              qric(i,k)=qric(i,k)+dum*qsic(i,k)
              nric(i,k)=nric(i,k)+dum*nsic(i,k)
              qsic(i,k)=(1._r8-dum)*qsic(i,k)
              nsic(i,k)=(1._r8-dum)*nsic(i,k)

              qrtot(i)=qrtot(i)+dum*qstot(i)
              nrtot(i)=nrtot(i)+dum*nstot(i)
              qstot(i)=(1._r8-dum)*qstot(i)
              nstot(i)=(1._r8-dum)*nstot(i)

              preci(i)=(1._r8-dum)*preci(i)

              ! Get heating tendency based on proportion of snow that
              ! actually melts.
              dum1 = dum * dum1/deltat

              meltsdttot(i,k)=meltsdttot(i,k) + dum1
              tlat(i,k)=tlat(i,k)+dum1

           end if

           ! freeze all rain at -5C for Arctic
           !=========================================================

           ! Difference in amount of heat between temperature at which
           ! all rain freezes, and the current state.
           dum = (rainfrze - t(i,k) - tlat(i,k)/cpp*deltat) * cpp

           ! Test if temperature is below threshold and snow is present.
           if (dum > 0._r8 .and. qrtot(i) > 0._r8) then

              ! heat produced if all rain freezes
              dum1 = xlf * qrtot(i)/(dz(i,k)*rho(i,k))*deltat

              ! ratio of heating needed to get to the threshold, to
              ! the total heating from freezing everything, is equal
              ! to the proportion of rain that actually freezes
              dum = min(1._r8, dum/dum1)
              dum = max(0._r8, dum)

              ! Freeze rain
              qsic(i,k)=qsic(i,k)+dum*qric(i,k)
              nsic(i,k)=nsic(i,k)+dum*nric(i,k)
              qric(i,k)=(1._r8-dum)*qric(i,k)
              nric(i,k)=(1._r8-dum)*nric(i,k)

              qstot(i)=qstot(i)+dum*qrtot(i)
              nstot(i)=nstot(i)+dum*nrtot(i)
              qrtot(i)=(1._r8-dum)*qrtot(i)
              nrtot(i)=(1._r8-dum)*nrtot(i)

              preci(i)=preci(i)+dum*(prect(i)-preci(i))

              ! Get heating tendency based on proportion of rain that
              ! actually freezes.
              dum1 = dum * dum1/deltat

              frzrdttot(i,k)=frzrdttot(i,k) + dum1
              tlat(i,k)=tlat(i,k)+dum1

           end if

           ! if rain/snow mix ratio is zero so should number concentration
           !=========================================================

           if (qsic(i,k) < qsmall) then
              qsic(i,k)=0._r8
              nsic(i,k)=0._r8
           end if

           if (qric(i,k) < qsmall) then
              qric(i,k)=0._r8
              nric(i,k)=0._r8
           end if

           ! make sure number concentration is a positive number to avoid
           ! taking root of negative

           nric(i,k)=max(nric(i,k),0._r8)
           nsic(i,k)=max(nsic(i,k),0._r8)

           ! get size distribution parameters for fallspeed calculations
           !=========================================================

           ! rain

           call size_dist_param_rain(qric(i,k), nric(i,k), lamr(i,k), n0r(i,k))

           if (lamr(i,k).ge.qsmall) then

              ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

              unr(i,k) = min(arn(i,k)*cons4/lamr(i,k)**br,9.1_r8*rhof(i,k))
              umr(i,k) = min(arn(i,k)*cons5/(6._r8*lamr(i,k)**br),9.1_r8*rhof(i,k))

           else
              umr(i,k)=0._r8
              unr(i,k)=0._r8
           end if

           !......................................................................
           ! snow

           call size_dist_param_snow(qsic(i,k), nsic(i,k), lams(i,k), n0s(i,k))

           if (lams(i,k) > 0._r8) then

              ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)

              ums(i,k) = min(asn(i,k)*cons8/(6._r8*lams(i,k)**bs),1.2_r8*rhof(i,k))
              uns(i,k) = min(asn(i,k)*cons7/lams(i,k)**bs,1.2_r8*rhof(i,k))

           else
              ums(i,k) = 0._r8
              uns(i,k) = 0._r8
           end if

        end do

        ! Done with vertical dependencies from precipitation.

     end do micro_vert_loop ! end k loop

     !-----------------------------------------------------
     ! convert rain/snow q and N for output to history, note,
     ! output is for gridbox average

     ! calculate precip fluxes
     ! calculate the precip flux (kg/m2/s) as mixingratio(kg/kg)*airdensity(kg/m3)*massweightedfallspeed(m/s)
     ! ---------------------------------------------------------------------

     rflx(:,2:) = rflx(:,2:) + (qric*rho*umr*cldmax)

     dumr(:,1)  = int_to_mid(:,1) * qric(:,1)
     qrout(:,1)  = qrout(:,1)  + int_to_mid(:,1)*qric(:,1)*cldmax(:,1)

     dumr(:,2:)  = interp_to_mid(qric,int_to_mid(:,2:))
     qrout(:,2:) = qrout(:,2:) + interp_to_mid(qric * cldmax, int_to_mid(:,2:))

     dumnr(:,1) = int_to_mid(:,1) * nric(:,1)
     nrout(:,1)  = nrout(:,1)  + &
          (int_to_mid(:,1)*nric(:,1)*cldmax(:,1)*rho(:,1))

     dumnr(:,2:) = interp_to_mid(nric,int_to_mid(:,2:))
     nrout(:,2:) = nrout(:,2:) + interp_to_mid(nric * cldmax * rho, int_to_mid(:,2:))
     ! Calculate rercld

     ! calculate mean size of combined rain and cloud water

     ! hm 11-22-11 modify to interpolate rain from interface to mid-point
     ! logic is to interpolate rain mass and number, then recalculate PSD
     ! parameters to get relevant parameters for mean size

     ! interpolate rain mass and number, store in dummy variables

     ! calculate n0r and lamr from interpolated mid-point rain mass and number
     ! divide by precip fraction to get in-precip (local) values of
     ! rain mass and number, divide by rhow to get rain number in kg^-1

     call size_dist_param_rain(dumr, dumnr, lamr, n0r)

     call calc_rercld(lamr, n0r, lamc, cdist1, pgam, dumr, qcic, &
          rercld)

     nsout(:,1)  = nsout(:,1)  + &
          (int_to_mid(:,1)*nsic(:,1)*cldmax(:,1)*rho(:,1))
     nsout(:,2:) = nsout(:,2:) + interp_to_mid(nsic * cldmax * rho, int_to_mid(:,2:))

     qsout(:,1)  = qsout(:,1)  + int_to_mid(:,1)*qsic(:,1)*cldmax(:,1)
     qsout(:,2:) = qsout(:,2:) + interp_to_mid(qsic * cldmax, int_to_mid(:,2:))

     sflx(:,2:) = sflx(:,2:) + (qsic*rho*ums*cldmax)

  ! assign variables back to start-of-timestep values
  !hm note: only nc is modified above (droplet activation tendency is added on)
  !hm       thus only nc needs to be assigned to start-of-timestep values
  !================================================================================

  nc = ncn

  !.............................................................................

  !================================================================================

  ! Re-apply droplet activation tendency
  nctend = nctend + npccn

  ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
  nevapr = nevapr + evapsnow
  prer_evap = nevapr2
  prain = prain + prodsnow

  sed_col_loop: do i=1,mgncol

     do k=1,nlev

        ! calculate sedimentation for cloud water and ice
        !================================================================================

        ! update in-cloud cloud mixing ratio and number concentration
        ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
        ! note: these are in-cloud values***, hence we divide by cloud fraction

        dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)/lcldm(i,k)
        dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)/icldm(i,k)
        dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k),0._r8)
        dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat)/icldm(i,k),0._r8)

        ! hm add 6/2/11 switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)
        end if

        ! hm add 6/2/11 switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)
        end if

        ! obtain new slope parameter to avoid possible singularity
!!== KZ_DCS 
        call size_dist_param_ice(dumi(i,k), dcst(i,k), dumni(i,k), lami(i,k), n0i(i,k)) 
!!== KZ_DCS 

        call size_dist_param_liq(dumc(i,k), dumnc(i,k), cdnl, rho(i,k), .true., &
             pgam(i,k), lamc(i,k))

        ! calculate number and mass weighted fall velocity for droplets and cloud ice
        !-------------------------------------------------------------------


        if (dumc(i,k).ge.qsmall) then

           vtrmc(i,k)=acn(i,k)*gamma(4._r8+bc+pgam(i,k))/ &
                (lamc(i,k)**bc*gamma(pgam(i,k)+4._r8))

           fc(k) = g*rho(i,k)*vtrmc(i,k)

           fnc(k) = g*rho(i,k)* &
                acn(i,k)*gamma(1._r8+bc+pgam(i,k))/ &
                (lamc(i,k)**bc*gamma(pgam(i,k)+1._r8))
        else
           fc(k) = 0._r8
           fnc(k)= 0._r8
        end if

        ! calculate number and mass weighted fall velocity for cloud ice

        if (dumi(i,k).ge.qsmall) then

           vtrmi(i,k)=min(ain(i,k)*cons17/(6._r8*lami(i,k)**bi), &
                1.2_r8*rhof(i,k))

           fi(k) = g*rho(i,k)*vtrmi(i,k)
           fni(k) = g*rho(i,k)* &
                min(ain(i,k)*cons16/lami(i,k)**bi,1.2_r8*rhof(i,k))
        else
           fi(k) = 0._r8
           fni(k)= 0._r8
        end if

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation

        dumc(i,k) = (qc(i,k)+qctend(i,k)*deltat)
        dumi(i,k) = (qi(i,k)+qitend(i,k)*deltat)
        dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat),0._r8)
        dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat),0._r8)

        if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
        if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8

     end do       !!! vertical loop

     ! initialize nstep for sedimentation sub-steps

     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + int(max( &
          maxval( fi/pdel(i,:)), &
          maxval( fc/pdel(i,:)), &
          maxval(fni/pdel(i,:)), &
          maxval(fnc/pdel(i,:))) &
          * deltat)

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     do n = 1,nstep

        if (do_cldice) then
           falouti  = fi  * dumi(i,:)
           faloutni = fni * dumni(i,:)
        else
           falouti  = 0._r8
           faloutni = 0._r8
        end if
        faloutc  = fc  * dumc(i,:)
        faloutnc = fnc * dumnc(i,:)

        ! top of model

        k = 1
        faltndi = falouti(k)/pdel(i,k)
        faltndni = faloutni(k)/pdel(i,k)
        faltndc = faloutc(k)/pdel(i,k)
        faltndnc = faloutnc(k)/pdel(i,k)

        ! add fallout terms to microphysical tendencies

        qitend(i,k) = qitend(i,k)-faltndi/nstep
        nitend(i,k) = nitend(i,k)-faltndni/nstep
        qctend(i,k) = qctend(i,k)-faltndc/nstep
        nctend(i,k) = nctend(i,k)-faltndnc/nstep

        ! sedimentation tendencies for output
        qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
        qisedten(i,k)=qisedten(i,k)-faltndi/nstep

        dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
        dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
        dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
        dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

        do k = 2,nlev

           ! for cloud liquid and ice, if cloud fraction increases with height
           ! then add flux from above to both vapor and cloud water of current level
           ! this means that flux entering clear portion of cell from above evaporates
           ! instantly

           dum=lcldm(i,k)/lcldm(i,k-1)
           dum=min(dum,1._r8)
           dum1=icldm(i,k)/icldm(i,k-1)
           dum1=min(dum1,1._r8)

           faltndqie=(falouti(k)-falouti(k-1))/pdel(i,k)
           faltndi=(falouti(k)-dum1*falouti(k-1))/pdel(i,k)
           faltndni=(faloutni(k)-dum1*faloutni(k-1))/pdel(i,k)
           faltndqce=(faloutc(k)-faloutc(k-1))/pdel(i,k)
           faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel(i,k)
           faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel(i,k)

           ! add fallout terms to eulerian tendencies

           qitend(i,k) = qitend(i,k)-faltndi/nstep
           nitend(i,k) = nitend(i,k)-faltndni/nstep
           qctend(i,k) = qctend(i,k)-faltndc/nstep
           nctend(i,k) = nctend(i,k)-faltndnc/nstep

           ! sedimentation tendencies for output
           qcsedten(i,k)=qcsedten(i,k)-faltndc/nstep
           qisedten(i,k)=qisedten(i,k)-faltndi/nstep

           ! add terms to to evap/sub of cloud water

           qvlat(i,k)=qvlat(i,k)-(faltndqie-faltndi)/nstep
           ! for output
           qisevap(i,k)=qisevap(i,k)-(faltndqie-faltndi)/nstep
           qvlat(i,k)=qvlat(i,k)-(faltndqce-faltndc)/nstep
           ! for output
           qcsevap(i,k)=qcsevap(i,k)-(faltndqce-faltndc)/nstep

           tlat(i,k)=tlat(i,k)+(faltndqie-faltndi)*xxls/nstep
           tlat(i,k)=tlat(i,k)+(faltndqce-faltndc)*xxlv/nstep

           dumi(i,k) = dumi(i,k)-faltndi*deltat/nstep
           dumni(i,k) = dumni(i,k)-faltndni*deltat/nstep
           dumc(i,k) = dumc(i,k)-faltndc*deltat/nstep
           dumnc(i,k) = dumnc(i,k)-faltndnc*deltat/nstep

        end do   !! k loop

        ! units below are m/s
        ! cloud water/ice sedimentation flux at surface
        ! is added to precip flux at surface to get total precip (cloud + precip water)
        ! rate

        prect(i) = prect(i)+(faloutc(nlev)+falouti(nlev))/g/nstep/1000._r8
        preci(i) = preci(i)+(              falouti(nlev))/g/nstep/1000._r8

     end do   !! nstep loop

     ! end sedimentation
     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     ! get new update for variables that includes sedimentation tendency
     ! note : here dum variables are grid-average, NOT in-cloud

     do k=1,nlev

        dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)
        dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)
        dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)
        dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)

        ! hm add 6/2/11 switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)*lcldm(i,k)
        end if

        ! hm add 6/2/11 switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)*icldm(i,k)
        end if

        if (dumc(i,k).lt.qsmall) dumnc(i,k)=0._r8
        if (dumi(i,k).lt.qsmall) dumni(i,k)=0._r8

        ! calculate instantaneous processes (melting, homogeneous freezing)
        !====================================================================

        if (do_cldice) then
           if (t(i,k)+tlat(i,k)/cpp*deltat > tmelt) then
              if (dumi(i,k) > 0._r8) then

                 ! limit so that melting does not push temperature below freezing
                 !-----------------------------------------------------------------
                 dum = -dumi(i,k)*xlf/cpp
                 if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt.tmelt) then
                    dum = (t(i,k)+tlat(i,k)/cpp*deltat-tmelt)*cpp/xlf
                    dum = dum/dumi(i,k)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qctend(i,k)=qctend(i,k)+dum*dumi(i,k)/deltat

                 ! for output
                 melttot(i,k)=dum*dumi(i,k)/deltat

                 ! assume melting ice produces droplet
                 ! mean volume radius of 8 micron

                 nctend(i,k)=nctend(i,k)+3._r8*dum*dumi(i,k)/deltat/ &
                      (4._r8*pi*5.12e-16_r8*rhow)

                 qitend(i,k)=((1._r8-dum)*dumi(i,k)-qi(i,k))/deltat
                 nitend(i,k)=((1._r8-dum)*dumni(i,k)-ni(i,k))/deltat
                 tlat(i,k)=tlat(i,k)-xlf*dum*dumi(i,k)/deltat
              end if
           end if

           ! homogeneously freeze droplets at -40 C
           !-----------------------------------------------------------------

           if (t(i,k)+tlat(i,k)/cpp*deltat < 233.15_r8) then
              if (dumc(i,k) > 0._r8) then

                 ! limit so that freezing does not push temperature above threshold
                 dum = dumc(i,k)*xlf/cpp
                 if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.233.15_r8) then
                    dum = -(t(i,k)+tlat(i,k)/cpp*deltat-233.15_r8)*cpp/xlf
                    dum = dum/dumc(i,k)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qitend(i,k)=qitend(i,k)+dum*dumc(i,k)/deltat
                 ! for output
                 homotot(i,k)=dum*dumc(i,k)/deltat

                 ! assume 25 micron mean volume radius of homogeneously frozen droplets
                 ! consistent with size of detrained ice in stratiform.F90
                 nitend(i,k)=nitend(i,k)+dum*3._r8*dumc(i,k)/(4._r8*3.14_r8*1.563e-14_r8* &
                      500._r8)/deltat
                 qctend(i,k)=((1._r8-dum)*dumc(i,k)-qc(i,k))/deltat
                 nctend(i,k)=((1._r8-dum)*dumnc(i,k)-nc(i,k))/deltat
                 tlat(i,k)=tlat(i,k)+xlf*dum*dumc(i,k)/deltat
              end if
           end if

           ! remove any excess over-saturation, which is possible due to non-linearity when adding
           ! together all microphysical processes
           !-----------------------------------------------------------------
           ! follow code similar to old CAM scheme

           qtmp=q(i,k)+qvlat(i,k)*deltat
           ttmp=t(i,k)+tlat(i,k)/cpp*deltat

           ! use rhw to allow ice supersaturation
           call qsat_water(ttmp, p(i,k), esn, qvn)
           qvn = min(qvn,1._r8)

           if (qtmp > qvn .and. qvn > 0) then
              ! expression below is approximate since there may be ice deposition
              dum = (qtmp-qvn)/(1._r8+cons27*qvn/(cpp*rv*ttmp**2))/deltat
              ! add to output cme
              cmeout(i,k) = cmeout(i,k)+dum
              ! now add to tendencies, partition between liquid and ice based on temperature
              if (ttmp > 268.15_r8) then
                 dum1=0.0_r8
                 ! now add to tendencies, partition between liquid and ice based on te
                 !-------------------------------------------------------
              else if (ttmp < 238.15_r8) then
                 dum1=1.0_r8
              else
                 dum1=(268.15_r8-ttmp)/30._r8
              end if

              dum = (qtmp-qvn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                   *qvn/(cpp*rv*ttmp**2))/deltat
              qctend(i,k)=qctend(i,k)+dum*(1._r8-dum1)
              ! for output
              qcrestot(i,k)=dum*(1._r8-dum1)
              qitend(i,k)=qitend(i,k)+dum*dum1
              qirestot(i,k)=dum*dum1
              qvlat(i,k)=qvlat(i,k)-dum
              ! for output
              qvres(i,k)=-dum
              tlat(i,k)=tlat(i,k)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
           end if
        end if

        ! calculate effective radius for pass to radiation code
        !=========================================================
        ! if no cloud water, default value is 10 micron for droplets,
        ! 25 micron for cloud ice

        ! update cloud variables after instantaneous processes to get effective radius
        ! variables are in-cloud to calculate size dist parameters

        dumc(i,k) = max(qc(i,k)+qctend(i,k)*deltat,0._r8)/lcldm(i,k)
        dumi(i,k) = max(qi(i,k)+qitend(i,k)*deltat,0._r8)/icldm(i,k)
        dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat,0._r8)/lcldm(i,k)
        dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat,0._r8)/icldm(i,k)

        ! hm add 6/2/11 switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k)=ncnst/rho(i,k)
        end if

        ! hm add 6/2/11 switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k)=ninst/rho(i,k)
        end if

        ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1

        dumc(i,k)=min(dumc(i,k),5.e-3_r8)
        dumi(i,k)=min(dumi(i,k),5.e-3_r8)


        ! cloud ice effective radius
        !-----------------------------------------------------------------

        if (do_cldice) then
           if (dumi(i,k).ge.qsmall) then

              dum_2D(i,k) = dumni(i,k)
!!== KZ_DCS 
              call size_dist_param_ice(dumi(i,k), dcst(i,k), dumni(i,k), lami(i,k), n0i(i,k))
!!== KZ_DCS 

              if (dumni(i,k) /=dum_2D(i,k)) then
                 ! adjust number conc if needed to keep mean size in reasonable range
                 nitend(i,k)=(dumni(i,k)*icldm(i,k)-ni(i,k))/deltat
              end if

              effi(i,k) = 1.5_r8/lami(i,k)*1.e6_r8

           else
              effi(i,k) = 25._r8
           end if

           ! ice effective diameter for david mitchell's optics
           deffi(i,k)=effi(i,k)*rhoi/917._r8*2._r8
        else
           ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
           ! radius has already been determined from the size distribution.
           effi(i,k) = re_ice(i,k) * 1.e6_r8      ! m -> um
           deffi(i,k)=effi(i,k) * 2._r8
        end if

        ! cloud droplet effective radius
        !-----------------------------------------------------------------
        if (dumc(i,k).ge.qsmall) then


           ! hm add 6/2/11 switch for specification of droplet and crystal number
           if (nccons) then
              ! make sure nc is consistence with the constant N by adjusting tendency, need
              ! to multiply by cloud fraction
              ! note that nctend may be further adjusted below if mean droplet size is
              ! out of bounds

              nctend(i,k)=(ncnst/rho(i,k)*lcldm(i,k)-nc(i,k))/deltat

           end if

           dum = dumnc(i,k)

           call size_dist_param_liq(dumc(i,k), dumnc(i,k), cdnl, rho(i,k), .true., &
                pgam(i,k), lamc(i,k))

           if (dum /= dumnc(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nctend(i,k)=(dumnc(i,k)*lcldm(i,k)-nc(i,k))/deltat
           end if

           effc(i,k) = gamma(pgam(i,k)+4._r8)/ &
                gamma(pgam(i,k)+3._r8)/lamc(i,k)/2._r8*1.e6_r8
           !assign output fields for shape here
           lamcrad(i,k)=lamc(i,k)
           pgamrad(i,k)=pgam(i,k)


           ! recalculate effective radius for constant number, in order to separate
           ! first and second indirect effects
           !======================================
           ! assume constant number of 10^8 kg-1

           dumnc(i,k)=1.e8_r8

           ! Pass in "false" adjust flag to prevent number from being changed within
           ! size distribution subroutine.
           call size_dist_param_liq(dumc(i,k), dumnc(i,k), cdnl, rho(i,k), .false., &
                pgam(i,k), lamc(i,k))

           effc_fn(i,k) = gamma(pgam(i,k)+4._r8)/ &
                gamma(pgam(i,k)+3._r8)/lamc(i,k)/2._r8*1.e6_r8

        else
           effc(i,k) = 10._r8
           lamcrad(i,k)=0._r8
           pgamrad(i,k)=0._r8
           effc_fn(i,k) = 10._r8
        end if

     end do ! vertical k loop

     do k=1,nlev
        ! if updated q (after microphysics) is zero, then ensure updated n is also zero
        !=================================================================================
        if (qc(i,k)+qctend(i,k)*deltat.lt.qsmall) nctend(i,k)=-nc(i,k)/deltat
        if (do_cldice .and. qi(i,k)+qitend(i,k)*deltat.lt.qsmall) nitend(i,k)=-ni(i,k)/deltat

     end do

  end do sed_col_loop! i loop

  ! DO STUFF FOR OUTPUT:
  !==================================================

  ! qc and qi are only used for output calculations past here,
  ! so add qctend and qitend back in one more time
  qc = qc + qctend*deltat
  qi = qi + qitend*deltat

  ! averaging for snow and rain number and diameter
  !--------------------------------------------------

  ! drout2/dsout2:
  ! diameter of rain and snow
  ! dsout:
  ! scaled diameter of snow (passed to radiation in CAM)
  ! reff_rain/reff_snow:
  ! calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual

  where (qrout .gt. 1.e-7_r8 &
       .and. nrout.gt.0._r8)
     qrout2 = qrout * cldmax
     nrout2 = nrout * cldmax
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just drout2 times constants.
     drout2 = avg_diameter(qrout, nrout, rho, rhow)
     freqr = cldmax

     reff_rain=1.5_r8*drout2*1.e6_r8
  elsewhere
     qrout2 = 0._r8
     nrout2 = 0._r8
     drout2 = 0._r8
     freqr = 0._r8
     reff_rain = 0._r8
  end where

  where (qsout .gt. 1.e-7_r8 &
       .and. nsout.gt.0._r8)
     qsout2 = qsout * cldmax
     nsout2 = nsout * cldmax
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just dsout2 times constants.
     dsout2 = avg_diameter(qsout, nsout, rho, rhosn)
     freqs = cldmax

     dsout=3._r8*rhosn/917._r8*dsout2

     reff_snow=1.5_r8*dsout2*1.e6_r8
  elsewhere
     dsout  = 0._r8
     qsout2 = 0._r8
     nsout2 = 0._r8
     dsout2 = 0._r8
     freqs  = 0._r8
     reff_snow=0._r8
  end where

  ! analytic radar reflectivity
  !--------------------------------------------------
  ! formulas from Matthew Shupe, NOAA/CERES
  ! *****note: radar reflectivity is local (in-precip average)
  ! units of mm^6/m^3

  do i = 1,mgncol
     do k=1,nlev
        if (qc(i,k).ge.qsmall  .and. (nc(i,k)+nctend(i,k)*deltat).gt.10._r8) then
           dum=(qc(i,k)/lcldm(i,k)*rho(i,k)*1000._r8)**2 &
                /(0.109_r8*(nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)/1.e6_r8)*lcldm(i,k)/cldmax(i,k)
        else
           dum=0._r8
        end if
        if (qi(i,k).ge.qsmall) then
           dum1=(qi(i,k)*rho(i,k)/icldm(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*icldm(i,k)/cldmax(i,k)
        else
           dum1=0._r8
        end if

        if (qsout(i,k).ge.qsmall) then
           dum1=dum1+(qsout(i,k)*rho(i,k)*1000._r8/0.1_r8)**(1._r8/0.63_r8)
        end if

        refl(i,k)=dum+dum1

        ! add rain rate, but for 37 GHz formulation instead of 94 GHz
        ! formula approximated from data of Matrasov (2007)
        ! rainrt is the rain rate in mm/hr
        ! reflectivity (dum) is in DBz

        if (rainrt(i,k).ge.0.001_r8) then
           dum=log10(rainrt(i,k)**6._r8)+16._r8

           ! convert from DBz to mm^6/m^3

           dum = 10._r8**(dum/10._r8)
        else
           ! don't include rain rate in R calculation for values less than 0.001 mm/hr
           dum=0._r8
        end if

        ! add to refl

        refl(i,k)=refl(i,k)+dum

        !output reflectivity in Z.
        areflz(i,k)=refl(i,k) * cldmax(i,k)

        ! convert back to DBz

        if (refl(i,k).gt.minrefl) then
           refl(i,k)=10._r8*log10(refl(i,k))
        else
           refl(i,k)=-9999._r8
        end if

        !set averaging flag
        if (refl(i,k).gt.mindbz) then
           arefl(i,k)=refl(i,k) * cldmax(i,k)
           frefl(i,k)=cldmax(i,k)
        else
           arefl(i,k)=0._r8
           areflz(i,k)=0._r8
           frefl(i,k)=0._r8
        end if

        ! bound cloudsat reflectivity

        csrfl(i,k)=min(csmax,refl(i,k))

        !set averaging flag
        if (csrfl(i,k).gt.csmin) then
           acsrfl(i,k)=refl(i,k) * cldmax(i,k)
           fcsrfl(i,k)=cldmax(i,k)
        else
           acsrfl(i,k)=0._r8
           fcsrfl(i,k)=0._r8
        end if

     end do
  end do

  !redefine fice here....
  dum_2D = qsout + qrout + qc + qi
  dumi = qsout + qi
  where (dumi .gt. qsmall .and. dum_2D .gt. qsmall)
     nfice=min(dumi/dum_2D,1._r8)
  elsewhere
     nfice=0._r8
  end where

  ! Avoid zero/near-zero division.
  qcsinksum_rate1ord = qcsinksum_rate1ord/max(qcsum_rate1ord,1.0e-30_r8)

end subroutine micro_mg_tend

!========================================================================
!FORMULAS
!========================================================================

! Calculate correction due to latent heat for evaporation/sublimation
elemental function calc_ab(t, qv, xxl) result(ab)
  real(r8), intent(in) :: t     ! Temperature
  real(r8), intent(in) :: qv    ! Saturation vapor pressure
  real(r8), intent(in) :: xxl   ! Latent heat

  real(r8) :: ab

  real(r8) :: dqsdt

  dqsdt = xxl*qv / (rv * t**2)
  ab = 1._r8 + dqsdt*xxl/cpp

end function calc_ab

! get cloud droplet size distribution parameters
elemental subroutine size_dist_param_liq(qcic, ncic, cdnl, rho, nadjflag, pgam, lamc)

  real(r8), intent(in) :: qcic
  real(r8), intent(inout) :: ncic
  real(r8), intent(in) :: cdnl
  real(r8), intent(in) :: rho
  logical,  intent(in) :: nadjflag ! Whether to adjust number concentration to fall
                                   ! within certain bounds

  real(r8), intent(out) :: pgam
  real(r8), intent(out) :: lamc

  real(r8) :: dumgam1
  real(r8) :: dumgam2
  real(r8) :: lammin
  real(r8) :: lammax

  if (qcic > qsmall) then

     if (nadjflag) then
        ! add upper limit to in-cloud number concentration to prevent numerical error
        ncic=min(ncic,qcic*1.e20_r8)
        ! add lower limit to in-cloud number concentration
        ncic=max(ncic,cdnl/rho) ! sghan minimum in #/cm
     end if

     ! get pgam from fit to observations of martin et al. 1994

     pgam=0.0005714_r8*(ncic/1.e6_r8*rho)+0.2714_r8
     pgam=1._r8/(pgam**2)-1._r8
     pgam=max(pgam,2._r8)
     pgam=min(pgam,15._r8)

     ! calculate lamc
     dumgam1 = gamma(pgam+1._r8)
     dumgam2 = gamma(pgam+4._r8)

     lamc = (pi/6._r8*rhow*ncic*dumgam2/ &
          (qcic*dumgam1))**(1._r8/3._r8)

     ! lammin, 50 micron diameter max mean size
     ! omsm fudge factors are to guaranteed that lamcrad falls within the
     ! table used by RRTMG.
     lammin = (pgam+1._r8)/(omsm*50.e-6_r8)
     lammax = (pgam+1._r8)*omsm/2.e-6_r8

     if (lamc < lammin) then
        lamc = lammin

        if (nadjflag) then
           ncic = 6._r8 * lamc**3 * qcic * dumgam1/ &
                (pi * rhow * dumgam2)
        end if
     else if (lamc > lammax) then
        lamc = lammax

        if (nadjflag) then
           ncic = 6._r8 * lamc**3 * qcic * dumgam1/ &
                (pi * rhow * dumgam2)
        end if
     end if

  else
     ! pgam not calculated in this case, so set it to a value likely to cause an error
     ! if it is accidentally used
     ! (gamma function undefined for negative integers)
     pgam = -100._r8
     lamc = 0._r8
  end if

end subroutine size_dist_param_liq

! get ice size distribution parameters
!!== KZ_DCS 
elemental subroutine size_dist_param_ice(qiic, dcst, niic, lami, n0i)
!!== KZ_DCS 
  real(r8), intent(in) :: qiic
!!== KZ_DCS 
  real(r8), intent(in) :: dcst 
!!== KZ_DCS 
  real(r8), intent(inout) :: niic

  real(r8), intent(out) :: lami
  real(r8), intent(out) :: n0i

  ! particle mass-diameter relationship
  ! currently we assume spherical particles for cloud ice/snow
  ! m = cD^d
  ! cloud ice mass-diameter relationship
  real(r8), parameter :: ci = rhoi*pi/6._r8

  ! local parameters
  real(r8), parameter :: lammaxi = 1._r8/10.e-6_r8
  real(r8) :: lammini

!!== KZ_DCS 
  if(dcs_tdep) then 
     lammini = 1._r8/(2._r8*dcst)
  else
     lammini = 1._r8/(2._r8*dcs)
  end if 
!!== KZ_DCS 

  if (qiic > qsmall) then

     ! add upper limit to in-cloud number concentration to prevent numerical error
     niic = min(niic, qiic * 1.e20_r8)

     lami = (cons1*ci*niic/qiic)**(1._r8/dsph)
     n0i = niic * lami

     ! check for slope
     ! adjust vars
     if (lami < lammini) then
        lami = lammini
        n0i = lami**(dsph+1._r8) * qiic/(ci*cons1)
        niic = n0i/lami
     else if (lami > lammaxi) then
        lami = lammaxi
        n0i = lami**(dsph+1._r8) * qiic/(ci*cons1)
        niic = n0i/lami
     end if
  else
     lami = 0._r8
     n0i  = 0._r8
  end if

end subroutine size_dist_param_ice

! get rain size distribution parameters
elemental subroutine size_dist_param_rain(qric, nric, lamr, n0r)
  real(r8), intent(in) :: qric
  real(r8), intent(inout) :: nric

  real(r8), intent(out) :: lamr
  real(r8), intent(out) :: n0r

  ! particle mass-diameter relationship
  ! currently we assume spherical particles for cloud ice/snow
  ! m = cD^d
  ! rain mass-diameter relationship
  ! Rain is hard-coded for spherical drops
  ! Therefore cr is rhow*pi, rather than
  ! using rhow*pi/6 and later multiplying
  ! by gamma(1+dsph) == 6.
  real(r8), parameter :: cr = rhow*pi

  ! local parameters
  real(r8), parameter :: lammaxr = 1._r8/20.e-6_r8
  real(r8), parameter :: lamminr = 1._r8/500.e-6_r8

  if (qric > qsmall) then

     lamr = (cr*nric/qric)**(1._r8/3._r8)
     n0r = nric * lamr

     ! check for slope
     ! adjust vars

     if (lamr < lamminr) then
        lamr = lamminr
        n0r = lamr**4 * qric/cr
        nric = n0r/lamr
     else if (lamr > lammaxr) then
        lamr = lammaxr
        n0r = lamr**4 * qric/cr
        nric = n0r/lamr
     end if

  else
     lamr = 0._r8
     n0r  = 0._r8
  end if

end subroutine size_dist_param_rain

! get snow size distribution parameters
elemental subroutine size_dist_param_snow(qsic, nsic, lams, n0s)
  real(r8), intent(in) :: qsic
  real(r8), intent(inout) :: nsic

  real(r8), intent(out) :: lams
  real(r8), intent(out) :: n0s

  ! particle mass-diameter relationship
  ! currently we assume spherical particles for cloud ice/snow
  ! m = cD^d
  ! cloud ice mass-diameter relationship
  real(r8), parameter :: cs = rhosn*pi/6._r8

  ! local parameters
  real(r8), parameter :: lammaxs = 1._r8/10.e-6_r8
  real(r8), parameter :: lammins = 1._r8/2000.e-6_r8

  if (qsic > qsmall) then

     lams = (cons1*cs*nsic/qsic)**(1._r8/dsph)
     n0s = nsic * lams

     ! check for slope
     ! adjust vars
     if (lams < lammins) then
        lams = lammins
        n0s = lams**(dsph+1._r8) * qsic/(cs*cons1)
        nsic = n0s/lams
     else if (lams > lammaxs) then
        lams = lammaxs
        n0s = lams**(dsph+1._r8) * qsic/(cs*cons1)
        nsic = n0s/lams
     end if

  else
     lams = 0._r8
     n0s  = 0._r8
  end if

end subroutine size_dist_param_snow

real(r8) elemental function avg_diameter(q, n, rho_air, rho_sub)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  ! Assumes that diameter follows an exponential distribution.
  real(r8), intent(in) :: q         ! mass mixing ratio
  real(r8), intent(in) :: n         ! number concentration (per volume)
  real(r8), intent(in) :: rho_air   ! local density of the air
  real(r8), intent(in) :: rho_sub   ! density of the particle substance

  avg_diameter = (pi * rho_sub * n/(q*rho_air))**(-1._r8/3._r8)

end function avg_diameter

real(r8) elemental function var_coef(relvar, a)
  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  real(r8), intent(in) :: relvar
  real(r8), intent(in) :: a

  var_coef = gamma(relvar + a) / (gamma(relvar) * relvar**a)

end function var_coef

!========================================================================
!MICROPHYSICAL PROCESS CALCULATIONS
!========================================================================

!========================================================================
! Initial ice deposition and sublimation loop.
! Run before the main loop
! This subroutine written by Peter Caldwell

elemental subroutine ice_deposition_sublimation(deltat, t, qv, qc, qi, ni, lcldm, &
!!== KZ_DCS 
                                                icldm, naai, rho, dv,qvl, qvi, dcst, &
!!== KZ_DCS 
                                                berg, vap_dep, ice_sublim)

  !INPUT VARS:
  !===============================================
  real(r8), intent(in) :: deltat
  real(r8), intent(in) :: t
  real(r8), intent(in) :: qv
  real(r8), intent(in) :: qc
  real(r8), intent(in) :: qi
  real(r8), intent(in) :: ni
  real(r8), intent(in) :: lcldm
  real(r8), intent(in) :: icldm
  real(r8), intent(in) :: naai
  real(r8), intent(in) :: rho
  real(r8), intent(in) :: dv
  real(r8), intent(in) :: qvl
  real(r8), intent(in) :: qvi
!!== KZ_DCS 
  real(r8), intent(in) :: dcst
!!== KZ_DCS 

  !OUTPUT VARS:
  !===============================================
  real(r8), intent(out) :: vap_dep !ice deposition (cell-ave value)
  real(r8), intent(out) :: ice_sublim !ice sublimation (cell-ave value)
  real(r8), intent(out) :: berg !bergeron enhancement (cell-ave value)

  !INTERNAL VARS:
  !===============================================
  real(r8) :: ab
  real(r8) :: epsi
  real(r8) :: qiic
  real(r8) :: niic
  real(r8) :: dum
  real(r8) :: lami
  real(r8) :: n0i

  if (qi>=qsmall) then

     !GET IN-CLOUD qi, ni
     !===============================================
     qiic = qi/icldm
     niic = ni/icldm

     !Compute linearized condensational heating correction
     ab=calc_ab(t, qvi, xxls)

     !Get slope and intercept of gamma distn for ice.
!!== KZ_DCS
     call size_dist_param_ice(qiic, dcst, niic, lami, n0i)
!!== KZ_DCS
     !Get depletion timescale=1/eps
     epsi = 2._r8*pi*n0i*rho*Dv/(lami*lami)

     !Compute deposition/sublimation
     vap_dep = epsi/ab*(qv - qvi)

     !Make this a grid-averaged quantity
     vap_dep=vap_dep*icldm

     !Split into deposition or sublimation.
     if (t<273.15_r8 .and. vap_dep>0._r8) then
        ice_sublim=0._r8
     else
     !hm, make ice_sublim negative for consistency with other evap/sub processes
        ice_sublim=min(vap_dep,0._r8)
        vap_dep=0._r8
     end if

     !sublimation occurs @ any T. Not so for berg.
     if (T<273.15_r8) then

        !Compute bergeron rate assuming cloud for whole step.
        berg = max(epsi/ab*(qvl - qvi), 0._r8)
     else !T>frz
        berg=0._r8
     end if !T<frz

  else !where qi<qsmall
     berg=0._r8
     vap_dep=0._r8
     ice_sublim=0._r8
  end if !qi>qsmall

end subroutine ice_deposition_sublimation

!========================================================================
! autoconversion of cloud liquid water to rain
! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
! minimum qc of 1 x 10^-8 prevents floating point error

elemental subroutine kk2000_liq_autoconversion(qcic, ncic, rho, relvar, &
                        prc, nprc, nprc1)

  real(r8), intent(in) :: qcic
  real(r8), intent(in) :: ncic
  real(r8), intent(in) :: rho
  real(r8), intent(in) :: relvar

  real(r8), intent(out) :: prc
  real(r8), intent(out) :: nprc
  real(r8), intent(out) :: nprc1

  real(r8) :: prc_coef
  real(r8) :: nprc_denom

  ! subcolumn modifications change coefficients
  if (microp_uniform) then
     prc_coef = 1._r8
     nprc_denom = (4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
  else
     prc_coef = var_coef(relvar, 2.47_r8)
     nprc_denom = cons22
  end if

  if (qcic .ge. icsmall) then

     ! nprc is increase in rain number conc due to autoconversion
     ! nprc1 is decrease in cloud droplet conc due to autoconversion

     ! assume exponential sub-grid distribution of qc, resulting in additional
     ! factor related to qcvar below
     ! hm switch for sub-columns, don't include sub-grid qc

     prc = prc_coef * &
          1350._r8 * qcic**2.47_r8 * (ncic/1.e6_r8*rho)**(-1.79_r8)
     nprc = prc/nprc_denom
     nprc1 = prc/(qcic/ncic)

  else
     prc=0._r8
     nprc=0._r8
     nprc1=0._r8
  end if

end subroutine kk2000_liq_autoconversion

!========================================================================
! Autoconversion of cloud ice to snow
! similar to Ferrier (1994)

!!== KZ_DCS 
elemental subroutine ice_autoconversion(t, qiic, lami, n0i, dcst, prci, nprci)
!!== KZ_DCS 

  real(r8), intent(in) :: t
  real(r8), intent(in) :: qiic
  real(r8), intent(in) :: lami
  real(r8), intent(in) :: n0i
!!== KZ_DCS 
  real(r8), intent(in) :: dcst
!!== KZ_DCS 

  real(r8), intent(out) :: prci
  real(r8), intent(out) :: nprci

  if (t .le. tmelt .and.qiic.ge.qsmall) then

     ! note: assumes autoconversion timescale of 180 sec

!!== KZ_DCS 
     if(dcs_tdep) then 
        nprci = n0i/(lami*180._r8)*exp(-lami*dcst)

        prci = pi*rhoi*n0i/(6._r8*180._r8)* &
             (dcst**3/lami+3._r8*dcst**2/lami**2+ &
             6._r8*dcst/lami**3+6._r8/lami**4)*exp(-lami*dcst)
!!== KZ_DCS 
     else 
        nprci = n0i/(lami*180._r8)*exp(-lami*dcs)

        prci = pi*rhoi*n0i/(6._r8*180._r8)* &
             (cons23/lami+3._r8*cons24/lami**2+ &
             6._r8*dcs/lami**3+6._r8/lami**4)*exp(-lami*dcs)
     end if 

  else
     prci=0._r8
     nprci=0._r8
  end if

end subroutine ice_autoconversion

! immersion freezing (Bigg, 1953)
!===================================

elemental subroutine immersion_freezing(t, pgam, lamc, cdist1, qcic, &
                        relvar, mnuccc, nnuccc)

  ! Temperature
  real(r8), intent(in) :: t

  ! Cloud droplet size distribution parameters
  real(r8), intent(in) :: pgam
  real(r8), intent(in) :: lamc
  real(r8), intent(in) :: cdist1

  ! MMR of in-cloud liquid water
  real(r8), intent(in) :: qcic

  ! Relative variance of cloud water
  real(r8), intent(in) :: relvar

  ! Output tendencies
  real(r8), intent(out) :: mnuccc ! MMR
  real(r8), intent(out) :: nnuccc ! Number

  ! Coefficients that will be omitted for sub-columns
  real(r8) :: dum, dum1


  if (microp_uniform) then
     dum = 1._r8
     dum1 = 1._r8
  else
     dum = var_coef(relvar, 2._r8)
     dum1 = var_coef(relvar, 1._r8)
  end if

  if (qcic >= qsmall .and. t < 269.15_r8) then

     mnuccc = dum * &
          pi*pi/36._r8*rhow* &
          cdist1*gamma(7._r8+pgam)* &
          bimm*(exp(aimm*(tmelt - t))-1._r8)/lamc**3/lamc**3

     nnuccc = dum1 * &
          pi/6._r8*cdist1*gamma(pgam+4._r8) &
          *bimm*(exp(aimm*(tmelt - t))-1._r8)/lamc**3

  else
     mnuccc = 0._r8
     nnuccc = 0._r8
  end if ! qcic > qsmall and t < 4 deg C

end subroutine immersion_freezing

! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
!===================================================================
! dust size and number in multiple bins are read in from companion routine

pure subroutine contact_freezing (t, p, rndst, nacon, pgam, lamc, cdist1, qcic, &
     relvar, mnucct, nnucct)

  real(r8), intent(in) :: t(:)            ! Temperature
  real(r8), intent(in) :: p(:)            ! Pressure
  real(r8), intent(in) :: rndst(:,:)      ! Radius (for multiple dust bins)
  real(r8), intent(in) :: nacon(:,:)      ! Number (for multiple dust bins)

  ! Size distribution parameters for cloud droplets
  real(r8), intent(in) :: pgam(:)
  real(r8), intent(in) :: lamc(:)
  real(r8), intent(in) :: cdist1(:)

  ! MMR of in-cloud liquid water
  real(r8), intent(in) :: qcic(:)

  ! Relative cloud water variance
  real(r8), intent(in) :: relvar(:)

  ! Output tendencies
  real(r8), intent(out) :: mnucct(:) ! MMR
  real(r8), intent(out) :: nnucct(:) ! Number

  real(r8) :: tcnt                  ! scaled relative temperature
  real(r8) :: viscosity             ! temperature-specific viscosity (kg/m/s)
  real(r8) :: mfp                   ! temperature-specific mean free path (m)

  ! Dimension these according to number of dust bins, inferred from rndst size
  real(r8) :: nslip(size(rndst,2))  ! slip correction factors
  real(r8) :: ndfaer(size(rndst,2)) ! aerosol diffusivities (m^2/sec)

  ! Coefficients not used for subcolumns
  real(r8) :: dum, dum1

  integer  :: i

  ! subcolumns

  do i = 1,size(t)

     if (qcic(i) >= qsmall .and. t(i) < 269.15_r8) then

        if (microp_uniform) then
           dum = 1._r8
           dum1 = 1._r8
        else
           dum = var_coef(relvar(i), 4._r8/3._r8)
           dum1 = var_coef(relvar(i), 1._r8/3._r8)
        endif

        tcnt=(270.16_r8-t(i))**1.3_r8
        viscosity = 1.8e-5_r8*(t(i)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
        mfp = 2.0_r8*viscosity/ &                         ! Mean free path (m)
                     (p(i)*sqrt( 8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i)) ))

        ! Note that these two are vectors.
        nslip = 1.0_r8+(mfp/rndst(i,:))*(1.257_r8+(0.4_r8*exp(-(1.1_r8*rndst(i,:)/mfp))))! Slip correction factor
        ndfaer = 1.381e-23_r8*t(i)*nslip/(6._r8*pi*viscosity*rndst(i,:))  ! aerosol diffusivity (m2/s)

        mnucct(i) = dum *  &
             dot_product(ndfaer,nacon(i,:)*tcnt)*pi*pi/3._r8*rhow* &
             cdist1(i)*gamma(pgam(i)+5._r8)/lamc(i)**4

        nnucct(i) =  dum1 *  &
             dot_product(ndfaer,nacon(i,:)*tcnt)*2._r8*pi*  &
             cdist1(i)*gamma(pgam(i)+2._r8)/lamc(i)

     else

        mnucct(i)=0._r8
        nnucct(i)=0._r8

     end if ! qcic > qsmall and t < 4 deg C
  end do

end subroutine contact_freezing

! snow self-aggregation from passarelli, 1978, used by reisner, 1998
!===================================================================
! this is hard-wired for bs = 0.4 for now
! ignore self-collection of cloud ice

elemental subroutine snow_self_aggregation(t, rho, asn, qsic, nsic, nsagg)

  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: rho  ! Density
  real(r8), intent(in) :: asn  ! fall speed parameter for snow

  ! In-cloud snow
  real(r8), intent(in) :: qsic ! MMR
  real(r8), intent(in) :: nsic ! Number

  ! Output number tendency
  real(r8), intent(out) :: nsagg

  if (qsic >= qsmall .and. t <= tmelt) then
     nsagg = -1108._r8*asn*eii* &
          pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)* &
          rho**((2._r8+bs)/3._r8)*qsic**((2._r8+bs)/3._r8)* &
          (nsic*rho)**((4._r8-bs)/3._r8) /(4._r8*720._r8*rho)
  else
     nsagg=0._r8
  end if

end subroutine snow_self_aggregation

! accretion of cloud droplets onto snow/graupel
!===================================================================
! here use continuous collection equation with
! simple gravitational collection kernel
! ignore collisions between droplets/cloud ice
! since minimum size ice particle for accretion is 50 - 150 micron

elemental subroutine accrete_cloud_water_snow(t, rho, asn, uns, mu, qcic, ncic, qsic, &
     pgam, lamc, lams, n0s, psacws, npsacws)

  real(r8), intent(in) :: t   ! Temperature
  real(r8), intent(in) :: rho ! Density
  real(r8), intent(in) :: asn ! Fallspeed parameter (snow)
  real(r8), intent(in) :: uns ! Current fallspeed   (snow)
  real(r8), intent(in) :: mu  ! Viscosity

  ! In-cloud liquid water
  real(r8), intent(in) :: qcic ! MMR
  real(r8), intent(in) :: ncic ! Number

  ! In-cloud snow
  real(r8), intent(in) :: qsic ! MMR

  ! Cloud droplet size parameters
  real(r8), intent(in) :: pgam
  real(r8), intent(in) :: lamc

  ! Snow size parameters
  real(r8), intent(in) :: lams
  real(r8), intent(in) :: n0s

  ! Output tendencies
  real(r8), intent(out) :: psacws  ! Mass mixing ratio
  real(r8), intent(out) :: npsacws ! Number concentration

  real(r8) :: dc0 ! Provisional mean droplet size
  real(r8) :: dum
  real(r8) :: eci ! collection efficiency for riming of snow by droplets

  ! ignore collision of snow with droplets above freezing

  if (qsic >= qsmall .and. t <= tmelt .and. qcic >= qsmall) then

     ! put in size dependent collection efficiency
     ! mean diameter of snow is area-weighted, since
     ! accretion is function of crystal geometric area
     ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)

     dc0 = (pgam+1._r8)/lamc
     dum = dc0*dc0*uns*rhow/(9._r8*mu*(1._r8/lams))
     eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))

     eci = max(eci,0._r8)
     eci = min(eci,1._r8)

     ! no impact of sub-grid distribution of qc since psacws
     ! is linear in qc

     psacws = pi/4._r8*asn*qcic*rho*n0s*eci*cons11 / lams**(bs+3._r8)
     npsacws = pi/4._r8*asn*ncic*rho*n0s*eci*cons11 / lams**(bs+3._r8)
  else
     psacws = 0._r8
     npsacws = 0._r8
  end if

end subroutine accrete_cloud_water_snow

! add secondary ice production due to accretion of droplets by snow
!===================================================================
! (Hallet-Mossop process) (from Cotton et al., 1986)

elemental subroutine secondary_ice_production(t, psacws, msacwi, nsacwi)
  real(r8), intent(in) :: t ! Temperature

  ! Accretion of cloud water to snow tendencies
  real(r8), intent(inout) :: psacws ! MMR

  ! Output (ice) tendencies
  real(r8), intent(out) :: msacwi ! MMR
  real(r8), intent(out) :: nsacwi ! Number

  if((t < 270.16_r8) .and. (t >= 268.16_r8)) then
     nsacwi = 3.5e8_r8*(270.16_r8-t)/2.0_r8*psacws
     msacwi = min(nsacwi*mi0, psacws)
  else if((t < 268.16_r8) .and. (t >= 265.16_r8)) then
     nsacwi = 3.5e8_r8*(t-265.16_r8)/3.0_r8*psacws
     msacwi = min(nsacwi*mi0, psacws)
  else
     nsacwi = 0.0_r8
     msacwi = 0.0_r8
  endif

  psacws = max(0.0_r8,psacws - nsacwi*mi0)

end subroutine secondary_ice_production

! accretion of rain water by snow
!===================================================================
! formula from ikawa and saito, 1991, used by reisner et al., 1998

elemental subroutine accrete_rain_snow(t, rho, umr, ums, unr, uns, qric, qsic, &
     lamr, n0r, lams, n0s, pracs, npracs )

  real(r8), intent(in) :: t   ! Temperature
  real(r8), intent(in) :: rho ! Density

  ! Fallspeeds
  ! mass-weighted
  real(r8), intent(in) :: umr ! rain
  real(r8), intent(in) :: ums ! snow
  ! number-weighted
  real(r8), intent(in) :: unr ! rain
  real(r8), intent(in) :: uns ! snow

  ! In cloud MMRs
  real(r8), intent(in) :: qric ! rain
  real(r8), intent(in) :: qsic ! snow

  ! Size distribution parameters
  ! rain
  real(r8), intent(in) :: lamr
  real(r8), intent(in) :: n0r
  ! snow
  real(r8), intent(in) :: lams
  real(r8), intent(in) :: n0s

  ! Output tendencies
  real(r8), intent(out) :: pracs  ! MMR
  real(r8), intent(out) :: npracs ! Number

  ! Collection efficiency for accretion of rain by snow
  real(r8), parameter :: ecr = 1.0_r8

  if (qric >= icsmall .and. qsic >= icsmall .and. t <= tmelt) then

     pracs = pi*pi*ecr*(((1.2_r8*umr-0.95_r8*ums)**2 + &
          0.08_r8*ums*umr)**0.5_r8 *  &
          rhow * rho * n0r * n0s * &
          (5._r8/(lamr**6 * lams)+ &
          2._r8/(lamr**5 * lams**2)+ &
          0.5_r8/(lamr**4 * lams**3)))

     npracs = pi/2._r8*rho*ecr* (1.7_r8*(unr-uns)**2 + &
          0.3_r8*unr*uns)**0.5_r8 * &
          n0r*n0s* &
          (1._r8/(lamr**3 * lams)+ &
          1._r8/(lamr**2 * lams**2)+ &
          1._r8/(lamr * lams**3))

  else
     pracs = 0._r8
     npracs = 0._r8
  end if

end subroutine accrete_rain_snow

! heterogeneous freezing of rain drops
!===================================================================
! follows from Bigg (1953)

elemental subroutine heterogeneous_rain_freezing(t, qric, nric, lamr, mnuccr, nnuccr)

  real(r8), intent(in) :: t    ! Temperature

  ! In-cloud rain
  real(r8), intent(in) :: qric ! MMR
  real(r8), intent(in) :: nric ! Number
  real(r8), intent(in) :: lamr ! size parameter

  ! Output tendencies
  real(r8), intent(out) :: mnuccr ! MMR
  real(r8), intent(out) :: nnuccr ! Number

  if (t < 269.15_r8 .and. qric >= qsmall) then

     ! Division by lamr**3 twice is old workaround to avoid overflow.
     ! Probably no longer necessary
     mnuccr = 20._r8*pi*pi*rhow*nric*bimm* &
          (exp(aimm*(tmelt - t))-1._r8)/lamr**3 &
          /lamr**3

     nnuccr = pi*nric*bimm* &
          (exp(aimm*(tmelt - t))-1._r8)/lamr**3
  else
     mnuccr = 0._r8
     nnuccr = 0._r8
  end if
end subroutine heterogeneous_rain_freezing

! accretion of cloud liquid water by rain
!===================================================================
! formula from Khrouditnov and Kogan (2000)
! gravitational collection kernel, droplet fall speed neglected

elemental subroutine accrete_cloud_water_rain(qric, qcic, ncic, &
                     relvar, accre_enhan, pra, npra)

  ! In-cloud rain
  real(r8), intent(in) :: qric ! MMR

  ! Cloud droplets
  real(r8), intent(in) :: qcic ! MMR
  real(r8), intent(in) :: ncic ! Number

  ! SGS variability
  real(r8), intent(in) :: relvar
  real(r8), intent(in) :: accre_enhan

  ! Output tendencies
  real(r8), intent(out) :: pra  ! MMR
  real(r8), intent(out) :: npra ! Number

  ! Coefficient that varies for subcolumns
  real(r8) :: pra_coef

  if (microp_uniform) then
     pra_coef = 1._r8
  else
     pra_coef = accre_enhan * var_coef(relvar, 1.15_r8)
  end if

  if (qric >= qsmall .and. qcic >= qsmall) then

     ! include sub-grid distribution of cloud water
     pra = pra_coef * 67._r8*(qcic*qric)**1.15_r8

     npra = pra/(qcic/ncic)

  else
     pra = 0._r8
     npra = 0._r8
  end if
end subroutine accrete_cloud_water_rain

! Self-collection of rain drops
!===================================================================
! from Beheng(1994)

elemental subroutine self_collection_rain(rho, qric, nric, nragg)

  real(r8), intent(in) :: rho  ! Air density

  ! Rain
  real(r8), intent(in) :: qric ! MMR
  real(r8), intent(in) :: nric ! Number

  ! Output number tendency
  real(r8), intent(out) :: nragg

  if (qric >= qsmall) then
     nragg = -8._r8*nric*qric*rho
  else
     nragg = 0._r8
  end if

end subroutine self_collection_rain

! Accretion of cloud ice by snow
!===================================================================
! For this calculation, it is assumed that the Vs >> Vi
! and Ds >> Di for continuous collection

elemental subroutine accrete_cloud_ice_snow(t, rho, asn, qiic, niic, qsic, &
     lams, n0s, prai, nprai)

  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: rho  ! Density

  real(r8), intent(in) :: asn  ! Snow fallspeed parameter

  ! Cloud ice
  real(r8), intent(in) :: qiic ! MMR
  real(r8), intent(in) :: niic ! Number

  real(r8), intent(in) :: qsic ! Snow MMR

  ! Snow size parameters
  real(r8), intent(in) :: lams
  real(r8), intent(in) :: n0s

  ! Output tendencies
  real(r8), intent(out) :: prai  ! MMR
  real(r8), intent(out) :: nprai ! Number

  if (qsic >= qsmall .and. qiic >= qsmall .and. t <= tmelt) then

     prai = pi/4._r8 * asn * qiic * rho * n0s * eii * cons11/ &
          lams**(bs+3._r8)

     nprai = pi/4._r8 * asn * niic * rho * n0s * eii * cons11/ &
          lams**(bs+3._r8)
  else
     prai = 0._r8
     nprai = 0._r8
  end if

end subroutine accrete_cloud_ice_snow

! calculate evaporation/sublimation of rain and snow
!===================================================================
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process

elemental subroutine evaporate_sublimate_precip(deltat, t, p, rho, dv, mu, sc, q, qvl, qvi, &
     lcldm, cldmax, arn, asn, qcic, qiic, qric, qsic, lamr, n0r, lams, n0s, &
     pre, prds)

  real(r8), intent(in) :: deltat ! timestep

  real(r8), intent(in) :: t    ! temperature
  real(r8), intent(in) :: p    ! pressure
  real(r8), intent(in) :: rho  ! air density
  real(r8), intent(in) :: dv   ! water vapor diffusivity
  real(r8), intent(in) :: mu   ! viscosity
  real(r8), intent(in) :: sc   ! schmidt number
  real(r8), intent(in) :: q    ! humidity
  real(r8), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), intent(in) :: qvi  ! saturation humidity (ice)
  real(r8), intent(in) :: lcldm  ! liquid cloud fraction
  real(r8), intent(in) :: cldmax ! precipitation fraction (maximum overlap)

  ! fallspeed parameters
  real(r8), intent(in) :: arn  ! rain
  real(r8), intent(in) :: asn  ! snow

  ! In-cloud MMRs
  real(r8), intent(in) :: qcic ! cloud liquid
  real(r8), intent(in) :: qiic ! cloud ice
  real(r8), intent(in) :: qric ! rain
  real(r8), intent(in) :: qsic ! snow

  ! Size parameters
  ! rain
  real(r8), intent(in) :: lamr
  real(r8), intent(in) :: n0r
  ! snow
  real(r8), intent(in) :: lams
  real(r8), intent(in) :: n0s

  ! Output tendencies
  real(r8), intent(out) :: pre
  real(r8), intent(out) :: prds

  ! checking for RH after rain evap
  real(r8) :: esn    ! saturation pressure
  real(r8) :: qvn    ! saturation humidity

  real(r8) :: qclr   ! water vapor mixing ratio in clear air
  real(r8) :: ab     ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale

  ! Temps/dummies
  real(r8) :: qtmp
  real(r8) :: ttmp

  real(r8) :: dum, dum1

  ! set temporary cloud fraction to zero if cloud water + ice is very small
  ! this will ensure that evaporation/sublimation of precip occurs over
  ! entire grid cell, since min cloud fraction is specified otherwise
  if (qcic+qiic < 1.e-6_r8) then
     dum = 0._r8
  else
     dum = lcldm
  end if

  ! only calculate if there is some precip fraction > cloud fraction

  if (cldmax > dum) then

     ! calculate q for out-of-cloud region
     qvn = min(qvl,1._r8)
     qclr=(q-dum*qvn)/(1._r8-dum)

     ! evaporation of rain
     if (qric.ge.qsmall) then

        ab = calc_ab(t, qvl, xxlv)
        eps = 2._r8*pi*n0r*rho*Dv* &
             (f1r/(lamr*lamr)+ &
             f2r*(arn*rho/mu)**0.5_r8* &
             sc**(1._r8/3._r8)*cons13/ &
             (lamr**(5._r8/2._r8+br/2._r8)))

        pre = eps*(qclr-qvl)/ab

        ! only evaporate in out-of-cloud region
        ! and distribute across cldmax
        pre=min(pre*(cldmax-dum),0._r8)
        pre=pre/cldmax
     else
        pre = 0._r8
     end if

     ! sublimation of snow
     if (qsic.ge.qsmall) then
        ab = calc_ab(t, qvi, xxls)
        eps = 2._r8*pi*n0s*rho*Dv* &
             (f1s/(lams*lams)+ &
             f2s*(asn*rho/mu)**0.5_r8* &
             sc**(1._r8/3._r8)*cons14/ &
             (lams**(5._r8/2._r8+bs/2._r8)))
        prds = eps*(qclr-qvi)/ab

        ! only sublimate in out-of-cloud region and distribute over cldmax
        prds=min(prds*(cldmax-dum),0._r8)
        prds=prds/cldmax
     else
        prds = 0._r8
     end if

  else
     prds = 0._r8
     pre = 0._r8
  end if

end subroutine evaporate_sublimate_precip

! bergeron process - evaporation of droplets and deposition onto snow
!===================================================================

elemental subroutine bergeron_process_snow(t, rho, dv, mu, sc, qvl, qvi, asn, &
     qcic, qsic, lams, n0s, bergs)

  real(r8), intent(in) :: t    ! temperature
  real(r8), intent(in) :: rho  ! air density
  real(r8), intent(in) :: dv   ! water vapor diffusivity
  real(r8), intent(in) :: mu   ! viscosity
  real(r8), intent(in) :: sc   ! schmidt number
  real(r8), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), intent(in) :: qvi  ! saturation humidity (ice)

  ! fallspeed parameter for snow
  real(r8), intent(in) :: asn

  ! In-cloud MMRs
  real(r8), intent(in) :: qcic ! cloud liquid
  real(r8), intent(in) :: qsic ! snow

  ! Size parameters for snow
  real(r8), intent(in) :: lams
  real(r8), intent(in) :: n0s

  ! Output tendencies
  real(r8), intent(out) :: bergs

  real(r8) :: ab     ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale

  if (qsic >= qsmall.and. qcic >= qsmall .and. t < tmelt) then
     ab = calc_ab(t, qvi, xxls)
     eps = 2._r8*pi*n0s*rho*Dv* &
          (f1s/(lams*lams)+ &
          f2s*(asn*rho/mu)**0.5_r8* &
          sc**(1._r8/3._r8)*cons14/ &
          (lams**(5._r8/2._r8+bs/2._r8)))
     bergs = eps*(qvl-qvi)/ab
  else
     bergs = 0._r8
  end if

end subroutine bergeron_process_snow

!========================================================================
!OUTPUT CALCULATIONS
!========================================================================

elemental subroutine calc_rercld(lamr, n0r, lamc, cdist1, pgam, dumr, qcic, &
     rercld)
  real(r8), intent(in) :: lamr          ! rain size parameter (slope)
  real(r8), intent(in) :: n0r           ! rain size parameter (intercept)
  real(r8), intent(in) :: lamc          ! size distribution parameter (slope)
  real(r8), intent(in) :: cdist1        ! for droplet freezing
  real(r8), intent(in) :: pgam          ! droplet size parameter
  real(r8), intent(in) :: dumr          ! in-cloud rain mass mixing ratio
  real(r8), intent(in) :: qcic          ! in-cloud cloud liquid

  real(r8), intent(inout) :: rercld     ! effective radius calculation for rain + cloud

  ! combined size of precip & cloud drops
  real(r8) :: Atmp

  ! Rain drops
  if (lamr > 0._r8) then
     Atmp = n0r * pi / (2._r8 * lamr**3._r8)
  else
     Atmp = 0._r8
  end if

  ! Add cloud drops
  if (lamc > 0._r8) then
     Atmp = Atmp + cdist1 * pi * gamma(pgam+3._r8)/(4._r8 * lamc**2._r8)
  end if

  if (Atmp > 0._r8) then
     rercld = rercld + 3._r8 *(dumr + qcic) / (4._r8 * rhow * Atmp)
  end if

end subroutine calc_rercld

!========================================================================
!UTILITIES
!========================================================================

pure subroutine micro_mg_get_cols(ncol, nlev, top_lev, qcn, qin, &
     mgncol, mgcols)

  ! Determines which columns microphysics should operate over by
  ! checking for non-zero cloud water/ice.

  integer, intent(in) :: ncol      ! Number of columns with meaningful data
  integer, intent(in) :: nlev      ! Number of levels to use
  integer, intent(in) :: top_lev   ! Top level for microphysics

  real(r8), intent(in) :: qcn(:,:) ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(:,:) ! cloud ice mixing ratio (kg/kg)

  integer, intent(out) :: mgncol   ! Number of columns MG will use
  integer, allocatable, intent(out) :: mgcols(:) ! column indices

  integer :: lev_offset  ! top_lev - 1 (defined here for consistency)
  logical :: ltrue(ncol) ! store tests for each column

  integer :: i, ii ! column indices

  if (allocated(mgcols)) deallocate(mgcols)

  lev_offset = top_lev - 1

  ! Using "any" along dimension 2 collapses across levels, but
  ! not columns, so we know if water is present at any level
  ! in each column.

  ltrue = any(qcn(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)
  ltrue = ltrue .or. any(qin(:ncol,top_lev:(nlev+lev_offset)) >= qsmall, 2)

  ! Scan for true values to get a usable list of indices.

  mgncol = count(ltrue)
  allocate(mgcols(mgncol))
  i = 0
  do ii = 1,ncol
     if (ltrue(ii)) then
        i = i + 1
        mgcols(i) = ii
     end if
  end do

end subroutine micro_mg_get_cols


!!== KZ_DCS
subroutine get_dcst(ncol,pver,temp,dcst)

implicit none

integer,  intent(in) :: ncol
integer,  intent(in) :: pver                 ! number of layers in columns
real(r8), intent(in) :: temp(ncol,pver)       ! input temperature (K)
real(r8), intent(out) :: dcst(ncol,pver)      ! temperature dependent dcs

integer :: i,k
real(r8) :: st


dcst = 400.e-6_r8

do k=1,pver
   do i=1,ncol
      st = temp(i,k) - 273.15
      if(st.le.-70.) then
         dcst(i,k) = 100.e-6_r8
      elseif(st.gt.-70. .and. st.le.-10.) then
         dcst(i,k) = 5.e-6_r8 * st  + 450.e-6_r8
      elseif(st.gt.-10.) then
         dcst(i,k) = 400.e-6_r8
      end if
   end do
end do

return

end subroutine get_dcst
!!== KZ_DCS





pure function interp_to_mid(orig_val, weights) result(new_val)
  ! Linear interpolation, here used to move from interfaces to midlevel
  real(r8), intent(in) :: orig_val(:,:)
  real(r8), intent(in) :: weights(:,:)

  ! New value is slightly smaller than the old value
  real(r8) :: new_val(size(orig_val,1),size(orig_val,2)-1)

  new_val = orig_val(:,:size(new_val,2))
  new_val = new_val + weights * (orig_val(:,2:) - new_val)

end function interp_to_mid

end module micro_mg1_5
