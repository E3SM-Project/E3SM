module clm_varcon

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing various model constants.
  !
  ! !USES:
  use shr_kind_mod  , only: r8 => shr_kind_r8
  use shr_const_mod , only: SHR_CONST_G,SHR_CONST_STEBOL,SHR_CONST_KARMAN
  use shr_const_mod , only: SHR_CONST_RWV,SHR_CONST_RDAIR,SHR_CONST_CPFW
  use shr_const_mod , only: SHR_CONST_CPICE,SHR_CONST_CPDAIR,SHR_CONST_LATVAP
  use shr_const_mod , only: SHR_CONST_LATSUB,SHR_CONST_LATICE,SHR_CONST_RHOFW
  use shr_const_mod , only: SHR_CONST_RHOICE,SHR_CONST_TKFRZ,SHR_CONST_REARTH
  use shr_const_mod , only: SHR_CONST_PDB, SHR_CONST_PI, SHR_CONST_CDAY
  use shr_const_mod , only: SHR_CONST_RGAS
  use clm_varpar    , only: numrad, nlevgrnd, nlevlak, nlevdecomp_full
  use clm_varpar    , only: ngases
  use clm_varpar    , only: nlayer
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !-----------------------------------------------------------------------
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: clm_varcon_init  ! initialize constants in clm_varcon
  !
  ! !REVISION HISTORY:
  ! Created by Mariana Vertenstein
  ! 27 February 2008: Keith Oleson; Add forcing height and aerodynamic parameters
  !-----------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Initialize mathmatical constants
  !------------------------------------------------------------------

  real(r8) :: rpi    = SHR_CONST_PI

  !------------------------------------------------------------------
  ! Initialize physical constants
  !------------------------------------------------------------------

  real(r8), parameter :: n_melt=0.7                         ! fsca shape parameter
  real(r8), parameter :: e_ice=6.0                          ! soil ice impedance factor
  real(r8), parameter :: pc = 0.4                           ! threshold probability
  real(r8), parameter :: mu = 0.13889                       ! connectivity exponent 
  real(r8) :: grav   = SHR_CONST_G                          ! gravity constant [m/s2]
  real(r8) :: sb     = SHR_CONST_STEBOL                     ! stefan-boltzmann constant  [W/m2/K4]
  real(r8) :: vkc    = SHR_CONST_KARMAN                     ! von Karman constant [-]
  real(r8) :: rwat   = SHR_CONST_RWV                        ! gas constant for water vapor [J/(kg K)]
  real(r8) :: rair   = SHR_CONST_RDAIR                      ! gas constant for dry air [J/kg/K]
  real(r8) :: roverg = SHR_CONST_RWV/SHR_CONST_G*1000._r8   ! Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  real(r8) :: cpliq  = SHR_CONST_CPFW                       ! Specific heat of water [J/kg-K]
  real(r8) :: cpice  = SHR_CONST_CPICE                      ! Specific heat of ice [J/kg-K]
  real(r8) :: cpair  = SHR_CONST_CPDAIR                     ! specific heat of dry air [J/kg/K]
  real(r8) :: hvap   = SHR_CONST_LATVAP                     ! Latent heat of evap for water [J/kg]
  real(r8) :: hsub   = SHR_CONST_LATSUB                     ! Latent heat of sublimation    [J/kg]
  real(r8) :: hfus   = SHR_CONST_LATICE                     ! Latent heat of fusion for ice [J/kg]
  real(r8) :: denh2o = SHR_CONST_RHOFW                      ! density of liquid water [kg/m3]
  real(r8) :: denice = SHR_CONST_RHOICE                     ! density of ice [kg/m3]
  real(r8) :: rgas   = SHR_CONST_RGAS                       ! universal gas constant [J/K/kmole]
  real(r8) :: tkair  = 0.023_r8                             ! thermal conductivity of air   [W/m/K]
  real(r8) :: tkice  = 2.290_r8                             ! thermal conductivity of ice   [W/m/K]
  real(r8) :: tkwat  = 0.57_r8                              ! thermal conductivity of water [W/m/K]
  real(r8) :: tfrz   = SHR_CONST_TKFRZ                      ! freezing temperature [K]
  real(r8), parameter :: tcrit  = 2.5_r8                    ! critical temperature to determine rain or snow
  real(r8) :: o2_molar_const = 0.209_r8                     ! constant atmospheric O2 molar ratio (mol/mol)

  real(r8) :: bdsno = 250._r8                               ! bulk density snow (kg/m**3)
  real(r8) :: alpha_aero = 1.0_r8                           ! constant for aerodynamic parameter weighting
  real(r8) :: tlsai_crit = 2.0_r8                           ! critical value of elai+esai for which aerodynamic parameters are maximum
  real(r8) :: watmin = 0.01_r8                              ! minimum soil moisture (mm)

  real(r8) :: re = SHR_CONST_REARTH*0.001_r8                ! radius of earth (km)

  real(r8), public, parameter :: degpsec = 15._r8/3600.0_r8 ! Degree's earth rotates per second
  real(r8), public, parameter ::  secspday= SHR_CONST_CDAY  ! Seconds per day
  integer,  public, parameter :: isecspday= secspday        ! Integer seconds per day
  real(r8), public, parameter ::  spval = 1.e36_r8          ! special value for real data
  integer , public, parameter :: ispval = -9999             ! special value for int data 
                                                            ! (keep this negative to avoid conflicts with possible valid values)

  ! These are tunable constants from clm2_3

  real(r8) :: zlnd = 0.01_r8        ! Roughness length for soil [m]
  real(r8) :: zsno = 0.0024_r8      ! Roughness length for snow [m]
  real(r8) :: csoilc = 0.004_r8     ! Drag coefficient for soil under canopy [-]
  real(r8) :: capr   = 0.34_r8      ! Tuning factor to turn first layer T into surface T
  real(r8) :: cnfac  = 0.5_r8       ! Crank Nicholson factor between 0 and 1
  real(r8) :: ssi    = 0.033_r8     ! Irreducible water saturation of snow
  real(r8) :: wimp   = 0.05_r8      ! Water impremeable if porosity less than wimp
  real(r8) :: pondmx = 0.0_r8       ! Ponding depth (mm)
  real(r8) :: pondmx_urban = 1.0_r8 ! Ponding depth for urban roof and impervious road (mm)

  real(r8) :: thk_bedrock = 3.0_r8  ! thermal conductivity of 'typical' saturated granitic rock 
                                    ! (Clauser and Huenges, 1995)(W/m/K)

  !!! C13
  real(r8), parameter :: preind_atm_del13c = -6.0   ! preindustrial value for atmospheric del13C
  real(r8), parameter :: preind_atm_ratio = SHR_CONST_PDB + (preind_atm_del13c * SHR_CONST_PDB)/1000.0  ! 13C/12C
  real(r8) :: c13ratio = preind_atm_ratio/(1.0+preind_atm_ratio) ! 13C/(12+13)C preind atmosphere

  !!! C14
  real(r8) :: c14ratio = 1.e-12_r8
  ! real(r8) :: c14ratio = 1._r8  ! debug lets set to 1 to try to avoid numerical errors

  ! Note that the wasteheat factors are currently set to zero until a better parameterization can be developed
  ! The prior parameterization appeared to be significantly overestimating wasteheat
  real(r8) :: ht_wasteheat_factor = 0.0_r8  !wasteheat factor for urban heating (-)
  real(r8) :: ac_wasteheat_factor = 0.0_r8  !wasteheat factor for urban air conditioning (-)
  real(r8) :: wasteheat_limit = 100._r8  !limit on wasteheat (W/m2)

  real(r8), parameter :: h2osno_max = 1000._r8    ! max allowed snow thickness (mm H2O)
  real(r8), parameter :: lapse_glcmec = 0.006_r8  ! surface temperature lapse rate (deg m-1)
                                                  ! Pritchard et al. (GRL, 35, 2008) use 0.006  
  real(r8), parameter :: glcmec_rain_snow_threshold = SHR_CONST_TKFRZ  ! temperature dividing rain & snow in downscaling (K)

  integer, private :: i  ! loop index

 !real(r8), parameter :: nitrif_n2o_loss_frac = 0.02_r8  ! fraction of N lost as N2O in nitrification (Parton et al., 2001)
  real(r8), parameter :: nitrif_n2o_loss_frac = 6.e-4_r8 ! fraction of N lost as N2O in nitrification (Li et al., 2000)
  real(r8), parameter :: frac_minrlztn_to_no3 = 0.2_r8   ! fraction of N mineralized that is dieverted to the nitrification stream (Parton et al., 2001)

  !------------------------------------------------------------------
  ! Set subgrid names
  !------------------------------------------------------------------

  character(len=16), parameter :: grlnd  = 'lndgrid'      ! name of lndgrid
  character(len=16), parameter :: namea  = 'gridcellatm'  ! name of atmgrid
  character(len=16), parameter :: nameg  = 'gridcell'     ! name of gridcells
  character(len=16), parameter :: namel  = 'landunit'     ! name of landunits
  character(len=16), parameter :: namec  = 'column'       ! name of columns
  character(len=16), parameter :: namep  = 'pft'          ! name of patches
  character(len=16), parameter :: nameCohort = 'cohort'   ! name of cohorts (ED specific)

  !------------------------------------------------------------------
  ! Initialize miscellaneous radiation constants
  !------------------------------------------------------------------

  real(r8) :: betads  = 0.5_r8            ! two-stream parameter betad for snow
  real(r8) :: betais  = 0.5_r8            ! two-stream parameter betai for snow
  real(r8) :: omegas(numrad)           ! two-stream parameter omega for snow by band
  data (omegas(i),i=1,numrad) /0.8_r8, 0.4_r8/

  ! Lake Model Constants will be defined in LakeCon.

  !------------------------------------------------------------------
  ! Soil depths are constants for now; lake depths can vary by gridcell
  ! zlak and dzlak correspond to the default 50 m lake depth.
  ! The values for the following arrays are set in routine iniTimeConst
  !------------------------------------------------------------------

  real(r8), allocatable :: zlak(:)         !lake z  (layers)
  real(r8), allocatable :: dzlak(:)        !lake dz (thickness)
  real(r8), allocatable :: zsoi(:)         !soil z  (layers)
  real(r8), allocatable :: dzsoi(:)        !soil dz (thickness)
  real(r8), allocatable :: zisoi(:)        !soil zi (interfaces)
  real(r8), allocatable :: dzsoi_decomp(:) !soil dz (thickness)
  integer , allocatable :: nlvic(:)        !number of CLM layers in each VIC layer (#)
  real(r8), allocatable :: dzvic(:)        !soil dz (thickness) of each VIC layer
  real(r8) ,allocatable :: zsoifl(:)       !original soil midpoint (used in interpolation of sand and clay)
  real(r8) ,allocatable :: zisoifl(:)      !original soil interface depth (used in interpolation of sand and clay)
  real(r8) ,allocatable :: dzsoifl(:)      !original soil thickness  (used in interpolation of sand and clay)

  !------------------------------------------------------------------
  ! (Non-tunable) Constants for the CH4 submodel (Tuneable constants in ch4varcon)
  !------------------------------------------------------------------
  ! Note some of these constants are also used in CNNitrifDenitrifMod

  real(r8), parameter :: catomw = 12.011_r8 ! molar mass of C atoms (g/mol)

  real(r8) :: s_con(ngases,4)    ! Schmidt # calculation constants (spp, #)
  data (s_con(1,i),i=1,4) /1898_r8, -110.1_r8, 2.834_r8, -0.02791_r8/ ! CH4
  data (s_con(2,i),i=1,4) /1801_r8, -120.1_r8, 3.7818_r8, -0.047608_r8/ ! O2
  data (s_con(3,i),i=1,4) /1911_r8, -113.7_r8, 2.967_r8, -0.02943_r8/ ! CO2

  real(r8) :: d_con_w(ngases,3)    ! water diffusivity constants (spp, #)  (mult. by 10^-4)
  data (d_con_w(1,i),i=1,3) /0.9798_r8, 0.02986_r8, 0.0004381_r8/ ! CH4
  data (d_con_w(2,i),i=1,3) /1.172_r8, 0.03443_r8, 0.0005048_r8/ ! O2
  data (d_con_w(3,i),i=1,3) /0.939_r8, 0.02671_r8, 0.0004095_r8/ ! CO2

  real(r8) :: d_con_g(ngases,2)    ! gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)
  data (d_con_g(1,i),i=1,2) /0.1875_r8, 0.0013_r8/ ! CH4
  data (d_con_g(2,i),i=1,2) /0.1759_r8, 0.00117_r8/ ! O2
  data (d_con_g(3,i),i=1,2) /0.1325_r8, 0.0009_r8/ ! CO2

  real(r8) :: c_h_inv(ngases)    ! constant (K) for Henry's law (4.12, Wania)
  data c_h_inv(1:3) /1600._r8, 1500._r8, 2400._r8/ ! CH4, O2, CO2

  real(r8) :: kh_theta(ngases)    ! Henry's constant (L.atm/mol) at standard temperature (298K)
  data kh_theta(1:3) /714.29_r8, 769.23_r8, 29.4_r8/ ! CH4, O2, CO2

  real(r8) :: kh_tbase = 298._r8 ! base temperature for calculation of Henry's constant (K)
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------------
  subroutine clm_varcon_init()
    !
    ! !DESCRIPTION:
    ! This subroutine initializes constant arrays in clm_varcon. 
    ! MUST be called  after clm_varpar_init.
    !
    ! USES
    use clm_varpar, only: nlevgrnd, nlevlak, nlevdecomp_full, nlevsoifl, nlayer
    !------------------------------------------------------------------------------

    allocate( zlak(1:nlevlak                 ))
    allocate( dzlak(1:nlevlak                ))
    allocate( zsoi(1:nlevgrnd                ))
    allocate( dzsoi(1:nlevgrnd               ))
    allocate( zisoi(0:nlevgrnd               ))
    allocate( dzsoi_decomp(1:nlevdecomp_full ))
    allocate( nlvic(1:nlayer                 ))
    allocate( dzvic(1:nlayer                 ))
    allocate( zsoifl(1:nlevsoifl             ))
    allocate( zisoifl(0:nlevsoifl            ))
    allocate( dzsoifl(1:nlevsoifl            ))

  end subroutine clm_varcon_init

end module clm_varcon
