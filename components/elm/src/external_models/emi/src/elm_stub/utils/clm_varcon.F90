module clm_varcon

!-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing various model constants. This includes physical constants, landunit
  ! and column indices, and others. This also contains functions for operating on these
  ! constants and associated variables.
  !
  ! TODO: Move some of the stuff here into column_varcon.F90 and landunit_varcon.F90 (bug 1928).
  !
  ! !USES:

  use shr_kind_mod , only: r8 => shr_kind_r8
  use shr_const_mod, only: SHR_CONST_G,SHR_CONST_STEBOL,SHR_CONST_KARMAN,     &
                           SHR_CONST_RWV,SHR_CONST_RDAIR,SHR_CONST_CPFW,      &
                           SHR_CONST_CPICE,SHR_CONST_CPDAIR,SHR_CONST_LATVAP, &
                           SHR_CONST_LATSUB,SHR_CONST_LATICE,SHR_CONST_RHOFW, &
                           SHR_CONST_RHOICE,SHR_CONST_TKFRZ,SHR_CONST_REARTH, &
                           SHR_CONST_PDB, SHR_CONST_PI, SHR_CONST_CDAY,       &
                           SHR_CONST_RGAS

  use clm_varpar, only : ngases

  !
  ! !PUBLIC TYPES:
  implicit none
  save
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

  real(r8), parameter :: n_melt=0.7     !fsca shape parameter
  real(r8), parameter :: e_ice=6.0      !soil ice impedance factor
  real(r8), parameter :: mu = 0.13889   !connectivity exponent
  real(r8) :: grav   = SHR_CONST_G      !gravity constant [m/s2]
  real(r8) :: sb     = SHR_CONST_STEBOL !stefan-boltzmann constant  [W/m2/K4]
  real(r8) :: vkc    = SHR_CONST_KARMAN !von Karman constant [-]
  real(r8) :: rwat   = SHR_CONST_RWV    !gas constant for water vapor [J/(kg K)]
  real(r8) :: rair   = SHR_CONST_RDAIR  !gas constant for dry air [J/kg/K]
  real(r8) :: roverg = SHR_CONST_RWV/SHR_CONST_G*1000._r8 !Rw/g constant = (8.3144/0.018)/(9.80616)*1000. mm/K
  real(r8) :: cpliq  = SHR_CONST_CPFW   !Specific heat of water [J/kg-K]
  real(r8) :: cpice  = SHR_CONST_CPICE  !Specific heat of ice [J/kg-K]
  real(r8) :: cpair  = SHR_CONST_CPDAIR !specific heat of dry air [J/kg/K]
  real(r8) :: hvap   = SHR_CONST_LATVAP !Latent heat of evap for water [J/kg]
  real(r8) :: hsub   = SHR_CONST_LATSUB !Latent heat of sublimation    [J/kg]
  real(r8) :: hfus   = SHR_CONST_LATICE !Latent heat of fusion for ice [J/kg]
  real(r8) :: denh2o = SHR_CONST_RHOFW  !density of liquid water [kg/m3]
  real(r8) :: denice = SHR_CONST_RHOICE !density of ice [kg/m3]
  real(r8) :: rgas   = SHR_CONST_RGAS   !universal gas constant [J/K/kmole]
  real(r8) :: tkair  = 0.023_r8     !thermal conductivity of air   [W/m/K]
  real(r8) :: tkice  = 2.290_r8     !thermal conductivity of ice   [W/m/K]
  real(r8) :: tkwat  = 0.57_r8       !thermal conductivity of water [W/m/K]
  real(r8) :: tfrz   = SHR_CONST_TKFRZ  !freezing temperature [K]
  real(r8), parameter :: tcrit  = 2.5_r8       !critical temperature to determine rain or snow
  real(r8) :: o2_molar_const = 0.209_r8   !constant atmospheric O2 molar ratio (mol/mol)

  real(r8) :: bdsno = 250._r8       !bulk density snow (kg/m**3)
  real(r8) :: alpha_aero = 1.0_r8   !constant for aerodynamic parameter weighting
  real(r8) :: tlsai_crit = 2.0_r8   !critical value of elai+esai for which aerodynamic parameters are maximum
  real(r8) :: watmin = 0.01_r8      !minimum soil moisture (mm)

  real(r8) :: re = SHR_CONST_REARTH*0.001_r8 !radius of earth (km)
  real(r8) :: oneatm = 1.01325e5_r8 !one standard atmospheric pressure
  real(r8), public, parameter :: degpsec = 15._r8/3600.0_r8 ! Degree's earth rotates per second

  real(r8), public, parameter ::  secspday= SHR_CONST_CDAY  ! Seconds per day
  integer,  public, parameter :: isecspday= secspday        ! Integer seconds per day
  real(r8), public, parameter ::  spval = 1.e36_r8  ! special value for real data
  integer , public, parameter :: ispval = -9999     ! special value for int data (keep this negative to avoid conflicts with possible valid values)

  ! These are tunable constants from clm2_3

  real(r8) :: zlnd = 0.01_r8      !Roughness length for soil [m]
  real(r8) :: zsno = 0.0024_r8    !Roughness length for snow [m]
  real(r8) :: csoilc = 0.004_r8   !Drag coefficient for soil under canopy [-]
  real(r8) :: capr   = 0.34_r8    !Tuning factor to turn first layer T into surface T
  real(r8) :: cnfac  = 0.5_r8     !Crank Nicholson factor between 0 and 1
  real(r8) :: ssi    = 0.033_r8   !Irreducible water saturation of snow
  real(r8) :: wimp   = 0.05_r8    !Water impremeable if porosity less than wimp
  real(r8) :: pondmx = 0.0_r8     !Ponding depth (mm)
  real(r8) :: pondmx_urban = 1.0_r8  !Ponding depth for urban roof and impervious road (mm)

  real(r8) :: thk_bedrock = 3.0_r8      ! thermal conductivity of 'typical' saturated granitic rock
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

  !  real(r8), parameter :: nitrif_n2o_loss_frac = 0.02_r8   !fraction of N lost as N2O in nitrification (Parton et al., 2001)
  real(r8), parameter :: nitrif_n2o_loss_frac = 6.e-4_r8   !fraction of N lost as N2O in nitrification (Li et al., 2000)
  real(r8), parameter :: frac_minrlztn_to_no3 = 0.2_r8   !fraction of N mineralized that is dieverted to the nitrification stream (Parton et al., 2001)

  !------------------------------------------------------------------
  ! Initialize landunit & column type constants
  !------------------------------------------------------------------

  integer, parameter :: istsoil    = 1  !soil         landunit type (natural vegetation)
  integer, parameter :: istcrop    = 2  !crop         landunit type
  integer, parameter :: istice     = 3  !land ice     landunit type (glacier)
  integer, parameter :: istice_mec = 4  !land ice (multiple elevation classes) landunit type
  integer, parameter :: istdlak    = 5  !deep lake    landunit type (now used for all lakes)
  integer, parameter :: istwet     = 6  !wetland      landunit type (swamp, marsh, etc.)

  integer, parameter :: isturb_MIN = 7  !minimum urban type index
  integer, parameter :: isturb_tbd = 7  !urban tbd    landunit type
  integer, parameter :: isturb_hd  = 8  !urban hd     landunit type
  integer, parameter :: isturb_md  = 9  !urban md     landunit type
  integer, parameter :: isturb_MAX = 9  !maximum urban type index

  integer, parameter :: max_lunit  = 9  !maximum value that lun%itype can have
                                        !(i.e., largest value in the above list)

  integer, parameter                   :: landunit_name_length = 12  ! max length of landunit names
  character(len=landunit_name_length)  :: landunit_names(max_lunit)  ! name of each landunit type


  real(r8), parameter :: catomw = 12.011_r8     ! molar mass of C atoms (g/mol)
  real(r8), parameter :: natomw = 14.007_r8     ! molar mass of N atoms (g/mol)
  real(r8), parameter :: patomw = 30.97_r8      ! molar mass of P atmos (g/mol)
  ! urban column types

  integer, parameter :: icol_roof        = 71
  integer, parameter :: icol_sunwall     = 72
  integer, parameter :: icol_shadewall   = 73
  integer, parameter :: icol_road_imperv = 74
  integer, parameter :: icol_road_perv   = 75

  ! parameters that depend on the above constants

  integer, parameter :: numurbl = isturb_MAX - isturb_MIN + 1   ! number of urban landunits


  character(len=16), parameter :: grlnd  = 'lndgrid'      ! name of lndgrid
  character(len=16), parameter :: nameg  = 'gridcell'     ! name of gridcells
  character(len=16), parameter :: namet  = 'topounit'     ! name of topographic units
  character(len=16), parameter :: namel  = 'landunit'     ! name of landunits
  character(len=16), parameter :: namec  = 'column'       ! name of columns
  character(len=16), parameter :: namep  = 'pft'          ! name of patches
  character(len=16), parameter :: nameCohort = 'cohort'   ! name of cohorts (ED specific)

  real(r8) :: d_con_g(ngases,2)    ! gas diffusivity constants (spp, #) (cm^2/s) (mult. by 10^-9)
  data (d_con_g(1,i),i=1,2) /0.1875_r8, 0.0013_r8/ ! CH4
  data (d_con_g(2,i),i=1,2) /0.1759_r8, 0.00117_r8/ ! O2
  data (d_con_g(3,i),i=1,2) /0.1325_r8, 0.0009_r8/ ! CO2

  real(r8) :: d_con_w(ngases,3)    ! water diffusivity constants (spp, #)  (mult. by 10^-4)
  data (d_con_w(1,i),i=1,3) /0.9798_r8, 0.02986_r8, 0.0004381_r8/ ! CH4
  data (d_con_w(2,i),i=1,3) /1.172_r8, 0.03443_r8, 0.0005048_r8/ ! O2
  data (d_con_w(3,i),i=1,3) /0.939_r8, 0.02671_r8, 0.0004095_r8/ ! CO2

  real(r8), allocatable :: zisoi(:)        !soil zi (interfaces)

contains


!------------------------------------------------------------------------------
  subroutine clm_varcon_init()
  use clm_varpar, only: nlevgrnd

  implicit none


    allocate( zisoi(0:nlevgrnd               ))
  end subroutine clm_varcon_init
end module clm_varcon
