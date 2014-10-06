module clm_varcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varcon
!
! !DESCRIPTION:
! Module containing various model constants
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
  use clm_varpar   , only: numrad, nlevgrnd, nlevlak
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 27 February 2008: Keith Oleson; Add forcing height and aerodynamic parameters
!
!EOP
!-----------------------------------------------------------------------

  !------------------------------------------------------------------
  ! Initialize mathmatical constants
  !------------------------------------------------------------------

  real(r8) :: rpi    = SHR_CONST_PI

  !------------------------------------------------------------------
  ! Initialize physical constants
  !------------------------------------------------------------------

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
  real(r8) :: tkwat  = 0.6_r8       !thermal conductivity of water [W/m/K]
  real(r8) :: tfrz   = SHR_CONST_TKFRZ  !freezing temperature [K]
  real(r8) :: tcrit  = 2.5_r8       !critical temperature to determine rain or snow
  real(r8) :: o2_molar_const = 0.209_r8   !constant atmospheric O2 molar ratio (mol/mol)

  real(r8) :: bdsno = 250._r8       !bulk density snow (kg/m**3)
  real(r8) :: alpha_aero = 1.0_r8   !constant for aerodynamic parameter weighting
  real(r8) :: tlsai_crit = 2.0_r8   !critical value of elai+esai for which aerodynamic parameters are maximum
  real(r8) :: watmin = 0.01_r8      !minimum soil moisture (mm)

  real(r8) :: re = SHR_CONST_REARTH*0.001_r8 !radius of earth (km)

  real(r8), public, parameter :: degpsec = 15._r8/3600.0_r8 ! Degree's earth rotates per second

  real(r8), public, parameter ::  secspday= SHR_CONST_CDAY  ! Seconds per day
  integer,  public, parameter :: isecspday= secspday        ! Integer seconds per day
  real(r8), public, parameter ::  spval = 1.e36_r8  ! special value for real data
  integer , public, parameter :: ispval = -9999     ! special value for int data

  ! These are tunable constants from clm2_3

  real(r8) :: zlnd = 0.01_r8      !Roughness length for soil [m]
  real(r8) :: zsno = 0.0024_r8    !Roughness length for snow [m]
  real(r8) :: csoilc = 0.004_r8   !Drag coefficient for soil under canopy [-]
  real(r8) :: capr   = 0.34_r8    !Tuning factor to turn first layer T into surface T
  real(r8) :: cnfac  = 0.5_r8     !Crank Nicholson factor between 0 and 1
  real(r8) :: ssi    = 0.033_r8   !Irreducible water saturation of snow
  real(r8) :: wimp   = 0.05_r8    !Water impremeable if porosity less than wimp
  real(r8) :: pondmx = 10.0_r8    !Ponding depth (mm)
  real(r8) :: pondmx_urban = 1.0_r8  !Ponding depth for urban roof and impervious road (mm)
  ! 4/14/05: PET
  ! Adding isotope code
  real(r8), parameter :: preind_atm_del13c = -6.0   ! preindustrial value for atmospheric del13C
  real(r8), parameter :: preind_atm_ratio = SHR_CONST_PDB + (preind_atm_del13c * SHR_CONST_PDB)/1000.0  ! 13C/12C
  real(r8) :: c13ratio = preind_atm_ratio/(1.0+preind_atm_ratio) ! 13C/(12+13)C preind atmosphere

  real(r8), parameter :: ht_efficiency_factor = 0.75_r8 !efficiency factor for urban heating (-)
  real(r8), parameter :: ac_efficiency_factor = 0.25_r8 !efficiency factor for urban air conditioning (-)
  real(r8) :: ht_wasteheat_factor = 1.0_r8/ht_efficiency_factor  !wasteheat factor for urban heating (-)
  real(r8) :: ac_wasteheat_factor = 1.0_r8/ac_efficiency_factor  !wasteheat factor for urban air conditioning (-)
  real(r8) :: wasteheat_limit = 100._r8  !limit on wasteheat (W/m2)
  
  real(r8), parameter :: h2osno_max = 1000._r8    ! max allowed snow thickness (mm H2O)
  real(r8), parameter :: lapse_glcmec = 0.006_r8  ! surface temperature lapse rate (deg m-1)
                                                  ! Pritchard et al. (GRL, 35, 2008) use 0.006  

  !------------------------------------------------------------------
  ! Initialize water type constants
  !------------------------------------------------------------------

  ! "land unit " types
  !   1     soil (includes vegetated landunits)
  !   2     land ice (glacier)
  !   3     deep lake
  !   4     shallow lake
  !   5     wetland (swamp, marsh, etc.)
  !   6     urban
  !   7     land ice (glacier) with multiple elevation classes
  !   8     crop

  integer :: istsoil    = 1  !soil         landunit type
  integer :: istice     = 2  !land ice     landunit type
  integer :: istdlak    = 3  !deep lake    landunit type
  integer :: istslak    = 4  !shallow lake landunit type
  integer :: istwet     = 5  !wetland      landunit type
  integer :: isturb     = 6  !urban        landunit type
  integer :: istice_mec = 7  !land ice (multiple elevation classes) landunit type
  integer :: istcrop    = 8  !crop         landunit type
  integer :: max_lunit  = 8  !maximum value that lun%itype can have
                             !(i.e., largest value in the above list)

  ! urban column types

  integer :: icol_roof        = 61
  integer :: icol_sunwall     = 62
  integer :: icol_shadewall   = 63
  integer :: icol_road_imperv = 64
  integer :: icol_road_perv   = 65

  !------------------------------------------------------------------
  ! Initialize miscellaneous radiation constants
  !------------------------------------------------------------------

  integer, private :: i  ! loop index

  real(r8), allocatable :: albsat(:,:) ! wet soil albedo by color class and waveband (1=vis,2=nir)
  real(r8), allocatable :: albdry(:,:) ! dry soil albedo by color class and waveband (1=vis,2=nir)

  real(r8) :: alblak(numrad)           ! albedo frozen lakes by waveband (1=vis, 2=nir)
  data (alblak(i),i=1,numrad) /0.60_r8, 0.40_r8/

  real(r8) :: betads  = 0.5_r8            ! two-stream parameter betad for snow
  real(r8) :: betais  = 0.5_r8            ! two-stream parameter betai for snow
  real(r8) :: omegas(numrad)           ! two-stream parameter omega for snow by band
  data (omegas(i),i=1,numrad) /0.8_r8, 0.4_r8/

  !------------------------------------------------------------------
  ! Soil and Lake depths are constants for now
  ! The values for the following arrays are set in routine iniTimeConst
  !------------------------------------------------------------------

  real(r8) :: zlak(1:nlevlak)     !lake z  (layers)
  real(r8) :: dzlak(1:nlevlak)    !lake dz (thickness)
  real(r8) :: zsoi(1:nlevgrnd)    !soil z  (layers)
  real(r8) :: dzsoi(1:nlevgrnd)   !soil dz (thickness)
  real(r8) :: zisoi(0:nlevgrnd)   !soil zi (interfaces)

end module clm_varcon
