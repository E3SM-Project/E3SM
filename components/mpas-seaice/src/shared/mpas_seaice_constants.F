!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_constants
!
!> \brief
!> \author Adrian K. Turner and Elizabeth Hunke, LANL
!> \date 2013-2014, 2023
!> \details
!
!-----------------------------------------------------------------------

module seaice_constants

  use mpas_derived_types
#ifdef CCSMCOUPLED
     use shr_const_mod, only : SHR_CONST_G
#endif
  ! Map shared physical references from PCD here; SeaIce-specific settings and
  ! quantities not yet represented in PCD remain defined locally below.
  use pcd, only : pcd_pi => pi, &
       pcd_seconds_per_day => seconds_per_day, &
       pcd_gravity => standard_acceleration_of_gravity, &
       pcd_ice_density => sea_ice_density_reference, &
       pcd_seawater_density => seawater_density_reference, &
       pcd_freshwater_density => pure_water_density_reference, &
       pcd_stefan_boltzmann => stefan_boltzmann_constant, &
       pcd_freezing_temperature => pure_water_freezing_temperature_reference, &
       pcd_ice_specific_heat => sea_ice_specific_heat_capacity_reference, &
       pcd_ice_thermal_conductivity => sea_ice_thermal_conductivity_reference, &
       pcd_air_specific_heat => dry_air_specific_heat_capacity_reference, &
       pcd_vapor_specific_heat => water_vapor_specific_heat_capacity_reference, &
       pcd_rwv => water_vapor_specific_gas_constant_reference, &
       pcd_rair => dry_air_specific_gas_constant_reference, &
       pcd_latent_heat_sublimation => latent_heat_of_sublimation_reference, &
       pcd_seawater_specific_heat => seawater_specific_heat_capacity_reference, &
       pcd_latent_heat_vaporization => latent_heat_of_vaporization_reference, &
       pcd_latent_heat_fusion => latent_heat_of_fusion_reference, &
       pcd_ice_reference_salinity => sea_ice_reference_salinity, &
       pcd_angular_velocity => angular_velocity, &
       pcd_carbon_molar_mass => carbon_molar_mass_reference, &
       pcd_iron_molar_mass => iron_molar_mass_reference

  private
  save

  real (kind=RKIND), parameter, public :: &
  ! fundamental constants
       pi                              = pcd_pi           ,&! pi
       seaiceSecondsPerDay             = pcd_seconds_per_day ! fixed day unit (s), NIST SP 330

  real (kind=RKIND), public :: &
#ifdef CCSMCOUPLED
     seaiceGravity                   = SHR_CONST_G
#else
     seaiceGravity                   = pcd_gravity
#endif

  ! physical constants
  real (kind=RKIND), public :: &
       seaiceDensityIce                = pcd_ice_density  ,&! density of ice (kg/m^3)
       seaiceDensitySeaWater           = pcd_seawater_density,&! density of seawater (kg/m^3)
       seaiceDensityFreshwater         = pcd_freshwater_density,&! density of fresh water (kg/m^3)

  ! thermodynamic constants
       seaiceStefanBoltzmann           = pcd_stefan_boltzmann,&! W/m^2/K^4
       seaiceIceSnowEmissivity         = 1.0_RKIND        ,&! emissivity of snow and ice
       seaiceFreshWaterFreezingPoint   = pcd_freezing_temperature,&! freezing temp of fresh ice (K)
       seaiceFreshIceSpecificHeat      = pcd_ice_specific_heat,&! specific heat of fresh ice (J/kg/K)
       seaiceAirSpecificHeat           = pcd_air_specific_heat,&! specific heat of air (J/kg/K)
       seaiceWaterVaporSpecificHeat    = pcd_vapor_specific_heat,&! specific heat of water vapor (J/kg/K)
       seaiceZvir                      = (pcd_rwv / pcd_rair) - 1.0_RKIND,&! rh2o/rair - 1.0
       seaiceLatentHeatSublimation     = pcd_latent_heat_sublimation,&! latent heat, sublimation freshwater (J/kg)
       seaiceIceSurfaceMeltingTemperature= pcd_freezing_temperature-pcd_freezing_temperature,&! melting temp. ice top surface  (C)
       seaiceSnowSurfaceMeltingTemperature= pcd_freezing_temperature-pcd_freezing_temperature,&! melting temp. snow top surface (C)
       seaiceVonKarmanConstant         = 0.4_RKIND        ,&! von Karman constant (not yet in PCD)
       seaiceSeaWaterSpecificHeat      = pcd_seawater_specific_heat,&! specific heat of ocn (J/kg/K)
                                                            ! freshwater value needed for enthalpy
       seaiceLatentHeatVaporization    = pcd_latent_heat_vaporization,&! latent heat, vaporization freshwater (J/kg)
       seaiceLatentHeatMelting         = pcd_latent_heat_fusion,&! latent heat of melting of fresh ice (J/kg)
       seaiceReferenceSalinity         = pcd_ice_reference_salinity,&! ice reference salinity (ppt)
       seaiceSnowPatchiness            = 0.005_RKIND      ,&! parameter for fractional snow area (m)

       gramsIronPerMolIron             = pcd_iron_molar_mass * 1.0e3_RKIND ! g iron per mol iron

  ! fundamental constants
  real (kind=RKIND), parameter, public :: &
       pii                                = pi, &    ! replace pii with seaicePi elsewhere (might need to be parameter)
       seaiceDegreesToRadians             = pi / 180.0_RKIND, &
       seaiceRadiansToDegrees             = 180.0_RKIND / pi, &
       seaiceSecondsPerYear               = seaiceSecondsPerDay * 365.0_RKIND, & ! fixed 365-day year; not calendar-aware
       seaiceDaysPerSecond                = 1.0_RKIND/seaiceSecondsPerDay

  real (kind=RKIND), public :: &
       seaicePi                           = pi ! pi

  character (len=*), public, parameter :: &
       coupleAlarmID                      = 'coupling'

  real(kind=RKIND), public :: &
       seaicePuny                         = 1.0e-11_RKIND   ! a small value

  ! physical constants
  real (kind=RKIND), parameter, public :: &
       omega                              = pcd_angular_velocity ! angular rotation rate of Earth [s-1]

  real(kind=RKIND), public :: &
       seaiceDensitySnow                  = 330.0_RKIND     ! density of snow (kg/m^3)

  ! thermodynamic constants
  real(kind=RKIND), public :: &
       seaiceStabilityReferenceHeight     = 10._RKIND   , & ! reference height for stability (m)
       seaiceSnowSurfaceScatteringLayer   = 0.040_RKIND , & ! snow surface scattering layer thickness (m)
       seaiceMeltingTemperatureDepression = 0.054_RKIND , & ! melting temp. depression factor (C/ppt)
       seaiceOceanAlbedo                  = 0.06_RKIND  , & ! ocean albedo
       seaiceIceSurfaceRoughness          = 0.0005_RKIND, & ! ice surface roughness (m)
       seaiceMaximumSalinity              = 3.2_RKIND   , & ! max salinity at ice base for BL99 (ppt)
                                                            ! for mushy thermo:
       seaiceFrazilSalinityReduction      = 3.0_RKIND   , & ! bulk salinity reduction of newly formed frazil (ppt)
       seaiceFrazilIcePorosity            = 0.75_RKIND  , & ! initial liquid fraction of frazil 
       seaiceFrazilMinimumThickness       = 0.05_RKIND      ! min thickness of new frazil i

  ! dynamics constants
  real(kind=RKIND), public :: &
       seaiceIceStrengthConstantHiblerP   = 2.75e4_RKIND , &! P* constant in Hibler strength formula
       seaiceIceStrengthConstantHiblerC   = 20._RKIND       ! C* constant in Hibler strength formula

  ! minimum sea ice area
  real(kind=RKIND), public :: &
       iceAreaMinimum, &
       iceThicknessMinimum, &
       snowThicknessMinimum

  real(kind=RKIND), public :: &
       seaiceBigNumber                   = 1.0e+30_RKIND, & ! a large number
       seaiceSnowMinimumDensity          = 100.0_RKIND  , & ! minimum snow density (kg/m^3)
       seaiceBrineDynamicViscosity       = 1.79e-3_RKIND, & ! dynamic viscosity of brine (kg/m/s)
       seaiceFreezingTemperatureConstant = -1.8_RKIND   , & ! freezing temp of seawater (C), used 
                                                            ! as Tsfcn for open water only when 
                                                            ! tfrz_option is 'minus1p8' or null
       seaiceExtinctionCoef              = 1.4_RKIND    , & ! vis extnctn coef in ice, wvlngth<700nm (1/m)
       seaiceFreshIceConductivity        = pcd_ice_thermal_conductivity, & ! thermal conductivity of fresh ice(W/m/deg)
                                                            ! (kice) is not used for mushy thermo
       seaiceSnowMinimumThickness        = 1.e-4_RKIND  , & ! min snow thickness for computing zTsn (m)
       ! weights for albedos for history and diagnostics
       ! 4 Jan 2007 BPB  Following are appropriate for complete cloud
       ! in a summer polar atmosphere with 1.5m bare sea ice surface:
       ! .636/.364 vis/nir with only 0.5% direct for each band.
       seaiceAlbedoWtVisibleDirect       = 0.00318_RKIND, & ! visible, direct
       seaiceAlbedoWtNearIRDirect        = 0.00182_RKIND, & ! near IR, direct
       seaiceAlbedoWtVisibleDiffuse      = 0.63282_RKIND, & ! visible, diffuse
       seaiceAlbedoWtNearIRDiffuse       = 0.36218_RKIND, & ! near IR, diffuse
       ! constants for saturation humidity over ice and ocean
       seaiceQsatQiceConstant            = 11637800._RKIND   ,&
       seaiceQsatTiceConstant            = 5897.8_RKIND      ,&
       seaiceQsatQocnConstant            = 627572.4_RKIND    ,&
       seaiceQsatTocnConstant            = 5107.4_RKIND

! Are these needed for standalone runs?
!      ! orbital parameters
!      integer, public :: iyear_AD  ! Year to calculate orbit for
!      real(kind=RKIND),public :: eccen  ! Earth's orbital eccentricity
!      real(kind=RKIND),public :: obliqr ! Earth's obliquity in radians
!      real(kind=RKIND),public :: lambm0 ! Mean longitude of perihelion at the
!                                           ! vernal equinox (radians)
!      real(kind=RKIND),public :: mvelpp ! Earth's moving vernal equinox longitude
!                                           ! of perihelion + pi (radians)
!      real(kind=RKIND),public :: obliq  ! obliquity in degrees
!      real(kind=RKIND),public :: mvelp  ! moving vernal equinox long
!      real(kind=RKIND),public :: decln  ! solar declination angle in radians
!      real(kind=RKIND),public :: eccf   ! earth orbit eccentricity factor
!      logical,public :: log_print ! Flags print of status/error

  ! biogeochemistry constants
  real(kind=RKIND), public :: &
       skeletalLayerThickness             = 0.03_RKIND    , & ! (m) skeletal layer thickness
       gramsCarbonPerMolCarbon            = pcd_carbon_molar_mass * 1.0e3_RKIND, & ! g carbon per mol carbon
       microgramsPerKilograms             = 1.0e9_RKIND       ! kg to ug conversion

  ! ocean biogeochemistry ISPOL values
  real(kind=RKIND), parameter, public :: &
       oceanAmmoniumISPOL        = 1.0_RKIND   , & ! mmol N m-3
       oceanDMSISPOL             = 0.1_RKIND   , & ! mmol S m-3
       oceanDMSPISPOL            = 0.1_RKIND   , & ! mmol S m-3
       oceanDiatomsISPOL         = 1.0_RKIND   , & ! mmol N m-3
       oceanSmallAlgaeISPOL      = 0.0057_RKIND, & ! mmol N m-3
       oceanPhaeocystisISPOL     = 0.0027_RKIND, & ! mmol N m-3
       oceanPolysaccharidsISPOL  = 16.2_RKIND  , & ! mmol C m-3
       oceanLipidsISPOL          = 9.0_RKIND   , & ! mmol C m-3
       oceanProteinsCarbonISPOL  = 9.0_RKIND   , & ! mmol C m-3
       oceanDICISPOL             = 1.0_RKIND   , & ! mmol C m-3
       oceanProteinsISPOL        = 12.9_RKIND  , & ! mmol N m-3
       oceanDissolvedIronISPOL   = 0.4_RKIND   , & ! mmol Fe m-3
       oceanParticulateIronISPOL = 2.0_RKIND   , & ! mmol Fe m-3
       oceanHumicsISPOL          = 1.0_RKIND       ! mmol C m-3

end module seaice_constants
