!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_constants
!
!> \brief
!> \author Adrian K. Turner and Elizabeth Hunke, LANL
!> \date 2013-2014, 2023
!> \details
!> Values in this module have been copied from the two ice_constants_colpkg.F90 
!> modules in the column physics package (in column/).  They will remain 
!> duplicates until the column physics package is deprecated.
!
!-----------------------------------------------------------------------

module seaice_constants

  use mpas_derived_types
#ifdef CCSMCOUPLED
  use shr_const_mod
#endif

  private
  save

#ifdef CCSMCOUPLED
      real (kind=RKIND), parameter, public :: &
  ! fundamental constants
       pi                              = SHR_CONST_PI     ,&! pi
       seaiceSecondsPerDay             = SHR_CONST_CDAY     ! seconds in calendar day

      real (kind=RKIND), public :: &
       seaiceGravity                   = SHR_CONST_G      ,&! gravitational acceleration (m/s^2)

  ! physical constants
       seaiceDensityIce                = SHR_CONST_RHOICE ,&! density of ice (kg/m^3)
       seaiceDensitySeaWater           = SHR_CONST_RHOSW  ,&! density of seawater (kg/m^3)
       seaiceDensityFreshwater         = SHR_CONST_RHOFW  ,&! density of fresh water (kg/m^3)

  ! thermodynamic constants
       seaiceStefanBoltzmann           = SHR_CONST_STEBOL ,&!  W/m^2/K^4
       seaiceIceSnowEmissivity         = 1.0_RKIND        ,&! emissivity of snow and ice
       seaiceFreshWaterFreezingPoint   = SHR_CONST_TKFRZ  ,&! freezing temp of fresh ice (K)
       seaiceFreshIceSpecificHeat      = SHR_CONST_CPICE  ,&! specific heat of fresh ice (J/kg/K)
       seaiceAirSpecificHeat           = SHR_CONST_CPDAIR ,&! specific heat of air (J/kg/K)
       seaiceWaterVaporSpecificHeat    = SHR_CONST_CPWV   ,&! specific heat of water vapor (J/kg/K)
       seaiceZvir                      = SHR_CONST_ZVIR   ,&! rh2o/rair - 1.0
       seaiceLatentHeatSublimation     = SHR_CONST_LATSUB ,&! latent heat, sublimation freshwater (J/kg)
       seaiceIceSurfaceMeltingTemperature= SHR_CONST_TKFRZ-SHR_CONST_TKFRZ ,&! melting temp. ice top surface  (C)
       seaiceSnowSurfaceMeltingTemperature= SHR_CONST_TKFRZ-SHR_CONST_TKFRZ ,&! melting temp. snow top surface (C)
       seaiceVonKarmanConstant         = SHR_CONST_KARMAN ,&! von Karman constant
       seaiceSeaWaterSpecificHeat      = SHR_CONST_CPSW   ,&! specific heat of ocn    (J/kg/K)
                                                            ! freshwater value needed for enthalpy
       seaiceLatentHeatVaporization    = SHR_CONST_LATVAP ,&! latent heat, vaporization freshwater (J/kg)
       seaiceLatentHeatMelting         = SHR_CONST_LATICE ,&! latent heat of melting of fresh ice (J/kg)
       seaiceReferenceSalinity         = SHR_CONST_ICE_REF_SAL,&! ice reference salinity (ppt)
       seaiceSnowPatchiness            = 0.005_RKIND        ! parameter for fractional snow area (m)

!         R_gC2molC = SHR_CONST_MWC    ,&! molar mass of carbon

#else
      real (kind=RKIND), parameter, public :: &
  ! fundamental constants
       pi                              = 3.14159265358979323846_RKIND, &
       seaiceSecondsPerDay             = 24.0_RKIND * 3600.0_RKIND
 
      real (kind=RKIND), public :: &
       seaiceGravity                   = 9.80616_RKIND ,&! gravitational acceleration (m/s^2)

  ! physical constants
       seaiceDensityIce                = 917.0_RKIND   ,&! density of ice (kg/m^3)
       seaiceDensitySeaWater           = 1026.0_RKIND  ,&! density of seawater (kg/m^3)
       seaiceDensityFreshwater         = 1000.0_RKIND  ,&! density of fresh water (kg/m^3)

  ! thermodynamic constants
       seaiceStefanBoltzmann           = 567.0e-10_RKIND,&! W/m^2/K^4
       seaiceIceSnowEmissivity         = 0.95_RKIND    ,&! emissivity of snow and ice ! 0.985_RKIND in Icepack
       seaiceFreshWaterFreezingPoint   = 273.15_RKIND  ,&! freezing temp of fresh ice (K)
       seaiceFreshIceSpecificHeat      = 2106._RKIND   ,&! specific heat of fresh ice (J/kg/K)
       seaiceAirSpecificHeat           = 1005.0_RKIND  ,&! specific heat of air (J/kg/K)
       seaiceWaterVaporSpecificHeat    = 1.81e3_RKIND  ,&! specific heat of water vapor (J/kg/K)
       seaiceZvir                      = 0.606_RKIND   ,&! rh2o/rair - 1.0
       seaiceLatentHeatSublimation     = 2.835e6_RKIND ,&! latent heat, sublimation freshwater (J/kg)
       seaiceIceSurfaceMeltingTemperature  = 0.0_RKIND ,&! melting temperature, ice top surface  (C)
       seaiceSnowSurfaceMeltingTemperature = 0.0_RKIND ,&! melting temperature, snow top surface  (C)
       seaiceVonKarmanConstant         = 0.4_RKIND     ,&! von Karman constant
       seaiceSeaWaterSpecificHeat      = 4218._RKIND   ,&! specific heat of ocn    (J/kg/K)
                                                         ! freshwater value needed for enthalpy
       seaiceLatentHeatVaporization    = 2.501e6_RKIND ,&! latent heat, vaporization freshwater (J/kg)
!       seaiceLatentHeatMelting         = seaiceLatentHeatSublimation & ! latent heat of melting of fresh ice (J/kg)
!                                       - seaiceLatentHeatVaporization, &
       seaiceLatentHeatMelting         = 3.34e5_RKIND  ,&! latent heat of melting of fresh ice (J/kg)
       seaiceReferenceSalinity         = 4._RKIND      ,&! ice reference salinity (ppt)
       seaiceSnowPatchiness            = 0.02_RKIND      ! parameter for fractional snow area (m)

!         R_gC2molC = 12.0107_RKIND, & ! molar mass of carbon
#endif

  ! fundamental constants
  real (kind=RKIND), parameter, public :: &
       pii                                = pi, &    !echmod - replace pii with seaicePi elsewhere
       seaiceDegreesToRadians             = pi / 180.0_RKIND, &
       seaiceRadiansToDegrees             = 180.0_RKIND / pi, &
       seaiceSecondsPerYear               = 24.0_RKIND * 3600.0_RKIND * 365.0_RKIND, & !echmod - incorrect for leap years
       seaiceDaysPerSecond                = 1.0_RKIND/seaiceSecondsPerDay

  real (kind=RKIND), public :: &
       seaicePi                           = pi ! pi

  character (len=*), public, parameter :: &
       coupleAlarmID                      = 'coupling'

  real(kind=RKIND), public :: &
       seaicePuny                         = 1.0e-11_RKIND   ! a small value

  ! physical constants
  real (kind=RKIND), parameter, public :: &
       omega                              = 7.29212e-5_RKIND ! angular rotation rate of the Earth [s-1]

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

!echmod - declare and define these constants explicitly in case Icepack changes its defaults
  real(kind=RKIND), public :: &
       seaiceBigNumber                   = 1.0e+30_RKIND, & ! a large number
       seaiceSnowMinimumDensity          = 100.0_RKIND  , & ! minimum snow density (kg/m^3)
       seaiceBrineDynamicViscosity       = 1.79e-3_RKIND, & ! dynamic viscosity of brine (kg/m/s)
       seaiceFreezingTemperatureConstant = -1.8_RKIND   , & ! freezing temp of seawater (C), used 
                                                            ! as Tsfcn for open water only when 
                                                            ! tfrz_option is 'minus1p8' or null
       seaiceExtinctionCoef              = 1.4_RKIND    , & ! vis extnctn coef in ice, wvlngth<700nm (1/m)
       seaiceFreshIceConductivity        = 2.03_RKIND   , & ! thermal conductivity of fresh ice(W/m/deg)
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

!echmod - Are these needed for standalone runs?
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
       skeletalLayerThickness             = 0.03_RKIND    ,&! (m) skeletal layer thickness
       gramsCarbonPerMolCarbon            = 12.0107_RKIND   ! g carbon per mol carbon

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
