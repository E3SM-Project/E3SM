!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_constants
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_constants

  use mpas_derived_types

  private
  save

  ! fundamental constants
  real (kind=RKIND), parameter, public :: &
       pii = 3.141592653589793_RKIND, &
       seaiceDegreesToRadians = pii / 180.0_RKIND, &
       seaiceRadiansToDegrees = 180.0_RKIND / pii, &
       seaiceSecondsPerYear = 24.0_RKIND * 3600.0_RKIND * 365.0_RKIND, &
       seaiceSecondsPerDay  = 24.0_RKIND * 3600.0_RKIND, &
       seaiceDaysPerSecond = 1.0_RKIND/seaiceSecondsPerDay

  ! Earth constants
  real (kind=RKIND), public :: &
       seaiceGravity                      ! gravitational acceleration (m/s^2)
  real (kind=RKIND), parameter, public :: &
       omega         = 7.29212e-5_RKIND   ! angular rotation rate of the Earth [s-1]

  character (len=*), public, parameter :: &
       coupleAlarmID = 'coupling'

  real(kind=RKIND), public :: &
       seaicePuny

  ! physical constants
  real(kind=RKIND), public :: &
       seaiceDensityIce, &      ! density of ice (kg/m^3)
       seaiceDensitySnow, &     ! density of snow (kg/m^3)
       seaiceDensitySeaWater, & ! density of seawater (kg/m^3)
       seaiceDensityFreshwater  ! density of freshwater (kg/m^3)

  ! thermodynamic constants
  real(kind=RKIND), public :: &
       seaiceStefanBoltzmann, &         ! J m-2 K-4 s-1
       seaiceIceSnowEmissivity, &       ! emissivity of snow and ice
       seaiceFreshWaterFreezingPoint, & ! freezing temp of fresh ice (K)
       seaiceAirSpecificHeat, &         ! specific heat of air (J/kg/K)
       seaiceLatentHeatSublimation, &   ! latent heat, sublimation freshwater (J/kg)
       seaiceLatentHeatMelting, &       ! latent heat of melting of fresh ice (J/kg)
       seaiceOceanAlbedo, &             ! Ocean albedo
       seaiceVonKarmanConstant, &       ! Von Karman constant
       seaiceIceSurfaceRoughness, &     ! ice surface roughness (m)
       seaiceStabilityReferenceHeight   ! stability reference height (m)

  ! dynamics constants
  real(kind=RKIND), public :: &
       seaiceIceStrengthConstantHiblerP, & ! P* constant in Hibler strength formulation
       seaiceIceStrengthConstantHiblerC, & ! C* constant in Hibler strength formulation
       seaiceIceOceanDragCoefficient       ! ice ocean drag coefficient

  ! minimum sea ice area
  real(kind=RKIND), public :: &
       iceAreaMinimum, &
       iceThicknessMinimum, &
       snowThicknessMinimum

   ! biogeochemistry constants
   real(kind=RKIND), public :: &
        skeletalLayerThickness, &
        gramsCarbonPerMolCarbon   ! g carbon per mol carbon

end module seaice_constants
