!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  cice_constants
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module cice_constants

  use mpas_derived_types

  use ice_constants_colpkg, only: &
       gravit, &
       rhoi, &
       rhos, &
       rhow, &
       puny, &
       stefan_boltzmann, &
       emissivity, &
       Tffresh, &
       cp_air, &
       Lsub, &
       Lfresh, &
       Pstar, &
       Cstar, &
       dragio, &
       albocn, &
       rhofresh, &
       vonkar, &
       iceruf, &
       zref

  private
  save

  ! fundamental constants
  real (kind=RKIND), parameter, public :: &
       pii = 3.141592653589793_RKIND, &
       ciceDegreesToRadians = pii / 180.0_RKIND, &
       ciceRadiansToDegrees = 180.0_RKIND / pii, &
       ciceSecondsPerYear = 24.0_RKIND * 3600.0_RKIND * 365.0_RKIND, &
       ciceSecondsPerDay  = 24.0_RKIND * 3600.0_RKIND, &
       ciceDaysPerSecond = 1.0_RKIND/ciceSecondsPerDay

  ! Earth constants
  real (kind=RKIND), parameter, public :: &
       ciceGravity = gravit, &    ! gravitational acceleration (m/s^2)
       earthRadius = 6371229.0_RKIND, & ! radius of Earth in [m]
       omega       = 7.29212e-5_RKIND   ! angular rotation rate of the Earth [s-1]

  character (len=*), public, parameter :: &
       coupleAlarmID = 'coupling'

  real(kind=RKIND), parameter, public :: &
       cicePuny = puny

  ! physical constants
  real(kind=RKIND), parameter, public :: &
       ciceDensityIce        = rhoi, & ! density of ice (kg/m^3)
       ciceDensitySnow       = rhos, & ! density of snow (kg/m^3)
       ciceDensitySeaWater   = rhow, & ! density of seawater (kg/m^3)
       ciceDensityFreshwater = rhofresh ! density of freshwater (kg/m^3)

  ! thermodynamic constants
  real(kind=RKIND), parameter, public :: &
       ciceStefanBoltzmann          = stefan_boltzmann, & ! J m-2 K-4 s-1
       ciceIceSnowEmissivity        = emissivity, &       ! emissivity of snow and ice
       ciceFreshWaterFreezingPoint  = Tffresh, &          ! freezing temp of fresh ice (K)
       ciceAirSpecificHeat          = cp_air, &           ! specific heat of air (J/kg/K)
       ciceLatentHeatSublimation    = Lsub, &             ! latent heat, sublimation freshwater (J/kg)
       ciceLatentHeatMelting        = Lfresh, &           ! latent heat of melting of fresh ice (J/kg)
       ciceOceanAlbedo              = albocn, &           ! Ocean albedo
       ciceVonKarmanConstant        = vonkar, &           ! Von Karman constant
       ciceIceSurfaceRoughness      = iceruf, &           ! ice surface roughness (m)
       ciceStabilityReferenceHeight = zref ! stability reference height (m)

  ! dynamics constants
  real(kind=RKIND), parameter, public :: &
       ciceIceStrengthConstantHiblerP = Pstar, & ! P* constant in Hibler strength formulation
       ciceIceStrengthConstantHiblerC = Cstar, & ! C* constant in Hibler strength formulation
       ciceIceOceanDragCoefficient    = dragio   ! ice ocean drag coefficient

  ! minimum sea ice area
  real(kind=RKIND), parameter, public :: &
       iceAreaMinimum       = cicePuny, &
       iceThicknessMinimum  = cicePuny, &
       snowThicknessMinimum = cicePuny


end module cice_constants
