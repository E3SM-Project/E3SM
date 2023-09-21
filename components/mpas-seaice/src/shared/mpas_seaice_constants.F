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
       pii = 3.141592653589793_RKIND, &          ! CESM uses SHR_CONST_PI
       seaiceDegreesToRadians = pii / 180.0_RKIND, &
       seaiceRadiansToDegrees = 180.0_RKIND / pii, &
       seaiceSecondsPerYear = 24.0_RKIND * 3600.0_RKIND * 365.0_RKIND, &
       seaiceSecondsPerDay  = 24.0_RKIND * 3600.0_RKIND, &
       seaiceDaysPerSecond = 1.0_RKIND/seaiceSecondsPerDay

  real (kind=RKIND), public :: &
       seaicePi                           ! pi

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
       seaiceFreshIceSpecificHeat, &    ! specific heat of fresh ice (J/kg/K)
       seaiceAirSpecificHeat, &         ! specific heat of air (J/kg/K)
       seaiceWaterVaporSpecificHeat, &  ! specific heat of water vapor (J/kg/K)
       seaiceZvir, &                    ! rh2o/rair - 1.0
       seaiceLatentHeatSublimation, &   ! latent heat, sublimation freshwater (J/kg)
       seaiceLatentHeatMelting, &       ! latent heat of melting of fresh ice (J/kg)
       seaiceIceSurfaceMeltingTemperature, &  ! melting temp. ice top surface  (C)
       seaiceSnowSurfaceMeltingTemperature, & ! melting temp. snow top surface  (C)
       seaiceMeltingTemperatureDepression, &  ! melting temp. depression factor (C/ppt)
       seaiceOceanAlbedo, &             ! Ocean albedo
       seaiceVonKarmanConstant, &       ! Von Karman constant
       seaiceIceSurfaceRoughness, &     ! ice surface roughness (m)
       seaiceSeaWaterSpecificHeat, &    ! specific heat of ocn (J/kg/K)
       seaiceLatentHeatVaporization, &  ! latent heat, vaporization freshwater (J/kg)
       seaiceReferenceSalinity, &       ! ice reference salinity (ppt)
       seaiceMaximumSalinity, &         ! ice maximum salinity (ppt)
       seaiceStabilityReferenceHeight, &! stability reference height (m)
       seaiceSnowPatchiness             ! snow patchiness parameter

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

   ! ocean biogeochemistry ISPOL values
   real(kind=RKIND), parameter, public :: &
         oceanAmmoniumISPOL        = 1.0_RKIND, & ! mmol N m-3
         oceanDMSISPOL             = 0.1_RKIND, & ! mmol S m-3
         oceanDMSPISPOL            = 0.1_RKIND, & ! mmol S m-3
         oceanDiatomsISPOL         = 1.0_RKIND, & ! mmol N m-3
         oceanSmallAlgaeISPOL      = 0.0057_RKIND, & ! mmol N m-3
         oceanPhaeocystisISPOL     = 0.0027_RKIND, & ! mmol N m-3
         oceanPolysaccharidsISPOL  = 16.2_RKIND, & ! mmol C m-3
         oceanLipidsISPOL          = 9.0_RKIND, & ! mmol C m-3
         oceanProteinsCarbonISPOL  = 9.0_RKIND, & ! mmol C m-3
         oceanDICISPOL             = 1.0_RKIND, & ! mmol C m-3
         oceanProteinsISPOL        = 12.9_RKIND, & ! mmol N m-3
         oceanDissolvedIronISPOL   = 0.4_RKIND, & ! mmol Fe m-3
         oceanParticulateIronISPOL = 2.0_RKIND,& ! mmol Fe m-3
         oceanHumicsISPOL          = 1.0_RKIND ! mmol C m-3

end module seaice_constants
