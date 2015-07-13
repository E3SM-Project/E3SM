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

  use mpas_constants, only: &
       pii

  private
  save
  
  real(kind=RKIND), parameter, public :: &
       ciceDegreesToRadians = pii / 180.0_RKIND, &
       ciceRadiansToDegrees = 180.0_RKIND / pii

  real(kind=RKIND), parameter, public :: &
       ciceDensityIce  = 917.0_RKIND, &
       ciceDensitySnow = 330.0_RKIND

  real(kind=RKIND), parameter, public :: &
       cicePuny = 1.0e-11_RKIND

  ! minimum sea ice area
  real(kind=RKIND), parameter, public :: &
       iceAreaMinimum       = cicePuny, &
       iceThicknessMinimum  = cicePuny, &
       snowThicknessMinimum = cicePuny


  real(kind=RKIND), parameter, public :: &
       ciceStefanBoltzmann = 5.67e-8_RKIND , & ! J m-2 K-4 s-1
       ciceIceSnowEmissivity = 0.95_RKIND, & ! emissivity of snow and ice
       ciceFreshWaterFreezingPoint = 273.15_RKIND, & ! freezing temp of fresh ice (K)
       ciceAirSpecificHeat = 1005.0_RKIND, & ! specific heat of air (J/kg/K)
       ciceLatentHeatSublimation = 2.835e6_RKIND ! latent heat, sublimation freshwater (J/kg)



 real(kind=RKIND), parameter, public :: &
      ciceSecondsPerYear = 365.0_RKIND * 24.0_RKIND * 3600.0_RKIND, &
      ciceSecondsPerDay = 24.0_RKIND * 3600.0_RKIND

end module cice_constants

