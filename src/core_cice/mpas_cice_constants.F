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

  use mpas_grid_types

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

  real(kind=RKIND), parameter :: &
       puny = 1.0e-11_RKIND

  ! minimum sea ice area
  real(kind=RKIND), parameter, public :: &
       iceAreaMinimum       = puny, &
       iceThicknessMinimum  = puny, &
       snowThicknessMinimum = puny

end module cice_constants

