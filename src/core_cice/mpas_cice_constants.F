module cice_constants

  use mpas_grid_types

  private
  
  real(kind=RKIND), parameter, public :: &
       pi = 3.14159265359_RKIND, &
       deg_to_rad = pi / 180.0_RKIND, &
       rad_to_deg = 180.0_RKIND / pi

  real(kind=RKIND), parameter, public :: &
       rhoi = 917.0_RKIND, &
       rhos = 330.0_RKIND


end module cice_constants

