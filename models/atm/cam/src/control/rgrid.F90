
module rgrid

  use pmgrid, only: plat, plon
  use pspect, only: pmmax, pmax, ptrm

  implicit none

  integer :: nlon(plat)        = plon ! num longitudes per latitude
  integer :: beglatpair(pmmax) = huge(1)
  integer :: nmmax(plat/2)     = huge(1)
  integer :: wnummax(plat)     = ptrm ! cutoff Fourier wavenumber

  logical :: fullgrid                   ! true => no grid reduction towards poles
end module rgrid
