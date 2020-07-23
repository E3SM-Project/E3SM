module shoc_iso_f
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains bridges from shoc fortran to scream c++.
!

interface

  !
  ! These are some routine math operations that are not BFB between
  ! fortran and C++ on all platforms, so fortran will need to use
  ! the C++ versions in order to stay BFB.
  !

end interface

end module shoc_iso_f
