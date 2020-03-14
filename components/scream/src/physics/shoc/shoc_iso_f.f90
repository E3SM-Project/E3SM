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

  function cxx_pow(base, exp) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base
    real(kind=c_real), value, intent(in)  :: exp

    ! return
    real(kind=c_real)               :: cxx_pow
  end function cxx_pow

  function cxx_sqrt(base) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base

    ! return
    real(kind=c_real)               :: cxx_sqrt
  end function cxx_sqrt

  function cxx_cbrt(base) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base

    ! return
    real(kind=c_real)               :: cxx_cbrt
  end function cxx_cbrt

  function cxx_gamma(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_gamma
  end function cxx_gamma

  function cxx_log(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_log
  end function cxx_log

  function cxx_log10(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_log10
  end function cxx_log10

  function cxx_exp(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_exp
  end function cxx_exp

end interface

end module shoc_iso_f
