module physics_share_f2c
  use iso_c_binding
  implicit none

#include "scream_config.f"
#ifdef SCREAM_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

!
! This file contains Fortran utilities common to all physics models in SCREAM.
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

  function cxx_expm1(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_expm1
  end function cxx_expm1
  
  function cxx_tanh(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: cxx_tanh
  end function cxx_tanh

  function cxx_erf(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: input

    ! return
    real(kind=c_real) :: cxx_erf
  end function cxx_erf

end interface

end module physics_share_f2c
