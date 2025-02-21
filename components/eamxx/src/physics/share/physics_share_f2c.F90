module physics_share_f2c
  use iso_c_binding
  implicit none

#include "eamxx_config.f"
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

  function scream_pow(base, exp) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base
    real(kind=c_real), value, intent(in)  :: exp

    ! return
    real(kind=c_real)               :: scream_pow
  end function scream_pow

  function scream_sqrt(base) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base

    ! return
    real(kind=c_real)               :: scream_sqrt
  end function scream_sqrt

  function scream_cbrt(base) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: base

    ! return
    real(kind=c_real)               :: scream_cbrt
  end function scream_cbrt

  function scream_gamma(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: scream_gamma
  end function scream_gamma

  function scream_log(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: scream_log
  end function scream_log

  function scream_log10(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: scream_log10
  end function scream_log10

  function scream_exp(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: scream_exp
  end function scream_exp

  function scream_expm1(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: scream_expm1
  end function scream_expm1

  function scream_tanh(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in) :: input

    ! return
    real(kind=c_real)            :: scream_tanh
  end function scream_tanh

  function scream_erf(input) bind(C)
    use iso_c_binding

    !arguments:
    real(kind=c_real), value, intent(in)  :: input

    ! return
    real(kind=c_real) :: scream_erf
  end function scream_erf

end interface

end module physics_share_f2c
