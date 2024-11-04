! These macros provide support for math functions that give bit-for-bit
! identical results when called via Fortran or C++ (essentially by calling
! C++ versions of the functions from Fortran) if the macro
! SCREAM_CONFIG_IS_CMAKE is defined.
!
! Make sure to place the following lines at the top of any modules
! that use these macros:
!
! use physics_share_f2c, only: scream_pow, scream_sqrt, scream_cbrt, scream_gamma, scream_log, &
!                              scream_log10, scream_exp, scream_tanh, scream_erf

#ifndef SCREAM_BFB_MATH_INC
#define SCREAM_BFB_MATH_INC

#define bfb_square(val) ((val)*(val))
#define bfb_cube(val)   ((val)*(val)*(val))
#define bfb_quad(val)   (bfb_square(bfb_square(val)))

! This conditional must match CPP logic for SCREAM_BFB_TESTING in scream_types.hpp
#if !defined(SCREAM_CONFIG_IS_CMAKE) || defined (NDEBUG) || defined (SCREAM_SHORT_TESTS)
#  define bfb_pow(base, exp) (base)**(exp)
#  define bfb_cbrt(base) (base)**(1.0D0/3.0D0)
#  define bfb_gamma(val) gamma(val)
#  define bfb_log(val) log(val)
#  define bfb_log10(val) log10(val)
#  define bfb_exp(val) exp(val)
#  define bfb_expm1(val) (exp(val) - 1)
#  define bfb_tanh(val) tanh(val)
#  define bfb_sqrt(val) sqrt(val)
#  define bfb_tanh(val) tanh(val)
#  define bfb_erf(val) erf(val)
#else
#  define bfb_pow(base, exp) scream_pow(base, exp)
#  define bfb_sqrt(base) scream_sqrt(base)
#  define bfb_cbrt(base) scream_cbrt(base)
#  define bfb_gamma(val) scream_gamma(val)
#  define bfb_log(val) scream_log(val)
#  define bfb_log10(val) scream_log10(val)
#  define bfb_exp(val) scream_exp(val)
#  define bfb_expm1(val) scream_expm1(val)
#  define bfb_tanh(val) scream_tanh(val)
#  define bfb_erf(val) scream_erf(val)
#endif

#endif
