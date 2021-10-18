! These macros provide support for math functions that give bit-for-bit
! identical results when called via Fortran or C++ (essentially by calling
! C++ versions of the functions from Fortran) if the macro
! SCREAM_CONFIG_IS_CMAKE is defined.
!
! Make sure to place the following lines at the top of any modules
! that use these macros:
!
! use physics_share_f2c, only: cxx_pow, cxx_sqrt, cxx_cbrt, cxx_gamma, cxx_log, &
!                              cxx_log10, cxx_exp, cxx_tanh, cxx_erf

#ifndef SCREAM_BFB_MATH_INC
#define SCREAM_BFB_MATH_INC

#define bfb_square(val) ((val)*(val))
#define bfb_cube(val)   ((val)*(val)*(val))
#define bfb_quad(val)   (bfb_square(bfb_square(val)))

#ifdef SCREAM_CONFIG_IS_CMAKE
#  define bfb_pow(base, exp) cxx_pow(base, exp)
#  define bfb_sqrt(base) cxx_sqrt(base)
#  define bfb_cbrt(base) cxx_cbrt(base)
#  define bfb_gamma(val) cxx_gamma(val)
#  define bfb_log(val) cxx_log(val)
#  define bfb_log10(val) cxx_log10(val)
#  define bfb_exp(val) cxx_exp(val)
#  define bfb_expm1(val) cxx_expm1(val)
#  define bfb_tanh(val) cxx_tanh(val)
#  define bfb_erf(val) cxx_erf(val)
#else
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
#endif

#endif
