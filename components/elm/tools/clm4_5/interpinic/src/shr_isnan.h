#undef COMMENT
#ifdef COMMENT
/*
! These compilers define the isnan function as a FORTRAN intrinsic
! (gfortran, g95, intel, pathscale, irix and hp)
*/
#endif

#if defined(LINUX) && (defined(__GFORTRAN__) || defined(__G95__) || defined(__INTEL_COMPILER) || (_LANGUAGE_FORTRAN90 == 1 && __unix == 1)) || defined(IRIX64) || defined(OSF1)
#define ISNAN_INTRINSIC
#endif

#ifdef COMMENT
/*
! These compilers have a FORTRAN intrinsic to detect NAN's
! (The above list+ IBM and SunOS)
*/
#endif

#if defined(AIX) || defined(ISNAN_INTRINSIC) || defined(SunOS)
#undef  NOFTN_INTRINSIC

#ifdef COMMENT
/*
! For other compilers -- link to the C isnan function
! (such as PGI or lahey)
*/
#endif

#else
#define NOFTN_INTRINSIC
#endif

#ifdef COMMENT
/*
! isnan is only defined for gfortran for version 4.3 and greater
*/
#endif

#ifdef __GFORTRAN__

#define GCC_VERSION (__GNUC__ * 100 + __GNUC_MINOR__)

#if GCC_VERSION < 403

#undef ISNAN_INTRINSIC
#define NOFTN_INTRINSIC

#endif

#endif
