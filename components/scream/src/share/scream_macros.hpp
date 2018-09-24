#ifndef SCREAM_MACROS_HPP
#define SCREAM_MACROS_HPP

#if defined __INTEL_COMPILER
# define vector_ivdep _Pragma("ivdep")
# ifdef _OPENMP
#  define vector_simd _Pragma("omp simd")
# else
#  define vector_simd _Pragma("simd")
# endif
# define vector_novec _Pragma("novector")
#elif defined __GNUG__
# define vector_ivdep _Pragma("GCC ivdep")
# define vector_simd _Pragma("GCC ivdep")
# define vector_novec
# define restrict __restrict__
#else
# define vector_ivdep
# define vector_simd
# define vector_novec
# define restrict
#endif

#endif // SCREAM_MACROS_HPP
