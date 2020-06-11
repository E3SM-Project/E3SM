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
#elif defined(__GNUG__) && !defined(__clang__)
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

// Annotate a loop with this symbol if vector_simd should work but
// currently does not due to a compiler issue. For example,
// compilation to reduce_min seems not to work in Intel 17.
#define vector_disabled vector_novec

#endif // SCREAM_MACROS_HPP
