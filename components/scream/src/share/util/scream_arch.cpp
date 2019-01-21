#ifdef FPE
# include <xmmintrin.h>
#endif

#ifdef _OPENMP
# include <omp.h>
#endif

#include "scream_arch.hpp"

/*
 * Implementations of scream_arch.hpp functions.
 */

namespace scream {
namespace util {

std::string active_avx_string () {
  std::string s;
#if defined __AVX512F__
  s += "-AVX512F";
#endif
#if defined __AVX2__
  s += "-AVX2";
#endif
#if defined __AVX__
  s += "-AVX";
#endif
  return s;
}

void dump_arch () {
  printf("ARCH: dp %d avx %s FPE %d nthread %d packn %d\n",
#ifdef DOUBLE_PRECISION
         1,
#else
         0,
#endif
         util::active_avx_string().c_str(),
#ifdef FPE
         1,
#else
         0,
#endif
#ifdef KOKKOS_ENABLE_OPENMP
         Kokkos::OpenMP::concurrency(),
#elif defined _OPENMP
         omp_get_max_threads(),
#else
         1,
#endif
         SCREAM_PACK_SIZE
         );
}

void initialize () {
#ifdef FPE
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &
                         ~( _MM_MASK_INVALID |
                            _MM_MASK_DIV_ZERO |
                            _MM_MASK_OVERFLOW ));
#endif
}

#ifdef SCREAM_FPE
static unsigned int constexpr exceptions =
  _MM_MASK_INVALID |
  _MM_MASK_DIV_ZERO |
  _MM_MASK_OVERFLOW;
#endif

void activate_floating_point_exceptions_if_enabled () {
#ifdef SCREAM_FPE
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~exceptions);
#endif
}

void deactivate_floating_point_exceptions_if_enabled () {
#ifdef SCREAM_FPE
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() | exceptions);
#endif
}

} // namespace util
} // namespace scream
