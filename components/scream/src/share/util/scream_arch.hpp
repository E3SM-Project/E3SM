#ifndef SCREAM_ARCH_HPP
#define SCREAM_ARCH_HPP

#include <string>

/*
 * Architecture-related calls
 */

#if defined __INTEL_COMPILER
# define vector_ivdep _Pragma("ivdep")
# ifdef _OPENMP
#  define vector_simd _Pragma("omp simd")
# else
#  define vector_simd _Pragma("simd")
# endif
#elif defined __GNUG__
# define vector_ivdep _Pragma("GCC ivdep")
# define vector_simd _Pragma("GCC ivdep")
# define restrict __restrict__
#else
# define vector_ivdep
# define vector_simd
# define restrict
#endif

namespace scream {
namespace util {

std::string active_avx_string();

void dump_arch();

void initialize();

void activate_floating_point_exceptions_if_enabled();

// Use only when the situation demands it.
void deactivate_floating_point_exceptions_if_enabled();

template <typename ExeSpace>
struct OnGpu { enum : bool { value = false }; };
#ifdef KOKKOS_ENABLE_CUDA
template <> struct OnGpu<Kokkos::Cuda> { enum : bool { value = true }; };
#endif

} // namespace util
} // namespace scream

#endif // SCREAM_ARCH_HPP
