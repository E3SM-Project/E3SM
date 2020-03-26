#ifndef SCREAM_ARCH_HPP
#define SCREAM_ARCH_HPP

#include <string>

/*
 * Architecture-related calls
 */

namespace scream {
namespace util {

std::string active_avx_string();

std::string config_string();

template <typename ExeSpace>
struct OnGpu { enum : bool { value = false }; };
#ifdef KOKKOS_ENABLE_CUDA
template <> struct OnGpu<Kokkos::Cuda> { enum : bool { value = true }; };
#endif

} // namespace util
} // namespace scream

#endif // SCREAM_ARCH_HPP
