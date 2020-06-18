#ifndef EKAT_ARCH_HPP
#define EKAT_ARCH_HPP

#include <string>
#include <Kokkos_Core.hpp>

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

#endif // EKAT_ARCH_HPP
