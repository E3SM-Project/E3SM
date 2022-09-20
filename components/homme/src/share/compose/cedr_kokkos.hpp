// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_KOKKOS_HPP
#define INCLUDE_CEDR_KOKKOS_HPP

#include <Kokkos_Core.hpp>

#if defined KOKKOS_ENABLE_CUDA || defined KOKKOS_ENABLE_HIP || defined KOKKOS_ENABLE_SYCL
# define CEDR_ENABLE_GPU
# if defined KOKKOS_ENABLE_CUDA
typedef Kokkos::Cuda CedrGpuExeSpace;
typedef Kokkos::CudaSpace CedrGpuSpace;
# endif
# if defined KOKKOS_ENABLE_HIP
typedef Kokkos::Experimental::HIP CedrGpuExeSpace;
typedef Kokkos::Experimental::HIPSpace CedrGpuSpace;
# endif
# if defined KOKKOS_ENABLE_SYCL
typedef Kokkos::Experimental::SYCL CedrGpuExeSpace;
typedef Kokkos::Experimental::SYCL> CedrGpuSpace;
# endif
#endif

#define KIF KOKKOS_INLINE_FUNCTION

// Clarify that a class member type is meant to be private but is
// marked public for Cuda visibility.
#define PRIVATE_CUDA public
#define PROTECTED_CUDA public

#if defined KOKKOS_COMPILER_GNU
// See https://github.com/kokkos/kokkos-kernels/issues/129 
# define ConstExceptGnu
#else
# define ConstExceptGnu const
#endif

namespace cedr {
namespace impl {

// Turn a View's MemoryTraits (traits::memory_traits) into the equivalent
// unsigned int mask.
template <typename View>
struct MemoryTraitsMask {
  enum : unsigned int {
#ifndef KOKKOS_VERSION
    // pre-v3
    value = ((View::traits::memory_traits::RandomAccess ? Kokkos::RandomAccess : 0) |
             (View::traits::memory_traits::Atomic ? Kokkos::Atomic : 0) |
             (View::traits::memory_traits::Restrict ? Kokkos::Restrict : 0) |
             (View::traits::memory_traits::Aligned ? Kokkos::Aligned : 0) |
             (View::traits::memory_traits::Unmanaged ? Kokkos::Unmanaged : 0))
#else
    // >= v3
    value = ((View::traits::memory_traits::is_random_access ? Kokkos::RandomAccess : 0) |
             (View::traits::memory_traits::is_atomic ? Kokkos::Atomic : 0) |
             (View::traits::memory_traits::is_restrict ? Kokkos::Restrict : 0) |
             (View::traits::memory_traits::is_aligned ? Kokkos::Aligned : 0) |
             (View::traits::memory_traits::is_unmanaged ? Kokkos::Unmanaged : 0))
#endif
  };
};

// Make the input View Unmanaged, whether or not it already is. One might
// imagine that View::unmanaged_type would provide this.
//   Use: Unmanaged<ViewType>
template <typename View>
using Unmanaged =
  // Provide a full View type specification, augmented with Unmanaged.
  Kokkos::View<typename View::traits::scalar_array_type,
               typename View::traits::array_layout,
               typename View::traits::device_type,
               Kokkos::MemoryTraits<
                 // All the current values...
                 MemoryTraitsMask<View>::value |
                 // ... |ed with the one we want, whether or not it's
                 // already there.
                 Kokkos::Unmanaged> >;

template <typename View>
using Const = typename View::const_type;

template <typename View>
using ConstUnmanaged = Const<Unmanaged<View> >;

template <typename ExeSpace>
struct DeviceType {
  typedef Kokkos::Device<typename ExeSpace::execution_space,
                         typename ExeSpace::memory_space> type;
};

#ifdef CEDR_ENABLE_GPU
typedef Kokkos::Device<CedrGpuSpace::execution_space,
                       CedrGpuSpace::memory_space> DefaultDeviceType;

template <> struct DeviceType<CedrGpuExeSpace> {
  typedef DefaultDeviceType type;
};
#else
typedef Kokkos::Device<Kokkos::DefaultExecutionSpace::execution_space,
                       Kokkos::DefaultExecutionSpace::memory_space> DefaultDeviceType;
#endif

template <typename ES> struct OnGpu {
  enum : bool { value =
#ifdef COMPOSE_MIMIC_GPU
                true
#else
                false
#endif
  };
};
#ifdef CEDR_ENABLE_GPU
template <> struct OnGpu<CedrGpuExeSpace> { enum : bool { value = true }; };
#endif

template <typename ExeSpace = Kokkos::DefaultExecutionSpace>
struct ExeSpaceUtils {
  using TeamPolicy = Kokkos::TeamPolicy<ExeSpace>;
  using Member = typename TeamPolicy::member_type;
  static TeamPolicy get_default_team_policy (int outer, int inner) {
#ifdef COMPOSE_MIMIC_GPU
    const int max_threads =
#ifdef KOKKOS_ENABLE_OPENMP
      ExeSpace::concurrency()
#else
      1
#endif
      ;
    const int team_size = max_threads < 7 ? max_threads : 7;
    return TeamPolicy(outer, team_size, 1);
#else
    return TeamPolicy(outer, 1, 1);
#endif
}
};

#ifdef CEDR_ENABLE_GPU
template <>
struct ExeSpaceUtils<CedrGpuExeSpace> {
  using TeamPolicy = Kokkos::TeamPolicy<CedrGpuExeSpace>;
  using Member = typename TeamPolicy::member_type;
  static TeamPolicy get_default_team_policy (int outer, int inner) {
    return TeamPolicy(outer, std::min(128, 32*((inner + 31)/32)), 1);
  }
};
#endif

// GPU-friendly replacements for std::*.
template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
void swap (T& a, T& b) { const T tmp = a; a = b; b = tmp; }

template <typename Real> struct TypeTraits {};

template <> struct TypeTraits<double> {
  static constexpr double epsilon = 2.2204460492503131e-16;
  static constexpr double infinity = 1e308;
};

template <> struct TypeTraits<float> {
  static constexpr double epsilon = 1.1920928955078125e-07;
  static constexpr double infinity = 1e38;
};

} // namespace impl
} // namespace cedr

#endif
