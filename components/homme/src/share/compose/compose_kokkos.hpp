#ifndef INCLUDE_COMPOSE_KOKKOS_HPP
#define INCLUDE_COMPOSE_KOKKOS_HPP

#include "compose.hpp"

#include <limits>

namespace Kokkos {

#ifndef ConstExceptGnu
# if defined KOKKOS_COMPILER_GNU
#  define ConstExceptGnu
# else
#  define ConstExceptGnu const
# endif
#endif

template <typename View>
using Const = typename View::const_type;

// GPU-friendly replacements for std::*.
#if KOKKOS_VERSION < 30700
template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }
#endif
template <typename T> KOKKOS_INLINE_FUNCTION
void swap (T& a, T& b) { const auto tmp = a; a = b; b = tmp; }

template <typename Real> struct NumericTraits;

template <> struct NumericTraits<double> {
  KOKKOS_INLINE_FUNCTION static double epsilon () {
    return
#ifdef COMPOSE_ENABLE_GPU
      2.2204460492503131e-16
#else
      std::numeric_limits<double>::epsilon()
#endif
      ;
  }
};

template <> struct NumericTraits<float> {
  KOKKOS_INLINE_FUNCTION static float epsilon () {
    return
#ifdef COMPOSE_ENABLE_GPU
      1.1920928955078125e-07
#else
      std::numeric_limits<float>::epsilon()
#endif
      ;
  }
};

template <typename ExeSpace>
struct DeviceType {
  typedef Kokkos::Device<typename ExeSpace::execution_space,
                         typename ExeSpace::memory_space> type;
};

#ifdef COMPOSE_ENABLE_GPU
typedef Kokkos::Device<ComposeGpuSpace::execution_space,
                       ComposeGpuSpace::memory_space> DefaultDeviceType;

template <> struct DeviceType<ComposeGpuExeSpace> {
  typedef DefaultDeviceType type;
};
#else
typedef Kokkos::Device<Kokkos::DefaultExecutionSpace::execution_space,
                       Kokkos::DefaultExecutionSpace::memory_space> DefaultDeviceType;
#endif

struct MachineTraits {
  // Host and device execution spaces.
#ifdef COMPOSE_PORT
  using HES = Kokkos::DefaultHostExecutionSpace;
  using DES = Kokkos::DefaultExecutionSpace;
#else
  using HES = Kokkos::Serial;
  using DES = Kokkos::Serial;
#endif
  using HDT = DeviceType<HES>::type;
  using DDT = DeviceType<DES>::type;
};

template <typename ES> struct OnGpu {
  enum : bool { value =
#ifdef COMPOSE_MIMIC_GPU
                true
#else
                false
#endif
  };
};
#ifdef COMPOSE_ENABLE_GPU
template <> struct OnGpu<ComposeGpuExeSpace> { enum : bool { value = true }; };
template <> struct OnGpu<MachineTraits> {}; // flag as an error at compile time
#endif

template <typename MT> using EnableIfOnGpu
  = typename std::enable_if<Kokkos::OnGpu<typename MT::DES>::value>::type;
template <typename MT> using EnableIfNotOnGpu
  = typename std::enable_if< ! Kokkos::OnGpu<typename MT::DES>::value>::type;

template <typename MT> struct SameSpace {
  enum { value = std::is_same<typename MT::HES, typename MT::DES>::value };
};
template <typename MT> using EnableIfSameSpace
  = typename std::enable_if<SameSpace<MT>::value>::type;
template <typename MT> using EnableIfDiffSpace
  = typename std::enable_if< ! SameSpace<MT>::value>::type;

// When Homme is running with HORIZ_OPENMP enabled, we can't do
//    const auto lv = gv; // gv is a managed view
// within the impl because everything in our impl is within the top-level
// HORIZ_OPENMP threaded region and the assignment above unsafely incr's and
// decr's the shared pointer count. Thus, we must carefully always obtain
// unmanaged Views to avoid the thread-unsafe ref-counting.

#ifdef COMPOSE_PORT
template <typename View> View unmanaged (
  const View& s, typename std::enable_if<View::rank_dynamic == 1>::type* = 0)
{ return View(s.data(), s.extent(0)); }
template <typename View> View unmanaged (
  const View& s, typename std::enable_if<View::rank_dynamic == 2>::type* = 0)
{ return View(s.data(), s.extent(0), s.extent(1)); }
template <typename View> View unmanaged (
  const View& s, typename std::enable_if<View::rank_dynamic == 3>::type* = 0)
{ return View(s.data(), s.extent(0), s.extent(1), s.extent(2)); }
template <typename View> View unmanaged (
  const View& s, typename std::enable_if<View::rank_dynamic == 4>::type* = 0)
{ return View(s.data(), s.extent(0), s.extent(1), s.extent(2), s.extent(3)); }
template <typename View> View unmanaged (
  const View& s, typename std::enable_if<View::rank_dynamic == 5>::type* = 0)
{ return View(s.data(), s.extent(0), s.extent(1), s.extent(2), s.extent(3), s.extent(4)); }
#else
template <typename View> const View& unmanaged (
  const View& s, typename std::enable_if<View::rank_dynamic>::type* = 0)
{ return s; }
#endif

// Copy by ref if not Cuda build.
#if defined COMPOSE_PORT && defined COMPOSE_ENABLE_GPU
# define COMPOSE_LAMBDA KOKKOS_LAMBDA
# define COMPOSE_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
# define COMPOSE_FORCEINLINE_FUNCTION KOKKOS_FORCEINLINE_FUNCTION
# define COMPOSE_FUNCTION KOKKOS_FUNCTION
#else
# define COMPOSE_LAMBDA [&]
# define COMPOSE_INLINE_FUNCTION inline
# define COMPOSE_FORCEINLINE_FUNCTION inline
# define COMPOSE_FUNCTION
#endif

template <typename V>
decltype(Kokkos::create_mirror_view(V())) cmvdc (const V& v) {
  const auto h = Kokkos::create_mirror_view(v);
  deep_copy(h, v);
  return h;
}

} // namespace Kokkos

#endif
