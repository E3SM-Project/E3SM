#ifndef OMEGA_KOKKOS_H
#define OMEGA_KOKKOS_H
//===-- base/OmegaKokkos.h - Omega extension of Kokkos ------*- C++ -*-===//
//
/// \file
/// \brief Extends Kokkos for Omega
///
/// This header extends Kokkos for Omega.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include <utility>

namespace OMEGA {

#define OMEGA_SCOPE(a, b) auto &a = b

using ExecSpace     = MemSpace::execution_space;
using HostExecSpace = HostMemSpace::execution_space;

template <typename V>
auto createHostMirrorCopy(const V &view)
    -> Kokkos::View<typename V::data_type, HostMemLayout, HostMemSpace> {
   return Kokkos::create_mirror_view_and_copy(HostExecSpace(), view);
}

template <typename V>
auto createDeviceMirrorCopy(const V &view)
    -> Kokkos::View<typename V::data_type, MemLayout, MemSpace> {
   return Kokkos::create_mirror_view_and_copy(ExecSpace(), view);
}

// function alias to follow Camel Naming Convention
template <typename D, typename S> void deepCopy(D &dst, const S &src) {
   Kokkos::deep_copy(dst, src);
}

template <typename E, typename D, typename S>
void deepCopy(E &space, D &dst, const S &src) {
   Kokkos::deep_copy(space, dst, src);
}

#if OMEGA_LAYOUT_RIGHT

template <int N, class... Args>
using Bounds = Kokkos::MDRangePolicy<
    ExecSpace, Kokkos::Rank<N, Kokkos::Iterate::Right, Kokkos::Iterate::Right>,
    Args...>;

#elif OMEGA_LAYOUT_LEFT

template <int N, class... Args>
using Bounds = Kokkos::MDRangePolicy<
    ExecSpace, Kokkos::Rank<N, Kokkos::Iterate::Left, Kokkos::Iterate::Left>,
    Args...>;

#else

#error "OMEGA Memory Layout is not defined."

#endif

// parallelFor: with label
template <int N, class F, class... Args>
inline void parallelFor(const std::string &label, const int (&upper_bounds)[N],
                        const F &f,
                        const int (&tile)[N] = DefaultTile<N>::value) {
   if constexpr (N == 1) {
      const auto policy = Kokkos::RangePolicy<Args...>(0, upper_bounds[0]);
      Kokkos::parallel_for(label, policy, f);

   } else {
      const int lower_bounds[N] = {0};
      const auto policy = Bounds<N, Args...>(lower_bounds, upper_bounds, tile);
      Kokkos::parallel_for(label, policy, f);
   }
}

// parallelFor: without label
template <int N, class F>
inline void parallelFor(const int (&upper_bounds)[N], const F &f,
                        const int (&tile)[N] = DefaultTile<N>::value) {
   parallelFor("", upper_bounds, f, tile);
}

// parallelReduce: with label
template <int N, class F, class R, class... Args>
inline void parallelReduce(const std::string &label,
                           const int (&upper_bounds)[N], const F &f,
                           R &&reducer,
                           const int (&tile)[N] = DefaultTile<N>::value) {
   if constexpr (N == 1) {
      const auto policy = Kokkos::RangePolicy<Args...>(0, upper_bounds[0]);
      Kokkos::parallel_reduce(label, policy, f, std::forward<R>(reducer));

   } else {
      const int lower_bounds[N] = {0};
      const auto policy = Bounds<N, Args...>(lower_bounds, upper_bounds, tile);
      Kokkos::parallel_reduce(label, policy, f, std::forward<R>(reducer));
   }
}

// parallelReduce: without label
template <int N, class F, class R, class... Args>
inline void parallelReduce(const int (&upper_bounds)[N], const F &f,
                           R &&reducer,
                           const int (&tile)[N] = DefaultTile<N>::value) {
   parallelReduce("", upper_bounds, f, std::forward<R>(reducer), tile);
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
