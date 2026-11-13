#ifndef OMEGA_KOKKOS_FLATPAR_H
#define OMEGA_KOKKOS_FLATPAR_H
//===-- infra/OmegaKokkosFlatPar.h - Omega flat parallelism wrappers ------*-
// C++ -*-===//
//
/// \file
/// \brief Omega flat parallelism wrappers
///
/// INTERNAL HEADER NOT MEANT TO BE INCLUDED DIRECTLY
//
//===--------------------------------------------------------------------------------===//

namespace OMEGA {

using Bounds1D = Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<int>>;

#if OMEGA_LAYOUT_RIGHT

template <int N>
using Bounds = Kokkos::MDRangePolicy<
    ExecSpace, Kokkos::Rank<N, Kokkos::Iterate::Right, Kokkos::Iterate::Right>,
    Kokkos::IndexType<int>>;

#elif OMEGA_LAYOUT_LEFT

template <int N>
using Bounds = Kokkos::MDRangePolicy<
    ExecSpace, Kokkos::Rank<N, Kokkos::Iterate::Left, Kokkos::Iterate::Left>,
    Kokkos::IndexType<int>>;

#else

#error "OMEGA Memory Layout is not defined."

#endif

// parallelFor: with label
template <int N, class F>
inline void parallelFor(const std::string &Label, const int (&UpperBounds)[N],
                        F &&Functor) {
   if constexpr (N == 1) {
      const auto Policy = Bounds1D(0, UpperBounds[0]);
      Kokkos::parallel_for(Label, Policy, std::forward<F>(Functor));

   } else {
#ifdef OMEGA_TARGET_DEVICE
      // On device convert the functor to use one dimensional indexing and use
      // 1D RangePolicy
      auto LinFunctor = LinearIdxWrapper{std::forward<F>(Functor), UpperBounds};
      int LinBound    = 1;
      for (int Rank = 0; Rank < N; ++Rank) {
         LinBound *= UpperBounds[Rank];
      }
      const auto Policy = Bounds1D(0, LinBound);
      Kokkos::parallel_for(Label, Policy, std::move(LinFunctor));
#else
      // On host use MDRangePolicy
      const int LowerBounds[N] = {0};
      const auto Policy        = Bounds<N>(LowerBounds, UpperBounds);
      Kokkos::parallel_for(Label, Policy, std::forward<F>(Functor));
#endif
   }
}

// parallelFor: without label
template <int N, class F>
inline void parallelFor(const int (&UpperBounds)[N], F &&Functor) {
   parallelFor("", UpperBounds, std::forward<F>(Functor));
}

// parallelReduce: with label
template <int N, class F, class... R>
inline void parallelReduce(const std::string &Label,
                           const int (&UpperBounds)[N], F &&Functor,
                           R &&...Reducers) {
   if constexpr (N == 1) {
      const auto Policy = Bounds1D(0, UpperBounds[0]);
      Kokkos::parallel_reduce(Label, Policy, std::forward<F>(Functor),
                              std::forward<R>(Reducers)...);

   } else {

#ifdef OMEGA_TARGET_DEVICE
      // On device convert the functor to use one dimensional indexing and use
      // 1D RangePolicy
      auto LinFunctor = LinearIdxWrapper{std::forward<F>(Functor), UpperBounds};
      int LinBound    = 1;
      for (int Rank = 0; Rank < N; ++Rank) {
         LinBound *= UpperBounds[Rank];
      }
      const auto Policy = Bounds1D(0, LinBound);
      Kokkos::parallel_reduce(Label, Policy, std::move(LinFunctor),
                              std::forward<R>(Reducers)...);
#else
      // On host use MDRangePolicy
      const int LowerBounds[N] = {0};
      const auto Policy        = Bounds<N>(LowerBounds, UpperBounds);
      Kokkos::parallel_reduce(Label, Policy, std::forward<F>(Functor),
                              std::forward<R>(Reducers)...);
#endif
   }
}

// parallelReduce: without label
template <int N, class F, class... R>
inline void parallelReduce(const int (&UpperBounds)[N], F &&Functor,
                           R &&...Reducers) {
   parallelReduce("", UpperBounds, std::forward<F>(Functor),
                  std::forward<R>(Reducers)...);
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
