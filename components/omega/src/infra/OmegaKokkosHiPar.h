#ifndef OMEGA_KOKKOS_HIPAR_H
#define OMEGA_KOKKOS_HIPAR_H
//===-- infra/OmegaKokkosHiPar.h - Omega hierarchical parallelism wrappers
//------*- C++ -*-===//
//
/// \file
/// \brief Omega hierarchical parallelism wrappers
///
/// INTERNAL HEADER NOT MEANT TO BE INCLUDED DIRECTLY
//
//===--------------------------------------------------------------------------------------===//

namespace OMEGA {

using TeamPolicy      = Kokkos::TeamPolicy<ExecSpace>;
using TeamMember      = TeamPolicy::member_type;
using ScratchMemSpace = ExecSpace::scratch_memory_space;
using Kokkos::PerTeam;

template <class T>
using ArrayScratch1D =
    Kokkos::View<T *, ScratchMemSpace, Kokkos::MemoryUnmanaged>;

using ArrayScratch1DReal = ArrayScratch1D<Real>;
using ArrayScratch1DI4   = ArrayScratch1D<I4>;

/// team_size for hierarchical parallelism
#ifdef OMEGA_TARGET_DEVICE
constexpr int OMEGA_TEAMSIZE = 64;
#else
constexpr int OMEGA_TEAMSIZE = 1;
#endif

#define INNER_LAMBDA [=]
// #define INNER_LAMBDA [&]

// Helper struct for providing information about scratch memory requirements
// TeamScratch<Real, I4>(4, 8) stores the number of bytes needed for
// 4 values of type Real and 8 vals of type I4
template <class... T> struct TeamScratch {
   size_t BytesPerTeam = 0;

   TeamScratch() = default;

   template <class... ArgT> TeamScratch(ArgT... Args) {
      static_assert(sizeof...(ArgT) == sizeof...(T));
      ((BytesPerTeam += ArrayScratch1D<T>::shmem_size(Args)), ...);
   }
};

template <int N> struct LaunchConfig {
   std::array<int, N> UpperBounds;
   int TeamSize;
   size_t ScratchBytesPerTeam;

   template <class... T>
   LaunchConfig(const int (&UpperBoundsIn)[N], int TeamSize,
                const TeamScratch<T...> &Scratch)
       : TeamSize(TeamSize), ScratchBytesPerTeam(Scratch.BytesPerTeam) {
      std::copy(std::begin(UpperBoundsIn), std::end(UpperBoundsIn),
                std::begin(UpperBounds));
   }

   template <class... T>
   LaunchConfig(const int (&UpperBounds)[N], const TeamScratch<T...> &Scratch)
       : LaunchConfig(UpperBounds, OMEGA_TEAMSIZE, Scratch) {}

   LaunchConfig(const int (&UpperBounds)[N], int TeamSize)
       : LaunchConfig(UpperBounds, TeamSize, TeamScratch<>{}) {}

   LaunchConfig(const int (&UpperBounds)[N])
       : LaunchConfig(UpperBounds, OMEGA_TEAMSIZE, TeamScratch<>{}) {}
};

KOKKOS_INLINE_FUNCTION void teamBarrier(const TeamMember &Team) {
   Team.team_barrier();
}

KOKKOS_INLINE_FUNCTION decltype(auto) teamScratch(const TeamMember &Team) {
   return Team.team_scratch(0);
}

// parallelForOuter: with label and with launch config
template <int N, class F>
inline void parallelForOuter(const std::string &Label,
                             const LaunchConfig<N> &LConfig, F &&Functor) {

   auto LinFunctor =
       LinearIdxWrapper{std::forward<F>(Functor), LConfig.UpperBounds};
   int LinBound = 1;
   for (int Rank = 0; Rank < N; ++Rank) {
      LinBound *= LConfig.UpperBounds[Rank];
   }

   auto Policy = TeamPolicy(LinBound, LConfig.TeamSize);

   if (LConfig.ScratchBytesPerTeam > 0) {
      Policy.set_scratch_size(0, Kokkos::PerTeam(LConfig.ScratchBytesPerTeam));
   }

   Kokkos::parallel_for(
       Label, Policy, KOKKOS_LAMBDA(const TeamMember &Team) {
          const int TeamId = Team.league_rank();
          LinFunctor(TeamId, Team);
       });
}

// parallelForOuter: without label and with launch config
template <int N, class F>
inline void parallelForOuter(const LaunchConfig<N> &LConfig, F &&Functor) {
   parallelForOuter("", LConfig, std::forward<F>(Functor));
}

// parallelForOuter: with label and with array bounds
template <int N, class F>
inline void parallelForOuter(const std::string &Label,
                             const int (&UpperBounds)[N], F &&Functor) {
   parallelForOuter(Label, LaunchConfig(UpperBounds), std::forward<F>(Functor));
}

// parallelForOuter: without label and with array bounds
template <int N, class F>
inline void parallelForOuter(const int (&UpperBounds)[N], F &&Functor) {
   parallelForOuter("", LaunchConfig(UpperBounds), std::forward<F>(Functor));
}

// This struct is used to get the right accumulator type to be used in
// the outer parallel lambda based on the final reduction variable type.
// The final reduction variable can be either a reference to
// an arithmetic type (int&, Real&) or a Kokkos reducer (Kokkos::Max<Real>).
// We need to know this type because nvcc does not allow generic lambdas.
template <class T, class Enable = void> struct AccumTypeHelper;

template <class T>
struct AccumTypeHelper<T, std::enable_if_t<std::is_arithmetic_v<T>>> {
   using Type = T;
};

template <class T>
struct AccumTypeHelper<T, std::enable_if_t<Kokkos::is_reducer_v<T>>> {
   using Type = typename T::value_type;
};

template <class T> using AccumType = typename AccumTypeHelper<T>::Type;

// parallelReduceOuter: with label and with launch config
template <int N, class F, class... R>
inline void parallelReduceOuter(const std::string &Label,
                                const LaunchConfig<N> &LConfig, F &&Functor,
                                R &&...Reducers) {

   auto LinFunctor =
       LinearIdxWrapper{std::forward<F>(Functor), LConfig.UpperBounds};
   int LinBound = 1;
   for (int Rank = 0; Rank < N; ++Rank) {
      LinBound *= LConfig.UpperBounds[Rank];
   }

   auto Policy = TeamPolicy(LinBound, LConfig.TeamSize);
   if (LConfig.ScratchBytesPerTeam > 0) {
      Policy.set_scratch_size(0, Kokkos::PerTeam(LConfig.ScratchBytesPerTeam));
   }

   Kokkos::parallel_reduce(
       Label, Policy,
       KOKKOS_LAMBDA(const TeamMember &Team,
                     AccumType<std::remove_reference_t<R>> &...Accums) {
          const int TeamId = Team.league_rank();
          LinFunctor(TeamId, Team, Accums...);
       },
       std::forward<R>(Reducers)...);
}

// parallelReduceOuter: without label and with launch config
template <int N, class F, class... R>
inline void parallelReduceOuter(const LaunchConfig<N> &LConfig, F &&Functor,
                                R &&...Reducers) {
   parallelReduceOuter("", LConfig, std::forward<F>(Functor),
                       std::forward<R>(Reducers)...);
}

// parallelReduceOuter: with label and with array bounds
template <int N, class F, class... R>
inline void parallelReduceOuter(const std::string &Label,
                                const int (&UpperBounds)[N], F &&Functor,
                                R &&...Reducers) {
   parallelReduceOuter(Label, LaunchConfig(UpperBounds),
                       std::forward<F>(Functor), std::forward<R>(Reducers)...);
}

// parallelReduceOuter: without label and with array bounds
template <int N, class F, class... R>
inline void parallelReduceOuter(const int (&UpperBounds)[N], F &&Functor,
                                R &&...Reducers) {
   parallelReduceOuter("", LaunchConfig(UpperBounds), std::forward<F>(Functor),
                       std::forward<R>(Reducers)...);
}

// Inclusive range of indices
struct Range {
   int First;
   int Last;
};

// parallelForInner

template <class F>
KOKKOS_FUNCTION void parallelForInner(const TeamMember &Team, Range Rng,
                                      F &&Functor) {
   const auto Policy = Kokkos::TeamThreadRange(Team, Rng.First, Rng.Last + 1);
   Kokkos::parallel_for(Policy, std::forward<F>(Functor));
}

template <class F>
KOKKOS_FUNCTION void parallelForInner(const TeamMember &Team, int UpperBound,
                                      F &&Functor) {
   parallelForInner(Team, Range{0, UpperBound - 1}, std::forward<F>(Functor));
}

// parallelReduceInner

template <class F, class... R>
KOKKOS_FUNCTION void parallelReduceInner(const TeamMember &Team, Range Rng,
                                         F &&Functor, R &&...Reducers) {
   const auto Policy = Kokkos::TeamThreadRange(Team, Rng.First, Rng.Last + 1);
   Kokkos::parallel_reduce(Policy, std::forward<F>(Functor),
                           std::forward<R>(Reducers)...);
}

template <class F, class... R>
KOKKOS_FUNCTION void parallelReduceInner(const TeamMember &Team, int UpperBound,
                                         F &&Functor, R &&...Reducers) {
   parallelReduceInner(Team, Range{0, UpperBound - 1}, std::forward<F>(Functor),
                       std::forward<R>(Reducers)...);
}

// parallelScanInner

template <class F, class... R>
KOKKOS_FUNCTION void parallelScanInner(const TeamMember &Team, Range Rng,
                                       F &&Functor, R &&...Reducers) {
   const auto Policy = Kokkos::TeamThreadRange(Team, Rng.First, Rng.Last + 1);
   Kokkos::parallel_scan(Policy, std::forward<F>(Functor),
                         std::forward<R>(Reducers)...);
}

template <class F, class... R>
KOKKOS_FUNCTION void parallelScanInner(const TeamMember &Team, int UpperBound,
                                       F &&Functor, R &&...Reducers) {
   parallelScanInner(Team, Range{0, UpperBound - 1}, std::forward<F>(Functor),
                     std::forward<R>(Reducers)...);
}

// parallelSearchInner
// Given a functor taking an index and returning a bool this function
// returns the first index in the range [0, UpperBound) for which the input
// functor returns true. If no such index is found it returns -1
template <class F>
KOKKOS_FUNCTION void parallelSearchInner(const TeamMember &Team, Range Rng,
                                         F &&Functor, int &Idx) {
   static_assert(std::is_same_v<std::invoke_result_t<F, int>, bool>,
                 "parallelSearchInner requires a functor that takes an int and "
                 "returns bool");

   // There are different implementations for host and device since the
   // parallel_reduce version doesn't return early leading to performance loss
   // on CPUs
#ifndef OMEGA_TARGET_DEVICE
   Idx = -1;
   for (int I = Rng.First; I <= Rng.Last; ++I) {
      if (Functor(I)) {
         Idx = I;
         break;
      }
   }
#else
   const auto Policy = Kokkos::TeamThreadRange(Team, Rng.First, Rng.Last + 1);
   Kokkos::parallel_reduce(
       Policy,
       INNER_LAMBDA(int I, int &Accum) {
          if (I <= Accum && Functor(I)) {
             Accum = I;
          }
       },
       Kokkos::Min<int>(Idx));
   if (Idx == Kokkos::reduction_identity<int>::min()) {
      Idx = -1;
   }
#endif
}

template <class F>
KOKKOS_FUNCTION void parallelSearchInner(const TeamMember &Team, int UpperBound,
                                         F &&Functor, int &Idx) {
   parallelSearchInner(Team, Range{0, UpperBound - 1}, std::forward<F>(Functor),
                       Idx);
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
