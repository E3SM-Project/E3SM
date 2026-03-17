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
using Kokkos::MemoryUnmanaged;
using Kokkos::PerTeam;
using Kokkos::TeamThreadRange;
using ArrayScratch1DReal =
    Kokkos::View<Real *, ScratchMemSpace, Kokkos::MemoryUnmanaged>;

/// team_size for hierarchical parallelism
#ifdef OMEGA_TARGET_DEVICE
constexpr int OMEGA_TEAMSIZE = 64;
#else
constexpr int OMEGA_TEAMSIZE = 1;
#endif

#define INNER_LAMBDA [=]
// #define INNER_LAMBDA [&]

KOKKOS_INLINE_FUNCTION void teamBarrier(const TeamMember &Team) {
   Team.team_barrier();
}

// parallelForOuter: with label
template <int N, class F>
inline void parallelForOuter(const std::string &Label,
                             const int (&UpperBounds)[N], F &&Functor,
                             int ScratchValsPerTeam = 0) {

   auto LinFunctor = LinearIdxWrapper{std::forward<F>(Functor), UpperBounds};
   int LinBound    = 1;
   for (int Rank = 0; Rank < N; ++Rank) {
      LinBound *= UpperBounds[Rank];
   }

   auto Policy = TeamPolicy(LinBound, OMEGA_TEAMSIZE);

   if (ScratchValsPerTeam > 0) {
      Policy.set_scratch_size(
          0, Kokkos::PerTeam(ScratchValsPerTeam * sizeof(Real)));
   }

   Kokkos::parallel_for(
       Label, Policy, KOKKOS_LAMBDA(const TeamMember &Team) {
          const int TeamId = Team.league_rank();
          LinFunctor(TeamId, Team);
       });
}

// parallelForOuter: without label
template <int N, class F>
inline void parallelForOuter(const int (&UpperBounds)[N], F &&Functor,
                             int ScratchValsPerTeam = 0) {
   parallelForOuter("", UpperBounds, std::forward<F>(Functor),
                    ScratchValsPerTeam);
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

// parallelReduceOuter: with label
template <int N, class F, class... R>
inline void parallelReduceOuter(const std::string &Label,
                                const int (&UpperBounds)[N], F &&Functor,
                                R &&...Reducers) {

   auto LinFunctor = LinearIdxWrapper{std::forward<F>(Functor), UpperBounds};
   int LinBound    = 1;
   for (int Rank = 0; Rank < N; ++Rank) {
      LinBound *= UpperBounds[Rank];
   }

   auto Policy = TeamPolicy(LinBound, OMEGA_TEAMSIZE);
   Kokkos::parallel_reduce(
       Label, Policy,
       KOKKOS_LAMBDA(const TeamMember &Team,
                     AccumType<std::remove_reference_t<R>> &...Accums) {
          const int TeamId = Team.league_rank();
          LinFunctor(TeamId, Team, Accums...);
       },
       std::forward<R>(Reducers)...);
}

// parallelReduceOuter: without label
template <int N, class F, class... R>
inline void parallelReduceOuter(const int (&UpperBounds)[N], F &&Functor,
                                R &&...Reducers) {
   parallelReduceOuter("", UpperBounds, std::forward<F>(Functor),
                       std::forward<R>(Reducers)...);
}

// parallelForInner

template <class F>
KOKKOS_FUNCTION void parallelForInner(const TeamMember &Team, int MinIndex,
                                      int MaxIndex, F &&Functor) {
   const auto Policy = TeamThreadRange(Team, MinIndex, MaxIndex + 1);
   Kokkos::parallel_for(Policy, std::forward<F>(Functor));
}

template <class F>
KOKKOS_FUNCTION void parallelForInner(const TeamMember &Team, int UpperBound,
                                      F &&Functor) {
   parallelForInner(Team, 0, UpperBound - 1, std::forward<F>(Functor));
}

// parallelReduceInner

template <class F, class... R>
KOKKOS_FUNCTION void parallelReduceInner(const TeamMember &Team, int MinIndex,
                                         int MaxIndex, F &&Functor,
                                         R &&...Reducers) {
   const auto Policy = TeamThreadRange(Team, MinIndex, MaxIndex + 1);
   Kokkos::parallel_reduce(Policy, std::forward<F>(Functor),
                           std::forward<R>(Reducers)...);
}

template <class F, class... R>
KOKKOS_FUNCTION void parallelReduceInner(const TeamMember &Team, int UpperBound,
                                         F &&Functor, R &&...Reducers) {
   parallelReduceInner(Team, 0, UpperBound - 1, std::forward<F>(Functor),
                       std::forward<R>(Reducers)...);
}

// parallelScanInner

template <class F, class... R>
KOKKOS_FUNCTION void parallelScanInner(const TeamMember &Team, int MinIndex,
                                       int MaxIndex, F &&Functor,
                                       R &&...Reducers) {
   const auto Policy = TeamThreadRange(Team, MinIndex, MaxIndex + 1);
   Kokkos::parallel_scan(Policy, std::forward<F>(Functor),
                         std::forward<R>(Reducers)...);
}

template <class F, class... R>
KOKKOS_FUNCTION void parallelScanInner(const TeamMember &Team, int UpperBound,
                                       F &&Functor, R &&...Reducers) {
   parallelScanInner(Team, 0, UpperBound - 1, std::forward<F>(Functor),
                     std::forward<R>(Reducers)...);
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
