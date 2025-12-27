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
#include "Error.h"
#include <functional>
#include <type_traits>
#include <utility>

namespace OMEGA {

#define OMEGA_SCOPE(a, b) auto &a = b

/// An enum is used to provide a shorthand for determining the type of
/// field. These correspond to the supported Omega data types (Real will be
/// identical to R4 or R8 depending on settings)
enum class ArrayDataType { Unknown, I4, I8, R4, R8 };

/// An enum is used to identify the location of the data - currently
/// either the device (the default) or explicitly on the host. Both refers
/// to the CPU-only case where the host and device are identical.
enum class ArrayMemLoc { Unknown, Device, Host, Both };

// determine ArrayDataType from Kokkos array type
template <class T> constexpr ArrayDataType checkArrayType() {
   if (std::is_same_v<typename T::non_const_value_type, I4>) {
      return ArrayDataType::I4;
   }

   if (std::is_same_v<typename T::non_const_value_type, I8>) {
      return ArrayDataType::I8;
   }

   if (std::is_same_v<typename T::non_const_value_type, R4>) {
      return ArrayDataType::R4;
   }

   if (std::is_same_v<typename T::non_const_value_type, R8>) {
      return ArrayDataType::R8;
   }

   return ArrayDataType::Unknown;
}

// determine ArrayMemLoc from Kokkos array type
template <class T> constexpr ArrayMemLoc findArrayMemLoc() {
   if (std::is_same_v<MemSpace, HostMemSpace>) {
      return ArrayMemLoc::Both;
   } else if (T::is_hostspace) {
      return ArrayMemLoc::Host;
   } else {
      return ArrayMemLoc::Device;
   }
}

/// Struct template to specify the rank of a supported Array
template <class T> struct ArrayRank {
   static constexpr bool Is1D = T::rank == 1;
   static constexpr bool Is2D = T::rank == 2;
   static constexpr bool Is3D = T::rank == 3;
   static constexpr bool Is4D = T::rank == 4;
   static constexpr bool Is5D = T::rank == 5;
};

using ExecSpace       = MemSpace::execution_space;
using HostExecSpace   = HostMemSpace::execution_space;
using TeamPolicy      = Kokkos::TeamPolicy<ExecSpace>;
using TeamMember      = TeamPolicy::member_type;
using ScratchMemSpace = ExecSpace::scratch_memory_space;
using Kokkos::MemoryUnmanaged;
using Kokkos::PerTeam;
using Kokkos::TeamThreadRange;
using RealScratchArray =
    Kokkos::View<Real *, ScratchMemSpace, Kokkos::MemoryUnmanaged>;

/// team_size for hierarchical parallelism
#ifdef OMEGA_TARGET_DEVICE
constexpr int OMEGA_TEAMSIZE = 64;
#else
constexpr int OMEGA_TEAMSIZE = 1;
#endif

// Takes a functor that uses multidimensional indexing
// and converts it into one that also accepts linear index
template <class F, int Rank> struct LinearIdxWrapper : F {
   static_assert(Rank >= 1 && Rank <= 5, "LinearIdxWrapper supports ranks 1-5");
   using F::operator();

   LinearIdxWrapper(F &&Functor, const int (&Bounds)[Rank])
       : F(std::move(Functor)) {
      computeStrides(Bounds);
   }

   LinearIdxWrapper(const F &Functor, const int (&Bounds)[Rank]) : F(Functor) {
      computeStrides(Bounds);
   }

   void computeStrides(const int (&Bounds)[Rank]) {
      if constexpr (Rank > 1) {
         Strides[Rank - 2] = Bounds[Rank - 1];
         for (int I = Rank - 3; I >= 0; --I) {
            Strides[I] = Bounds[I + 1] * Strides[I + 1];
         }
      }
   }

   template <int N = Rank, class... Args>
   KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<N == 2>
   operator()(int Idx, Args &&...OtherArgs) const {
      const int I1 = Idx / Strides[0];
      const int I2 = Idx - I1 * Strides[0];

      (*this)(I1, I2, std::forward<Args>(OtherArgs)...);
   }

   template <int N = Rank, class... Args>
   KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<N == 3>
   operator()(int Idx, Args &&...OtherArgs) const {
      const int I1 = Idx / Strides[0];
      Idx -= I1 * Strides[0];
      const int I2 = Idx / Strides[1];
      const int I3 = Idx - I2 * Strides[1];

      (*this)(I1, I2, I3, std::forward<Args>(OtherArgs)...);
   }

   template <int N = Rank, class... Args>
   KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<N == 4>
   operator()(int Idx, Args &&...OtherArgs) const {
      const int I1 = Idx / Strides[0];
      Idx -= I1 * Strides[0];
      const int I2 = Idx / Strides[1];
      Idx -= I2 * Strides[1];
      const int I3 = Idx / Strides[2];
      const int I4 = Idx - I3 * Strides[2];

      (*this)(I1, I2, I3, I4, std::forward<Args>(OtherArgs)...);
   }

   template <int N = Rank, class... Args>
   KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<N == 5>
   operator()(int Idx, Args &&...OtherArgs) const {
      const int I1 = Idx / Strides[0];
      Idx -= I1 * Strides[0];
      const int I2 = Idx / Strides[1];
      Idx -= I2 * Strides[1];
      const int I3 = Idx / Strides[2];
      Idx -= I3 * Strides[2];
      const int I4 = Idx / Strides[3];
      const int I5 = Idx - I4 * Strides[3];

      (*this)(I1, I2, I3, I4, I5, std::forward<Args>(OtherArgs)...);
   }

// SYCL doesn't allow 0-length arrays so add one extra element even though
// it is not needed
#ifdef KOKKOS_ENABLE_SYCL
   int Strides[Rank];
#else
   int Strides[Rank - 1];
#endif
};

template <typename V>
auto createHostMirrorCopy(const V &View)
    -> Kokkos::View<typename V::data_type, HostMemLayout, HostMemSpace> {
   return Kokkos::create_mirror_view_and_copy(HostExecSpace(), View);
}

template <typename V>
auto createDeviceMirrorCopy(const V &View)
    -> Kokkos::View<typename V::data_type, MemLayout, MemSpace> {
   return Kokkos::create_mirror_view_and_copy(ExecSpace(), View);
}

// function alias to follow Camel Naming Convention
template <typename D, typename S> void deepCopy(D &&Dst, S &&Src) {
   Kokkos::deep_copy(std::forward<D>(Dst), std::forward<S>(Src));
}

template <typename E, typename D, typename S>
void deepCopy(E &Space, D &Dst, const S &Src) {
   Kokkos::deep_copy(Space, Dst, Src);
}

// Check if two arrays are identical
template <class ArrayTypeA, class ArrayTypeB>
bool arraysEqual(const ArrayTypeA &A, const ArrayTypeB &B) {
   OMEGA_REQUIRE(A.span_is_contiguous() && B.span_is_contiguous(),
                 "arraysEqual works only for contiguous arrays");
   OMEGA_REQUIRE(A.size() == B.size(),
                 "arrayEqual can only compare arrays of equal size");

   // This is a debug utility and not performance critical
   // so just copy to the host and compare there
   const auto AH = createHostMirrorCopy(A);
   const auto BH = createHostMirrorCopy(B);

   bool Equal = true;
   for (size_t I = 0; I < AH.size(); I++) {
      if (AH.data()[I] != BH.data()[I]) {
         Equal = false;
         break;
      }
   }
   return Equal;
}

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

/// Hierarchical parallelism wrappers

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

// parallelForInner
template <class F>
KOKKOS_FUNCTION void parallelForInner(const TeamMember &Team, int UpperBound,
                                      F &&Functor) {
   const auto Policy = TeamThreadRange(Team, UpperBound);
   Kokkos::parallel_for(Policy, std::forward<F>(Functor));
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

// parallelReduceInner
template <class F, class... R>
KOKKOS_FUNCTION void parallelReduceInner(const TeamMember &Team, int UpperBound,
                                         F &&Functor, R &&...Reducers) {
   const auto Policy = TeamThreadRange(Team, UpperBound);
   Kokkos::parallel_reduce(Policy, std::forward<F>(Functor),
                           std::forward<R>(Reducers)...);
}

// parallelScanInner
template <class F, class... R>
KOKKOS_FUNCTION void parallelScanInner(const TeamMember &Team, int UpperBound,
                                       F &&Functor, R &&...Reducers) {
   const auto Policy = TeamThreadRange(Team, UpperBound);
   Kokkos::parallel_scan(Policy, std::forward<F>(Functor),
                         std::forward<R>(Reducers)...);
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
