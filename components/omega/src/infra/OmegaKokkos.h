#ifndef OMEGA_KOKKOS_H
#define OMEGA_KOKKOS_H
//===-- infra/OmegaKokkos.h - Omega extension of Kokkos ------*- C++ -*-===//
//
/// \file
/// \brief Extends Kokkos for Omega
///
/// This header extends Kokkos for Omega.
//
//===-------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Error.h"
#include <algorithm>
#include <array>
#include <functional>
#include <type_traits>
#include <utility>

namespace OMEGA {

#define OMEGA_SCOPE(a, b) auto &a = b

using ExecSpace     = MemSpace::execution_space;
using HostExecSpace = HostMemSpace::execution_space;

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
   if constexpr (std::is_same_v<MemSpace, HostMemSpace>) {
      return ArrayMemLoc::Both;
   } else if constexpr (std::is_same_v<typename T::memory_space,
                                       HostMemSpace>) {
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

// Array type traits: Map scalar types to concrete OMEGA array types
template <typename T> struct Array1D;
template <> struct Array1D<I4> {
   using type = Array1DI4;
};
template <> struct Array1D<I8> {
   using type = Array1DI8;
};
template <> struct Array1D<R4> {
   using type = Array1DR4;
};
template <> struct Array1D<R8> {
   using type = Array1DR8;
};

template <typename T> struct Array2D;
template <> struct Array2D<I4> {
   using type = Array2DI4;
};
template <> struct Array2D<I8> {
   using type = Array2DI8;
};
template <> struct Array2D<R4> {
   using type = Array2DR4;
};
template <> struct Array2D<R8> {
   using type = Array2DR8;
};

template <typename T> struct Array3D;
template <> struct Array3D<I4> {
   using type = Array3DI4;
};
template <> struct Array3D<I8> {
   using type = Array3DI8;
};
template <> struct Array3D<R4> {
   using type = Array3DR4;
};
template <> struct Array3D<R8> {
   using type = Array3DR8;
};

template <typename T> struct Array4D;
template <> struct Array4D<I4> {
   using type = Array4DI4;
};
template <> struct Array4D<I8> {
   using type = Array4DI8;
};
template <> struct Array4D<R4> {
   using type = Array4DR4;
};
template <> struct Array4D<R8> {
   using type = Array4DR8;
};

template <typename T> struct Array5D;
template <> struct Array5D<I4> {
   using type = Array5DI4;
};
template <> struct Array5D<I8> {
   using type = Array5DI8;
};
template <> struct Array5D<R4> {
   using type = Array5DR4;
};
template <> struct Array5D<R8> {
   using type = Array5DR8;
};

// Helper aliases for cleaner syntax
template <typename T> using Array1D_t = typename Array1D<T>::type;

template <typename T> using Array2D_t = typename Array2D<T>::type;

template <typename T> using Array3D_t = typename Array3D<T>::type;

template <typename T> using Array4D_t = typename Array4D<T>::type;

template <typename T> using Array5D_t = typename Array5D<T>::type;

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
                 "arraysEqual can only compare arrays of equal size");

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

// Works the same as Kokkos::subview, but returns an unmanaged view
// This is safe when we know that the subview won't outlive the input view
// and can be slightly faster
template <class View, class... Args>
KOKKOS_FUNCTION auto subviewUnmanaged(const View &A, Args... args) {
   return Kokkos::subview(View(A.data(), A.layout()), args...);
}

// Takes a functor that uses multidimensional indexing
// and converts it into one that also accepts linear index
template <class F, int Rank> class LinearIdxWrapper : private F {
   static_assert(Rank >= 1 && Rank <= 5, "LinearIdxWrapper supports ranks 1-5");

 public:
   template <class Array>
   LinearIdxWrapper(F &&Functor, Array &&Bounds) : F(std::move(Functor)) {
      computeStrides(std::forward<Array>(Bounds));
   }

   template <class Array>
   LinearIdxWrapper(const F &Functor, Array &&Bounds) : F(Functor) {
      computeStrides(std::forward<Array>(Bounds));
   }

   template <class Array> void computeStrides(Array &&Bounds) {
      if constexpr (Rank > 1) {
         Strides[Rank - 2] = Bounds[Rank - 1];
         for (int I = Rank - 3; I >= 0; --I) {
            Strides[I] = Bounds[I + 1] * Strides[I + 1];
         }
      }
   }

   template <int N = Rank, class... Args>
   KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<N == 1>
   operator()(int Idx, Args &&...OtherArgs) const {
      F::operator()(Idx, std::forward<Args>(OtherArgs)...);
   }

   template <int N = Rank, class... Args>
   KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<N == 2>
   operator()(int Idx, Args &&...OtherArgs) const {
      const int I1 = Idx / Strides[0];
      const int I2 = Idx - I1 * Strides[0];

      F::operator()(I1, I2, std::forward<Args>(OtherArgs)...);
   }

   template <int N = Rank, class... Args>
   KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<N == 3>
   operator()(int Idx, Args &&...OtherArgs) const {
      const int I1 = Idx / Strides[0];
      Idx -= I1 * Strides[0];
      const int I2 = Idx / Strides[1];
      const int I3 = Idx - I2 * Strides[1];

      F::operator()(I1, I2, I3, std::forward<Args>(OtherArgs)...);
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

      F::operator()(I1, I2, I3, I4, std::forward<Args>(OtherArgs)...);
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

      F::operator()(I1, I2, I3, I4, I5, std::forward<Args>(OtherArgs)...);
   }

 private:
// SYCL doesn't allow 0-length arrays so add one extra element even though
// it is not needed
#ifdef KOKKOS_ENABLE_SYCL
   int Strides[Rank];
#else
   int Strides[Rank - 1];
#endif
};

// Deduction guides for deducing Rank
template <class F, int Rank>
LinearIdxWrapper(F, const int (&)[Rank]) -> LinearIdxWrapper<F, Rank>;
template <class F, size_t Rank>
LinearIdxWrapper(F, std::array<int, Rank>) -> LinearIdxWrapper<F, Rank>;

} // end namespace OMEGA

// Flat parallelism wrappers
#include "OmegaKokkosFlatPar.h"

// Hierarchical parallelism wrappers
#include "OmegaKokkosHiPar.h"

//===----------------------------------------------------------------------===//
#endif
