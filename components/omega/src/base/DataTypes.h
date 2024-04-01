#ifndef OMEGA_DATA_TYPES_H
#define OMEGA_DATA_TYPES_H
//===-- base/DataTypes.h - data type and array definitions ------*- C++ -*-===//
//
/// \file
/// \brief Defines standard data types and Kokkos array aliases
///
/// This header defines fixed-length data types to enforce levels of precision
/// where needed. In addition, it supplies a generic real type that is double
/// precision by default but can be switched throughout using a preprocessor
/// definition SINGLE_PRECISION. Finally, all arrays in OMEGA are defined
/// as Kokkos arrays to enable allocation and kernel launching on accelerator
/// devices. Because the Kokkos definitions can be lengthy, this header defines
/// useful aliases for up to 5D arrays in all supported types on either the
/// host or device.
//
//===----------------------------------------------------------------------===//

#include "Kokkos_Core.hpp"

namespace OMEGA {

// Standard integer and floating point types
using I4 = std::int32_t; ///< alias for 32-bit integer
using I8 = std::int64_t; ///< alias for 64-bit integer
using R4 = float;        ///< alias for 32-bit (single prec) real
using R8 = double;       ///< alias for 64-bit (double prec) real

/// generic real 64-bit (default) or 32-bit (if -DSINGLE_PRECISION used)
#ifdef SINGLE_PRECISION
using Real = float;
#else
using Real = double;
#endif

// user-defined literal for generic reals
KOKKOS_INLINE_FUNCTION constexpr Real operator""_Real(long double x) {
   return x;
}

// Aliases for Kokkos memory spaces
#ifdef OMEGA_ENABLE_CUDA
using MemSpace = Kokkos::CudaSpace;
#elif OMEGA_ENABLE_HIP
using MemSpace = Kokkos::Experimental::HIPSpace;
#elif OMEGA_ENABLE_OPENMP
using MemSpace = Kokkos::HostSpace;
#else
#error \
    "None of OMEGA_ENABLE_CUDA, OMEGA_ENABLE_HIP, and OMEGA_ENABLE_OPENMP is defined."
#endif

// Set default tile length
#ifndef OMEGA_TILE_LENGTH
#define OMEGA_TILE_LENGTH 64
#endif

// Aliases for Kokkos memory layouts
#ifdef OMEGA_LAYOUT_RIGHT

using MemLayout    = Kokkos::LayoutRight;
using MemInvLayout = Kokkos::LayoutLeft;

// Default tile configurations
template <int N> struct DefaultTile;

template <> struct DefaultTile<1> {
   static constexpr int value[] = {OMEGA_TILE_LENGTH};
};

template <> struct DefaultTile<2> {
   static constexpr int value[] = {1, OMEGA_TILE_LENGTH};
};

template <> struct DefaultTile<3> {
   static constexpr int value[] = {1, 1, OMEGA_TILE_LENGTH};
};

template <> struct DefaultTile<4> {
   static constexpr int value[] = {1, 1, 1, OMEGA_TILE_LENGTH};
};

template <> struct DefaultTile<5> {
   static constexpr int value[] = {1, 1, 1, 1, OMEGA_TILE_LENGTH};
};

#elif OMEGA_LAYOUT_LEFT

using MemLayout    = Kokkos::LayoutLeft;
using MemInvLayout = Kokkos::LayoutRight;

// Default tile configurations
template <int N> struct DefaultTile;

template <> struct DefaultTile<1> {
   static constexpr int value[] = {OMEGA_TILE_LENGTH};
};

template <> struct DefaultTile<2> {
   static constexpr int value[] = {OMEGA_TILE_LENGTH, 1};
};

template <> struct DefaultTile<3> {
   static constexpr int value[] = {OMEGA_TILE_LENGTH, 1, 1};
};

template <> struct DefaultTile<4> {
   static constexpr int value[] = {OMEGA_TILE_LENGTH, 1, 1, 1};
};

template <> struct DefaultTile<5> {
   static constexpr int value[] = {OMEGA_TILE_LENGTH, 1, 1, 1, 1};
};

#else
#error "OMEGA Memory Layout is not defined."
#endif

using HostMemSpace     = Kokkos::HostSpace;
using HostMemLayout    = MemLayout;
using HostMemInvLayout = MemInvLayout;

#define MAKE_OMEGA_VIEW_DIMS(N, V, T, ML, MS)  \
   using N##1D##T = Kokkos::V<T *, ML, MS>;    \
   using N##2D##T = Kokkos::V<T **, ML, MS>;   \
   using N##3D##T = Kokkos::V<T ***, ML, MS>;  \
   using N##4D##T = Kokkos::V<T ****, ML, MS>; \
   using N##5D##T = Kokkos::V<T *****, ML, MS>;

#define MAKE_OMEGA_VIEW_TYPES(N, V, ML, MS) \
   MAKE_OMEGA_VIEW_DIMS(N, V, I4, ML, MS)   \
   MAKE_OMEGA_VIEW_DIMS(N, V, I8, ML, MS)   \
   MAKE_OMEGA_VIEW_DIMS(N, V, R4, ML, MS)   \
   MAKE_OMEGA_VIEW_DIMS(N, V, R8, ML, MS)   \
   MAKE_OMEGA_VIEW_DIMS(N, V, Real, ML, MS)

// Aliases for Kokkos device arrays of various dimensions and types
MAKE_OMEGA_VIEW_TYPES(Array, View, MemLayout, MemSpace)

// Aliases for Kokkos host arrays of various dimensions and types
MAKE_OMEGA_VIEW_TYPES(HostArray, View, HostMemLayout, HostMemSpace)
} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
