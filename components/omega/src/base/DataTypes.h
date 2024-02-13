#ifndef OMEGA_DATA_TYPES_H
#define OMEGA_DATA_TYPES_H
//===-- base/DataTypes.h - data type and array definitions ------*- C++ -*-===//
//
/// \file
/// \brief Defines standard data types and YAKL array aliases
///
/// This header defines fixed-length data types to enforce levels of precision
/// where needed. In addition, it supplies a generic real type that is double
/// precision by default but can be switched throughout using a preprocessor
/// definition SINGLE_PRECISION. Finally, all arrays in OMEGA are defined
/// as YAKL arrays to enable allocation and kernel launching on accelerator
/// devices. Because the YAKL definitions can be lengthy, this header defines
/// useful aliases for up to 5D arrays in all supported types on either the
/// host or device.
//
//===----------------------------------------------------------------------===//

#include "YAKL.h"
#include <cstdint>

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
YAKL_INLINE constexpr Real operator""_Real(long double x) { return x; }

// Aliases for YAKL arrays - by default, all arrays are on the device and
// use C-ordering.
/// Aliases for YAKL device arrays of various dimensions and types
using Array1DI4   = yakl::Array<I4, 1, yakl::memDevice, yakl::styleC>;
using Array1DI8   = yakl::Array<I8, 1, yakl::memDevice, yakl::styleC>;
using Array1DR4   = yakl::Array<R4, 1, yakl::memDevice, yakl::styleC>;
using Array1DR8   = yakl::Array<R8, 1, yakl::memDevice, yakl::styleC>;
using Array1DReal = yakl::Array<Real, 1, yakl::memDevice, yakl::styleC>;
using Array2DI4   = yakl::Array<I4, 2, yakl::memDevice, yakl::styleC>;
using Array2DI8   = yakl::Array<I8, 2, yakl::memDevice, yakl::styleC>;
using Array2DR4   = yakl::Array<R4, 2, yakl::memDevice, yakl::styleC>;
using Array2DR8   = yakl::Array<R8, 2, yakl::memDevice, yakl::styleC>;
using Array2DReal = yakl::Array<Real, 2, yakl::memDevice, yakl::styleC>;
using Array3DI4   = yakl::Array<I4, 3, yakl::memDevice, yakl::styleC>;
using Array3DI8   = yakl::Array<I8, 3, yakl::memDevice, yakl::styleC>;
using Array3DR4   = yakl::Array<R4, 3, yakl::memDevice, yakl::styleC>;
using Array3DR8   = yakl::Array<R8, 3, yakl::memDevice, yakl::styleC>;
using Array3DReal = yakl::Array<Real, 3, yakl::memDevice, yakl::styleC>;
using Array4DI4   = yakl::Array<I4, 4, yakl::memDevice, yakl::styleC>;
using Array4DI8   = yakl::Array<I8, 4, yakl::memDevice, yakl::styleC>;
using Array4DR4   = yakl::Array<R4, 4, yakl::memDevice, yakl::styleC>;
using Array4DR8   = yakl::Array<R8, 4, yakl::memDevice, yakl::styleC>;
using Array4DReal = yakl::Array<Real, 4, yakl::memDevice, yakl::styleC>;
using Array5DI4   = yakl::Array<I4, 5, yakl::memDevice, yakl::styleC>;
using Array5DI8   = yakl::Array<I8, 5, yakl::memDevice, yakl::styleC>;
using Array5DR4   = yakl::Array<R4, 5, yakl::memDevice, yakl::styleC>;
using Array5DR8   = yakl::Array<R8, 5, yakl::memDevice, yakl::styleC>;
using Array5DReal = yakl::Array<Real, 5, yakl::memDevice, yakl::styleC>;

// Also need similar aliases for arrays on the host
/// Aliases for YAKL host arrays of various dimensions and types
using ArrayHost1DI4   = yakl::Array<I4, 1, yakl::memHost, yakl::styleC>;
using ArrayHost1DI8   = yakl::Array<I8, 1, yakl::memHost, yakl::styleC>;
using ArrayHost1DR4   = yakl::Array<R4, 1, yakl::memHost, yakl::styleC>;
using ArrayHost1DR8   = yakl::Array<R8, 1, yakl::memHost, yakl::styleC>;
using ArrayHost1DReal = yakl::Array<Real, 1, yakl::memHost, yakl::styleC>;
using ArrayHost2DI4   = yakl::Array<I4, 2, yakl::memHost, yakl::styleC>;
using ArrayHost2DI8   = yakl::Array<I8, 2, yakl::memHost, yakl::styleC>;
using ArrayHost2DR4   = yakl::Array<R4, 2, yakl::memHost, yakl::styleC>;
using ArrayHost2DR8   = yakl::Array<R8, 2, yakl::memHost, yakl::styleC>;
using ArrayHost2DReal = yakl::Array<Real, 2, yakl::memHost, yakl::styleC>;
using ArrayHost3DI4   = yakl::Array<I4, 3, yakl::memHost, yakl::styleC>;
using ArrayHost3DI8   = yakl::Array<I8, 3, yakl::memHost, yakl::styleC>;
using ArrayHost3DR4   = yakl::Array<R4, 3, yakl::memHost, yakl::styleC>;
using ArrayHost3DR8   = yakl::Array<R8, 3, yakl::memHost, yakl::styleC>;
using ArrayHost3DReal = yakl::Array<Real, 3, yakl::memHost, yakl::styleC>;
using ArrayHost4DI4   = yakl::Array<I4, 4, yakl::memHost, yakl::styleC>;
using ArrayHost4DI8   = yakl::Array<I8, 4, yakl::memHost, yakl::styleC>;
using ArrayHost4DR4   = yakl::Array<R4, 4, yakl::memHost, yakl::styleC>;
using ArrayHost4DR8   = yakl::Array<R8, 4, yakl::memHost, yakl::styleC>;
using ArrayHost4DReal = yakl::Array<Real, 4, yakl::memHost, yakl::styleC>;
using ArrayHost5DI4   = yakl::Array<I4, 5, yakl::memHost, yakl::styleC>;
using ArrayHost5DI8   = yakl::Array<I8, 5, yakl::memHost, yakl::styleC>;
using ArrayHost5DR4   = yakl::Array<R4, 5, yakl::memHost, yakl::styleC>;
using ArrayHost5DR8   = yakl::Array<R8, 5, yakl::memHost, yakl::styleC>;
using ArrayHost5DReal = yakl::Array<Real, 5, yakl::memHost, yakl::styleC>;

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif
