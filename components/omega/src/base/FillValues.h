#ifndef OMEGA_FILL_VALUES_H
#define OMEGA_FILL_VALUES_H
//===-- base/FillValues.h - standard fill value constants -------*- C++ -*-===//
///
/// \file
/// \brief Standard fill value constants matching NetCDF-C NC_FILL_* values.
///
/// These values exactly match the SCORPIO PIO_FILL_* constants (which wrap
/// NC_FILL_* from netcdf.h). Defined here as literals so this header can be
/// included without pulling in pio.h or mpi.h, which have strict include-order
/// requirements that make them unsafe to include from generic headers.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"

namespace OMEGA {

/// Type-indexed fill value. Instantiating FillValue<T> for an unsupported
/// type T produces a link error (primary template is intentionally undefined).
template <typename T> constexpr T FillValue;

template <> inline constexpr I4 FillValue<I4> = -2147483647; ///< NC_FILL_INT
template <>
inline constexpr I8 FillValue<I8> = -9223372036854775806LL; ///< NC_FILL_INT64
template <>
inline constexpr R4 FillValue<R4> = 9.9692099683868690e+36f; ///< NC_FILL_FLOAT
template <>
inline constexpr R8 FillValue<R8> = 9.9692099683868690e+36; ///< NC_FILL_DOUBLE

/// Named aliases for readability in comparison and test code.
constexpr I4 FillValueI4 = FillValue<I4>;
constexpr I8 FillValueI8 = FillValue<I8>;
constexpr R4 FillValueR4 = FillValue<R4>;
constexpr R8 FillValueR8 = FillValue<R8>;
#if defined(SINGLE_PRECISION)
constexpr Real FillValueReal = FillValue<R4>;
#else
constexpr Real FillValueReal = FillValue<R8>;
#endif

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // OMEGA_FILL_VALUES_H
