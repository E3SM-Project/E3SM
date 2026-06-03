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

constexpr I4 FillValueI4 = -2147483647;             ///< NC_FILL_INT
constexpr I8 FillValueI8 = -9223372036854775806LL;  ///< NC_FILL_INT64
constexpr R4 FillValueR4 = 9.9692099683868690e+36f; ///< NC_FILL_FLOAT
constexpr R8 FillValueR8 = 9.9692099683868690e+36;  ///< NC_FILL_DOUBLE
#if defined(SINGLE_PRECISION)
constexpr Real FillValueReal = FillValueR4;
#else
constexpr Real FillValueReal = FillValueR8;
#endif

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // OMEGA_FILL_VALUES_H
