#ifndef SCREAM_TYPES_HPP
#define SCREAM_TYPES_HPP

#include "ekat/ekat.hpp"
#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "scream_config.h"

namespace scream
{

// Scalar types
using ekat::Int;
#ifdef SCREAM_DOUBLE_PRECISION
using Real = double;
#else
using Real = float;
#endif

// Kokkos types
using ekat::KokkosTypes;
using ekat::DefaultDevice;
using ekat::HostDevice;
using ekat::Unmanaged;

// Miscellanea

// An enum to be used with object that have 'repository'-like behavior
enum class RepoState {
  Clean,
  Open,
  Closed
};

#ifdef NDEBUG
static constexpr bool SCREAM_BFB_TESTING = true;
#else
static constexpr bool SCREAM_BFB_TESTING = false;
#endif

/*
 * Utility function for handling floating point literals,
 * so that they match the scream precision. This is
 * especially useful for bfb tests agaisnt fortran,
 * to ensure that literals are not a source of round-off differences.
 */
template<typename T> KOKKOS_INLINE_FUNCTION
constexpr typename std::enable_if<std::is_arithmetic<T>::value,Real>::type
sp (const T val) {
  return Real(val);
}

// Micro-utility, that given an enum returns the underlying int.
// The only use of this is if you need to sort scoped enums.
template<typename EnumT>
KOKKOS_FUNCTION
constexpr typename
std::enable_if<std::is_enum<EnumT>::value,
               typename std::underlying_type<EnumT>::type
              >::type
etoi (const EnumT e) {
  return static_cast<typename std::underlying_type<EnumT>::type>(e);
}

} // namespace scream

#endif // SCREAM_TYPES_HPP
