#ifndef SHOC_PHYSICS_CONSTANTS_HPP
#define SHOC_PHYSICS_CONSTANTS_HPP

#include "share/scream_types.hpp"

#include <vector>

namespace scream {
namespace physics {
namespace shoc {

/*
 * Mathematical constants used by SHOC.
 */

template <typename Scalar>
struct Constants
{
  // Upper limit for mixing length [m]
  static constexpr Scalar maxlen        = 20000.0;
  // Mixing length scaling parameter
  static constexpr Scalar length_fac    = 0.5;
};

} // namespace shoc
} // namespace physics
} // namespace scream

#endif
