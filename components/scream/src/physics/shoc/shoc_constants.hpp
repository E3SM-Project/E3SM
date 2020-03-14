#ifndef SHOC_CONSTANTS_HPP
#define SHOC_CONSTANTS_HPP

#include "share/scream_types.hpp"

#include <vector>

namespace scream {
namespace shoc {

/*
 * Mathematical constants used by SHOC.
 *
 * Note that a potential optimization could be to change the type of
 * Scalar constants that have integer values to int.
 */

template <typename Scalar>
struct Constants
{
  static constexpr Scalar Ggr     = 1004.64;
  static constexpr Scalar RGas    = 1004.64;
  static constexpr Scalar Rv      = 1004.64;
  static constexpr Scalar Cp      = 1004.64;
  static constexpr Scalar LatCond = 1004.64;
  static constexpr Scalar LatIce  = 1004.64;
  static constexpr Scalar Eps     = 1004.64;
  static constexpr Scalar Vk      = 1004.64;
};

} // namespace shoc
} // namespace scream

#endif
