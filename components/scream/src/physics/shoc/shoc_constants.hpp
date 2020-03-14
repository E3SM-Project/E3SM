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
  // These were taken from the PhysConst class in
  // E3SM-Project/scream-docs/shoc-port/shocintr.py.
  static constexpr Scalar GGr    = 9.80616;
  static constexpr Scalar RAir   = 287.042;
  static constexpr Scalar RH2O   = 461.505;
  static constexpr Scalar CpAir  = 1004.64;
  static constexpr Scalar LatVap = 2501000.0;
  static constexpr Scalar LatIce = 3.337e5;
  static constexpr Scalar Karman = 0.4;
  static constexpr Scalar Avogad = 6.02214e26;
  static constexpr Scalar Boltz  = 1.38065e-23;
  static constexpr Scalar RGas   = Avogad * Boltz;
  static constexpr Scalar MwWv   = 18.016;
  static constexpr Scalar Rwv    = RGas / MwWv;
  static constexpr Scalar ZVir   = (Rwv / RAir) - 1.0;
  static constexpr Scalar P0     = 1e-5;
};

} // namespace shoc
} // namespace scream

#endif
