/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_PHYSICAL_CONSTANTS_HPP
#define HOMMEXX_PHYSICAL_CONSTANTS_HPP

#include "Types.hpp"

namespace Homme
{

namespace PhysicalConstants
{
  constexpr Real Rwater_vapor  = 461.5;
  constexpr Real Cpwater_vapor = 1870.0;
  constexpr Real Rgas          = 287.04;
  constexpr Real cp            = 1005.0;
  constexpr Real kappa         = Rgas / cp;
// real Earth
  constexpr Real rearth0       = 6.376e6;
  constexpr Real rrearth0      = 1.0 / rearth0;
  constexpr Real g             = 9.80616;
  constexpr Real p0            = 100000;         // [mbar]

  constexpr Real Tref          = 288;

  constexpr Real Tref_lapse_rate = 0.0065;
};

} // namespace Homme

#endif // HOMMEXX_PHYSICAL_CONSTANTS_HPP
