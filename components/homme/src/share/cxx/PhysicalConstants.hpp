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

struct PhysicalConstants
{
  static constexpr Real Rwater_vapor  = 461.5;
  static constexpr Real Cpwater_vapor = 1870.0;
  static constexpr Real Rgas          = 287.04;
  static constexpr Real cp            = 1005.0;
  static constexpr Real kappa         = Rgas / cp;
  static constexpr Real rrearth       = 1.0 / 6.376e6;
};

} // namespace Homme

#endif // HOMMEXX_PHYSICAL_CONSTANTS_HPP
