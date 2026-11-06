/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_PHYSICAL_CONSTANTS_HPP
#define HOMMEXX_PHYSICAL_CONSTANTS_HPP

#include "Types.hpp"

// Numerical values come from the Physical Constants Dictionary (share/util_cxx/pcd_const.h)
#include "pcd_const.h"

namespace Homme
{

namespace PhysicalConstants
{
  constexpr Real Rwater_vapor  = pcd::water_vapor_specific_gas_constant_reference;
  constexpr Real Cpwater_vapor = pcd::water_vapor_specific_heat_capacity_reference;
  constexpr Real Rgas          = pcd::dry_air_specific_gas_constant_reference;
  constexpr Real cp            = pcd::dry_air_specific_heat_capacity_reference;
  constexpr Real rearth0       = pcd::mean_radius;
  constexpr Real g             = pcd::standard_acceleration_of_gravity;
  constexpr Real kappa         = Rgas / cp;
// real Earth
  constexpr Real rrearth0      = 1.0 / rearth0;
  constexpr Real p0            = 100000;         // [Pa] (no PCD equivalent: standard_atmosphere is 101325 Pa)

  constexpr Real Tref          = 288;

  constexpr Real Tref_lapse_rate = 0.0065;
};

} // namespace Homme

#endif // HOMMEXX_PHYSICAL_CONSTANTS_HPP
