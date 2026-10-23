#ifndef SCREAM_COSP_CONSTANTS_HPP
#define SCREAM_COSP_CONSTANTS_HPP

#include "share/core/eamxx_types.hpp"

namespace scream {
namespace cosp {

/*
 * Constants for the C++ port of COSP.
 *
 * A note on literals: the Fortran COSP library is compiled WITHOUT -fdefault-real-8
 * (the EAMxx Fortran flags only reach targets linking scream_share), so any literal
 * without a _wp suffix in COSP is single precision, even when it is assigned to, or
 * compared against, a real(wp) variable. To stay BFB with the Fortran, such literals
 * are written here as float literals cast to Scalar, e.g. Scalar(0.7490f).
 */

template <typename Scalar>
struct Constants
{
  // ---------------------------------------------------------------------------------------
  // Hydrometeor flags used in frac_out (subcolumn cloud type), see scops.F90
  static constexpr int I_LSC = 1; // Large-scale (stratiform) cloud
  static constexpr int I_CVC = 2; // Convective cloud

  // ---------------------------------------------------------------------------------------
  // MODIS simulator parameters, set in cosp_modis_init (cosp_modis_interface.F90)
  static constexpr int    num_trial_res = 15;
  static constexpr Scalar re_water_min  = 4.0;
  static constexpr Scalar re_water_max  = 30.0;
  static constexpr Scalar re_ice_min    = 5.0;
  static constexpr Scalar re_ice_max    = 90.0;
  static constexpr int phaseIsNone         = 0;
  static constexpr int phaseIsLiquid       = 1;
  static constexpr int phaseIsIce          = 2;
  static constexpr int phaseIsUndetermined = 3;

  // Ice density used in MODIS phase partitioning (cosp_optics.F90, 0.93_wp)
  static constexpr Scalar ice_density = 0.93;
};

} // namespace cosp
} // namespace scream

#endif // SCREAM_COSP_CONSTANTS_HPP
