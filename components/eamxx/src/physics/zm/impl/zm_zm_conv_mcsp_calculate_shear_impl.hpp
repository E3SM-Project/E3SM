#ifndef ZM_ZM_CONV_MCSP_CALCULATE_SHEAR_IMPL_HPP
#define ZM_ZM_CONV_MCSP_CALCULATE_SHEAR_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_conv_mcsp_calculate_shear. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_conv_mcsp_calculate_shear(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const uview_1d<const Real>& state_pmid, // physics state mid-point pressure
  const uview_1d<const Real>& state_u, // physics state u momentum
  const uview_1d<const Real>& state_v, // physics state v momentum
  // Outputs
  Real& mcsp_shear)
{
  //----------------------------------------------------------------------------
  // Purpose: calculate shear for MCSP
  //----------------------------------------------------------------------------

  // Local variables
  Real storm_u = 0.0;         // u wind at storm reference level set by MCSP_storm_speed_pref
  Real storm_v = 0.0;         // v wind at storm reference level set by MCSP_storm_speed_pref
  Real storm_u_shear = 0.0;   // u shear at storm reference level set by MCSP_storm_speed_pref
  Real storm_v_shear = 0.0;   // v shear at storm reference level set by MCSP_storm_speed_pref

  //----------------------------------------------------------------------------
  // Interpolate wind to pressure level specified by MCSP_storm_speed_pref

  // Find the interpolation indices using parallel reduction
  Int k_below = pver - 1; // default to lowest level
  Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, pver - 1),
    [&] (const Int& k, Int& min_k) {
      if (state_pmid(k) >= ZMC::MCSP_storm_speed_pref && state_pmid(k + 1) < ZMC::MCSP_storm_speed_pref) {
        if (k < min_k) {
          min_k = k;
        }
      }
    },
    Kokkos::Min<Int>(k_below));
  team.team_barrier();

  // Linear interpolation
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    if (k_below < pver - 1) {
      const Real p_above = state_pmid(k_below);
      const Real p_below = state_pmid(k_below + 1);
      const Real weight = (ZMC::MCSP_storm_speed_pref - p_below) / (p_above - p_below);

      storm_u = state_u(k_below + 1) + weight * (state_u(k_below) - state_u(k_below + 1));
      storm_v = state_v(k_below + 1) + weight * (state_v(k_below) - state_v(k_below + 1));
    } else {
      // Target pressure is below all model levels or exact match at lowest level
      storm_u = state_u(pver - 1);
      storm_v = state_v(pver - 1);
    }

    //----------------------------------------------------------------------------
    // calculate low-level shear
    if (state_pmid(pver - 1) > ZMC::MCSP_storm_speed_pref) {
      storm_u_shear = storm_u - state_u(pver - 1);
      storm_v_shear = storm_v - state_v(pver - 1);
    } else {
      storm_u_shear = -999.0;
      storm_v_shear = -999.0;
    }
    mcsp_shear = storm_u_shear;
  });
}

} // namespace zm
} // namespace scream

#endif
