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
  // Outputs
  Real& mcsp_shear)
{
  //----------------------------------------------------------------------------
  // Purpose: calculate shear for MCSP
  //----------------------------------------------------------------------------

  // Local variables
  Real storm_u = 0;         // u wind at storm reference level set by MCSP_storm_speed_pref

  //----------------------------------------------------------------------------
  // Interpolate wind to pressure level specified by MCSP_storm_speed_pref

  // Find the interpolation indices using parallel reduction
  Int k_below = pver - 1; // default to lowest level
  Kokkos::parallel_reduce(Kokkos::TeamVectorRange(team, pver - 1),
    [&] (const Int& k, Int& min_k) {
      if (state_pmid(k) < ZMC::MCSP_storm_speed_pref && state_pmid(k+1) >= ZMC::MCSP_storm_speed_pref) {
        if (k < min_k) {
          min_k = k;
        }
      }
    },
    Kokkos::Min<Int>(k_below));
  team.team_barrier();

  // Linear interpolation
  Kokkos::single(Kokkos::PerTeam(team), [&] () {
    if (state_pmid(0) >= ZMC::MCSP_storm_speed_pref) {
      storm_u = state_u(0);
    }
    else if (state_pmid(pver - 1) < ZMC::MCSP_storm_speed_pref) {
      storm_u = state_u(pver - 1);
    }
    else {
      EKAT_KERNEL_ASSERT(k_below < pver-1);
      const Real dpu = ZMC::MCSP_storm_speed_pref - state_pmid(k_below);
      const Real dpl = state_pmid(k_below+1) - ZMC::MCSP_storm_speed_pref;
      storm_u = (state_u(k_below)*dpl + state_u(k_below+1)*dpu) / (dpl + dpu);
    }

    //----------------------------------------------------------------------------
    // calculate low-level shear
    if (state_pmid(pver - 1) > ZMC::MCSP_storm_speed_pref) {
      mcsp_shear = storm_u - state_u(pver - 1);
    } else {
      mcsp_shear = -999;
    }
  });
}

} // namespace zm
} // namespace scream

#endif
