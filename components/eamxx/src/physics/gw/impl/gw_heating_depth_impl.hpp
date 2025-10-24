#ifndef GW_GW_HEATING_DEPTH_IMPL_HPP
#define GW_GW_HEATING_DEPTH_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_heating_depth. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_heating_depth(
  // Inputs
  const MemberType& team,
  const GwConvectInit& init,
  const Int& pver,
  const Real& maxq0_conversion_factor,
  const Real& hdepth_scaling_factor,
  const bool& use_gw_convect_old,
  const uview_1d<const Real>& zm,
  const uview_1d<const Real>& netdt,
  // Outputs
  Int& mini,
  Int& maxi,
  Real& hdepth,
  Real& maxq0_out,
  Real& maxq0)
{
  // Find indices for the top and bottom of the heating range.
  EKAT_KERNEL_REQUIRE_MSG(!use_gw_convect_old, "The old gw convect algorithm is not supported");

  //---------------------------------------------------------------------
  // cleaner version that addresses bug in original where heating max and
  // depth were too low whenever heating <=0 occurred in the middle of
  // the heating profile (ex. at the melting level)
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    mini = 0;
    maxi = 0;
    for (Int k = pver-1; k >= 0; --k) {
      if ( zm(k) < GWC::heating_altitude_max ) {
        if ( netdt(k) > 0 ) {
          // Set mini as first spot where heating rate is positive
          if (mini == 0) mini = k;
          maxi = k;
        }
      }
      else {
        if (mini == 0) mini = k;
        if (maxi == 0) maxi = k;
      }
    }

    // apply tunable scaling factor for the heating depth
    // Heating depth in km.
    hdepth = (zm(maxi) -zm(mini))/1000;

    // Confine depth to table range.
    hdepth = ekat::impl::min(hdepth, static_cast<Real>(init.maxh));

    hdepth *= hdepth_scaling_factor;
  });

  team.team_barrier();

  // Maximum heating rate.
  Kokkos::parallel_reduce(
    Kokkos::TeamVectorRange(team, maxi, mini+1), [&] (const int k, Real& lmax) {
      if (netdt(k) > lmax) {
        lmax = netdt(k);
      }
    }, Kokkos::Max<Real>(maxq0));

  team.team_barrier();

  Kokkos::single(Kokkos::PerTeam(team), [&] {
    //output max heating rate in K/day
    maxq0_out = maxq0 * 24 * 3600;

    // Multipy by conversion factor
    maxq0 *= maxq0_conversion_factor;
  });
}

} // namespace gw
} // namespace scream

#endif
