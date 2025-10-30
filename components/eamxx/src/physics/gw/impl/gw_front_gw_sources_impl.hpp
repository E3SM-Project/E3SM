#ifndef GW_GW_FRONT_GW_SOURCES_IMPL_HPP
#define GW_GW_FRONT_GW_SOURCES_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_front_gw_sources. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_front_gw_sources(
  // Inputs
  const MemberType& team,
  const GwFrontInit& finit,
  const Int& pgwv,
  const Int& pver,
  const Int& kbot,
  const uview_1d<const Real>& frontgf,
  // Outputs
  const uview_2d<Real>& tau)
{
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, tau.size()), [&] (const int k) {
      tau.data()[k] = 0;
    });

  team.team_barrier();

  // GW generation depends on frontogenesis at specified level (may be below
  // actual source level).
  const bool launch_wave = frontgf(finit.kfront) > finit.frontgfc;

  const int num_pgwv = 2*pgwv + 1;
  if (launch_wave) {
    Kokkos::parallel_for(
      Kokkos::TeamVectorRange(team, num_pgwv), [&] (const int l) {
        tau(l, kbot+1)  = finit.fav(l);
      });
  }
}

} // namespace gw
} // namespace scream

#endif
