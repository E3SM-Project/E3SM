#ifndef GW_GW_CONVECT_PROJECT_WINDS_IMPL_HPP
#define GW_GW_CONVECT_PROJECT_WINDS_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

#include <ekat_subview_utils.hpp>

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_convect_project_winds. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_convect_project_winds(
  // Inputs
  const MemberType& team,
  const GwConvectInit& init,
  const Int& pver,
  const uview_1d<const Real>& u,
  const uview_1d<const Real>& v,
  // Outputs
  Real& xv,
  Real& yv,
  const uview_1d<Real>& ubm,
  const uview_1d<Real>& ubi)
{
  // source wind speed and direction
  const Real u_src = u(init.k_src_wind);
  const Real v_src = v(init.k_src_wind);

  // Get the unit vector components and magnitude of the source wind.
  //
  // xv and yv are THREAD-PRIVATE outputs: every team thread computes its own
  // copy here, because every thread needs them (for ubm below, and for
  // utgw/vtgw in the caller). Callers must pass per-thread storage, not a
  // shared view element; a caller that wants to store them in a view must
  // publish them once itself (e.g. from a Kokkos::single). Computing them
  // inside a Kokkos::single instead would leave every other thread with an
  // uninitialized xv/yv -- a barrier synchronizes memory, it cannot broadcast
  // a register.
  //
  // The magnitude goes to the shared ubi array, so it is published once.
  Real src_wind_mag;
  get_unit_vector(u_src, v_src, xv, yv, src_wind_mag);
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    ubi(init.k_src_wind + 1) = src_wind_mag;
  });

  team.team_barrier();

  // Project the local wind at midpoints onto the source wind.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, pver), [&] (const int k) {
      ubm(k) = dot_2d(u(k), v(k), xv, yv);
    });
  team.team_barrier();

  // Compute the interface wind projection by averaging the midpoint winds.
  // Use the top level wind at the top interface.
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    ubi(0) = ubm(0);
  });

  midpoint_interp(team, ubm, ekat::subview(ubi, Kokkos::pair<int, int>{1, pver}));
}

} // namespace gw
} // namespace scream

#endif
