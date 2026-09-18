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

  // Get the unit vector components and magnitude at the surface.
  //
  // NOTE: this must run on *every* team thread, not inside a Kokkos::single.
  // xv and yv are thread-private scalars owned by the caller (locals in the
  // team lambda in GWDrag::run_impl), so a single-thread write leaves every
  // other thread in the team with an uninitialized xv/yv, and no team_barrier
  // can fix that -- a barrier synchronizes memory, it cannot broadcast a
  // register. Those threads then go on to compute ubm(k) = dot_2d(...,xv,yv)
  // and utgw(k) = ubt*xv for whichever levels they own, poisoning the result.
  // Running it redundantly is cheap and gives every thread the same value;
  // this matches gw_front_project_winds and gw_oro_src. The ubi() store below
  // becomes a same-value write from all threads, which is benign.
  get_unit_vector(u_src, v_src, xv, yv, ubi(init.k_src_wind + 1));

  team.team_barrier();

  // Project the local wind at midpoints onto the source wind.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, pver), [&] (const int k) {
      ubm(k) = dot_2d(u(k), v(k), xv, yv);
    });
  team.team_barrier();

  // Compute the interface wind projection by averaging the midpoint winds.
  // Use the top level wind at the top interface.
  ubi(0) = ubm(0);

  midpoint_interp(team, ubm, ekat::subview(ubi, Kokkos::pair<int, int>{1, pver}));
}

} // namespace gw
} // namespace scream

#endif
