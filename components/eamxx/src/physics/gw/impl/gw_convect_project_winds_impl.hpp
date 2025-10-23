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
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    get_unit_vector(u_src, v_src, xv, yv, ubi(init.k_src_wind + 1));
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
  ubi(0) = ubm(0);

  midpoint_interp(team, ubm, ekat::subview(ubi, Kokkos::pair<int, int>{1, pver}));
}

} // namespace gw
} // namespace scream

#endif
