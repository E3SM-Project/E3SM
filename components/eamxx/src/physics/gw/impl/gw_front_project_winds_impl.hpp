#ifndef GW_GW_FRONT_PROJECT_WINDS_IMPL_HPP
#define GW_GW_FRONT_PROJECT_WINDS_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

#include <ekat_subview_utils.hpp>

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_front_project_winds. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_front_project_winds(
  // Inputs
  const MemberType& team,
  const Int& pver,
  const Int& kbot,
  const uview_1d<const Real>& u,
  const uview_1d<const Real>& v,
  // Outputs
  Real& xv,
  Real& yv,
  const uview_1d<Real>& ubm,
  const uview_1d<Real>& ubi)
{
  // Just use the source level interface values for the source wind speed
  // and direction (unit vector).
  const Real usrc = GWC::half*(u(kbot+1)+u(kbot));
  const Real vsrc = GWC::half*(v(kbot+1)+v(kbot));

  // Get the unit vector components and magnitude at the surface.
  get_unit_vector(usrc, vsrc, xv, yv, ubi(kbot+1));

  // Project the local wind at midpoints onto the source wind.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, kbot+1), [&] (const int k) {
      ubm(k) = dot_2d(u(k), v(k), xv, yv);
    });
  team.team_barrier();

  // Compute the interface wind projection by averaging the midpoint winds.
  // Use the top level wind at the top interface.
  ubi(0) = ubm(0);

  midpoint_interp(team,
                  ekat::subview(ubm, Kokkos::pair<int, int>{0, kbot+1}),
                  ekat::subview(ubi, Kokkos::pair<int, int>{1, kbot+1}));
}

} // namespace gw
} // namespace scream

#endif
