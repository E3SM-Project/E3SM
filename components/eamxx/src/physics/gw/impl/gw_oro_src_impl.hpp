#ifndef GW_GW_ORO_SRC_IMPL_HPP
#define GW_GW_ORO_SRC_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU
#include "share/util/eamxx_utils.hpp"
#include <ekat_subview_utils.hpp>

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_oro_src. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_oro_src(
  // Inputs
  const MemberType& team,
  const GwCommonInit& init,
  const Int& pver,
  const Int& pgwv,
  const uview_1d<const Real>& u,
  const uview_1d<const Real>& v,
  const uview_1d<const Real>& t,
  const Real& sgh,
  const uview_1d<const Real>& pmid,
  const uview_1d<const Real>& pint,
  const uview_1d<const Real>& dpm,
  const uview_1d<const Real>& zm,
  const uview_1d<const Real>& nm,
  // Outputs
  Int& src_level,
  Int& tend_level,
  const uview_2d<Real>& tau,
  const uview_1d<Real>& ubm,
  const uview_1d<Real>& ubi,
  Real& xv,
  Real& yv,
  const uview_1d<Real>& c)
{
  //--------------------------------------------------------------------------
  // Average the basic state variables for the wave source over the depth of
  // the orographic standard deviation. Here we assume that the appropiate
  // values of wind, stability, etc. for determining the wave source are
  // averages over the depth of the atmosphere penterated by the typical
  // mountain.
  // Reduces to the bottom midpoint values when sgh=0, such as over ocean.
  //--------------------------------------------------------------------------

  const Real hdsp = 2 * sgh; // Surface streamline displacement height (2*sgh).

  Int k = pver-1;

  // Averages over source region.
  Real rsrc = pmid(k)/(C::Rair*t(k)) * dpm(k); // Density.
  Real usrc = u(k) * dpm(k); // Zonal wind.
  Real vsrc = v(k) * dpm(k); // Meridional wind.
  Real nsrc = nm(k)* dpm(k); // B-V frequency.

  Real rsrc_sum, usrc_sum, vsrc_sum, nsrc_sum;
  Kokkos::parallel_reduce(
    Kokkos::TeamVectorRange(team, pver/2 - 1, pver-1), [&] (const int k, Real& lrsrc, Real& lusrc, Real& lvsrc, Real& lnsrc) {
      if (hdsp > std::sqrt(zm(k)*zm(k+1))) {
        lrsrc += pmid(k) / (C::Rair*t(k))* dpm(k);
        lusrc += u(k) * dpm(k);
        lvsrc += v(k) * dpm(k);
        lnsrc += (nm(k)* dpm(k));
      }
    }, rsrc_sum, usrc_sum, vsrc_sum, nsrc_sum);

  Kokkos::parallel_reduce(
    Kokkos::TeamVectorRange(team, pver/2 - 1, pver-1), [&] (const int k, Int& lmin) {
      if (lmin > pver - 1) {
        lmin = pver - 1;
      }
      if (hdsp > std::sqrt(zm(k)*zm(k+1)) && lmin > k) {
        lmin = k;
      }
    }, Kokkos::Min<Int>(src_level));

  team.team_barrier();

  rsrc += rsrc_sum;
  usrc += usrc_sum;
  vsrc += vsrc_sum;
  nsrc += nsrc_sum;

  // Difference in interface pressure across source region.
  const Real dpsrc = pint(pver) - pint(src_level);

  rsrc = rsrc / dpsrc;
  usrc = usrc / dpsrc;
  vsrc = vsrc / dpsrc;
  nsrc = nsrc / dpsrc;

  // Get the unit vector components and magnitude at the surface.
  get_unit_vector(usrc, vsrc, xv, yv, ubi(pver));

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

  team.team_barrier();

  // Determine the orographic c=0 source term following McFarlane (1987).
  // Set the source top interface index to pver, if the orographic term is
  // zero.
  Real tauoro = 0;
  if ((ubi(pver) > GWC::orovmin) &&  (hdsp > GWC::orohmin)) {
    // Max orographic standard deviation to use.
    const Real sghmax = init.fcrit2 * bfb_square(ubi(pver) / nsrc);
    // c=0 stress from orography.
    tauoro = init.oroko2 * ekat::impl::min(bfb_square(hdsp), sghmax) * rsrc * nsrc * ubi(pver);
  }
  else {
    src_level = pver;
  }

  // Set the phase speeds and wave numbers in the direction of the source
  // wind. Set the source stress magnitude (positive only, note that the
  // sign of the stress is the same as (c-u).
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, src_level, pver+1), [&] (const int k) {
      tau(pgwv, k) = tauoro;
    });

  // Allow wind tendencies all the way to the model bottom.
  tend_level = pver - 1;

  // adjust to c indexing. Up to this point, src_level was used to index into 0:pver arrays
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    --src_level;
  });

  // No spectrum; phase speed is just 0.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, c.size()), [&] (const int k) {
      c(k) = 0;
    });
}

} // namespace gw
} // namespace scream

#endif
