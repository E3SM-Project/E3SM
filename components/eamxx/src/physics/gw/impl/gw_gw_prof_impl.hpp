#ifndef GW_GW_PROF_IMPL_HPP
#define GW_GW_PROF_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

#include <ekat_subview_utils.hpp>

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_prof. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_prof(
  // Inputs
  const MemberType& team,
  const Int& pver,
  const Real& cpair,
  const uview_1d<const Real>& t,
  const uview_1d<const Real>& pmid,
  const uview_1d<const Real>& pint,
  // Outputs
  const uview_1d<Real>& rhoi,
  const uview_1d<Real>& ti,
  const uview_1d<Real>& nm,
  const uview_1d<Real>& ni)
{
  // Minimum value of Brunt-Vaisalla frequency squared.
  static constexpr Real n2min = 1.e-8;

  //-----------------------------------------------------------------------
  // Determine the interface densities and Brunt-Vaisala frequencies.
  //-----------------------------------------------------------------------

  // The top interface values are calculated assuming an isothermal
  // atmosphere above the top level.
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    ti(0) = t(0);
    rhoi(0) = pint(0) / (C::Rair*ti(0));
    ni(0) = sqrt(C::gravit*C::gravit / (cpair*ti(0)));
  });

  // Interior points use centered differences.
  midpoint_interp(team, t, ekat::subview(ti, Kokkos::pair<int, int>{1, pver}));
  team.team_barrier();
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, 1, pver), [&] (const int k) {
    rhoi(k) = pint(k) / (C::Rair*ti(k));
    const Real dtdp = (t(k)-t(k-1)) / (pmid(k)-pmid(k-1));
    const Real n2 = C::gravit*C::gravit/ti(k) * (1/cpair - rhoi(k)*dtdp);
    ni(k) = std::sqrt(ekat::impl::max(n2min, n2));
  });

  // Bottom interface uses bottom level temperature, density; next interface
  // B-V frequency.
  team.team_barrier();
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    ti(pver) = t(pver-1);
    rhoi(pver) = pint(pver) / (C::Rair*ti(pver));
    ni(pver) = ni(pver-1);
  });

  //------------------------------------------------------------------------
  // Determine the midpoint Brunt-Vaisala frequencies.
  //------------------------------------------------------------------------
  team.team_barrier();
  midpoint_interp(team, ni, nm);
}

} // namespace gw
} // namespace scream

#endif
