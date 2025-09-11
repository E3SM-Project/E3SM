#ifndef GW_MOMENTUM_ENERGY_CONSERVATION_IMPL_HPP
#define GW_MOMENTUM_ENERGY_CONSERVATION_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw momentum_energy_conservation. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::momentum_energy_conservation(
  // Inputs
  const MemberType& team,
  const Int& pver,
  const Int& tend_level,
  const Real& dt,
  const uview_2d<const Real>& taucd,
  const uview_1d<const Real>& pint,
  const uview_1d<const Real>& pdel,
  const uview_1d<const Real>& u,
  const uview_1d<const Real>& v,
  // Inputs/Outputs
  const uview_1d<Real>& dudt,
  const uview_1d<Real>& dvdt,
  const uview_1d<Real>& dsdt,
  const uview_1d<Real>& utgw,
  const uview_1d<Real>& vtgw,
  const uview_1d<Real>& ttgw)
{
  static constexpr Real half=0.5;

  // Total mass from ground to source level: rho*dz = dp/gravit
  Real dz = 0;
  Kokkos::parallel_reduce(
    Kokkos::TeamVectorRange(team, tend_level+1, pver), [&] (const int k, Real& lsum) {
    lsum += pdel(k) / C::gravit;
  }, Kokkos::Sum<Real>(dz));

  // Tendency for U & V below source level.
  const Real ut_dz = -(taucd(tend_level+1, GWC::east) +
                       taucd(tend_level+1, GWC::west))/dz;
  const Real vt_dz = -(taucd(tend_level+1, GWC::north) +
                       taucd(tend_level+1, GWC::south))/dz;

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, tend_level+1, pver), [&] (const int k) {
    dudt(k) += ut_dz;
    dvdt(k) += vt_dz;
    utgw(k) += ut_dz;
    vtgw(k) += vt_dz;
  });

  team.team_barrier();

  // Net gain/loss of total energy in the column.
  Real dE = 0;
  Kokkos::parallel_reduce(
    Kokkos::TeamVectorRange(team, 0, pver), [&] (const int k, Real& lsum) {
      lsum += pdel(k) * (dsdt(k) +
                         dudt(k)*(u(k)+dudt(k)*half*dt) +
                         dvdt(k)*(v(k)+dvdt(k)*half*dt) );
    }, Kokkos::Sum<Real>(dE));

  dE = dE/(pint(pver)-pint(tend_level+1));

  // Subtract net gain/loss of total energy below source level.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, tend_level+1, pver), [&] (const int k) {
    dsdt(k) -= dE;
    ttgw(k) -= dE;
  });
}

} // namespace gw
} // namespace scream

#endif
