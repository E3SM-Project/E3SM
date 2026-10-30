#ifndef RAYLEIGH_FRICTION_IMPL_HPP
#define RAYLEIGH_FRICTION_IMPL_HPP

#include "rayleigh_friction_functions.hpp" // for ETI only but harmless for GPU

#include <ekat_team_policy_utils.hpp>
#include <ekat_assert.hpp>

namespace scream {
namespace rayleigh_friction {

/*
 * Implementation of the Rayleigh friction functions. Clients should NOT
 * #include this file, but include rayleigh_friction_functions.hpp instead.
 */

template<typename S, typename D>
void Functions<S,D>::rayleigh_friction_init(
  const int&             nlevs,
  const int&             rayk0,
  const Scalar&          raykrange,
  const Scalar&          raytau0,
  const view_1d<Scalar>& otau)
{
  EKAT_REQUIRE_MSG(raytau0 >= 0,
    "Error! Rayleigh friction decay time must be non-negative.\n"
    "  - raytau0: " + std::to_string(raytau0) + "\n");
  EKAT_REQUIRE_MSG(static_cast<int>(otau.extent(0)) >= nlevs,
    "Error! Rayleigh friction otau view is too small.\n");

  // If raytau0 == 0, no Rayleigh friction is applied
  if (raytau0 == 0) {
    Kokkos::deep_copy(otau, 0);
    return;
  }

  // Width of the profile. The default (raykrange=0) gives x=2 at the model top.
  const Scalar krange = (raykrange == 0) ? (rayk0 - 1) / Scalar(2) : raykrange;
  EKAT_REQUIRE_MSG(krange > 0,
    "Error! Invalid Rayleigh friction profile width.\n"
    "  - rayk0:     " + std::to_string(rayk0) + "\n"
    "  - raykrange: " + std::to_string(raykrange) + "\n"
    "  If raykrange = 0, rayk0 must be > 1. Otherwise raykrange must be > 0.\n");

  // Convert decay time from days to seconds and invert
  const Scalar tau0  = 86400 * raytau0;
  const Scalar otau0 = 1 / tau0;

  const typename KT::RangePolicy policy (0, nlevs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int& k) {
    // rayk0 follows the EAM convention of 1-based level indices
    const Scalar x = (rayk0 - (k+1)) / krange;
    otau(k) = otau0 * (1 + Kokkos::tanh(x)) / 2;
  });
}

template<typename S, typename D>
void Functions<S,D>::rayleigh_friction_tend(
  const int&                   ncols,
  const int&                   nlevs,
  const Scalar&                dt,
  const view_1d<const Scalar>& otau,
  const view_3d<Scalar>&       horiz_winds,
  const view_2d<Scalar>&       T_mid)
{
  using C   = physics::Constants<Scalar>;
  using TPF = ekat::TeamPolicyFactory<typename KT::ExeSpace>;

  const Scalar cpair  = C::Cpair.value;
  const Scalar rztodt = 1 / dt;

  const auto policy = TPF::get_default_team_policy(ncols, nlevs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const typename KT::MemberType& team) {
    const int i = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&] (const int& k) {
      // Euler backward: u_new = u / (1 + otau*dt), so the tendency is c1*u.
      // c3*(u^2+v^2) is the kinetic energy dissipation rate, applied as heating.
      const Scalar c2 = 1 / (1 + otau(k)*dt);
      const Scalar c1 = -otau(k) * c2;
      const Scalar c3 = Scalar(0.5) * (1 - c2*c2) * rztodt;

      const Scalar u = horiz_winds(i,0,k);
      const Scalar v = horiz_winds(i,1,k);

      horiz_winds(i,0,k) = u + c1*u*dt;
      horiz_winds(i,1,k) = v + c1*v*dt;
      T_mid(i,k)        += c3*(u*u + v*v)*dt / cpair;
    });
  });
}

} // namespace rayleigh_friction
} // namespace scream

#endif // RAYLEIGH_FRICTION_IMPL_HPP
