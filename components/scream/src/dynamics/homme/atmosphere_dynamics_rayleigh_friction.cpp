#include "atmosphere_dynamics.hpp"

// Scream includes
#include "share/util/scream_common_physics_functions.hpp"

namespace scream {

void HommeDynamics::rayleigh_friction_init()
{
  constexpr int N = HOMMEXX_PACK_SIZE;
  using KT = KokkosTypes<DefaultDevice>;
  using Pack = ekat::Pack<Real,N>;

  // Rayleigh friction paramaters
  m_rayk0     = m_params.get<int>("Rayleigh Friction Vertical Level", 2);
  m_raykrange = m_params.get<Real>("Rayleigh Friction Range", 0.0);
  m_raytau0   = m_params.get<Real>("Rayleigh Friction Decay Time", 5.0);

  // Calculate decay rate profile, otau.
  // If m_raytau0 == 0 no Rayleigh friction is applied.
  if (m_raytau0 != 0) {
    Real krange; // range of rayleigh friction profile
    Real tau0;   // approximate value of decay time at model top
    Real otau0;  // inverse of tau0

    krange = m_raykrange;
    if (m_raykrange == 0) krange = (m_rayk0 - 1.0)/2.0;

    tau0 = 86400.0*m_raytau0; // convert to seconds
    otau0 = 0;
    if (tau0 != 0) otau0 = 1.0/tau0;

    auto otau = m_buffer.otau;
    const auto nlevs = m_phys_grid->get_num_vertical_levels();
    const int npacks = ekat::PackInfo<Pack::n>::num_packs(nlevs);

    Kokkos::parallel_for(KT::RangePolicy(0, npacks),
                         KOKKOS_LAMBDA (const int ilev) {
      const auto range_pack = ekat::range<Pack>(ilev*Pack::n+1);
      const Pack x = (m_rayk0 - range_pack)/krange;
      otau(ilev) = otau0*(1.0 + ekat::tanh(x))/2.0;
    });
  }
}

void HommeDynamics::rayleigh_friction_apply(const Real dt) const
{
  constexpr int N = HOMMEXX_PACK_SIZE;
  using PF = PhysicsFunctions<DefaultDevice>;
  using KT = KokkosTypes<DefaultDevice>;
  using Pack = ekat::Pack<Real,N>;
  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;

  const auto ncols = m_phys_grid->get_num_local_dofs();
  const auto nlevs = m_phys_grid->get_num_vertical_levels();
  const auto npacks= ekat::PackInfo<N>::num_packs(nlevs);

  const auto horiz_winds_view = get_field_out("horiz_winds").get_view<Pack***>();
  const auto T_mid_view       = get_field_out("T_mid").get_view<Pack**>();

  const auto policy = ESU::get_default_team_policy(ncols, npacks);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int& icol = team.league_rank();

    // If m_raytau0 == 0 no Rayleigh friction is applied.
    if (m_raytau0 != 0) {
      const auto otau = m_buffer.otau;
      auto u_wind = ekat::subview(horiz_winds_view, icol, 0);
      auto v_wind = ekat::subview(horiz_winds_view, icol, 1);
      auto T_mid  = ekat::subview(T_mid_view, icol);
      PF::apply_rayleigh_friction(team, dt, otau, u_wind, v_wind, T_mid);
    }
  });
}

} // namespace scream
