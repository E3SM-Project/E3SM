#include "atmosphere_dynamics.hpp"

// Scream includes
#include "share/util/scream_common_physics_functions.hpp"

// Homme includes
#include "dynamics/homme/homme_dimensions.hpp"

namespace scream {

void HommeDynamics::rayleigh_friction_init()
{
  // Rayleigh friction paramaters
  m_rayk0     = m_params.get<int>("rayleigh_friction_vertical_level", 2);
  m_raykrange = m_params.get<Real>("rayleigh_friction_range", 0.0);
  m_raytau0   = m_params.get<Real>("rayleigh_friction_decay_time", 5.0);

  // If m_raytau0==0, then no Rayleigh friction is applied. Return.
  if (m_raytau0 == 0) return;

  // Input file is read in 1-based indexing, convert 
  // to 0-based for computations
  m_rayk0 -= 1;

  constexpr int N = SCREAM_PACK_SIZE;
  using Pack = RPack<N>;

  // Calculate decay rate profile, otau.
  const int nlevs = m_dyn_grid->get_num_vertical_levels();
  const auto npacks= ekat::PackInfo<N>::num_packs(nlevs);
  m_otau = decltype(m_otau)("otau", npacks);

  // Local paramters
  Real krange; // range of rayleigh friction profile
  Real tau0;   // approximate value of decay time at model top
  Real otau0;  // inverse of tau0

  krange = m_raykrange;
  if (m_raykrange == 0) krange = m_rayk0/2.0;

  tau0 = 86400.0*m_raytau0; // convert to seconds
  otau0 = 0;
  if (tau0 != 0) otau0 = 1.0/tau0;

  // locals for lambda captures to avoid issues on GPU
  auto otau = m_otau;
  auto rayk0 = m_rayk0;
    
  Kokkos::parallel_for(KT::RangePolicy(0, npacks),
                       KOKKOS_LAMBDA (const int ilev) {
    const auto range_pack = ekat::range<Pack>(ilev*N);
    const Pack x = (rayk0 - range_pack)/krange;
    otau(ilev) = otau0*(1.0 + ekat::tanh(x))/2.0;
  });
}

void HommeDynamics::rayleigh_friction_apply(const Real dt) const
{
  constexpr int N = HOMMEXX_PACK_SIZE;
  using PF = PhysicsFunctions<DefaultDevice>;
  using Pack = RPack<N>;
  using ESU = ekat::ExeSpaceUtils<KT::ExeSpace>;

  // If m_raytau0==0, then no Rayleigh friction is applied. Return.
  if (m_raytau0 == 0) return;

  const auto ncols = m_phys_grid->get_num_local_dofs();
  const auto nlevs = m_phys_grid->get_num_vertical_levels();
  const auto npacks= ekat::PackInfo<N>::num_packs(nlevs);

  const auto horiz_winds_view = get_field_out("horiz_winds").get_view<Pack***>();
  const auto T_mid_view       = get_field_out("T_mid").get_view<Pack**>();

  // local for lambda captures to avoid issues on GPU
  auto otau = m_otau;

  const auto policy = ESU::get_default_team_policy(ncols, npacks);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int& icol = team.league_rank();

    auto u_wind = ekat::subview(horiz_winds_view, icol, 0);
    auto v_wind = ekat::subview(horiz_winds_view, icol, 1);
    auto T_mid  = ekat::subview(T_mid_view, icol);
    PF::apply_rayleigh_friction(team, dt, otau, u_wind, v_wind, T_mid);
  });
}

} // namespace scream
