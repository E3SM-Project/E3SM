#include "share/property_checks/mass_and_energy_conservation_check.hpp"
#include "physics/share/physics_constants.hpp"
#include <iomanip>

namespace scream
{

MassAndEnergyConservationCheck::MassAndEnergyConservationCheck (const std::shared_ptr<const AbstractGrid>& grid,
                                                                const std::shared_ptr<const Field>&        pseudo_density_ptr,
                                                                const std::shared_ptr<const Field>&        ps_ptr,
                                                                const std::shared_ptr<const Field>&        phis_ptr,
                                                                const std::shared_ptr<const Field>&        horiz_winds_ptr,
                                                                const std::shared_ptr<const Field>&        T_mid_ptr,
                                                                const std::shared_ptr<const Field>&        qv_ptr,
                                                                const std::shared_ptr<const Field>&        qc_ptr,
                                                                const std::shared_ptr<const Field>&        qr_ptr,
                                                                const std::shared_ptr<const Field>&        qi_ptr,
                                                                const std::shared_ptr<const Field>&        vapor_flux_ptr,
                                                                const std::shared_ptr<const Field>&        water_flux_ptr,
                                                                const std::shared_ptr<const Field>&        ice_flux_ptr,
                                                                const std::shared_ptr<const Field>&        heat_flux_ptr)
  : m_grid (grid)
  , m_dt (std::nan(""))
  , m_tol (std::numeric_limits<Real>::max())
{  
  m_num_cols = m_grid->get_num_local_dofs();
  m_num_levs = m_grid->get_num_vertical_levels();

  m_current_mass   = view_1d<Real> ("current_total_water",  m_num_cols);
  m_current_energy = view_1d<Real> ("current_total_energy", m_num_cols);

  m_fields["pseudo_density"] = pseudo_density_ptr;
  m_fields["ps"]             = ps_ptr;
  m_fields["phis"]           = phis_ptr;
  m_fields["horiz_winds"]    = horiz_winds_ptr;
  m_fields["T_mid"]          = T_mid_ptr;
  m_fields["qv"]             = qv_ptr;
  m_fields["qc"]             = qc_ptr;
  m_fields["qr"]             = qr_ptr;
  m_fields["qi"]             = qi_ptr;
  m_fields["vapor_flux"]     = vapor_flux_ptr;
  m_fields["water_flux"]     = water_flux_ptr;
  m_fields["ice_flux"]       = ice_flux_ptr;
  m_fields["heat_flux"]      = heat_flux_ptr;

  // Require that all fields needed for mass and energy computation are not null.
  const bool all_computation_fields_exist =
      pseudo_density_ptr && ps_ptr          &&
      phis_ptr           && horiz_winds_ptr &&
      T_mid_ptr          && qv_ptr          &&
      qc_ptr             && qr_ptr          &&
      qi_ptr;
  EKAT_REQUIRE_MSG(all_computation_fields_exist,
                   "Error! Currently we require mass and energy conservation "
                   "check to contain all fields related to the mass and energy "
                   "computation.\n");

  // Require any process that add this checker to define all fluxes.
  // Fluxes which are not relevant to a certain process should be set to 0.
  EKAT_REQUIRE_MSG(vapor_flux_ptr && water_flux_ptr &&
                   ice_flux_ptr   && heat_flux_ptr,
                   "Error! If a process adds this check, it must define all "
                   "boundary fluxes. Fluxes which are not relevant to a "
                   "certain process should be set to 0.\n");
}

void MassAndEnergyConservationCheck::compute_current_mass ()
{
  auto mass = m_current_mass;

  const auto pseudo_density = m_fields.at("pseudo_density")->get_view<const Real**>();
  const auto qv = m_fields.at("qv")->get_view<const Real**>();
  const auto qc = m_fields.at("qc")->get_view<const Real**>();
  const auto qi = m_fields.at("qi")->get_view<const Real**>();
  const auto qr = m_fields.at("qr")->get_view<const Real**>();

  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, m_num_levs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();

    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto qv_i             = ekat::subview(qv, i);
    const auto qc_i             = ekat::subview(qc, i);
    const auto qi_i             = ekat::subview(qi, i);
    const auto qr_i             = ekat::subview(qr, i);

    mass(i) = compute_total_mass_on_column(team, pseudo_density_i, qv_i, qc_i, qi_i, qr_i);
  });
}

void MassAndEnergyConservationCheck::compute_current_energy ()
{
  auto energy = m_current_energy;

  const auto pseudo_density = m_fields.at("pseudo_density")->get_view<const Real**>();
  const auto T_mid = m_fields.at("T_mid")->get_view<const Real**>();
  const auto horiz_winds = m_fields.at("horiz_winds")->get_view<const Real***>();
  const auto qv = m_fields.at("qv")->get_view<const Real**>();
  const auto qc = m_fields.at("qc")->get_view<const Real**>();
  const auto qr = m_fields.at("qr")->get_view<const Real**>();
  const auto ps = m_fields.at("ps")->get_view<const Real*>();
  const auto phis = m_fields.at("phis")->get_view<const Real*>();

  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, m_num_levs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();

    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto T_mid_i          = ekat::subview(T_mid, i);
    const auto horiz_winds_i    = ekat::subview(horiz_winds, i);
    const auto qv_i             = ekat::subview(qv, i);
    const auto qc_i             = ekat::subview(qc, i);
    const auto qr_i             = ekat::subview(qr, i);

    energy(i) = compute_total_energy_on_column(team, pseudo_density_i, T_mid_i, horiz_winds_i,
                                               qv_i, qc_i, qr_i, ps(i), phis(i));
  });
}

PropertyCheck::ResultAndMsg MassAndEnergyConservationCheck::check() const
{
  auto mass   = m_current_mass;
  auto energy = m_current_energy;

  view_1d<Real> rel_err_mass  ("rel_err_mass",   m_num_cols);
  view_1d<Real> rel_err_energy("rel_err_energy", m_num_cols);

  EKAT_REQUIRE_MSG(!std::isnan(m_dt), "Error! Timestep dt must be set in MassAndEnergyConservationCheck "
                                      "before running check().");
  auto dt = m_dt;

  const auto pseudo_density = m_fields.at("pseudo_density")->get_view<const Real**> ();
  const auto T_mid          = m_fields.at("T_mid"         )->get_view<const Real**> ();
  const auto horiz_winds    = m_fields.at("horiz_winds"   )->get_view<const Real***>();
  const auto qv             = m_fields.at("qv"            )->get_view<const Real**> ();
  const auto qc             = m_fields.at("qc"            )->get_view<const Real**> ();
  const auto qi             = m_fields.at("qi"            )->get_view<const Real**> ();
  const auto qr             = m_fields.at("qr"            )->get_view<const Real**> ();
  const auto ps             = m_fields.at("ps"            )->get_view<const Real*>  ();
  const auto phis           = m_fields.at("phis"          )->get_view<const Real*>  ();

  const auto vapor_flux = m_fields.at("vapor_flux")->get_view<const Real*>();
  const auto water_flux = m_fields.at("water_flux")->get_view<const Real*>();
  const auto ice_flux   = m_fields.at("ice_flux"  )->get_view<const Real*>();
  const auto heat_flux  = m_fields.at("heat_flux" )->get_view<const Real*>();

  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, m_num_levs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {

    const int i = team.league_rank();

    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto T_mid_i          = ekat::subview(T_mid, i);
    const auto horiz_winds_i    = ekat::subview(horiz_winds, i);
    const auto qv_i             = ekat::subview(qv, i);
    const auto qc_i             = ekat::subview(qc, i);
    const auto qi_i             = ekat::subview(qi, i);
    const auto qr_i             = ekat::subview(qr, i);

    // Calculate total mass
    const Real tm = compute_total_mass_on_column(team, pseudo_density_i, qv_i, qc_i, qi_i, qr_i);
    const Real previous_tm = mass(i);

    // Calculate expected total mass. Here, dt should be set to the timestep of the 
    // subcycle for the process call this check. This effectively scales the boundary
    // fluxes by 1/num_subcycles (dt = model_dt/num_subcycles) so that we only include
    // the expected change after one substep (not a full timestep).
    const Real tm_exp = previous_tm +
                        compute_mass_boundary_flux_on_column(vapor_flux(i), water_flux(i))*dt;

    // Calculate relative error of total mass
    const Real err_m = std::abs(tm-tm_exp);
    rel_err_mass(i) = err_m/previous_tm;

    // Calculate total energy
    const Real te = compute_total_energy_on_column(team, pseudo_density_i, T_mid_i, horiz_winds_i,
                                                   qv_i, qc_i, qr_i, ps(i), phis(i));
    const Real previous_te = energy(i);

    // Calculate expected total energy. See the comment above for an explanation of dt. 
    const Real te_exp = previous_te +
                        compute_energy_boundary_flux_on_column(vapor_flux(i), water_flux(i), ice_flux(i), heat_flux(i))*dt;

    // Calculate relative error of total energy
    const Real err_e = std::abs(te-te_exp);
    rel_err_energy(i) = err_e/previous_te;
  });

  // Use Kokkos::MaxLoc to find the largest error for both mass and energy
  // and test that they are below the tolerance.
  using maxloc_t = Kokkos::MaxLoc<Real, int>;
  using maxloc_value_t = typename maxloc_t::value_type;
  maxloc_value_t maxloc_mass;
  maxloc_value_t maxloc_energy;

  Kokkos::parallel_reduce("mass_reduce", m_num_cols, KOKKOS_LAMBDA(int i, maxloc_value_t& result) {
    if (rel_err_mass(i) > result.val) {
      result.val = rel_err_mass(i);
      result.loc = i;
    }
  }, maxloc_t(maxloc_mass));

  Kokkos::parallel_reduce("energy_reduce", m_num_cols, KOKKOS_LAMBDA(int i, maxloc_value_t& result) {
    if (rel_err_energy(i) > result.val) {
      result.val = rel_err_energy(i);
      result.loc = i;
    }
  }, maxloc_t(maxloc_energy));

  PropertyCheck::ResultAndMsg res_and_msg;
  if (maxloc_mass.val < m_tol && maxloc_energy.val < m_tol) {
    res_and_msg.result = CheckResult::Pass;
  } else {
    res_and_msg.result = CheckResult::Fail;
  }

  std::stringstream msg;
  if (res_and_msg.result == CheckResult::Pass) {
    msg << "Check passed.\n";
  } else {
    msg << "Check failed.\n";
  }
  msg << "  - check name: " << this->name() << "\n"
      << "  - error tolerance: " << m_tol << "\n";

  // Output relative errors with lat/lon information (if available)
  AbstractGrid::dofs_list_h_type gids = m_grid->get_dofs_gids_host();
  AbstractGrid::geo_view_h_type lat, lon;
  const bool has_latlon = m_grid->has_geometry_data("lat") && m_grid->has_geometry_data("lon");
  if (has_latlon) {
    lat = m_grid->get_geometry_data_host("lat");
    lon = m_grid->get_geometry_data_host("lon");
  }
  msg << "  - mass relative error: " << maxloc_mass.val << "\n"
      << "    - global dof: " << gids(maxloc_mass.loc) << "\n";
  if (has_latlon) {
    msg << "    - (lat, lon): (" << lat(maxloc_mass.loc) << ", " << lon(maxloc_mass.loc) << ")\n";
  }
  msg << "  - energy relative error: " << maxloc_energy.val << "\n"
      << "    - global dof: " << gids(maxloc_energy.loc) << "\n";
  if (has_latlon) {
    msg << "    - (lat, lon): (" << lat(maxloc_energy.loc) << ", " << lon(maxloc_energy.loc) << ")\n";
  }
  res_and_msg.msg = msg.str();

  return res_and_msg;
}


KOKKOS_INLINE_FUNCTION
Real MassAndEnergyConservationCheck::compute_total_mass_on_column (const KT::MemberType&       team,
                                                                   const uview_1d<const Real>& pseudo_density,
                                                                   const uview_1d<const Real>& qv,
                                                                   const uview_1d<const Real>& qc,
                                                                   const uview_1d<const Real>& qi,
                                                                   const uview_1d<const Real>& qr) const
{
  using PC = scream::physics::Constants<Real>;
  const Real gravit = PC::gravit;

  Real total_mass(0);
  ekat::ExeSpaceUtils<typename KT::ExeSpace>::parallel_reduce(team, 0, m_num_levs,
                                                              [&] (const int lev, Real& local_mass) {
    local_mass += (qv(lev)+
                   qc(lev)+
                   qi(lev)+
                   qr(lev))*pseudo_density(lev)/gravit;
  }, total_mass);

  return total_mass;
}

KOKKOS_INLINE_FUNCTION
Real MassAndEnergyConservationCheck::compute_mass_boundary_flux_on_column (const Real vapor_flux,
                                                                           const Real water_flux) const
{
  using PC = scream::physics::Constants<Real>;
  const Real RHO_H2O  = PC::RHO_H2O;

  return vapor_flux - water_flux*RHO_H2O;
}

KOKKOS_INLINE_FUNCTION
Real MassAndEnergyConservationCheck::compute_total_energy_on_column (const KT::MemberType&       team,
                                                                     const uview_1d<const Real>& pseudo_density,
                                                                     const uview_1d<const Real>& T_mid,
                                                                     const uview_2d<const Real>& horiz_winds,
                                                                     const uview_1d<const Real>& qv,
                                                                     const uview_1d<const Real>& qc,
                                                                     const uview_1d<const Real>& qr,
                                                                     const Real                  ps,
                                                                     const Real                  phis) const
{
  using PC = scream::physics::Constants<Real>;
  const Real LatVap = PC::LatVap;
  const Real LatIce = PC::LatIce;
  const Real gravit = PC::gravit;
  const Real Cpair  = PC::Cpair;

  Real total_energy(0);
  ekat::ExeSpaceUtils<typename KT::ExeSpace>::parallel_reduce(team, 0, m_num_levs,
                                                              [&] (const int lev, Real& local_energy) {
    const auto u2 = horiz_winds(0,lev)*horiz_winds(0,lev);
    const auto v2 = horiz_winds(1,lev)*horiz_winds(1,lev);

    local_energy += (T_mid(lev)*Cpair +
                     0.5*(u2+v2) +
                     (LatVap+LatIce)*qv(lev) +
                     LatIce*(qc(lev)+qr(lev)))*pseudo_density(lev)/gravit;
  }, total_energy);
  total_energy += phis*ps/gravit;

  return total_energy;
}

KOKKOS_INLINE_FUNCTION
Real MassAndEnergyConservationCheck::compute_energy_boundary_flux_on_column (const Real vapor_flux,
                                                                             const Real water_flux,
                                                                             const Real ice_flux,
                                                                             const Real heat_flux) const
{
  using PC = scream::physics::Constants<Real>;
  const Real LatVap = PC::LatVap;
  const Real LatIce = PC::LatIce;
  const Real RHO_H2O  = PC::RHO_H2O;

  return vapor_flux*(LatVap+LatIce) - (water_flux-ice_flux)*RHO_H2O*LatIce + heat_flux;
}

} // namespace scream
