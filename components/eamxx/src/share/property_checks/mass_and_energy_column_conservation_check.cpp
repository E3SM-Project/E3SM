#include "share/property_checks/mass_and_energy_column_conservation_check.hpp"
#include "physics/share/physics_constants.hpp"
#include <iomanip>

namespace scream
{

MassAndEnergyColumnConservationCheck::
MassAndEnergyColumnConservationCheck (const std::shared_ptr<const AbstractGrid>& grid,
                                      const Real                                 mass_error_tolerance,
                                      const Real                                 energy_error_tolerance,
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
  , m_mass_tol (mass_error_tolerance)
  , m_energy_tol (energy_error_tolerance)
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

void MassAndEnergyColumnConservationCheck::compute_current_mass ()
{
  auto mass = m_current_mass;
  const auto ncols = m_num_cols;
  const auto nlevs = m_num_levs;

  const auto pseudo_density = m_fields.at("pseudo_density")->get_view<const Real**>();
  const auto qv = m_fields.at("qv")->get_view<const Real**>();
  const auto qc = m_fields.at("qc")->get_view<const Real**>();
  const auto qi = m_fields.at("qi")->get_view<const Real**>();
  const auto qr = m_fields.at("qr")->get_view<const Real**>();

  const auto policy = ExeSpaceUtils::get_default_team_policy(ncols, nlevs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();

    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto qv_i             = ekat::subview(qv, i);
    const auto qc_i             = ekat::subview(qc, i);
    const auto qi_i             = ekat::subview(qi, i);
    const auto qr_i             = ekat::subview(qr, i);

    mass(i) = compute_total_mass_on_column(team, nlevs, pseudo_density_i, qv_i, qc_i, qi_i, qr_i);
  });
}

void MassAndEnergyColumnConservationCheck::compute_current_energy ()
{
  auto energy = m_current_energy;
  const auto ncols = m_num_cols;
  const auto nlevs = m_num_levs;

  const auto pseudo_density = m_fields.at("pseudo_density")->get_view<const Real**>();
  const auto T_mid = m_fields.at("T_mid")->get_view<const Real**>();
  const auto horiz_winds = m_fields.at("horiz_winds")->get_view<const Real***>();
  const auto qv = m_fields.at("qv")->get_view<const Real**>();
  const auto qc = m_fields.at("qc")->get_view<const Real**>();
  const auto qr = m_fields.at("qr")->get_view<const Real**>();
  const auto ps = m_fields.at("ps")->get_view<const Real*>();
  const auto phis = m_fields.at("phis")->get_view<const Real*>();

  const auto policy = ExeSpaceUtils::get_default_team_policy(ncols, nlevs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();

    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto T_mid_i          = ekat::subview(T_mid, i);
    const auto horiz_winds_i    = ekat::subview(horiz_winds, i);
    const auto qv_i             = ekat::subview(qv, i);
    const auto qc_i             = ekat::subview(qc, i);
    const auto qr_i             = ekat::subview(qr, i);

    energy(i) = compute_total_energy_on_column(team, nlevs, pseudo_density_i, T_mid_i, horiz_winds_i,
                                               qv_i, qc_i, qr_i, ps(i), phis(i));
  });
}

PropertyCheck::ResultAndMsg MassAndEnergyColumnConservationCheck::check() const
{
  auto mass   = m_current_mass;
  auto energy = m_current_energy;
  const auto ncols = m_num_cols;
  const auto nlevs = m_num_levs;

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

  // Use Kokkos::MaxLoc to find the largest error for both mass and energy
  using maxloc_t = Kokkos::MaxLoc<Real, int>;
  using maxloc_value_t = typename maxloc_t::value_type;
  maxloc_value_t maxloc_mass;
  maxloc_value_t maxloc_energy;

  // Mass error calculation
  const auto policy = ExeSpaceUtils::get_default_team_policy(ncols, nlevs);
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA (const KT::MemberType& team,
                                                 maxloc_value_t&       result) {
    const int i = team.league_rank();

    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto qv_i             = ekat::subview(qv, i);
    const auto qc_i             = ekat::subview(qc, i);
    const auto qi_i             = ekat::subview(qi, i);
    const auto qr_i             = ekat::subview(qr, i);

    // Calculate total mass
    const Real tm = compute_total_mass_on_column(team, nlevs, pseudo_density_i, qv_i, qc_i, qi_i, qr_i);
    const Real previous_tm = mass(i);

    // Calculate expected total mass. Here, dt should be set to the timestep of the
    // subcycle for the process that called this check. This effectively scales the boundary
    // fluxes by 1/num_subcycles (dt = model_dt/num_subcycles) so that we only include
    // the expected change after one substep (not a full timestep).
    const Real tm_exp = previous_tm +
                        compute_mass_boundary_flux_on_column(vapor_flux(i), water_flux(i))*dt;

    // Calculate relative error of total mass
    const Real rel_err_mass = std::abs(tm-tm_exp)/previous_tm;

    // Test relative error against current max value
    if (rel_err_mass > result.val) {
      result.val = rel_err_mass;
      result.loc = i;
    }
  }, maxloc_t(maxloc_mass));

  // Energy error calculation
  Kokkos::parallel_reduce(policy, KOKKOS_LAMBDA (const KT::MemberType& team,
                                                 maxloc_value_t&       result) {

    const int i = team.league_rank();

    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto T_mid_i          = ekat::subview(T_mid, i);
    const auto horiz_winds_i    = ekat::subview(horiz_winds, i);
    const auto qv_i             = ekat::subview(qv, i);
    const auto qc_i             = ekat::subview(qc, i);
    const auto qr_i             = ekat::subview(qr, i);

    // Calculate total energy
    const Real te = compute_total_energy_on_column(team, nlevs, pseudo_density_i, T_mid_i, horiz_winds_i,
                                                   qv_i, qc_i, qr_i, ps(i), phis(i));
    const Real previous_te = energy(i);

    // Calculate expected total energy. See the comment above for an explanation of dt.
    const Real te_exp = previous_te +
                        compute_energy_boundary_flux_on_column(vapor_flux(i), water_flux(i), ice_flux(i), heat_flux(i))*dt;

    // Calculate relative error of total energy
    const Real rel_err_energy = std::abs(te-te_exp)/previous_te;

    // Test relative error against current max value
    if (rel_err_energy > result.val) {
      result.val = rel_err_energy;
      result.loc = i;
    }
  }, maxloc_t(maxloc_energy));

  // Check if mass and/or energy values were below tolerance.
  const bool mass_below_tol   = (maxloc_mass.val   < m_mass_tol);
  const bool energy_below_tol = (maxloc_energy.val < m_energy_tol);

  PropertyCheck::ResultAndMsg res_and_msg;
  if (mass_below_tol && energy_below_tol) {
    // If both vals are below the tolerance, the check passes.
    res_and_msg.result = CheckResult::Pass;
    return res_and_msg;
  }

  // If one or more values is above the tolerance, the check fails.
  res_and_msg.result = CheckResult::Fail;

  // We output relative errors with lat/lon information (if available)
  using gid_t = AbstractGrid::gid_type;
  auto gids = m_grid->get_dofs_gids().get_view<const gid_t*,Host>();
  typename Field::view_host_t<const Real*> lat, lon;
  const bool has_latlon = m_grid->has_geometry_data("lat") && m_grid->has_geometry_data("lon");
  if (has_latlon) {
    lat = m_grid->get_geometry_data("lat").get_view<const Real*, Host>();
    lon = m_grid->get_geometry_data("lon").get_view<const Real*, Host>();
  }

  std::stringstream msg;
  msg << "Check failed.\n"
      << "  - check name: " << this->name() << "\n";
  if (not mass_below_tol) {
    msg << "  - mass error tolerance: " << m_mass_tol << "\n";
    msg << "  - mass relative error: " << maxloc_mass.val << "\n"
        << "    - global dof: " << gids(maxloc_mass.loc) << "\n";
    if (has_latlon) {
      msg << "    - (lat, lon): (" << lat(maxloc_mass.loc) << ", " << lon(maxloc_mass.loc) << ")\n";
    }
    res_and_msg.fail_loc_indices.resize(1,maxloc_mass.loc);
    res_and_msg.fail_loc_tags = m_fields.at("phis")->get_header().get_identifier().get_layout().tags();
  }
  if (not energy_below_tol) {
    msg << "  - energy error tolerance: " << m_energy_tol << "\n";
    msg << "  - energy relative error: " << maxloc_energy.val << "\n"
        << "    - global dof: " << gids(maxloc_energy.loc) << "\n";
    if (has_latlon) {
      msg << "    - (lat, lon): (" << lat(maxloc_energy.loc) << ", " << lon(maxloc_energy.loc) << ")\n";
    }
    res_and_msg.fail_loc_indices.resize(1,maxloc_energy.loc);
    res_and_msg.fail_loc_tags = m_fields.at("phis")->get_header().get_identifier().get_layout().tags();
  }

  res_and_msg.msg = msg.str();

  return res_and_msg;
}


KOKKOS_INLINE_FUNCTION
Real MassAndEnergyColumnConservationCheck::
compute_total_mass_on_column (const KT::MemberType&       team,
                              const int                   nlevs,
                              const uview_1d<const Real>& pseudo_density,
                              const uview_1d<const Real>& qv,
                              const uview_1d<const Real>& qc,
                              const uview_1d<const Real>& qi,
                              const uview_1d<const Real>& qr)
{
  using PC = scream::physics::Constants<Real>;

  const Real gravit = PC::gravit;

  return ExeSpaceUtils::parallel_reduce<Real>(team, 0, nlevs,
                                              [&] (const int lev, Real& local_mass) {
    local_mass += (qv(lev)+
                   qc(lev)+
                   qi(lev)+
                   qr(lev))*pseudo_density(lev)/gravit;
  });
}

KOKKOS_INLINE_FUNCTION
Real MassAndEnergyColumnConservationCheck::
compute_mass_boundary_flux_on_column (const Real vapor_flux,
                                      const Real water_flux)
{
  using PC = scream::physics::Constants<Real>;
  const Real RHO_H2O  = PC::RHO_H2O;

  return vapor_flux - water_flux*RHO_H2O;
}

KOKKOS_INLINE_FUNCTION
Real MassAndEnergyColumnConservationCheck::
compute_total_energy_on_column (const KT::MemberType&       team,
                                const int                   nlevs,
                                const uview_1d<const Real>& pseudo_density,
                                const uview_1d<const Real>& T_mid,
                                const uview_2d<const Real>& horiz_winds,
                                const uview_1d<const Real>& qv,
                                const uview_1d<const Real>& qc,
                                const uview_1d<const Real>& qr,
                                const Real                  ps,
                                const Real                  phis)
{
  using PC = scream::physics::Constants<Real>;
  const Real LatVap = PC::LatVap;
  const Real LatIce = PC::LatIce;
  const Real gravit = PC::gravit;
  const Real Cpair  = PC::Cpair;

  Real total_energy =
    ExeSpaceUtils::parallel_reduce<Real>(team, 0, nlevs,
                                         [&] (const int lev, Real& local_energy) {
    const auto u2 = horiz_winds(0,lev)*horiz_winds(0,lev);
    const auto v2 = horiz_winds(1,lev)*horiz_winds(1,lev);

    local_energy += (T_mid(lev)*Cpair +
                     0.5*(u2+v2) +
                     (LatVap+LatIce)*qv(lev) +
                     LatIce*(qc(lev)+qr(lev)))*pseudo_density(lev)/gravit;
  });
  total_energy += phis*ps/gravit;

  return total_energy;
}

KOKKOS_INLINE_FUNCTION
Real MassAndEnergyColumnConservationCheck::
compute_energy_boundary_flux_on_column (const Real vapor_flux,
                                        const Real water_flux,
                                        const Real ice_flux,
                                        const Real heat_flux)
{
  using PC = scream::physics::Constants<Real>;
  const Real LatVap = PC::LatVap;
  const Real LatIce = PC::LatIce;
  const Real RHO_H2O  = PC::RHO_H2O;

  return vapor_flux*(LatVap+LatIce) - (water_flux-ice_flux)*RHO_H2O*LatIce + heat_flux;
}

} // namespace scream
