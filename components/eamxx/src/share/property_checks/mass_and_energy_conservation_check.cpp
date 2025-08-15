#include "share/property_checks/mass_and_energy_conservation_check.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/field/field_utils.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_reduction_utils.hpp>

#include <iomanip>

namespace scream
{

MassAndEnergyConservationCheck::
MassAndEnergyConservationCheck (const ekat::Comm& comm,
                                const std::shared_ptr<const AbstractGrid>& grid,
                                const Real          mass_error_tolerance,
                                const Real          energy_error_tolerance,
                                const Field&        pseudo_density,
                                const Field&        ps,
                                const Field&        phis,
                                const Field&        horiz_winds,
                                const Field&        T_mid,
                                const Field&        qv,
                                const Field&        qc,
                                const Field&        qr,
                                const Field&        qi,
                                const Field&        vapor_flux,
                                const Field&        water_flux,
                                const Field&        ice_flux,
                                const Field&        heat_flux)
  : m_comm (comm)
  , m_grid (grid)
  , m_dt (std::nan(""))
  , m_mass_tol (mass_error_tolerance)
  , m_energy_tol (energy_error_tolerance)
{
  m_num_cols = m_grid->get_num_local_dofs();
  m_num_levs = m_grid->get_num_vertical_levels();

  m_current_mass   = view_1d<Real> ("current_total_water",  m_num_cols);
  m_current_energy = view_1d<Real> ("current_total_energy", m_num_cols);

  m_energy_change = view_1d<Real> ("energy change fixer", m_num_cols);

  m_fields["pseudo_density"] = pseudo_density;
  m_fields["ps"]             = ps;
  m_fields["phis"]           = phis;
  m_fields["horiz_winds"]    = horiz_winds;
  m_fields["T_mid"]          = T_mid;
  m_fields["qv"]             = qv;
  m_fields["qc"]             = qc;
  m_fields["qr"]             = qr;
  m_fields["qi"]             = qi;
  m_fields["vapor_flux"]     = vapor_flux;
  m_fields["water_flux"]     = water_flux;
  m_fields["ice_flux"]       = ice_flux;
  m_fields["heat_flux"]      = heat_flux;

  //allocate Fields for fixer reductions
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  //this does not exist, why? cause it can be a LEV or a COL array?
  //FieldLayout scalar1d = m_grid->get_1d_scalar_layout();
  FieldLayout scalar1d{{COL}, {m_num_cols}};
  FieldIdentifier s1_fid("s1", scalar1d, kg / kg, m_grid->name());
  field_version_s1 = Field(s1_fid);

  field_version_s1.allocate_view();

}

void MassAndEnergyConservationCheck::compute_current_mass ()
{
  using TPF = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;

  auto mass = m_current_mass;
  const auto ncols = m_num_cols;
  const auto nlevs = m_num_levs;

  const auto pseudo_density = m_fields.at("pseudo_density").get_view<const Real**>();
  const auto qv = m_fields.at("qv").get_view<const Real**>();
  const auto qc = m_fields.at("qc").get_view<const Real**>();
  const auto qi = m_fields.at("qi").get_view<const Real**>();
  const auto qr = m_fields.at("qr").get_view<const Real**>();

  const auto policy = TPF::get_default_team_policy(ncols, nlevs);
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

void MassAndEnergyConservationCheck::compute_current_energy ()
{
  using TPF = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;

  auto energy = m_current_energy;
  const auto ncols = m_num_cols;
  const auto nlevs = m_num_levs;

  const auto pseudo_density = m_fields.at("pseudo_density").get_view<const Real**>();
  const auto T_mid = m_fields.at("T_mid").get_view<const Real**>();
  const auto horiz_winds = m_fields.at("horiz_winds").get_view<const Real***>();
  const auto qv = m_fields.at("qv").get_view<const Real**>();
  const auto qc = m_fields.at("qc").get_view<const Real**>();
  const auto qr = m_fields.at("qr").get_view<const Real**>();
  const auto ps = m_fields.at("ps").get_view<const Real*>();
  const auto phis = m_fields.at("phis").get_view<const Real*>();

  const auto policy = TPF::get_default_team_policy(ncols, nlevs);
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

PropertyCheck::ResultAndMsg MassAndEnergyConservationCheck::check() const
{
  auto mass   = m_current_mass;
  auto energy = m_current_energy;
  const auto ncols = m_num_cols;
  const auto nlevs = m_num_levs;

  EKAT_REQUIRE_MSG(!std::isnan(m_dt), "Error! Timestep dt must be set in MassAndEnergyConservationCheck "
                                      "before running check().");
  auto dt = m_dt;

  const auto pseudo_density = m_fields.at("pseudo_density").get_view<const Real**> ();
  const auto T_mid          = m_fields.at("T_mid"         ).get_view<const Real**> ();
  const auto horiz_winds    = m_fields.at("horiz_winds"   ).get_view<const Real***>();
  const auto qv             = m_fields.at("qv"            ).get_view<const Real**> ();
  const auto qc             = m_fields.at("qc"            ).get_view<const Real**> ();
  const auto qi             = m_fields.at("qi"            ).get_view<const Real**> ();
  const auto qr             = m_fields.at("qr"            ).get_view<const Real**> ();
  const auto ps             = m_fields.at("ps"            ).get_view<const Real*>  ();
  const auto phis           = m_fields.at("phis"          ).get_view<const Real*>  ();

  const auto vapor_flux = m_fields.at("vapor_flux").get_view<const Real*>();
  const auto water_flux = m_fields.at("water_flux").get_view<const Real*>();
  const auto ice_flux   = m_fields.at("ice_flux"  ).get_view<const Real*>();
  const auto heat_flux  = m_fields.at("heat_flux" ).get_view<const Real*>();

  // Use Kokkos::MaxLoc to find the largest error for both mass and energy
  using maxloc_t = Kokkos::MaxLoc<Real, int>;
  using maxloc_value_t = typename maxloc_t::value_type;
  maxloc_value_t maxloc_mass;
  maxloc_value_t maxloc_energy;

  // Mass error calculation
  using TPF = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;
  const auto policy = TPF::get_default_team_policy(ncols, nlevs);
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
  using gid_type = AbstractGrid::gid_type;
  auto gids = m_grid->get_dofs_gids().get_view<const gid_type*,Host>();
  typename Field::view_host_t<const Real*> lat, lon;
  const bool has_latlon = m_grid->has_geometry_data("lat") && m_grid->has_geometry_data("lon");
  if (has_latlon) {
    lat = m_grid->get_geometry_data("lat").get_view<const Real*, Host>();
    lon = m_grid->get_geometry_data("lon").get_view<const Real*, Host>();
  }
  const bool has_additional_col_info = not additional_data_fields().empty();
  using namespace ShortFieldTagsNames;

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
    if (has_additional_col_info) {
      msg << "    - additional data (w/ local column index):\n";
      for (auto& f : additional_data_fields()) {
        f.sync_to_host();
        msg << "\n";
        print_field_hyperslab(f, {COL}, {maxloc_mass.loc}, msg);
      }
      msg << "\n    END OF ADDITIONAL DATA\n";
    }
    res_and_msg.fail_loc_indices.resize(1,maxloc_mass.loc);
    res_and_msg.fail_loc_tags = m_fields.at("phis").get_header().get_identifier().get_layout().tags();
  }
  if (not energy_below_tol) {
    msg << "  - energy error tolerance: " << m_energy_tol << "\n";
    msg << "  - energy relative error: " << maxloc_energy.val << "\n"
        << "    - global dof: " << gids(maxloc_energy.loc) << "\n";
    if (has_latlon) {
      msg << "    - (lat, lon): (" << lat(maxloc_energy.loc) << ", " << lon(maxloc_energy.loc) << ")\n";
    }
    if (has_additional_col_info) {
      msg << "    - additional data (w/ local column index):\n";
      for (auto& f : additional_data_fields()) {
        f.sync_to_host();
        msg << "\n";
        print_field_hyperslab(f, {COL}, {maxloc_energy.loc}, msg);
      }
      msg << "\n    END OF ADDITIONAL DATA\n";
    }
    res_and_msg.fail_loc_indices.resize(1,maxloc_energy.loc);
    res_and_msg.fail_loc_tags = m_fields.at("phis").get_header().get_identifier().get_layout().tags();
  }

  res_and_msg.msg = msg.str();

  return res_and_msg;
}

void MassAndEnergyConservationCheck::global_fixer(const bool & print_debug_info)
{
  using TPF = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;
  const auto ncols = m_num_cols;
  const auto nlevs = m_num_levs;

  //keep dt for fixers with fluxes.
  //we do not plan to use fluxes in dycore fixer, but the code is already there
  EKAT_REQUIRE_MSG(!std::isnan(m_dt), "Error! Timestep dt must be set in MassAndEnergyConservationCheck "
                                      "before running check().");
  auto dt = m_dt;

  const auto pseudo_density = m_fields.at("pseudo_density").get_view<const Real**> ();
  const auto T_mid          = m_fields.at("T_mid"         ).get_view<      Real**> ();
  const auto horiz_winds    = m_fields.at("horiz_winds"   ).get_view<const Real***>();
  const auto qv             = m_fields.at("qv"            ).get_view<const Real**> ();
  const auto qc             = m_fields.at("qc"            ).get_view<const Real**> ();
  const auto qi             = m_fields.at("qi"            ).get_view<const Real**> ();
  const auto qr             = m_fields.at("qr"            ).get_view<const Real**> ();
  const auto ps             = m_fields.at("ps"            ).get_view<const Real*>  ();
  const auto phis           = m_fields.at("phis"          ).get_view<const Real*>  ();

  const auto vapor_flux = m_fields.at("vapor_flux").get_view<const Real*>();
  const auto water_flux = m_fields.at("water_flux").get_view<const Real*>();
  const auto ice_flux   = m_fields.at("ice_flux"  ).get_view<const Real*>();
  const auto heat_flux  = m_fields.at("heat_flux" ).get_view<const Real*>();

  const auto policy = TPF::get_default_team_policy(ncols, nlevs);

  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto field_view_s1 = field_version_s1.get_view<Real*>();

  auto area = m_grid->get_geometry_data("area").clone();
  auto area_view = area.get_view<const Real*>();

  //reprosum vars
  const Real* send;
  Real recv;
  Int nlocal = ncols;
  Int ncount = 1;

//unite 1st and 2nd ||4 and repro calls  
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();
    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    // Calculate total gas mass (sum dp, no water loading)
    field_view_s1(i) = compute_gas_mass_on_column(team, nlevs, pseudo_density_i) * area_view(i);
  });
  Kokkos::fence();

  auto s1_host = field_version_s1.get_view<const Real*,Host>();
  send = s1_host.data();

  field_version_s1.sync_to_host(); 
  eamxx_repro_sum(send, &recv, nlocal, ncount, MPI_Comm_c2f(m_comm.mpi_comm()));
  field_version_s1.sync_to_dev();

  m_total_gas_mass_after = recv;

//this ||4 needs to be 2 for-loops, with one summing each 4 cols first, serially
//for pg2 grids (if on np4 grids, it would require much more work?)  
  auto energy_change = m_energy_change;
  auto current_energy = m_current_energy;
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {

    const int i = team.league_rank();
    const auto pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto T_mid_i          = ekat::subview(T_mid, i);
    const auto horiz_winds_i    = ekat::subview(horiz_winds, i);
    const auto qv_i             = ekat::subview(qv, i);
    const auto qc_i             = ekat::subview(qc, i);
    const auto qr_i             = ekat::subview(qr, i);
    const auto qi_i             = ekat::subview(qi, i);

    // Calculate total energy
    const auto new_energy_for_fixer = compute_total_energy_on_column(team, nlevs, pseudo_density_i, T_mid_i, horiz_winds_i,
                                                   qv_i, qc_i, qr_i, ps(i), phis(i));
    Kokkos::single(Kokkos::PerTeam(team),[&]() {
      energy_change(i) = compute_energy_boundary_flux_on_column(vapor_flux(i), water_flux(i), ice_flux(i), heat_flux(i))*dt;
      field_view_s1(i) = (current_energy(i)-new_energy_for_fixer-energy_change(i)) * area_view(i);
    });
  });
  Kokkos::fence();

  field_version_s1.sync_to_host(); 
  eamxx_repro_sum(send, &recv, nlocal, ncount, MPI_Comm_c2f(m_comm.mpi_comm()));
  field_version_s1.sync_to_dev();
  m_pb_fixer = recv;

  if(print_debug_info) {
    //total energy needed for relative error
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
      const int i = team.league_rank();
      Kokkos::single(Kokkos::PerTeam(team),[&]() {
        field_view_s1(i) = current_energy(i) * area_view(i);
      });
    });
    Kokkos::fence();

    field_version_s1.sync_to_host(); 
    eamxx_repro_sum(send, &recv, nlocal, ncount, MPI_Comm_c2f(m_comm.mpi_comm()));
    field_version_s1.sync_to_dev();

    m_total_energy_before = recv;
  }

  using PC = scream::physics::Constants<Real>;
  const Real cpdry = PC::Cpair;
  m_pb_fixer /= (cpdry*m_total_gas_mass_after); // T change due to fixer

  //add the fixer to temperature
  const auto pb_fixer=m_pb_fixer;
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
    const int i = team.league_rank();
    const auto T_mid_i = ekat::subview(T_mid, i);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevs), [&] (const int k){ 
       T_mid_i(k) += pb_fixer;
    });
  });//adding fix to T

  if(print_debug_info){
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA (const KT::MemberType& team) {

      const int i = team.league_rank();
      const auto pseudo_density_i = ekat::subview(pseudo_density, i);
      const auto T_mid_i          = ekat::subview(T_mid, i);
      const auto horiz_winds_i    = ekat::subview(horiz_winds, i);
      const auto qv_i             = ekat::subview(qv, i);
      const auto qc_i             = ekat::subview(qc, i);
      const auto qr_i             = ekat::subview(qr, i);
      const auto qi_i             = ekat::subview(qi, i);

      // Calculate total energy
      const auto new_energy_for_fixer = compute_total_energy_on_column(team, nlevs, pseudo_density_i, T_mid_i, horiz_winds_i,
                                                   qv_i, qc_i, qr_i, ps(i), phis(i));
      //overwrite the "new" fields with relative change

      Kokkos::single(Kokkos::PerTeam(team),[&]() {
        field_view_s1(i) = (current_energy(i)-new_energy_for_fixer-energy_change(i))*area_view(i);
      });
    });
    Kokkos::fence();

    field_version_s1.sync_to_host(); 
    eamxx_repro_sum(send, &recv, nlocal, ncount, MPI_Comm_c2f(m_comm.mpi_comm()));
    field_version_s1.sync_to_dev();

    m_echeck = recv/m_total_energy_before;
  }

};//global_fixer


KOKKOS_INLINE_FUNCTION
Real MassAndEnergyConservationCheck::
compute_total_mass_on_column (const KT::MemberType&       team,
                              const int                   nlevs,
                              const uview_1d<const Real>& pseudo_density,
                              const uview_1d<const Real>& qv,
                              const uview_1d<const Real>& qc,
                              const uview_1d<const Real>& qi,
                              const uview_1d<const Real>& qr)
{
  using PC = scream::physics::Constants<Real>;
  using RU = ekat::ReductionUtils<DefaultDevice::execution_space>;

  const Real gravit = PC::gravit;

  return RU::parallel_reduce<Real>(team, 0, nlevs,
                                   [&] (const int lev, Real& local_mass) {
    local_mass += (qv(lev)+
                   qc(lev)+
                   qi(lev)+
                   qr(lev))*pseudo_density(lev)/gravit;
  });
}

KOKKOS_INLINE_FUNCTION
Real MassAndEnergyConservationCheck::
compute_gas_mass_on_column (const KT::MemberType&       team,
                            const int                   nlevs,
                            const uview_1d<const Real>& pseudo_density)
{
  using RU = ekat::ReductionUtils<DefaultDevice::execution_space>;
  using PC = scream::physics::Constants<Real>;
  const Real gravit = PC::gravit;

  return RU::parallel_reduce<Real>(team, 0, nlevs,
                                              [&] (const int lev, Real& local_mass) {
    local_mass += pseudo_density(lev)/gravit;
  });
}

KOKKOS_INLINE_FUNCTION
Real MassAndEnergyConservationCheck::
compute_mass_boundary_flux_on_column (const Real vapor_flux,
                                      const Real water_flux)
{
  using PC = scream::physics::Constants<Real>;
  const Real RHO_H2O  = PC::RHO_H2O;

  return vapor_flux - water_flux*RHO_H2O;
}

KOKKOS_INLINE_FUNCTION
Real MassAndEnergyConservationCheck::
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
  using RU = ekat::ReductionUtils<DefaultDevice::execution_space>;

  const Real LatVap = PC::LatVap;
  const Real LatIce = PC::LatIce;
  const Real gravit = PC::gravit;
  const Real Cpair  = PC::Cpair;

  Real total_energy =
    RU::parallel_reduce<Real>(team, 0, nlevs,
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
Real MassAndEnergyConservationCheck::
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

Real MassAndEnergyConservationCheck::get_echeck() const{
  return m_echeck;
}

Real MassAndEnergyConservationCheck::get_total_energy_before() const{
  return m_total_energy_before;
}

Real MassAndEnergyConservationCheck::get_pb_fixer() const{
  return m_pb_fixer;
}


} // namespace scream
