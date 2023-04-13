#include <physics/mam/eamxx_mam_microphysics.hpp>
#include <share/property_checks/field_lower_bound_check.hpp>
#include <share/property_checks/field_within_interval_check.hpp>

#include "scream_config.h" // for SCREAM_CIME_BUILD

#include <ekat/ekat_assert.hpp>

namespace scream
{

MAMMicrophysics::MAMMicrophysics(
    const ekat::Comm& comm,
    const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params),
    nucleation_(nullptr) {
}

AtmosphereProcessType MAMMicrophysics::type() const {
  return AtmosphereProcessType::Physics;
}

std::string MAMMicrophysics::name() const {
  return "MAMMicrophysics";
}

void MAMMicrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  // The units of mixing ratio q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto q_unit = kg/kg;
  q_unit.set_string("kg/kg");
  Units nondim(0,0,0,0,0,0,0);

  grid_ = grids_manager->get_grid("Physics");
  const auto& grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs(); // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels(); // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout_col{ {COL}, {ncol_} };

  // Layout for 3D (2d horiz X 1d vertical) variables
  FieldLayout scalar3d_layout_mid{ {COL, LEV}, {ncol_, nlev_} };

  // Define fields needed in mam4xx.
  const auto m2 = m*m;
  const auto s2 = s*s;

  // atmospheric quantities
  add_field<Required>("T_mid", scalar3d_layout_mid, K, grid_name);
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name); // total pressure
  add_field<Required>("qv", scalar3d_layout_mid, q_unit, grid_name, "tracers");
  add_field<Required>("pbl_height", scalar2d_layout_col, m, grid_name);
  add_field<Required>("pseudo_density", scalar3d_layout_mid, q_unit, grid_name); // pdel
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name);

  // tracer group (stores all aerosol prognostics)
  add_group<Updated>("tracers", grid_name, Bundling::Required);
}

void MAMMicrophysics::
set_computed_group_impl(const FieldGroup& group) {
  const auto& name = group.m_info->m_group_name;
  EKAT_REQUIRE_MSG(name=="tracers",
    "Error! MAM4 expects a 'tracers' field group (got '" << name << "')\n");

  EKAT_REQUIRE_MSG(group.m_info->m_bundled,
    "Error! MAM4 expects bundled fields for tracers.\n");

  // How many aerosol/gas tracers do we expect? Recall that we maintain
  // both cloudborne and interstitial aerosol tracers.
  int num_aero_tracers =
    aero_config_.num_gas_ids() +  // gas tracers
    2 * aero_config_.num_modes(); // modal number mixing ratio tracers
  for (int m = 0; m < aero_config_.num_modes(); ++m) {
    auto m_index = static_cast<mam4::ModeIndex>(m);
    for (int a = 0; a < aero_config_.num_aerosol_ids(); ++a) {
      auto a_id = static_cast<mam4::AeroId>(a);
      if (mam4::aerosol_index_for_mode(m_index, a_id) != -1) {
        num_aero_tracers += 2; // aerosol mass mixing ratios (interstitial, cloudborne)
      }
    }
  }

  EKAT_REQUIRE_MSG(group.m_info->size() >= num_aero_tracers,
    "Error! MAM4 requires at least " << num_aero_tracers << " aerosol tracers.");
}

void MAMMicrophysics::initialize_impl(const RunType run_type) {
  const auto& T_mid = get_field_in("T_mid").get_view<Real**>();
  const auto& p_mid = get_field_in("p_mid").get_view<const Real**>();
  const auto& qv = get_field_in("qv").get_view<Real**>();
  const auto& pblh = get_field_in("pbl_height").get_view<Real*>();
  const auto& p_del = get_field_in("pseudo_density").get_view<const Real**>();
  const auto& cldfrac = get_field_in("cldfrac_tot").get_view<Real**>(); // FIXME: tot or liq?
  const auto& tracers = get_group_out("tracers");
  const auto& tracers_info = tracers.m_info;

  // Alias local variables from temporary buffer
  // e.g. auto z_mid       = buffer_.z_mid;

  // Perform any initialization work.
  if (run_type==RunType::Initial){
    /* e.g.
    Kokkos::deep_copy(sgs_buoy_flux,0.0);
    Kokkos::deep_copy(tk,0.0);
    Kokkos::deep_copy(tke,0.0004);
    Kokkos::deep_copy(tke_copy,0.0004);
    Kokkos::deep_copy(cldfrac_liq,0.0);
    */
  }

  // Find indices of qv and aerosol-related quantities
  auto qv_index  = tracers_info->m_subview_idx.at("qv");
  // FIXME

  //preprocess_.set_variables(ncol_, nlev_, T_mid, p_mid, qv, height,
  //                          p_del, pblh, q_soag, q_h2so4, q_nh3, q_aitken_so4);

  // FIXME: here we run the aerosol microphysics parameterizations

  //postprocess_.set_variables(ncol_, nlev_, qv, q_soag, q_h2so4,
  //                           q_nh3, q_aitken_so4);

  // Set field property checks for the fields in this process
  /* e.g.
  using Interval = FieldWithinIntervalCheck;
  using LowerBound = FieldLowerBoundCheck;
  add_postcondition_check<Interval>(get_field_out("T_mid"),m_grid,130.0,500.0,false);
  add_postcondition_check<LowerBound>(get_field_out("pbl_height"),m_grid,0);
  add_postcondition_check<Interval>(get_field_out("cldfrac_liq"),m_grid,0.0,1.0,false);
  add_postcondition_check<LowerBound>(get_field_out("tke"),m_grid,0);
  */

  // Setup WSM for internal local variables
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);
  // FIXME
  //workspace_mgr_.setup(buffer_.wsm_data, nlev_+1, 13+(n_wind_slots+n_trac_slots), default_policy);

  // FIXME: aerosol process initialization goes here!
}

void MAMMicrophysics::run_impl(const double dt) {

  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);

  // preprocess input
  Kokkos::parallel_for("preprocess", default_policy, preprocess_);
  Kokkos::fence();

  // Reset internal WSM variables.
  //workspace_mgr_.reset_internals();

  // FIXME: Aerosol stuff goes here!

  // postprocess output
  Kokkos::parallel_for("postprocess", default_policy, postprocess_);
  Kokkos::fence();
}

void MAMMicrophysics::finalize_impl()
{
}

} // namespace scream
