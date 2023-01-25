#include "ekat/ekat_assert.hpp"
#include "physics/mam4/aerosol_microphysics.hpp"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "scream_config.h" // for SCREAM_CIME_BUILD

namespace scream
{

MAM4AerosolMicrophysics::MAM4AerosolMicrophysics(
    const ekat::Comm& comm,
    const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params),
    nucleation_(nullptr) {
}

AtmosphereProcessType type() const {
  return AtmosphereProcessType::Physics;
}

std::string MAM4AerosolMicrophysics::name() const {
  return "MAM4AerosolMicrophysics";
}

void MAM4AerosolMicrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  // The units of mixing ratio q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto q_unit = kg/kg;
  q_unit.set_string("kg/kg");
  Units nondim(0,0,0,0,0,0,0);

  grid_ = grids_manager->get_grid("Physics");
  const auto& grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs(); // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout_col{ {COL}, {ncol_} };

  // Define fields needed in mam4xx.
  const auto m2 = m*m;
  const auto s2 = s*s;

  // atmospheric quantities
  add_field<Required>("T_mid", scalar3d_layout_mid, K, grid_name, ps);
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name, ps); // total pressure
  add_field<Required>("qv", scalar3d_layout_mid, Qunit, grid_name, "tracers", ps);
  add_field<Required>("pbl_height", scalar2d_layout_col, m, grid_name);
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Qunit, grid_name, ps); // pdel
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name, ps);

  // tracer group (stores all aerosol prognostics)
  add_group<Updated>("tracers",grid_name,ps,Bundling::Required);
}

void MAM4AerosolMicrophysics::
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
  for (int m = 0; m < aero_config_.num_modes(); ++m) { // aerosol tracers
    for (int a = 0; a < aero_config_.num_aerosol_ids(); ++a) {
      if (aerosol_index_for_mode(m, a) != -1) {
        num_aero_tracers += 2;
      }
    }
  }

  EKAT_REQUIRE_MSG(group.m_info->size() >= num_aero_tracers,
    "Error! MAM4 requires at least " << num_aero_tracers << " aerosol tracers.");
}

size_t MAM4AerosolMicrophysics::requested_buffer_size_in_bytes() const override {
  // Number of Reals needed by local views in the interface
  const size_t interface_request = Buffer::num_1d_scalar_ncol*ncol_*sizeof(Real) +
                                   Buffer::num_1d_scalar_nlev*nlev_packs*sizeof(Real) +
                                   Buffer::num_2d_vector_mid*ncol_*nlev_packs*sizeof(Real) +
                                   Buffer::num_2d_vector_int*ncol_*nlevi_packs*sizeof(Real) +
                                   Buffer::num_2d_vector_tr*ncol_*num_tracer_packs*sizeof(Real);

  // Number of Reals needed by the WorkspaceManager passed to aerosol processes
//  const auto policy       = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_packs);
//  const size_t wsm_request= WSM::get_total_bytes_needed(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy);

  return interface_request;// + wsm_request;
}

void MAM4AerosolMicrophysics::init_buffers(const ATMBufferManager &buffer_manager) override {
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d scalar views
  using scalar_view_t = decltype(buffer_.cell_length);
  scalar_view_t* _1d_scalar_view_ptrs[Buffer::num_1d_scalar_ncol] =
    {&buffer_.cell_length, &buffer_.wpthlp_sfc, &buffer_.wprtp_sfc, &buffer_.upwp_sfc, &buffer_.vpwp_sfc
    };
  for (int i = 0; i < Buffer::num_1d_scalar_ncol; ++i) {
    *_1d_scalar_view_ptrs[i] = scalar_view_t(mem, ncol_);
    mem += _1d_scalar_view_ptrs[i]->size();
  }

  buffer_.pref_mid = decltype(buffer_.pref_mid)(s_mem, nlev_packs);
  s_mem += buffer_.pref_mid.size();

  using spack_2d_view_t = decltype(buffer_.z_mid);
  spack_2d_view_t* _2d_spack_mid_view_ptrs[Buffer::num_2d_vector_mid] = {
    &buffer_.z_mid, &buffer_.rrho, &buffer_.thv, &buffer_.dz, &buffer_.zt_grid, &buffer_.wm_zt,
    &buffer_.inv_exner, &buffer_.thlm, &buffer_.qw, &buffer_.dse, &buffer_.tke_copy, &buffer_.qc_copy,
    &buffer_.shoc_ql2, &buffer_.shoc_mix, &buffer_.isotropy, &buffer_.w_sec, &buffer_.wqls_sec, &buffer_.brunt
  };

  spack_2d_view_t* _2d_spack_int_view_ptrs[Buffer::num_2d_vector_int] = {
    &buffer_.z_int, &buffer_.rrho_i, &buffer_.zi_grid, &buffer_.thl_sec, &buffer_.qw_sec,
    &buffer_.qwthl_sec, &buffer_.wthl_sec, &buffer_.wqw_sec, &buffer_.wtke_sec, &buffer_.uw_sec,
    &buffer_.vw_sec, &buffer_.w3
  };

  for (int i = 0; i < Buffer::num_2d_vector_mid; ++i) {
    *_2d_spack_mid_view_ptrs[i] = spack_2d_view_t(s_mem, ncol_, nlev_packs);
    s_mem += _2d_spack_mid_view_ptrs[i]->size();
  }

  for (int i = 0; i < Buffer::num_2d_vector_int; ++i) {
    *_2d_spack_int_view_ptrs[i] = spack_2d_view_t(s_mem, ncol_, nlevi_packs);
    s_mem += _2d_spack_int_view_ptrs[i]->size();
  }
  buffer_.wtracer_sfc = decltype(buffer_.wtracer_sfc)(s_mem, ncol_, num_tracer_packs);
  s_mem += buffer_.wtracer_sfc.size();

  // WSM data
  buffer_.wsm_data = s_mem;

  // Compute workspace manager size to check used memory
  // vs. requested memory
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_packs);
  const int wsm_size = WSM::get_total_bytes_needed(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy)/sizeof(Spack);
  s_mem += wsm_size;

  size_t used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for SHOCMacrophysics.");
}

void MAM4AerosolMicrophysics::initialize_impl(const RunType run_type) {
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
  auto qv_index  = tracer_info->m_subview_idx.at("qv");
  // FIXME

  preprocess_.set_variables(ncol_, nlev_, T_mid, p_mid, qv, height,
                            p_del, pblh, q_soag, q_h2so4, q_nh3, q_aitken_so4);

  // FIXME: stuff goes here!

  postprocess_.set_variables(ncol_, nlev_, qv, q_soag, q_h2so4,
                             q_nh3, q_aitken_so4);

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

void MAM4AerosolMicrophysics::run_impl(const int dt) {

  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);

  // preprocess input
  Kokkos::parallel_for("preprocess", default_policy, preprocess_);
  Kokkos::fence();

  // Reset internal WSM variables.
  workspace_mgr_.reset_internals();

  // FIXME: Aerosol stuff goes here!

  // postprocess output
  Kokkos::parallel_for("postprocess", default_policy, postprocess_);
  Kokkos::fence();
}

void MAM4AerosolMicrophysics::finalize_impl()
{
}

} // namespace scream
