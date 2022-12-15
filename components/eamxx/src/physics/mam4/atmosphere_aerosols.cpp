#include "ekat/ekat_assert.hpp"
#include "physics/shoc/atmosphere_macrophysics.hpp"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "scream_config.h" // for SCREAM_CIME_BUILD

namespace scream
{

MAM4Aerosols::MAM4Aerosols(const ekat::Comm& comm,
                           const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params) {
}

MAM4Aerosols::AtmosphereProcessType type() const {
  return AtmosphereProcessType::Physics;
}

std::string MAM4Aerosols::name() const {
  return "MAM4Aerosols";
}

void MAM4Aerosols::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  // The units of mixing ratio q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto q_unit = kg/kg;
  q_unit.set_string("kg/kg");
  Units nondim(0,0,0,0,0,0,0);

  grid_ = grids_manager->get_grid("Physics");
  const auto& grid_name = grid_->name();

  num_cols_ = grid_->get_num_local_dofs(); // Number of columns on this rank
  num_levs_ = grid_->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout_col{ {COL}, {num_cols_} };

  // Layout for surf_mom_flux
  FieldLayout  surf_mom_flux_layout { {COL, CMP}, {num_cols_, 2} };

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {num_cols_,num_levs_} };
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {num_cols_,num_levs_+1} };

  // Layout for horiz_wind field
  FieldLayout horiz_wind_layout { {COL,CMP,LEV}, {num_cols_,2,num_levs_} };

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

void MAM4Aerosols::
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

size_t MAM4Aerosols::requested_buffer_size_in_bytes() const override {
  // Number of Reals needed by local views in the interface
  const size_t interface_request = Buffer::num_1d_scalar_ncol*num_cols_*sizeof(Real) +
                                   Buffer::num_1d_scalar_nlev*nlev_packs*sizeof(Spack) +
                                   Buffer::num_2d_vector_mid*num_cols_*nlev_packs*sizeof(Spack) +
                                   Buffer::num_2d_vector_int*num_cols_*nlevi_packs*sizeof(Spack) +
                                   Buffer::num_2d_vector_tr*num_cols_*num_tracer_packs*sizeof(Spack);

  // Number of Reals needed by the WorkspaceManager passed to shoc_main
  const auto policy       = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(num_cols_, nlev_packs);
  const int n_wind_slots  = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots  = ekat::npack<Spack>(num_tracers_+3)*Spack::n;
  const size_t wsm_request= WSM::get_total_bytes_needed(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy);

  return interface_request + wsm_request;
}

void MAM4Aerosols::init_buffers(const ATMBufferManager &buffer_manager) override {
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d scalar views
  using scalar_view_t = decltype(m_buffer.cell_length);
  scalar_view_t* _1d_scalar_view_ptrs[Buffer::num_1d_scalar_ncol] =
    {&m_buffer.cell_length, &m_buffer.wpthlp_sfc, &m_buffer.wprtp_sfc, &m_buffer.upwp_sfc, &m_buffer.vpwp_sfc
    };
  for (int i = 0; i < Buffer::num_1d_scalar_ncol; ++i) {
    *_1d_scalar_view_ptrs[i] = scalar_view_t(mem, num_cols_);
    mem += _1d_scalar_view_ptrs[i]->size();
  }

  m_buffer.pref_mid = decltype(m_buffer.pref_mid)(s_mem, nlev_packs);
  s_mem += m_buffer.pref_mid.size();

  using spack_2d_view_t = decltype(m_buffer.z_mid);
  spack_2d_view_t* _2d_spack_mid_view_ptrs[Buffer::num_2d_vector_mid] = {
    &m_buffer.z_mid, &m_buffer.rrho, &m_buffer.thv, &m_buffer.dz, &m_buffer.zt_grid, &m_buffer.wm_zt,
    &m_buffer.inv_exner, &m_buffer.thlm, &m_buffer.qw, &m_buffer.dse, &m_buffer.tke_copy, &m_buffer.qc_copy,
    &m_buffer.shoc_ql2, &m_buffer.shoc_mix, &m_buffer.isotropy, &m_buffer.w_sec, &m_buffer.wqls_sec, &m_buffer.brunt
  };

  spack_2d_view_t* _2d_spack_int_view_ptrs[Buffer::num_2d_vector_int] = {
    &m_buffer.z_int, &m_buffer.rrho_i, &m_buffer.zi_grid, &m_buffer.thl_sec, &m_buffer.qw_sec,
    &m_buffer.qwthl_sec, &m_buffer.wthl_sec, &m_buffer.wqw_sec, &m_buffer.wtke_sec, &m_buffer.uw_sec,
    &m_buffer.vw_sec, &m_buffer.w3
  };

  for (int i = 0; i < Buffer::num_2d_vector_mid; ++i) {
    *_2d_spack_mid_view_ptrs[i] = spack_2d_view_t(s_mem, num_cols_, nlev_packs);
    s_mem += _2d_spack_mid_view_ptrs[i]->size();
  }

  for (int i = 0; i < Buffer::num_2d_vector_int; ++i) {
    *_2d_spack_int_view_ptrs[i] = spack_2d_view_t(s_mem, num_cols_, nlevi_packs);
    s_mem += _2d_spack_int_view_ptrs[i]->size();
  }
  m_buffer.wtracer_sfc = decltype(m_buffer.wtracer_sfc)(s_mem, num_cols_, num_tracer_packs);
  s_mem += m_buffer.wtracer_sfc.size();

  // WSM data
  m_buffer.wsm_data = s_mem;

  // Compute workspace manager size to check used memory
  // vs. requested memory
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(num_cols_, nlev_packs);
  const int wsm_size = WSM::get_total_bytes_needed(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy)/sizeof(Spack);
  s_mem += wsm_size;

  size_t used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for SHOCMacrophysics.");
}

void MAM4Aerosols::initialize_impl(const RunType run_type) {
  const auto& T_mid = get_field_in("T_mid").get_view<Spack**>();
  const auto& p_mid = get_field_in("p_mid").get_view<const Spack**>();
  const auto& qv = get_field_in("qv").get_view<Spack**>();
  const auto& pblh = get_field_in("pbl_height").get_view<Real*>();
  const auto& p_del = get_field_in("pseudo_density").get_view<const Spack**>();
  const auto& cldfrac = get_field_in("cldfrac_tot").get_view<Spack**>(); // FIXME: tot or liq?
  const auto& tracers = get_group_out("tracers");
  const auto& tracers_info = tracers.m_info;

  // Alias local variables from temporary buffer
  // e.g. auto z_mid       = m_buffer.z_mid;

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

  mam4_preprocess.set_variables(num_cols_, num_levs_, T_mid, p_mid, qv, height,
                                pdel, pblh, q_soag, q_h2so4, q_nh3, q_aitken_so4);

  // input
  input.pres        = p_mid;
  input.pdel        = pseudo_density;
  input.temp        = mam4_preprocess.t_mid;

  // input/output
  // e.g. input_output.host_dse     = mam4_preprocess.shoc_s;

  // output (prognostic)
  // e.g. output.pblh     = get_field_out("pbl_height").get_view<Real*>();

  // output (diagnostic)
  // e.g. history_output.shoc_mix  = m_buffer.shoc_mix;

  mam4_postprocess.set_variables(num_cols_, num_levs_, qv, q_soag, q_h2so4,
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
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(num_cols_, nlev_);
  // FIXME
  //workspace_mgr.setup(m_buffer.wsm_data, nlev_+1, 13+(n_wind_slots+n_trac_slots), default_policy);

  // FIXME: aerosol process initialization goes here!
}

void MAM4Aerosols::run_impl(const int dt) {

  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(num_cols_, nlev_);

  // preprocess input
  Kokkos::parallel_for("mam4_preprocess", default_policy, mam4_preprocess);
  Kokkos::fence();

  // For now set the host timestep to the MAM4 timestep. This forces
  // number of MAM4 timesteps (nadv) to be 1.
  // TODO: input parameter?
  hdtime = dt;
  m_nadv = std::max(hdtime/dt,1);

  // Reset internal WSM variables.
  workspace_mgr.reset_internals();

  // FIXME: Aerosol stuff goes here!

  // postprocess output
  Kokkos::parallel_for("mam4_postprocess", default_policy, mam4_postprocess);
  Kokkos::fence();
}

void MAM4Aerosols::finalize_impl()
{
}

} // namespace scream
