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
    aero_config_(), nucleation_(new mam4::NucleationProcess(aero_config_)) {
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
  auto n_unit = 1/kg;
  n_unit.set_string("#/kg");
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

  // aerosol tracers of interest: mass (q) and number (n) mixing ratios
  add_field<Updated>("q_aitken_so4", scalar3d_layout_mid, q_unit, grid_name, "tracers");
  add_field<Updated>("n_aitken_so4", scalar3d_layout_mid, n_unit, grid_name, "tracers");

  // aerosol-related gases: mass mixing ratios
  add_field<Updated>("q_soag", scalar3d_layout_mid, q_unit, grid_name, "tracers");
  add_field<Updated>("q_h2so4", scalar3d_layout_mid, q_unit, grid_name, "tracers");

  // Tracer group -- do we need this in addition to the tracers above?
  add_group<Updated>("tracers", grid_name, 1, Bundling::Required);
}

// this checks whether we have the tracers we expect
void MAMMicrophysics::
set_computed_group_impl(const FieldGroup& group) {
  const auto& name = group.m_info->m_group_name;
  EKAT_REQUIRE_MSG(name=="tracers",
    "Error! MAM4 expects a 'tracers' field group (got '" << name << "')\n");

  EKAT_REQUIRE_MSG(group.m_info->m_bundled,
    "Error! MAM4 expects bundled fields for tracers.\n");

  // How many aerosol/gas tracers do we expect? Recall that we maintain
  // both cloudborne and interstitial aerosol tracers.
  int num_aero_tracers = 4; // for now, just 2 gases + aitken so4 n,q
  /*
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
  */

  EKAT_REQUIRE_MSG(group.m_info->size() >= num_aero_tracers,
    "Error! MAM4 requires at least " << num_aero_tracers << " aerosol tracers.");
}

size_t MAMMicrophysics::requested_buffer_size_in_bytes() const
{
  // number of Reals needed by local views in interface
  const size_t iface_request = sizeof(Real) *
    (Buffer::num_2d_mid * ncol_ * nlev_ +
     Buffer::num_2d_iface * ncol_ * (nlev_+1));

  // FIXME: Need to figure out whether we need this stuff
  /*
  // Number of Reals needed by the WorkspaceManager passed to shoc_main
  const auto policy       = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);
  const int n_wind_slots  = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots  = ekat::npack<Spack>(m_num_tracers+3)*Spack::n;
  const size_t wsm_request= WSM::get_total_bytes_needed(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy);
  */

  return iface_request;// + wsm_request;
}

void MAMMicrophysics::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error! Insufficient buffer size.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // set view pointers for midpoint fields
  using view_2d_t = decltype(buffer_.z_mid);
  view_2d_t* view_2d_mid_ptrs[Buffer::num_2d_mid] = {&buffer_.z_mid, &buffer_.dz};
  for (int i = 0; i < Buffer::num_2d_mid; ++i) {
    *view_2d_mid_ptrs[i] = view_2d_t(mem, ncol_, nlev_);
    mem += view_2d_mid_ptrs[i]->size();
  }

  // set view pointers for interface fields
  view_2d_t* view_2d_iface_ptrs[Buffer::num_2d_iface] = {&buffer_.z_iface};
  for (int i = 0; i < Buffer::num_2d_iface; ++i) {
    *view_2d_iface_ptrs[i] = view_2d_t(mem, ncol_, nlev_+1);
    mem += view_2d_iface_ptrs[i]->size();
  }

  // WSM data
  buffer_.wsm_data = mem;

  /* FIXME: this corresponds to the FIXME in the above function
  // Compute workspace manager size to check used memory
  // vs. requested memory
  const auto policy      = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);
  const int n_wind_slots = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots = ekat::npack<Spack>(m_num_tracers+3)*Spack::n;
  const int wsm_size     = WSM::get_total_bytes_needed(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy)/sizeof(Spack);
  mem += wsm_size;
  */

  size_t used_mem = (mem - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for MAMMicrophysics.");
}

void MAMMicrophysics::initialize_impl(const RunType run_type) {
  const auto& T_mid = get_field_in("T_mid").get_view<Real**>();
  const auto& p_mid = get_field_in("p_mid").get_view<Real**>();
  const auto& qv = get_field_in("qv").get_view<Real**>();
  const auto& pblh = get_field_in("pbl_height").get_view<Real*>();
  const auto& p_del = get_field_in("pseudo_density").get_view<Real**>();
  const auto& cldfrac = get_field_in("cldfrac_tot").get_view<Real**>(); // FIXME: tot or liq?

  const auto& tracers = get_group_out("tracers");
  const auto& tracers_info = tracers.m_info;
  int num_tracers = tracers_info->size();

  // Alias local variables from temporary buffer
  auto z_mid = buffer_.z_mid;
  auto dz    = buffer_.dz;
  auto z_iface = buffer_.z_iface;

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

  // set atmosphere state data
  T_mid_ = T_mid;
  p_mid_ = p_mid;
  qv_ = qv;
  pdel_ = p_del;
//  cloud_f_ = cloud_f; // cloud fraction
//  uv_ = uv; // updraft velocity

  // For now, set z_surf to zero.
  const Real z_surf = 0.0;

  // Determine indices of aerosol/gas tracers for wet<->dry conversion
  auto q_aitken_so4_index  = tracers_info->m_subview_idx.at("q_aitken_so4");
  auto q_soag_index        = tracers_info->m_subview_idx.at("q_soag");
  auto q_h2so4_index       = tracers_info->m_subview_idx.at("q_h2so4");
  int num_aero_tracers = 3; // for now, just 2 gases + aitken so4
  view_1d_int convert_wet_dry_idx_d("convert_wet_dry_idx_d", num_aero_tracers);
  auto convert_wet_dry_idx_h = Kokkos::create_mirror_view(convert_wet_dry_idx_d);
  for (int it=0, iq=0; it < num_tracers; ++it) {
    if ((it == q_aitken_so4_index) || (it == q_soag_index) || (it == q_h2so4_index)) {
      convert_wet_dry_idx_h(iq) = it;
      ++iq;
    }
  }
  Kokkos::deep_copy(convert_wet_dry_idx_d, convert_wet_dry_idx_h);

  preprocess_.set_variables(ncol_, nlev_, z_surf, convert_wet_dry_idx_d, T_mid,
                            p_mid, qv, z_mid, z_iface, dz, pdel_, pblh, q_soag_,
                            q_h2so4_, q_aitken_so4_);

  // FIXME: here we run the aerosol microphysics parameterizations

  //postprocess_.set_variables(ncol_, nlev_, qv, q_soag, q_h2so4, q_aitken_so4);

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
  // FIXME: do we need this?
  //const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);
  //workspace_mgr_.setup(buffer_.wsm_data, nlev_+1, 13+(n_wind_slots+n_trac_slots), default_policy);

  // FIXME: aerosol process initialization goes here!
}

void MAMMicrophysics::run_impl(const double dt) {

  const auto scan_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);
  const auto policy      = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  // Reset internal WSM variables.
  //workspace_mgr_.reset_internals();

  // nothing depends on simulation time (yet), so we can just use zero for now
  double t = 0.0;

  // FIXME: for now, we set cloud fraction and updraft velocity to zero.
  view_1d cloud_f("cloud fraction", nlev_);
  view_1d uv("updraft velocity", nlev_);
  Kokkos::deep_copy(cloud_f, 0.0);
  Kokkos::deep_copy(uv, 0.0);

  // Compute nucleation tendencies on all local columns and accumulate them
  // into our tracer state.
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ThreadTeam& team) {
    const Int icol = team.league_rank(); // column index

    // extract column-specific atmosphere state data
    haero::Atmosphere atm(nlev_, pblh_(icol));
    atm.temperature = ekat::subview(T_mid_, icol);
    atm.pressure = ekat::subview(p_mid_, icol);
    atm.vapor_mixing_ratio = ekat::subview(qv_, icol);
    atm.height = ekat::subview(height_, icol);
    atm.hydrostatic_dp = ekat::subview(pdel_, icol);
    atm.cloud_fraction = cloud_f;
    atm.updraft_vel_ice_nucleation = uv;

    // extract column-specific subviews into aerosol prognostics
    using AeroConfig = mam4::AeroConfig;
    using ModeIndex = mam4::ModeIndex;
    using AeroId = mam4::AeroId;
    using GasId = mam4::GasId;

    mam4::Prognostics progs(nlev_);
    int iait = static_cast<int>(ModeIndex::Aitken);
    progs.n_mode_i[iait] = ekat::subview(n_aitken_, icol);
    int iso4 = mam4::aerosol_index_for_mode(ModeIndex::Aitken, AeroId::SO4);
    progs.q_aero_i[iait][iso4] = ekat::subview(q_aitken_so4_, icol);
    int ih2so4 = static_cast<int>(GasId::H2SO4);
    progs.q_gas[ih2so4] = ekat::subview(q_h2so4_, icol);

    mam4::Diagnostics diags(nlev_);

    // run the nucleation process to obtain tendencies
    mam4::Tendencies tends(nlev_);
    nucleation_->compute_tendencies(team, t, dt, atm, progs, diags, tends);

    // accumulate tendencies into tracers state
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_), [&](const int klev) {
      progs.n_mode_i[iait](klev) += dt * tends.n_mode_i[iait](klev);
      progs.q_aero_i[iait][iso4](klev) += dt * tends.q_aero_i[iait][iso4](klev);
      progs.q_gas[ih2so4](klev) += dt * tends.q_gas[ih2so4](klev);
    });
  });

  // postprocess output
  Kokkos::parallel_for("postprocess", policy, postprocess_);
  Kokkos::fence();
}

void MAMMicrophysics::finalize_impl()
{
}

} // namespace scream
