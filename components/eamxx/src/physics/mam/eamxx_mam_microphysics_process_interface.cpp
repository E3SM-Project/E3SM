#include <physics/mam/eamxx_mam_microphysics_process_interface.hpp>
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
  return "mam4_micro";
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

  // atmospheric quantities
  add_field<Required>("omega", scalar3d_layout_mid, Pa/s, grid_name); // vertical pressure velocity
  add_field<Required>("T_mid", scalar3d_layout_mid, K, grid_name); // Temperature
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name); // total pressure
  add_field<Required>("qv", scalar3d_layout_mid, q_unit, grid_name, "tracers"); // specific humidity
  add_field<Required>("qi", scalar3d_layout_mid, q_unit, grid_name, "tracers"); // ice wet mixing ratio
  add_field<Required>("ni", scalar3d_layout_mid, n_unit, grid_name, "tracers"); // ice number mixing ratio
  add_field<Required>("pbl_height", scalar2d_layout_col, m, grid_name); // planetary boundary layer height
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid_name); // pdel, hydrostatic pressure
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name); // cloud fraction

  // droplet activation can alter cloud liquid and number mixing ratios
  add_field<Updated>("qc", scalar3d_layout_mid, q_unit, grid_name, "tracers"); // cloud liquid wet mixing ratio
  add_field<Updated>("nc", scalar3d_layout_mid, n_unit, grid_name, "tracers"); // cloud liquid wet number mixing ratio

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing ratios
  for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    const char* int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    add_field<Updated>(int_nmr_field_name, scalar3d_layout_mid, n_unit, grid_name, "tracers");
    for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char* int_mmr_field_name = mam_coupling::int_aero_mmr_field_name(m, a);
      if (strlen(int_mmr_field_name) > 0) {
        add_field<Updated>(int_mmr_field_name, scalar3d_layout_mid, q_unit, grid_name, "tracers");
      }
    }
  }

  // aerosol-related gases: mass mixing ratios
  for (int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char* gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_field<Updated>(gas_mmr_field_name, scalar3d_layout_mid, q_unit, grid_name, "tracers");
  }

  // Tracers group -- do we need this in addition to the tracers above? In any
  // case, this call should be idempotent, so it can't hurt.
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

  // How many aerosol/gas tracers do we expect?
  int num_tracers = mam_coupling::num_aero_modes() +
                2 * mam_coupling::num_aero_tracers() + // interstitial + cloudborne aerosol mmrs
                    mam_coupling::num_aero_gases();
  EKAT_REQUIRE_MSG(group.m_info->size() >= num_tracers,
    "Error! MAM4 requires at least " << num_tracers << " aerosol tracers.");
}

size_t MAMMicrophysics::requested_buffer_size_in_bytes() const
{
  // number of Reals needed by local views in interface
  const size_t request = sizeof(Real) *
    (mam_coupling::Buffer::num_2d_mid * ncol_ * nlev_ +
     mam_coupling::Buffer::num_2d_iface * ncol_ * (nlev_+1));

  // FIXME: Need to figure out whether we need this stuff
  /*
  // Number of Reals needed by the WorkspaceManager passed to shoc_main
  const auto policy       = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);
  const int n_wind_slots  = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots  = ekat::npack<Spack>(m_num_tracers+3)*Spack::n;
  const size_t wsm_request= WSM::get_total_bytes_needed(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy);
  */

  return request;// + wsm_request;
}

void MAMMicrophysics::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
                   "Error! Insufficient buffer size.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // set view pointers for midpoint fields
  using view_2d_t = decltype(buffer_.z_mid);
  view_2d_t* view_2d_mid_ptrs[mam_coupling::Buffer::num_2d_mid] = {
    &buffer_.z_mid,
    &buffer_.dz,
    &buffer_.qv_dry,
    &buffer_.qc_dry,
    &buffer_.nc_dry,
    &buffer_.qi_dry,
    &buffer_.ni_dry,
    &buffer_.w_updraft,
  };
  for (int i = 0; i < mam_coupling::Buffer::num_2d_mid; ++i) {
    *view_2d_mid_ptrs[i] = view_2d_t(mem, ncol_, nlev_);
    mem += view_2d_mid_ptrs[i]->size();
  }

  // set view pointers for interface fields
  view_2d_t* view_2d_iface_ptrs[mam_coupling::Buffer::num_2d_iface] = {&buffer_.z_iface};
  for (int i = 0; i < mam_coupling::Buffer::num_2d_iface; ++i) {
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

  const auto& T_mid = get_field_in("T_mid").get_view<const Real**>();
  const auto& p_mid = get_field_in("p_mid").get_view<const Real**>();
  const auto& qv_wet = get_field_in("qv").get_view<const Real**>();
  const auto& pblh = get_field_in("pbl_height").get_view<const Real*>();
  const auto& p_del = get_field_in("pseudo_density").get_view<const Real**>();
  const auto& cldfrac = get_field_in("cldfrac_tot").get_view<const Real**>(); // FIXME: tot or liq?
  const auto& qc_wet = get_field_out("qc").get_view<Real**>();
  const auto& nc_wet = get_field_out("nc").get_view<Real**>();
  const auto& qi_wet = get_field_in("qi").get_view<const Real**>();
  const auto& ni_wet = get_field_in("ni").get_view<const Real**>();
  const auto& omega = get_field_in("omega").get_view<const Real**>();

  const auto& tracers = get_group_out("tracers");
  const auto& tracers_info = tracers.m_info;

  // Alias local variables from temporary buffer
  auto z_mid = buffer_.z_mid;
  auto dz    = buffer_.dz;
  auto z_iface = buffer_.z_iface;
  auto qv_dry = buffer_.qv_dry;
  auto qc_dry = buffer_.qc_dry;
  auto nc_dry = buffer_.nc_dry;
  auto qi_dry = buffer_.qi_dry;
  auto ni_dry = buffer_.ni_dry;
  auto w_updraft = buffer_.w_updraft;

  // Perform any initialization work.
  if (run_type==RunType::Initial){
  }

  // set atmosphere state data
  atm_.T_mid = T_mid;
  atm_.p_mid = p_mid;
  atm_.qv_wet = qv_wet;
  atm_.qv_dry = qv_dry;
  atm_.qc_wet = qc_wet;
  atm_.qc_dry = qc_wet;
  atm_.nc_wet = nc_wet;
  atm_.nc_dry = nc_dry;
  atm_.qi_wet = qi_wet;
  atm_.qi_dry = qi_dry;
  atm_.ni_wet = ni_wet;
  atm_.ni_dry = ni_dry;
  atm_.pdel = p_del;
  atm_.cldfrac = cldfrac;
  atm_.pblh = pblh;

  // set aerosol state data (interstitial aerosols only)
  for (int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    const char* int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    aero_.wet_int_aero_nmr[m] = get_field_out(int_nmr_field_name).get_view<Real**>();
    aero_.dry_int_aero_nmr[m] = buffer_.dry_int_aero_nmr[m];
    for (int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      const char* int_mmr_field_name = mam_coupling::int_aero_mmr_field_name(m, a);
      if (strlen(int_mmr_field_name) > 0) {
        aero_.wet_int_aero_mmr[m][a] = get_field_out(int_mmr_field_name).get_view<Real**>();
        aero_.dry_int_aero_mmr[m][a] = buffer_.dry_int_aero_mmr[m][a];
      }
    }
  }

  // set aerosol-related gas state data
  for (int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char* mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    aero_.wet_gas_mmr[g] = get_field_out(mmr_field_name).get_view<Real**>();
    aero_.dry_gas_mmr[g] = buffer_.dry_gas_mmr[g];
  }

  // FIXME: For now, set z_surf to zero.
  const Real z_surf = 0.0;

  // set up our preprocess/postprocess functors
  preprocess_.initialize(ncol_, nlev_, z_surf, atm_, aero_);
  postprocess_.initialize(ncol_, nlev_, atm_, aero_);

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

  // configure the nucleation parameterization
  mam4::NucleationProcess::ProcessConfig nuc_config{};
  nuc_config.dens_so4a_host = 1770.0;
  nuc_config.mw_so4a_host = 115.0;
  nuc_config.newnuc_method_user_choice = 2;
  nuc_config.pbl_nuc_wang2008_user_choice = 1;
  nuc_config.adjust_factor_pbl_ratenucl = 1.0;
  nuc_config.accom_coef_h2so4 = 1.0;
  nuc_config.newnuc_adjust_factor_dnaitdt = 1.0;
  nucleation_->init(nuc_config);
}

void MAMMicrophysics::run_impl(const double dt) {

  const auto scan_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);
  const auto policy      = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  // Reset internal WSM variables.
  //workspace_mgr_.reset_internals();

  // FIXME: nothing depends on simulation time (yet), so we can just use zero for now
  double t = 0.0;

  // Compute nucleation tendencies on all local columns and accumulate them
  // into our tracer state.
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const ThreadTeam& team) {
    const Int icol = team.league_rank(); // column index

    // extract column-specific atmosphere state data
    auto atm = atmosphere_for_column(atm_, icol);

    // set surface state data
    haero::Surface sfc{};

    // extract column-specific subviews into aerosol prognostics
    auto progs = mam_coupling::interstitial_aerosols_for_column(aero_, icol);

    // set up diagnostics
    mam4::Diagnostics diags(nlev_);

    // grab views from the buffer to store tendencies
    mam4::Tendencies tends(nlev_);
    /*
    tends.q_gas[ih2so4] = ekat::subview(buffer_.q_h2so4_tend, icol);
    tends.n_mode_i[iait] = ekat::subview(buffer_.n_aitken_tend, icol);
    tends.q_aero_i[iait][iso4] = ekat::subview(buffer_.q_aitken_so4_tend, icol);
    */

    // run the nucleation process to obtain tendencies
    nucleation_->compute_tendencies(team, t, dt, atm, sfc, progs, diags, tends);
    /*
#ifndef NDEBUG
    const int lev_idx = 0;
    if (icol == 0) {
    m_atm_logger->debug("tends.q_gas[ih2so4] = {}, tends.n_mode_i[iait] = {}, tends.q_aero_i[iait][iso4] = {}",
      tends.q_gas[ih2so4](lev_idx), tends.n_mode_i[iait](lev_idx), tends.q_aero_i[iait][iso4](lev_idx));
    }
#endif

    // accumulate tendencies into prognostics
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_), [&](const int klev) {
      progs.q_gas[ih2so4](klev) += dt * tends.q_gas[ih2so4](klev);
      progs.n_mode_i[iait](klev) += dt * tends.n_mode_i[iait](klev);
      progs.q_aero_i[iait][iso4](klev) += dt * tends.q_aero_i[iait][iso4](klev);
    });
    */
  });

  // postprocess output
  Kokkos::parallel_for("postprocess", policy, postprocess_);
  Kokkos::fence();
}

void MAMMicrophysics::finalize_impl()
{
}

} // namespace scream
