#include "ekat/ekat_assert.hpp"
#include "physics/shoc/eamxx_shoc_process_interface.hpp"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "eamxx_config.h" // for SCREAM_CIME_BUILD

namespace scream
{

// =========================================================================================
SHOCMacrophysics::SHOCMacrophysics (const ekat::Comm& comm,const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  /* Anything that can be initialized without grid information can be initialized here.
   * Like universal constants, shoc options.
   */
}

// =========================================================================================
void SHOCMacrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();

  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d = m_grid->get_2d_scalar_layout();

  // Layout for surf_mom_flux
  FieldLayout vector2d = m_grid->get_2d_vector_layout(2);

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(true);
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(false);

  // Layout for horiz_wind field
  FieldLayout vector3d_mid = m_grid->get_3d_vector_layout(true,2);

  // Define fields needed in SHOC.
  // Note: shoc_main is organized by a set of 5 structures, variables below are organized
  //       using the same approach to make it easier to follow.

  constexpr int ps = Spack::n;

  const auto nondim = Units::nondimensional();
  const auto m2 = pow(m,2);
  const auto s2 = pow(s,2);

  // These variables are needed by the interface, but not actually passed to shoc_main.
  add_field<Required>("omega",          scalar3d_mid, Pa/s, grid_name, ps);
  add_field<Required>("surf_sens_flux", scalar2d    , W/m2, grid_name);
  add_field<Required>("surf_mom_flux",  vector2d    , N/m2, grid_name);

  add_field<Updated>("surf_evap",       scalar2d    , kg/(m2*s), grid_name);
  add_field<Updated> ("T_mid",          scalar3d_mid, K,         grid_name, ps);
  add_tracer<Updated>("qv", m_grid, kg/kg, ps);

  // If TMS is a process, add surface drag coefficient to required fields
  if (m_params.get<bool>("apply_tms", false)) {
    add_field<Required>("surf_drag_coeff_tms", scalar2d,  kg/(m2*s), grid_name);
  }

  // Input variables
  add_field<Required>("p_mid",          scalar3d_mid, Pa,    grid_name, ps);
  add_field<Required>("p_int",          scalar3d_int, Pa,    grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_mid, Pa,    grid_name, ps);
  add_field<Required>("phis",           scalar2d    , m2/s2, grid_name, ps);

  // Input/Output variables
  add_field<Updated>("horiz_winds",   vector3d_mid,   m/s,     grid_name, ps);
  add_field<Updated>("sgs_buoy_flux", scalar3d_mid, K*(m/s), grid_name, ps);
  add_field<Updated>("eddy_diff_mom", scalar3d_mid, m2/s,    grid_name, ps);
  add_field<Updated>("cldfrac_liq",   scalar3d_mid, nondim,  grid_name, ps);
  add_tracer<Updated>("tke", m_grid, m2/s2, ps);
  add_tracer<Updated>("qc",  m_grid, kg/kg, ps);

  // Output variables
  add_field<Computed>("pbl_height",    scalar2d    , m,            grid_name);
  add_field<Computed>("inv_qc_relvar", scalar3d_mid, pow(kg/kg,2), grid_name, ps);
  add_field<Computed>("eddy_diff_heat",   scalar3d_mid, m2/s,        grid_name, ps);
  add_field<Computed>("w_variance",       scalar3d_mid, m2/s2,       grid_name, ps);
  add_field<Computed>("cldfrac_liq_prev", scalar3d_mid, nondim,      grid_name, ps);
  add_field<Computed>("ustar",            scalar2d,     m/s,         grid_name, ps);
  add_field<Computed>("obklen",           scalar2d,     m,           grid_name, ps);

  // Extra SHOC output diagnostics
  if (m_params.get<bool>("extra_shoc_diags", false)) {

    // Diagnostic output - mid point grid
    add_field<Computed>("brunt", scalar3d_mid, pow(s,-1), grid_name, ps);
    add_field<Computed>("shoc_mix", scalar3d_mid, m, grid_name, ps);
    add_field<Computed>("isotropy", scalar3d_mid, s, grid_name, ps);

    // Diagnostic output - interface grid
    add_field<Computed>("wthl_sec", scalar3d_int, K*(m/s), grid_name, ps);
    add_field<Computed>("thl_sec", scalar3d_int, pow(K,2), grid_name, ps);
    add_field<Computed>("wqw_sec", scalar3d_int, (kg/kg)*(m/s), grid_name, ps);
    add_field<Computed>("qw_sec", scalar3d_int, pow(kg/kg,2), grid_name, ps);
    add_field<Computed>("uw_sec", scalar3d_int, pow(m/s,2), grid_name, ps);
    add_field<Computed>("vw_sec", scalar3d_int, pow(m/s,2), grid_name, ps);
    add_field<Computed>("w3", scalar3d_int, pow(m/s,3), grid_name, ps);

  } // Extra SHOC output diagnostics

  // Tracer group
  add_group<Updated>("tracers", grid_name, ps, Bundling::Required);

  // Boundary flux fields for energy and mass conservation checks
  if (has_column_conservation_check()) {
    add_field<Computed>("vapor_flux", scalar2d, kg/(m2*s), grid_name);
    add_field<Computed>("water_flux", scalar2d, m/s,     grid_name);
    add_field<Computed>("ice_flux",   scalar2d, m/s,     grid_name);
    add_field<Computed>("heat_flux",  scalar2d, W/m2,    grid_name);
  }
}

// =========================================================================================
void SHOCMacrophysics::
set_computed_group_impl (const FieldGroup& group)
{
  EKAT_REQUIRE_MSG(group.m_info->size() >= 3,
                   "Error! Shoc requires at least 3 tracers (tke, qv, qc) as inputs.");

  const auto& name = group.m_info->m_group_name;

  EKAT_REQUIRE_MSG(name=="tracers",
    "Error! We were not expecting a field group called '" << name << "\n");

  EKAT_REQUIRE_MSG(group.m_info->m_bundled,
      "Error! Shoc expects bundled fields for tracers.\n");

  // Calculate number of advected tracers
  m_num_tracers = group.m_info->size();
}

// =========================================================================================
size_t SHOCMacrophysics::requested_buffer_size_in_bytes() const
{
  const int nlev_packs       = ekat::npack<Spack>(m_num_levs);
  const int nlevi_packs      = ekat::npack<Spack>(m_num_levs+1);
  const int num_tracer_packs = ekat::npack<Spack>(m_num_tracers);

  // Number of Reals needed by local views in the interface
  const size_t interface_request = Buffer::num_1d_scalar_ncol*m_num_cols*sizeof(Real) +
                                   Buffer::num_1d_scalar_nlev*nlev_packs*sizeof(Spack) +
                                   Buffer::num_2d_vector_mid*m_num_cols*nlev_packs*sizeof(Spack) +
                                   Buffer::num_2d_vector_int*m_num_cols*nlevi_packs*sizeof(Spack) +
                                   Buffer::num_2d_vector_tr*m_num_cols*num_tracer_packs*sizeof(Spack);

  // Number of Reals needed by the WorkspaceManager passed to shoc_main
  const auto policy       = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);
  const int n_wind_slots  = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots  = ekat::npack<Spack>(m_num_tracers+3)*Spack::n;
  const size_t wsm_request= WSM::get_total_bytes_needed(nlevi_packs, 14+(n_wind_slots+n_trac_slots), policy);

  return interface_request + wsm_request;
}

// =========================================================================================
void SHOCMacrophysics::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d scalar views
  using scalar_view_t = decltype(m_buffer.wpthlp_sfc);
  scalar_view_t* _1d_scalar_view_ptrs[Buffer::num_1d_scalar_ncol] =
    {&m_buffer.wpthlp_sfc, &m_buffer.wprtp_sfc, &m_buffer.upwp_sfc, &m_buffer.vpwp_sfc
#ifdef SCREAM_SHOC_SMALL_KERNELS
     , &m_buffer.se_b, &m_buffer.ke_b, &m_buffer.wv_b, &m_buffer.wl_b
     , &m_buffer.se_a, &m_buffer.ke_a, &m_buffer.wv_a, &m_buffer.wl_a
     , &m_buffer.kbfs, &m_buffer.ustar2, &m_buffer.wstar
#endif
    };
  for (int i = 0; i < Buffer::num_1d_scalar_ncol; ++i) {
    *_1d_scalar_view_ptrs[i] = scalar_view_t(mem, m_num_cols);
    mem += _1d_scalar_view_ptrs[i]->size();
  }

  Spack* s_mem = reinterpret_cast<Spack*>(mem);

  // 2d packed views
  const int nlev_packs       = ekat::npack<Spack>(m_num_levs);
  const int nlevi_packs      = ekat::npack<Spack>(m_num_levs+1);
  const int num_tracer_packs = ekat::npack<Spack>(m_num_tracers);

  m_buffer.pref_mid = decltype(m_buffer.pref_mid)(s_mem, nlev_packs);
  s_mem += m_buffer.pref_mid.size();

  using spack_2d_view_t = decltype(m_buffer.z_mid);
  spack_2d_view_t* _2d_spack_mid_view_ptrs[Buffer::num_2d_vector_mid] = {
    &m_buffer.z_mid, &m_buffer.rrho, &m_buffer.thv, &m_buffer.dz, &m_buffer.zt_grid, &m_buffer.wm_zt,
    &m_buffer.inv_exner, &m_buffer.thlm, &m_buffer.qw, &m_buffer.dse, &m_buffer.tke_copy, &m_buffer.qc_copy,
    &m_buffer.shoc_ql2, &m_buffer.shoc_mix, &m_buffer.isotropy, &m_buffer.w_sec, &m_buffer.wqls_sec, &m_buffer.brunt
#ifdef SCREAM_SHOC_SMALL_KERNELS
    , &m_buffer.rho_zt, &m_buffer.shoc_qv, &m_buffer.tabs, &m_buffer.dz_zt
#endif
  };

  spack_2d_view_t* _2d_spack_int_view_ptrs[Buffer::num_2d_vector_int] = {
    &m_buffer.z_int, &m_buffer.rrho_i, &m_buffer.zi_grid, &m_buffer.thl_sec, &m_buffer.qw_sec,
    &m_buffer.qwthl_sec, &m_buffer.wthl_sec, &m_buffer.wqw_sec, &m_buffer.wtke_sec, &m_buffer.uw_sec,
    &m_buffer.vw_sec, &m_buffer.w3
#ifdef SCREAM_SHOC_SMALL_KERNELS
    , &m_buffer.dz_zi
#endif
  };

  for (int i = 0; i < Buffer::num_2d_vector_mid; ++i) {
    *_2d_spack_mid_view_ptrs[i] = spack_2d_view_t(s_mem, m_num_cols, nlev_packs);
    s_mem += _2d_spack_mid_view_ptrs[i]->size();
  }

  for (int i = 0; i < Buffer::num_2d_vector_int; ++i) {
    *_2d_spack_int_view_ptrs[i] = spack_2d_view_t(s_mem, m_num_cols, nlevi_packs);
    s_mem += _2d_spack_int_view_ptrs[i]->size();
  }
  m_buffer.wtracer_sfc = decltype(m_buffer.wtracer_sfc)(s_mem, m_num_cols, num_tracer_packs);
  s_mem += m_buffer.wtracer_sfc.size();

  // WSM data
  m_buffer.wsm_data = s_mem;

  // Compute workspace manager size to check used memory
  // vs. requested memory
  const auto policy      = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);
  const int n_wind_slots = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots = ekat::npack<Spack>(m_num_tracers+3)*Spack::n;
  const int wsm_size     = WSM::get_total_bytes_needed(nlevi_packs, 14+(n_wind_slots+n_trac_slots), policy)/sizeof(Spack);
  s_mem += wsm_size;

  size_t used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for SHOCMacrophysics.");
}

// =========================================================================================
void SHOCMacrophysics::initialize_impl (const RunType run_type)
{
  // Gather runtime options
  runtime_options.lambda_low    = m_params.get<double>("lambda_low");
  runtime_options.lambda_high   = m_params.get<double>("lambda_high");
  runtime_options.lambda_slope  = m_params.get<double>("lambda_slope");
  runtime_options.lambda_thresh = m_params.get<double>("lambda_thresh");
  runtime_options.thl2tune      = m_params.get<double>("thl2tune");
  runtime_options.qw2tune       = m_params.get<double>("qw2tune");
  runtime_options.qwthl2tune    = m_params.get<double>("qwthl2tune");
  runtime_options.w2tune        = m_params.get<double>("w2tune");
  runtime_options.length_fac    = m_params.get<double>("length_fac");
  runtime_options.c_diag_3rd_mom = m_params.get<double>("c_diag_3rd_mom");
  runtime_options.Ckh           = m_params.get<double>("Ckh");
  runtime_options.Ckm           = m_params.get<double>("Ckm");
  // Initialize all of the structures that are passed to shoc_main in run_impl.
  // Note: Some variables in the structures are not stored in the field manager.  For these
  //       variables a local view is constructed.
  const auto& T_mid               = get_field_out("T_mid").get_view<Spack**>();
  const auto& p_mid               = get_field_in("p_mid").get_view<const Spack**>();
  const auto& p_int               = get_field_in("p_int").get_view<const Spack**>();
  const auto& pseudo_density      = get_field_in("pseudo_density").get_view<const Spack**>();
  const auto& omega               = get_field_in("omega").get_view<const Spack**>();
  const auto& surf_sens_flux      = get_field_in("surf_sens_flux").get_view<const Real*>();
  const auto& surf_evap           = get_field_in("surf_evap").get_view<const Real*>();
  const auto& surf_mom_flux       = get_field_in("surf_mom_flux").get_view<const Real**>();
  const auto& qtracers            = get_group_out("tracers").m_bundle->get_view<Spack***>();
  const auto& qc                  = get_field_out("qc").get_view<Spack**>();
  const auto& qv                  = get_field_out("qv").get_view<Spack**>();
  const auto& tke                 = get_field_out("tke").get_view<Spack**>();
  const auto& cldfrac_liq         = get_field_out("cldfrac_liq").get_view<Spack**>();
  const auto& cldfrac_liq_prev    = get_field_out("cldfrac_liq_prev").get_view<Spack**>();
  const auto& sgs_buoy_flux       = get_field_out("sgs_buoy_flux").get_view<Spack**>();
  const auto& tk                  = get_field_out("eddy_diff_mom").get_view<Spack**>();
  const auto& inv_qc_relvar       = get_field_out("inv_qc_relvar").get_view<Spack**>();
  const auto& phis                = get_field_in("phis").get_view<const Real*>();

  // Alias local variables from temporary buffer
  auto z_mid       = m_buffer.z_mid;
  auto z_int       = m_buffer.z_int;
  auto wpthlp_sfc  = m_buffer.wpthlp_sfc;
  auto wprtp_sfc   = m_buffer.wprtp_sfc;
  auto upwp_sfc    = m_buffer.upwp_sfc;
  auto vpwp_sfc    = m_buffer.vpwp_sfc;
  auto rrho        = m_buffer.rrho;
  auto rrho_i      = m_buffer.rrho_i;
  auto thv         = m_buffer.thv;
  auto dz          = m_buffer.dz;
  auto zt_grid     = m_buffer.zt_grid;
  auto zi_grid     = m_buffer.zi_grid;
  auto wtracer_sfc = m_buffer.wtracer_sfc;
  auto wm_zt       = m_buffer.wm_zt;
  auto inv_exner   = m_buffer.inv_exner;
  auto thlm        = m_buffer.thlm;
  auto qw          = m_buffer.qw;
  auto dse         = m_buffer.dse;
  auto tke_copy    = m_buffer.tke_copy;
  auto qc_copy     = m_buffer.qc_copy;
  auto shoc_ql2    = m_buffer.shoc_ql2;

  // For now, set z_int(i,nlevs) = z_surf = 0
  const Real z_surf = 0.0;

  // Some SHOC variables should be initialized uniformly if an Initial run
  if (run_type==RunType::Initial){
    Kokkos::deep_copy(sgs_buoy_flux,0.0);
    Kokkos::deep_copy(tk,0.0);
    Kokkos::deep_copy(tke,0.0004);
    Kokkos::deep_copy(tke_copy,0.0004);
    Kokkos::deep_copy(cldfrac_liq,0.0);
  }

  shoc_preprocess.set_variables(m_num_cols,m_num_levs,m_num_tracers,z_surf,
                                T_mid,p_mid,p_int,pseudo_density,omega,phis,surf_sens_flux,surf_evap,
                                surf_mom_flux,qtracers,qv,qc,qc_copy,tke,tke_copy,z_mid,z_int,
                                dse,rrho,rrho_i,thv,dz,zt_grid,zi_grid,wpthlp_sfc,wprtp_sfc,upwp_sfc,vpwp_sfc,
                                wtracer_sfc,wm_zt,inv_exner,thlm,qw, cldfrac_liq, cldfrac_liq_prev);

  // Input Variables:
  input.zt_grid     = shoc_preprocess.zt_grid;
  input.zi_grid     = shoc_preprocess.zi_grid;
  input.pres        = p_mid;
  input.presi       = p_int;
  input.pdel        = pseudo_density;
  input.thv         = shoc_preprocess.thv;
  input.w_field     = shoc_preprocess.wm_zt;
  input.wthl_sfc    = shoc_preprocess.wpthlp_sfc;
  input.wqw_sfc     = shoc_preprocess.wprtp_sfc;
  input.uw_sfc      = shoc_preprocess.upwp_sfc;
  input.vw_sfc      = shoc_preprocess.vpwp_sfc;
  input.wtracer_sfc = shoc_preprocess.wtracer_sfc;
  input.inv_exner   = shoc_preprocess.inv_exner;
  input.phis        = phis;

  // Input/Output Variables
  input_output.host_dse     = shoc_preprocess.shoc_s;
  input_output.tke          = shoc_preprocess.tke_copy;
  input_output.thetal       = shoc_preprocess.thlm;
  input_output.qw           = shoc_preprocess.qw;
  input_output.horiz_wind   = get_field_out("horiz_winds").get_view<Spack***>();
  input_output.wthv_sec     = sgs_buoy_flux;
  input_output.qtracers     = shoc_preprocess.qtracers;
  input_output.tk           = tk;
  input_output.shoc_cldfrac = cldfrac_liq;
  input_output.shoc_ql      = qc_copy;

  // Output Variables
  output.pblh     = get_field_out("pbl_height").get_view<Real*>();
  output.shoc_ql2 = shoc_ql2;
  output.tkh      = get_field_out("eddy_diff_heat").get_view<Spack**>();
  output.ustar    = get_field_out("ustar").get_view<Real*>();
  output.obklen   = get_field_out("obklen").get_view<Real*>();

  // Ouput (diagnostic)
  history_output.shoc_mix  = m_buffer.shoc_mix;
  history_output.isotropy  = m_buffer.isotropy;
  history_output.w_sec     = get_field_out("w_variance").get_view<Spack**>();
  history_output.thl_sec   = m_buffer.thl_sec;
  history_output.qw_sec    = m_buffer.qw_sec;
  history_output.qwthl_sec = m_buffer.qwthl_sec;
  history_output.wthl_sec  = m_buffer.wthl_sec;
  history_output.wqw_sec   = m_buffer.wqw_sec;
  history_output.wtke_sec  = m_buffer.wtke_sec;
  history_output.uw_sec    = m_buffer.uw_sec;
  history_output.vw_sec    = m_buffer.vw_sec;
  history_output.w3        = m_buffer.w3;
  history_output.wqls_sec  = m_buffer.wqls_sec;
  history_output.brunt     = m_buffer.brunt;

#ifdef SCREAM_SHOC_SMALL_KERNELS
  temporaries.se_b = m_buffer.se_b;
  temporaries.ke_b = m_buffer.ke_b;
  temporaries.wv_b = m_buffer.wv_b;
  temporaries.wl_b = m_buffer.wl_b;
  temporaries.se_a = m_buffer.se_a;
  temporaries.ke_a = m_buffer.ke_a;
  temporaries.wv_a = m_buffer.wv_a;
  temporaries.wl_a = m_buffer.wl_a;
  temporaries.kbfs = m_buffer.kbfs;
  temporaries.ustar2 = m_buffer.ustar2;
  temporaries.wstar = m_buffer.wstar;

  temporaries.rho_zt = m_buffer.rho_zt;
  temporaries.shoc_qv = m_buffer.shoc_qv;
  temporaries.tabs = m_buffer.tabs;
  temporaries.dz_zt = m_buffer.dz_zt;
  temporaries.dz_zi = m_buffer.dz_zi;
#endif

  shoc_postprocess.set_variables(m_num_cols,m_num_levs,m_num_tracers,
                                 rrho,qv,qw,qc,qc_copy,tke,tke_copy,qtracers,shoc_ql2,
                                 cldfrac_liq,inv_qc_relvar,
                                 T_mid, dse, z_mid, phis);

  if (has_column_conservation_check()) {
    const auto& vapor_flux = get_field_out("vapor_flux").get_view<Real*>();
    const auto& water_flux = get_field_out("water_flux").get_view<Real*>();
    const auto& ice_flux   = get_field_out("ice_flux").get_view<Real*>();
    const auto& heat_flux  = get_field_out("heat_flux").get_view<Real*>();
    shoc_postprocess.set_mass_and_energy_fluxes (surf_evap, surf_sens_flux,
                                                 vapor_flux, water_flux,
                                                 ice_flux, heat_flux);
  }

  // Set field property checks for the fields in this process
  using Interval = FieldWithinIntervalCheck;
  using LowerBound = FieldLowerBoundCheck;
  add_postcondition_check<Interval>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_postcondition_check<Interval>(get_field_out("qc"),m_grid,0.0,0.1,false);
  add_postcondition_check<Interval>(get_field_out("horiz_winds"),m_grid,-400.0,400.0,false);
  add_postcondition_check<LowerBound>(get_field_out("pbl_height"),m_grid,0);
  add_postcondition_check<Interval>(get_field_out("cldfrac_liq"),m_grid,0.0,1.0,false);
  add_postcondition_check<LowerBound>(get_field_out("tke"),m_grid,0);
  // For qv, ensure it doesn't get negative, by allowing repair of any neg value.
  // TODO: use a repairable lb that clips only "small" negative values
  add_postcondition_check<Interval>(get_field_out("qv"),m_grid,0,0.2,true);

  // Setup WSM for internal local variables
  const auto nlev_packs  = ekat::npack<Spack>(m_num_levs);
  const auto nlevi_packs = ekat::npack<Spack>(m_num_levs+1);
  const int n_wind_slots = ekat::npack<Spack>(2)*Spack::n;
  const int n_trac_slots = ekat::npack<Spack>(m_num_tracers+3)*Spack::n;
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);
  workspace_mgr.setup(m_buffer.wsm_data, nlevi_packs, 14+(n_wind_slots+n_trac_slots), default_policy);

  // Calculate pref_mid, and use that to calculate
  // maximum number of levels in pbl from surface
  const auto pref_mid = m_buffer.pref_mid;
  const auto s_pref_mid = ekat::scalarize(pref_mid);
  const auto hyam = m_grid->get_geometry_data("hyam").get_view<const Real*>();
  const auto hybm = m_grid->get_geometry_data("hybm").get_view<const Real*>();
  const auto ps0 = C::P0;
  const auto psref = ps0;
  Kokkos::parallel_for(Kokkos::RangePolicy<>(0, m_num_levs), KOKKOS_LAMBDA (const int lev) {
    s_pref_mid(lev) = ps0*hyam(lev) + psref*hybm(lev);
  });
  Kokkos::fence();

  const int ntop_shoc = 0;
  const int nbot_shoc = m_num_levs;
  m_npbl = SHF::shoc_init(nbot_shoc,ntop_shoc,pref_mid);

  // Compute cell length for input dx and dy.
  const auto ncols = m_num_cols;
  view_1d cell_length("cell_length", ncols);
  if (m_grid->has_geometry_data("dx_short")) {
    // In this case IOP is running with a planar geometry
    auto dx = m_grid->get_geometry_data("dx_short").get_view<const Real,Host>()();
    Kokkos::deep_copy(cell_length, dx*1000); // convert km -> m
  } else {
    const auto area = m_grid->get_geometry_data("area").get_view<const Real*>();
    const auto lat  = m_grid->get_geometry_data("lat").get_view<const Real*>();
    Kokkos::parallel_for(ncols, KOKKOS_LAMBDA (const int icol) {
      // For now, we are considering dy=dx. Here, we
      // will need to compute dx/dy instead of cell_length
      // if we have dy!=dx.
      cell_length(icol) = PF::calculate_dx_from_area(area(icol),lat(icol));;
    });
  }
  input.dx = cell_length;
  input.dy = cell_length;
}

// =========================================================================================
void SHOCMacrophysics::run_impl (const double dt)
{
  EKAT_REQUIRE_MSG (dt<=300,
      "Error! SHOC is intended to run with a timestep no longer than 5 minutes.\n"
      "       Please, reduce timestep (perhaps increasing subcycling iterations).\n");

  const auto nlev_packs  = ekat::npack<Spack>(m_num_levs);
  const auto scan_policy    = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_num_cols, nlev_packs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);

  // Preprocessing of SHOC inputs. Kernel contains a parallel_scan,
  // so a special TeamPolicy is required.
  Kokkos::parallel_for("shoc_preprocess",
                       scan_policy,
                       shoc_preprocess);
  Kokkos::fence();

  if (m_params.get<bool>("apply_tms", false)) {
    apply_turbulent_mountain_stress();
  }

  if (m_params.get<bool>("check_flux_state_consistency", false)) {
    check_flux_state_consistency(dt);
  }

  // For now set the host timestep to the shoc timestep. This forces
  // number of SHOC timesteps (nadv) to be 1.
  // TODO: input parameter?
  hdtime = dt;
  m_nadv = std::max(static_cast<int>(round(hdtime/dt)),1);

  // Reset internal WSM variables.
  workspace_mgr.reset_internals();

  // Run shoc main
  SHF::shoc_main(m_num_cols, m_num_levs, m_num_levs+1, m_npbl, m_nadv, m_num_tracers, dt,
                 workspace_mgr,runtime_options,input,input_output,output,history_output
#ifdef SCREAM_SHOC_SMALL_KERNELS
                 , temporaries
#endif
                 );

  // Postprocessing of SHOC outputs
  Kokkos::parallel_for("shoc_postprocess",
                       default_policy,
                       shoc_postprocess);
  Kokkos::fence();

  // Extra SHOC output diagnostics
  if (m_params.get<bool>("extra_shoc_diags", false)) {

    const auto& shoc_mix = get_field_out("shoc_mix").get_view<Spack**>();
    Kokkos::deep_copy(shoc_mix,history_output.shoc_mix);

    const auto& brunt = get_field_out("brunt").get_view<Spack**>();
    Kokkos::deep_copy(brunt,history_output.brunt);

    const auto& w3 = get_field_out("w3").get_view<Spack**>();
    Kokkos::deep_copy(w3,history_output.w3);

    const auto& isotropy = get_field_out("isotropy").get_view<Spack**>();
    Kokkos::deep_copy(isotropy,history_output.isotropy);

    const auto& wthl_sec = get_field_out("wthl_sec").get_view<Spack**>();
    Kokkos::deep_copy(wthl_sec,history_output.wthl_sec);

    const auto& wqw_sec = get_field_out("wqw_sec").get_view<Spack**>();
    Kokkos::deep_copy(wqw_sec,history_output.wqw_sec);

    const auto& uw_sec = get_field_out("uw_sec").get_view<Spack**>();
    Kokkos::deep_copy(uw_sec,history_output.uw_sec);

    const auto& vw_sec = get_field_out("vw_sec").get_view<Spack**>();
    Kokkos::deep_copy(vw_sec,history_output.vw_sec);

    const auto& qw_sec = get_field_out("qw_sec").get_view<Spack**>();
    Kokkos::deep_copy(qw_sec,history_output.qw_sec);

    const auto& thl_sec = get_field_out("thl_sec").get_view<Spack**>();
    Kokkos::deep_copy(thl_sec,history_output.thl_sec);

  } // Extra SHOC output diagnostics
}
// =========================================================================================
void SHOCMacrophysics::finalize_impl()
{
  // Do nothing
}
// =========================================================================================
void SHOCMacrophysics::apply_turbulent_mountain_stress()
{
  auto surf_drag_coeff_tms = get_field_in("surf_drag_coeff_tms").get_view<const Real*>();
  auto horiz_winds         = get_field_in("horiz_winds").get_view<const Spack***>();

  auto rrho_i   = m_buffer.rrho_i;
  auto upwp_sfc = m_buffer.upwp_sfc;
  auto vpwp_sfc = m_buffer.vpwp_sfc;

  const int nlev_v  = (m_num_levs-1)/Spack::n;
  const int nlev_p  = (m_num_levs-1)%Spack::n;
  const int nlevi_v = m_num_levs/Spack::n;
  const int nlevi_p = m_num_levs%Spack::n;

  Kokkos::parallel_for("apply_tms", KT::RangePolicy(0, m_num_cols), KOKKOS_LAMBDA (const int i) {
    upwp_sfc(i) -= surf_drag_coeff_tms(i)*horiz_winds(i,0,nlev_v)[nlev_p]/rrho_i(i,nlevi_v)[nlevi_p];
    vpwp_sfc(i) -= surf_drag_coeff_tms(i)*horiz_winds(i,1,nlev_v)[nlev_p]/rrho_i(i,nlevi_v)[nlevi_p];
  });
}
// =========================================================================================
void SHOCMacrophysics::check_flux_state_consistency(const double dt)
{
  using PC = scream::physics::Constants<Real>;
  const Real gravit = PC::gravit;
  const Real qmin   = 1e-12; // minimum permitted constituent concentration (kg/kg)

  const auto& pseudo_density = get_field_in ("pseudo_density").get_view<const Spack**>();
  const auto& surf_evap      = get_field_out("surf_evap").get_view<Real*>();
  const auto& qv             = get_field_out("qv").get_view<Spack**>();

  const auto nlevs           = m_num_levs;
  const auto nlev_packs      = ekat::npack<Spack>(nlevs);
  const auto last_pack_idx   = (nlevs-1)/Spack::n;
  const auto last_pack_entry = (nlevs-1)%Spack::n;
  const auto policy          = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);
  Kokkos::parallel_for("check_flux_state_consistency",
                       policy,
                       KOKKOS_LAMBDA (const KT::MemberType& team) {
    const auto i = team.league_rank();

    const auto& pseudo_density_i = ekat::subview(pseudo_density, i);
    const auto& qv_i             = ekat::subview(qv, i);

    // reciprocal of pseudo_density at the bottom layer
    const auto rpdel = 1.0/pseudo_density_i(last_pack_idx)[last_pack_entry];

    // Check if the negative surface latent heat flux can exhaust
    // the moisture in the lowest model level. If so, apply fixer.
    const auto condition = surf_evap(i) - (qmin - qv_i(last_pack_idx)[last_pack_entry])/(dt*gravit*rpdel);
    if (condition < 0) {
      const auto cc = abs(surf_evap(i)*dt*gravit);

      auto tracer_mass = [&](const int k) {
        return qv_i(k)*pseudo_density_i(k);
      };
      Real mm = ekat::ExeSpaceUtils<KT::ExeSpace>::view_reduction(team, 0, nlevs, tracer_mass);

      EKAT_KERNEL_ASSERT_MSG(mm >= cc, "Error! Total mass of column vapor should be greater than mass of surf_evap.\n");

      Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_packs), [&](const int& k) {
        const auto adjust = cc*qv_i(k)*pseudo_density_i(k)/mm;
        qv_i(k) = (qv_i(k)*pseudo_density_i(k) - adjust)/pseudo_density_i(k);
      });

      surf_evap(i) = 0;
    }
  });
}
// =========================================================================================
} // namespace scream
