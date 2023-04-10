#include "ekat/ekat_assert.hpp"
#include "physics/shoc/atmosphere_macrophysics.hpp"

#include "share/property_checks/field_lower_bound_check.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "scream_config.h" // for SCREAM_CIME_BUILD

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

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Qunit = kg/kg;
  Qunit.set_string("kg/kg");
  auto nondim = Units::nondimensional();

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();

  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  m_cell_area = m_grid->get_geometry_data("area").get_view<const Real*>(); // area of each cell
  m_cell_lat  = m_grid->get_geometry_data("lat").get_view<const Real*>(); // area of each cell

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout_col{ {COL}, {m_num_cols} };

  // Layout for surf_mom_flux
  FieldLayout  surf_mom_flux_layout { {COL, CMP}, {m_num_cols, 2} };

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {m_num_cols,m_num_levs+1} };

  // Layout for horiz_wind field
  FieldLayout horiz_wind_layout { {COL,CMP,LEV}, {m_num_cols,2,m_num_levs} };

  // Define fields needed in SHOC.
  // Note: shoc_main is organized by a set of 5 structures, variables below are organized
  //       using the same approach to make it easier to follow.

  constexpr int ps = Spack::n;

  const auto m2 = m*m;
  const auto s2 = s*s;

  // These variables are needed by the interface, but not actually passed to shoc_main.
  add_field<Required>("omega",            scalar3d_layout_mid,  Pa/s,    grid_name, ps);
  add_field<Required>("surf_sens_flux",   scalar2d_layout_col,  W/m2,    grid_name);
  add_field<Required>("surf_evap",        scalar2d_layout_col,  kg/m2/s, grid_name);
  add_field<Required>("surf_mom_flux",    surf_mom_flux_layout, N/m2,    grid_name);

  add_field<Updated> ("T_mid",            scalar3d_layout_mid, K,       grid_name, ps);
  add_field<Updated> ("qv",               scalar3d_layout_mid, Qunit,   grid_name, "tracers", ps);

  // Input variables
  add_field<Required>("p_mid",          scalar3d_layout_mid, Pa,    grid_name, ps);
  add_field<Required>("p_int",          scalar3d_layout_int, Pa,    grid_name, ps);
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa,    grid_name, ps);
  add_field<Required>("phis",           scalar2d_layout_col, m2/s2, grid_name, ps);

  // Input/Output variables
  add_field<Updated>("tke",           scalar3d_layout_mid, m2/s2,   grid_name, "tracers", ps);
  add_field<Updated>("horiz_winds",   horiz_wind_layout,   m/s,     grid_name, ps);
  add_field<Updated>("sgs_buoy_flux", scalar3d_layout_mid, K*(m/s), grid_name, ps);
  add_field<Updated>("eddy_diff_mom", scalar3d_layout_mid, m2/s,    grid_name, ps);
  add_field<Updated>("qc",            scalar3d_layout_mid, Qunit,   grid_name, "tracers", ps);
  add_field<Updated>("cldfrac_liq",   scalar3d_layout_mid, nondim,  grid_name, ps);

  // Output variables
  add_field<Computed>("pbl_height",    scalar2d_layout_col, m,           grid_name);
  add_field<Computed>("inv_qc_relvar", scalar3d_layout_mid, Qunit*Qunit, grid_name, ps);

  // Tracer group
  add_group<Updated>("tracers",grid_name,ps,Bundling::Required);

  // Boundary flux fields for energy and mass conservation checks
  if (has_column_conservation_check()) {
    add_field<Computed>("vapor_flux", scalar2d_layout_col, kg/m2/s, grid_name);
    add_field<Computed>("water_flux", scalar2d_layout_col, m/s,     grid_name);
    add_field<Computed>("ice_flux",   scalar2d_layout_col, m/s,     grid_name);
    add_field<Computed>("heat_flux",  scalar2d_layout_col, W/m2,    grid_name);
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
  const size_t wsm_request= WSM::get_total_bytes_needed(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy);

  return interface_request + wsm_request;
}

// =========================================================================================
void SHOCMacrophysics::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d scalar views
  using scalar_view_t = decltype(m_buffer.cell_length);
  scalar_view_t* _1d_scalar_view_ptrs[Buffer::num_1d_scalar_ncol] =
    {&m_buffer.cell_length, &m_buffer.wpthlp_sfc, &m_buffer.wprtp_sfc, &m_buffer.upwp_sfc, &m_buffer.vpwp_sfc
#ifdef SCREAM_SMALL_KERNELS
     , &m_buffer.se_b, &m_buffer.ke_b, &m_buffer.wv_b, &m_buffer.wl_b
     , &m_buffer.se_a, &m_buffer.ke_a, &m_buffer.wv_a, &m_buffer.wl_a
     , &m_buffer.ustar, &m_buffer.kbfs, &m_buffer.obklen, &m_buffer.ustar2, &m_buffer.wstar
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
#ifdef SCREAM_SMALL_KERNELS
    , &m_buffer.rho_zt, &m_buffer.shoc_qv, &m_buffer.dz_zt, &m_buffer.tkh
#endif
  };

  spack_2d_view_t* _2d_spack_int_view_ptrs[Buffer::num_2d_vector_int] = {
    &m_buffer.z_int, &m_buffer.rrho_i, &m_buffer.zi_grid, &m_buffer.thl_sec, &m_buffer.qw_sec,
    &m_buffer.qwthl_sec, &m_buffer.wthl_sec, &m_buffer.wqw_sec, &m_buffer.wtke_sec, &m_buffer.uw_sec,
    &m_buffer.vw_sec, &m_buffer.w3
#ifdef SCREAM_SMALL_KERNELS
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
  const int wsm_size     = WSM::get_total_bytes_needed(nlevi_packs, 13+(n_wind_slots+n_trac_slots), policy)/sizeof(Spack);
  s_mem += wsm_size;

  size_t used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for SHOCMacrophysics.");
}

// =========================================================================================
void SHOCMacrophysics::initialize_impl (const RunType run_type)
{
  // Initialize all of the structures that are passed to shoc_main in run_impl.
  // Note: Some variables in the structures are not stored in the field manager.  For these
  //       variables a local view is constructed.
  const auto& T_mid          = get_field_out("T_mid").get_view<Spack**>();
  const auto& p_mid          = get_field_in("p_mid").get_view<const Spack**>();
  const auto& p_int          = get_field_in("p_int").get_view<const Spack**>();
  const auto& pseudo_density = get_field_in("pseudo_density").get_view<const Spack**>();
  const auto& omega          = get_field_in("omega").get_view<const Spack**>();
  const auto& surf_sens_flux = get_field_in("surf_sens_flux").get_view<const Real*>();
  const auto& surf_evap      = get_field_in("surf_evap").get_view<const Real*>();
  const auto& surf_mom_flux  = get_field_in("surf_mom_flux").get_view<const Real**>();
  const auto& qtracers       = get_group_out("tracers").m_bundle->get_view<Spack***>();
  const auto& qc             = get_field_out("qc").get_view<Spack**>();
  const auto& qv             = get_field_out("qv").get_view<Spack**>();
  const auto& tke            = get_field_out("tke").get_view<Spack**>();
  const auto& cldfrac_liq    = get_field_out("cldfrac_liq").get_view<Spack**>();
  const auto& sgs_buoy_flux  = get_field_out("sgs_buoy_flux").get_view<Spack**>();
  const auto& tk             = get_field_out("eddy_diff_mom").get_view<Spack**>();
  const auto& inv_qc_relvar  = get_field_out("inv_qc_relvar").get_view<Spack**>();
  const auto& phis           = get_field_in("phis").get_view<const Real*>();
  const auto& tracer_info    = get_group_out("tracers").m_info; // obtain tracer info structure

  // Alias local variables from temporary buffer
  auto z_mid       = m_buffer.z_mid;
  auto z_int       = m_buffer.z_int;
  auto cell_length = m_buffer.cell_length;
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
  // Find index of qv (water vapor, kg/kg(wet-air) and tke (J/kg(wet-air)) in the qtracer 3d view
  // These indices are later used for converting tracers from wet mmr to dry mmr and vice-versa
  auto qv_index  = tracer_info->m_subview_idx.at("qv");
  auto tke_index = tracer_info->m_subview_idx.at("tke");

  //Device view to store indices of tracers which will participate in wet<->dry conversion; we are excluding
  //"tke" [as it is not "water based" tracer] and "qv"[as "qv" (before conversion) is needed for
  //computing conversion for all other tracers] from qtracers view

  view_1d_int convert_wet_dry_idx_d("convert_wet_dry_idx_d",m_num_tracers-2);  //2 tracers, qv and tke, are excluded

  //mirror view on host
  auto convert_wet_dry_idx_h = Kokkos::create_mirror_view(convert_wet_dry_idx_d);
  //loop over all tracers to store of all tracer indices except for tke and qv
  for (int it=0,iq=0; it<m_num_tracers; ++it) {
    if (it!=qv_index && it!= tke_index) { //skip if "it" is a tke or qv index
      convert_wet_dry_idx_h(iq) = it;
      ++iq;
    }
  }

  // copy to device
  Kokkos::deep_copy(convert_wet_dry_idx_d,convert_wet_dry_idx_h);


  shoc_preprocess.set_variables(m_num_cols,m_num_levs,m_num_tracers,convert_wet_dry_idx_d,z_surf,m_cell_area,m_cell_lat,
                                T_mid,p_mid,p_int,pseudo_density,omega,phis,surf_sens_flux,surf_evap,
                                surf_mom_flux,qtracers,qv,qc,qc_copy,tke,tke_copy,z_mid,z_int,cell_length,
                                dse,rrho,rrho_i,thv,dz,zt_grid,zi_grid,wpthlp_sfc,wprtp_sfc,upwp_sfc,vpwp_sfc,
                                wtracer_sfc,wm_zt,inv_exner,thlm,qw);

  // Input Variables:
  input.dx          = shoc_preprocess.cell_length;
  input.dy          = shoc_preprocess.cell_length;
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

  // Ouput (diagnostic)
  history_output.shoc_mix  = m_buffer.shoc_mix;
  history_output.isotropy  = m_buffer.isotropy;
  history_output.w_sec     = m_buffer.w_sec;
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

#ifdef SCREAM_SMALL_KERNELS
  temporaries.se_b = m_buffer.se_b;
  temporaries.ke_b = m_buffer.ke_b;
  temporaries.wv_b = m_buffer.wv_b;
  temporaries.wl_b = m_buffer.wl_b;
  temporaries.se_a = m_buffer.se_a;
  temporaries.ke_a = m_buffer.ke_a;
  temporaries.wv_a = m_buffer.wv_a;
  temporaries.wl_a = m_buffer.wl_a;
  temporaries.ustar = m_buffer.ustar;
  temporaries.kbfs = m_buffer.kbfs;
  temporaries.obklen = m_buffer.obklen;
  temporaries.ustar2 = m_buffer.ustar2;
  temporaries.wstar = m_buffer.wstar;

  temporaries.rho_zt = m_buffer.rho_zt;
  temporaries.shoc_qv = m_buffer.shoc_qv;
  temporaries.dz_zt = m_buffer.dz_zt;
  temporaries.dz_zi = m_buffer.dz_zi;
  temporaries.tkh = m_buffer.tkh;
#endif

  shoc_postprocess.set_variables(m_num_cols,m_num_levs,m_num_tracers,convert_wet_dry_idx_d,
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
  workspace_mgr.setup(m_buffer.wsm_data, nlevi_packs, 13+(n_wind_slots+n_trac_slots), default_policy);

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
}

// =========================================================================================
void SHOCMacrophysics::run_impl (const double dt)
{
  EKAT_REQUIRE_MSG (dt<=300,
      "Error! SHOC is intended to run with a timestep no longer than 5 minutes.\n"
      "       Please, reduce timestep (perhaps increasing subcycling iteratinos).\n");

  const auto nlev_packs  = ekat::npack<Spack>(m_num_levs);
  const auto scan_policy    = ekat::ExeSpaceUtils<KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(m_num_cols, nlev_packs);
  const auto default_policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);

  // Preprocessing of SHOC inputs. Kernel contains a parallel_scan,
  // so a special TeamPolicy is required.
  Kokkos::parallel_for("shoc_preprocess",
                       scan_policy,
                       shoc_preprocess);
  Kokkos::fence();

  // For now set the host timestep to the shoc timestep. This forces
  // number of SHOC timesteps (nadv) to be 1.
  // TODO: input parameter?
  hdtime = dt;
  m_nadv = std::max(static_cast<int>(round(hdtime/dt)),1);

  // Reset internal WSM variables.
  workspace_mgr.reset_internals();

  // Run shoc main
  SHF::shoc_main(m_num_cols, m_num_levs, m_num_levs+1, m_npbl, m_nadv, m_num_tracers, dt,
                 workspace_mgr,input,input_output,output,history_output
#ifdef SCREAM_SMALL_KERNELS
                 , temporaries
#endif
                 );

  // Postprocessing of SHOC outputs
  Kokkos::parallel_for("shoc_postprocess",
                       default_policy,
                       shoc_postprocess);
  Kokkos::fence();
}
// =========================================================================================
void SHOCMacrophysics::finalize_impl()
{
  // Do nothing
}
// =========================================================================================

} // namespace scream
