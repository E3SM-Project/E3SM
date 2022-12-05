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

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Qunit = kg/kg;
  Qunit.set_string("kg/kg");
  Units nondim(0,0,0,0,0,0,0);

  grid_ = grids_manager->get_grid("Physics");
  const auto& grid_name = grid_->name();

  num_cols_ = grid_->get_num_local_dofs(); // Number of columns on this rank
  num_levs_ = grid_->get_num_vertical_levels();  // Number of levels per column

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
  auto z_mid       = m_buffer.z_mid;
  auto z_int       = m_buffer.z_int;
  auto cell_length = m_buffer.cell_length;
  auto dz          = m_buffer.dz;
  auto zt_grid     = m_buffer.zt_grid;
  auto zi_grid     = m_buffer.zi_grid;
  auto wtracer_sfc = m_buffer.wtracer_sfc;
  auto wm_zt       = m_buffer.wm_zt;

  // Perform any initialization work.
  if (run_type==RunType::Initial){
    /*
    Kokkos::deep_copy(sgs_buoy_flux,0.0);
    Kokkos::deep_copy(tk,0.0);
    Kokkos::deep_copy(tke,0.0004);
    Kokkos::deep_copy(tke_copy,0.0004);
    Kokkos::deep_copy(cldfrac_liq,0.0);
    */
  }

  // Find index of qv (water vapor, kg/kg(wet-air) and tke (J/kg(wet-air)) in the qtracer 3d view
  // These indices are later used for converting tracers from wet mmr to dry mmr and vice-versa
  auto qv_index  = tracer_info->m_subview_idx.at("qv");

  //Device view to store indices of tracers which will participate in wet<->dry conversion; we are excluding
  //"tke" [as it is not "water based" tracer] and "qv"[as "qv" (before conversion) is needed for
  //computing conversion for all other tracers] from qtracers view
  view_1d_int convert_wet_dry_idx_d("convert_wet_dry_idx_d",m_num_tracers-1);  // qv is excluded

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

  shoc_postprocess.set_variables(m_num_cols,m_num_levs,m_num_tracers,convert_wet_dry_idx_d,
                                 rrho,qv,qw,qc,qc_copy,tke,tke_copy,qtracers,shoc_ql2,
                                 cldfrac_liq,inv_qc_relvar,
                                 T_mid, dse, z_mid, phis);

  // Set field property checks for the fields in this process
  using Interval = FieldWithinIntervalCheck;
  using LowerBound = FieldLowerBoundCheck;
  add_postcondition_check<Interval>(get_field_out("T_mid"),m_grid,130.0,500.0,false);
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
  const auto hyam = m_grid->get_geometry_data("hyam");
  const auto hybm = m_grid->get_geometry_data("hybm");
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

void MAM4Aerosols::run_impl(const int dt) {

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
  m_nadv = std::max(hdtime/dt,1);

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

void MAM4Aerosols::finalize_impl()
{
}

} // namespace scream
