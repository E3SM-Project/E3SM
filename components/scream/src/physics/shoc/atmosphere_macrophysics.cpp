#include "ekat/ekat_assert.hpp"
#include "physics/shoc/atmosphere_macrophysics.hpp"

namespace scream
{

// =========================================================================================
SHOCMacrophysics::SHOCMacrophysics (const ekat::Comm& comm,const ekat::ParameterList& params)
  : m_shoc_comm   (comm)
  , m_shoc_params (params)
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
  Units nondim(0,0,0,0,0,0,0);

  const auto& grid_name = m_shoc_params.get<std::string>("Grid");
  auto grid = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for pref_mid_field
  FieldLayout pref_mid_layout{ {LEV}, {m_num_levs} };

  // Layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout_col{ {COL}, {m_num_cols} };

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {m_num_cols,m_num_levs+1} };

  // Layout for horiz_wind field
  FieldLayout horiz_wind_layout { {COL,CMP,LEV}, {m_num_cols,2,m_num_levs} };

  // Define fields needed in SHOC.
  // Note: shoc_main is organized by a set of 5 structures, variables below are organized
  //       using the same approach to make it easier to follow.

  // These variables are needed by the interface, but not actually passed to shoc_main.
  add_required_field("pref_mid",         pref_mid_layout,     Pa,     grid_name);
  add_required_field("T_mid",            scalar3d_layout_mid, nondim, grid_name);
  add_required_field("zi",               scalar3d_layout_int, m,      grid_name);
  add_required_field("zm",               scalar3d_layout_mid, K,      grid_name);
  add_required_field("omega",            scalar3d_layout_mid, K,      grid_name);
  add_required_field("surf_sens_flux",   scalar2d_layout_col, K,      grid_name);
  add_required_field("surf_latent_flux", scalar2d_layout_col, K,      grid_name);
  add_required_field("surf_u_mom_flux",  scalar2d_layout_col, K,      grid_name);
  add_required_field("surf_v_mom_flux",  scalar2d_layout_col, K,      grid_name);
  add_required_field("qv",               scalar3d_layout_mid, Qunit,  grid_name);

  add_computed_field("T_mid",            scalar3d_layout_mid, nondim, grid_name);
  add_computed_field("qv",               scalar3d_layout_mid, Qunit,  grid_name);

  // Input variables
  add_required_field("host_dx",        scalar2d_layout_col, m,  grid_name);
  add_required_field("host_dy",        scalar2d_layout_col, m,  grid_name);
  add_required_field("p_mid",          scalar3d_layout_mid, Pa, grid_name);
  add_required_field("pint",           scalar3d_layout_int, Pa, grid_name);
  add_required_field("pseudo_density", scalar3d_layout_mid, Pa, grid_name);
  add_required_field("phis",           scalar2d_layout_col, m,  grid_name);

  // Input/Output variables
  add_required_field("tke",           scalar3d_layout_mid, (m*m)/(s*s), grid_name);
  add_required_field("horiz_winds",   horiz_wind_layout,   m/s,         grid_name);
  add_required_field("sgs_buoy_flux", scalar3d_layout_mid, K*(m/s),     grid_name);
  add_required_field("eddy_diff_mom", scalar3d_layout_mid, (m*m)/s,     grid_name);
  add_required_field("qc",            scalar3d_layout_mid, Qunit,           grid_name);
  add_required_field("cldfrac_liq",   scalar3d_layout_mid, Pa,          grid_name);

  add_computed_field("tke",           scalar3d_layout_mid, (m*m)/(s*s), grid_name);
  add_computed_field("horiz_winds",   horiz_wind_layout,   m/s,         grid_name);
  add_computed_field("sgs_buoy_flux", scalar3d_layout_mid, K*(m/s),     grid_name);
  add_computed_field("eddy_diff_mom", scalar3d_layout_mid, (m*m)/s,     grid_name);
  add_computed_field("qc",            scalar3d_layout_mid, Qunit,       grid_name);
  add_computed_field("cldfrac_liq",   scalar3d_layout_mid, Pa,          grid_name);

  // Output variables
  add_computed_field("pbl_height",    scalar2d_layout_col, m, grid_name);
  add_computed_field("inv_qc_relvar", scalar3d_layout_mid, Qunit*Qunit, grid_name);

  // Tracer group
  m_inout_groups_req.emplace("tracers",grid->name());
}
// =========================================================================================
void SHOCMacrophysics::
set_updated_group (const FieldGroup<Real>& group)
{
  EKAT_REQUIRE_MSG(group.m_info->size() >= 3,
                   "Error! Shoc requires at least 3 tracers (tke, qv, qc) as inputs.");

  const auto& name = group.m_info->m_group_name;

  EKAT_REQUIRE_MSG(name=="tracers",
    "Error! We were not expecting a field group called '" << name << "\n");

  EKAT_REQUIRE_MSG(group.m_info->m_bundled,
      "Error! Shoc expects bundled fields for tracers.\n");

  // Add Q bundle as in/out field
  m_shoc_fields_in["Q"]  = *group.m_bundle;
  m_shoc_fields_out["Q"] = *group.m_bundle;

  // Calculate number of advected tracers
  m_num_tracers = group.m_info->size();
}

// =========================================================================================
void SHOCMacrophysics::initialize_impl (const util::TimeStamp& t0)
{
  m_current_ts = t0;

  // Initialize all of the structures that are passed to shoc_main in run_impl.
  // Note: Some variables in the structures are not stored in the field manager.  For these
  //       variables a local view is constructed.
  auto T_mid            = m_shoc_fields_out["T_mid"].get_reshaped_view<Spack**>();
  auto z_int            = m_shoc_fields_in["zi"].get_reshaped_view<const Spack**>();
  auto z_mid            = m_shoc_fields_in["zm"].get_reshaped_view<const Spack**>();
  auto p_mid            = m_shoc_fields_in["p_mid"].get_reshaped_view<const Spack**>();
  auto pseudo_density   = m_shoc_fields_in["pseudo_density"].get_reshaped_view<const Spack**>();
  auto omega            = m_shoc_fields_in["omega"].get_reshaped_view<const Spack**>();
  auto surf_sens_flux   = m_shoc_fields_in["surf_sens_flux"].get_reshaped_view<const Pack1d*>();
  auto surf_latent_flux = m_shoc_fields_in["surf_latent_flux"].get_reshaped_view<const Pack1d*>();
  auto surf_u_mom_flux  = m_shoc_fields_in["surf_u_mom_flux"].get_reshaped_view<const Pack1d*>();
  auto surf_v_mom_flux  = m_shoc_fields_in["surf_v_mom_flux"].get_reshaped_view<const Pack1d*>();
  auto qc               = m_shoc_fields_out["qc"].get_reshaped_view<Spack**>();
  auto qv               = m_shoc_fields_out["qv"].get_reshaped_view<Spack**>();
  auto tke              = m_shoc_fields_out["tke"].get_reshaped_view<Spack**>();
  auto Q                = m_shoc_fields_out["Q"].get_reshaped_view<Spack***>();
  auto cldfrac_liq      = m_shoc_fields_out["cldfrac_liq"].get_reshaped_view<Spack**>();
  auto sgs_buoy_flux    = m_shoc_fields_out["sgs_buoy_flux"].get_reshaped_view<Spack**>();
  auto inv_qc_relvar    = m_shoc_fields_out["inv_qc_relvar"].get_reshaped_view<Spack**>();
  auto phis             = m_shoc_fields_in["phis"].get_reshaped_view<const Pack1d*>();

  const int nlev_packs = ekat::npack<Spack>(m_num_levs);
  const int nlevi_packs = ekat::npack<Spack>(m_num_levs+1);
  const int num_tracer_packs = ekat::npack<Spack>(m_num_tracers);

  view_1d wpthlp_sfc("wpthlp_sfc",m_num_cols),
          wprtp_sfc("wprtp_sfc",m_num_cols),
          upwp_sfc("upwp_sfc",m_num_cols),
          vpwp_sfc("vpwp_sfc",m_num_cols);

  view_2d rrho("rrho",m_num_cols,nlev_packs),
          rrho_i("rrhoi",m_num_cols,nlevi_packs),
          thv("thv",m_num_cols,nlev_packs),
          dz("dz",m_num_cols,nlev_packs),
          zt_grid("zt_grid",m_num_cols,nlev_packs),
          zi_grid("zi_grid",m_num_cols,nlevi_packs),
          wtracer_sfc("wtracer_sfc",m_num_cols,num_tracer_packs),
          wm_zt("wm_zt",m_num_cols,nlev_packs),
          exner("exner",m_num_cols,nlev_packs),
          thlm("thlm",m_num_cols,nlev_packs),
          qw("qw",m_num_cols,nlev_packs),
          s("s",m_num_cols,nlev_packs),
          qv_copy("qv_copy",m_num_cols,nlev_packs),
          qc_copy("qc_copy",m_num_cols,nlev_packs),
          tke_copy("tke_copy",m_num_cols,nlev_packs);

  // TODO: Transpose of the tracers should be handled internally in shoc,
  //       removing this allocation.
  view_3d tracers("tracers",m_num_cols,m_num_levs,num_tracer_packs);

  shoc_preprocess.set_variables(m_num_cols,m_num_levs,m_num_tracers,nlev_packs,num_tracer_packs,T_mid,
                                z_int,z_mid,p_mid,pseudo_density,omega,phis,surf_sens_flux,surf_latent_flux,
                                surf_u_mom_flux,surf_v_mom_flux,qv,qv_copy,Q,qc,qc_copy,tke,tke_copy,
                                s,rrho,rrho_i,thv,dz,zt_grid,zi_grid,wpthlp_sfc,wprtp_sfc,upwp_sfc,vpwp_sfc,
                                wtracer_sfc,wm_zt,exner,thlm,qw,tracers);

  // Input Variables:
  input.host_dx     = m_shoc_fields_in["host_dx"].get_reshaped_view<const Pack1d*>();
  input.host_dy     = m_shoc_fields_in["host_dy"].get_reshaped_view<const Pack1d*>();
  input.zt_grid     = shoc_preprocess.zt_grid;
  input.zi_grid     = shoc_preprocess.zi_grid;
  input.pres        = p_mid;
  input.presi       = m_shoc_fields_in["pint"].get_reshaped_view<const Spack**>();
  input.pdel        = pseudo_density;
  input.thv         = shoc_preprocess.thv;
  input.w_field     = shoc_preprocess.wm_zt;
  input.wthl_sfc    = shoc_preprocess.wpthlp_sfc;
  input.wqw_sfc     = shoc_preprocess.wprtp_sfc;
  input.uw_sfc      = shoc_preprocess.upwp_sfc;
  input.vw_sfc      = shoc_preprocess.vpwp_sfc;
  input.wtracer_sfc = shoc_preprocess.wtracer_sfc;
  input.exner       = shoc_preprocess.exner;
  input.phis        = phis;

  // Input/Output Variables
  input_output.host_dse     = shoc_preprocess.shoc_s;
  input_output.tke          = shoc_preprocess.tke_copy;
  input_output.thetal       = shoc_preprocess.thlm;
  input_output.qw           = shoc_preprocess.qw;
  input_output.horiz_wind   = m_shoc_fields_out["horiz_winds"].get_reshaped_view<Spack***>();
  input_output.wthv_sec     = sgs_buoy_flux;
  input_output.qtracers     = shoc_preprocess.tracers;
  input_output.tk           = m_shoc_fields_out["eddy_diff_mom"].get_reshaped_view<Spack**>();
  input_output.shoc_cldfrac = cldfrac_liq;
  input_output.shoc_ql      = shoc_preprocess.qc_copy;

  // Output Variables
  output.pblh = m_shoc_fields_out["pbl_height"].get_reshaped_view<Pack1d*>();

  view_2d shoc_ql2("shoc_ql2",m_num_cols,nlev_packs);
  output.shoc_ql2 = shoc_ql2;

  // Ouput (diagnostic)
  // TODO: A temporary buffer should be added to the AD for these local
  //       views.
  view_2d shoc_mix("shoc_mix",m_num_cols,nlev_packs);
  view_2d isotropy("isotropy",m_num_cols,nlev_packs);
  view_2d w_sec("w_sec",m_num_cols,nlev_packs);
  view_2d thl_sec("thl_sec",m_num_cols,nlevi_packs);
  view_2d qw_sec("qw_sec",m_num_cols,nlevi_packs);
  view_2d qwthl_sec("qwthl_sec",m_num_cols,nlevi_packs);
  view_2d wthl_sec("wthl_sec",m_num_cols,nlevi_packs);
  view_2d wqw_sec("wqw_sec",m_num_cols,nlevi_packs);
  view_2d wtke_sec("wtke_sec",m_num_cols,nlevi_packs);
  view_2d uw_sec("uw_sec",m_num_cols,nlevi_packs);
  view_2d vw_sec("vw_sec",m_num_cols,nlevi_packs);
  view_2d w3("w3",m_num_cols,nlevi_packs);
  view_2d wqls_sec("wqls_sec",m_num_cols,nlev_packs);
  view_2d brunt("brunt",m_num_cols,nlev_packs);

  history_output.shoc_mix  = shoc_mix;
  history_output.isotropy  = isotropy;
  history_output.w_sec     = w_sec;
  history_output.thl_sec   = thl_sec;
  history_output.qw_sec    = qw_sec;
  history_output.qwthl_sec = qwthl_sec;
  history_output.wthl_sec  = wthl_sec;
  history_output.wqw_sec   = wqw_sec;
  history_output.wtke_sec  = wtke_sec;
  history_output.uw_sec    = uw_sec;
  history_output.vw_sec    = vw_sec;
  history_output.w3        = w3;
  history_output.wqls_sec  = wqls_sec;
  history_output.brunt     = brunt;

  shoc_postprocess.set_variables(m_num_cols,m_num_levs,m_num_tracers,nlev_packs,num_tracer_packs,
                                 rrho,qv,qv_copy,qc,qc_copy,tke,tke_copy,shoc_ql2,Q,tracers,
                                 cldfrac_liq,sgs_buoy_flux,inv_qc_relvar);
}

// =========================================================================================
void SHOCMacrophysics::run_impl (const Real dt)
{
  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_shoc_fields_in) {
    it.second.sync_to_host();
  }
  for (auto& it : m_shoc_fields_out) {
    it.second.sync_to_host();
  }

  {
    const auto nlev_packs = ekat::npack<Spack>(m_num_cols);
    const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);
    Kokkos::parallel_for("shoc_preprocess",
                         policy,
                         shoc_preprocess);
  }
  Kokkos::fence();


  // Calculate maximum number of levels in pbl from surface
  const auto pref_mid = m_shoc_fields_in["pref_mid"].get_reshaped_view<const Spack*>();
  const int ntop_shoc = 0;
  const int nbot_shoc = m_num_levs;
  m_npbl = SHF::shoc_init(nbot_shoc,ntop_shoc,pref_mid);

  // For now set the host timestep to the shoc timestep. This forces
  // number of SHOC timesteps (nadv) to be 1.
  // TODO: input parameter?
  hdtime = dt;
  m_nadv = ekat::impl::max(hdtime/dt,sp(1));

  // shoc_main() expects 3 extra slots in qtracer array
  // used for solving
  // TODO: This should be handled internally (with tracer transpose)
  Kokkos::resize(input_output.qtracers,
                 m_num_cols,m_num_levs,ekat::npack<Spack>(m_num_tracers+3));

  // Run shoc main
  SHF::shoc_main(m_num_cols, m_num_levs, m_num_levs+1, m_npbl, m_nadv, m_num_tracers, dt,
                 input,input_output,output,history_output);

  // Remove extra slots
  // TODO: This should be handled internally (with tracer transpose)
  Kokkos::resize(input_output.qtracers,
                 m_num_cols,m_num_levs,ekat::npack<Spack>(m_num_tracers));

  {
    const auto nlev_packs = ekat::npack<Spack>(m_num_cols);
    const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nlev_packs);
    Kokkos::parallel_for("shoc_postprocess",
                         policy,
                         shoc_postprocess);
  }
  Kokkos::fence();

  // Get a copy of the current timestamp (at the beginning of the step) and
  // advance it, updating the shoc fields.
  auto ts = timestamp();
  ts += dt;
  //Q->get_header().get_tracking().update_time_stamp(ts);
  for (auto& f : m_shoc_fields_out) {
    f.second.get_header().get_tracking().update_time_stamp(ts);
  }
}
// =========================================================================================
void SHOCMacrophysics::finalize_impl()
{
  // Do nothing
}
// =========================================================================================

void SHOCMacrophysics::register_fields (FieldRepository<Real>& field_repo) const {
  std::set<ci_string> q_names =
    { "qc", "qv", "tke"};

  for (const auto& fid : get_required_fields()) {
    const auto& name = fid.name();
    if (q_names.count(name)>0) {
      field_repo.register_field<Spack>(fid,"tracers");
    } else {
      field_repo.register_field<Spack>(fid);
    }
  }
  for (const auto& fid : get_computed_fields()) {
    const auto& name = fid.name();
    if (q_names.count(name)>0) {
      field_repo.register_field<Spack>(fid,"tracers");
    } else {
      field_repo.register_field<Spack>(fid);
    }
  }
}

void SHOCMacrophysics::set_required_field_impl (const Field<const Real>& f) {

  const auto& name = f.get_header().get_identifier().name();
  m_shoc_fields_in.emplace(name,f);

  // Add myself as customer to the field
  add_me_as_customer(f);
}

void SHOCMacrophysics::set_computed_field_impl (const Field<Real>& f) {

  const auto& name = f.get_header().get_identifier().name();
  m_shoc_fields_out.emplace(name,f);

  // Add myself as provider for the field
  add_me_as_provider(f);
}

} // namespace scream
