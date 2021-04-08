#include "physics/p3/atmosphere_microphysics.hpp"
// Needed for p3_init, the only F90 code still used.
#include "physics/p3/p3_f90.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
/*
 * P3 Microphysics routines
*/

  using namespace p3;

  using view_1d  = typename P3F::view_1d<Real>;
  using view_2d  = typename P3F::view_2d<Spack>;
  using sview_2d = typename KokkosTypes<DefaultDevice>::template view_2d<Real>;

// =========================================================================================
P3Microphysics::P3Microphysics (const ekat::Comm& comm, const ekat::ParameterList& params)
 : m_p3_comm (comm)
 , m_p3_params (params)
{
  // Nothing to do here
}

// =========================================================================================
void P3Microphysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  Units nondim(0,0,0,0,0,0,0);
  auto mm = m/1000;

  const auto& grid_name = m_p3_params.get<std::string>("Grid");
  auto grid = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  FieldLayout scalar3d_layout_int { {COL,ILEV}, {m_num_cols,m_num_levs+1} };

  // Define fields needed in P3.
  // Note: p3_main is organized by a set of 5 structures, variables below are organized
  //       using the same approach to make it easier to follow.

  // These variables are needed by the interface, but not actually passed to p3_main. 
  add_required_field("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name);
  add_required_field("p_mid",       scalar3d_layout_mid, Pa,     grid_name);
  add_required_field("zi",          scalar3d_layout_int, m,      grid_name);
  add_required_field("T_mid",       scalar3d_layout_mid, K,      grid_name);
  add_computed_field("T_mid",       scalar3d_layout_mid, K,      grid_name);  // T_mid is the only one of these variables that is also updated.

  // Prognostic State:  (all fields are both input and output)
  add_required_field("qv",     scalar3d_layout_mid, Q,    grid_name);
  add_required_field("qc",     scalar3d_layout_mid, Q,    grid_name);
  add_required_field("qr",     scalar3d_layout_mid, Q,    grid_name);
  add_required_field("qi",     scalar3d_layout_mid, Q,    grid_name);
  add_required_field("qm",     scalar3d_layout_mid, Q,    grid_name);
  add_required_field("nc",     scalar3d_layout_mid, 1/kg, grid_name);
  add_required_field("nr",     scalar3d_layout_mid, 1/kg, grid_name);
  add_required_field("ni",     scalar3d_layout_mid, 1/kg, grid_name);
  add_required_field("bm",     scalar3d_layout_mid, 1/kg, grid_name);
  //
  add_computed_field("qv",     scalar3d_layout_mid, Q,    grid_name);
  add_computed_field("qc",     scalar3d_layout_mid, Q,    grid_name);
  add_computed_field("qr",     scalar3d_layout_mid, Q,    grid_name);
  add_computed_field("qi",     scalar3d_layout_mid, Q,    grid_name);
  add_computed_field("qm",     scalar3d_layout_mid, Q,    grid_name);
  add_computed_field("nc",     scalar3d_layout_mid, 1/kg, grid_name);
  add_computed_field("nr",     scalar3d_layout_mid, 1/kg, grid_name);
  add_computed_field("ni",     scalar3d_layout_mid, 1/kg, grid_name);
  add_computed_field("bm",     scalar3d_layout_mid, 1/kg, grid_name);
  // Diagnostic Inputs: (only the X_prev fields are both input and output, all others are just inputs)
  add_required_field("nc_nuceat_tend",     scalar3d_layout_mid, 1/(kg*s), grid_name);
  add_required_field("nc_activated",       scalar3d_layout_mid, nondim,   grid_name);
  add_required_field("ni_activated",       scalar3d_layout_mid, 1/kg,     grid_name);
  add_required_field("inv_qc_relvar",      scalar3d_layout_mid, nondim,   grid_name);
  add_required_field("pseudo_density",     scalar3d_layout_mid, Pa,       grid_name);
  add_required_field("qv_prev_micro_step", scalar3d_layout_mid, Q,        grid_name);
  add_required_field("T_prev_micro_step",  scalar3d_layout_mid, K,        grid_name); 
  //
  add_computed_field("qv_prev_micro_step", scalar3d_layout_mid, Q,        grid_name);
  add_computed_field("T_prev_micro_step",  scalar3d_layout_mid, K,        grid_name);
  // Diagnostic Outputs: (all fields are just outputs w.r.t. P3)
  m_computed_fields.emplace("eff_radius_qc",     scalar3d_layout_mid, m,      grid_name);
  m_computed_fields.emplace("eff_radius_qi",     scalar3d_layout_mid, m,      grid_name);
  // History Only: (all fields are just outputs and are really only meant for I/O purposes)
  m_computed_fields.emplace("micro_liq_ice_exchange", scalar3d_layout_mid, nondim, grid_name);
  m_computed_fields.emplace("micro_vap_liq_exchange", scalar3d_layout_mid, nondim, grid_name);
  m_computed_fields.emplace("micro_vap_ice_exchange", scalar3d_layout_mid, nondim, grid_name);

}

// =========================================================================================
void P3Microphysics::initialize_impl (const util::TimeStamp& t0)
{
  m_current_ts = t0;

  // Initialize p3
  p3_init();

  // Initialize all of the structures that are passed to p3_main in run_impl.
  // Note: Some variables in the structures are not stored in the field manager.  For these
  //       variables a local view is constructed.
  const Int nk_pack = ekat::npack<Spack>(m_num_levs);
  auto pmid  = m_p3_fields_in["p_mid"].get_reshaped_view<const Pack**>();
  auto T_atm = m_p3_fields_out["T_mid"].get_reshaped_view<Pack**>();
  auto ast   = m_p3_fields_in["cldfrac_tot"].get_reshaped_view<const Pack**>();
  auto zi    = m_p3_fields_in["zi"].get_reshaped_view<const Pack**>();
  view_2d exner("exner",m_num_cols,nk_pack);
  view_2d th_atm("th_atm",m_num_cols,nk_pack);
  view_2d cld_frac_l("cld_frac_l",m_num_cols,nk_pack);
  view_2d cld_frac_i("cld_frac_i",m_num_cols,nk_pack);
  view_2d cld_frac_r("cld_frac_r",m_num_cols,nk_pack);
  view_2d dz("dz",m_num_cols,nk_pack);
  // -- Set values for the post-amble structure
  p3_preproc.set_variables(m_num_cols,nk_pack,pmid,T_atm,ast,zi,
                        exner, th_atm, cld_frac_l, cld_frac_i, cld_frac_r, dz);
  // --Prognostic State Variables:
  prog_state.qc     = m_p3_fields_out["qc"].get_reshaped_view<Pack**>();
  prog_state.nc     = m_p3_fields_out["nc"].get_reshaped_view<Pack**>();
  prog_state.qr     = m_p3_fields_out["qr"].get_reshaped_view<Pack**>();
  prog_state.nr     = m_p3_fields_out["nr"].get_reshaped_view<Pack**>();
  prog_state.qi     = m_p3_fields_out["qi"].get_reshaped_view<Pack**>();
  prog_state.qm     = m_p3_fields_out["qm"].get_reshaped_view<Pack**>();
  prog_state.ni     = m_p3_fields_out["ni"].get_reshaped_view<Pack**>();
  prog_state.bm     = m_p3_fields_out["bm"].get_reshaped_view<Pack**>();
  prog_state.qv     = m_p3_fields_out["qv"].get_reshaped_view<Pack**>();
  prog_state.th     = p3_preproc.th_atm;
  // --Diagnostic Input Variables:
  diag_inputs.nc_nuceat_tend  = m_p3_fields_in["nc_nuceat_tend"].get_reshaped_view<const Pack**>();
  diag_inputs.nccn            = m_p3_fields_in["nc_activated"].get_reshaped_view<const Pack**>();
  diag_inputs.ni_activated    = m_p3_fields_in["ni_activated"].get_reshaped_view<const Pack**>();
  diag_inputs.inv_qc_relvar   = m_p3_fields_in["inv_qc_relvar"].get_reshaped_view<const Pack**>();
  diag_inputs.pres            = m_p3_fields_in["p_mid"].get_reshaped_view<const Pack**>();
  diag_inputs.dpres           = m_p3_fields_in["pseudo_density"].get_reshaped_view<const Pack**>();
  auto qv_prev                = m_p3_fields_out["qv_prev_micro_step"].get_reshaped_view<Pack**>();
  diag_inputs.qv_prev         = qv_prev;
  auto t_prev                 = m_p3_fields_out["T_prev_micro_step"].get_reshaped_view<Pack**>();
  diag_inputs.t_prev          = t_prev;
  diag_inputs.cld_frac_l      = p3_preproc.cld_frac_l;
  diag_inputs.cld_frac_i      = p3_preproc.cld_frac_i;
  diag_inputs.cld_frac_r      = p3_preproc.cld_frac_r;
  diag_inputs.dz              = p3_preproc.dz;
  diag_inputs.exner           = p3_preproc.exner;
  // --Diagnostic Outputs
  view_1d precip_liq_surf("precip_liq_surf",m_num_cols);
  view_1d precip_ice_surf("precip_ice_surf",m_num_cols);
  view_2d qv2qi_depos_tend("qv2qi_depos_tend",m_num_cols,m_num_levs);
  view_2d rho_qi("rho_qi",m_num_cols,m_num_levs);
  view_2d precip_liq_flux("precip_liq_flux",m_num_cols,m_num_levs);
  view_2d precip_ice_flux("precip_ice_flux",m_num_cols,m_num_levs);

  diag_outputs.diag_eff_radius_qc = m_p3_fields_out["eff_radius_qc"].get_reshaped_view<Pack**>();
  diag_outputs.diag_eff_radius_qi = m_p3_fields_out["eff_radius_qi"].get_reshaped_view<Pack**>();

  diag_outputs.precip_liq_surf  = precip_liq_surf; 
  diag_outputs.precip_ice_surf  = precip_ice_surf; 
  diag_outputs.qv2qi_depos_tend = qv2qi_depos_tend; 
  diag_outputs.rho_qi           = rho_qi;
  diag_outputs.precip_liq_flux  = precip_liq_flux;
  diag_outputs.precip_ice_flux  = precip_ice_flux;
  // --Infrastructure
  // dt is passed as an argument to run_impl
  infrastructure.it  = 0;
  infrastructure.its = 0;
  infrastructure.ite = m_num_cols-1;
  infrastructure.kts = 0;
  infrastructure.kte = m_num_levs-1;
  infrastructure.predictNc = true;     // Hard-coded for now, TODO: make this a runtime option 
  infrastructure.prescribedCCN = true; // Hard-coded for now, TODO: make this a runtime option
  sview_2d col_location("col_location", m_num_cols, 3);
  infrastructure.col_location = col_location; // TODO: Initialize this here and now when P3 has access to lat/lon for each column.
  // --History Only
  history_only.liq_ice_exchange = m_p3_fields_out["micro_liq_ice_exchange"].get_reshaped_view<Pack**>();
  history_only.vap_liq_exchange = m_p3_fields_out["micro_vap_liq_exchange"].get_reshaped_view<Pack**>();
  history_only.vap_ice_exchange = m_p3_fields_out["micro_vap_ice_exchange"].get_reshaped_view<Pack**>();
  // Deprecated -- These are fields actively being deleted, but are still needed for F90 BFB tests.
  view_2d nevapr("nevapr",m_num_cols,nk_pack);
  view_2d qr_evap_tend("qr_evap_tend",m_num_cols,nk_pack);
  view_2d precip_total_tend("precip_total_tend",m_num_cols,nk_pack);
  view_2d mu_c("mu_c",m_num_cols,nk_pack);
  view_2d lamc("lamc",m_num_cols,nk_pack);
  deprecated.nevapr            = nevapr; 
  deprecated.qr_evap_tend      = qr_evap_tend; 
  deprecated.precip_total_tend = precip_total_tend; 
  deprecated.mu_c              = mu_c; 
  deprecated.lamc              = lamc; 
  // -- Set values for the post-amble structure
  p3_postproc.set_variables(m_num_cols,nk_pack,prog_state.th,p3_preproc.exner,T_atm,t_prev,prog_state.qv,qv_prev,
      diag_outputs.diag_eff_radius_qc,diag_outputs.diag_eff_radius_qi,ast);
}

// =========================================================================================
void P3Microphysics::run_impl (const Real dt)
{

  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_p3_fields_in) {
    it.second.sync_to_host();
  }
  for (auto& it : m_p3_fields_out) {
    it.second.sync_to_host();
  }

  // Copy outputs back to device
  // LB: why?!?
  for (auto& it : m_p3_fields_out) {
    it.second.sync_to_dev();
  }

  // Gather views needed to pre-process local variables.
  auto T_atm  = m_p3_fields_out["T_mid"].get_reshaped_view<Pack**>();
  auto ast    = m_p3_fields_in["cldfrac_tot"].get_reshaped_view<const Pack**>();
  auto zi     = m_p3_fields_in["zi"].get_reshaped_view<const Pack**>();
  auto pmid   = m_p3_fields_in["p_mid"].get_reshaped_view<const Pack**>();

  // Assign values to local arrays used by P3, these are now stored in p3_loc.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_preproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();

  // Update the variables in the p3 input structures with local values.

  infrastructure.dt = dt;
  infrastructure.it++;

  // Run p3 main
  P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
                                       history_only, deprecated, m_num_cols, m_num_levs);

  // Conduct the post-processing of the p3_main output.
  Kokkos::parallel_for(
    "p3_main_local_vals",
    Kokkos::RangePolicy<>(0,m_num_cols),
    p3_postproc
  ); // Kokkos::parallel_for(p3_main_local_vals)
  Kokkos::fence();

  // Get a copy of the current timestamp (at the beginning of the step) and
  // advance it, updating the p3 fields.
  auto ts = timestamp();
  ts += dt;
  for (auto& f : m_p3_fields_out) {
    f.second.get_header().get_tracking().update_time_stamp(ts);
  }

}

// =========================================================================================
void P3Microphysics::finalize_impl()
{
  // Do nothing
}

// =========================================================================================
void P3Microphysics::
register_fields (const std::map<std::string,std::shared_ptr<FieldManager<Real>>>& field_mgrs) const {
  std::set<ci_string> q_names =
    { "qv","qc","qr","qi","qm","nc","nr","ni","bm"};

  const auto& grid_name = m_p3_params.get<std::string>("Grid");
  auto& field_mgr = *field_mgrs.at(grid_name);
  for (const auto& fid : get_required_fields()) {
    const auto& name = fid.name();
    if (q_names.count(name)>0) {
      field_mgr.register_field<Pack>(fid,"TRACERS");
    } else {
      field_mgr.register_field<Pack>(fid);
    }
  }
  for (const auto& fid : get_computed_fields()) {
    const auto& name = fid.name();
    if (q_names.count(name)>0) {
      field_mgr.register_field<Pack>(fid,"TRACERS");
    } else {
      field_mgr.register_field<Pack>(fid);
    }
  }
}

void P3Microphysics::set_required_field_impl (const Field<const Real>& f) {

  const auto& name = f.get_header().get_identifier().name();
  m_p3_fields_in.emplace(name,f);
  m_p3_host_views_in[name] = f.get_view<Host>();
  m_raw_ptrs_in[name] = m_p3_host_views_in[name].data();

  // Add myself as customer to the field
  add_me_as_customer(f);
}

void P3Microphysics::set_computed_field_impl (const Field<      Real>& f) {

  const auto& name = f.get_header().get_identifier().name();
  m_p3_fields_out.emplace(name,f);
  m_p3_host_views_out[name] = f.get_view<Host>();
  m_raw_ptrs_out[name] = m_p3_host_views_out[name].data();

  // Add myself as provider for the field
  add_me_as_provider(f);
}

} // namespace scream
