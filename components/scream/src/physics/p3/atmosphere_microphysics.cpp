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
  auto micron = m / 1000000;

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

  constexpr int ps = Pack::n;

  // These variables are needed by the interface, but not actually passed to p3_main. 
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name, ps);
  add_field<Required>("p_mid",       scalar3d_layout_mid, Pa,     grid_name, ps);
  add_field<Required>("z_int",       scalar3d_layout_int, m,      grid_name, ps);
  add_field<Updated> ("T_mid",       scalar3d_layout_mid, K,      grid_name, ps);  // T_mid is the only one of these variables that is also updated.

  // Prognostic State:  (all fields are both input and output)
  add_field<Updated>("qv",     scalar3d_layout_mid, Q,    grid_name, "tracers", ps);
  add_field<Updated>("qc",     scalar3d_layout_mid, Q,    grid_name, "tracers", ps);
  add_field<Updated>("qr",     scalar3d_layout_mid, Q,    grid_name, "tracers", ps);
  add_field<Updated>("qi",     scalar3d_layout_mid, Q,    grid_name, "tracers", ps);
  add_field<Updated>("qm",     scalar3d_layout_mid, Q,    grid_name, "tracers", ps);
  add_field<Updated>("nc",     scalar3d_layout_mid, 1/kg, grid_name, "tracers", ps);
  add_field<Updated>("nr",     scalar3d_layout_mid, 1/kg, grid_name, "tracers", ps);
  add_field<Updated>("ni",     scalar3d_layout_mid, 1/kg, grid_name, "tracers", ps);
  add_field<Updated>("bm",     scalar3d_layout_mid, 1/kg, grid_name, "tracers", ps);

  // Diagnostic Inputs: (only the X_prev fields are both input and output, all others are just inputs)
  add_field<Required>("nc_nuceat_tend",     scalar3d_layout_mid, 1/(kg*s), grid_name, ps);
  add_field<Required>("nc_activated",       scalar3d_layout_mid, 1/kg,     grid_name, ps);
  add_field<Required>("ni_activated",       scalar3d_layout_mid, 1/kg,     grid_name, ps);
  add_field<Required>("inv_qc_relvar",      scalar3d_layout_mid, Q*Q,      grid_name, ps);
  add_field<Required>("pseudo_density",     scalar3d_layout_mid, Pa,       grid_name, ps);
  add_field<Updated> ("qv_prev_micro_step", scalar3d_layout_mid, Q,        grid_name, ps);
  add_field<Updated> ("T_prev_micro_step",  scalar3d_layout_mid, K,        grid_name, ps);

  // Diagnostic Outputs: (all fields are just outputs w.r.t. P3)
  add_field<Computed>("eff_radius_qc",      scalar3d_layout_mid, micron, grid_name, ps);
  add_field<Computed>("eff_radius_qi",      scalar3d_layout_mid, micron, grid_name, ps);

  // History Only: (all fields are just outputs and are really only meant for I/O purposes)
  add_field<Computed>("micro_liq_ice_exchange", scalar3d_layout_mid, Q, grid_name, ps);
  add_field<Computed>("micro_vap_liq_exchange", scalar3d_layout_mid, Q, grid_name, ps);
  add_field<Computed>("micro_vap_ice_exchange", scalar3d_layout_mid, Q, grid_name, ps);

}

// =========================================================================================
int P3Microphysics::requested_buffer_size_in_bytes() const
{
  const Int nk_pack = ekat::npack<Spack>(m_num_levs);

  // Number of Reals needed by local views in the interface
  const int interface_request =
      // 1d view scalar, size (ncol)
      Buffer::num_1d_scalar*m_num_cols*sizeof(Real) +
      // 2d view packed, size (ncol, nlev_packs)
      Buffer::num_2d_vector*m_num_cols*nk_pack*sizeof(Spack) +
      // 2d view scalar, size (ncol, 3)
      m_num_cols*3*sizeof(Real);

  // Number of Reals needed by the WorkspaceManager passed to p3_main
  const auto policy       = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nk_pack);
  const int wsm_request   = WSM::get_total_bytes_needed(nk_pack, 52, policy);

  return interface_request + wsm_request;
}

// =========================================================================================
void P3Microphysics::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d scalar views
  m_buffer.precip_liq_surf = decltype(m_buffer.precip_liq_surf)(mem, m_num_cols);
  mem += m_buffer.precip_liq_surf.size();
  m_buffer.precip_ice_surf = decltype(m_buffer.precip_ice_surf)(mem, m_num_cols);
  mem += m_buffer.precip_ice_surf.size();

  // 2d scalar views
  m_buffer.col_location = decltype(m_buffer.col_location)(mem, m_num_cols, 3);
  mem += m_buffer.col_location.size();

  Spack* s_mem = reinterpret_cast<Spack*>(mem);

  // 2d packed views
  const Int nk_pack = ekat::npack<Spack>(m_num_levs);

  m_buffer.inv_exner = decltype(m_buffer.inv_exner)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.inv_exner.size();
  m_buffer.th_atm = decltype(m_buffer.th_atm)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.th_atm.size();
  m_buffer.cld_frac_l = decltype(m_buffer.cld_frac_l)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.cld_frac_l.size();
  m_buffer.cld_frac_i = decltype(m_buffer.cld_frac_i)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.cld_frac_i.size();
  m_buffer.cld_frac_r = decltype(m_buffer.cld_frac_r)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.cld_frac_r.size();
  m_buffer.dz = decltype(m_buffer.dz)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.dz.size();
  m_buffer.qv2qi_depos_tend = decltype(m_buffer.qv2qi_depos_tend)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.qv2qi_depos_tend.size();
  m_buffer.rho_qi = decltype(m_buffer.rho_qi)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.rho_qi.size();
  m_buffer.precip_liq_flux = decltype(m_buffer.precip_liq_flux)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.precip_liq_flux.size();
  m_buffer.precip_ice_flux = decltype(m_buffer.precip_ice_flux)(s_mem, m_num_cols, nk_pack);
  s_mem += m_buffer.precip_ice_flux.size();

  // WSM data
  m_buffer.wsm_data = s_mem;

  // Compute workspace manager size to check used memory
  // vs. requested memory
  const auto policy  = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nk_pack);
  const int wsm_size = WSM::get_total_bytes_needed(nk_pack, 52, policy)/sizeof(Spack);
  s_mem += wsm_size;

  int used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for P3Microphysics.");
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
  const  auto& pmid           = m_p3_fields_in["p_mid"].get_reshaped_view<const Pack**>();
  const  auto& pseudo_density = m_p3_fields_in["pseudo_density"].get_reshaped_view<const Pack**>();
  const  auto& T_atm          = m_p3_fields_out["T_mid"].get_reshaped_view<Pack**>();
  const  auto& cld_frac_t     = m_p3_fields_in["cldfrac_tot"].get_reshaped_view<const Pack**>();
  const  auto& zi             = m_p3_fields_in["z_int"].get_reshaped_view<const Pack**>();
  const  auto& qv             = m_p3_fields_out["qv"].get_reshaped_view<Pack**>();

  // Alias local variables from temporary buffer
  auto inv_exner  = m_buffer.inv_exner;
  auto th_atm     = m_buffer.th_atm;
  auto cld_frac_l = m_buffer.cld_frac_l;
  auto cld_frac_i = m_buffer.cld_frac_i;
  auto cld_frac_r = m_buffer.cld_frac_r;
  auto dz         = m_buffer.dz;

  // -- Set values for the post-amble structure
  p3_preproc.set_variables(m_num_cols,nk_pack,pmid,pseudo_density,T_atm,cld_frac_t,qv,
                        inv_exner, th_atm, cld_frac_l, cld_frac_i, cld_frac_r, dz);
  // --Prognostic State Variables:
  prog_state.qc     = m_p3_fields_out["qc"].get_reshaped_view<Pack**>();
  prog_state.nc     = m_p3_fields_out["nc"].get_reshaped_view<Pack**>();
  prog_state.qr     = m_p3_fields_out["qr"].get_reshaped_view<Pack**>();
  prog_state.nr     = m_p3_fields_out["nr"].get_reshaped_view<Pack**>();
  prog_state.qi     = m_p3_fields_out["qi"].get_reshaped_view<Pack**>();
  prog_state.qm     = m_p3_fields_out["qm"].get_reshaped_view<Pack**>();
  prog_state.ni     = m_p3_fields_out["ni"].get_reshaped_view<Pack**>();
  prog_state.bm     = m_p3_fields_out["bm"].get_reshaped_view<Pack**>();
  prog_state.th     = p3_preproc.th_atm;
  prog_state.qv     = p3_preproc.qv;
  // --Diagnostic Input Variables:
  diag_inputs.nc_nuceat_tend  = m_p3_fields_in["nc_nuceat_tend"].get_reshaped_view<const Pack**>();
  diag_inputs.nccn            = m_p3_fields_in["nc_activated"].get_reshaped_view<const Pack**>();
  diag_inputs.ni_activated    = m_p3_fields_in["ni_activated"].get_reshaped_view<const Pack**>();
  diag_inputs.inv_qc_relvar   = m_p3_fields_in["inv_qc_relvar"].get_reshaped_view<const Pack**>();
  diag_inputs.pres            = m_p3_fields_in["p_mid"].get_reshaped_view<const Pack**>();
  diag_inputs.dpres           = p3_preproc.pseudo_density;
  auto qv_prev                = m_p3_fields_out["qv_prev_micro_step"].get_reshaped_view<Pack**>();
  diag_inputs.qv_prev         = qv_prev;
  auto t_prev                 = m_p3_fields_out["T_prev_micro_step"].get_reshaped_view<Pack**>();
  diag_inputs.t_prev          = t_prev;
  diag_inputs.cld_frac_l      = p3_preproc.cld_frac_l;
  diag_inputs.cld_frac_i      = p3_preproc.cld_frac_i;
  diag_inputs.cld_frac_r      = p3_preproc.cld_frac_r;
  diag_inputs.dz              = p3_preproc.dz;
  diag_inputs.inv_exner       = p3_preproc.inv_exner;
  // --Diagnostic Outputs
  diag_outputs.diag_eff_radius_qc = m_p3_fields_out["eff_radius_qc"].get_reshaped_view<Pack**>();
  diag_outputs.diag_eff_radius_qi = m_p3_fields_out["eff_radius_qi"].get_reshaped_view<Pack**>();

  diag_outputs.precip_liq_surf  = m_buffer.precip_liq_surf;
  diag_outputs.precip_ice_surf  = m_buffer.precip_ice_surf;
  diag_outputs.qv2qi_depos_tend = m_buffer.qv2qi_depos_tend;
  diag_outputs.rho_qi           = m_buffer.rho_qi;
  diag_outputs.precip_liq_flux  = m_buffer.precip_liq_flux;
  diag_outputs.precip_ice_flux  = m_buffer.precip_ice_flux;
  // --Infrastructure
  // dt is passed as an argument to run_impl
  infrastructure.it  = 0;
  infrastructure.its = 0;
  infrastructure.ite = m_num_cols-1;
  infrastructure.kts = 0;
  infrastructure.kte = m_num_levs-1;
  infrastructure.predictNc = true;     // Hard-coded for now, TODO: make this a runtime option 
  infrastructure.prescribedCCN = true; // Hard-coded for now, TODO: make this a runtime option
  infrastructure.col_location = m_buffer.col_location; // TODO: Initialize this here and now when P3 has access to lat/lon for each column.
  // --History Only
  history_only.liq_ice_exchange = m_p3_fields_out["micro_liq_ice_exchange"].get_reshaped_view<Pack**>();
  history_only.vap_liq_exchange = m_p3_fields_out["micro_vap_liq_exchange"].get_reshaped_view<Pack**>();
  history_only.vap_ice_exchange = m_p3_fields_out["micro_vap_ice_exchange"].get_reshaped_view<Pack**>();
  // -- Set values for the post-amble structure
  p3_postproc.set_variables(m_num_cols,nk_pack,prog_state.th,pmid,T_atm,t_prev,prog_state.qv,qv_prev,
      diag_outputs.diag_eff_radius_qc,diag_outputs.diag_eff_radius_qi);
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

  // WorkspaceManager for internal local variables
  const Int nk_pack = ekat::npack<Spack>(m_num_levs);
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nk_pack);
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr(m_buffer.wsm_data, nk_pack, 52, policy);

  // Run p3 main
  P3F::p3_main(prog_state, diag_inputs, diag_outputs, infrastructure,
               history_only, workspace_mgr, m_num_cols, m_num_levs);

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
