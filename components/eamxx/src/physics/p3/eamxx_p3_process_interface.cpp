#include "share/property_checks/field_within_interval_check.hpp"
#include "share/property_checks/field_lower_bound_check.hpp"
#include "p3_functions.hpp"
#include "eamxx_p3_process_interface.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{

// =========================================================================================
P3Microphysics::P3Microphysics (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

// =========================================================================================
void P3Microphysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ekat::prefixes;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto nondim = Units::nondimensional();
  auto micron = micro*m;
  auto m2 = pow(m,2);

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

  // --Infrastructure
  // dt is passed as an argument to run_impl
  infrastructure.it  = 0;
  infrastructure.its = 0;
  infrastructure.ite = m_num_cols-1;
  infrastructure.kts = 0;
  infrastructure.kte = m_num_levs-1;
  infrastructure.predictNc = m_params.get<bool>("do_predict_nc",true);
  infrastructure.prescribedCCN = m_params.get<bool>("do_prescribed_ccn",true);

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout { {COL}, {m_num_cols} };

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  // Define fields needed in P3.
  // Note: p3_main is organized by a set of 5 structures, variables below are organized
  //       using the same approach to make it easier to follow.

  constexpr int ps = Pack::n;

  // These variables are needed by the interface, but not actually passed to p3_main.
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid_name, ps);

//should we use one pressure only, wet/full?
  add_field<Required>("p_mid",       scalar3d_layout_mid, Pa,     grid_name, ps);
  add_field<Required>("p_dry_mid",   scalar3d_layout_mid, Pa,     grid_name, ps);
  add_field<Updated> ("T_mid",       scalar3d_layout_mid, K,      grid_name, ps);  // T_mid is the only one of these variables that is also updated.

  // Prognostic State:  (all fields are both input and output)
  add_tracer<Updated>("qv", m_grid, kg/kg, ps);
  add_tracer<Updated>("qc", m_grid, kg/kg, ps);
  add_tracer<Updated>("qr", m_grid, kg/kg, ps);
  add_tracer<Updated>("qi", m_grid, kg/kg, ps);
  add_tracer<Updated>("qm", m_grid, kg/kg, ps);
  add_tracer<Updated>("nc", m_grid, 1/kg,  ps);
  add_tracer<Updated>("nr", m_grid, 1/kg,  ps);
  add_tracer<Updated>("ni", m_grid, 1/kg,  ps);
  add_tracer<Updated>("bm", m_grid, 1/kg,  ps);

  // Diagnostic Inputs: (only the X_prev fields are both input and output, all others are just inputs)
  add_field<Required>("nc_nuceat_tend",     scalar3d_layout_mid, 1/(kg*s), grid_name, ps);
  if (infrastructure.prescribedCCN) {
    add_field<Required>("nccn",               scalar3d_layout_mid, 1/kg,     grid_name, ps);
  }
  add_field<Required>("ni_activated",       scalar3d_layout_mid, 1/kg,         grid_name, ps);
  add_field<Required>("inv_qc_relvar",      scalar3d_layout_mid, pow(kg/kg,2), grid_name, ps);
  add_field<Required>("pseudo_density",     scalar3d_layout_mid, Pa,           grid_name, ps);
  add_field<Required>("pseudo_density_dry", scalar3d_layout_mid, Pa,           grid_name, ps);
  add_field<Updated> ("qv_prev_micro_step", scalar3d_layout_mid, kg/kg,        grid_name, ps);
  add_field<Updated> ("T_prev_micro_step",  scalar3d_layout_mid, K,            grid_name, ps);

  // Diagnostic Outputs: (all fields are just outputs w.r.t. P3)
  add_field<Updated>("precip_liq_surf_mass", scalar2d_layout,     kg/m2,     grid_name, "ACCUMULATED");
  add_field<Updated>("precip_ice_surf_mass", scalar2d_layout,     kg/m2,     grid_name, "ACCUMULATED");
  add_field<Computed>("eff_radius_qc",       scalar3d_layout_mid, micron,    grid_name, ps);
  add_field<Computed>("eff_radius_qi",       scalar3d_layout_mid, micron,    grid_name, ps);
  add_field<Computed>("eff_radius_qr",       scalar3d_layout_mid, micron,    grid_name, ps);
  add_field<Computed>("precip_total_tend",   scalar3d_layout_mid, kg/(kg*s), grid_name, ps);
  add_field<Computed>("nevapr",              scalar3d_layout_mid, kg/(kg*s), grid_name, ps);

  // History Only: (all fields are just outputs and are really only meant for I/O purposes)
  // TODO: These should be averaged over subcycle as well.  But there is no simple mechanism
  //       yet to reset these values at the beginning of the atmosphere timestep.  When this
  //       mechanism is developed we should add these variables to the accumulated variables.
  add_field<Computed>("micro_liq_ice_exchange", scalar3d_layout_mid, kg/kg,  grid_name, ps);
  add_field<Computed>("micro_vap_liq_exchange", scalar3d_layout_mid, kg/kg,  grid_name, ps);
  add_field<Computed>("micro_vap_ice_exchange", scalar3d_layout_mid, kg/kg,  grid_name, ps);
  add_field<Computed>("rainfrac",               scalar3d_layout_mid, nondim, grid_name, ps);

  // Boundary flux fields for energy and mass conservation checks
  if (has_column_conservation_check()) {
    add_field<Computed>("vapor_flux", scalar2d_layout, kg/m2/s, grid_name);
    add_field<Computed>("water_flux", scalar2d_layout, m/s,     grid_name);
    add_field<Computed>("ice_flux",   scalar2d_layout, m/s,     grid_name);
    add_field<Computed>("heat_flux",  scalar2d_layout, W/m2,    grid_name);
  }
}

// =========================================================================================
size_t P3Microphysics::requested_buffer_size_in_bytes() const
{
  const Int nk_pack    = ekat::npack<Spack>(m_num_levs);
  const Int nk_pack_p1 = ekat::npack<Spack>(m_num_levs+1);

  // Number of Reals needed by local views in the interface
  const size_t interface_request =
      // 1d view scalar, size (ncol)
      Buffer::num_1d_scalar*m_num_cols*sizeof(Real) +
      // 2d view packed, size (ncol, nlev_packs)
      Buffer::num_2d_vector*m_num_cols*nk_pack*sizeof(Spack) +
      Buffer::num_2dp1_vector*m_num_cols*nk_pack_p1*sizeof(Spack) +
      // 2d view scalar, size (ncol, 3)
      m_num_cols*3*sizeof(Real);

  // Number of Reals needed by the WorkspaceManager passed to p3_main
  const auto policy       = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nk_pack);
  const size_t wsm_request   = WSM::get_total_bytes_needed(nk_pack_p1, 52, policy);

  return interface_request + wsm_request;
}

// =========================================================================================
void P3Microphysics::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d scalar views
  using scalar_1d_view_t = decltype(m_buffer.precip_liq_surf_flux);
  scalar_1d_view_t* _1d_scalar_view_ptrs[Buffer::num_1d_scalar] = {
    &m_buffer.precip_liq_surf_flux, &m_buffer.precip_ice_surf_flux
  };
  for (int i=0; i<Buffer::num_1d_scalar; ++i) {
    *_1d_scalar_view_ptrs[i] = scalar_1d_view_t(mem, m_num_cols);
    mem += _1d_scalar_view_ptrs[i]->size();
  }

  // 2d scalar views
  m_buffer.col_location = decltype(m_buffer.col_location)(mem, m_num_cols, 3);
  mem += m_buffer.col_location.size();

  Spack* s_mem = reinterpret_cast<Spack*>(mem);

  // 2d packed views
  const Int nk_pack    = ekat::npack<Spack>(m_num_levs);
  const Int nk_pack_p1 = ekat::npack<Spack>(m_num_levs+1);

  using spack_2d_view_t = decltype(m_buffer.inv_exner);
  spack_2d_view_t* _2d_spack_mid_view_ptrs[Buffer::num_2d_vector] = {
    &m_buffer.inv_exner, &m_buffer.th_atm, &m_buffer.cld_frac_l, &m_buffer.cld_frac_i,
    &m_buffer.dz, &m_buffer.qv2qi_depos_tend, &m_buffer.rho_qi, &m_buffer.unused
#ifdef SCREAM_P3_SMALL_KERNELS
    , &m_buffer.mu_r, &m_buffer.T_atm, &m_buffer.lamr, &m_buffer.logn0r, &m_buffer.nu,
    &m_buffer.cdist, &m_buffer.cdist1, &m_buffer.cdistr, &m_buffer.inv_cld_frac_i,
    &m_buffer.inv_cld_frac_l, &m_buffer.inv_cld_frac_r, &m_buffer.qc_incld, &m_buffer.qr_incld,
    &m_buffer.qi_incld, &m_buffer.qm_incld, &m_buffer.nc_incld, &m_buffer.nr_incld,
    &m_buffer.ni_incld, &m_buffer.bm_incld, &m_buffer.inv_dz, &m_buffer.inv_rho, &m_buffer.ze_ice,
    &m_buffer.ze_rain, &m_buffer.prec, &m_buffer.rho, &m_buffer.rhofacr, &m_buffer.rhofaci,
    &m_buffer.acn, &m_buffer.qv_sat_l, &m_buffer.qv_sat_i, &m_buffer.sup, &m_buffer.qv_supersat_i,
    &m_buffer.tmparr2, &m_buffer.exner, &m_buffer.diag_equiv_reflectivity, &m_buffer.diag_vm_qi,
    &m_buffer.diag_diam_qi, &m_buffer.pratot, &m_buffer.prctot, &m_buffer.qtend_ignore,
    &m_buffer.ntend_ignore, &m_buffer.mu_c, &m_buffer.lamc, &m_buffer.qr_evap_tend, &m_buffer.v_qc,
    &m_buffer.v_nc, &m_buffer.flux_qx, &m_buffer.flux_nx, &m_buffer.v_qit, &m_buffer.v_nit,
    &m_buffer.flux_nit, &m_buffer.flux_bir, &m_buffer.flux_qir, &m_buffer.flux_qit, &m_buffer.v_qr,
    &m_buffer.v_nr
#endif
  };
  for (int i=0; i<Buffer::num_2d_vector; ++i) {
    *_2d_spack_mid_view_ptrs[i] = spack_2d_view_t(s_mem, m_num_cols, nk_pack);
    s_mem += _2d_spack_mid_view_ptrs[i]->size();
  }

  spack_2d_view_t* _2d_spack_int_view_ptrs[Buffer::num_2dp1_vector] = {
    &m_buffer.precip_liq_flux, &m_buffer.precip_ice_flux
  };
  for (int i=0; i<Buffer::num_2dp1_vector; ++i) {
    *_2d_spack_int_view_ptrs[i] = spack_2d_view_t(s_mem, m_num_cols, nk_pack_p1);
    s_mem += _2d_spack_int_view_ptrs[i]->size();
  }

  // WSM data
  m_buffer.wsm_data = s_mem;

  // Compute workspace manager size to check used memory
  // vs. requested memory
  const auto policy  = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nk_pack);
  const int wsm_size = WSM::get_total_bytes_needed(nk_pack_p1, 52, policy)/sizeof(Spack);
  s_mem += wsm_size;

  size_t used_mem = (reinterpret_cast<Real*>(s_mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for P3Microphysics.");
}

// =========================================================================================
void P3Microphysics::initialize_impl (const RunType /* run_type */)
{
  // Gather runtime options from file
  runtime_options.load_runtime_options_from_file(m_params);

  // Set property checks for fields in this process
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("T_mid"),m_grid,100.0,500.0,false);
  add_invariant_check<FieldWithinIntervalCheck>(get_field_out("qv"),m_grid,1e-13,0.2,true);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("qc"),m_grid,0.0,0.1,false);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("qi"),m_grid,0.0,0.1,false);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("qr"),m_grid,0.0,0.1,false);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("qm"),m_grid,0.0,0.1,false);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("nc"),m_grid,0.0,1.e11,false);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("nr"),m_grid,0.0,1.e10,false);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("ni"),m_grid,0.0,1.e10,false);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("bm"),m_grid,0.0,1.0,false);
  // The following checks on precip have been changed to lower bound checks, from an interval check.
  // TODO: Change back to interval check when it is possible to pass dt_atm for the check.  Because
  //       precip is now an accumulated mass, the upper bound is dependent on the timestep.
  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precip_liq_surf_mass"),m_grid,0.0,false);
  add_postcondition_check<FieldLowerBoundCheck>(get_field_out("precip_ice_surf_mass"),m_grid,0.0,false);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("eff_radius_qc"),m_grid,0.0,1.0e2,false);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("eff_radius_qi"),m_grid,0.0,5.0e3,false);
  add_postcondition_check<FieldWithinIntervalCheck>(get_field_out("eff_radius_qr"),m_grid,0.0,5.0e3,false);

  // Initialize p3
  P3F::p3_init(/* write_tables = */ false,
               this->get_comm().am_i_root());

  // Initialize all of the structures that are passed to p3_main in run_impl.
  // Note: Some variables in the structures are not stored in the field manager.  For these
  //       variables a local view is constructed.
  const Int nk_pack = ekat::npack<Spack>(m_num_levs);
  const Int nk_pack_p1 = ekat::npack<Spack>(m_num_levs+1);
  const  auto& pmid           = get_field_in("p_mid").get_view<const Pack**>();
  const  auto& pmid_dry       = get_field_in("p_dry_mid").get_view<const Pack**>();
  const  auto& pseudo_density = get_field_in("pseudo_density").get_view<const Pack**>();
  const  auto& pseudo_density_dry = get_field_in("pseudo_density_dry").get_view<const Pack**>();
  const  auto& T_atm          = get_field_out("T_mid").get_view<Pack**>();
  const  auto& cld_frac_t     = get_field_in("cldfrac_tot").get_view<const Pack**>();
  const  auto& qv             = get_field_out("qv").get_view<Pack**>();
  const  auto& qc             = get_field_out("qc").get_view<Pack**>();
  const  auto& nc             = get_field_out("nc").get_view<Pack**>();
  const  auto& qr             = get_field_out("qr").get_view<Pack**>();
  const  auto& nr             = get_field_out("nr").get_view<Pack**>();
  const  auto& qi             = get_field_out("qi").get_view<Pack**>();
  const  auto& qm             = get_field_out("qm").get_view<Pack**>();
  const  auto& ni             = get_field_out("ni").get_view<Pack**>();
  const  auto& bm             = get_field_out("bm").get_view<Pack**>();
  auto qv_prev                = get_field_out("qv_prev_micro_step").get_view<Pack**>();
  const auto& precip_liq_surf_mass = get_field_out("precip_liq_surf_mass").get_view<Real*>();
  const auto& precip_ice_surf_mass = get_field_out("precip_ice_surf_mass").get_view<Real*>();
  auto cld_frac_r             = get_field_out("rainfrac").get_view<Pack**>();

  // Alias local variables from temporary buffer
  auto inv_exner  = m_buffer.inv_exner;
  auto th_atm     = m_buffer.th_atm;
  auto cld_frac_l = m_buffer.cld_frac_l;
  auto cld_frac_i = m_buffer.cld_frac_i;
  auto dz         = m_buffer.dz;

  // -- Set values for the pre-amble structure
  p3_preproc.set_variables(m_num_cols,nk_pack,pmid,pmid_dry,pseudo_density,pseudo_density_dry,
                        T_atm,cld_frac_t,
                        qv, qc, nc, qr, nr, qi, qm, ni, bm, qv_prev,
                        inv_exner, th_atm, cld_frac_l, cld_frac_i, cld_frac_r, dz);
  // --Prognostic State Variables:
  prog_state.qc     = p3_preproc.qc;
  prog_state.nc     = p3_preproc.nc;
  prog_state.qr     = p3_preproc.qr;
  prog_state.nr     = p3_preproc.nr;
  prog_state.qi     = p3_preproc.qi;
  prog_state.qm     = p3_preproc.qm;
  prog_state.ni     = p3_preproc.ni;
  prog_state.bm     = p3_preproc.bm;
  prog_state.th     = p3_preproc.th_atm;
  prog_state.qv     = p3_preproc.qv;
  // --Diagnostic Input Variables:
  diag_inputs.nc_nuceat_tend  = get_field_in("nc_nuceat_tend").get_view<const Pack**>();
  if (infrastructure.prescribedCCN) {
    diag_inputs.nccn          = get_field_in("nccn").get_view<const Pack**>();
  } else {
    diag_inputs.nccn          = m_buffer.unused; //TODO set value of unused to something like 0.0 or nan as a layer of protection that it isn't being used.
  }
  diag_inputs.ni_activated    = get_field_in("ni_activated").get_view<const Pack**>();
  diag_inputs.inv_qc_relvar   = get_field_in("inv_qc_relvar").get_view<const Pack**>();

  // P3 will use dry pressure for dry qv_sat
  diag_inputs.pres            = get_field_in("p_dry_mid").get_view<const Pack**>();
  diag_inputs.dpres           = p3_preproc.pseudo_density_dry; //give dry density as input
  diag_inputs.qv_prev         = p3_preproc.qv_prev;
  auto t_prev                 = get_field_out("T_prev_micro_step").get_view<Pack**>();
  diag_inputs.t_prev          = t_prev;
  diag_inputs.cld_frac_l      = p3_preproc.cld_frac_l;
  diag_inputs.cld_frac_i      = p3_preproc.cld_frac_i;
  diag_inputs.cld_frac_r      = p3_preproc.cld_frac_r;
  diag_inputs.dz              = p3_preproc.dz;
  diag_inputs.inv_exner       = p3_preproc.inv_exner;
  // --Diagnostic Outputs
  diag_outputs.diag_eff_radius_qc = get_field_out("eff_radius_qc").get_view<Pack**>();
  diag_outputs.diag_eff_radius_qi = get_field_out("eff_radius_qi").get_view<Pack**>();
  diag_outputs.diag_eff_radius_qr = get_field_out("eff_radius_qr").get_view<Pack**>();
  diag_outputs.precip_total_tend  = get_field_out("precip_total_tend").get_view<Pack**>();
  diag_outputs.nevapr             = get_field_out("nevapr").get_view<Pack**>();

  diag_outputs.precip_liq_surf  = m_buffer.precip_liq_surf_flux;
  diag_outputs.precip_ice_surf  = m_buffer.precip_ice_surf_flux;
  diag_outputs.qv2qi_depos_tend = m_buffer.qv2qi_depos_tend;
  diag_outputs.rho_qi           = m_buffer.rho_qi;
  diag_outputs.precip_liq_flux  = m_buffer.precip_liq_flux;
  diag_outputs.precip_ice_flux  = m_buffer.precip_ice_flux;
  // -- Infrastructure, what is left to assign
  infrastructure.col_location = m_buffer.col_location; // TODO: Initialize this here and now when P3 has access to lat/lon for each column.
  // --History Only
  history_only.liq_ice_exchange = get_field_out("micro_liq_ice_exchange").get_view<Pack**>();
  history_only.vap_liq_exchange = get_field_out("micro_vap_liq_exchange").get_view<Pack**>();
  history_only.vap_ice_exchange = get_field_out("micro_vap_ice_exchange").get_view<Pack**>();
#ifdef SCREAM_P3_SMALL_KERNELS
  // Temporaries
  temporaries.mu_r                    = m_buffer.mu_r;
  temporaries.T_atm                   = m_buffer.T_atm;
  temporaries.lamr                    = m_buffer.lamr;
  temporaries.logn0r                  = m_buffer.logn0r;
  temporaries.nu                      = m_buffer.nu;
  temporaries.cdist                   = m_buffer.cdist;
  temporaries.cdist1                  = m_buffer.cdist1;
  temporaries.cdistr                  = m_buffer.cdistr;
  temporaries.inv_cld_frac_i          = m_buffer.inv_cld_frac_i;
  temporaries.inv_cld_frac_l          = m_buffer.inv_cld_frac_l;
  temporaries.inv_cld_frac_r          = m_buffer.inv_cld_frac_r;
  temporaries.qc_incld                = m_buffer.qc_incld;
  temporaries.qr_incld                = m_buffer.qr_incld;
  temporaries.qi_incld                = m_buffer.qi_incld;
  temporaries.qm_incld                = m_buffer.qm_incld;
  temporaries.nc_incld                = m_buffer.nc_incld;
  temporaries.nr_incld                = m_buffer.nr_incld;
  temporaries.ni_incld                = m_buffer.ni_incld;
  temporaries.bm_incld                = m_buffer.bm_incld;
  temporaries.inv_dz                  = m_buffer.inv_dz;
  temporaries.inv_rho                 = m_buffer.inv_rho;
  temporaries.ze_ice                  = m_buffer.ze_ice;
  temporaries.ze_rain                 = m_buffer.ze_rain;
  temporaries.prec                    = m_buffer.prec;
  temporaries.rho                     = m_buffer.rho;
  temporaries.rhofacr                 = m_buffer.rhofacr;
  temporaries.rhofaci                 = m_buffer.rhofaci;
  temporaries.acn                     = m_buffer.acn;
  temporaries.qv_sat_l                = m_buffer.qv_sat_l;
  temporaries.qv_sat_i                = m_buffer.qv_sat_i;
  temporaries.sup                     = m_buffer.sup;
  temporaries.qv_supersat_i           = m_buffer.qv_supersat_i;
  temporaries.tmparr2                 = m_buffer.tmparr2;
  temporaries.exner                   = m_buffer.exner;
  temporaries.diag_equiv_reflectivity = m_buffer.diag_equiv_reflectivity;
  temporaries.diag_vm_qi              = m_buffer.diag_vm_qi;
  temporaries.diag_diam_qi            = m_buffer.diag_diam_qi;
  temporaries.pratot                  = m_buffer.pratot;
  temporaries.prctot                  = m_buffer.prctot;
  temporaries.qtend_ignore            = m_buffer.qtend_ignore;
  temporaries.ntend_ignore            = m_buffer.ntend_ignore;
  temporaries.mu_c                    = m_buffer.mu_c;
  temporaries.lamc                    = m_buffer.lamc;
  temporaries.qr_evap_tend            = m_buffer.qr_evap_tend;
  temporaries.v_qc                    = m_buffer.v_qc;
  temporaries.v_nc                    = m_buffer.v_nc;
  temporaries.flux_qx                 = m_buffer.flux_qx;
  temporaries.flux_nx                 = m_buffer.flux_nx;
  temporaries.v_qit                   = m_buffer.v_qit;
  temporaries.v_nit                   = m_buffer.v_nit;
  temporaries.flux_nit                = m_buffer.flux_nit;
  temporaries.flux_bir                = m_buffer.flux_bir;
  temporaries.flux_qir                = m_buffer.flux_qir;
  temporaries.flux_qit                = m_buffer.flux_qit;
  temporaries.v_qr                    = m_buffer.v_qr;
  temporaries.v_nr                    = m_buffer.v_nr;
#endif

  // -- Set values for the post-amble structure
  p3_postproc.set_variables(m_num_cols,nk_pack,
                            prog_state.th,pmid,pmid_dry,T_atm,t_prev,
                            pseudo_density,pseudo_density_dry,
                            prog_state.qv, prog_state.qc, prog_state.nc, prog_state.qr,prog_state.nr,
                            prog_state.qi, prog_state.qm, prog_state.ni,prog_state.bm,qv_prev,
                            diag_outputs.diag_eff_radius_qc,diag_outputs.diag_eff_radius_qi,
                            diag_outputs.diag_eff_radius_qr,
                            diag_outputs.precip_liq_surf,diag_outputs.precip_ice_surf,
                            precip_liq_surf_mass,precip_ice_surf_mass);

  if (has_column_conservation_check()) {
    const auto& vapor_flux = get_field_out("vapor_flux").get_view<Real*>();
    const auto& water_flux = get_field_out("water_flux").get_view<Real*>();
    const auto& ice_flux   = get_field_out("ice_flux").get_view<Real*>();
    const auto& heat_flux  = get_field_out("heat_flux").get_view<Real*>();
    p3_postproc.set_mass_and_energy_fluxes(vapor_flux, water_flux, ice_flux, heat_flux);
  }

  // Load tables
  P3F::init_kokkos_ice_lookup_tables(lookup_tables.ice_table_vals, lookup_tables.collect_table_vals);
  P3F::init_kokkos_tables(lookup_tables.vn_table_vals, lookup_tables.vm_table_vals,
                          lookup_tables.revap_table_vals, lookup_tables.mu_r_table_vals,
                          lookup_tables.dnu_table_vals);

  // Setup WSM for internal local variables
  const auto policy = ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(m_num_cols, nk_pack);
  workspace_mgr.setup(m_buffer.wsm_data, nk_pack_p1, 52, policy);
}

// =========================================================================================
void P3Microphysics::finalize_impl()
{
  // Do nothing
}
// =========================================================================================
} // namespace scream
