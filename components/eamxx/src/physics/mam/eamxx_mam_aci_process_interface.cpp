#include <physics/mam/eamxx_mam_aci_process_interface.hpp>

// ACI functions are stored in the following hpp file
#include <physics/mam/eamxx_mam_aci_functions.hpp>

// For EKAT units package
#include <physics/mam/physical_limits.hpp>

#include "ekat/util/ekat_units.hpp"
/*
-----------------------------------------------------------------
NOTES:
1. w_variance assumes that we are using SE dycore. If we use EUL
dycore, we need to get its value from previous dynamic time step

2. w_variance and eddy_diff_heat are at midpoints in EAMxx, we
derived their interface values using boundary conditions to be
the top and botton values of the mid point arrays. We assume that
this assumption should not cause any issues.

FUTURE WORK:
1. MAM4xx submodule should point to MAM4xx main branch
2. Link hetrozenous freezing outputs to microphysics
3. Add postcondition and invariant checks
5. Resolve comments about top_lev
6. Avoid using c-style static arrays in ptend and other arrays
7. Use std::string rather than c-strings
8. Remove a Kokkos:fence and combine two kernels while computing w_sec_int_
9. Fix double counting of tracer advection by modifying SHOC.
10. A git issue for computing top_lev, moving liq_cldfrac in ACI and using TKE
directly
14. Merge kernels so that we are not calling one from another
16. improve the way qqcw is populated
17.delete fence mentioned by Luca
-----------------------------------------------------------------
*/

namespace scream {
MAMAci::MAMAci(const ekat::Comm &comm, const ekat::ParameterList &params)
    : MAMGenericInterface(comm, params) {
  // Asserts for the runtime or namelist options
  EKAT_REQUIRE_MSG(m_params.isParameter("wsubmin"),
                   "ERROR: wsubmin is missing from mam_aci parameter list.");
  EKAT_REQUIRE_MSG(m_params.isParameter("enable_aero_vertical_mix"),
                   "ERROR: enable_aero_vertical_mix is missing from mam_aci "
                   "parameter list.");
  EKAT_REQUIRE_MSG(
      m_params.isParameter("top_level_mam4xx"),
      "ERROR: top_level_mam4xx is missing from mam_aci parameter list.");
  check_fields_intervals_ =
      m_params.get<bool>("create_fields_interval_checks", false);
}

// ================================================================
//  SET_GRIDS
// ================================================================
void MAMAci::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  // set grid for all the inputs and outputs
  // use physics grid
  grid_ = grids_manager->get_grid("Physics");

  // Name of the grid
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 2D (2d horiz) variable
  const FieldLayout scalar2d = grid_->get_2d_scalar_layout();

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  const FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);

  using namespace ekat::units;
  constexpr auto n_unit = 1 / kg;  // units of number mixing ratios of tracers

  constexpr auto nondim = ekat::units::Units::nondimensional();

  add_tracers_wet_atm();
  add_fields_dry_atm();

  constexpr auto m2 = pow(m, 2);
  constexpr auto s2 = pow(s, 2);

  // NOTE: w_variance im microp_aero.F90 in EAM is at "itim_old" dynamics time
  // step. Since, we are using SE dycore, itim_old is 1 which is equivalent to
  // the current time step. For other dycores (such as EUL), it may be different
  // and we might need to revisit this.

  // Vertical velocity variance at midpoints
  add_field<Required>("w_variance", scalar3d_mid, m2 / s2, grid_name);

  // NOTE: "cldfrac_liq" is updated in SHOC. "cldfrac_liq" in C++ code is
  // equivalent to "alst" in the shoc_intr.F90. In the C++ code, it is used as
  // "shoc_cldfrac" and in the F90 code it is called "cloud_frac"

  // Liquid stratiform cloud fraction at midpoints
  add_field<Required>("cldfrac_liq", scalar3d_mid, nondim, grid_name);

  // Previous value of liquid stratiform cloud fraction at midpoints
  add_field<Required>("cldfrac_liq_prev", scalar3d_mid, nondim, grid_name);

  // Eddy diffusivity for heat
  add_field<Required>("eddy_diff_heat", scalar3d_mid, m2 / s, grid_name);

  // Number of modes
  constexpr int nmodes = mam4::AeroConfig::num_modes();

  // layout for 3D (ncol, nmodes, nlevs)
  FieldLayout scalar3d_mid_nmodes =
      grid_->get_3d_vector_layout(true, nmodes, "nmodes");

  // dry diameter of aerosols [m]
  add_field<Required>("dgnum", scalar3d_mid_nmodes, m, grid_name);

  // ========================================================================
  // Output from this whole process
  // ========================================================================

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  // add tracers, e.g., num_a1, soa_a1
  add_tracers_interstitial_aerosol();
  // add tracer gases, e.g., O3
  add_tracers_gases();
  // add fields e.g., num_c1, soa_c1
  add_fields_cloudborne_aerosol();

  // ------------------------------------------------------------------------
  // Output from ice nucleation process
  // ------------------------------------------------------------------------

  // number of activated aerosol for ice nucleation[#/kg]
  add_field<Computed>("ni_activated", scalar3d_mid, n_unit, grid_name);

  // ------------------------------------------------------------------------
  // Output from droplet activation process (dropmixnuc)
  // ------------------------------------------------------------------------

  // tendency in droplet number mixing ratio [#/kg/s]
  add_field<Computed>("nc_nuceat_tend", scalar3d_mid, n_unit / s, grid_name);

  // FIXME: [TEMPORARY]droplet number mixing ratio source tendency [#/kg/s]
  add_field<Computed>("nsource", scalar3d_mid, n_unit / s, grid_name);

  // FIXME: [TEMPORARY]droplet number mixing ratio tendency due to mixing
  // [#/kg/s]
  add_field<Computed>("ndropmix", scalar3d_mid, n_unit / s, grid_name);

  // FIXME: [TEMPORARY]droplet number as seen by ACI [#/kg]
  add_field<Computed>("nc_inp_to_aci", scalar3d_mid, n_unit / s, grid_name);
  constexpr auto cm_tmp = m / 100;         // FIXME: [TEMPORARY] remove this
  constexpr auto cm3    = pow(cm_tmp, 3);  // FIXME: [TEMPORARY] remove this
  // FIXME: [TEMPORARY] remove the following ccn outputs
  add_field<Computed>("ccn_0p02", scalar3d_mid, cm3, grid_name);
  add_field<Computed>("ccn_0p05", scalar3d_mid, cm3, grid_name);
  add_field<Computed>("ccn_0p1", scalar3d_mid, cm3, grid_name);
  add_field<Computed>("ccn_0p2", scalar3d_mid, cm3, grid_name);
  add_field<Computed>("ccn_0p5", scalar3d_mid, cm3, grid_name);
  add_field<Computed>("ccn_1p0", scalar3d_mid, cm3, grid_name);

  // ------------------------------------------------------------------------
  // Output from hetrozenous freezing
  // ------------------------------------------------------------------------

  constexpr auto cm = m / 100;

  // units of number mixing ratios of tracers
  constexpr auto frz_unit = 1 / (cm * cm * cm * s);
  //  heterogeneous freezing by immersion nucleation [cm^-3 s^-1]
  add_field<Computed>("hetfrz_immersion_nucleation_tend", scalar3d_mid,
                      frz_unit, grid_name);

  // heterogeneous freezing by contact nucleation [cm^-3 s^-1]
  add_field<Computed>("hetfrz_contact_nucleation_tend", scalar3d_mid, frz_unit,
                      grid_name);

  // heterogeneous freezing by deposition nucleation [cm^-3 s^-1]
  add_field<Computed>("hetfrz_deposition_nucleation_tend", scalar3d_mid,
                      frz_unit, grid_name);
}  // function set_grids ends

// ================================================================
//  INIT_BUFFERS
// ================================================================

void MAMAci::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");
  buffer_.set_num_scratch(num_2d_scratch_);
  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(
      used_mem == requested_buffer_size_in_bytes(),
      "Error! Used memory != requested memory for MAMMicrophysics.");
}  // function init_buffers ends

// ================================================================
//  INITIALIZE_IMPL
// ================================================================
void MAMAci::initialize_impl(const RunType run_type) {
  // ------------------------------------------------------------------------
  // ## Runtime options
  // ------------------------------------------------------------------------
  // Check the interval values for the following fields used by this interface.
  // NOTE: We do not include aerosol and gas species, e.g., soa_a1, num_a1,
  // because we automatically added these fields.
  const std::map<std::string, std::pair<Real, Real>> ranges_aci = {
      // aci compute
      {"ni_activated", {-1e100, 1e100}},                     // FIXME
      {"nc_nuceat_tend", {-1e100, 1e100}},                   // FIXME
      {"nsource", {-1e10, 1e10}},                            // FIXME
      {"ndropmix", {-1e10, 1e10}},                           // FIXME
      {"nc_inp_to_aci", {-1e10, 1e10}},                      // FIXME
      {"ccn_0p02", {-1e10, 1e10}},                           // FIXME
      {"ccn_0p05", {-1e10, 1e10}},                           // FIXME
      {"ccn_0p1", {-1e10, 1e10}},                            // FIXME
      {"ccn_0p2", {-1e10, 1e10}},                            // FIXME
      {"ccn_0p5", {-1e10, 1e10}},                            // FIXME
      {"ccn_1p0", {-1e10, 1e10}},                            // FIXME
      {"hetfrz_immersion_nucleation_tend", {-1e10, 1e10}},   // FIXME
      {"hetfrz_contact_nucleation_tend", {-1e10, 1e10}},     // FIXME
      {"hetfrz_deposition_nucleation_tend", {-1e10, 1e10}},  // FIXME
      // aci required
      {"w_variance", {-1e10, 1e10}},        // FIXME
      {"cldfrac_liq", {-1e10, 1e10}},       // FIXME
      {"cldfrac_liq_prev", {-1e10, 1e10}},  // FIXME
      {"eddy_diff_heat", {-1e10, 1e10}},    // FIXME
      {"dgnum", {-1e10, 1e10}},             // FIXME
  };

  set_ranges_process(ranges_aci);
  add_interval_checks();
  wsubmin_                  = m_params.get<double>("wsubmin");
  enable_aero_vertical_mix_ = m_params.get<bool>("enable_aero_vertical_mix");
  top_lev_                  = m_params.get<int>("top_level_mam4xx");

  // ------------------------------------------------------------------------
  // Input fields read in from IC file, namelist or other processes
  // ------------------------------------------------------------------------
  w_sec_mid_    = get_field_in("w_variance").get_view<const Real **>();
  dgnum_        = get_field_in("dgnum").get_view<const Real ***>();
  liqcldf_      = get_field_in("cldfrac_liq").get_view<const Real **>();
  liqcldf_prev_ = get_field_in("cldfrac_liq_prev").get_view<const Real **>();
  kvh_mid_      = get_field_in("eddy_diff_heat").get_view<const Real **>();

  populate_wet_atm(wet_atm_);
  populate_dry_atm(dry_atm_, buffer_);

  // ------------------------------------------------------------------------
  // Output fields to be used by other processes
  // ------------------------------------------------------------------------
  // ice nucleation output
  naai_ = get_field_out("ni_activated").get_view<Real **>();

  // droplet activation output
  tendnd_ = get_field_out("nc_nuceat_tend").get_view<Real **>();

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  // It populates wet_aero struct (wet_aero_) with:
  // interstitial aerosol, e.g., soa_a_1
  populate_interstitial_wet_aero(wet_aero_);
  // gases, e.g., O3
  populate_gases_wet_aero(wet_aero_);
  // cloudborne aerosol, e.g., soa_c_1
  populate_cloudborne_wet_aero(wet_aero_);
  // It populates dry_aero struct (dry_aero_) with:
  // interstitial aerosol, e.g., soa_a_1
  populate_interstitial_dry_aero(dry_aero_, buffer_);
  // gases, e.g., O3
  populate_gases_dry_aero(dry_aero_, buffer_);
  // cloudborne aerosol, e.g., soa_c_1
  populate_cloudborne_dry_aero(dry_aero_, buffer_);

  // hetrozenous freezing outputs
  hetfrz_immersion_nucleation_tend_ =
      get_field_out("hetfrz_immersion_nucleation_tend").get_view<Real **>();
  hetfrz_contact_nucleation_tend_ =
      get_field_out("hetfrz_contact_nucleation_tend").get_view<Real **>();
  hetfrz_deposition_nucleation_tend_ =
      get_field_out("hetfrz_deposition_nucleation_tend").get_view<Real **>();

  //---------------------------------------------------------------------------------
  // Allocate memory for the class members
  // (Kokkos::resize only works on host to allocates memory)
  //---------------------------------------------------------------------------------

  set_field_w_scratch_buffer(rho_, buffer_, true);
  set_field_w_scratch_buffer(w0_, buffer_, true);
  Kokkos::resize(tke_, ncol_, nlev_ + 1);
  set_field_w_scratch_buffer(wsub_, buffer_, true);
  set_field_w_scratch_buffer(wsubice_, buffer_, true);
  set_field_w_scratch_buffer(wsig_, buffer_, true);
  set_field_w_scratch_buffer(w2_, buffer_, true);
  set_field_w_scratch_buffer(cloud_frac_, buffer_, true);
  set_field_w_scratch_buffer(cloud_frac_prev_, buffer_, true);
  set_field_w_scratch_buffer(aitken_dry_dia_, buffer_, true);
  set_field_w_scratch_buffer(rpdel_, buffer_, true);

  //---------------------------------------------------------------------------------
  // Diagnotics variables from the ice nucleation scheme
  //---------------------------------------------------------------------------------

  // number conc of ice nuclei due to heterogeneous freezing [1/m3]
  set_field_w_scratch_buffer(nihf_, buffer_, true);

  // number conc of ice nuclei due to immersionfreezing (hetero nuc) [1/m3]
  set_field_w_scratch_buffer(niim_, buffer_, true);

  // number conc of ice nuclei due to deposition nucleation (hetero nuc)[1/m3]
  set_field_w_scratch_buffer(nidep_, buffer_, true);

  // number conc of ice nuclei due to meyers deposition [1/m3]
  set_field_w_scratch_buffer(nimey_, buffer_, true);

  // number of activated aerosol for ice nucleation(homogeneous frz only)[#/kg]
  set_field_w_scratch_buffer(naai_hom_, buffer_, true);

  //---------------------------------------------------------------------------------
  // Diagnotics variables from the droplet activation scheme
  //---------------------------------------------------------------------------------

  // activation fraction for aerosol number [fraction]
  const int num_aero_modes = mam_coupling::num_aero_modes();
  Kokkos::resize(factnum_, ncol_, num_aero_modes, nlev_);

  // cloud droplet number mixing ratio [#/kg]
  set_field_w_scratch_buffer(qcld_, buffer_, true);

  // number conc of aerosols activated at supersat [#/m^3]
  // NOTE:  activation fraction fluxes are defined as
  // fluxn = [flux of activated aero. number into cloud[#/m^2/s]]
  //        / [aero. number conc. in updraft, just below cloudbase [#/m^3]]
  Kokkos::resize(ccn_, ncol_, nlev_, mam4::ndrop::psat);

  // column-integrated droplet number [#/m2]
  set_field_w_scratch_buffer(ndropcol_, buffer_, true);

  // droplet number mixing ratio tendency due to mixing [#/kg/s]
  // Kokkos::resize(ndropmix_, ncol_, nlev_);
  // Temporarily output ndropmix_
  ndropmix_ = get_field_out("ndropmix").get_view<Real **>();

  // droplet number mixing ratio source tendency [#/kg/s]
  // Kokkos::resize(nsource_, ncol_, nlev_);
  // Temporarily output nsource_
  nsource_ = get_field_out("nsource").get_view<Real **>();

  // Temporarily output nc_inp_to_aci_
  nc_inp_to_aci_ = get_field_out("nc_inp_to_aci").get_view<Real **>();

  // FIXME: [TEMPORARY] remove the following ccn outputs
  ccn_0p02_ = get_field_out("ccn_0p02").get_view<Real **>();
  ccn_0p05_ = get_field_out("ccn_0p05").get_view<Real **>();
  ccn_0p1_  = get_field_out("ccn_0p1").get_view<Real **>();
  ccn_0p2_  = get_field_out("ccn_0p2").get_view<Real **>();
  ccn_0p5_  = get_field_out("ccn_0p5").get_view<Real **>();
  ccn_1p0_  = get_field_out("ccn_1p0").get_view<Real **>();

  // subgrid vertical velocity [m/s]
  set_field_w_scratch_buffer(wtke_, buffer_, true);

  for(int i = 0; i < dropmix_scratch_; ++i) {
    set_field_w_scratch_buffer(dropmixnuc_scratch_mem_[i], buffer_, true);
  }
  for(int i = 0; i < mam4::ndrop::ncnst_tot; ++i) {
    // column tendency for diagnostic output
    set_field_w_scratch_buffer(coltend_[i], buffer_, true);
    // column tendency
    set_field_w_scratch_buffer(coltend_cw_[i], buffer_, true);
  }
  for(int i = 0; i < mam4::aero_model::pcnst; ++i) {
    set_field_w_scratch_buffer(ptend_q_[i], buffer_, true);
  }
  for(int i = 0; i < mam4::ndrop::pver; ++i) {
    for(int j = 0; j < 2; ++j) {
      Kokkos::resize(raercol_cw_[i][j], ncol_, mam4::ndrop::ncnst_tot);
      Kokkos::resize(raercol_[i][j], ncol_, mam4::ndrop::ncnst_tot);
    }
  }

  // nact : fractional aero. number activation rate [/s]
  Kokkos::resize(nact_, ncol_, nlev_, mam_coupling::num_aero_modes());

  // mact : fractional aero. mass activation rate [/s]
  Kokkos::resize(mact_, ncol_, nlev_, mam_coupling::num_aero_modes());

  // Eddy diffusivity of heat at the interfaces
  Kokkos::resize(kvh_int_, ncol_, nlev_ + 1);

  // Vertical velocity variance at the interfaces
  Kokkos::resize(w_sec_int_, ncol_, nlev_ + 1);
  // Allocate work arrays
  for(int icnst = 0; icnst < mam4::ndrop::ncnst_tot; ++icnst) {
    qqcw_fld_work_[icnst] = view_2d("qqcw_fld_work_", ncol_, nlev_);
    set_field_w_scratch_buffer(qqcw_fld_work_[icnst], buffer_, true);
  }
  state_q_work_ =
      view_3d("state_q_work_", ncol_, nlev_, mam4::aero_model::pcnst);

  //---------------------------------------------------------------------------------
  // Diagnotics variables from the hetrozenous ice nucleation scheme
  //---------------------------------------------------------------------------------

  for(int i = 0; i < hetro_scratch_; ++i)
    set_field_w_scratch_buffer(diagnostic_scratch_[i], buffer_, true);

  //---------------------------------------------------------------------------------
  // Initialize the processes
  //---------------------------------------------------------------------------------

  mam4::AeroConfig aero_config;
  // configure the nucleation parameterization
  mam4::NucleateIce::Config nucleate_ice_config;
  nucleate_ice_.init(aero_config, nucleate_ice_config);

  // configure the heterogeneous freezing parameterization
  mam4::Hetfrz::Config hetfrz_config;
  hetfrz_.init(aero_config, hetfrz_config);

  //---------------------------------------------------------------------------------
  // Setup preprocessing and post processing
  //---------------------------------------------------------------------------------
  // set up our preprocess  and postprocess functors
  preprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_,
                         dry_aero_);

}  // end function initialize_impl

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMAci::run_impl(const double dt) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of local derivied
  // quantities
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  haero::ThreadTeamPolicy team_policy(ncol_, Kokkos::AUTO);

  // FIXME: Temporary assignment of nc
  mam_coupling::copy_view_lev_slice(team_policy, wet_atm_.nc, nlev_,  // inputs
                                    nc_inp_to_aci_);                  // output

  compute_w0_and_rho(team_policy, dry_atm_, top_lev_, nlev_,
                     // output
                     w0_, rho_);

  compute_tke_at_interfaces(team_policy, w_sec_mid_, dry_atm_.dz, nlev_,
                            w_sec_int_,
                            // output
                            tke_);

  Kokkos::fence();  // wait for tke_ to be computed.
  compute_subgrid_scale_velocities(team_policy, tke_, wsubmin_, top_lev_, nlev_,
                                   // output
                                   wsub_, wsubice_, wsig_);

  // We need dry diameter for only aitken mode
  Kokkos::deep_copy(
      aitken_dry_dia_,
      ekat::subview_1(dgnum_, static_cast<int>(mam4::ModeIndex::Aitken)));

  Kokkos::fence();  // wait for aitken_dry_dia_ to be copied.

  //  Compute Ice nucleation
  //  NOTE: The Fortran version uses "ast" for cloud fraction which is
  //  equivalent to "cldfrac_tot" in FM. It is part of the "dry_atm_" struct
  compute_nucleate_ice_tendencies(
      nucleate_ice_, team_policy, dry_atm_, dry_aero_, wsubice_,
      aitken_dry_dia_, nlev_, dt,
      // output
      nihf_, niim_, nidep_, nimey_, naai_hom_,
      // ## output to be used by the other processes ##
      naai_);

  // Compute cloud fractions based on cloud threshold
  store_liquid_cloud_fraction(team_policy, dry_atm_, liqcldf_, liqcldf_prev_,
                              top_lev_, nlev_,
                              // output
                              cloud_frac_, cloud_frac_prev_);

  mam_coupling::compute_recipical_pseudo_density(team_policy, dry_atm_.p_del,
                                                 nlev_,
                                                 // output
                                                 rpdel_);

  Kokkos::fence();  // wait for rpdel_ to be computed.

  //  Compute activated CCN number tendency (tendnd_) and updated
  //  cloud borne aerosols (stored in a work array) and interstitial
  //  aerosols tendencies
  call_function_dropmixnuc(
      team_policy, dt, dry_atm_, rpdel_, kvh_mid_, kvh_int_, wsub_, cloud_frac_,
      cloud_frac_prev_, dry_aero_, nlev_, enable_aero_vertical_mix_,
      // output
      coltend_, coltend_cw_, qcld_, ndropcol_, ndropmix_, nsource_, wtke_, ccn_,
      // ## output to be used by the other processes ##
      qqcw_fld_work_, ptend_q_, factnum_, tendnd_,
      // work arrays
      raercol_cw_, raercol_, state_q_work_, nact_, mact_,
      dropmixnuc_scratch_mem_);
  Kokkos::fence();  // wait for ptend_q_ to be computed.

  Kokkos::deep_copy(ccn_0p02_,
                    Kokkos::subview(ccn_, Kokkos::ALL(), Kokkos::ALL(), 0));
  Kokkos::deep_copy(ccn_0p05_,
                    Kokkos::subview(ccn_, Kokkos::ALL(), Kokkos::ALL(), 1));
  Kokkos::deep_copy(ccn_0p1_,
                    Kokkos::subview(ccn_, Kokkos::ALL(), Kokkos::ALL(), 2));
  Kokkos::deep_copy(ccn_0p2_,
                    Kokkos::subview(ccn_, Kokkos::ALL(), Kokkos::ALL(), 3));
  Kokkos::deep_copy(ccn_0p5_,
                    Kokkos::subview(ccn_, Kokkos::ALL(), Kokkos::ALL(), 4));
  Kokkos::deep_copy(ccn_1p0_,
                    Kokkos::subview(ccn_, Kokkos::ALL(), Kokkos::ALL(), 5));

  //---------------------------------------------------------------------------
  //  NOTE: DO NOT UPDATE cloud borne aerosols using the qqcw_fld_work_ array
  //  at this point as heterozenous freezing needs to use cloud borne aerosols
  //  before they are changed by the droplet activation (dropmixnuc) process.
  //---------------------------------------------------------------------------

  // Compute hetrozenous freezing
  call_hetfrz_compute_tendencies(
      team_policy, hetfrz_, dry_atm_, dry_aero_, factnum_, dt, nlev_,
      // ## output to be used by the other processes ##
      hetfrz_immersion_nucleation_tend_, hetfrz_contact_nucleation_tend_,
      hetfrz_deposition_nucleation_tend_,
      // work arrays
      diagnostic_scratch_);

  //---------------------------------------------------------------
  // Now update interstitial and cloud borne aerosols
  //---------------------------------------------------------------

  // Update cloud borne aerosols
  update_cloud_borne_aerosols(qqcw_fld_work_, nlev_,
                              // output
                              dry_aero_);

  // Update interstitial aerosols using tendencies
  update_interstitial_aerosols(team_policy, ptend_q_, nlev_, dt,
                               // output
                               dry_aero_);

  // call post processing to convert dry mixing ratios to wet mixing ratios

  post_process(wet_aero_, dry_aero_, dry_atm_);
  Kokkos::fence();  // wait before returning to calling function
}

}  // namespace scream
