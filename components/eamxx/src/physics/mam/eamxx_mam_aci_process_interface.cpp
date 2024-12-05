#include <physics/mam/eamxx_mam_aci_process_interface.hpp>

// ACI functions are stored in the following hpp file
#include <physics/mam/eamxx_mam_aci_functions.hpp>

// For EKAT units package
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
    : AtmosphereProcess(comm, params) {
  // Asserts for the runtime or namelist options
  EKAT_REQUIRE_MSG(m_params.isParameter("wsubmin"),
                   "ERROR: wsubmin is missing from mam_aci parameter list.");
  EKAT_REQUIRE_MSG(
      m_params.isParameter("top_level_mam4xx"),
      "ERROR: top_level_mam4xx is missing from mam_aci parameter list.");
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
  const FieldLayout scalar3d_int = grid_->get_3d_scalar_layout(false);

  // layout for 2D (1d horiz X 1d vertical) variable
  FieldLayout scalar2d_layout_col{{COL}, {ncol_}};

  using namespace ekat::units;
  constexpr auto q_unit = kg / kg;  // units of mass mixing ratios of tracers
  constexpr auto n_unit = 1 / kg;   // units of number mixing ratios of tracers

  constexpr auto nondim = ekat::units::Units::nondimensional();

  // atmospheric quantities
  // specific humidity [kg/kg]
  add_tracer<Required>("qv", grid_, q_unit);

  // cloud liquid mass mixing ratio [kg/kg]
  add_tracer<Required>("qc", grid_, q_unit);

  // cloud ice mass mixing ratio [kg/kg]
  add_tracer<Required>("qi", grid_, q_unit);

  // cloud liquid number mixing ratio [1/kg]
  add_tracer<Required>("nc", grid_, n_unit);

  // cloud ice number mixing ratio [1/kg]
  add_tracer<Required>("ni", grid_, n_unit);

  // Temperature[K] at midpoints
  add_field<Required>("T_mid", scalar3d_mid, K, grid_name);

  // Vertical pressure velocity [Pa/s] at midpoints
  add_field<Required>("omega", scalar3d_mid, Pa / s, grid_name);

  // Total pressure [Pa] at midpoints
  add_field<Required>("p_mid", scalar3d_mid, Pa, grid_name);

  // Total pressure [Pa] at interfaces
  add_field<Required>("p_int", scalar3d_int, Pa, grid_name);

  // Layer thickness(pdel) [Pa] at midpoints
  add_field<Required>("pseudo_density", scalar3d_mid, Pa, grid_name);

  // planetary boundary layer height
  add_field<Required>("pbl_height", scalar2d_layout_col, m, grid_name);

  // cloud fraction [nondimensional] computed by eamxx_cld_fraction_process
  add_field<Required>("cldfrac_tot", scalar3d_mid, nondim, grid_name);

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
  for(int mode = 0; mode < mam_coupling::num_aero_modes(); ++mode) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const char *int_nmr_field_name =
        mam_coupling::int_aero_nmr_field_name(mode);
    add_tracer<Updated>(int_nmr_field_name, grid_, n_unit);

    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    // NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are
    // NOT advected
    const char *cld_nmr_field_name =
        mam_coupling::cld_aero_nmr_field_name(mode);
    add_field<Updated>(cld_nmr_field_name, scalar3d_mid, n_unit, grid_name);

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const char *int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(mode, a);
      if(strlen(int_mmr_field_name) > 0) {
        add_tracer<Updated>(int_mmr_field_name, grid_, q_unit);
      }
      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      // NOTE: DO NOT add cld borne aerosols to the "tracer" group as these are
      // NOT advected
      const char *cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(mode, a);
      if(strlen(cld_mmr_field_name) > 0) {
        add_field<Updated>(cld_mmr_field_name, scalar3d_mid, q_unit, grid_name);
      }
    }  // end for loop num species
  }    // end for loop for num modes

  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    add_tracer<Updated>(gas_mmr_field_name, grid_, q_unit);
  }  // end for loop num gases

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

  wsubmin_ = m_params.get<double>("wsubmin");
  top_lev_ = m_params.get<int>("top_level_mam4xx");

  // ------------------------------------------------------------------------
  // Input fields read in from IC file, namelist or other processes
  // ------------------------------------------------------------------------
  w_sec_mid_    = get_field_in("w_variance").get_view<const Real **>();
  dgnum_        = get_field_in("dgnum").get_view<const Real ***>();
  liqcldf_      = get_field_in("cldfrac_liq").get_view<const Real **>();
  liqcldf_prev_ = get_field_in("cldfrac_liq_prev").get_view<const Real **>();
  kvh_mid_      = get_field_in("eddy_diff_heat").get_view<const Real **>();

  // store fields only to be converted to dry mmrs in wet_atm_
  wet_atm_.qv = get_field_in("qv").get_view<const Real **>();
  wet_atm_.qc = get_field_in("qc").get_view<const Real **>();
  wet_atm_.nc = get_field_in("nc").get_view<const Real **>();
  wet_atm_.qi = get_field_in("qi").get_view<const Real **>();
  wet_atm_.ni = get_field_in("ni").get_view<const Real **>();

  // store rest fo the atm fields in dry_atm_in
  dry_atm_.z_surf = 0;
  dry_atm_.T_mid  = get_field_in("T_mid").get_view<const Real **>();
  dry_atm_.p_mid  = get_field_in("p_mid").get_view<const Real **>();
  dry_atm_.p_int  = get_field_in("p_int").get_view<const Real **>();
  dry_atm_.p_del  = get_field_in("pseudo_density").get_view<const Real **>();
  dry_atm_.omega  = get_field_in("omega").get_view<const Real **>();

  // store fields converted to dry mmr from wet mmr in dry_atm_
  dry_atm_.qv = buffer_.qv_dry;
  dry_atm_.qc = buffer_.qc_dry;
  dry_atm_.nc = buffer_.nc_dry;
  dry_atm_.qi = buffer_.qi_dry;
  dry_atm_.ni = buffer_.ni_dry;

  // pbl_height
  dry_atm_.pblh = get_field_in("pbl_height").get_view<const Real *>();

  // geometric thickness of layers (m)
  dry_atm_.dz = buffer_.dz;

  // geopotential height above surface at interface levels (m)
  dry_atm_.z_iface = buffer_.z_iface;

  // geopotential height above surface at mid levels (m)
  dry_atm_.z_mid = buffer_.z_mid;

  // total cloud fraction
  dry_atm_.cldfrac = get_field_in("cldfrac_tot").get_view<const Real **>();

  // computed updraft velocity
  dry_atm_.w_updraft = buffer_.w_updraft;

  // ------------------------------------------------------------------------
  // Output fields to be used by other processes
  // ------------------------------------------------------------------------
  // ice nucleation output
  naai_ = get_field_out("ni_activated").get_view<Real **>();

  // droplet activation output
  tendnd_ = get_field_out("nc_nuceat_tend").get_view<Real **>();

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    // interstitial aerosol tracers of interest: number (n) mixing ratios
    const char *int_nmr_field_name = mam_coupling::int_aero_nmr_field_name(m);
    wet_aero_.int_aero_nmr[m] =
        get_field_out(int_nmr_field_name).get_view<Real **>();
    dry_aero_.int_aero_nmr[m] = buffer_.dry_int_aero_nmr[m];

    // cloudborne aerosol tracers of interest: number (n) mixing ratios
    const char *cld_nmr_field_name = mam_coupling::cld_aero_nmr_field_name(m);
    wet_aero_.cld_aero_nmr[m] =
        get_field_out(cld_nmr_field_name).get_view<Real **>();
    dry_aero_.cld_aero_nmr[m] = buffer_.dry_cld_aero_nmr[m];

    for(int a = 0; a < mam_coupling::num_aero_species(); ++a) {
      // (interstitial) aerosol tracers of interest: mass (q) mixing ratios
      const char *int_mmr_field_name =
          mam_coupling::int_aero_mmr_field_name(m, a);

      if(strlen(int_mmr_field_name) > 0) {
        wet_aero_.int_aero_mmr[m][a] =
            get_field_out(int_mmr_field_name).get_view<Real **>();
        dry_aero_.int_aero_mmr[m][a] = buffer_.dry_int_aero_mmr[m][a];
      }

      // (cloudborne) aerosol tracers of interest: mass (q) mixing ratios
      const char *cld_mmr_field_name =
          mam_coupling::cld_aero_mmr_field_name(m, a);
      if(strlen(cld_mmr_field_name) > 0) {
        wet_aero_.cld_aero_mmr[m][a] =
            get_field_out(cld_mmr_field_name).get_view<Real **>();
        dry_aero_.cld_aero_mmr[m][a] = buffer_.dry_cld_aero_mmr[m][a];
      }
    }
  }
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    const char *gas_mmr_field_name = mam_coupling::gas_mmr_field_name(g);
    wet_aero_.gas_mmr[g] =
        get_field_out(gas_mmr_field_name).get_view<Real **>();
    dry_aero_.gas_mmr[g] = buffer_.dry_gas_mmr[g];
  }

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

  Kokkos::resize(rho_, ncol_, nlev_);
  Kokkos::resize(w0_, ncol_, nlev_);
  Kokkos::resize(tke_, ncol_, nlev_ + 1);
  Kokkos::resize(wsub_, ncol_, nlev_);
  Kokkos::resize(wsubice_, ncol_, nlev_);
  Kokkos::resize(wsig_, ncol_, nlev_);
  Kokkos::resize(w2_, ncol_, nlev_);
  Kokkos::resize(cloud_frac_, ncol_, nlev_);
  Kokkos::resize(cloud_frac_prev_, ncol_, nlev_);
  Kokkos::resize(aitken_dry_dia_, ncol_, nlev_);
  Kokkos::resize(rpdel_, ncol_, nlev_);

  //---------------------------------------------------------------------------------
  // Diagnotics variables from the ice nucleation scheme
  //---------------------------------------------------------------------------------

  // number conc of ice nuclei due to heterogeneous freezing [1/m3]
  Kokkos::resize(nihf_, ncol_, nlev_);

  // number conc of ice nuclei due to immersionfreezing (hetero nuc) [1/m3]
  Kokkos::resize(niim_, ncol_, nlev_);

  // number conc of ice nuclei due to deposition nucleation (hetero nuc)[1/m3]
  Kokkos::resize(nidep_, ncol_, nlev_);

  // number conc of ice nuclei due to meyers deposition [1/m3]
  Kokkos::resize(nimey_, ncol_, nlev_);

  // number of activated aerosol for ice nucleation(homogeneous frz only)[#/kg]
  Kokkos::resize(naai_hom_, ncol_, nlev_);

  //---------------------------------------------------------------------------------
  // Diagnotics variables from the droplet activation scheme
  //---------------------------------------------------------------------------------

  // activation fraction for aerosol number [fraction]
  const int num_aero_modes = mam_coupling::num_aero_modes();
  Kokkos::resize(factnum_, ncol_, num_aero_modes, nlev_);

  // cloud droplet number mixing ratio [#/kg]
  Kokkos::resize(qcld_, ncol_, nlev_);

  // number conc of aerosols activated at supersat [#/m^3]
  // NOTE:  activation fraction fluxes are defined as
  // fluxn = [flux of activated aero. number into cloud[#/m^2/s]]
  //        / [aero. number conc. in updraft, just below cloudbase [#/m^3]]
  Kokkos::resize(ccn_, ncol_, nlev_, mam4::ndrop::psat);

  // column-integrated droplet number [#/m2]
  Kokkos::resize(ndropcol_, ncol_, nlev_);

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
  Kokkos::resize(wtke_, ncol_, nlev_);

  for(int i = 0; i < dropmix_scratch_; ++i) {
    Kokkos::resize(dropmixnuc_scratch_mem_[i], ncol_, nlev_);
  }
  for(int i = 0; i < mam4::ndrop::ncnst_tot; ++i) {
    // column tendency for diagnostic output
    Kokkos::resize(coltend_[i], ncol_, nlev_);
    // column tendency
    Kokkos::resize(coltend_cw_[i], ncol_, nlev_);
  }
  for(int i = 0; i < mam4::aero_model::pcnst; ++i) {
    Kokkos::resize(ptend_q_[i], ncol_, nlev_);
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
  }
  state_q_work_ =
      view_3d("state_q_work_", ncol_, nlev_, mam4::aero_model::pcnst);

  //---------------------------------------------------------------------------------
  // Diagnotics variables from the hetrozenous ice nucleation scheme
  //---------------------------------------------------------------------------------

  for(int i = 0; i < hetro_scratch_; ++i)
    Kokkos::resize(diagnostic_scratch_[i], ncol_, nlev_);

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

  postprocess_.initialize(ncol_, nlev_, wet_atm_, wet_aero_, dry_atm_,
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
      cloud_frac_prev_, dry_aero_, nlev_,
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
  Kokkos::parallel_for("postprocess", scan_policy, postprocess_);
  Kokkos::fence();  // wait before returning to calling function
}

}  // namespace scream
