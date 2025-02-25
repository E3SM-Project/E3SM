#include "physics/mam/eamxx_mam_wetscav_process_interface.hpp"

/*
-----------------------------------------------------------------
NOTES:
1. We should connect surface fluxes and add code to update the fluxes
2. Identify diagnostic variables and remove them from FM
3. Add assert statements to check output ranges
*/

namespace scream {

// =========================================================================================
MAMWetscav::MAMWetscav(const ekat::Comm &comm,
                       const ekat::ParameterList &params)
    : MAMGenericInterface(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
    check_fields_intervals_   = m_params.get<bool>("mam4_check_fields_intervals", false);

}

// ================================================================
//  SET_GRIDS
// ================================================================
void MAMWetscav::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

  grid_                = grids_manager->get_grid("Physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column
  const int nmodes    = mam4::AeroConfig::num_modes();  // Number of modes
  constexpr int pcnst = mam4::aero_model::pcnst;

  // layout for 3D (2d horiz X 1d vertical) variables at level
  // midpoints/interfaces
  FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);
  FieldLayout scalar3d_int = grid_->get_3d_scalar_layout(false);

  // layout for 2D (1d horiz X 1d vertical) variables
  FieldLayout scalar2d = grid_->get_2d_scalar_layout();

  // layout for 3D (ncol, nmodes, nlevs)
  FieldLayout scalar3d_mid_nmodes =
      grid_->get_3d_vector_layout(true, nmodes, mam_coupling::num_modes_tag_name());

  // layout for 2D (ncol, pcnst)
  FieldLayout scalar2d_pconst =
      grid_->get_2d_vector_layout(pcnst, "num_phys_constants");

  // --------------------------------------------------------------------------
  // These variables are "required" or pure inputs for the process
  // --------------------------------------------------------------------------

  // ----------- Atmospheric quantities -------------
  // Specific humidity [kg/kg]
  add_tracers_wet_and_dry_atm();

  static constexpr auto m2 = m * m;
  static constexpr auto s2 = s * s;

  // Surface geopotential [m2/s2]
  add_field<Required>("phis", scalar2d, m2 / s2, grid_name);

  //----------- Variables from microphysics scheme -------------

  // Evaporation from stratiform rain [kg/kg/s]
  add_field<Required>("nevapr", scalar3d_mid, kg / kg / s, grid_name);

  // Stratiform rain production rate [kg/kg/s]
  add_field<Required>("precip_total_tend", scalar3d_mid, kg / kg / s,
                      grid_name);

  // For variables that are non dimensional (e.g., fractions etc.)
  static constexpr auto nondim = Units::nondimensional();

  //----------- Variables from macrophysics scheme -------------

  // Total cloud fraction [fraction]
  add_field<Required>("cldfrac_liq", scalar3d_mid, nondim, grid_name);

  // ---------------------------------------------------------------------
  // These variables are "updated" or inputs/outputs for the process
  // ---------------------------------------------------------------------
  // FIXME: we have not added code to update the surface fluxes.
  // -- surface fluxes (input/outpts) for the coupler's cam_out data struture
  // for the land model

  // Wet deposition of hydrophilic black carbon [kg/m2/s]
  add_field<Updated>("wetdep_hydrophilic_bc", scalar3d_mid, kg / m2 / s,
                     grid_name);

  // Dry deposition of hydrophilic black carbon [kg/m2/s]
  add_field<Updated>("drydep_hydrophilic_bc", scalar3d_mid, kg / m2 / s,
                     grid_name);

  // Wet deposition of hydrophilic organic carbon [kg/m2/s]
  add_field<Updated>("wetdep_hydrophilic_oc", scalar3d_mid, kg / m2 / s,
                     grid_name);

  // Dry deposition of hydrophilic organic carbon [kg/m2/s]
  add_field<Updated>("drydep_hydrophilic_oc", scalar3d_mid, kg / m2 / s,
                     grid_name);

  // Wet deposition of dust (bin1) [kg/m2/s]
  add_field<Updated>("wetdep_dust_bin1", scalar3d_mid, kg / m2 / s, grid_name);

  // Wet deposition of dust (bin2) [kg/m2/s]
  add_field<Updated>("wetdep_dust_bin2", scalar3d_mid, kg / m2 / s, grid_name);

  // Wet deposition of dust (bin3) [kg/m2/s]
  add_field<Updated>("wetdep_dust_bin3", scalar3d_mid, kg / m2 / s, grid_name);

  // Wet deposition of dust (bin4) [kg/m2/s]
  add_field<Updated>("wetdep_dust_bin4", scalar3d_mid, kg / m2 / s, grid_name);

  // Interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios

  // NOTE: Interstitial aerosols are updated in the interface using the
  // "tendencies" from the wetscavenging process
  add_tracers_aerosol_and_gases();

  // -------------------------------------------------------------
  // These variables are "Computed" or outputs for the process
  // -------------------------------------------------------------
  static constexpr auto m3 = m * m * m;

  // Aerosol dry particle diameter [m]
  add_field<Computed>("dgnum", scalar3d_mid_nmodes, m, grid_name);

  // Wet aerosol density [kg/m3]
  add_field<Computed>("wetdens", scalar3d_mid_nmodes, kg / m3, grid_name);

  // Aerosol water [kg/kg]
  add_field<Computed>("qaerwat", scalar3d_mid_nmodes, kg / kg, grid_name);

  // Wet aerosol diameter [m]
  add_field<Computed>("dgnumwet", scalar3d_mid_nmodes, m, grid_name);

  // Fraction of transported species that are insoluble [fraction]
  add_field<Computed>("fracis", scalar3d_mid, nondim, grid_name);

  // Aerosol wet deposition (interstitial) [kg/m2/s]
  add_field<Computed>("aerdepwetis", scalar2d_pconst, kg / m2 / s, grid_name);

  // Aerosol wet deposition (cloud water) [kg/m2/s]
  add_field<Computed>("aerdepwetcw", scalar2d_pconst, kg / m2 / s, grid_name);
}

// ================================================================
//  INIT_BUFFERS
// ================================================================
// ON HOST, initializeѕ the Buffer type with sufficient memory to store
// intermediate (dry) quantities on the given number of columns with the given
// number of vertical levels. Returns the number of bytes allocated.
void MAMWetscav::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");

  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for MAMWetscav.");
}

// ================================================================
//  INITIALIZE_IMPL
// ================================================================
void MAMWetscav::initialize_impl(const RunType run_type) {
  // ---------------------------------------------------------------
  // Input fields read in from IC file, namelist or other processes
  // ---------------------------------------------------------------
  // print_fields_names();
  add_interval_checks();

  populate_wet_and_dry_atm();
  dry_atm_.phis  = get_field_in("phis").get_view<const Real *>();

  // store fields converted to dry mmr from wet mmr in dry_atm_
  dry_atm_.qv        = buffer_.qv_dry;
  dry_atm_.qc        = buffer_.qc_dry;
  dry_atm_.nc        = buffer_.nc_dry;
  dry_atm_.qi        = buffer_.qi_dry;
  dry_atm_.ni        = buffer_.ni_dry;
  dry_atm_.z_mid     = buffer_.z_mid;
  dry_atm_.dz        = buffer_.dz;
  dry_atm_.z_iface   = buffer_.z_iface;
  dry_atm_.w_updraft = buffer_.w_updraft;

  // interstitial and cloudborne aerosol tracers of interest: mass (q) and
  // number (n) mixing ratios
  populate_wet_and_dry_aero();

  //---------------------------------------------------------------------------------
  // Allocate memory
  //---------------------------------------------------------------------------------
  // Alllocate aerosol-related gas tendencies
  for(int g = 0; g < mam_coupling::num_aero_gases(); ++g) {
    dry_aero_tends_.gas_mmr[g] = view_2d("gas_mmr", ncol_, nlev_);
  }

  // Allocate aerosol state tendencies (interstitial aerosols only)
  for(int imode = 0; imode < mam_coupling::num_aero_modes(); ++imode) {
    dry_aero_tends_.int_aero_nmr[imode] = view_2d("int_aero_nmr", ncol_, nlev_);

    for(int ispec = 0; ispec < mam_coupling::num_aero_species(); ++ispec) {
      dry_aero_tends_.int_aero_mmr[imode][ispec] =
          view_2d("int_aero_mmr", ncol_, nlev_);
    }
  }

  // Allocate work array
  const int work_len = mam4::wetdep::get_aero_model_wetdep_work_len();
  work_              = view_2d("work", ncol_, work_len);

  // TODO: Following variables are from convective parameterization (not
  // implemented yet in EAMxx), so should be zero for now

  sh_frac_ = view_2d("sh_frac", ncol_, nlev_);
  Kokkos::deep_copy(sh_frac_, 0);

  // Deep convective cloud fraction [fraction]
  dp_frac_ = view_2d("dp_frac", ncol_, nlev_);
  Kokkos::deep_copy(dp_frac_, 0);

  // Evaporation rate of shallow convective precipitation >=0. [kg/kg/s]
  evapcsh_ = view_2d("evapcsh", ncol_, nlev_);
  Kokkos::deep_copy(evapcsh_, 0);

  // Evaporation rate of deep convective precipitation >=0. [kg/kg/s]
  evapcdp_ = view_2d("evapcdp", ncol_, nlev_);
  Kokkos::deep_copy(evapcdp_, 0);

  // Rain production, shallow convection [kg/kg/s]
  rprdsh_ = view_2d("rprdsh", ncol_, nlev_);
  Kokkos::deep_copy(rprdsh_, 0);

  // Rain production, deep convection [kg/kg/s]
  rprddp_ = view_2d("rprddp", ncol_, nlev_);
  Kokkos::deep_copy(rprddp_, 0);

  // In cloud water mixing ratio, deep convection
  icwmrdp_ = view_2d("icwmrdp", ncol_, nlev_);
  Kokkos::deep_copy(icwmrdp_, 0);

  // In cloud water mixing ratio, shallow convection
  icwmrsh_ = view_2d("icwmrsh", ncol_, nlev_);
  Kokkos::deep_copy(icwmrsh_, 0);

  // Detraining cld H20 from deep convection [kg/kg/s]
  dlf_ = view_2d("dlf", ncol_, nlev_);
  Kokkos::deep_copy(dlf_, 0);
  }

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMWetscav::run_impl(const double dt) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of all variables
  // needed by this process or setting up MAM4xx classes and their objects
  pre_process();
  Kokkos::fence();

  const mam_coupling::DryAtmosphere &dry_atm = dry_atm_;
  const auto &dry_aero                       = dry_aero_;
  const auto &work                           = work_;
  const auto &dry_aero_tends                 = dry_aero_tends_;

  // ---------------------------------------------------------------
  // These variables are "required" or pure inputs for the process
  // ---------------------------------------------------------------

  //----------- Variables from convective scheme -------------

  // TODO: Following variables are from convective parameterization (not
  // implemented yet in EAMxx), so should be zero for now

  auto sh_frac = sh_frac_;

  // Deep convective cloud fraction [fraction]
  auto dp_frac = dp_frac_;

  // Evaporation rate of shallow convective precipitation >=0. [kg/kg/s]
  auto evapcsh = evapcsh_;

  // Evaporation rate of deep convective precipitation >=0. [kg/kg/s]
  auto evapcdp = evapcdp_;

  // Rain production, shallow convection [kg/kg/s]
  auto rprdsh = rprdsh_;

  // Rain production, deep convection [kg/kg/s]
  auto rprddp = rprddp_;

  // In cloud water mixing ratio, deep convection
  auto icwmrdp = icwmrdp_;

  // In cloud water mixing ratio, shallow convection
  auto icwmrsh = icwmrsh_;

  // Detraining cld H20 from deep convection [kg/kg/s]
  auto dlf = dlf_;

  //----------- Variables from macrophysics scheme -------------
  // Total cloud fraction
  auto cldt = get_field_in("cldfrac_liq").get_view<const Real **>();

  //----------- Variables from microphysics scheme -------------

  // Evaporation from stratiform rain [kg/kg/s]
  auto nevapr = get_field_in("nevapr").get_view<const Real **>();

  // Stratiform rain production rate [kg/kg/s]
  auto prain = get_field_in("precip_total_tend").get_view<const Real **>();
  // ------------------------------------------------------------------
  // These variables are "Computed" or pure outputs for the process
  // ------------------------------------------------------------------

  const auto aerdepwetis = get_field_out("aerdepwetis").get_view<Real **>();
  const auto aerdepwetcw = get_field_out("aerdepwetcw").get_view<Real **>();

  const auto wet_geometric_mean_diameter_i =
      get_field_out("dgnumwet").get_view<Real ***>();
  const auto dry_geometric_mean_diameter_i =
      get_field_out("dgnum").get_view<Real ***>();
  const auto qaerwat = get_field_out("qaerwat").get_view<Real ***>();
  const auto wetdens = get_field_out("wetdens").get_view<Real ***>();

  const auto policy =
      ekat::ExeSpaceUtils<KT::ExeSpace>::get_default_team_policy(ncol_, nlev_);

  // Making a local copy of 'nlev_' because we cannot use a member of a class
  // inside a parallel_for.
  const int nlev = nlev_;

  // Zero out tendencies otherwise, they are initialized to junk values
  for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
    Kokkos::deep_copy(dry_aero_tends_.int_aero_nmr[m], 0);
    for(int a = 0; a < mam4::num_species_mode(m); ++a) {
      Kokkos::deep_copy(dry_aero_tends_.int_aero_mmr[m][a], 0);
    }
  }

  // Loop over atmosphere columns
  Kokkos::parallel_for(
      policy, KOKKOS_LAMBDA(const ThreadTeam &team) {
        const int icol = team.league_rank();  // column index

        const auto atm = mam_coupling::atmosphere_for_column(dry_atm, icol);
        // set surface state data
        // fetch column-specific subviews into aerosol prognostics
        mam4::Prognostics progs =
            mam_coupling::aerosols_for_column(dry_aero, icol);
        // fetch column-specific subviews into aerosol tendencies
        // Note: we are only updating interstitial aerosols.
        mam4::Tendencies tends =
            mam_coupling::interstitial_aerosols_tendencies_for_column(
                dry_aero_tends, icol);
        /// shallow_convective_precipitation_production
        const auto rprdsh_icol = ekat::subview(rprdsh, icol);
        // deep_convective_precipitation_production
        const auto rprddp_icol = ekat::subview(rprddp, icol);
        // deep_convective_precipitation_evaporation
        const auto evapcdp_icol = ekat::subview(evapcdp, icol);
        // shallow_convective_precipitation_evaporation =
        const auto evapcsh_icol = ekat::subview(evapcsh, icol);
        // deep_convective_cloud_fraction
        const auto dp_frac_icol = ekat::subview(dp_frac, icol);
        // shallow_convective_cloud_fraction    =
        const auto sh_frac_icol = ekat::subview(sh_frac, icol);

        const auto icwmrdp_col  = ekat::subview(icwmrdp, icol);
        const auto icwmrsh_icol = ekat::subview(icwmrsh, icol);
        const auto nevapr_icol  = ekat::subview(nevapr, icol);
        const auto cldt_icol    = ekat::subview(cldt, icol);

        const auto dlf_icol   = ekat::subview(dlf, icol);
        auto aerdepwetis_icol = ekat::subview(aerdepwetis, icol);
        auto aerdepwetcw_icol = ekat::subview(aerdepwetcw, icol);
        auto work_icol        = ekat::subview(work, icol);
        auto wet_diameter_icol =
            ekat::subview(wet_geometric_mean_diameter_i, icol);
        auto dry_diameter_icol =
            ekat::subview(dry_geometric_mean_diameter_i, icol);
        auto qaerwat_icol     = ekat::subview(qaerwat, icol);
        auto wetdens_icol     = ekat::subview(wetdens, icol);
        const auto prain_icol = ekat::subview(prain, icol);

        Real scavimptblnum[mam4::aero_model::nimptblgrow_total]
                          [mam4::AeroConfig::num_modes()];
        Real scavimptblvol[mam4::aero_model::nimptblgrow_total]
                          [mam4::AeroConfig::num_modes()];

        mam4::wetdep::init_scavimptbl(scavimptblvol, scavimptblnum);

        mam4::wetdep::aero_model_wetdep(
            team, atm, progs, tends, dt,
            // inputs
            cldt_icol, rprdsh_icol, rprddp_icol, evapcdp_icol, evapcsh_icol,
            dp_frac_icol, sh_frac_icol, icwmrdp_col, icwmrsh_icol, nevapr_icol,
            dlf_icol, prain_icol, scavimptblnum, scavimptblvol,
            // outputs
            wet_diameter_icol, dry_diameter_icol, qaerwat_icol, wetdens_icol,
            aerdepwetis_icol, aerdepwetcw_icol, work_icol);
        team.team_barrier();
        // update interstitial aerosol state
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev), [&](int kk) {
          for(int m = 0; m < mam_coupling::num_aero_modes(); ++m) {
            const auto n_mode_i       = progs.n_mode_i[m];
            const auto tends_n_mode_i = tends.n_mode_i[m];
            n_mode_i(kk) += tends_n_mode_i(kk) * dt;
            for(int a = 0; a < mam4::num_species_mode(m); ++a) {
              const auto q_aero_i       = progs.q_aero_i[m][a];
              const auto tends_q_aero_i = tends.q_aero_i[m][a];
              q_aero_i(kk) += tends_q_aero_i(kk) * dt;
            }
          }
        });  // parallel_for for update interstitial aerosol state
      });    // icol parallel_for loop

  // call post processing to convert dry mixing ratios to wet mixing ratios
  // and update the state
  post_process();
  Kokkos::fence();  // wait before returning to calling function
}

// =========================================================================================
}  // namespace scream
