#include "physics/mam/eamxx_mam_dry_deposition_process_interface.hpp"

// Drydep functions are stored in the following hpp file
#include <physics/mam/eamxx_mam_dry_deposition_functions.hpp>

// For reading fractional land use file
#include <physics/mam/readfiles/fractional_land_use.hpp>

namespace scream {

using FracLandUseFunc = frac_landuse::fracLandUseFunctions<Real, DefaultDevice>;

MAMDryDep::MAMDryDep(const ekat::Comm &comm, const ekat::ParameterList &params)
    : MAMGenericInterface(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
  check_fields_intervals_ =
      m_params.get<bool>("create_fields_interval_checks", false);
}

// ================================================================
//  SET_GRIDS
// ================================================================
void MAMDryDep::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;

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

  // layout for 2D (ncol, pcnst)
  constexpr int pcnst = mam4::aero_model::pcnst;
  const FieldLayout vector2d_pcnst =
      grid_->get_2d_vector_layout(pcnst, "num_phys_constants");
  const FieldLayout vector2d_class =
      grid_->get_2d_vector_layout(n_land_type, "class");

  // Layout for 4D (2d horiz X 1d vertical x number of modes) variables
  // at mid points
  const int num_aero_modes       = mam_coupling::num_aero_modes();
  const FieldLayout vector3d_mid = grid_->get_3d_vector_layout(
      true, num_aero_modes, mam_coupling::num_modes_tag_name());

  using namespace ekat::units;
  auto nondim = ekat::units::Units::nondimensional();

  auto m3 = m * m * m;  // meter cubed

  // --------------------------------------------------------------------------
  // These variables are "Required" or pure inputs for the process
  // --------------------------------------------------------------------------

  add_tracers_wet_atm();
  add_fields_dry_atm();

  static constexpr auto m2 = m * m;
  static constexpr auto s2 = s * s;

  // Surface geopotential [m2/s2] (Require only for building DS)
  add_field<Required>("phis", scalar2d, m2 / s2, grid_name);

  //----------- Variables from microphysics scheme -------------

  // Total cloud fraction [fraction] (Require only for building DS)
  add_field<Required>("cldfrac_tot", scalar3d_mid, nondim, grid_name);

  //----------- Variables from coupler (land component)---------
  // Obukhov length [m]
  add_field<Required>("obklen", scalar2d, m, grid_name);

  // Surface friction velocty or ustar[m/s]
  add_field<Required>("ustar", scalar2d, m / s, grid_name);

  // Land fraction [fraction]
  add_field<Required>("landfrac", scalar2d, nondim, grid_name);

  // Friction velocity from land model [m/s]
  add_field<Required>("fv", scalar2d, m / s, grid_name);

  // Aerodynamical resistance from land model [s/m]
  add_field<Required>("ram1", scalar2d, s / m, grid_name);

  //----------- Variables from coupler (ice component)---------

  // Ice fraction [unitless]
  add_field<Required>("icefrac", scalar2d, nondim, grid_name);

  //----------- Variables from coupler (ocean component)---------
  // Ocean fraction [unitless]
  add_field<Required>("ocnfrac", scalar2d, nondim, grid_name);

  //----------- Variables from other mam4xx processes ------------
  // Geometric mean wet diameter for number distribution [m]
  add_field<Required>("dgnumwet", vector3d_mid, m, grid_name);

  // Wet density of interstitial aerosol [kg/m3]
  add_field<Required>("wetdens", vector3d_mid, kg / m3, grid_name);

  // ---------------------------------------------------------------------
  // These variables are "updated" or inputs/outputs for the process
  // ---------------------------------------------------------------------

  // (interstitial) aerosol tracers of interest: mass (q) and number (n) mixing
  // ratios
  // add tracers, e.g., num_a1, soa_a1
  add_tracers_interstitial_aerosol();
  // add tracer gases, e.g., O3
  add_tracers_gases();
  // add fields e.g., num_c1, soa_c1
  add_fields_cloudborne_aerosol();

  // -------------------------------------------------------------
  // These variables are "Computed" or outputs for the process
  // -------------------------------------------------------------
  // FIXME: These are diagnostics, remove them from FM after initial evaluation
  // surface deposition flux of cloud-borne  aerosols, [kg/m2/s] or [1/m2/s]
  add_field<Computed>("deposition_flux_of_cloud_borne_aerosols", vector2d_pcnst,
                      1 / m2 / s, grid_name);
  // surface deposition flux of interstitial aerosols, [kg/m2/s] or [1/m2/s]
  add_field<Computed>("deposition_flux_of_interstitial_aerosols",
                      vector2d_pcnst, 1 / m2 / s, grid_name);

  // Fractional land use [fraction]
  add_field<Computed>("fraction_landuse", vector2d_class, nondim, grid_name);
  // -------------------------------------------------------------
  // setup to enable reading fractional land use file
  // -------------------------------------------------------------

  const auto mapping_file = m_params.get<std::string>("drydep_remap_file", "");
  const std::string frac_landuse_data_file =
      m_params.get<std::string>("fractional_land_use_file");

  // Field to be read from file
  const std::string field_name = "fraction_landuse";

  // Dimensions of the filed
  const std::string dim_name1 = "ncol";
  const std::string dim_name2 = "class";

  // initialize the file read
  FracLandUseFunc::init_frac_landuse_file_read(
      ncol_, field_name, dim_name1, dim_name2, grid_, frac_landuse_data_file,
      mapping_file, horizInterp_, dataReader_);  // output

}  // set_grids

// ================================================================
//  REQUEST_BUFFER_SIZE_IN_BYTES
// ================================================================
// ON HOST, returns the number of bytes of device memory needed by
// the above. Buffer type given the number of columns and vertical
// levels
size_t MAMDryDep::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_,0, 0);
}  // requested_buffer_size_in_bytes

// ================================================================
//  INIT_BUFFERS
// ================================================================
// ON HOST, initializeѕ the Buffer type with sufficient memory to
// store intermediate (dry) quantities on the given number of
// columns with the given number of vertical levels. Returns the
// number of bytes allocated.
void MAMDryDep::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");

  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_, 0);
  EKAT_REQUIRE_MSG(used_mem == requested_buffer_size_in_bytes(),
                   "Error! Used memory != requested memory for MAMDryDep.");
}  // init_buffers

// ================================================================
//  INITIALIZE_IMPL
// ================================================================
void MAMDryDep::initialize_impl(const RunType run_type) {
  // ---------------------------------------------------------------
  // Input fields read in from IC file, namelist or other processes
  // ---------------------------------------------------------------
  // Check the interval values for the following fields used by this interface.
  // NOTE: We do not include aerosol and gas species, e.g., soa_a1, num_a1,
  // because we automatically added these fields.
  const std::map<std::string, std::pair<Real, Real>> ranges_dry_deposition = {
      // dry deposition
      {"dgnumwet", {-1e10, 1e10}},                                    // FIXME
      {"fv", {-1e10, 1e10}},                                          // FIXME
      {"icefrac", {-1e10, 1e10}},                                     // FIXME
      {"landfrac", {-1e10, 1e10}},                                    // FIXME
      {"obklen", {-1e10, 1e10}},                                      // FIXME
      {"ocnfrac", {-1e10, 1e10}},                                     // FIXME
      {"ram1", {-1e10, 1e10}},                                        // FIXME
      {"ustar", {-1e10, 1e10}},                                       // FIXME
      {"wetdens", {-1e10, 1e10}},                                     // FIXME
      {"deposition_flux_of_cloud_borne_aerosols", {-1e100, 1e100}},   // FIXME
      {"deposition_flux_of_interstitial_aerosols", {-1e100, 1e100}},  // FIXME
      {"fraction_landuse", {-1e100, 1e100}},                          // FIXME
  };
  set_ranges_process(ranges_dry_deposition);
  // Check pre/post condition interval values for all fields employed by this
  // interface
  add_interval_checks();

  populate_wet_atm(wet_atm_);
  populate_dry_atm(dry_atm_, buffer_);
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

  //-----------------------------------------------------------------
  // Allocate memory
  //-----------------------------------------------------------------
  const int pcnst = mam4::aero_model::pcnst;

  // Output of the the mixing ratio tendencies [kg/kg/s or 1/kg/s]
  ptend_q_ = view_3d("ptend_q_", ncol_, nlev_, pcnst);

  // Deposition velocity of turbulent dry deposition [m/s]
  vlc_trb_ = view_3d("vlc_trb_", mam4::AeroConfig::num_modes(),
                     aerosol_categories_, ncol_);
  // Deposition velocity of gravitational settling [m/s]
  vlc_grv_ = view_4d("vlc_grv_", mam4::AeroConfig::num_modes(),
                     aerosol_categories_, ncol_, nlev_);
  // Deposition velocity, [m/s]
  // Fraction landuse weighted sum of vlc_grv and vlc_trb
  vlc_dry_ = view_4d("vlc_dry_", mam4::AeroConfig::num_modes(),
                     aerosol_categories_, ncol_, nlev_);

  // Work array to hold the mixing ratios [kg/kg or 1/kg]
  // Packs AerosolState::int_aero_nmr and AerosolState::int_aero_nmr
  // into one array.
  qtracers_ = view_3d("qtracers_", ncol_, nlev_, pcnst);

  // Work array to hold the air density [kg/m3]
  rho_ = view_2d("rho", ncol_, nlev_);

  // Work array to hold cloud borne aerosols mixing ratios [kg/kg or 1/kg]
  // Filled with Prognostics::n_mode_c and Prognostics::q_aero_c
  qqcw_ = view_3d("qqcw_", pcnst, ncol_, nlev_);

  // Work array to hold tendency for 1 species [kg/kg/s] or [1/kg/s]
  dqdt_tmp_ = view_3d("dqdt_tmp_", pcnst, ncol_, nlev_);

  //-----------------------------------------------------------------
  // Read fractional land use data
  //-----------------------------------------------------------------
  frac_landuse_fm_ = get_field_out("fraction_landuse").get_view<Real **>();
  // This data is time-independent, we read all data here for the
  // entire simulation
  FracLandUseFunc::update_frac_land_use_data_from_file(
      dataReader_, *horizInterp_,
      frac_landuse_);  // output

  // Copy fractional landuse values to a FM array to be used by other processes
  Kokkos::deep_copy(frac_landuse_fm_, frac_landuse_);
}  // initialize_impl

// =========================================================================================
void MAMDryDep::run_impl(const double dt) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  pre_process(wet_aero_, dry_aero_, wet_atm_, dry_atm_);
  Kokkos::fence();

  // -------------------------------------------------------------
  // Inputs fields for the process
  // -------------------------------------------------------------

  // Geometric mean wet diameter for number distribution [m]
  auto dgncur_awet_ = get_field_in("dgnumwet").get_view<const Real ***>();
  // Wet density of interstitial aerosol [kg/m3]
  auto wet_dens_ = get_field_in("wetdens").get_view<const Real ***>();
  // Obukhov length [m]
  auto obukhov_length_ = get_field_in("obklen").get_view<const Real *>();
  // Land fraction [unitless]
  auto land_fraction_ = get_field_in("landfrac").get_view<const Real *>();
  // Ice fraction [unitless]
  auto ice_fraction_ = get_field_in("icefrac").get_view<const Real *>();
  // Ocean fraction [unitless]
  auto ocean_fraction_ = get_field_in("ocnfrac").get_view<const Real *>();
  // Friction velocity from land model [m/s]
  auto friction_velocity_ = get_field_in("fv").get_view<const Real *>();
  // Aerodynamical resistance from land model [s/m]
  auto aerodynamical_resistance_ =
      get_field_in("ram1").get_view<const Real *>();
  //  Sfc friction velocity or ustar [m/s]
  auto surface_friction_velocty_ =
      get_field_in("ustar").get_view<const Real *>();

  // -------------------------------------------------------------
  // Output fields for the process
  // -------------------------------------------------------------
  // Surface deposition flux of cloud-borne  aerosols, [kg/m2/s] or [1/m2/s]
  auto aerdepdrycw_ = get_field_out("deposition_flux_of_cloud_borne_aerosols")
                          .get_view<Real **>();
  // Surface deposition flux of interstitial aerosols, [kg/m2/s] or [1/m2/s]
  auto aerdepdryis_ = get_field_out("deposition_flux_of_interstitial_aerosols")
                          .get_view<Real **>();

  //--------------------------------------------------------------------
  // Call drydeposition and get tendencies
  //--------------------------------------------------------------------
  compute_tendencies(ncol_, nlev_, dt, obukhov_length_,
                     surface_friction_velocty_, land_fraction_, ice_fraction_,
                     ocean_fraction_, friction_velocity_,
                     aerodynamical_resistance_, frac_landuse_, dgncur_awet_,
                     wet_dens_, dry_atm_, dry_aero_,
                     // Inouts-outputs
                     qqcw_,
                     // Outputs
                     ptend_q_, aerdepdrycw_, aerdepdryis_,
                     // work arrays
                     rho_, vlc_dry_, vlc_trb_, vlc_grv_, dqdt_tmp_, qtracers_);
  Kokkos::fence();

  // Update the interstitial aerosols using ptend.
  update_interstitial_mmrs(ptend_q_, dt, ncol_, nlev_,  // inputs
                           dry_aero_);                  // output

  // Update the interstitial aerosols
  update_cloudborne_mmrs(qqcw_, dt, nlev_,  // inputs
                         dry_aero_);        // output

  // call post processing to convert dry mixing ratios to wet mixing ratios
  // and update the state
  post_process(wet_aero_, dry_aero_, dry_atm_);
  Kokkos::fence();  // wait before returning to calling function
}  // run_impl
}  // namespace scream
