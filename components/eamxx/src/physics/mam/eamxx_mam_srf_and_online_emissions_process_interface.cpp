#include <physics/mam/eamxx_mam_srf_and_online_emissions_process_interface.hpp>

// For surface and online emission functions
#include <physics/mam/eamxx_mam_srf_and_online_emissions_functions.hpp>

// For reading soil erodibility file
#include <physics/mam/readfiles/soil_erodibility.hpp>

namespace scream {

// For reading soil erodibility file
using soilErodibilityFunc =
    soil_erodibility::soilErodibilityFunctions<Real, DefaultDevice>;

// ================================================================
//  Constructor
// ================================================================
MAMSrfOnlineEmiss::MAMSrfOnlineEmiss(const ekat::Comm &comm,
                                     const ekat::ParameterList &params)
    : MAMGenericInterface(comm, params) {
  // FIXME: Do we want to read dust emiss factor from the namelist??
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants.
   */
  check_fields_intervals_ =
      m_params.get<bool>("create_fields_interval_checks", false);
}

// ================================================================
//  SET_GRIDS
// ================================================================
void MAMSrfOnlineEmiss::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  grid_                 = grids_manager->get_grid("Physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column

  using namespace ekat::units;
  constexpr auto m2     = pow(m, 2);
  constexpr auto s2     = pow(s, 2);
  constexpr auto nondim = ekat::units::Units::nondimensional();

  const FieldLayout scalar2d   = grid_->get_2d_scalar_layout();
  const FieldLayout scalar3d_m = grid_->get_3d_scalar_layout(true);   // mid
  const FieldLayout scalar3d_i = grid_->get_3d_scalar_layout(false);  // int

  // For U and V components of wind
  const FieldLayout vector3d = grid_->get_3d_vector_layout(true, 2);

  // For components of dust flux
  const FieldLayout vector4d = grid_->get_2d_vector_layout(4);

  // --------------------------------------------------------------------------
  // These variables are "Required" or pure inputs for the process
  // --------------------------------------------------------------------------

  // ----------- Atmospheric quantities -------------

  // -- Variables required for building DS to compute z_mid --
  // Specific humidity [kg/kg]
  add_tracers_wet_atm();
  add_fields_dry_atm();

  //----------- Variables from microphysics scheme -------------

  // Surface geopotential [m2/s2]
  add_field<Required>("phis", scalar2d, m2 / s2, grid_name);

  // -- Variables required for online dust and sea salt emissions --

  // U and V components of the wind[m/s]
  add_field<Required>("horiz_winds", vector3d, m / s, grid_name);

  //----------- Variables from coupler (ocean component)---------
  // Ocean fraction [unitless]
  add_field<Required>("ocnfrac", scalar2d, nondim, grid_name);

  // Sea surface temperature [K]
  add_field<Required>("sst", scalar2d, K, grid_name);

  // dust fluxes [kg/m^2/s]: Four flux values for each column
  add_field<Required>("dstflx", vector4d, kg / m2 / s, grid_name);

  // -------------------------------------------------------------
  // These variables are "Updated" or input-outputs for the process
  // -------------------------------------------------------------

  constexpr int pcnst = mam4::aero_model::pcnst;
  const FieldLayout vector2d_pcnst =
      grid_->get_2d_vector_layout(pcnst, "num_phys_constituents");

  // Constituent fluxes of species in [kg/m2/s]
  // FIXME: confirm if it is Updated or Computed
  add_field<Updated>("constituent_fluxes", vector2d_pcnst, kg / m2 / s,
                     grid_name);

  // Surface emissions remapping file
  auto srf_map_file = m_params.get<std::string>("srf_remap_file", "");

  // FIXME: We can extract the following info about each species
  // in a separate hpp file
  //--------------------------------------------------------------------
  // Init dms srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ dms;
  // File name, name and sectors
  dms.data_file    = m_params.get<std::string>("srf_emis_specifier_for_DMS");
  dms.species_name = "dms";
  dms.sectors      = {"DMS"};
  srf_emiss_species_.push_back(dms);  // add to the vector

  //--------------------------------------------------------------------
  // Init so2 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ so2;
  // File name, name and sectors
  so2.data_file    = m_params.get<std::string>("srf_emis_specifier_for_SO2");
  so2.species_name = "so2";
  so2.sectors      = {"AGR", "RCO", "SHP", "SLV", "TRA", "WST"};
  srf_emiss_species_.push_back(so2);  // add to the vector

  //--------------------------------------------------------------------
  // Init bc_a4 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ bc_a4;
  // File name, name and sectors
  bc_a4.data_file = m_params.get<std::string>("srf_emis_specifier_for_bc_a4");
  bc_a4.species_name = "bc_a4";
  bc_a4.sectors      = {"AGR", "ENE", "IND", "RCO", "SHP", "SLV", "TRA", "WST"};
  srf_emiss_species_.push_back(bc_a4);  // add to the vector

  //--------------------------------------------------------------------
  // Init num_a1 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ num_a1;
  // File name, name and sectors
  num_a1.data_file = m_params.get<std::string>("srf_emis_specifier_for_num_a1");
  num_a1.species_name = "num_a1";
  num_a1.sectors      = {"num_a1_SO4_AGR", "num_a1_SO4_SHP", "num_a1_SO4_SLV",
                         "num_a1_SO4_WST"};
  srf_emiss_species_.push_back(num_a1);  // add to the vector

  //--------------------------------------------------------------------
  // Init num_a2 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ num_a2;
  // File name, name and sectors
  num_a2.data_file = m_params.get<std::string>("srf_emis_specifier_for_num_a2");
  num_a2.species_name = "num_a2";
  num_a2.sectors      = {"num_a2_SO4_RCO", "num_a2_SO4_TRA"};
  srf_emiss_species_.push_back(num_a2);  // add to the vector

  //--------------------------------------------------------------------
  // Init num_a4 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ num_a4;
  // File name, name and sectors
  num_a4.data_file = m_params.get<std::string>("srf_emis_specifier_for_num_a4");
  num_a4.species_name = "num_a4";
  num_a4.sectors      = {
           "num_a1_BC_AGR",  "num_a1_BC_ENE",  "num_a1_BC_IND",  "num_a1_BC_RCO",
           "num_a1_BC_SHP",  "num_a1_BC_SLV",  "num_a1_BC_TRA",  "num_a1_BC_WST",
           "num_a1_POM_AGR", "num_a1_POM_ENE", "num_a1_POM_IND", "num_a1_POM_RCO",
           "num_a1_POM_SHP", "num_a1_POM_SLV", "num_a1_POM_TRA", "num_a1_POM_WST"};
  srf_emiss_species_.push_back(num_a4);  // add to the vector

  //--------------------------------------------------------------------
  // Init pom_a4 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ pom_a4;
  // File name, name and sectors
  pom_a4.data_file = m_params.get<std::string>("srf_emis_specifier_for_pom_a4");
  pom_a4.species_name = "pom_a4";
  pom_a4.sectors = {"AGR", "ENE", "IND", "RCO", "SHP", "SLV", "TRA", "WST"};
  srf_emiss_species_.push_back(pom_a4);  // add to the vector

  //--------------------------------------------------------------------
  // Init so4_a1 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ so4_a1;
  // File name, name and sectors
  so4_a1.data_file = m_params.get<std::string>("srf_emis_specifier_for_so4_a1");
  so4_a1.species_name = "so4_a1";
  so4_a1.sectors      = {"AGR", "SHP", "SLV", "WST"};
  srf_emiss_species_.push_back(so4_a1);

  //--------------------------------------------------------------------
  // Init so4_a2 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ so4_a2;
  // File name, name and sectors
  so4_a2.data_file = m_params.get<std::string>("srf_emis_specifier_for_so4_a2");
  so4_a2.species_name = "so4_a2";
  so4_a2.sectors      = {"RCO", "TRA"};
  srf_emiss_species_.push_back(so4_a2);

  //--------------------------------------------------------------------
  // Init data structures to read and interpolate
  //--------------------------------------------------------------------
  for(srf_emiss_ &ispec_srf : srf_emiss_species_) {
    srfEmissFunc::init_srf_emiss_objects(
        ncol_, grid_, ispec_srf.data_file, ispec_srf.sectors, srf_map_file,
        // output
        ispec_srf.horizInterp_, ispec_srf.data_start_, ispec_srf.data_end_,
        ispec_srf.data_out_, ispec_srf.dataReader_);
  }  // srf emissions file read init

  // -------------------------------------------------------------
  // Setup to enable reading soil erodibility file
  // -------------------------------------------------------------

  const std::string soil_erodibility_data_file =
      m_params.get<std::string>("soil_erodibility_file");

  // Field to be read from file
  const std::string soil_erod_fld_name = "mbl_bsn_fct_geo";

  // Dimensions of the field
  const std::string soil_erod_dname = "ncol";

  // initialize the file read
  soilErodibilityFunc::init_soil_erodibility_file_read(
      ncol_, soil_erod_fld_name, soil_erod_dname, grid_,
      soil_erodibility_data_file, srf_map_file, serod_horizInterp_,
      serod_dataReader_);  // output

  // -------------------------------------------------------------
  // Setup to enable reading marine organics file
  // -------------------------------------------------------------
  const std::string marine_organics_data_file =
      m_params.get<std::string>("marine_organics_file");

  // Fields to be read from file (order matters as they are read in the same
  // order)
  const std::vector<std::string> marine_org_fld_name = {
      "TRUEPOLYC", "TRUEPROTC", "TRUELIPC"};

  // Dimensions of the field
  const std::string marine_org_dname = "ncol";

  // initialize the file read
  marineOrganicsFunc::init_marine_organics_file_read(
      ncol_, marine_org_fld_name, marine_org_dname, grid_,
      marine_organics_data_file, srf_map_file,
      // output
      morg_horizInterp_, morg_data_start_, morg_data_end_, morg_data_out_,
      morg_dataReader_);

}  // set_grid ends

// ================================================================
//  REQUEST_BUFFER_SIZE_IN_BYTES
// ================================================================
// ON HOST, returns the number of bytes of device memory needed by the above
// Buffer type given the number of columns and vertical levels
size_t MAMSrfOnlineEmiss::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_);
}

// ================================================================
//  INIT_BUFFERS
// ================================================================
// ON HOST, initializes the Buffer type with sufficient memory to store
// intermediate (dry) quantities on the given number of columns with the given
// number of vertical levels. Returns the number of bytes allocated.
void MAMSrfOnlineEmiss::init_buffers(const ATMBufferManager &buffer_manager) {
  EKAT_REQUIRE_MSG(
      buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(),
      "Error! Insufficient buffer size.\n");

  size_t used_mem =
      mam_coupling::init_buffer(buffer_manager, ncol_, nlev_, buffer_);
  EKAT_REQUIRE_MSG(
      used_mem == requested_buffer_size_in_bytes(),
      "Error! Used memory != requested memory for MAMSrfOnlineEmiss.");
}

// ================================================================
//  INITIALIZE_IMPL
// ================================================================
void MAMSrfOnlineEmiss::initialize_impl(const RunType run_type) {
  // Check the interval values for the following fields used by this interface.
  // NOTE: We do not include aerosol and gas species, e.g., soa_a1, num_a1,
  // because we automatically added these fields.
  const std::map<std::string, std::pair<Real, Real>> ranges_emissions = {
      {"sst", {-1e10, 1e10}},  // FIXME
      {"dstflx", {-1e10, 1e10}}};
  set_ranges_process(ranges_emissions);
  add_interval_checks();

  // ---------------------------------------------------------------
  // Input fields read in from IC file, namelist or other processes
  // ---------------------------------------------------------------
  populate_wet_atm(wet_atm_);
  populate_dry_atm(dry_atm_, buffer_);

  // ---------------------------------------------------------------
  // Output fields
  // ---------------------------------------------------------------
  // Constituent fluxes of species in [kg/m2/s]
  constituent_fluxes_ = get_field_out("constituent_fluxes").get_view<Real **>();

  // ---------------------------------------------------------------
  // Allocate memory for local and work arrays
  // ---------------------------------------------------------------

  // Work array to store fluxes after unit conversions to kg/m2/s
  fluxes_in_mks_units_ = view_1d("fluxes_in_mks_units", ncol_);

  // Current month ( 0-based)
  const int curr_month = start_of_step_ts().get_month() - 1;

  // Load the first month into data_end.

  // Note: At the first time step, the data will be moved into data_beg,
  // and data_end will be reloaded from file with the new month.

  //--------------------------------------------------------------------
  // Update surface emissions from file
  //--------------------------------------------------------------------
  for(srf_emiss_ &ispec_srf : srf_emiss_species_) {
    srfEmissFunc::update_srfEmiss_data_from_file(
        ispec_srf.dataReader_, start_of_step_ts(), curr_month, *ispec_srf.horizInterp_,
        ispec_srf.data_end_);  // output
  }

  //-----------------------------------------------------------------
  // Read Soil erodibility data
  //-----------------------------------------------------------------
  // This data is time-independent, we read all data here for the
  // entire simulation
  soilErodibilityFunc::update_soil_erodibility_data_from_file(
      serod_dataReader_, *serod_horizInterp_,
      soil_erodibility_);  // output

  //--------------------------------------------------------------------
  // Update marine orgaincs from file
  //--------------------------------------------------------------------
  // Time dependent data
  marineOrganicsFunc::update_marine_organics_data_from_file(
      morg_dataReader_, start_of_step_ts(), curr_month, *morg_horizInterp_,
      morg_data_end_);  // output

  //-----------------------------------------------------------------
  // Setup preprocessing and post processing
  //-----------------------------------------------------------------
  preprocess_.initialize(ncol_, nlev_, wet_atm_, dry_atm_);

}  // end initialize_impl()

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMSrfOnlineEmiss::run_impl(const double dt) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  // Constituent fluxes [kg/m^2/s]
  auto constituent_fluxes = this->constituent_fluxes_;

  // Zero out constituent fluxes only for gasses and aerosols
  init_fluxes(ncol_,                // in
              constituent_fluxes);  // in-out
  Kokkos::fence();
  // Gather time and state information for interpolation
  const auto ts = end_of_step_ts();

  //--------------------------------------------------------------------
  // Online emissions from dust and sea salt
  //--------------------------------------------------------------------

  // --- Interpolate marine organics data --

  // Update TimeState, note the addition of dt
  morg_timeState_.t_now = ts.frac_of_year_in_days();

  // Update time state and if the month has changed, update the data.
  marineOrganicsFunc::update_marine_organics_timestate(
      morg_dataReader_, ts, *morg_horizInterp_,
      // output
      morg_timeState_, morg_data_start_, morg_data_end_);

  // Call the main marine organics routine to get interpolated forcings.
  marineOrganicsFunc::marineOrganics_main(morg_timeState_, morg_data_start_,
                                          morg_data_end_, morg_data_out_);

  // Marine organics emission data read from the file (order is important here)
  const const_view_1d mpoly = ekat::subview(morg_data_out_.emiss_sectors, 0);
  const const_view_1d mprot = ekat::subview(morg_data_out_.emiss_sectors, 1);
  const const_view_1d mlip  = ekat::subview(morg_data_out_.emiss_sectors, 2);

  // Ocean fraction [unitless]
  const const_view_1d ocnfrac =
      get_field_in("ocnfrac").get_view<const Real *>();

  // Sea surface temperature [K]
  const const_view_1d sst = get_field_in("sst").get_view<const Real *>();

  // U wind component [m/s]
  const const_view_2d u_wind =
      get_field_in("horiz_winds").get_component(0).get_view<const Real **>();

  // V wind component [m/s]
  const const_view_2d v_wind =
      get_field_in("horiz_winds").get_component(1).get_view<const Real **>();

  // Dust fluxes [kg/m^2/s]: Four flux values for each column
  const const_view_2d dstflx = get_field_in("dstflx").get_view<const Real **>();

  // Soil edodibility [fraction]
  const const_view_1d soil_erodibility = this->soil_erodibility_;

  // Vertical layer height at midpoints
  const const_view_2d z_mid = dry_atm_.z_mid;

  compute_online_dust_nacl_emiss(ncol_, nlev_, ocnfrac, sst, u_wind, v_wind,
                                 dstflx, mpoly, mprot, mlip, soil_erodibility,
                                 z_mid,
                                 // output
                                 constituent_fluxes);
  Kokkos::fence();
  //--------------------------------------------------------------------
  // Interpolate srf emiss data read in from emissions files
  //--------------------------------------------------------------------

  for(srf_emiss_ &ispec_srf : srf_emiss_species_) {
    // Update TimeState, note the addition of dt
    ispec_srf.timeState_.t_now = ts.frac_of_year_in_days();

    // Update time state and if the month has changed, update the data.
    srfEmissFunc::update_srfEmiss_timestate(
        ispec_srf.dataReader_, ts, *ispec_srf.horizInterp_,
        // output
        ispec_srf.timeState_, ispec_srf.data_start_, ispec_srf.data_end_);

    // Call the main srfEmiss routine to get interpolated aerosol forcings.
    srfEmissFunc::srfEmiss_main(ispec_srf.timeState_, ispec_srf.data_start_,
                                ispec_srf.data_end_, ispec_srf.data_out_);

    //--------------------------------------------------------------------
    // Modify units to MKS units (from molecules/cm2/s to kg/m2/s)
    //--------------------------------------------------------------------
    // Get species index in array with pcnst dimension (e.g., state_q or
    // constituent_fluxes_)
    const int species_index = spcIndex_in_pcnst_.at(ispec_srf.species_name);

    // modify units from molecules/cm2/s to kg/m2/s
    auto fluxes_in_mks_units = this->fluxes_in_mks_units_;
    const Real mfactor =
        amufac * mam4::gas_chemistry::adv_mass[species_index - offset_];
    const view_1d ispec_outdata0 =
        ekat::subview(ispec_srf.data_out_.emiss_sectors, 0);
    // Parallel loop over all the columns to update units
    Kokkos::parallel_for(
        "srf_emis_fluxes", ncol_, KOKKOS_LAMBDA(int icol) {
          fluxes_in_mks_units(icol) = ispec_outdata0(icol) * mfactor;
          constituent_fluxes(icol, species_index) = fluxes_in_mks_units(icol);
        });
  }  // for loop for species
  Kokkos::fence();
}  // run_impl ends
// =============================================================================
}  // namespace scream
