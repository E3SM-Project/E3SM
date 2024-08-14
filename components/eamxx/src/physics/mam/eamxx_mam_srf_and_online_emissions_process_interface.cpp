//#include <ekat/ekat_assert.hpp>
#include <physics/mam/eamxx_mam_srf_and_online_emissions_process_interface.hpp>

namespace scream {

// ================================================================
//  Constructor
// ================================================================
MAMSrfOnlineEmiss::MAMSrfOnlineEmiss(const ekat::Comm &comm,
                                     const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants.
   */
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

  static constexpr int pcnst = mam4::aero_model::pcnst;
  const FieldLayout scalar2d_pcnct =
      grid_->get_2d_vector_layout(pcnst, "num_phys_constituents");

  // -------------------------------------------------------------
  // These variables are "Computed" or outputs for the process
  // -------------------------------------------------------------
  static constexpr Units m2(m * m, "m2");
  // Constituent fluxes of species in [kg/m2/s]
  add_field<Computed>("constituent_fluxes", scalar2d_pcnct, kg / m2 / s,
                      grid_name);

  // Surface emissions remapping file
  auto srf_map_file = m_params.get<std::string>("srf_remap_file", "");

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
  }

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
  const int curr_month = timestamp().get_month() - 1;

  // Load the first month into data_end.

  // Note: At the first time step, the data will be moved into data_beg,
  // and data_end will be reloaded from file with the new month.

  //--------------------------------------------------------------------
  // Update surface emissions from file
  //--------------------------------------------------------------------
  for(srf_emiss_ &ispec_srf : srf_emiss_species_) {
    srfEmissFunc::update_srfEmiss_data_from_file(
        ispec_srf.dataReader_, timestamp(), curr_month, *ispec_srf.horizInterp_,
        ispec_srf.data_end_);  // output
  }

  //-----------------------------------------------------------------
  // Setup preprocessing and post processing
  //-----------------------------------------------------------------
  preprocess_.initialize(constituent_fluxes_);

}  // end initialize_impl()

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMSrfOnlineEmiss::run_impl(const double dt) {
  // Zero output
  Kokkos::deep_copy(preprocess_.constituent_fluxes_pre_, 0);

  // Gather time and state information for interpolation
  auto ts = timestamp() + dt;

  //--------------------------------------------------------------------
  // Interpolate srf emiss data
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
    auto constituent_fluxes  = this->constituent_fluxes_;
    const Real mfactor =
        amufac * mam4::gas_chemistry::adv_mass[species_index - offset_];
    // Parallel loop over all the columns to update units
    Kokkos::parallel_for(
        "fluxes", ncol_, KOKKOS_LAMBDA(int icol) {
          fluxes_in_mks_units(icol) =
              ispec_srf.data_out_.emiss_sectors(0, icol) * mfactor;
          constituent_fluxes(icol, species_index) = fluxes_in_mks_units(icol);
        });

  }  // for loop for species
  Kokkos::fence();
}  // run_imple ends

// =============================================================================
}  // namespace scream
