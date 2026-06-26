#include <physics/mam/eamxx_mam_srf_and_online_emissions_process_interface.hpp>

// For surface and online emission functions
#include <physics/mam/eamxx_mam_srf_and_online_emissions_functions.hpp>

#include "share/algorithm/eamxx_data_interpolation.hpp"

#include <ekat_team_policy_utils.hpp>
#include <iostream>

namespace scream {

// ================================================================
//  Constructor
// ================================================================
MAMSrfOnlineEmiss::MAMSrfOnlineEmiss(const ekat::Comm &comm,
                                     const ekat::ParameterList &params)
    : MAMGenericInterface(comm, params) {

  // FIXME: temporary solution to fix the sp test fails for the PR
  double dbl_dust_emis_scale_factor = m_params.get<double>("srf_emis_scale_factor_for_dust", 1.0);
  dust_emis_scale_factor = static_cast<Real>(dbl_dust_emis_scale_factor);

  double dbl_seasalt_emis_scale_factor = m_params.get<double>("srf_emis_scale_factor_for_seasalt", 1.0);
  seasalt_emis_scale_factor = static_cast<Real>(dbl_seasalt_emis_scale_factor);

  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants.
   */
  check_fields_intervals_ =
      m_params.get<bool>("create_fields_interval_checks", false);
}

// ================================================================
//  SET_GRIDS
// ================================================================
void MAMSrfOnlineEmiss::create_requests() {
  grid_                 = m_grids_manager->get_grid("physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();       // Number of columns on this rank
  nlev_ = grid_->get_num_vertical_levels();  // Number of levels per column

  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  constexpr auto m2     = pow(m, 2);
  constexpr auto s2     = pow(s, 2);

  const FieldLayout scalar2d   = grid_->get_2d_scalar_layout();
  const FieldLayout scalar3d_m = grid_->get_3d_scalar_layout(LEV);   // mid
  const FieldLayout scalar3d_i = grid_->get_3d_scalar_layout(ILEV);  // int

  // For U and V components of wind
  const FieldLayout vector3d = grid_->get_3d_vector_layout(LEV, 2);

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
  
  // cloud liquid number mixing ratio [1/kg]
  auto n_unit           = 1 / kg;   // units of number mixing ratios of tracers
  add_tracer<Required>("nc", grid_, n_unit);
  
  //----------- Variables from microphysics scheme -------------

  // Surface geopotential [m2/s2]
  add_field<Required>("phis", scalar2d, m2 / s2, grid_name);

  // -- Variables required for online dust and sea salt emissions --

  // U and V components of the wind[m/s]
  add_field<Required>("horiz_winds", vector3d, m / s, grid_name);

  //----------- Variables from coupler (ocean component)---------
  // Ocean fraction [unitless]
  add_field<Required>("ocnfrac", scalar2d, none, grid_name);

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

  // Soil erodibility [fraction]
  add_field<Computed>("soil_erodibility", scalar2d, none, grid_name);

  // Surface emissions remapping file
  auto srf_map_file = m_params.get<std::string>("srf_remap_file", "");

  // FIXME: We can extract the following info about each species
  // in a separate hpp file
  //--------------------------------------------------------------------
  // Init dms srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ dms;
  // File name, name and sectors
  dms.data_file    = m_params.get<std::string>("srf_emis_specifier_for_dms");
  dms.species_name = "dms";
  dms.sectors      = {"DMS"};
  dms.scale_factor = m_params.get<Real>("srf_emis_scale_factor_for_dms", 1.0);
  srf_emiss_species_.push_back(dms);  // add to the vector

  //--------------------------------------------------------------------
  // Init so2 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss_ so2;
  // File name, name and sectors
  so2.data_file    = m_params.get<std::string>("srf_emis_specifier_for_so2");
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
  // Register sector fields in FM for surface emissions.
  // DataInterpolation is set up in initialize_impl.
  //--------------------------------------------------------------------
  // -------------------------------------------------------------
  // Setup to enable reading marine organics file
  // -------------------------------------------------------------


  // Fields to be read from file (order matters as they are read in the same
  // order)
  const std::vector<std::string> marine_org_fld_name = {
      "TRUEPOLYC", "TRUEPROTC", "TRUELIPC"};
  for (const auto &field_name : marine_org_fld_name) {
    add_field<Computed>("morg_" + field_name, scalar2d, none, grid_name);
  }


}  // set_grid ends

// ================================================================
//  REQUEST_BUFFER_SIZE_IN_BYTES
// ================================================================
// ON HOST, returns the number of bytes of device memory needed by the above
// Buffer type given the number of columns and vertical levels
size_t MAMSrfOnlineEmiss::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_, 0, 0);
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

  // Get dust emission scheme from namelist
  auto dust_emis_scheme = m_params.get<int>("dust_emis_scheme", 1);

  // ---------------------------------------------------------------
  // Input fields read in from IC file, namelist or other processes
  // ---------------------------------------------------------------
  populate_wet_atm(wet_atm_);
  populate_dry_atm(dry_atm_, buffer_);

  // ---------------------------------------------------------------
  // Output fields
  // ---------------------------------------------------------------
  // Constituent fluxes of species in [kg/m2/s]
  constituent_fluxes_ = get_field_out("constituent_fluxes");

    const std::string marine_organics_data_file =
      m_params.get<std::string>("marine_organics_file");
  const auto marine_map_file = m_params.get<std::string>("srf_remap_file", "");

     std::cout << "MAMSrfOnlineEmiss: get_field_out marine_organics_data_files "
             << std::endl;
   const std::vector<std::string> marine_org_fld_name = {
      "TRUEPOLYC", "TRUEPROTC", "TRUELIPC"};           
  morg_fields_={};
  morg_fields_.reserve(marine_org_fld_name.size());
  for (const auto &field_name : marine_org_fld_name) {
    morg_fields_.push_back(
        get_field_out("morg_" + field_name).alias(field_name));
  }

  std::cout << "MAMSrfOnlineEmiss: before DataInterpolation ctor for marine organics"
            << std::endl;
  morg_data_interp_ = std::make_shared<DataInterpolation>(grid_, morg_fields_);
  std::cout << "MAMSrfOnlineEmiss: after DataInterpolation ctor for marine organics"
            << std::endl;
  std::cout << "MAMSrfOnlineEmiss: before set_logger (marine organics)"
            << std::endl;
  morg_data_interp_->set_logger(m_atm_logger);
  std::cout << "MAMSrfOnlineEmiss: after set_logger (marine organics)"
            << std::endl;
  std::cout << "MAMSrfOnlineEmiss: before setup_periodic_time_database (marine organics)"
            << std::endl;
  morg_data_interp_->setup_periodic_time_database({marine_organics_data_file});
  std::cout << "MAMSrfOnlineEmiss: after setup_periodic_time_database (marine organics)"
            << std::endl;
  std::cout << "MAMSrfOnlineEmiss: before create_horiz_remappers (marine organics)"
            << std::endl;
  morg_data_interp_->create_horiz_remappers(
      marine_map_file == "none" ? "" : marine_map_file);
  std::cout << "MAMSrfOnlineEmiss: after create_horiz_remappers (marine organics)"
            << std::endl;
  DataInterpolation::VertRemapData remap_data;
  remap_data.vr_type = DataInterpolation::None;
  std::cout << "MAMSrfOnlineEmiss: before create_vert_remapper (marine organics)"
            << std::endl;
  morg_data_interp_->create_vert_remapper(remap_data);
  std::cout << "MAMSrfOnlineEmiss: after create_vert_remapper (marine organics)"
            << std::endl;
  //--------------------------------------------------------------------
  // Setup data interpolation for surface emissions.
  //--------------------------------------------------------------------
  {
    using namespace ekat::units;
    using namespace ShortFieldTagsNames;
    const FieldLayout scalar2d = grid_->get_2d_scalar_layout();
    const auto srf_map_file    = m_params.get<std::string>("srf_remap_file", "");
    const auto srf_time_interp = DataInterpolation::Linear;
    for(srf_emiss_ &ispec_srf : srf_emiss_species_) {
      std::vector<Field> srf_fields;
      srf_fields.reserve(ispec_srf.sectors.size());
      for(const auto &sector_name : ispec_srf.sectors) {
        Field field(FieldIdentifier(sector_name, scalar2d, none, grid_->name()));
        field.allocate_view();
        srf_fields.push_back(field);
      }
      ispec_srf.emiss_sector_fields_ = srf_fields;

      ispec_srf.data_interp_ = std::make_shared<DataInterpolation>(grid_, srf_fields);
      ispec_srf.data_interp_->set_logger(m_atm_logger);
      ispec_srf.data_interp_->setup_periodic_time_database(
          {ispec_srf.data_file});
      ispec_srf.data_interp_->create_horiz_remappers(
          srf_map_file == "none" ? "" : srf_map_file);

      DataInterpolation::VertRemapData remap_data;
      remap_data.vr_type = DataInterpolation::None;
      ispec_srf.data_interp_->create_vert_remapper(remap_data);

      ispec_srf.data_interp_->init_time_interpolation(start_of_step_ts(), srf_time_interp);
    }
  }

    // Current month ( 0-based)
    const int curr_month = start_of_step_ts().get_month() - 1;

  //-----------------------------------------------------------------
  // Read Soil erodibility data
  //-----------------------------------------------------------------
  if (dust_emis_scheme == 1) {
    // This data is time-independent, we read all data here for the
    // entire simulation   
    const auto soil_erodibility_data_file =
        m_params.get<std::string>("soil_erodibility_file");
    const auto srf_map_file = m_params.get<std::string>("srf_remap_file", "");
    const std::string soil_erod_fld_name = "mbl_bsn_fct_geo";

    std::vector<Field> soil_erod_fields = {
        get_field_out("soil_erodibility").alias(soil_erod_fld_name)};
    std::cout << "MAMSrfOnlineEmiss: before DataInterpolation ctor for soil erodibility"
              << std::endl;
    auto soil_erod_data_interp =
        std::make_shared<DataInterpolation>(grid_, soil_erod_fields);
    std::cout << "MAMSrfOnlineEmiss: after DataInterpolation ctor for soil erodibility"
              << std::endl;
    std::cout << "MAMSrfOnlineEmiss: before set_logger (soil erodibility)"
              << std::endl;
    soil_erod_data_interp->set_logger(m_atm_logger);
    std::cout << "MAMSrfOnlineEmiss: after set_logger (soil erodibility)"
              << std::endl;
    std::cout << "MAMSrfOnlineEmiss: before setup_static_database (soil erodibility)"
              << std::endl;
    soil_erod_data_interp->setup_static_database({soil_erodibility_data_file});
    std::cout << "MAMSrfOnlineEmiss: after setup_static_database (soil erodibility)"
              << std::endl;
    std::cout << "MAMSrfOnlineEmiss: before create_horiz_remappers (soil erodibility)"
              << std::endl;
    soil_erod_data_interp->create_horiz_remappers(
        srf_map_file == "none" ? "" : srf_map_file);
    std::cout << "MAMSrfOnlineEmiss: after create_horiz_remappers (soil erodibility)"
              << std::endl;
    DataInterpolation::VertRemapData remap_data;
    remap_data.vr_type = DataInterpolation::None;
    std::cout << "MAMSrfOnlineEmiss: before create_vert_remapper (soil erodibility)"
              << std::endl;
    soil_erod_data_interp->create_vert_remapper(remap_data);
    std::cout << "MAMSrfOnlineEmiss: after create_vert_remapper (soil erodibility)"
              << std::endl;
    std::cout << "MAMSrfOnlineEmiss: before run static interpolation (soil erodibility)"
              << std::endl;
    soil_erod_data_interp->run();
    std::cout << "MAMSrfOnlineEmiss: after run static interpolation (soil erodibility)"
              << std::endl;

    soil_erodibility_ = soil_erodibility_field_.get_view<const Real *>();
  } else if (dust_emis_scheme == 2) {
    // For dust emission scheme 2, override soil erodibility to 1
    auto soil_erod_ones = view_1d("soil_erod_ones", ncol_);
    Kokkos::deep_copy(soil_erod_ones, 1.0);
    soil_erodibility_ = soil_erod_ones;
  }

  //--------------------------------------------------------------------
  // Update marine orgaincs from file
  //--------------------------------------------------------------------
  std::cout << "MAMSrfOnlineEmiss: before init_time_interpolation (marine organics)"
            << std::endl;
  morg_data_interp_->init_time_interpolation(start_of_step_ts(),
                                             DataInterpolation::Linear);
  std::cout << "MAMSrfOnlineEmiss: after init_time_interpolation (marine organics)"
            << std::endl;

  //-----------------------------------------------------------------
  // Setup preprocessing and post processing
  //-----------------------------------------------------------------
  preprocess_.initialize(ncol_, nlev_, wet_atm_, dry_atm_);

}  // end initialize_impl()

// ================================================================
//  RUN_IMPL
// ================================================================
void MAMSrfOnlineEmiss::run_impl(const double dt) {
  using TPF = ekat::TeamPolicyFactory<KT::ExeSpace>;
  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  // Constituent fluxes [kg/m^2/s]
  auto constituent_fluxes = constituent_fluxes_.get_view<Real **>();

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
  std::cout << "MAMSrfOnlineEmiss: before run(ts) (marine organics)" << std::endl;
  morg_data_interp_->run(ts);
  std::cout << "MAMSrfOnlineEmiss: after run(ts) (marine organics)" << std::endl;

  // Marine organics emission data read from the file (order is important here)
  const const_view_1d mpoly = morg_fields_[0].get_view<const Real *>();
  const const_view_1d mprot = morg_fields_[1].get_view<const Real *>();
  const const_view_1d mlip  = morg_fields_[2].get_view<const Real *>();

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
                                 dust_emis_scale_factor,
                                 seasalt_emis_scale_factor,
                                 // output
                                 constituent_fluxes);
  Kokkos::fence();
  //--------------------------------------------------------------------
  // Interpolate srf emiss data read in from emissions files
  //--------------------------------------------------------------------

  for(srf_emiss_ &ispec_srf : srf_emiss_species_) {
        ispec_srf.data_interp_->run(ts);

    //--------------------------------------------------------------------
    // Modify units to MKS units (from molecules/cm2/s to kg/m2/s)
    //--------------------------------------------------------------------
    // Get species index in array with pcnst dimension (e.g., state_q or
    // constituent_fluxes_)
    const int species_index = spcIndex_in_pcnst_.at(ispec_srf.species_name);

    auto constituent_fluxes_ispe_srf = constituent_fluxes_.get_component(species_index);
    // modify units from molecules/cm2/s to kg/m2/s
    constituent_fluxes_ispe_srf.deep_copy(0.0);

    for(const auto &sector_field : ispec_srf.emiss_sector_fields_) {
        constituent_fluxes_ispe_srf.update(sector_field, 1, 1);
    }

        const Real mfactor = amufac * ispec_srf.scale_factor *
                                                 mam4::gas_chemistry::adv_mass[species_index - offset_];
    constituent_fluxes_ispe_srf.scale(mfactor);
  }  // for loop for species
  Kokkos::fence();
}  // run_impl ends
// =============================================================================
}  // namespace scream
