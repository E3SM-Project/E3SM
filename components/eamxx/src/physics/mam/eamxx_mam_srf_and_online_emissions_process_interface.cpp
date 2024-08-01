#include <ekat/ekat_assert.hpp>
#include <physics/mam/eamxx_mam_srf_and_online_emissions_process_interface.hpp>

namespace scream {

// ================================================================
//  Constructor
// ================================================================
MAMSrfOnlineEmiss::MAMSrfOnlineEmiss(const ekat::Comm &comm,
                                     const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}

// ================================================================
//  SET_GRIDS
// ================================================================
void MAMSrfOnlineEmiss::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  grid_                 = grids_manager->get_grid("Physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();  // Number of columns on this rank

  using namespace ekat::units;
  FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);

  static constexpr int pcnst = mam4::aero_model::pcnst;

  const FieldLayout scalar2d_pcnct =
      grid_->get_2d_vector_layout(pcnst, "num_phys_constituents");

  // --------------------------------------------------------------------------
  // These variables are "Required" or pure inputs for the process
  // --------------------------------------------------------------------------
  // Temperature[K] at midpoints
  add_field<Required>("T_mid", scalar3d_mid, K, grid_name);

  // -------------------------------------------------------------
  // These variables are "Computed" or outputs for the process
  // -------------------------------------------------------------
  static constexpr Units m2(m * m, "m2");
  add_field<Computed>("constituent_fluxes", scalar2d_pcnct, kg / m2 / s,
                      grid_name);

  // Surface emissions remapping file
  std::string srf_map_file = m_params.get<std::string>("srf_remap_file");

  //--------------------------------------------------------------------
  // Init dms srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss dms;

  // File name and sectors
  dms.data_file = m_params.get<std::string>("srf_emis_specifier_for_DMS");
  dms.sectors   = {"DMS"};
  srf_emiss_species_.push_back(dms);  // add to the vector

  //--------------------------------------------------------------------
  // Init so2 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss so2;
  // File name and sectors
  so2.data_file = m_params.get<std::string>("srf_emis_specifier_for_SO2");
  so2.sectors   = {"AGR", "RCO", "SHP", "SLV", "TRA", "WST"};
  srf_emiss_species_.push_back(so2);  // add to the vector
  //--------------------------------------------------------------------
  // Init bc_a4 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss bc_a4;
  // File name and sectors
  bc_a4.data_file = m_params.get<std::string>("srf_emis_specifier_for_bc_a4");
  bc_a4.sectors   = {"AGR", "ENE", "IND", "RCO", "SHP", "SLV", "TRA", "WST"};
  srf_emiss_species_.push_back(bc_a4);  // add to the vector

  //--------------------------------------------------------------------
  // Init num_a1 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss num_a1;
  // File name and sectors
  num_a1.data_file = m_params.get<std::string>("srf_emis_specifier_for_num_a1");
  num_a1.sectors   = {"num_a1_SO4_AGR", "num_a1_SO4_SHP", "num_a1_SO4_SLV",
                      "num_a1_SO4_WST"};
  srf_emiss_species_.push_back(num_a1);  // add to the vector

  //--------------------------------------------------------------------
  // Init num_a2 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss num_a2;
  // File name and sectors
  num_a2.data_file = m_params.get<std::string>("srf_emis_specifier_for_num_a2");
  num_a2.sectors   = {"num_a2_SO4_RCO", "num_a2_SO4_TRA"};
  srf_emiss_species_.push_back(num_a2);  // add to the vector

  //--------------------------------------------------------------------
  // Init num_a4 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss num_a4;
  // File name and sectors
  num_a4.data_file = m_params.get<std::string>("srf_emis_specifier_for_num_a4");
  num_a4.sectors   = {
        "num_a1_BC_AGR",  "num_a1_BC_ENE",  "num_a1_BC_IND",  "num_a1_BC_RCO",
        "num_a1_BC_SHP",  "num_a1_BC_SLV",  "num_a1_BC_TRA",  "num_a1_BC_WST",
        "num_a1_POM_AGR", "num_a1_POM_ENE", "num_a1_POM_IND", "num_a1_POM_RCO",
        "num_a1_POM_SHP", "num_a1_POM_SLV", "num_a1_POM_TRA", "num_a1_POM_WST"};
  srf_emiss_species_.push_back(num_a4);  // add to the vector

  //--------------------------------------------------------------------
  // Init pom_a4 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss pom_a4;
  // File name and sectors
  pom_a4.data_file = m_params.get<std::string>("srf_emis_specifier_for_pom_a4");
  pom_a4.sectors   = {"AGR", "ENE", "IND", "RCO", "SHP", "SLV", "TRA", "WST"};
  srf_emiss_species_.push_back(pom_a4);  // add to the vector

  //--------------------------------------------------------------------
  // Init so4_a1 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss so4_a1;
  // File name and sectors
  so4_a1.data_file = m_params.get<std::string>("srf_emis_specifier_for_so4_a1");
  so4_a1.sectors   = {"AGR", "SHP", "SLV", "WST"};
  srf_emiss_species_.push_back(so4_a1);

  //--------------------------------------------------------------------
  // Init so4_a2 srf emiss data structures
  //--------------------------------------------------------------------
  srf_emiss so4_a2;
  // File name and sectors
  so4_a2.data_file = m_params.get<std::string>("srf_emis_specifier_for_so4_a2");
  so4_a2.sectors   = {"RCO", "TRA"};
  srf_emiss_species_.push_back(so4_a2);

  //--------------------------------------------------------------------
  // Init data structures to read and interpolate
  //--------------------------------------------------------------------

  for(const sr_emiss &ispec_srf : srf_emiss_species_) {
    srfEmissFunc::init_srf_emiss_objects(
        ncol_, grid_, ispec_srf.data_file, ispec_srf.sectors, srf_map_file,
        // output
        ispec_srf.HorizInterp_, ispec_srf.Data_start_, ispec_srf.Data_end_,
        ispec_srf.Data_out_, ispec_srf.DataReader_);
  }

}  // set_grid

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
  constituent_fluxes_ = get_field_out("constituent_fluxes").get_view<Real **>();

  // Current month ( 0-based)
  const int curr_month = timestamp().get_month() - 1;

  // Load the first month into srfEmiss_end.

  // Note: At the first time step, the data will be moved into srfEmiss_beg,
  // and srfEmiss_end will be reloaded from file with the new month.

  //--------------------------------------------------------------------
  // Update surface emissions from file
  //--------------------------------------------------------------------
  for(auto &ispec_srf : srf_emiss_species_) {
    srfEmissFunc::update_srfEmiss_data_from_file(
        ispec_srf.DataReader_, timestamp(), curr_month, *ispec_srf.HorizInterp_,
        ispec_srf.Data_end_);
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
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  Kokkos::fence();

  // Gather time and state information for interpolation
  auto ts = timestamp() + dt;

  //--------------------------------------------------------------------
  // Interpolate dms srf emiss data
  //--------------------------------------------------------------------

  for(auto &ispec_srf : srf_emiss_species_) {
    // Update srfEmissTimeState, note the addition of dt
    ispec_srf.TimeState_.t_now = ts.frac_of_year_in_days();

    // Update time state and if the month has changed, update the data.
    srfEmissFunc::update_srfEmiss_timestate(
        ispec_srf.DataReader_, ts, *ispec_srf.HorizInterp_,
        ispec_srf.TimeState_, ispec_srf.Data_start_, ispec_srf.Data_end_);

    // Call the main srfEmiss routine to get interpolated aerosol forcings.
    srfEmissFunc::srfEmiss_main(ispec_srf.TimeState_, ispec_srf.Data_start_,
                                ispec_srf.Data_end_, ispec_srf.Data_out_);
    // update flux
    auto clfx = Kokkos::subview(constituent_fluxes_, Kokkos::ALL(),
                                static_cast<int>(spcIndex_in_pcnst::DMS));
    Kokkos::deep_copy(cflx, dmsSrfEmissData_out_.emiss_sectors[0]);
  }
  /*
// update flux
auto constituent_fluxes_DMS =
    Kokkos::subview(constituent_fluxes_, Kokkos::ALL(),
                    static_cast<int>(spcIndex_in_pcnst::DMS));
Kokkos::deep_copy(constituent_fluxes_DMS,
                  dmsSrfEmissData_out_.emiss_sectors[0]);

//--------------------------------------------------------------------
// Interpolate so2 srf emiss data
//--------------------------------------------------------------------
// Update srfEmissTimeState, note the addition of dt
so2SrfEmissTimeState_.t_now = ts.frac_of_year_in_days();

// Update time state and if the month has changed, update the data.
srfEmissFunc::update_srfEmiss_timestate(
    so2SrfEmissDataReader_, ts, *so2SrfEmissHorizInterp_,
    so2SrfEmissTimeState_, so2SrfEmissData_start_, so2SrfEmissData_end_);

// Call the main srfEmiss routine to get interpolated aerosol forcings.
srfEmissFunc::srfEmiss_main(so2SrfEmissTimeState_, so2SrfEmissData_start_,
                            so2SrfEmissData_end_, so2SrfEmissData_out_);
// update flux
auto constituent_fluxes_SO2 =
    Kokkos::subview(constituent_fluxes_, Kokkos::ALL(),
                    static_cast<int>(spcIndex_in_pcnst::SO2));
Kokkos::deep_copy(constituent_fluxes_SO2,
                  so2SrfEmissData_out_.emiss_sectors[0]);

//--------------------------------------------------------------------
// Interpolate bc_a4 srf emiss data
//--------------------------------------------------------------------
// Update srfEmissTimeState, note the addition of dt
bc_a4SrfEmissTimeState_.t_now = ts.frac_of_year_in_days();

// Update time state and if the month has changed, update the data.
srfEmissFunc::update_srfEmiss_timestate(
    bc_a4SrfEmissDataReader_, ts, *bc_a4SrfEmissHorizInterp_,
    bc_a4SrfEmissTimeState_, bc_a4SrfEmissData_start_,
    bc_a4SrfEmissData_end_);

// Call the main srfEmiss routine to get interpolated aerosol forcings.
srfEmissFunc::srfEmiss_main(bc_a4SrfEmissTimeState_, bc_a4SrfEmissData_start_,
                            bc_a4SrfEmissData_end_, bc_a4SrfEmissData_out_);
// update flux
auto constituent_fluxes_bc_a4 =
    Kokkos::subview(constituent_fluxes_, Kokkos::ALL(),
                    static_cast<int>(spcIndex_in_pcnst::bc_a4));
Kokkos::deep_copy(constituent_fluxes_bc_a4,
                  bc_a4SrfEmissData_out_.emiss_sectors[0]);

//--------------------------------------------------------------------
// Interpolate num_a1 srf emiss data
//--------------------------------------------------------------------
// Update srfEmissTimeState, note the addition of dt
num_a1SrfEmissTimeState_.t_now = ts.frac_of_year_in_days();

// Update time state and if the month has changed, update the data.
srfEmissFunc::update_srfEmiss_timestate(
    num_a1SrfEmissDataReader_, ts, *num_a1SrfEmissHorizInterp_,
    num_a1SrfEmissTimeState_, num_a1SrfEmissData_start_,
    num_a1SrfEmissData_end_);

// Call the main srfEmiss routine to get interpolated aerosol forcings.
srfEmissFunc::srfEmiss_main(num_a1SrfEmissTimeState_,
                            num_a1SrfEmissData_start_,
                            num_a1SrfEmissData_end_, num_a1SrfEmissData_out_);
// update flux
auto constituent_fluxes_num_a1 =
    Kokkos::subview(constituent_fluxes_, Kokkos::ALL(),
                    static_cast<int>(spcIndex_in_pcnst::num_a1));
Kokkos::deep_copy(constituent_fluxes_num_a1,
                  num_a1SrfEmissData_out_.emiss_sectors[0]);

//--------------------------------------------------------------------
// Interpolate num_a2 srf emiss data
//--------------------------------------------------------------------
// Update srfEmissTimeState, note the addition of dt
num_a2SrfEmissTimeState_.t_now = ts.frac_of_year_in_days();

// Update time state and if the month has changed, update the data.
srfEmissFunc::update_srfEmiss_timestate(
    num_a2SrfEmissDataReader_, ts, *num_a2SrfEmissHorizInterp_,
    num_a2SrfEmissTimeState_, num_a2SrfEmissData_start_,
    num_a2SrfEmissData_end_);

// Call the main srfEmiss routine to get interpolated aerosol forcings.
srfEmissFunc::srfEmiss_main(num_a2SrfEmissTimeState_,
                            num_a2SrfEmissData_start_,
                            num_a2SrfEmissData_end_, num_a2SrfEmissData_out_);
// update flux
auto constituent_fluxes_num_a2 =
    Kokkos::subview(constituent_fluxes_, Kokkos::ALL(),
                    static_cast<int>(spcIndex_in_pcnst::num_a2));
Kokkos::deep_copy(constituent_fluxes_num_a2,
                  num_a2SrfEmissData_out_.emiss_sectors[0]);

//--------------------------------------------------------------------
// Interpolate num_a4 srf emiss data
//--------------------------------------------------------------------
// Update srfEmissTimeState, note the addition of dt
num_a4SrfEmissTimeState_.t_now = ts.frac_of_year_in_days();

// Update time state and if the month has changed, update the data.
srfEmissFunc::update_srfEmiss_timestate(
    num_a4SrfEmissDataReader_, ts, *num_a4SrfEmissHorizInterp_,
    num_a4SrfEmissTimeState_, num_a4SrfEmissData_start_,
    num_a4SrfEmissData_end_);

// Call the main srfEmiss routine to get interpolated aerosol forcings.
srfEmissFunc::srfEmiss_main(num_a4SrfEmissTimeState_,
                            num_a4SrfEmissData_start_,
                            num_a4SrfEmissData_end_, num_a4SrfEmissData_out_);
// update flux
auto constituent_fluxes_num_a4 =
    Kokkos::subview(constituent_fluxes_, Kokkos::ALL(),
                    static_cast<int>(spcIndex_in_pcnst::num_a4));
Kokkos::deep_copy(constituent_fluxes_num_a4,
                  num_a4SrfEmissData_out_.emiss_sectors[0]);

//--------------------------------------------------------------------
// Interpolate pom_a4 srf emiss data
//--------------------------------------------------------------------
// Update srfEmissTimeState, note the addition of dt
pom_a4SrfEmissTimeState_.t_now = ts.frac_of_year_in_days();

// Update time state and if the month has changed, update the data.
srfEmissFunc::update_srfEmiss_timestate(
    pom_a4SrfEmissDataReader_, ts, *pom_a4SrfEmissHorizInterp_,
    pom_a4SrfEmissTimeState_, pom_a4SrfEmissData_start_,
    pom_a4SrfEmissData_end_);

// Call the main srfEmiss routine to get interpolated aerosol forcings.
srfEmissFunc::srfEmiss_main(pom_a4SrfEmissTimeState_,
                            pom_a4SrfEmissData_start_,
                            pom_a4SrfEmissData_end_, pom_a4SrfEmissData_out_);
// update flux
auto constituent_fluxes_pom_a4 =
    Kokkos::subview(constituent_fluxes_, Kokkos::ALL(),
                    static_cast<int>(spcIndex_in_pcnst::pom_a4));
Kokkos::deep_copy(constituent_fluxes_pom_a4,
                  pom_a4SrfEmissData_out_.emiss_sectors[0]);

//--------------------------------------------------------------------
// Interpolate so4_a1 srf emiss data
//--------------------------------------------------------------------
// Update srfEmissTimeState, note the addition of dt
so4_a1SrfEmissTimeState_.t_now = ts.frac_of_year_in_days();

// Update time state and if the month has changed, update the data.
srfEmissFunc::update_srfEmiss_timestate(
    so4_a1SrfEmissDataReader_, ts, *so4_a1SrfEmissHorizInterp_,
    so4_a1SrfEmissTimeState_, so4_a1SrfEmissData_start_,
    so4_a1SrfEmissData_end_);

// Call the main srfEmiss routine to get interpolated aerosol forcings.
srfEmissFunc::srfEmiss_main(so4_a1SrfEmissTimeState_,
                            so4_a1SrfEmissData_start_,
                            so4_a1SrfEmissData_end_, so4_a1SrfEmissData_out_);
// update flux
auto constituent_fluxes_so4_a1 =
    Kokkos::subview(constituent_fluxes_, Kokkos::ALL(),
                    static_cast<int>(spcIndex_in_pcnst::so4_a1));
Kokkos::deep_copy(constituent_fluxes_so4_a1,
                  so4_a1SrfEmissData_out_.emiss_sectors[0]);

//--------------------------------------------------------------------
// Interpolate so4_a2 srf emiss data
//--------------------------------------------------------------------
// Update srfEmissTimeState, note the addition of dt
so4_a2SrfEmissTimeState_.t_now = ts.frac_of_year_in_days();

// Update time state and if the month has changed, update the data.
srfEmissFunc::update_srfEmiss_timestate(
    so4_a2SrfEmissDataReader_, ts, *so4_a2SrfEmissHorizInterp_,
    so4_a2SrfEmissTimeState_, so4_a2SrfEmissData_start_,
    so4_a2SrfEmissData_end_);

// Call the main srfEmiss routine to get interpolated aerosol forcings.
srfEmissFunc::srfEmiss_main(so4_a2SrfEmissTimeState_,
                            so4_a2SrfEmissData_start_,
                            so4_a2SrfEmissData_end_, so4_a2SrfEmissData_out_);
// update flux
auto constituent_fluxes_so4_a2 =
    Kokkos::subview(constituent_fluxes_, Kokkos::ALL(),
                    static_cast<int>(spcIndex_in_pcnst::so4_a2));
Kokkos::deep_copy(constituent_fluxes_so4_a2,
                  so4_a2SrfEmissData_out_.emiss_sectors[0]);

for(int i = 19; i < 30; ++i) {
  std::cout << "BALLI:" << so4_a2SrfEmissData_out_.emiss_sectors[0](i) << ":"
            << i << ":" << mam4::gas_chemistry::adv_mass[21 - offset_]
            << std::endl;
}*/

  /* Rough notes:

  Here we should implement or port the chem_emissions subroutine in
  chemistry.F90. Basically call two subroutines, aero_model_emissions and
  set_srf_emissions.

  Here is the code:

    ! initialize chemistry constituent surface fluxes to zero
    do m = 2,pcnst
       n = map2chm(m)
       if (n>0) cam_in%cflx(:,m) = 0._r8
    enddo

    ! aerosol emissions ...
    call aero_model_emissions( state, & ! in
                               cam_in ) ! out

    ! prescribed emissions from file ...

    !-----------------------------------------------------------------------
    !        ... Set surface emissions
    !-----------------------------------------------------------------------
    call set_srf_emissions( lchnk, ncol, sflx(:,:) )

    do m = 1,pcnst
       n = map2chm(m)
       if ( n /= h2o_ndx .and. n > 0 ) then
          cam_in%cflx(:ncol,m) = cam_in%cflx(:ncol,m) + sflx(:ncol,n)
          call outfld( sflxnam(m), cam_in%cflx(:ncol,m), ncol,lchnk )
       endif
    enddo


  */

  std::cout << "End of surface emissions run" << std::endl;
}

// =============================================================================
}  // namespace scream
