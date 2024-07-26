#include <ekat/ekat_assert.hpp>
#include <physics/mam/eamxx_mam_srf_and_online_emissions_process_interface.hpp>

namespace scream {

// =========================================================================================
MAMSrfOnlineEmiss::MAMSrfOnlineEmiss(const ekat::Comm &comm,
                                     const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  /* Anything that can be initialized without grid information can be
   * initialized here. Like universal constants, mam wetscav options.
   */
}

// =========================================================================================
void MAMSrfOnlineEmiss::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  grid_                 = grids_manager->get_grid("Physics");
  const auto &grid_name = grid_->name();

  ncol_ = grid_->get_num_local_dofs();  // Number of columns on this rank

  using namespace ekat::units;
  FieldLayout scalar3d_mid = grid_->get_3d_scalar_layout(true);

  // Temperature[K] at midpoints
  add_field<Required>("T_mid", scalar3d_mid, K, grid_name);

  // Surface emissions remapping file
  std::string srf_map_file = m_params.get<std::string>("srf_remap_file");

  //--------------------------------------------------------------------
  // Init so2 srf emiss data structures
  //--------------------------------------------------------------------
  // File name
  std::string so2_data_file =
      m_params.get<std::string>("srf_emis_specifier_for_SO2");

  // Number of sectors
  static constexpr int so2_num_sectors = 6;

  // Sector names in file
  std::array<std::string, so2_num_sectors> so2_sectors = {"AGR", "RCO", "SHP",
                                                          "SLV", "TRA", "WST"};

  srfEmissFunc::init_srf_emiss_objects(
      ncol_, so2_num_sectors, grid_, so2_data_file, so2_sectors, srf_map_file,
      // output
      so2SrfEmissHorizInterp_, so2SrfEmissData_start_, so2SrfEmissData_end_,
      so2SrfEmissData_out_, so2SrfEmissDataReader_);

}  // set_grid

// =========================================================================================
// ON HOST, returns the number of bytes of device memory needed by the above
// Buffer type given the number of columns and vertical levels
size_t MAMSrfOnlineEmiss::requested_buffer_size_in_bytes() const {
  return mam_coupling::buffer_size(ncol_, nlev_);
}

// =========================================================================================
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
// =========================================================================================
// // TODO: comments!
// void MAMSrfOnlineEmiss::set_emissions_names() {} \\ end set_emissions_names()
// =========================================================================================
// inline void set_emissions_layouts() {

// } // end set_emissions_layouts()
// =========================================================================================
void MAMSrfOnlineEmiss::initialize_impl(const RunType run_type) {
  // Load the first month into srfEmiss_end.
  // Note: At the first time step, the data will be moved into srfEmiss_beg,
  //       and srfEmiss_end will be reloaded from file with the new month.
  const int curr_month = timestamp().get_month() - 1;  // 0-based
  for(int i = 19; i < 30; ++i) {
    std::cout << "BALLI-bef:"
              << so2SrfEmissData_end_.data.emiss_sectors.at(1)(i) << std::endl;
  }
  srfEmissFunc::update_srfEmiss_data_from_file(
      so2SrfEmissDataReader_, timestamp(), curr_month, *so2SrfEmissHorizInterp_,
      so2SrfEmissData_end_);
  for(int i = 19; i < 30; ++i) {
    std::cout << "BALLI:" << so2SrfEmissData_end_.data.emiss_sectors[2](i)
              << ":" << i << std::endl;
  }

}  // end initialize_impl()

// =============================================================================
void MAMSrfOnlineEmiss::run_impl(const double dt) {
  const auto scan_policy = ekat::ExeSpaceUtils<
      KT::ExeSpace>::get_thread_range_parallel_scan_team_policy(ncol_, nlev_);

  // preprocess input -- needs a scan for the calculation of atm height
  // Kokkos::parallel_for("preprocess", scan_policy, preprocess_);
  // Kokkos::fence();

  // Gather time and state information for interpolation
  auto ts = timestamp() + dt;
  // Update the srfEmissTimeState to reflect the current time, note the addition
  // of dt
  so2SrfEmissTimeState_.t_now = ts.frac_of_year_in_days();
  // Update time state and if the month has changed, update the data.
  srfEmissFunc::update_srfEmiss_timestate(
      so2SrfEmissDataReader_, ts, *so2SrfEmissHorizInterp_,
      so2SrfEmissTimeState_, so2SrfEmissData_start_, so2SrfEmissData_end_);

  // Call the main srfEmiss routine to get interpolated aerosol forcings.
  srfEmissFunc::srfEmiss_main(so2SrfEmissTimeState_, so2SrfEmissData_start_,
                              so2SrfEmissData_end_, so2SrfEmissData_out_);
  for(int i = 19; i < 30; ++i) {
    std::cout << "BALLI:" << so2SrfEmissData_out_.emiss_sectors[2](i) << ":"
              << i << std::endl;
  }

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
