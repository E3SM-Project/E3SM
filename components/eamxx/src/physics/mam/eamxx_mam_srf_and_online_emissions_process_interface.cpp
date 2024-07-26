#include <ekat/ekat_assert.hpp>
#include <physics/mam/eamxx_mam_srf_and_online_emissions_process_interface.hpp>

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"

// for SCREAM_CIME_BUILD
#include "scream_config.h"

/*
Future work:
Write comments
write in/outs for all variables clearly
*/

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

  // Reading so2 srf emiss data
  std::string srf_map_file = "";
  std::string so2_data_file =
      "/compyfs/inputdata/atm/scream/mam4xx/emissions/test_DECK_ne30/"
      "cmip6_mam4_so2_surf_ne2np4_2010_clim_c20240723.nc";
  static constexpr int num_sectors_so2                = 6;
  std::array<std::string, num_sectors_so2> so2_fields = {"AGR", "RCO", "SHP",
                                                         "SLV", "TRA", "WST"};

  // Init horizontal remap
  so2SrfEmissHorizInterp_ = srfEmissFunc::create_horiz_remapper(
      grid_, so2_data_file, so2_fields, srf_map_file);

  // 2. Initialize the size of the SPAData structures.
  srfEmissData_start_ = srfEmissFunc::srfEmissInput(ncol_, num_sectors_so2);
  srfEmissData_end_   = srfEmissFunc::srfEmissInput(ncol_, num_sectors_so2);
  srfEmissData_out_.init(ncol_, num_sectors_so2,
                         true);  // FIXME: should it be true or false???

  // 3. Create reader for srfEmiss data. The reader is an
  //    AtmosphereInput object
  srfEmissDataReader_ = srfEmissFunc::create_srfEmiss_data_reader(
      so2SrfEmissHorizInterp_, so2_data_file);
}

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
    std::cout << "BALLI-bef:" << srfEmissData_end_.data.emiss_sectors.at(1)(i)
              << std::endl;
  }
  srfEmissFunc::update_srfEmiss_data_from_file(
      srfEmissDataReader_, timestamp(), curr_month, *so2SrfEmissHorizInterp_,
      srfEmissData_end_);
  for(int i = 19; i < 30; ++i) {
    std::cout << "BALLI:" << srfEmissData_end_.data.emiss_sectors[2](i) << ":"
              << i << std::endl;
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
  srfEmissTimeState_.t_now = ts.frac_of_year_in_days();
  // Update time state and if the month has changed, update the data.
  srfEmissFunc::update_srfEmiss_timestate(
      srfEmissDataReader_, ts, *so2SrfEmissHorizInterp_, srfEmissTimeState_,
      srfEmissData_start_, srfEmissData_end_);

  // Call the main srfEmiss routine to get interpolated aerosol forcings.
  srfEmissFunc::srfEmiss_main(srfEmissTimeState_, srfEmissData_start_,
                              srfEmissData_end_, srfEmiss_temp_,
                              srfEmissData_out_);
  for(int i = 19; i < 30; ++i) {
    std::cout << "BALLI:" << srfEmissData_out_.emiss_sectors[2](i) << ":" << i
              << std::endl;
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
