#ifndef MARINE_ORGANICS_HPP
#define MARINE_ORGANICS_HPP

// For AtmosphereInput
#include "share/io/scorpio_input.hpp"

namespace scream {
namespace marine_organics {

template <typename ScalarType, typename DeviceType>
struct marineOrganicsFunctions {
  using Device = DeviceType;

  using KT      = KokkosTypes<Device>;
  using view_2d = typename KT::template view_2d<Real>;

  struct marineOrganicsTimeState {
    marineOrganicsTimeState() = default;
    // Whether the timestate has been initialized.
    // The current month
    int current_month = -1;
    // Julian Date for the beginning of the month, as defined in
    //           /src/share/util/scream_time_stamp.hpp
    // See this file for definition of Julian Date.
    Real t_beg_month;
    // Current simulation Julian Date
    Real t_now;
    // Number of days in the current month, cast as a Real
    Real days_this_month;
  };  // marineOrganicsTimeState

  struct marineOrganicsData {
    marineOrganicsData() = default;
    marineOrganicsData(const int ncol_, const int nsectors_)
        : ncols(ncol_), nsectors(nsectors_) {
      init(ncols, nsectors, true);
    }

    void init(const int ncol, const int nsector, const bool allocate) {
      ncols    = ncol;
      nsectors = nsector;
      if(allocate) emiss_sectors = view_2d("morgAllSectors", nsectors, ncols);
    }  // marineOrganicsData init

    // Basic spatial dimensions of the data
    int ncols, nsectors;
    view_2d emiss_sectors;
  };  // marineOrganicsData

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------
  struct marineOrganicsInput {
    marineOrganicsInput() = default;
    marineOrganicsInput(const int ncols_, const int nsectors_) {
      init(ncols_, nsectors_);
    }

    void init(const int ncols_, const int nsectors_) {
      data.init(ncols_, nsectors_, true);
    }
    marineOrganicsData data;  // All marineOrganics fields
  };                          // marineOrganicsInput

  // The output is really just marineOrganicsData, but for clarity it might
  // help to see a marineOrganicsOutput along a marineOrganicsInput in functions
  // signatures
  using marineOrganicsOutput = marineOrganicsData;

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  static std::shared_ptr<AbstractRemapper> create_horiz_remapper(
      const std::shared_ptr<const AbstractGrid> &model_grid,
      const std::string &marineOrganics_data_file, const std::string &map_file,
      const std::vector<std::string> &field_name, const std::string &dim_name1);

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  static std::shared_ptr<AtmosphereInput> create_data_reader(
      const std::shared_ptr<AbstractRemapper> &horiz_remapper,
      const std::string &data_file);

  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------
#if 0
  static void update_marine_organics_data_from_file(
      std::shared_ptr<AtmosphereInput> &scorpio_reader,
      AbstractRemapper &horiz_interp, const_view_1d &input);
#endif
  // -------------------------------------------------------------------------------------------
  // -------------------------------------------------------------------------------------------

  static void init_marine_organics_file_read(
      const int ncol, const std::vector<std::string> field_name,
      const std::string dim_name1,
      const std::shared_ptr<const AbstractGrid> &grid,
      const std::string &data_file, const std::string &mapping_file,
      // output
      std::shared_ptr<AbstractRemapper> &marineOrganicsHorizInterp,
      marineOrganicsInput morg_data_start_, marineOrganicsInput morg_data_end_,
      marineOrganicsData morg_data_out_,
      std::shared_ptr<AtmosphereInput> &marineOrganicsDataReader);

};  // struct marineOrganicsFunctions

}  // namespace marine_organics
}  // namespace scream
#endif  // MARINE_ORGANICS_HPP

#include "marine_organics_impl.hpp"