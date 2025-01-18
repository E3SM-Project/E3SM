#ifndef MARINE_ORGANICS_HPP
#define MARINE_ORGANICS_HPP

// For AtmosphereInput
#include "share/io/scorpio_input.hpp"

namespace scream {
namespace marine_organics {

template <typename ScalarType, typename DeviceType>
struct marineOrganicsFunctions {
  using Device = DeviceType;

  using KT         = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;
  using view_2d    = typename KT::template view_2d<Real>;

  // -------------------------------------------------------------------------------------------
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
    marineOrganicsData(const int &ncol_, const int &nfields_)
        : ncols(ncol_), nsectors(nfields_) {
      init(ncols, nsectors, true);
    }

    void init(const int &ncol, const int &nsector, const bool allocate) {
      ncols    = ncol;
      nsectors = nsector;
      if(allocate) emiss_sectors = view_2d("morgAllSectors", nsectors, ncols);
    }  // marineOrganicsData init

    // Basic spatial dimensions of the data
    int ncols, nsectors;
    view_2d emiss_sectors;
  };  // marineOrganicsData

  // -------------------------------------------------------------------------------------------
  struct marineOrganicsInput {
    marineOrganicsInput() = default;
    marineOrganicsInput(const int &ncols_, const int &nfields_) {
      init(ncols_, nfields_);
    }

    void init(const int &ncols_, const int &nfields_) {
      data.init(ncols_, nfields_, true);
    }
    marineOrganicsData data;  // All marineOrganics fields
  };                          // marineOrganicsInput

  // The output is really just marineOrganicsData, but for clarity it might
  // help to see a marineOrganicsOutput along a marineOrganicsInput in functions
  // signatures
  using marineOrganicsOutput = marineOrganicsData;

  // -------------------------------------------------------------------------------------------
  static std::shared_ptr<AbstractRemapper> create_horiz_remapper(
      const std::shared_ptr<const AbstractGrid> &model_grid,
      const std::string &marineOrganics_data_file, const std::string &map_file,
      const std::vector<std::string> &field_name, const std::string &dim_name1);

  // -------------------------------------------------------------------------------------------
  static std::shared_ptr<AtmosphereInput> create_data_reader(
      const std::shared_ptr<AbstractRemapper> &horiz_remapper,
      const std::string &data_file);

  // -------------------------------------------------------------------------------------------
  static void update_marine_organics_data_from_file(
      std::shared_ptr<AtmosphereInput> &scorpio_reader,
      const util::TimeStamp &ts,
      const int &time_index,  // zero-based
      AbstractRemapper &horiz_interp,
      marineOrganicsInput &marineOrganics_input);

  // -------------------------------------------------------------------------------------------
  static void update_marine_organics_timestate(
      std::shared_ptr<AtmosphereInput> &scorpio_reader,
      const util::TimeStamp &ts, AbstractRemapper &horiz_interp,
      marineOrganicsTimeState &time_state, marineOrganicsInput &beg,
      marineOrganicsInput &end);

  // -------------------------------------------------------------------------------------------
  static void marineOrganics_main(const marineOrganicsTimeState &time_state,
                                  const marineOrganicsInput &data_beg,
                                  const marineOrganicsInput &data_end,
                                  const marineOrganicsOutput &data_out);

  // -------------------------------------------------------------------------------------------
  static void perform_time_interpolation(
      const marineOrganicsTimeState &time_state,
      const marineOrganicsInput &data_beg, const marineOrganicsInput &data_end,
      const marineOrganicsOutput &data_out);

  // -------------------------------------------------------------------------------------------
  // Performs convex interpolation of x0 and x1 at point t
  template <typename ScalarX, typename ScalarT>
  KOKKOS_INLINE_FUNCTION static ScalarX linear_interp(const ScalarX &x0,
                                                      const ScalarX &x1,
                                                      const ScalarT &t);

  // -------------------------------------------------------------------------------------------
  static void init_marine_organics_file_read(
      const int &ncol, const std::vector<std::string> &field_name,
      const std::string &dim_name1,
      const std::shared_ptr<const AbstractGrid> &grid,
      const std::string &data_file, const std::string &mapping_file,
      // output
      std::shared_ptr<AbstractRemapper> &marineOrganicsHorizInterp,
      marineOrganicsInput &morg_data_start_,
      marineOrganicsInput &morg_data_end_, marineOrganicsData &morg_data_out_,
      std::shared_ptr<AtmosphereInput> &marineOrganicsDataReader);

};  // struct marineOrganicsFunctions

}  // namespace marine_organics
}  // namespace scream
#endif  // MARINE_ORGANICS_HPP

#include "marine_organics_impl.hpp"