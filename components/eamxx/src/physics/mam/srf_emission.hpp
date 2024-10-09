#ifndef SRF_EMISSION_HPP
#define SRF_EMISSION_HPP

#include "share/util/scream_timing.hpp"

namespace scream::mam_coupling {
template <typename ScalarType, typename DeviceType>
struct srfEmissFunctions {
  using Device = DeviceType;

  using KT         = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  struct srfEmissTimeState {
    srfEmissTimeState() = default;
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
  };  // srfEmissTimeState

  struct srfEmissData {
    srfEmissData() = default;
    srfEmissData(const int ncol_, const int nsectors_)
        : ncols(ncol_), nsectors(nsectors_) {
      init(ncols, nsectors, true);
    }

    void init(const int ncol, const int nsector, const bool allocate) {
      ncols    = ncol;
      nsectors = nsector;
      if(allocate) emiss_sectors = view_2d("AllSectors", nsectors, ncols);
    }  // srfEmissData init

    // Basic spatial dimensions of the data
    int ncols, nsectors;
    view_2d emiss_sectors;
  };  // srfEmissData

  struct srfEmissInput {
    srfEmissInput() = default;
    srfEmissInput(const int ncols_, const int nsectors_) {
      init(ncols_, nsectors_);
    }

    void init(const int ncols_, const int nsectors_) {
      data.init(ncols_, nsectors_, true);
    }
    srfEmissData data;  // All srfEmiss fields
  };                    // srfEmissInput

  // The output is really just srfEmissData, but for clarity it might
  // help to see a srfEmissOutput along a srfEmissInput in functions signatures
  using srfEmissOutput = srfEmissData;

  /* -------------------------------------------------------------------------------------------
   */
  // Surface emissions routines
  template <std::size_t numSectors>
  static std::shared_ptr<AbstractRemapper> create_horiz_remapper(
      const std::shared_ptr<const AbstractGrid> &model_grid,
      const std::string &srfEmiss_data_file,
      const std::array<std::string, numSectors> &field_names,
      const std::string &map_file);

  static std::shared_ptr<AbstractRemapper> create_horiz_remapper(
      const std::shared_ptr<const AbstractGrid> &model_grid,
      const std::string &srfEmiss_data_file,
      const std::vector<std::string> &field_names, const std::string &map_file);

  static std::shared_ptr<AtmosphereInput> create_srfEmiss_data_reader(
      const std::shared_ptr<AbstractRemapper> &horiz_remapper,
      const std::string &srfEmiss_data_file);

  static void srfEmiss_main(const srfEmissTimeState &time_state,
                            const srfEmissInput &data_beg,
                            const srfEmissInput &data_end,
                            const srfEmissOutput &data_out);

  static void update_srfEmiss_data_from_file(
      std::shared_ptr<AtmosphereInput> &scorpio_reader,
      const util::TimeStamp &ts,
      const int time_index,  // zero-based
      AbstractRemapper &srfEmiss_horiz_interp, srfEmissInput &srfEmiss_input);
  static void update_srfEmiss_timestate(
      std::shared_ptr<AtmosphereInput> &scorpio_reader,
      const util::TimeStamp &ts, AbstractRemapper &srfEmiss_horiz_interp,
      srfEmissTimeState &time_state, srfEmissInput &srfEmiss_beg,
      srfEmissInput &srfEmiss_end);

  // The following three are called during srfEmiss_main
  static void perform_time_interpolation(const srfEmissTimeState &time_state,
                                         const srfEmissInput &data_beg,
                                         const srfEmissInput &data_end,
                                         const srfEmissOutput &data_out);

  // Performs convex interpolation of x0 and x1 at point t
  template <typename ScalarX, typename ScalarT>
  KOKKOS_INLINE_FUNCTION static ScalarX linear_interp(const ScalarX &x0,
                                                      const ScalarX &x1,
                                                      const ScalarT &t);
  template <std::size_t numSectors>
  static void init_srf_emiss_objects(
      const int ncol, const std::shared_ptr<const AbstractGrid> &grid,
      const std::string &data_file,
      const std::array<std::string, numSectors> &sectors,
      const std::string &srf_map_file,
      // output
      std::shared_ptr<AbstractRemapper> &SrfEmissHorizInterp,
      srfEmissInput &SrfEmissData_start, srfEmissInput &SrfEmissData_end,
      srfEmissOutput &SrfEmissData_out,
      std::shared_ptr<AtmosphereInput> &SrfEmissDataReader);

  static void init_srf_emiss_objects(
      const int ncol, const std::shared_ptr<const AbstractGrid> &grid,
      const std::string &data_file, const std::vector<std::string> &sectors,
      const std::string &srf_map_file,
      // output
      std::shared_ptr<AbstractRemapper> &SrfEmissHorizInterp,
      srfEmissInput &SrfEmissData_start, srfEmissInput &SrfEmissData_end,
      srfEmissOutput &SrfEmissData_out,
      std::shared_ptr<AtmosphereInput> &SrfEmissDataReader);

};  // struct srfEmissFunctions
}  // namespace scream::mam_coupling
#endif  // SRF_EMISSION_HPP

#include "srf_emission_impl.hpp"
