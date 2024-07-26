#ifndef SRF_EMISSION_HPP
#define SRF_EMISSION_HPP

#include "share/util/scream_timing.hpp"

namespace scream::mam_coupling {
namespace {

template <typename ScalarType, typename DeviceType>
struct srfEmissFunctions {
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
    srfEmissData(const int ncol_) { init(ncol_, true); }

    void init(const int ncol_, const bool allocate) {
      ncols = ncol_;

      if(allocate) {
        /*AGR = view_1d("", ncols);
        RCO = view_1d("", ncols);
        SHP = view_1d("", ncols);
        SLV = view_1d("", ncols);
        TRA = view_1d("", ncols);
        WST = view_1d("", ncols);*/
        for(auto &component : emiss_components) {
          component = view_1d("", ncols);
        }
      }
    }  // srfEmissData init

    // Basic spatial dimensions of the data
    int ncols;
    std::array<view_1d, 6> emiss_components;

    /*view_1d AGR;
    view_1d RCO;
    view_1d SHP;
    view_1d SLV;
    view_1d TRA;
    view_1d WST;*/
  };  // srfEmissData

  struct srfEmissInput {
    srfEmissInput() = default;
    srfEmissInput(const int ncols_) { init(ncols_); }

    void init(const int ncols_) { data.init(ncols_, true); }
    srfEmissData data;  // All srfEmiss fields
  };                    // srfEmissInput

  // The output is really just srfEmissData, but for clarity it might
  // help to see a srfEmissOutput along a srfEmissInput in functions signatures
  using srfEmissOutput = srfEmissData;

  /* -------------------------------------------------------------------------------------------
   */
  // Surface emissions routines
  template <std::size_t FN>
  static std::shared_ptr<AbstractRemapper> create_horiz_remapper(
      const std::shared_ptr<const AbstractGrid> &model_grid,
      const std::string &srfEmiss_data_file,
      const std::array<std::string, FN> &field_names,
      const std::string &map_file);

  static std::shared_ptr<AtmosphereInput> create_srfEmiss_data_reader(
      const std::shared_ptr<AbstractRemapper> &horiz_remapper,
      const std::string &srfEmiss_data_file);

  static void srfEmiss_main(const srfEmissTimeState &time_state,
                            const srfEmissInput &data_beg,
                            const srfEmissInput &data_end,
                            const srfEmissInput &data_tmp,  // Temporary
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

};  // struct srfEmissFunctions
}  // namespace
}  // namespace scream::mam_coupling
#endif  // SRF_EMISSION_HPP

#include "srf_emission_impl.hpp"