#ifndef ONLINE_EMISSION_HPP
#define ONLINE_EMISSION_HPP

#include "share/util/scream_timing.hpp"

namespace scream::mam_coupling {
namespace {

template <typename ScalarType, typename DeviceType> struct onlineEmissions {
  using Device = DeviceType;

  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  struct onlineEmissTimeState {
    onlineEmissTimeState() = default;
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
  }; // onlineEmissTimeState

  struct onlineEmissData {
    // Basic spatial dimensions of the data
    int ncols;
    view_2d flux_data;
    // local copy of main fluxes array
    view_2d constituent_fluxes;
    // FIXME: read this from input or get from mam4xx::aero_model_emissions?
    const std::vector<std::string> spec_names = {"so2",    "soag",   "bc_a4",
                                                 "num_a1", "num_a2", "num_a4",
                                                 "pom_a4", "so4_a1", "so4_a2"};
    // FIXME: change this when the above is dynamically-determined
    int nspec = 9;
    const std::string root_IC_str = "initial_condition_";

    onlineEmissData() = default;
    onlineEmissData(const int ncol_, const int nspec_)
        : ncols(ncol_), nspec(nspec_) {
      init(ncols, nspec, true);
    }

    onlineEmissData(const int ncol_) : ncols(ncol_) {
      init(ncols, nspec, true);
    }

    void init(const int ncol_, const int nspec_, const bool allocate_) {
      ncols = ncol_;
      nspec = nspec_;
      if (allocate_)
        flux_data = view_2d("onlineEmissData", nspec, ncols);
    } // onlineEmissData init
    void init(const int ncol_, const bool allocate_) {
      ncols = ncol_;
      if (allocate_)
        flux_data = view_2d("onlineEmissData", nspec, ncols);
    } // onlineEmissData init
    void init(const int ncol_) {
      ncols = ncol_;
      flux_data = view_2d("onlineEmissData", nspec, ncols);
    } // onlineEmissData init
  };  // onlineEmissData

  // The output is really just onlineEmissData, but for clarity it might
  // help to see a onlineEmissOutput along a onlineEmissInput in functions
  // signatures
  // using onlineEmissOutput = onlineEmissData;

  // ---------------------------------------------------------------------------
  // Online emissions routines
  // ---------------------------------------------------------------------------
  // FIXME:
  // static void onlineEmiss_main(const onlineEmissTimeState &time_state,
  //                           const onlineEmissInput &data_beg,
  //                           const onlineEmissInput &data_end,
  //                           const onlineEmissOutput &data_out);

  static void init_from_input_file(const ekat::ParameterList &m_params,
                                   onlineEmissData &data);
  static void
  transfer_to_cflux(const onlineEmissData &data,
                       const std::map<std::string, int> idx_map,
                       view_2d &fluxes);
  // static void update_onlineEmiss_timestate(
  //     std::shared_ptr<AtmosphereInput> &scorpio_reader,
  //     const util::TimeStamp &ts, AbstractRemapper &onlineEmiss_horiz_interp,
  //     onlineEmissTimeState &time_state, onlineEmissInput &onlineEmiss_beg,
  //     onlineEmissInput &onlineEmiss_end);

}; // struct onlineEmissions
} // namespace
} // namespace scream::mam_coupling
#endif // ONLINE_EMISSION_HPP

#include "online_emission_impl.hpp"
