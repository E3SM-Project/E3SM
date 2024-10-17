#ifndef COMMON_FILE_READ_FUNC_HPP
#define COMMON_FILE_READ_FUNC_HPP

namespace scream {
namespace fileRead {

template <typename ScalarType, typename DeviceType>
struct commonFileReadFunc {
  using Device = DeviceType;

  using KT         = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  struct timeState {
    timeState() = default;
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
  };  // timeState
};

}  // namespace fileRead
}  // namespace scream

#endif  // COMMON_FILE_READ_FUNC_HPP