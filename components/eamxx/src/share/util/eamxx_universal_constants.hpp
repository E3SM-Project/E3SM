#ifndef SCREAM_UNIVERSAL_CONSTANTS_HPP
#define SCREAM_UNIVERSAL_CONSTANTS_HPP

#include <limits>
#include <type_traits>

namespace scream {

namespace constants {

constexpr int seconds_per_day       = 86400;
constexpr int days_per_nonleap_year = 365;

// Universal fill value for variables
// NOTE: for floating point numbers, use the SAME numerical value, so that
//       we don't need to be aware of the precision of variables when checking
//       against fill_value (e.g., when reading in double data that was saved
//       in single precision)
template<typename T>
constexpr T fill_value =
  std::is_integral_v<T> ? std::numeric_limits<int>::max() / 2
                        : std::is_floating_point_v<T> ? std::numeric_limits<float>::max() / static_cast<float>(1e5)
                                                      : std::numeric_limits<char>::max();

} // namespace constants

} // namespace scream

#endif // SCREAM_UNIVERSAL_CONSTANTS_HPP
