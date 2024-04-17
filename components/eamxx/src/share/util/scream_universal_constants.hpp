#ifndef SCREAM_UNIVERSAL_CONSTANTS_HPP
#define SCREAM_UNIVERSAL_CONSTANTS_HPP

#include <limits>
#include <type_traits>

namespace scream {

namespace constants {

constexpr int seconds_per_day       = 86400;
constexpr int days_per_nonleap_year = 365;

// Universal fill value for variables
// TODO: When we switch to supporting C++17 we can use a simple `inline constexpr` rather than a struct
template<typename T>
struct DefaultFillValue {
  static constexpr bool is_float = std::is_floating_point<T>::value;
  static constexpr bool is_int   = std::is_integral<T>::value;
  static constexpr T value = is_int ? std::numeric_limits<int>::max() / 2 :
	  is_float ? std::numeric_limits<float>::max() / 1e5 : std::numeric_limits<char>::max();

};

} // namespace constants

} // namespace scream

#endif // SCREAM_UNIVERSAL_CONSTANTS_HPP
