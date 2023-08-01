#ifndef SCREAM_UNIVERSAL_CONSTANTS_HPP
#define SCREAM_UNIVERSAL_CONSTANTS_HPP

#include <limits>

namespace scream {

namespace constants {

constexpr int seconds_per_day       = 86400;
constexpr int days_per_nonleap_year = 365;

// Universal fill value for variables
constexpr float DEFAULT_FILL_VALUE = std::numeric_limits<float>::max() / 1e5;

} // namespace constants

} // namespace scream

#endif // SCREAM_UNIVERSAL_CONSTANTS_HPP
