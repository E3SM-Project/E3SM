#ifndef SCREAM_STRING_UTILS_HPP
#define SCREAM_STRING_UTILS_HPP

#include <string>
#include <sstream>
#include <algorithm>

namespace scream {
namespace util {

inline std::string upper_case (const std::string& s) {
  std::string s_up = s;
  std::transform(s_up.begin(), s_up.end(), s_up.begin(),
                 [](unsigned char c)->char { return std::toupper(c); }
                );
  return s_up;
}

inline std::string lower_case (const std::string& s) {
  std::string s_lo = s;
  std::transform(s_lo.begin(), s_lo.end(), s_lo.begin(),
                 [](unsigned char c)->char { return std::tolower(c); }
                );
  return s_lo;
}

inline std::string strint (const std::string& s, const int i) {
  std::stringstream ss;
  ss << s << " " << i;
  return ss.str();
}

} // namespace util
} // namespace scream

#endif // SCREAM_STRING_UTILS_HPP
