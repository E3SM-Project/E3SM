#ifndef SCREAM_STRING_UTILS_HPP
#define SCREAM_STRING_UTILS_HPP

#include <string>
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

} // namespace util
} // namespace scream

#endif // SCREAM_STRING_UTILS_HPP
