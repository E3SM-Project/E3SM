#ifndef SCREAM_UTILS_HPP
#define SCREAM_UTILS_HPP

#include <string>
#include <typeinfo>

namespace scream {
namespace util {

template <typename Real> struct is_single_precision {};
template <> struct is_single_precision<float> { enum : bool { value = true }; };
template <> struct is_single_precision<double> { enum : bool { value = false }; };

bool eq(const std::string& a, const char* const b1, const char* const b2 = 0);

void activate_floating_point_exceptions_if_enabled();

// Use only when the situation demands it.
void deactivate_floating_point_exceptions_if_enabled();

std::string active_avx_string();
std::string config_string();

// A templated class to return a human-readable name for a type (defaults to type_info implementation)
template<typename T>
struct TypeName {
  static std::string name () {
    return typeid(T).name();
  }
};

} // namespace util
} // namespace scream

#endif // SCREAM_UTILS_HPP
