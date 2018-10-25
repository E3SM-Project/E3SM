/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_MATH_UTILS_HPP
#define HOMMEXX_MATH_UTILS_HPP

#include <Kokkos_Core.hpp>
#include <cmath>

namespace Homme
{

template <typename FPType>
KOKKOS_INLINE_FUNCTION constexpr FPType min(const FPType val_1,
                                            const FPType val_2) {
  return val_1 < val_2 ? val_1 : val_2;
}

template <typename FPType, typename... FPPack>
KOKKOS_INLINE_FUNCTION constexpr FPType min(const FPType val, FPPack... pack) {
  return val < min(pack...) ? val : min(pack...);
}

template <typename FPType>
KOKKOS_INLINE_FUNCTION constexpr FPType max(const FPType val_1,
                                            const FPType val_2) {
  return val_1 > val_2 ? val_1 : val_2;
}

template <typename FPType, typename... FPPack>
KOKKOS_INLINE_FUNCTION constexpr FPType max(const FPType val, FPPack... pack) {
  return val > max(pack...) ? val : max(pack...);
}

// Computes the greatest common denominator of a and b with Euclid's algorithm
KOKKOS_INLINE_FUNCTION constexpr int gcd(const int a, const int b) {
	return (a % b == 0) ? b : gcd(b, a % b);
}

static_assert(gcd(1, 6) == 1, "gcd is broken");
static_assert(gcd(25, 20) == 5, "gcd is broken");
static_assert(gcd(29, 20) == 1, "gcd is broken");
static_assert(gcd(24, 16) == 8, "gcd is broken");

template <typename... int_pack>
KOKKOS_INLINE_FUNCTION constexpr int gcd(const int a, const int b, int_pack... pack) {
	return gcd(gcd(a, b), pack...);
}

static_assert(gcd(16, 24, 28) == 4, "gcd is broken");
static_assert(gcd(24, 16, 28) == 4, "gcd is broken");

// Computes the least common multiple of a and b
// Divide b by gcd(a, b) before multiplying to prevent overflows
KOKKOS_INLINE_FUNCTION constexpr int lcm(const int a, const int b) {
	return a * (b / gcd(a, b));
}

static_assert(lcm(1, 6) == 6, "lcm is broken");
static_assert(lcm(25, 20) == 100, "lcm is broken");
static_assert(lcm(29, 20) == 29 * 20, "lcm is broken");
static_assert(lcm(24, 16) == 48, "lcm is broken");

template <typename... int_pack>
KOKKOS_INLINE_FUNCTION constexpr int lcm(const int a, const int b, int_pack... pack) {
	return lcm(a, lcm(b, pack...));
}

static_assert(lcm(16, 24, 28) == 336, "lcm is broken");
static_assert(lcm(24, 16, 28) == 336, "lcm is broken");

template <typename ViewType>
typename std::enable_if<
    !std::is_same<typename ViewType::non_const_value_type, Scalar>::value, Real>::type
frobenius_norm(const ViewType view, bool ignore_nans = false) {
  typename ViewType::pointer_type data = view.data();

  size_t length = view.size();

  // Note: use Kahan algorithm to increase accuracy
  Real norm = 0;
  Real c = 0;
  Real temp, y;
  for (size_t i = 0; i < length; ++i) {
    if (std::isnan(data[i]) && ignore_nans) {
      continue;
    }
    y = data[i] * data[i] - c;
    temp = norm + y;
    c = (temp - norm) - y;
    norm = temp;
  }

  return std::sqrt(norm);
}

template <typename ViewType>
typename std::enable_if<
    std::is_same<typename ViewType::non_const_value_type, Scalar>::value, Real>::type
frobenius_norm(const ViewType view, bool ignore_nans = false) {
  typename ViewType::pointer_type data = view.data();

  size_t length = view.size();

  // Note: use Kahan algorithm to increase accuracy
  Real norm = 0;
  Real c = 0;
  Real temp, y;
  for (size_t i = 0; i < length; ++i) {
    for (int v = 0; v < VECTOR_SIZE; ++v) {
      if (std::isnan(data[i][v]) && ignore_nans) {
        continue;
      }
      y = data[i][v] * data[i][v] - c;
      temp = norm + y;
      c = (temp - norm) - y;
      norm = temp;
    }
  }

  return std::sqrt(norm);
}

} // namespace Homme

#endif // HOMMEXX_MATH_UTILS_HPP
