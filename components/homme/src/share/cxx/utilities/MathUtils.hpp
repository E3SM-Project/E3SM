/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_MATH_UTILS_HPP
#define HOMMEXX_MATH_UTILS_HPP

#include <Kokkos_Core.hpp>
#include <cmath>
#include "Types.hpp"

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

template <typename FPType>
KOKKOS_INLINE_FUNCTION constexpr FPType square(const FPType& x) {
  return x*x;
}

// Set the least significant n bits of the fraction to 0. Use this for
// BFB-testing between GPU and host. Zero the n least significant bits of y in
//     y = pow(x, e),
// for example.
//   Table 8 of Appendix E of the Cuda programming guide gives the accuracy of
// various functions. However, we also have to account for carries in the
// calculations. Thus, even if an answer is accurate to 2 ULP, carries can make
// more bits differ. Thus, zeroing (or any sort of rounding) is not assured to
// work. However, you can reduce the probability of a diff to a small number by
// choosing n to be large. Full double precision range is 53 bits; choose
// something like n = 20 to get a carry-based diff every 1 in a million calls,
// still leaving 33 bits.
//   To remain sensitive to the full bits in the argument to the non-BFB
// function, such as x in pow above, pass replace = x. The n least significant
// bits of the fraction part of x is used to fill the 0s introduced in this
// function. Thus, the final result z of
//     y = pow(x, e);
//     z = zeroulpn(y, n, x)
// is sensitive to all bits of x.
KOKKOS_INLINE_FUNCTION
double zeroulpn (double a, const int n, double replace) {
  const auto f = static_cast<double>(1LL << 54);
  int e;
  a = std::frexp(a, &e);
  const int sign = a >= 0 ? 1 : -1;
  a = std::abs(a);
  long long int x = f*a;
  x >>= n;
  x <<= n;
  { // Replace the n least significant bits with those from the fraction part of
    // 'replace'.
    const long long int ones = 0xffffffffffffffff;
    auto mask = ones;
    mask >>= n;
    mask <<= n;
    mask = ones & (mask ^ ones);
    int er;
    replace = std::abs(std::frexp(replace, &er));
    long long int r = f*replace;
    x = x | (mask & r);
  }
  a = x/f;
  return std::ldexp(sign*a, e);
}

template <typename Pack>
KOKKOS_INLINE_FUNCTION
void zero_ulp_n (Pack& a, const int n, const Pack& replace) {
  for (int s = 0; s < VECTOR_SIZE; ++s)
    a[s] = zeroulpn(a[s], n, replace[s]);
}

template <typename Pack>
KOKKOS_INLINE_FUNCTION
Pack bfb_pow (const Pack& a, const Real& e) {
  Pack b;
  for (int s = 0; s < VECTOR_SIZE; ++s)
    b[s] = std::pow(a[s], e);
  if (OnGpu<ExecSpace>::value)
    zero_ulp_n(b, 20, a);
  return b;
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

template<typename T, int N, typename... Properties>
KOKKOS_INLINE_FUNCTION
void binary_search (ViewType<T[N],Properties...> array,
                    const T& pivot, int& k) {
  int lo,hi;

  // Initialize lower and upper bounds
  if (array(k)>pivot) {
    lo = 0;
    hi = k;
  } else {
    lo = k;
    hi = N-1;
  }

  // Binary search until hi=lo+1
  while (hi>lo+1) {
    k = (lo+hi)/2;
    if (array(k)>pivot) {
      hi = k;
    } else {
      lo = k;
    }
  }
  k = lo;
}

} // namespace Homme

#endif // HOMMEXX_MATH_UTILS_HPP
