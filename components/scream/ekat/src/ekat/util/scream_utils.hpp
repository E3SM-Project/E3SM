#ifndef SCREAM_UTILS_HPP
#define SCREAM_UTILS_HPP

#include "ekat/scream_types.hpp"
#include "ekat/scream_assert.hpp"

#include <cstdio>
#include <cstring>
#include <memory>
#include <map>
#include <string>
#include <typeinfo>

#ifndef KOKKOS_ENABLE_CUDA
# include <cmath>
# include <algorithm>
#endif

namespace scream {

/*
 * Utility function for handling floating point literals,
 * so that they match the scream precision. This is
 * especially useful for bfb tests agaisnt fortran,
 * to ensure that literals are not a source of round-off differences.
 */
template<typename T> KOKKOS_INLINE_FUNCTION
constexpr typename std::enable_if<std::is_arithmetic<T>::value,Real>::type
sp (const T val) {
  return Real(val);
}

namespace util {

template <typename Real> struct is_single_precision {};
template <> struct is_single_precision<float> { enum : bool { value = true }; };
template <> struct is_single_precision<double> { enum : bool { value = false }; };

template <typename Real> struct is_double_precision {};
template <> struct is_double_precision<float> { enum : bool { value = false }; };
template <> struct is_double_precision<double> { enum : bool { value = true }; };

bool eq(const std::string& a, const char* const b1, const char* const b2 = 0);

// A templated class to return a human-readable name for a type (defaults to type_info implementation)
template<typename T>
struct TypeName {
  static std::string name () {
    return typeid(T).name();
  }
};

#ifdef KOKKOS_ENABLE_CUDA
// Replacements for namespace std functions that don't run on the GPU.
template <typename T> KOKKOS_INLINE_FUNCTION
const T& min (const T& a, const T& b) { return a < b ? a : b; }
template <typename T> KOKKOS_INLINE_FUNCTION
const T& max (const T& a, const T& b) { return a > b ? a : b; }
KOKKOS_INLINE_FUNCTION bool isfinite (const Real& a) {
  return a == a && a != INFINITY && a != -INFINITY;
}

template <typename T> KOKKOS_INLINE_FUNCTION
const T* max_element (const T* const begin, const T* const end) {
  const T* me = begin;
  for (const T* it = begin + 1; it < end; ++it)
    if ( ! (*it < *me)) // use operator<
      me = it;
  return me;
}
KOKKOS_INLINE_FUNCTION
size_t strlen(const char* str)
{
  scream_kassert(str != NULL);
  const char *char_ptr;
  for (char_ptr = str; ; ++char_ptr)  {
    if (*char_ptr == '\0') return char_ptr - str;
  }
}
KOKKOS_INLINE_FUNCTION
void strcpy(char* dst, const char* src)
{
  scream_kassert(dst != NULL && src != NULL);
  while(*dst++ = *src++);
}
KOKKOS_INLINE_FUNCTION
int strcmp(const char* first, const char* second)
{
  while(*first && (*first == *second))
  {
    first++;
    second++;
  }
  return *(const unsigned char*)first - *(const unsigned char*)second;
}
#else
using std::min;
using std::max;
using std::isfinite;
using std::max_element;
using std::strlen;
using std::strcpy;
using std::strcmp;
#endif

KOKKOS_INLINE_FUNCTION
bool is_nan (const Real& a) {
#ifdef __CUDA_ARCH__
  return isnan(a);
#else
  return std::isnan(a);
#endif
}

template <typename Integer> KOKKOS_INLINE_FUNCTION
void set_min_max (const Integer& lim0, const Integer& lim1,
                  Integer& min, Integer& max) {
  min = util::min(lim0, lim1);
  max = util::max(lim0, lim1);
}

template <typename Integer, typename Integer1> KOKKOS_INLINE_FUNCTION
void set_min_max (const Integer& lim0, const Integer& lim1,
                  Integer& min, Integer& max, const Integer1& vector_size) {
  min = util::min(lim0, lim1) / vector_size;
  max = util::max(lim0, lim1) / vector_size;
}

template <typename Real> KOKKOS_INLINE_FUNCTION
Real reldif (const Real& a, const Real& b) {
  return std::abs(b - a)/std::abs(a);
}

struct TransposeDirection {
  enum Enum { c2f, f2c };
};

// Switch whether i (column index) or k (level index) is the fast
// index. TransposeDirection::c2f makes i faster; f2c makes k faster.
template <TransposeDirection::Enum direction, typename Scalar>
void transpose(const Scalar* sv, Scalar* dv, Int ni, Int nk) {
  for (Int k = 0; k < nk; ++k)
    for (Int i = 0; i < ni; ++i)
      if (direction == TransposeDirection::c2f)
        dv[ni*k + i] = sv[nk*i + k];
      else
        dv[nk*i + k] = sv[ni*k + i];
}

namespace check_overloads
{
// Note: the trick used here is taken from https://stackoverflow.com/a/50631844/1093346

// Equality operator
template <typename LHS, typename RHS>
struct EqualExists {
  template<typename T, typename U>
  static auto test(T &&t, U &&u) -> decltype(t == u, void(), std::true_type{});
  static auto test(...) -> std::false_type;
  using type = decltype(test(std::declval<LHS>(),std::declval<RHS>()));

  static constexpr bool value = type::value;
};

// Streaming operator
template<typename T>
struct StreamExists {
  template<typename U>
  static auto test(U*)
    -> decltype(
        std::declval<std::ostream&>() << std::declval<U>(),
        std::true_type());

  template<typename>
  static std::false_type test(...);

  using type = decltype(test<T>(0));

  static constexpr bool value = type::value;
};

} // namespace check_overloads

} // namespace util
} // namespace scream

#endif // SCREAM_UTILS_HPP
