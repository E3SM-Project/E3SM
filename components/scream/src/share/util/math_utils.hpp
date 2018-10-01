#ifndef SCREAM_MATH_UTILS_HPP
#define SCREAM_MATH_UTILS_HPP

#include <Kokkos_Core.hpp>

#include "share/scream_types.hpp"

namespace scream {

namespace util {

// Replacements foro namespace std functions that don't run on the GPU.
#ifdef KOKKOS_ENABLE_CUDA
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
#else
using std::min;
using std::max;
using std::isfinite;
using std::max_element;
#endif // KOKKOS_ENABLE_CUDA

template <typename Real>
Real reldif (const Real& a, const Real& b) {
  return std::abs(b - a)/std::abs(a);
}

} // namespace util
} // namespace scream

#endif // SCREAM_MATH_UTILS_HPP
