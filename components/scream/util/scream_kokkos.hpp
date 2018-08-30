#ifndef INCLUDE_SCREAM_KOKKOS
#define INCLUDE_SCREAM_KOKKOS

#include <Kokkos_Core.hpp>

namespace scream {
// Replacements for namespace std functions that don't run on the GPU.
namespace ko {
#ifdef KOKKOS_ENABLE_CUDA
template <typename T>
KOKKOS_INLINE_FUNCTION const T &min(const T &a, const T &b) {
  return a < b ? a : b;
}
template <typename T>
KOKKOS_INLINE_FUNCTION const T &max(const T &a, const T &b) {
  return a > b ? a : b;
}
KOKKOS_INLINE_FUNCTION bool isfinite(const Real &a) {
  return a == a && a != INFINITY && a != -INFINITY;
}
template <typename T>
KOKKOS_INLINE_FUNCTION const T *max_element(const T *const begin,
                                            const T *const end) {
  const T *me = begin;
  for(const T *it = begin + 1; it < end; ++it)
    if(!(*it < *me))  // use operator<
      me = it;
  return me;
}
#else
using std::isfinite;
using std::max;
using std::max_element;
using std::min;
#endif
}  // namespace ko
}  // namespace scream

#endif
