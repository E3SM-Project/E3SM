#ifndef EAMXX_SIMPLE_LINEAR_INTERP_HPP
#define EAMXX_SIMPLE_LINEAR_INTERP_HPP

#include <array>
#include <algorithm>
#include <cstddef>

namespace scream {

/*
 * Simple 1D linear interpolation from coordinate grid x1 to x2.
 *
 * Assumes x1 and x2 are monotonically increasing.
 *
 * Template parameters:
 *   N1 - number of points in the source grid
 *   N2 - number of points in the target grid
 *
 * Arguments:
 *   x1        - source coordinate array (size N1)
 *   x2        - target coordinate array (size N2)
 *   y1        - source values (size N1)
 *   y2        - interpolated output values (size N2)  [output]
 */
template <std::size_t N1, std::size_t N2>
void simple_linear_interp(
  // const Kokkos::View<const Real[N1], Kokkos::MemoryTraits<Kokkos::Unmanaged>>& x1,
  // const Kokkos::View<const Real[N2], Kokkos::MemoryTraits<Kokkos::Unmanaged>>& x2,
  // const Kokkos::View<const Real[N1], Kokkos::MemoryTraits<Kokkos::Unmanaged>>& y1,
  // const Kokkos::View<      Real[N2], Kokkos::MemoryTraits<Kokkos::Unmanaged>>& y2);
  const Kokkos::View<const Real[N1], Kokkos::HostSpace>& x1,
  const Kokkos::View<const Real[N1], Kokkos::HostSpace>& y1,
  const Kokkos::View<const Real[N2], Kokkos::HostSpace>& x2,
  const Kokkos::View<      Real[N2], Kokkos::HostSpace>& y2);

} // namespace scream

#endif