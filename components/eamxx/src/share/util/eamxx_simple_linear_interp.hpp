#ifndef EAMXX_SIMPLE_LINEAR_INTERP_HPP
#define EAMXX_SIMPLE_LINEAR_INTERP_HPP

#include "share/core/eamxx_types.hpp"

#include <Kokkos_Core.hpp>
#include <ekat_assert.hpp>
#include <cstddef>

namespace scream {

/*
 * Simple 1D linear interpolation from coordinate grid x1 to x2.
 *
 * Assumes x1 and x2 are monotonically increasing.
 *
 *   N1 - number of points in the source grid
 *   N2 - number of points in the target grid
 *
 * Arguments:
 *   x1       - source coordinate array (size N1)
 *   y1       - source values (size N1)
 *   x2       - target coordinate array (size N2)
 *   y2       - interpolated output values (size N2)  [output]
 */

template<typename V1D>
KOKKOS_INLINE_FUNCTION
void simple_linear_interp(
  const V1D& x1,
  const V1D& y1,
  const V1D& x2,
        V1D& y2)
{
  const std::size_t N1 = x1.size();
  const std::size_t N2 = x2.size();

  // Runtime checks only valid on host
  KOKKOS_IF_ON_HOST((
    EKAT_REQUIRE_MSG(N1 == (int)y1.extent(0), "Variable sizes do not match");
    EKAT_REQUIRE_MSG(N2 == (int)y2.extent(0), "Variable sizes do not match");
    EKAT_REQUIRE_MSG(N1 >= 2, "Source grid must have at least 2 points");
    EKAT_REQUIRE_MSG(N2 >= 1, "Target grid must have at least 1 point");
  ))

  for (std::size_t k2 = 0; k2 < N2; ++k2) {
    const Real xq = x2(k2);

    std::size_t lo = 0, hi_search = N1;
    while (lo < hi_search) {
      const std::size_t mid = lo + (hi_search - lo) / 2;
      if (x1(mid) < xq) lo = mid + 1;
      else              hi_search = mid;
    }
    std::size_t hi = lo;
    hi = hi < 1      ? 1      : hi;
    hi = hi > N1 - 1 ? N1 - 1 : hi;
    lo = hi - 1;

    const Real x_lo = x1(lo), x_hi = x1(hi);
    const Real y_lo = y1(lo), y_hi = y1(hi);
    const Real dx = x_hi - x_lo;

    y2(k2) = (dx != 0.0)
      ? y_lo + (y_hi - y_lo) * (xq - x_lo) / dx
      : y_lo;
  }
}

} // namespace scream

#endif