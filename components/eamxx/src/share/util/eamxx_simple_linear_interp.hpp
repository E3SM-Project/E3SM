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
 *
 * Host-side dispatcher: deduces the execution space from the output view and
 * launches a parallel_for over the (independent) target points. All four
 * views must share a memory space.
 */

template<typename Vx1, typename Vy1, typename Vx2, typename Vy2>
void simple_linear_interp(
  const Vx1& x1,
  const Vy1& y1,
  const Vx2& x2,
  const Vy2& y2)
{
  using exe_space = typename Vy2::execution_space;
  const std::size_t N1 = x1.extent(0);
  const std::size_t N2 = x2.extent(0);

  EKAT_REQUIRE_MSG(N1 == y1.extent(0), "Variable sizes do not match");
  EKAT_REQUIRE_MSG(N2 == y2.extent(0), "Variable sizes do not match");
  EKAT_REQUIRE_MSG(N1 >= 2, "Source grid must have at least 2 points");
  EKAT_REQUIRE_MSG(N2 >= 1, "Target grid must have at least 1 point");

  Kokkos::parallel_for("simple_linear_interp",
    Kokkos::RangePolicy<exe_space>(0, N2),
    KOKKOS_LAMBDA(const std::size_t k2) {
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
    });
}

} // namespace scream

#endif
