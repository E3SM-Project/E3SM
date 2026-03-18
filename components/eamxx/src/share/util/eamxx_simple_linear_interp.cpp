#include "eamxx_simple_linear_interp.hpp" // for ETI only but harmless for GPU

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
 *   x1       - source coordinate array (size N1)
 *   x2       - target coordinate array (size N2)
 *   y1       - source values (size N1)
 *   y2       - interpolated output values (size N2)  [output]
 */
template <std::size_t N1, std::size_t N2>
KOKKOS_INLINE_FUNCTION
void simple_linear_interp(
  // const Kokkos::View<const Real[N1], Kokkos::MemoryTraits<Kokkos::Unmanaged>>& x1,
  // const Kokkos::View<const Real[N2], Kokkos::MemoryTraits<Kokkos::Unmanaged>>& x2,
  // const Kokkos::View<const Real[N1], Kokkos::MemoryTraits<Kokkos::Unmanaged>>& y1,
  // const Kokkos::View<      Real[N2], Kokkos::MemoryTraits<Kokkos::Unmanaged>>& y2)
  const Kokkos::View<const Real[N1], Kokkos::HostSpace>& x1,
  const Kokkos::View<const Real[N1], Kokkos::HostSpace>& y1,
  const Kokkos::View<const Real[N2], Kokkos::HostSpace>& x2,
  const Kokkos::View<      Real[N2], Kokkos::HostSpace>& y2)
{
  static_assert(N1 >= 2, "Source grid must have at least 2 points");
  static_assert(N2 >= 1, "Target grid must have at least 1 point");

  for (std::size_t k2 = 0; k2 < N2; ++k2) {
    const double xq = x2(k2);

    // Binary search for bracketing interval in x1.
    // Finds the first index where x1(idx) >= xq, then step back
    // one to get the left bracket (clamped to valid range).
    std::size_t lo = 0, hi_search = N1;
    while (lo < hi_search) {
      const std::size_t mid = lo + (hi_search - lo) / 2;
      if (x1(mid) < xq) lo = mid + 1;
      else              hi_search = mid;
    }
    // lo is now the first index with x1(lo) >= xq (i.e., lower_bound)
    std::size_t hi = lo;
    hi = hi < 1      ? 1      : hi;  // clamp lo
    hi = hi > N1 - 1 ? N1 - 1 : hi;  // clamp hi
    lo = hi - 1;

    const double x_lo = x1(lo), x_hi = x1(hi);
    const double y_lo = y1(lo), y_hi = y1(hi);

    const double dx = x_hi - x_lo;

    // If source points coincide (degenerate interval), hold left value
    y2(k2) = (dx != 0.0)
      ? y_lo + (y_hi - y_lo) * (xq - x_lo) / dx
      : y_lo;
  }
}

// /*
//  * Simple 1D linear interpolation from coordinate grid x1 to x2.
//  *
//  * Assumes x1 and x2 are monotonically increasing.
//  *
//  * Template parameters:
//  *   N1 - number of points in the source grid
//  *   N2 - number of points in the target grid
//  *
//  * Arguments:
//  *   x1       - source coordinate array (size N1)
//  *   x2       - target coordinate array (size N2)
//  *   y1       - source values (size N1)
//  *   y2       - interpolated output values (size N2)  [output]
//  */
// template <std::size_t N1, std::size_t N2>
// void simple_linear_interp(
//   std::array<const double, N1>& x1,
//   std::array<const double, N2>& x2,
//   std::array<const double, N1>& y1,
//   std::array<      double, N2>& y2)
// {
//   static_assert(N1 >= 2, "Source grid must have at least 2 points");
//   static_assert(N2 >= 1, "Target grid must have at least 1 point");

//   for (std::size_t k2 = 0; k2 < N2; ++k2) {
//     const double xq = x2[k2];

//     // Find the bracketing interval in x1 via binary search.
//     // lower_bound gives the first element >= xq; we want the
//     // left bracket, so step back one (clamped to valid range).
//     auto it = std::lower_bound(x1.begin(), x1.end(), xq);

//     std::size_t hi = std::distance(x1.begin(), it);
//     hi = std::clamp(hi, std::size_t{1}, N1 - 1); // clamp to [1, N1-1]
//     const std::size_t lo = hi - 1;

//     const double x_lo = x1[lo], x_hi = x1[hi];
//     const double y_lo = y1[lo], y_hi = y1[hi];

//     const double dx = x_hi - x_lo;

//     // If source points coincide (degenerate interval), hold left value
//     y2[k2] = (dx != 0.0)
//       ? y_lo + (y_hi - y_lo) * (xq - x_lo) / dx
//       : y_lo;

//   }
// }

} // namespace scream
