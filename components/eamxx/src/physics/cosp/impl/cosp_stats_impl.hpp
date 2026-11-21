#ifndef SCREAM_COSP_STATS_IMPL_HPP
#define SCREAM_COSP_STATS_IMPL_HPP

#include "physics/cosp/cosp_cxx_functions.hpp"

namespace scream {
namespace cosp {

// hist2D (cosp_stats.F90). Bins are closed on the left and open on the right.
template <typename S, typename D>
template <typename V1, typename V2, typename B1, typename B2, typename H>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::hist2d(const int npts, const V1& var1, const V2& var2,
                            const int nbin1, const B1& bin1, const int nbin2, const B2& bin2,
                            const H& hist)
{
  for (int i = 0; i < nbin1; ++i) {
    for (int k = 0; k < nbin2; ++k) {
      int count = 0;
      for (int n = 0; n < npts; ++n) {
        if (var1(n) >= bin1(i) && var1(n) < bin1(i+1) &&
            var2(n) >= bin2(k) && var2(n) < bin2(k+1)) {
          ++count;
        }
      }
      hist(i,k) = count;
    }
  }
}

} // namespace cosp
} // namespace scream

#endif // SCREAM_COSP_STATS_IMPL_HPP
