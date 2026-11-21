#ifndef SCREAM_COSP_SCOPS_IMPL_HPP
#define SCREAM_COSP_SCOPS_IMPL_HPP

#include "physics/cosp/cosp_cxx_functions.hpp"

#include <climits>

namespace scream {
namespace cosp {

/*
 * Subcolumn generator (scops.F90, and the seeding in subsample_and_optics in
 * cosp_c2f.F90).
 */

// The RNG seed of one point, from its surface pressure and the min/max of the surface
// pressure over the points of this call (so it depends on the domain decomposition).
// With equal min and max the Fortran converts a NaN to integer, which gives INT_MIN on
// x86 before the +1; that is mimicked here.
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
typename Functions<S,D>::RngState
Functions<S,D>::scops_seed(const int npoints, const S psfc, const S psfc_min, const S psfc_max)
{
  if (npoints == 1) return static_cast<RngState>(psfc);
  const S x = (psfc - psfc_min)/(psfc_max - psfc_min)*100000;
  if (!(x >= S(INT_MIN) && x < S(INT_MAX))) return RngState(INT_MIN) + 1;
  return static_cast<RngState>(x) + 1;
}

// SCOPS for one point: frac_out(i,k) is 0 (clear), 1 (stratiform) or 2 (convective) for
// subcolumn i and level k. seed is this point's RNG stream, and seed_last that of the last
// point of the call: at the first level the Fortran sets the thresholds from ran(npoints)
// for every point. Those thresholds are always overwritten before being used, but the
// draws still advance the stream. threshold is scratch space of size ncol.
template <typename S, typename D>
template <typename VC, typename VV, typename VF, typename VT>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::scops(const int nlev, const int ncol, RngState seed, RngState seed_last,
                           const VC& cc, const VV& conv, int overlap, const VF& frac_out,
                           const VT& threshold)
{
  // Invalid overlap assumptions fall back to maximum-random
  if (overlap != 1 && overlap != 2 && overlap != 3) overlap = 3;

  // Position of each subcolumn in the gridbox; computed in single precision in the Fortran
  auto boxpos = [&] (const int i) -> S {
    return S((static_cast<float>(i+1) - 0.5f)/static_cast<float>(ncol));
  };

  for (int k = 0; k < nlev; ++k) {
    if (k == 0) {
      for (int i = 0; i < ncol; ++i) {
        if (overlap == 1) {
          threshold(i) = boxpos(i);
        } else {
          get_rng(seed);
          const S ran_last = get_rng(seed_last);
          threshold(i) = conv(0) + (1 - conv(0))*ran_last;
        }
      }
    }

    for (int i = 0; i < ncol; ++i) {
      // Flag for maximum overlapped convective cloud
      const S maxocc = boxpos(i) <= conv(k) ? 1 : 0;

      // Flag for maximum overlapped stratiform cloud, and the lower bound of the new
      // threshold
      S threshold_min = 0, maxosc = 0;
      if (overlap == 1) {
        threshold_min = conv(k);
        maxosc        = 1;
      } else if (overlap == 2) {
        threshold_min = conv(k);
        maxosc        = 0;
      } else {
        const S cc_min = k > 0 ? Kokkos::min(cc(k-1), cc(k)) : Kokkos::min(S(0), cc(k));
        threshold_min = Kokkos::max(conv(k), cc_min);
        maxosc = (threshold(i) < cc_min && threshold(i) > conv(k)) ? 1 : 0;
      }

      const S ran = get_rng(seed);
      threshold(i) = maxocc*boxpos(i) +
                     (1 - maxocc)*(maxosc*threshold(i) +
                                   (1 - maxosc)*(threshold_min + (1 - threshold_min)*ran));

      frac_out(i,k) = cc(k) > threshold(i) ? 1 : 0;
      if (threshold(i) <= conv(k) && conv(k) > 0) frac_out(i,k) = C::I_CVC;
    }
  }
}

} // namespace cosp
} // namespace scream

#endif // SCREAM_COSP_SCOPS_IMPL_HPP
