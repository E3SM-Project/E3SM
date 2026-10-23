#ifndef SCREAM_COSP_MISR_IMPL_HPP
#define SCREAM_COSP_MISR_IMPL_HPP

#include "physics/cosp/cosp_cxx_functions.hpp"

namespace scream {
namespace cosp {

/*
 * MISR_SUBCOLUMN (MISR_simulator.F90) for one point and one subcolumn. Levels go from
 * the top of the atmosphere down, as in the Fortran. Outputs the column optical depth
 * and the MISR cloud-top height. The dist_model_layertops diagnostic is not computed,
 * since EAMxx does not use it.
 */
template <typename S, typename D>
template <typename VT, typename VZ, typename VA>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::misr_subcolumn(const int nlev, const VT& dtau, const VZ& zfull,
                                    const VA& at, const int sunlit, S& tau_out, S& ztop)
{
  if (sunlit != 1) {
    tau_out = C::R_UNDEF;
    ztop    = C::R_UNDEF;
    return;
  }

  ztop = 0;

  // Column optical depth
  tau_out = 0;
  for (int k = 0; k < nlev; ++k) tau_out += dtau(k);

  int thres_crossed = 0;
  S cloud_dtau = 0;
  S tau_cum = 0; // sum(dtau(1:ilev)), accumulated in the same order as the Fortran
  for (int k = 0; k < nlev; ++k) {
    const S dt = dtau(k);
    tau_cum += dt;

    // If there is a cloud, start the counter and store this height
    if (thres_crossed == 0 && dt > 0) {
      thres_crossed = 1;
      cloud_dtau    = 0;
    }

    if (thres_crossed < 99 && thres_crossed > 0) {
      if (dt == 0) {
        // End of the current cloud layer without a CTH boundary: restart the counter
        cloud_dtau = 0;
      } else {
        cloud_dtau = cloud_dtau + dt;
      }

      // If the cloud is continuous but optically thin (< 1) from above the current
      // layer's cloud top down to this level, MISR will see a top below the layer top
      if (dt > 0 && (cloud_dtau - dt) < 1) {
        if (dt < 1 || k == 0 || k == nlev-1) {
          // MISR likely penetrates to the middle of this layer
          ztop = zfull(k);
        } else {
          // Take the OD = 1 level into this layer
          ztop = S(0.5)*(zfull(k) + zfull(k-1)) - S(0.5)*(zfull(k-1) - zfull(k+1))/dt;
        }
      }

      // A distinctive water layer: take this as the CTH level
      if (dt > 1 && at(k) > 273) {
        thres_crossed = 99;
      }

      // If the column is optically thick down to here, MISR can't see anything below
      if (tau_cum > 5) {
        thres_crossed = 99;
      }
    }
  }

  // A cloud was found but no MISR cloud top was set
  if (thres_crossed == 1) {
    if (tau_out > S(0.5f)) {
      // keep the MISR detected CTH
    } else if (tau_out > S(0.2f)) {
      // MISR may detect it, but won't have a good height
      ztop = -1;
    } else {
      // MISR not likely to detect it at all
      ztop = 0;
    }
  }
}

/*
 * MISR_COLUMN (MISR_simulator.F90) for one point: the tau/CTH joint histogram, in
 * percent, indexed (itau, ihgt). The column cloud area and mean cloud-top height are
 * not computed, since EAMxx does not use them.
 */
template <typename S, typename D>
template <typename VZ, typename VT, typename H>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::misr_column(const int ncol, const VZ& ztop, const int sunlit,
                                 const VT& tau, const H& fq)
{
  constexpr int ntau = C::numMISRTauBins;
  constexpr int nhgt = C::numMISRHgtBins;

  if (sunlit != 1) {
    for (int i = 0; i < ntau; ++i)
      for (int k = 0; k < nhgt; ++k)
        fq(i,k) = C::R_UNDEF;
    return;
  }

  // Bin edges: optical depth, and height converted from km to m as in the Fortran
  S tau_edges[ntau+1], hgt_edges[nhgt+1];
  for (int i = 0; i <= ntau; ++i) tau_edges[i] = C::tau_binBounds(i);
  for (int k = 0; k <= nhgt; ++k) hgt_edges[k] = 1000*C::hgt_binBounds(k);

  // Histogram over the subcolumns, with optically thin (tau_wrk) and cloud-free
  // (ztop_wrk) subcolumns moved out of range, as the Fortran does with -999
  for (int i = 0; i < ntau; ++i) {
    for (int k = 0; k < nhgt; ++k) {
      int count = 0;
      for (int n = 0; n < ncol; ++n) {
        const S tau_wrk  = tau(n) > C::misr_tauchk ? tau(n) : S(-999);
        const S ztop_wrk = ztop(n) == 0 ? S(-999) : ztop(n);
        if (tau_wrk  >= tau_edges[i] && tau_wrk  < tau_edges[i+1] &&
            ztop_wrk >= hgt_edges[k] && ztop_wrk < hgt_edges[k+1]) {
          ++count;
        }
      }
      fq(i,k) = S(100)*S(count)/ncol;
    }
  }
}

} // namespace cosp
} // namespace scream

#endif // SCREAM_COSP_MISR_IMPL_HPP
