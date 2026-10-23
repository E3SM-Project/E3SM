#ifndef SCREAM_COSP_MODIS_IMPL_HPP
#define SCREAM_COSP_MODIS_IMPL_HPP

#include "physics/cosp/cosp_cxx_functions.hpp"

namespace scream {
namespace cosp {

/*
 * MODIS simulator (modis_simulator.F90). Only the pieces needed for the tau/CTP joint
 * histogram are ported. Levels go from the top down.
 *
 * Note that xmu = 0.866 and minConservativeW0 = 0.9999999 are real(wp) parameters
 * initialized from single precision literals in the Fortran.
 */

namespace impl {

// Shared first part of two_stream and two_stream_reflectance (delta-Eddington, Joseph
// et al. 1976). The exponentials are only computed for non-conservative scattering.
template <typename S>
struct TwoStreamCoefs {
  S tau, gamma1, gamma2, rk, beta, ef1, ef2;
  bool conservative;
};

template <typename S>
KOKKOS_INLINE_FUNCTION
TwoStreamCoefs<S> two_stream_coefs(const S tauint, const S gint, const S w0int)
{
  constexpr S xmu = S(0.866f);
  constexpr S minConservativeW0 = S(0.9999999f);

  TwoStreamCoefs<S> c;
  const S f  = gint*gint;
  c.tau      = (1 - w0int*f)*tauint;
  const S w0 = (1 - f)*w0int/(1 - w0int*f);
  const S g  = (gint - f)/(1 - f);

  // gamma3 and gamma4 are only needed for beam = 1
  c.gamma1 = (7 - w0*(4 + 3*g))/4;
  c.gamma2 = -(1 - w0*(4 - 3*g))/4;

  c.conservative = w0int > minConservativeW0;
  if (!c.conservative) {
    c.rk = Kokkos::sqrt(c.gamma1*c.gamma1 - c.gamma2*c.gamma2);
    const S rkx = c.rk*xmu;
    const S r4 = (1 - rkx*rkx)*(c.rk + c.gamma1);
    const S r5 = (1 - rkx*rkx)*(c.rk - c.gamma1);
    c.beta = -r5/r4;
    const S e1 = Kokkos::min(c.rk*c.tau, S(500));
    c.ef1 = Kokkos::exp(-e1);
    c.ef2 = Kokkos::exp(-2*e1);
  }
  return c;
}

} // namespace impl

// two_stream_reflectance (beam = 2)
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::two_stream_reflectance(const S tauint, const S gint, const S w0int)
{
  const auto c = impl::two_stream_coefs(tauint, gint, w0int);
  if (c.conservative) {
    return c.gamma1*c.tau/(1 + c.gamma1*c.tau);
  }
  return (c.gamma2*(1 - c.ef2))/((c.rk + c.gamma1)*(1 - c.beta*c.ef2));
}

// two_stream (beam = 2): reflectance and transmittance of one layer
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::two_stream(const S tauint, const S gint, const S w0int, S& ref, S& tra)
{
  const auto c = impl::two_stream_coefs(tauint, gint, w0int);
  if (c.conservative) {
    ref = c.gamma1*c.tau/(1 + c.gamma1*c.tau);
    tra = 1 - ref;
  } else {
    ref = (c.gamma2*(1 - c.ef2))/((c.rk + c.gamma1)*(1 - c.beta*c.ef2));
    tra = (2*c.rk*c.ef1)/((c.rk + c.gamma1)*(1 - c.beta*c.ef2));
  }
}

// compute_toa_reflectace + adding_doubling: TOA reflectance of the cloudy (tau > 0)
// layers, combined from the top down over a black surface
template <typename S, typename D>
template <typename VT, typename VG, typename VW>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::compute_toa_reflectance(const int nlev, const VT& tau, const VG& g,
                                          const VW& w0)
{
  bool first = true;
  S refl_cum = 0, tran_cum = 0;
  for (int k = 0; k < nlev; ++k) {
    if (!(tau(k) > 0)) continue;
    S refl, tran;
    two_stream(tau(k), g(k), w0(k), refl, tran);
    if (first) {
      refl_cum = refl;
      tran_cum = tran;
      first = false;
    } else {
      const S denom = 1 - refl_cum*refl;
      const S r = refl_cum + refl*(tran_cum*tran_cum)/denom;
      const S t = (tran_cum*tran)/denom;
      refl_cum = r;
      tran_cum = t;
    }
  }
  return refl_cum;
}

// cloud_top_pressure: extinction-weighted pressure down to tau_limit, trapezoidal in p.
// The Fortran passes (/ 0, tau(1:nlev) /) and phalf(1:nlev) to a routine declared with
// size nlev, so the bottom layer is never used; reproduced here.
template <typename S, typename D>
template <typename VT, typename VPH>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::modis_cloud_top_pressure(const int nlev, const VT& tau, const VPH& phalf,
                                           const S tau_limit)
{
  S total_tau = 0, total_product = 0;
  for (int k = 0; k < nlev-1; ++k) {
    const S dtau = tau(k), ptop = phalf(k), pbot = phalf(k+1);
    if (total_tau + dtau > tau_limit) {
      const S dx = tau_limit - total_tau;
      total_tau = total_tau + dx;
      // Trapezoidal rule when taking less than a full step
      total_product = total_product + ptop*dx + (pbot - ptop)*(dx*dx)/(2*dtau);
    } else {
      total_tau     = total_tau + dtau;
      total_product = total_product + dtau*(pbot + ptop)/2;
    }
    if (total_tau >= tau_limit) break;
  }
  return total_tau > 0 ? total_product/total_tau : S(0);
}

// weight_by_extinction: extinction-weighted value of f down to tau_limit
template <typename S, typename D>
template <typename VT, typename VF>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::modis_weight_by_extinction(const int nlev, const VT& tau, const VF& f,
                                             const S tau_limit)
{
  S total_tau = 0, total_product = 0;
  for (int k = 0; k < nlev; ++k) {
    if (total_tau + tau(k) > tau_limit) {
      const S dx = tau_limit - total_tau;
      total_tau     = total_tau + dx;
      total_product = total_product + dx*f(k);
    } else {
      total_tau     = total_tau + tau(k);
      total_product = total_product + tau(k)*f(k);
    }
    if (total_tau >= tau_limit) break;
  }
  return total_tau > 0 ? total_product/total_tau : S(0);
}

// The size retrieval: interpolate_to_min over the trial radii, with the predicted
// reflectances from two_stream_reflectance. trial_re, g and w0 for the trial radii are
// what cosp_modis_init precomputes.
template <typename S, typename D>
KOKKOS_INLINE_FUNCTION
S Functions<S,D>::modis_retrieve_size(const int phase, const S tau, const S obs_refl)
{
  constexpr int n = C::num_trial_res;
  const bool liq = phase == C::phaseIsLiquid || phase == C::phaseIsUndetermined;
  const S re_min = liq ? C::re_water_min : C::re_ice_min;
  const S re_max = liq ? C::re_water_max : C::re_ice_max;
  const int trial_phase = liq ? C::phaseIsLiquid : C::phaseIsIce;

  S x[n], diff[n];
  for (int i = 0; i < n; ++i) {
    x[i] = re_min + (re_max - re_min)/(n-1)*S(i);
    const S refl = two_stream_reflectance(tau, get_g_nir(trial_phase, x[i]),
                                          get_ssa_nir(trial_phase, x[i]));
    diff[i] = refl - obs_refl;
  }

  // interpolate_to_min: first location of the minimum of |diff|, then a bracketing pair
  int loc = 0;
  for (int i = 1; i < n; ++i) {
    if (Kokkos::abs(diff[i]) < Kokkos::abs(diff[loc])) loc = i;
  }
  int lb, ub;
  if (loc == 0) {
    lb = 0; ub = 1;
  } else if (loc == n-1) {
    lb = n-2; ub = n-1;
  } else if (diff[loc-1]*diff[loc] < 0) {
    lb = loc-1; ub = loc;
  } else {
    lb = loc; ub = loc+1;
  }
  S re;
  if (diff[lb]*diff[ub] < 0) {
    // Interpolate the root position linearly if we bracket the root
    re = x[ub] - diff[ub]*(x[ub] - x[lb])/(diff[ub] - diff[lb]);
  } else {
    re = C::re_fill;
  }
  return S(1.0e-06)*re;
}

// modis_subcolumn for one subcolumn of a sunlit point. phalf is the grid-mean interface
// pressure (Pa) and isccp_ctp the ISCCP cloud-top pressure of this subcolumn (Pa).
template <typename S, typename D>
template <typename VPH, typename VT, typename VF, typename VG, typename VW>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::modis_subcolumn(const int nlev, const VPH& phalf, const VT& tau,
                                     const VF& liq_frac, const VG& g, const VW& w0,
                                     const S isccp_ctp, int& phase, S& ctp, S& tau_ret,
                                     S& size_ret)
{
  // Optical depth retrieval, and cloud detection
  S tau_sum = 0;
  for (int k = 0; k < nlev; ++k) tau_sum += tau(k);
  const bool cloudy = tau_sum >= C::min_OpticalThickness;

  if (!cloudy) {
    ctp      = C::R_UNDEF;
    phase    = C::phaseIsNone;
    size_ret = C::R_UNDEF;
    tau_ret  = C::R_UNDEF;
    return;
  }
  tau_ret = tau_sum;

  // Cloud-top pressure (CO2 slicing)
  ctp = modis_cloud_top_pressure(nlev, tau, phalf, C::CO2Slicing_TauLimit);

  // Phase
  const S liq = modis_weight_by_extinction(nlev, tau, liq_frac, C::phase_TauLimit);
  if (liq >= C::phaseDiscrimination_Threshold) {
    phase = C::phaseIsLiquid;
  } else if (liq <= 1 - C::phaseDiscrimination_Threshold) {
    phase = C::phaseIsIce;
  } else {
    phase = C::phaseIsUndetermined;
  }

  // Size, from the observed near-IR reflectance
  const S obs_refl = compute_toa_reflectance(nlev, tau, g, w0);
  size_ret = modis_retrieve_size(phase, tau_ret, obs_refl);
  if (size_ret < 0 && size_ret != C::R_UNDEF) {
    size_ret = S(1.0e-06)*C::re_fill;
  }

  // Low clouds: use the ISCCP cloud-top pressure, as MODIS reverts to a thermal method
  if (ctp > C::CO2Slicing_PressureLimit) {
    ctp = isccp_ctp;
  }
}

// modis_column for one sunlit point: the tau/CTP joint histogram, in percent, indexed
// (itau, ip) with pressure bins in the Fortran order (low to high pressure)
template <typename S, typename D>
template <typename VPh, typename VC, typename VT, typename VS, typename H>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::modis_column(const int ncol, const VPh& phase, const VC& ctp,
                                  const VT& tau, const VS& size, const H& fq)
{
  constexpr int ntau  = C::numMODISTauBins;
  constexpr int npres = C::numMODISPresBins;
  for (int i = 0; i < ntau; ++i) {
    for (int k = 0; k < npres; ++k) {
      int count = 0;
      for (int n = 0; n < ncol; ++n) {
        // Only successful retrievals are counted
        const bool cloudy = phase(n) != C::phaseIsNone && size(n) > 0;
        const S tau_wrk = cloudy ? tau(n) : S(-999);
        const S ctp_wrk = cloudy ? ctp(n) : S(-999);
        if (tau_wrk >= C::tau_binBounds(i) && tau_wrk < C::tau_binBounds(i+1) &&
            ctp_wrk >= 100*C::pres_binBounds(k) && ctp_wrk < 100*C::pres_binBounds(k+1)) {
          ++count;
        }
      }
      fq(i,k) = (S(count)/ncol)*S(100);
    }
  }
}

} // namespace cosp
} // namespace scream

#endif // SCREAM_COSP_MODIS_IMPL_HPP
