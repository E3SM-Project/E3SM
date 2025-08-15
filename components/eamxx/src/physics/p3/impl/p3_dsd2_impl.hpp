#ifndef P3_DSD2_IMPL_HPP
#define P3_DSD2_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 dsd2 functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::
get_cloud_dsd2(
  const Spack& qc, Spack& nc, Spack& mu_c, const Spack& rho, Spack& nu,
  const view_dnu_table& dnu, Spack& lamc, Spack& cdist, Spack& cdist1, 
  const Smask& context)
{
  lamc.set(context   , 0);
  cdist.set(context  , 0);
  cdist1.set(context , 0);
  nu.set(context     , 0);
  mu_c.set(context   , 0);

  constexpr Scalar qsmall = C::QSMALL;
  const auto qc_gt_small = qc >= qsmall && context;

  if (qc_gt_small.any()) {
    constexpr Scalar nsmall = C::NSMALL;
    constexpr Scalar cons1  = C::CONS1;
    // set minimum nc to prevent floating point error
    {
      Spack mu_c_local;
      nc.set(qc_gt_small, max(nc, nsmall));
      mu_c_local = sp(0.0005714)*(nc * sp(1.e-6) * rho) + sp(0.2714);
      mu_c_local = 1/(mu_c_local * mu_c_local) - 1;
      mu_c_local = max(mu_c_local, 2);
      mu_c_local = min(mu_c_local, 15);

      mu_c.set(qc_gt_small, mu_c_local);
    }

    // interpolate for mass distribution spectral shape parameter (for SB warm processes)
    if (P3C::iparam == 1) {
      IntSmallPack dumi = IntSmallPack(mu_c) - 1;
      Spack dnu0, dnu1;
      ekat::index_and_shift<1>(dnu, dumi, dnu0, dnu1);
      nu.set(qc_gt_small, dnu0 + (dnu1 - dnu0) * (mu_c - Spack(dumi) - 1));
    }

    // calculate lamc
    lamc.set(qc_gt_small, cbrt(cons1 * nc * (mu_c + 3) * (mu_c + 2) * (mu_c + 1) / qc));

    // apply lambda limiters
    Spack lammin = (mu_c + 1)*sp(2.5e+4); // min: 40 micron mean diameter
    Spack lammax = (mu_c + 1)*sp(1.e+6);   // max:  1 micron mean diameter

    Smask lamc_lt_min = lamc < lammin && qc_gt_small;
    Smask lamc_gt_max = lamc > lammax && qc_gt_small;
    Smask min_or_max = lamc_lt_min || lamc_gt_max;
    lamc.set(lamc_lt_min, lammin);
    lamc.set(lamc_gt_max, lammax);

    nc.set(min_or_max, 6 * (lamc * lamc * lamc) * qc / (C::Pi * C::RHO_H2O * (mu_c + 3) * (mu_c + 2) * (mu_c + 1)));

    cdist.set(qc_gt_small, nc * (mu_c+1) / lamc);
    cdist1.set(qc_gt_small, nc / tgamma(mu_c + 1));
  }
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::
get_rain_dsd2 (
  const Spack& qr, Spack& nr, Spack& mu_r,
  Spack& lamr,
  const P3Runtime& runtime_options,
  const Smask& context)
{
  constexpr auto nsmall = C::NSMALL;
  constexpr auto qsmall = C::QSMALL;
  constexpr auto cons1  = C::CONS1;

  lamr.set(context  , 0);

  const auto qr_gt_small = qr >= qsmall && context;

  const Scalar mu_r_const = runtime_options.constant_mu_rain;

  if (qr_gt_small.any()) {
    // use lookup table to get mu
    // mu-lambda relationship is from Cao et al. (2008), eq. (7)

    // find spot in lookup table
    // (scaled N/q for lookup table parameter space)
    const auto nr_lim = max(nr, nsmall);

    // Apply constant mu_r:  Recall the switch to v4 tables means constant mu_r
    mu_r.set(qr_gt_small, mu_r_const);
    // recalculate slope based on mu_r
    const auto mass_to_d3_factor = cons1 * (mu_r + 3) * (mu_r + 2) * (mu_r + 1);
    lamr.set(qr_gt_small, cbrt(mass_to_d3_factor * nr_lim / qr));

    // check for slope
    const auto lammax = (mu_r+1.)*sp(1.e+5);
    //Below, 500 is inverse of max allowable number-weighted mean raindrop size=2mm
    //Since breakup is explicitly included, mean raindrop size can be relatively small
    const auto lammin = (mu_r+1.)*500; 

    // apply lambda limiters for rain
    const auto lt = qr_gt_small && (lamr < lammin);
    const auto gt = qr_gt_small && (lamr > lammax);
    const auto either = lt || gt;
    nr.set(qr_gt_small, nr_lim);
    if (either.any()) {
      lamr.set(lt, lammin);
      lamr.set(gt, lammax);
      ekat_masked_loop(either, s) {
        nr[s] = lamr[s]*lamr[s]*lamr[s] * qr[s] / mass_to_d3_factor[s];
      }
    }
  }
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::
get_cdistr_logn0r (
  const Spack& qr, const Spack& nr, const Spack& mu_r,
  const Spack& lamr, Spack& cdistr, Spack& logn0r,
  const Smask& context)
{
  constexpr auto qsmall = C::QSMALL;

  cdistr.set(context, 0);
  logn0r.set(context, 0);

  const auto qr_gt_small = qr >= qsmall && context;

  if (qr_gt_small.any()) {
    cdistr.set(qr_gt_small, nr/tgamma(mu_r + 1));
    // note: logn0r is calculated as log10(n0r)
    logn0r.set(qr_gt_small, log10(cdistr) + (mu_r + 1) * log10(lamr));
  }
}

} // namespace p3
} // namespace scream

#endif // P3_DSD2_IMPL_HPP
