#ifndef P3_FUNCTIONS_DSD2_IMPL_HPP
#define P3_FUNCTIONS_DSD2_IMPL_HPP

#include <fstream>

#include "p3_functions.hpp"
#include "p3_constants.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 dsd2 functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
void Functions<S,D>::
get_cloud_dsd2(const Smask& qc_gt_small, const Spack& qc, Spack& nc, Spack& mu_c, const Spack& rho, Spack& nu,
               const view_1d<const Scalar>& dnu, Spack& lamc, Spack& cdist, Spack& cdist1, const Spack& lcldm)
{
  if (qc_gt_small.any()) {
    // set minimum nc to prevent floating point error
    nc   = pack::max(nc, C::NSMALL);
    mu_c = 0.0005714*(nc * 1.e-6 * rho) + 0.2714;
    mu_c = 1./(pack::pow(mu_c, 2)) - 1.;
    mu_c = pack::max(mu_c, 2.);
    mu_c = pack::min(mu_c, 15.);

    // interpolate for mass distribution spectral shape parameter (for SB warm processes)
    if (P3C::iparam == 1) {
      IntSmallPack dumi = IntSmallPack(mu_c) - 1;
      Spack dnu0, dnu1;
      pack::index_and_shift<1>(dnu, dumi, dnu0, dnu1);
      nu   = dnu0 + (dnu1 - dnu0) * (mu_c - Spack(dumi) - 1);
    }

    // calculate lamc
    lamc = pack::pow(C::CONS1 * nc * (mu_c+3.) * (mu_c + 2.) * (mu_c + 1.) / qc, C::THIRD);

    // apply lambda limiters
    Spack lammin = (mu_c + 1.)*2.5e+4; // min: 40 micron mean diameter
    Spack lammax = (mu_c + 1.)*1.e+6;   // max:  1 micron mean diameter

    Smask lamc_lt_min = lamc < lammin;
    Smask lamc_gt_max = lamc > lammax;
    Smask min_or_max = lamc_lt_min || lamc_gt_max;
    lamc.set(lamc_lt_min, lammin);
    lamc.set(lamc_gt_max, lammax);
    nc.set(min_or_max, 6. * pack::pow(lamc, 3) * qc / (C::Pi * C::RHOW * (mu_c + 3.) * (mu_c + 2.) * (mu_c + 1.)));

    cdist  = nc * (mu_c+1.) / lamc;
    cdist1 = nc * lcldm / pack::tgamma(mu_c + 1.);
  }

  lamc.set(!qc_gt_small, 0.);
  cdist.set(!qc_gt_small, 0.);
  cdist1.set(!qc_gt_small, 0.);
}

template <typename S, typename D>
void Functions<S,D>::
get_rain_dsd2 (
    const view_1d_table& mu_r_table,
    const Smask& qr_gt_small, const Spack& qr, Spack& nr, Spack& mu_r,
    Spack& rdumii, IntSmallPack& dumii, Spack& lamr,
    Spack& cdistr, Spack& logn0r)
{
  constexpr auto nsmall = Constants<Scalar>::NSMALL;
  constexpr auto thrd = Constants<Scalar>::THIRD;
  constexpr auto cons1 = Constants<Scalar>::CONS1;

  lamr = 0;
  cdistr = 0;
  logn0r = 0;

  // use lookup table to get mu
  // mu-lambda relationship is from Cao et al. (2008), eq. (7)

  // find spot in lookup table
  // (scaled N/q for lookup table parameter space)
  const auto nr_lim = max(nr, nsmall);
  Spack inv_dum(0);
  inv_dum.set(qr_gt_small,
              pow(qr / (cons1 * nr_lim * 6.0), thrd));

  mu_r = 0;
  {
    const auto m1 = qr_gt_small && (inv_dum < 282.e-6);
    mu_r.set(m1, 8.282);
  }
  {
    const auto m2 = qr_gt_small && (inv_dum >= 282.e-6) && (inv_dum < 502.e-6);
    if (m2.any()) {
      scream_masked_loop(m2, s) {
        // Linearly interpolate mu_r.
        Scalar rdumiis = (inv_dum[s] - 250.e-6)*0.5e6;
        rdumiis = util::max<Scalar>(rdumiis, 1.0);
        rdumiis = util::min<Scalar>(rdumiis, 150.0);
        rdumii[s] = rdumiis;
        Int dumiis = rdumiis;
        dumiis = util::min(dumiis, 149);
        dumii[s] = dumiis;
        const auto mu_r_im1 = mu_r_table(dumiis-1);
        mu_r[s] = mu_r_im1 + (mu_r_table(dumiis) - mu_r_im1) * (rdumiis - dumiis);
      }
    }
  }

  // recalculate slope based on mu_r
  lamr.set(qr_gt_small,
           pow(cons1 * nr_lim * (mu_r + 3) *
               (mu_r + 2) * (mu_r + 1)/qr,
               thrd));

  // check for slope
  const auto lammax = (mu_r+1.)*1.e+5;
  // set to small value since breakup is explicitly included (mean size 0.8 mm)
  const auto lammin = (mu_r+1.)*1250.0;
  // apply lambda limiters for rain
  const auto lt = qr_gt_small && (lamr < lammin);
  const auto gt = qr_gt_small && (lamr > lammax);
  const auto either = lt || gt;
  nr.set(qr_gt_small, nr_lim);
  if (either.any()) {
    lamr.set(lt, lammin);
    lamr.set(gt, lammax);
    scream_masked_loop(either, s) {
      nr[s] = std::exp(3*std::log(lamr[s]) + std::log(qr[s]) +
                       std::log(std::tgamma(mu_r[s] + 1)) - std::log(std::tgamma(mu_r[s] + 4)))
        / cons1;
    }
  }
}

} // namespace p3
} // namespace scream

#endif
