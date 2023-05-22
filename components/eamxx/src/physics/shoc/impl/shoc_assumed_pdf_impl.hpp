#ifndef SHOC_SHOC_ASSUMED_PDF_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp"

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_assumed_pdf. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 *
 *  Purpose of this function is to calculate the
 *  double Gaussian PDF of SHOC, which is the centerpiece
 *  of the scheme.  The main outputs are the SGS cloud
 *  fraction and liquid water amount, in addition to the
 *  SGS buoyancy flux which is needed to close the SGS
 *  TKE equation.  This code follows the appendix of
 *  Larson et al. (2002) for Analytic Double Gaussian 1.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_assumed_pdf(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const uview_1d<const Spack>& thetal,
  const uview_1d<const Spack>& qw,
  const uview_1d<const Spack>& w_field,
  const uview_1d<const Spack>& thl_sec,
  const uview_1d<const Spack>& qw_sec,
  const uview_1d<const Spack>& wthl_sec,
  const uview_1d<const Spack>& w_sec,
  const uview_1d<const Spack>& wqw_sec,
  const uview_1d<const Spack>& qwthl_sec,
  const uview_1d<const Spack>& w3,
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const Workspace&             workspace,
  const uview_1d<Spack>&       shoc_cldfrac,
  const uview_1d<Spack>&       shoc_ql,
  const uview_1d<Spack>&       wqls,
  const uview_1d<Spack>&       wthv_sec,
  const uview_1d<Spack>&       shoc_ql2)
{
  // Define temporary variables
  uview_1d<Spack> wthl_sec_zt, wqw_sec_zt, w3_zt,
                  thl_sec_zt, qwthl_sec_zt, qw_sec_zt;
  workspace.template take_many_contiguous_unsafe<6>(
    {"wthl_sec_zt", "wqw_sec_zt", "w3_zt",
     "thl_sec_zt", "qwthl_sec_zt", "qw_sec_zt"},
    {&wthl_sec_zt, &wqw_sec_zt, &w3_zt,
     &thl_sec_zt, &qwthl_sec_zt, &qw_sec_zt});

  // Tolerances, thresholds, and constants
  const Scalar thl_tol = 1e-2;
  const Scalar rt_tol = 1e-4;
  const Scalar w_tol_sqd = 4e-4;
  const Scalar w_thresh = 0;
  const Scalar largeneg = SC::largeneg;
  const Scalar rair = C::Rair;
  const Scalar rv = C::RV;
  const bool dothetal_skew = SC::dothetal_skew;
  const Scalar basepres = C::P0;
  const Scalar cp = C::CP;
  const Scalar lcond = C::LatVap;
  const Scalar pi = C::Pi;
  const Scalar basetemp = C::basetemp;
  const Scalar epsterm = rair/rv;
  const Scalar Tl_min = 100;

  // Interpolate many variables from interface grid to thermo grid
  linear_interp(team,zi_grid,zt_grid,w3,w3_zt,nlevi,nlev,largeneg);
  linear_interp(team,zi_grid,zt_grid,thl_sec,thl_sec_zt,nlevi,nlev,0);
  linear_interp(team,zi_grid,zt_grid,wthl_sec,wthl_sec_zt,nlevi,nlev,largeneg);
  linear_interp(team,zi_grid,zt_grid,qwthl_sec,qwthl_sec_zt,nlevi,nlev,largeneg);
  linear_interp(team,zi_grid,zt_grid,wqw_sec,wqw_sec_zt,nlevi,nlev,largeneg);
  linear_interp(team,zi_grid,zt_grid,qw_sec,qw_sec_zt,nlevi,nlev,0);

  // The following is morally a const var, but there are issues with
  // gnu and std=c++14. The macro ConstExceptGnu is defined in ekat_kokkos_types.hpp.
  ConstExceptGnu Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    // Initialize cloud variables to zero
    shoc_cldfrac(k) = 0;
    if (k==0) { shoc_ql(k)[0] = 0; }
    shoc_ql2(k) = 0;

    const auto pval = pres(k);

    // Get all needed input moments for the PDF at this particular point
    const auto thl_first = thetal(k);
    const auto w_first = w_field(k);
    const auto qw_first = qw(k);
    const auto w3var = w3_zt(k);
    const auto thlsec = thl_sec_zt(k);
    const auto qwsec = qw_sec_zt(k);
    const auto qwthlsec = qwthl_sec_zt(k);
    const auto wqwsec = wqw_sec_zt(k);
    const auto wthlsec = wthl_sec_zt(k);
    const auto w2sec = w_sec(k);

    // Compute square roots of some variables so we don't have to compute these again
    const auto sqrtw2 = ekat::sqrt(w_sec(k));
    const auto sqrtthl = ekat::max(thl_tol,ekat::sqrt(thlsec));
    const auto sqrtqt = ekat::max(rt_tol,ekat::sqrt(qwsec));

    // Find parameters for vertical velocity
    Spack Skew_w(0), w1_1(w_first), w1_2(w_first), w2_1(0), w2_2(0), a(0.5);
    {
      const Smask condition = w_sec(k) > w_tol_sqd;

      const Scalar tmp_val(0.4);
      const Scalar one_m_tmp_val(1 - tmp_val);
      const Scalar sqrtw2t(std::sqrt(1-tmp_val));

      Skew_w.set(condition, w3var/ekat::sqrt(ekat::cube(w_sec(k))));
      a.set(condition,
            ekat::max(sp(0.01),
                      ekat::min(sp(0.99),
                                sp(0.5)*(1 - Skew_w*ekat::sqrt(1/(4*(one_m_tmp_val*one_m_tmp_val*one_m_tmp_val)
                                                                  + ekat::square(Skew_w)))))));

      w1_1.set(condition, ekat::sqrt((1 - a)/a)*sqrtw2t);
      w1_2.set(condition, -1*ekat::sqrt(a/(1 - a))*sqrtw2t);
      w2_1.set(condition, tmp_val*w_sec(k));
      w2_2.set(condition, tmp_val*w_sec(k));
    }

    // Find parameters for thetal
    Spack thl1_1(thl_first), thl1_2(thl_first), thl2_1(0), thl2_2(0),
          sqrtthl2_1(0), sqrtthl2_2(0);
    {
      const Smask condition =  thlsec > (thl_tol*thl_tol) && ekat::abs(w1_2 - w1_1) > w_thresh;

      const Spack corrtest1 = ekat::max(-1, ekat::min(1, wthlsec/(sqrtw2*sqrtthl)));
      const Spack tmp_val_1(-corrtest1/w1_1), tmp_val_2(-corrtest1/w1_2);

      Spack Skew_thl(0);
      if (dothetal_skew == true) {
        const auto tsign = ekat::abs(tmp_val_1 - tmp_val_2);
        Skew_thl.set(tsign>sp(0.4), sp(1.2)*Skew_w);
        Skew_thl.set(tsign>sp(0.2) && tsign<=sp(0.4), (((sp(1.2)*Skew_w)/sp(0.2))*(tsign-sp(0.2))));
      }

      if (condition.any()) {
        thl2_1.set(condition,
                   ekat::min(100,
                             ekat::max(0, (3*tmp_val_1*(1 - a*ekat::square(tmp_val_2) - (1-a)*ekat::square(tmp_val_1))
                                           - (Skew_thl - a*ekat::cube(tmp_val_2) - (1 - a)*ekat::cube(tmp_val_1)))
                                       /(3*a*(tmp_val_1 - tmp_val_2))))*thlsec);
        thl2_2.set(condition,
                   ekat::min(100,
                           ekat::max(0, (-3*tmp_val_2*(1 - a*ekat::square(tmp_val_2)
                                         - (1 - a)*ekat::square(tmp_val_1))
                                         + (Skew_thl - a*ekat::cube(tmp_val_2) - (1 - a)*ekat::cube(tmp_val_1)))
                                         /(3*(1 - a)*(tmp_val_1 - tmp_val_2))))*thlsec);

        thl1_1.set(condition, tmp_val_2*sqrtthl+thl_first);
        thl1_2.set(condition, tmp_val_1*sqrtthl+thl_first);

        sqrtthl2_1.set(condition, ekat::sqrt(thl2_1));
        sqrtthl2_2.set(condition, ekat::sqrt(thl2_2));
      }

    }

    // Find parameters for total water mixing ratio
    Spack qw1_1(qw_first), qw1_2(qw_first), qw2_1(0), qw2_2(0),
          sqrtqw2_1(0), sqrtqw2_2(0);
    {
      const Smask condition = qwsec > (rt_tol*rt_tol) && ekat::abs(w1_2 - w1_1) > w_thresh;

      const Spack corrtest2 = ekat::max(-1, ekat::min(1, wqwsec/(sqrtw2*sqrtqt)));
      const Spack tmp_val_1(-corrtest2/w1_1), tmp_val_2(-corrtest2/w1_2);

      const auto tsign = ekat::abs(tmp_val_1 - tmp_val_2);
      Spack Skew_qw(0);
      Skew_qw.set(tsign>sp(0.4), sp(1.2)*Skew_w);
      Skew_qw.set(tsign>sp(0.2) && tsign<=sp(0.4), (((sp(1.2)*Skew_w)/sp(0.2))*(tsign-sp(0.2))));

      if (condition.any()) {
        qw2_1.set(condition,
                ekat::min(100,
                          ekat::max(0, (3*tmp_val_1*(1 - a*ekat::square(tmp_val_2) - (1 - a)*ekat::square(tmp_val_1))
                                        - (Skew_qw - a*ekat::cube(tmp_val_2) - (1 - a)*ekat::cube(tmp_val_1)))
                                        /(3*a*(tmp_val_1 - tmp_val_2))))*qwsec);
        qw2_2.set(condition,
                ekat::min(100,
                          ekat::max(0, (-3*tmp_val_2*(1 - a*ekat::square(tmp_val_2) - (1 - a)*ekat::square(tmp_val_1))
                                        + (Skew_qw - a*ekat::cube(tmp_val_2) - (1 - a)*ekat::cube(tmp_val_1)))
                                        /(3*(1 - a)*(tmp_val_1 - tmp_val_2))))*qwsec);
      }

      qw1_1.set(condition, tmp_val_2*sqrtqt+qw_first);
      qw1_2.set(condition, tmp_val_1*sqrtqt+qw_first);

      sqrtqw2_1.set(condition, ekat::sqrt(qw2_1));
      sqrtqw2_2.set(condition, ekat::sqrt(qw2_2));
    }

    // Convert from tilde variables to "real" variables
     w1_1 *= sqrtw2;
     w1_1 += w_first;
     w1_2 *= sqrtw2;
     w1_2 += w_first;

     // Find within-plume correlations.
      Spack r_qwthl_1(0);
      {
        const Spack testvar = a*sqrtqw2_1*sqrtthl2_1 + (1 - a)*sqrtqw2_2*sqrtthl2_2;
        const auto testvar_ne_zero = testvar != 0;
        if (testvar_ne_zero.any()) {
          r_qwthl_1.set(testvar_ne_zero,
                        ekat::max(-1,
                                  ekat::min(1, (qwthlsec - a*(qw1_1 - qw_first)*(thl1_1 - thl_first)
                                                - (1 - a)*(qw1_2 - qw_first)*(thl1_2 - thl_first))/testvar)));
        }
      }

      // Begin to compute cloud property statistics
      Spack Tl1_1 = thl1_1/(ekat::pow(basepres/pval,(rair/cp)));
      Spack Tl1_2 = thl1_2/(ekat::pow(basepres/pval,(rair/cp)));

      const auto index_range = ekat::range<IntSmallPack>(k*Spack::n);
      const Smask active_entries = (index_range < nlev);

      // Do NaN Checking here
      const Smask is_nan_Tl1_1 = isnan(Tl1_1) && active_entries;
      const Smask is_nan_Tl1_2 = isnan(Tl1_2) && active_entries;
      if (is_nan_Tl1_1.any() || is_nan_Tl1_2.any()) {
        printf("WARNING: NaN Detected in Tl1_1 or Tl1_2!\n");
        for (int i=0; i<is_nan_Tl1_1.n; i++) {
          if (is_nan_Tl1_1[i] || is_nan_Tl1_2[i]) {
            printf(
              "Tl1 NaN Detected: lev, Tl1_1, Tl1_2: %d, %16.9e, %16.9e\n"
              "  thetal, qw, pressure, thl_sec, qw_sec, w2sec:"
              " %16.9e, %16.9e, %16.9e, %16.9e, %16.9e, %16.9e, %16.9e\n"
              "  w_first, w3, qwthlsec, wqwsec, wthlsec, a:"
              " %16.9e, %16.9e, %16.9e, %16.9e, %16.9e\n"
              "  w1_1, w1_2, w2_1, w2_2, thl1_1, thl1_2, thl2_1, thl2_2:"
              " %16.9e, %16.9e, %16.9e, %16.9e, %16.9e, %16.9e, %16.9e, %16.9e\n"
              "  qw1_1, qw1_2, qw2_1, qw2_2:"
              " %16.9e, %16.9e, %16.9e, %16.9e\n",
              index_range[i], Tl1_1[i], Tl1_2[i],
              thl_first[i], qw_first[i], pval[i], thlsec[i], qwsec[i], w2sec[i],
              w_first[i], w3var[i], qwthlsec[i], wqwsec[i], wthlsec[i], a[i],
              w1_1[i], w1_2[i], w2_1[i], w2_2[i], thl1_1[i], thl1_2[i], thl2_1[i], thl2_2[i],
              qw1_1[i], qw1_2[i], qw2_1[i], qw2_2[i]);
          }
        }
      }

      // Check to ensure Tl1_1 and Tl1_2 are not excessively small, set to threshold value if so.
      const Smask is_small_Tl1_1 = (Tl1_1 <= Tl_min) && active_entries;
      const Smask is_small_Tl1_2 = (Tl1_2 <= Tl_min) && active_entries;
      if( is_small_Tl1_1.any() ) {
        Tl1_1.set(is_small_Tl1_1,Tl_min);
        int n_mask = 0;
        for (int i=0; i<is_small_Tl1_1.n; i++) {
          if (is_small_Tl1_1[i]) {
            n_mask++;
          }
        }
        printf("WARNING: Tl1_1 has %d values <= allowable value.  Resetting to minimum value.\n",n_mask);
      }
      if( is_small_Tl1_2.any() ) {
        Tl1_2.set(is_small_Tl1_2,Tl_min);
        int n_mask = 0;
        for (int i=0; i<is_small_Tl1_2.n; i++) {
          if (is_small_Tl1_2[i]) {
            n_mask++;
          }
        }
        printf("WARNING: Tl1_2 has %d values <= allowable value.  Resetting to minimum value.\n",n_mask);
      }

      // Compute qs and beta
      Spack qs1(0), qs2(0), beta1(0), beta2(0);
      {
        // Compute MurphyKoop_svp
        const int liquid = 0;
        const Spack esval1_1 = scream::physics::Functions<S,D>::MurphyKoop_svp(Tl1_1,liquid,active_entries,"shoc::shoc_assumed_pdf (Tl1_1)");
        const Spack esval1_2 = scream::physics::Functions<S,D>::MurphyKoop_svp(Tl1_2,liquid,active_entries,"shoc::shoc_assumed_pdf (Tl1_2)");
        const Spack lstarn(lcond);

        qs1 = sp(0.622)*esval1_1/ekat::max(esval1_1, pval - esval1_1);
        beta1 = (rair/rv)*(lstarn/(rair*Tl1_1))*(lstarn/(cp*Tl1_1));

        // Only compute qs2 and beta2 if the two plumes are not equal
        const Smask condition = (Tl1_1 != Tl1_2);
        qs2 = qs1;
        beta2 = beta1;

        qs2.set(condition, sp(0.622)*esval1_2/ekat::max(esval1_2, pval - esval1_2));
        beta2.set(condition, (rair/rv)*(lstarn/(rair*Tl1_2))*(lstarn/(cp*Tl1_2)));
      }

      // Cloud computations

      // Compute s terms.
      Spack s1(0), std_s1(0), qn1(0), C1(0), ql1(0),
            s2(0), std_s2(0), qn2(0), C2(0), ql2(0);
      {
        const Scalar sqrt2(std::sqrt(Scalar(2.0))), sqrt2pi(std::sqrt(2*pi));

        // First plume
        const Spack cthl1=((1 + beta1*qw1_1)/ekat::square(1 + beta1*qs1))*(cp/lcond)*
                                            beta1*qs1*ekat::pow(pval/basepres, (rair/cp));
        const Spack cqt1 = 1/(1 + beta1*qs1);

        std_s1 = ekat::sqrt(ekat::max(0,
                                      ekat::square(cthl1)*thl2_1
                                      + ekat::square(cqt1)*qw2_1 - 2*cthl1*sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1));
        const auto std_s1_not_small = std_s1 > std::sqrt(std::numeric_limits<Scalar>::min()) * 100;
        s1 = qw1_1-qs1*((1 + beta1*qw1_1)/(1 + beta1*qs1));
        if (std_s1_not_small.any()) {
          C1.set(std_s1_not_small, sp(0.5)*(1 + ekat::erf(s1/(sqrt2*std_s1))));
        }
        C1.set(!std_s1_not_small && s1 > 0, 1);
        const auto std_s1_C1_not_small = std_s1_not_small && C1 != 0;
        if (std_s1_C1_not_small.any()) {
          qn1.set(std_s1_C1_not_small, s1*C1+(std_s1/sqrt2pi)*ekat::exp(-sp(0.5)*ekat::square(s1/std_s1)));
        }
        qn1.set(!std_s1_not_small && s1 > 0, s1);
        
        // Checking to prevent empty clouds
        const auto qn1_le_zero = qn1 <= 0;
        C1.set(qn1_le_zero,0);
        qn1.set(qn1_le_zero,0);
        
        ql1 = ekat::min(qn1, qw1_1);

        // Second plume
        // Only compute variables of the second plume if the two plumes are not equal
        const Smask equal(qw1_1==qw1_2 && thl2_1==thl2_2 && qs1==qs2);
        std_s2.set(equal, std_s1);
        s2.set(equal, s1);
        C2.set(equal, C1);
        qn2.set(equal, qn1);

        const Spack cthl2(!equal,
                          ((1 + beta2*qw1_2)/ekat::square(1 + beta2*qs2))*(cp/lcond)*
                           beta2*qs2*ekat::pow(pval/basepres, (rair/cp)));
        const Spack cqt2(!equal,
                         1/(1 + beta2*qs2));

        const auto nequal = !equal;
        if (nequal.any()) {
          std_s2.set(nequal,
                     ekat::sqrt(ekat::max(0,
                                          ekat::square(cthl2)*thl2_2
                                          + ekat::square(cqt2)*qw2_2 - 2*cthl2*sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1)));
          s2.set(nequal, qw1_2-qs2*((1 + beta2*qw1_2)/(1 + beta2*qs2)));
          const auto std_s2_not_small = std_s2 > std::sqrt(std::numeric_limits<Scalar>::min()) * 100;
          const auto nequal_std_s2_not_small = nequal && std_s2_not_small;
          if (nequal_std_s2_not_small.any()) {
            C2.set(nequal_std_s2_not_small, sp(0.5)*(1 + ekat::erf(s2/(sqrt2*std_s2))));
          }
          C2.set(nequal && !std_s2_not_small && s2 > 0, 1);
          const auto nequal_std_s2_C2_not_small = nequal_std_s2_not_small && C2 != 0;
          if (nequal_std_s2_C2_not_small.any()) {
            qn2.set(nequal_std_s2_C2_not_small, s2*C2+(std_s2/sqrt2pi)*ekat::exp(-sp(0.5)*ekat::square(s2/std_s2)));
          }
          qn2.set(nequal && !std_s2_not_small && s2 > 0, s2);
        }
        
        // Checking to prevent empty clouds
        const auto qn2_le_zero = qn2 <= 0;
        C2.set(qn2_le_zero,0);
        qn2.set(qn2_le_zero,0);

        ql2 = ekat::min(qn2, qw1_2);
      }

      // Compute SGS cloud fraction
      shoc_cldfrac(k) = ekat::min(1, a*C1 + (1 - a)*C2);

      // Compute SGS liquid water mixing ratio
      shoc_ql(k) = ekat::max(0, a*ql1 + (1 - a)*ql2);

      // Compute cloud liquid variance (CLUBB formulation, adjusted to SHOC parameters based)
      shoc_ql2(k) = ekat::max(0, a*(s1*ql1 + C1*ekat::square(std_s1))
                                 + (1 - a)*(s2*ql2 + C2*ekat::square(std_s2))
                                 - ekat::square(shoc_ql(k)));

      // Compute liquid water flux
      wqls(k) = a*((w1_1 - w_first)*ql1) + (1 - a)*((w1_2 - w_first)*ql2);

      // Compute the SGS buoyancy flux
      wthv_sec(k) = wthlsec + ((1 - epsterm)/epsterm)*basetemp*wqwsec
                   + ((lcond/cp)*ekat::pow(basepres/pval, (rair/cp))
                   - (1/epsterm)*basetemp)*wqls(k);
  });

  // Release temporary variables from the workspace
  workspace.template release_many_contiguous<6>(
    {&wthl_sec_zt, &wqw_sec_zt, &w3_zt,
     &thl_sec_zt, &qwthl_sec_zt, &qw_sec_zt});
}

} // namespace shoc
} // namespace scream

#endif
