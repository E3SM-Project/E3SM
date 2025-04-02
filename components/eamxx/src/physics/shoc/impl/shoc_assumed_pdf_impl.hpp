#ifndef SHOC_SHOC_ASSUMED_PDF_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp"

#include <iomanip>

namespace scream {
namespace shoc {

/*
Add some functions to avoid cuda compilation warnings
TODO: move this to ekat or something
*/
#ifdef KOKKOS_ENABLE_CUDA
KOKKOS_INLINE_FUNCTION constexpr Real safe_min() {
  return Kokkos::Experimental::norm_min_v<Real>;
}
#else
KOKKOS_INLINE_FUNCTION constexpr Real safe_min() {
  return std::numeric_limits<Real>::min();
}
#endif

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
  const Scalar&                dtime,
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
  const uview_1d<Spack>&       shoc_ql2,
  const uview_1d<Spack>&       shoc_cond,
  const uview_1d<Spack>&       shoc_evap)
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

    // Initialize shoc_cond and shoc_evap to zero
    shoc_cond(k) = 0;
    shoc_evap(k) = 0;

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
    Spack Skew_w, w1_1, w1_2, w2_1, w2_2, a;
    shoc_assumed_pdf_vv_parameters(w_first, w2sec, w3var, w_tol_sqd, Skew_w, w1_1, w1_2, w2_1, w2_2, a);

    // Find parameters for thetal
    Spack thl1_1, thl1_2, thl2_1, thl2_2, sqrtthl2_1, sqrtthl2_2;
    shoc_assumed_pdf_thl_parameters(wthlsec,sqrtw2,sqrtthl,thlsec,thl_first,w1_1,w1_2,Skew_w,a,
                                    thl_tol, w_thresh,
                                    thl1_1, thl1_2, thl2_1, thl2_2, sqrtthl2_1, sqrtthl2_2);

    // Find parameters for total water mixing ratio
    Spack qw1_1, qw1_2, qw2_1, qw2_2, sqrtqw2_1, sqrtqw2_2;
    shoc_assumed_pdf_qw_parameters(wqwsec, sqrtw2, Skew_w, sqrtqt, qwsec, w1_2, w1_1, qw_first, a,
                                   rt_tol, w_thresh,
                                   qw1_1, qw1_2, qw2_1, qw2_2, sqrtqw2_1, sqrtqw2_2);

    // Convert from tilde variables to "real" variables
    shoc_assumed_pdf_tilde_to_real(w_first, sqrtw2, w1_1);
    shoc_assumed_pdf_tilde_to_real(w_first, sqrtw2, w1_2);

     // Find within-plume correlations.
      Spack r_qwthl_1;
      shoc_assumed_pdf_inplume_correlations(sqrtqw2_1, sqrtthl2_1, a, sqrtqw2_2, sqrtthl2_2, qwthlsec,
                                            qw1_1, qw_first, thl1_1, thl_first, qw1_2, thl1_2,
                                            r_qwthl_1);

      // Begin to compute cloud property statistics
      Spack Tl1_1, Tl1_2;
      shoc_assumed_pdf_compute_temperature(thl1_1,pval, Tl1_1);
      shoc_assumed_pdf_compute_temperature(thl1_2,pval, Tl1_2);

      const auto index_range = ekat::range<IntSmallPack>(k*Spack::n);
      const Smask active_entries = (index_range < nlev);

      // Do NaN Checking here
      const Smask is_nan_Tl1_1 = isnan(Tl1_1) && active_entries;
      const Smask is_nan_Tl1_2 = isnan(Tl1_2) && active_entries;
      if (is_nan_Tl1_1.any() || is_nan_Tl1_2.any()) {
        Kokkos::printf("WARNING: NaN Detected in Tl1_1 or Tl1_2!\n");
        for (int i=0; i<is_nan_Tl1_1.n; i++) {
          if (is_nan_Tl1_1[i] || is_nan_Tl1_2[i]) {
            Kokkos::printf(
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
	Kokkos::printf("WARNING: Tl1_1 has %d values <= allowable value.  Resetting to minimum value.\n",n_mask);
      }
      if( is_small_Tl1_2.any() ) {
        Tl1_2.set(is_small_Tl1_2,Tl_min);
        int n_mask = 0;
        for (int i=0; i<is_small_Tl1_2.n; i++) {
          if (is_small_Tl1_2[i]) {
            n_mask++;
          }
        }
	Kokkos::printf("WARNING: Tl1_2 has %d values <= allowable value.  Resetting to minimum value.\n",n_mask);
      }

      // Compute qs and beta
      Spack qs1, qs2, beta1, beta2;
      shoc_assumed_pdf_compute_qs(Tl1_1, Tl1_2, pval, active_entries, qs1, beta1, qs2, beta2);

      // Cloud computations

      // Compute s terms.
      Spack s1, std_s1, qn1, C1, ql1, s2, std_s2, qn2, C2, ql2;
      Spack dum;

      // First plume
      shoc_assumed_pdf_compute_s(qw1_1, qs1, beta1, pval, thl2_1, qw2_1,
                                 sqrtthl2_1, sqrtqw2_1, r_qwthl_1,
                                 s1, std_s1, qn1, C1);

      // Second plume
      // Only compute variables of the second plume if the two plumes are not equal
      const Smask equal(qw1_1==qw1_2 && thl2_1==thl2_2 && qs1==qs2);
      std_s2.set(equal, std_s1);
      s2.set(equal, s1);
      C2.set(equal, C1);
      qn2.set(equal, qn1);

      const auto nequal = !equal;
      if (nequal.any()) {
        shoc_assumed_pdf_compute_s(qw1_2, qs2, beta2, pval, thl2_2, qw2_2,
                                   sqrtthl2_2, sqrtqw2_2, r_qwthl_1,
                                   s2, std_s2, qn2, C2);
      }

      ql1 = ekat::min(qn1, qw1_1);
      ql2 = ekat::min(qn2, qw1_2);

      // Compute SGS cloud fraction
      shoc_cldfrac(k) = ekat::min(1, a*C1 + (1 - a)*C2);

      //Compute cond and evap tendencies 
      dum = ekat::max(0, a*ql1 + (1 - a)*ql2);
      shoc_cond(k) = ekat::max(0,(dum - shoc_ql(k))/dtime);
      shoc_evap(k) = ekat::max(0,(shoc_ql(k) - dum)/dtime);

      // Compute SGS liquid water mixing ratio
      shoc_assumed_pdf_compute_sgs_liquid(a, ql1, ql2, shoc_ql(k));

      // Compute cloud liquid variance (CLUBB formulation, adjusted to SHOC parameters based)
      shoc_assumed_pdf_compute_cloud_liquid_variance(a, s1, ql1, C1, std_s1,
                                                     s2, ql2, C2, std_s2, shoc_ql(k),
                                                     shoc_ql2(k));

      // Compute liquid water flux
      shoc_assumed_pdf_compute_liquid_water_flux(a, w1_1, w_first, ql1, w1_2, ql2, wqls(k));

      // Compute the SGS buoyancy flux
      shoc_assumed_pdf_compute_buoyancy_flux(wthlsec, wqwsec, pval, wqls(k),
                                             wthv_sec(k));
  });

  // Release temporary variables from the workspace
  workspace.template release_many_contiguous<6>(
    {&wthl_sec_zt, &wqw_sec_zt, &w3_zt,
     &thl_sec_zt, &qwthl_sec_zt, &qw_sec_zt});
}

} // namespace shoc
} // namespace scream

#endif
