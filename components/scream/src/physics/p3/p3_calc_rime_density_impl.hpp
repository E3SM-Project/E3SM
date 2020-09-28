#ifndef P3_CALC_RIME_DENSITY_IMPL_HPP
#define P3_CALC_RIME_DENSITY_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of rime density calculation function.
 * Clients should NOT #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::calc_rime_density(
  const Spack& T_atm, const Spack& rhofaci,
  const Spack& table_val_qi_fallspd, const Spack& acn,
  const Spack& lamc, const Spack& mu_c,
  const Spack& qc_incld, const Spack& qc2qi_collect_tend,
  Spack& vtrmi1, Spack& rho_qm_cloud,
  const Smask& context)
{
  constexpr Scalar qsmall   = C::QSMALL;
  constexpr Scalar T_zerodegc = C::T_zerodegc;
  constexpr Scalar bcn      = C::bcn;

  const auto qc2qi_collect_tend_not_small_and_t_freezing = (qc2qi_collect_tend >= qsmall) &&
                                              (T_atm < T_zerodegc) && context;

  // NOTE: applicable for cloud only; modify when rain is added back
  if (qc2qi_collect_tend_not_small_and_t_freezing.any()) {
    // Mass-weighted mean ice fallspeed
    vtrmi1.set(qc2qi_collect_tend_not_small_and_t_freezing, table_val_qi_fallspd * rhofaci);

    // If qc_incld isn't small, compute the rime density.
    const auto qc2qi_collect_tend_and_qc_not_small_and_t_freezing =
      qc2qi_collect_tend_not_small_and_t_freezing && (qc_incld >= qsmall);
    if (qc2qi_collect_tend_and_qc_not_small_and_t_freezing.any()) {

      // Droplet fall speed (using Stokes' formulation, with analytic soln).
      Spack Vt_qc = acn * tgamma(4+bcn+mu_c) /
        (pow(lamc, bcn) * tgamma(4+mu_c));

      // Use mass-weighted mean size
      Spack D_c = (4 + mu_c) / lamc;
      Spack V_impact = abs(vtrmi1-Vt_qc);
      Spack inv_Tc = 1/min(sp(-0.001), T_atm-T_zerodegc);
      Spack Ri = max(1, min(sp(-0.5e+6) * D_c * V_impact * inv_Tc, 12));

      const auto Ri_le_8 = (Ri <= sp(8.0));
      rho_qm_cloud.set(qc2qi_collect_tend_and_qc_not_small_and_t_freezing and Ri_le_8,
                    (sp(0.051) + sp(0.114)*Ri - sp(0.0055)*square(Ri))*1000);

      // For Ri > 8, assume a linear fit between 8 and 12.
      // rhorime = 900 kg m-3 at Ri = 12
      // This is somewhat ad-hoc but allows a smoother transition
      // in rime density up to wet growth.
      rho_qm_cloud.set(qc2qi_collect_tend_and_qc_not_small_and_t_freezing and !Ri_le_8,
                    sp(611.) + sp(72.25) * (Ri-8));
    }

    // If qc_incld is small, just set rho_qm_cloud to 400.
    const auto qc2qi_collect_tend_not_small_and_qc_small_and_t_freezing =
      qc2qi_collect_tend_not_small_and_t_freezing && (qc_incld < qsmall);
    rho_qm_cloud.set(qc2qi_collect_tend_not_small_and_qc_small_and_t_freezing, 400);
  }

  // Handle the cases we haven't handled above.
  vtrmi1.set(!qc2qi_collect_tend_not_small_and_t_freezing && context, 0); // no velocity if no ice
  rho_qm_cloud.set(!qc2qi_collect_tend_not_small_and_t_freezing && context, 400);
}

} // namespace p3
} // namespace scream

#endif
