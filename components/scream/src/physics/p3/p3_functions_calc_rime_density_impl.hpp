#ifndef P3_FUNCTIONS_CALC_RIME_DENSITY_IMPL_HPP
#define P3_FUNCTIONS_CALC_RIME_DENSITY_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 rime density calculation function.
 * Clients should NOT #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::calc_rime_density(const Spack& t, const Spack& rhofaci,
                    const Spack& f1pr02, const Spack& acn,
                    const Spack& lamc, const Spack& mu_c,
                    const Spack& qc_incld, const Spack& qccol,
                    Spack& vtrmi1, Spack& rhorime_c)
{
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar RainFrze = C::RainFrze;
  constexpr Scalar ZeroDegC = C::ZeroDegC;
  constexpr Scalar bcn = C::bcn;

  const auto qccol_not_small_and_t_freezing = (qccol >= qsmall) &&
                                              (t < RainFrze);
  const auto qccol_and_qc_not_small_and_t_freezing = (qccol >= qsmall) &&
                                                     (qc_incld >= qsmall) &&
                                                     (t < RainFrze);

  // NOTE: applicable for cloud only; modify when rain is added back
  if (qccol_not_small_and_t_freezing.any()) {
    // Mass-weighted mean ice fallspeed
    vtrmi1.set(qccol_not_small_and_t_freezing, f1pr02 * rhofaci);

    if (qccol_and_qc_not_small_and_t_freezing.any()) {

      // Droplet fall speed (using Stokes' formulation, with analytic soln).
      Spack Vt_qc;
      Vt_qc.set(qccol_and_qc_not_small_and_t_freezing,
                acn * tgamma(sp(4.0)+bcn+mu_c) /
                (pow(lamc, bcn) * tgamma(sp(4.0)+mu_c)));

      // Use mass-weighted mean size
      Spack D_c, V_impact, inv_Tc, Ri;
      D_c.set(qccol_and_qc_not_small_and_t_freezing, (sp(4.0) + mu_c) / lamc);
      V_impact.set(qccol_and_qc_not_small_and_t_freezing, abs(vtrmi1-Vt_qc));
      inv_Tc.set(qccol_and_qc_not_small_and_t_freezing,
                 1/min(-0.001, t-ZeroDegC));
      Ri.set(qccol_and_qc_not_small_and_t_freezing,
             max(1, min(-0.5e+6 * D_c * V_impact * inv_Tc, sp(12.0))));
      const auto Ri_le_8 = (Ri <= sp(8.0));
      rhorime_c.set(qccol_and_qc_not_small_and_t_freezing and Ri_le_8,
                    (sp(0.051) + sp(0.114)*Ri - sp(0.0055)*square(Ri))*sp(1000.0));

      // For Ri > 8, assume a linear fit between 8 and 12,
      // rhorime = 900 kg m-3 at Ri = 12
      // this is somewhat ad-hoc but allows a smoother transition
      // in rime density up to wet growth
      rhorime_c.set(qccol_and_qc_not_small_and_t_freezing and not Ri_le_8,
                    sp(611.)+sp(72.25)*(Ri-sp(8.)));
    }
    rhorime_c.set(not qccol_and_qc_not_small_and_t_freezing, 0);
  }
  // Handle the cases we haven't handled above.
  rhorime_c.set(not qccol_not_small_and_t_freezing, sp(400.0));
}

} // namespace p3
} // namespace scream

#endif
