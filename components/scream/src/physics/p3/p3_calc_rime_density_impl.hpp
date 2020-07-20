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
  const Spack& t, const Spack& rhofaci,
  const Spack& f1pr02, const Spack& acn,
  const Spack& lamc, const Spack& mu_c,
  const Spack& qc_incld, const Spack& qccol,
  Spack& vtrmi1, Spack& rhorime_c,
  const Smask& context)
{
  constexpr Scalar qsmall   = C::QSMALL;
  constexpr Scalar ZeroDegC = C::ZeroDegC;
  constexpr Scalar bcn      = C::bcn;

  const auto qccol_not_small_and_t_freezing = (qccol >= qsmall) &&
                                              (t < ZeroDegC) && context;

  // NOTE: applicable for cloud only; modify when rain is added back
  if (qccol_not_small_and_t_freezing.any()) {
    // Mass-weighted mean ice fallspeed
    vtrmi1.set(qccol_not_small_and_t_freezing, f1pr02 * rhofaci);

    // If qc_incld isn't small, compute the rime density.
    const auto qccol_and_qc_not_small_and_t_freezing =
      qccol_not_small_and_t_freezing && (qc_incld >= qsmall);
    if (qccol_and_qc_not_small_and_t_freezing.any()) {

      // Droplet fall speed (using Stokes' formulation, with analytic soln).
      Spack Vt_qc = acn * tgamma(4+bcn+mu_c) /
        (pow(lamc, bcn) * tgamma(4+mu_c));

      // Use mass-weighted mean size
      Spack D_c = (4 + mu_c) / lamc;
      Spack V_impact = abs(vtrmi1-Vt_qc);
      Spack inv_Tc = 1/min(sp(-0.001), t-ZeroDegC);
      Spack Ri = max(1, min(sp(-0.5e+6) * D_c * V_impact * inv_Tc, 12));

      const auto Ri_le_8 = (Ri <= sp(8.0));
      rhorime_c.set(qccol_and_qc_not_small_and_t_freezing and Ri_le_8,
                    (sp(0.051) + sp(0.114)*Ri - sp(0.0055)*square(Ri))*1000);

      // For Ri > 8, assume a linear fit between 8 and 12.
      // rhorime = 900 kg m-3 at Ri = 12
      // This is somewhat ad-hoc but allows a smoother transition
      // in rime density up to wet growth.
      rhorime_c.set(qccol_and_qc_not_small_and_t_freezing and !Ri_le_8,
                    sp(611.) + sp(72.25) * (Ri-8));
    }

    // If qc_incld is small, just set rhorime_c to 400.
    const auto qccol_not_small_and_qc_small_and_t_freezing =
      qccol_not_small_and_t_freezing && (qc_incld < qsmall);
    rhorime_c.set(qccol_not_small_and_qc_small_and_t_freezing, 400);
  }

  // Handle the cases we haven't handled above.
  vtrmi1.set(!qccol_not_small_and_t_freezing && context, 0); // no velocity if no ice
  rhorime_c.set(!qccol_not_small_and_t_freezing && context, 400);
}

} // namespace p3
} // namespace scream

#endif
