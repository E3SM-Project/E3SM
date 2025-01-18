#ifndef P3_CLDLIQ_IMM_FREEZING_IMPL_HPP
#define P3_CLDLIQ_IMM_FREEZING_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "p3_subgrid_variance_scaling_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 contact and immersion freezing droplets function.
 * Clients should NOT #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cldliq_immersion_freezing(
  const Spack& T_atm, const Spack& lamc,
  const Spack& mu_c, const Spack& cdist1,
  const Spack& qc_incld, const Spack& inv_qc_relvar,
  Spack& qc2qi_hetero_freeze_tend, Spack& nc2ni_immers_freeze_tend,
  const P3Runtime& runtime_options,
  const Smask& context)
{
  constexpr Scalar qsmall   = C::QSMALL;
  constexpr Scalar T_rainfrz = C::T_rainfrz;
  constexpr Scalar T_zerodegc = C::T_zerodegc;
  constexpr Scalar CONS5    = C::CONS5;
  constexpr Scalar CONS6    = C::CONS6;
  const Scalar immersion_freezing_exponent =
      runtime_options.immersion_freezing_exponent;

  const auto qc_not_small_and_t_freezing = (qc_incld >= qsmall) &&
                                           (T_atm <= T_rainfrz) && context;
  if (qc_not_small_and_t_freezing.any()) {
    Spack expAimmDt, inv_lamc3;
    expAimmDt.set(qc_not_small_and_t_freezing,
                  exp(immersion_freezing_exponent * (T_zerodegc - T_atm)));
    inv_lamc3.set(qc_not_small_and_t_freezing, cube(1/lamc));

    Spack sgs_var_coef;
    // sgs_var_coef = subgrid_variance_scaling(inv_qc_relvar, 2);
    sgs_var_coef = 1;

    qc2qi_hetero_freeze_tend.set(qc_not_small_and_t_freezing,
               sgs_var_coef * CONS6 * cdist1 * tgamma(7+mu_c) * expAimmDt *
               square(inv_lamc3));
    nc2ni_immers_freeze_tend.set(qc_not_small_and_t_freezing,
               CONS5 * cdist1 * tgamma(sp(4.0)+mu_c) * expAimmDt * inv_lamc3);
  }
}

} // namespace p3
} // namespace scream

#endif
