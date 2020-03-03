#ifndef P3_FUNCTIONS_CLDLIQ_IMM_FREEZING_IMPL_HPP
#define P3_FUNCTIONS_CLDLIQ_IMM_FREEZING_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 contact and immersion freezing droplets function.
 * Clients should NOT #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cldliq_immersion_freezing(const Spack& t, const Spack& lamc,
                            const Spack& mu_c, const Spack& cdist1,
                            const Spack& qc_incld, Spack& qcheti, Spack& ncheti)
{
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar AIMM = C::AIMM;
  constexpr Scalar RainFrze = C::RainFrze;
  constexpr Scalar ZeroDegC = C::ZeroDegC;
  constexpr Scalar CONS5 = C::CONS5;
  constexpr Scalar CONS6 = C::CONS6;

  const auto qc_not_small_and_t_freezing = (qc_incld >= qsmall) &&
                                           (t <= RainFrze);
  if (qc_not_small_and_t_freezing.any()) {
    Spack expAimmDt, inv_lamc3;
    expAimmDt.set(qc_not_small_and_t_freezing, exp(AIMM * (ZeroDegC-t)));
    inv_lamc3.set(qc_not_small_and_t_freezing, cube(1/lamc));
    qcheti.set(qc_not_small_and_t_freezing,
               CONS6 * cdist1 * tgamma(sp(7.0)+mu_c) * expAimmDt *
               square(inv_lamc3));
    ncheti.set(qc_not_small_and_t_freezing,
               CONS5 * cdist1 * tgamma(sp(4.0)+mu_c) * expAimmDt * inv_lamc3);
  }
}

} // namespace p3
} // namespace scream

#endif
