#ifndef P3_RAIN_IMM_FREEZING_IMPL_HPP
#define P3_RAIN_IMM_FREEZING_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 rain immersion freezing function. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_immersion_freezing(const Spack& T_atm, const Spack& lamr,
                          const Spack& mu_r, const Spack& cdistr,
                          const Spack& qr_incld, Spack& qr2qi_immers_freeze_tend, Spack& nr2ni_immers_freeze_tend,
                          const Smask& context)
{
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar T_rainfrz = C::T_rainfrz;
  constexpr Scalar T_zerodegc = C::T_zerodegc;
  constexpr Scalar AIMM = C::AIMM;
  constexpr Scalar CONS5 = C::CONS5;
  constexpr Scalar CONS6 = C::CONS6;

  const auto qr_not_small_and_t_freezing = (qr_incld >= qsmall) &&
                                           (T_atm <= T_rainfrz) && context;
  if (qr_not_small_and_t_freezing.any()) {
    qr2qi_immers_freeze_tend.set(qr_not_small_and_t_freezing,
               CONS6 *
               exp(log(cdistr) + log(tgamma(sp(7.)+mu_r)) - sp(6.)*log(lamr)) *
               exp(AIMM*(T_zerodegc-T_atm)));
    nr2ni_immers_freeze_tend.set(qr_not_small_and_t_freezing,
               CONS5 *
               exp(log(cdistr) + log(tgamma(sp(4.)+mu_r)) - sp(3.)*log(lamr)) *
               exp(AIMM*(T_zerodegc-T_atm)));
  }
}

} // namespace p3
} // namespace scream

#endif
