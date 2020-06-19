#ifndef P3_FUNCTIONS_ICE_RELAXATION_TIMESCALE_IMPL_HPP
#define P3_FUNCTIONS_ICE_RELAXATION_TIMESCALE_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_relaxation_timescale(
  const Spack& rho, const Spack& temp, const Spack& rhofaci, const Spack& f1pr05,
  const Spack& f1pr14, const Spack& dv, const Spack& mu, const Spack& sc,
  const Spack& qitot_incld, const Spack& nitot_incld,
  Spack& epsi, Spack& epsi_tot,
  const Smask& context)
{
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar tmelt = C::Tmelt;
  constexpr Scalar zero = C::ZERO;
  constexpr Scalar pi = C::Pi;

  const auto t_is_negative = temp < tmelt;
  const auto qitot_incld_ge_small = qitot_incld >= qsmall;

  const auto any_if = qitot_incld_ge_small && t_is_negative && context;

  /*!-----------------------------
   * calculate total inverse ice relaxation timescale combined for all ice categories
   * note 'f1pr' values are normalized, so we need to multiply by N
   */
  epsi.set(any_if,
           ((f1pr05+f1pr14*pack::cbrt(sc)*sqrt(rhofaci*rho/mu))*
            sp(2.0)*pi*rho*dv)*nitot_incld);

  epsi_tot.set(any_if, epsi_tot+epsi);

  epsi.set(!any_if && context, zero);
}

} // namespace p3
} // namespace scream

#endif
