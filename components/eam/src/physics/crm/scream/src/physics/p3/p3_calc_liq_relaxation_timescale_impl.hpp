#ifndef P3_CALC_LIQ_RELAXATION_TIMESCALE_IMPL_HPP
#define P3_CALC_LIQ_RELAXATION_TIMESCALE_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 liquid relaxation timescale calculation.
 * Clients should NOT #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::calc_liq_relaxation_timescale(
  const view_2d_table& revap_table_vals,
  const Spack& rho, const Scalar& f1r, const Scalar& f2r,
  const Spack& dv, const Spack& mu, const Spack& sc,
  const Spack& mu_r, const Spack& lamr, const Spack& cdistr,
  const Spack& cdist, const Spack& qr_incld, const Spack& qc_incld,
  Spack& epsr, Spack& epsc,
  const Smask& context)
{
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar pi = C::Pi;

  const auto qr_not_small = (qr_incld >= qsmall) && context;
  epsr = 0;
  if (qr_not_small.any()) {
    Table3 table;
    lookup(mu_r, lamr, table, qr_not_small);
    epsr.set(qr_not_small, 2 * pi * cdistr * rho * dv *
             (f1r * tgamma(mu_r + 2)/lamr +
              f2r * sqrt(rho/mu) * cbrt(sc) *
              apply_table(revap_table_vals, table)));
  }

  const auto qc_not_small = (qc_incld >= qsmall) && context;
  epsc.set(qc_not_small, 2 * pi * rho * dv * cdist);
  const auto qc_small = !qc_not_small;
  epsc.set(qc_small && context, 0);
}

} // namespace p3
} // namespace scream

#endif
