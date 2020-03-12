#ifndef P3_FUNCTIONS_CALC_LIQ_RELAXATION_TIMESCALE_IMPL_HPP
#define P3_FUNCTIONS_CALC_LIQ_RELAXATION_TIMESCALE_IMPL_HPP

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
::calc_liq_relaxation_timescale(const view_2d_table& revap_table,
                                const Spack& rho, const Spack& f1r, const Spack& f2r,
                                const Spack& dv, const Spack& mu, const Spack& sc,
                                const Spack& mu_r, const Spack& lamr, const Spack& cdistr,
                                const Spack& cdist, const Spack& qr_incld, const Spack& qc_incld,
                                Spack& epsr, Spack& epsc)
{
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar pi = C::Pi;

  const auto qr_not_small = (qr_incld >= qsmall);
  Table3 table;
  lookup(qr_not_small, mu_r, lamr, table);
  epsr.set(qr_not_small, sp(2.0) * pi * cdistr * rho * dv *
           (f1r * tgamma(mu_r + sp(2.0))/lamr +
            f2r * sqrt(rho/mu) * cbrt(sc) *
            apply_table(qr_not_small, revap_table, table)));
  const auto qr_small = not qr_not_small;
  epsr.set(qr_small, 0);

  const auto qc_not_small = (qc_incld >= qsmall);
  epsc.set(qc_not_small, sp(2.0) * pi * rho * dv * cdist);
  const auto qc_small = not qc_not_small;
  epsc.set(qc_small, 0);
}

} // namespace p3
} // namespace scream

#endif
