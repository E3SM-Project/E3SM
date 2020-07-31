#ifndef P3_CLOUD_RAIN_ACC_IMPL_HPP
#define P3_CLOUD_RAIN_ACC_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "p3_subgrid_variance_scaling_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 cloud rain accretion function. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cloud_rain_accretion(
  const Spack& rho, const Spack& inv_rho,
  const Spack& qc_incld, const Spack& nc_incld,
  const Spack& qr_incld, const Spack& qc_relvar,
  Spack& qcacc, Spack& ncacc,
  const Smask& context)
{
  constexpr Scalar qsmall = C::QSMALL;

  Spack sgs_var_coef;
  sgs_var_coef = subgrid_variance_scaling(qc_relvar, sp(1.15) );

  const auto qr_and_qc_not_small = (qr_incld >= qsmall) && (qc_incld >= qsmall) && context;
  if (qr_and_qc_not_small.any()) {
    // Khroutdinov and Kogan (2000)
    qcacc.set(qr_and_qc_not_small,
              sgs_var_coef * sp(67.0) * pow(qc_incld * qr_incld, sp(1.15)));
    ncacc.set(qr_and_qc_not_small, qcacc * nc_incld / qc_incld);

    qcacc.set(ncacc == 0 && context, 0);
    ncacc.set(qcacc == 0 && context, 0);
  }
}

} // namespace p3
} // namespace scream

#endif
