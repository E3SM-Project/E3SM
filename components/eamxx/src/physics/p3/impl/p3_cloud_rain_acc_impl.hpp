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
  const Spack& qr_incld, const Spack& inv_qc_relvar,
  Spack& qc2qr_accret_tend, Spack& nc_accret_tend,
  const P3Runtime& runtime_options,
  const Smask& context)
{
  constexpr Scalar qsmall = C::QSMALL;

  const Scalar accretion_prefactor   = runtime_options.accretion_prefactor;
  const Scalar accretion_qc_exponent = runtime_options.accretion_qc_exponent;
  const Scalar accretion_qr_exponent = runtime_options.accretion_qr_exponent;

  Spack sgs_var_coef;
  // sgs_var_coef = subgrid_variance_scaling(inv_qc_relvar, sp(1.15) );
  sgs_var_coef = 1;

  const auto qr_and_qc_not_small = (qr_incld >= qsmall) && (qc_incld >= qsmall) && context;
  if(qr_and_qc_not_small.any()) {
    // Khroutdinov and Kogan (2000)
    // TODO: always default to second branch after BFB stuff is addressed
    if(accretion_qc_exponent == accretion_qr_exponent) {
      qc2qr_accret_tend.set(qr_and_qc_not_small,
                            sgs_var_coef * accretion_prefactor *
                                pow(qc_incld * qr_incld, accretion_qr_exponent));
    } else {
      qc2qr_accret_tend.set(qr_and_qc_not_small,
                            sgs_var_coef * accretion_prefactor *
                                pow(qc_incld, accretion_qc_exponent) *
                                pow(qr_incld, accretion_qr_exponent));
    }
    nc_accret_tend.set(qr_and_qc_not_small,
                       qc2qr_accret_tend * nc_incld / qc_incld);

    qc2qr_accret_tend.set(nc_accret_tend == 0 && context, 0);
    nc_accret_tend.set(qc2qr_accret_tend == 0 && context, 0);
  }
}

} // namespace p3
} // namespace scream

#endif
