#ifndef P3_FUNCTIONS_CLOUD_RAIN_ACC_IMPL_HPP
#define P3_FUNCTIONS_CLOUD_RAIN_ACC_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 cloud rain accretion function. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::cloud_rain_accretion(const Spack& rho, const Spack& inv_rho,
                       const Spack& qc_incld, const Spack& nc_incld,
                       const Spack& qr_incld, Spack& qcacc, Spack& ncacc)
{
  constexpr Scalar qsmall = C::QSMALL;

  const auto qr_and_qc_not_small = (qr_incld >= qsmall) && (qc_incld >= qsmall);
  if (qr_and_qc_not_small.any()) {
    // Khroutdinov and Kogan (2000)
    qcacc.set(qr_and_qc_not_small,
              sp(67.0) * pow(qc_incld * qr_incld, sp(1.15)));
    ncacc.set(qr_and_qc_not_small, qcacc * nc_incld / qc_incld);

    qcacc.set(ncacc == 0, 0);
    ncacc.set(qcacc == 0, 0);
  }
}

} // namespace p3
} // namespace scream

#endif
