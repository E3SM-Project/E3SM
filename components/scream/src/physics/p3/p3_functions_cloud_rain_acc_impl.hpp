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
  constexpr int iparam = P3C::iparam;
}

} // namespace p3
} // namespace scream

#endif
