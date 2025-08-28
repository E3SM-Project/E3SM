#ifndef P3_DROPLET_SELF_COLL_IMPL_HPP
#define P3_DROPLET_SELF_COLL_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 droplet self collection function. Clients should NOT
 * #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::droplet_self_collection(
  const Spack& rho, const Spack& inv_rho,
  const Spack& qc_incld, const Spack& mu_c,
  const Spack& nu, const Spack& nc2nr_autoconv_tend, Spack& nc_selfcollect_tend,
  const Smask& context)
{
  constexpr Scalar qsmall = C::QSMALL;

  const auto qc_not_small = (qc_incld >= qsmall) && context;
  if (qc_not_small.any()) {
    if(use_KK){
    // Khroutdinov and Kogan (2000)
      nc_selfcollect_tend.set(qc_not_small, 0);}
    else {
      nc_selfcollect_tend.set(qc_not_small,-kc*pow(1.e-3*rho*qc_incld,2)*(nu+2)/(nu+1)*         
              1.e+6*inv_rho+nc2nr_autoconv_tend;}
  }
}

} // namespace p3
} // namespace scream

#endif
