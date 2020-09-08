#ifndef P3_IMPOSE_MAX_TOTAL_NI_IMPL_HPP
#define P3_IMPOSE_MAX_TOTAL_NI_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::impose_max_total_ni(
  Spack& ni_local, const Scalar& max_total_ni, const Spack& inv_rho_local,
  const Smask& context)
{
  //--------------------------------------------------------------------------------
  // Impose maximum total ice number concentration (total of all ice categories).
  // If the sum of all ni(:) exceeds maximum allowable, each category to preserve
  // ratio of number between categories.
  //--------------------------------------------------------------------------------

  const auto ni_not_small = ni_local >= 1e-20 && context;
  if (ni_not_small.any()){
    ni_local.set(ni_not_small, ni_local*min(max_total_ni * inv_rho_local/ni_local, 1));
  }
}

} // namespace p3
} // namespace scream

#endif
