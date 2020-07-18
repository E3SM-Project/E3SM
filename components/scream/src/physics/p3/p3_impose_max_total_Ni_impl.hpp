#ifndef P3_IMPOSE_MAX_TOTAL_NI_IMPL_HPP
#define P3_IMPOSE_MAX_TOTAL_NI_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::impose_max_total_Ni(
  Spack& nitot_local, const Scalar& max_total_Ni, const Spack& inv_rho_local,
  const Smask& context)
{
  //--------------------------------------------------------------------------------
  // Impose maximum total ice number concentration (total of all ice categories).
  // If the sum of all nitot(:) exceeds maximum allowable, each category to preserve
  // ratio of number between categories.
  //--------------------------------------------------------------------------------

  const auto nitot_not_small = nitot_local >= 1e-20 && context;
  if (nitot_not_small.any()){
    nitot_local.set(nitot_not_small, nitot_local*min(max_total_Ni * inv_rho_local/nitot_local, 1));
  }
}

} // namespace p3
} // namespace scream

#endif
