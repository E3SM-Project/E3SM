#ifndef P3_FUNCTIONS_ICE_SED_IMPL_HPP
#define P3_FUNCTIONS_ICE_SED_IMPL_HPP

#include "p3_functions.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 ice sedimentation function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>
::calc_bulk_rho_rime(
  const Smask& qi_gt_small, const Spack& qi_tot, Spack& qi_rim, Spack& bi_rim)
{
}

} // namespace p3
} // namespace scream

#endif
