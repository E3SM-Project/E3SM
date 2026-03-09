#ifndef ZM_ZM_CONV_MCSP_CALCULATE_SHEAR_IMPL_HPP
#define ZM_ZM_CONV_MCSP_CALCULATE_SHEAR_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_conv_mcsp_calculate_shear. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_conv_mcsp_calculate_shear(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const uview_1d<const Real>& state_pmid, // physics state mid-point pressure
  const uview_1d<const Real>& state_u, // physics state u momentum
  const uview_1d<const Real>& state_v, // physics state v momentum
  // Outputs
  Real& mcsp_shear)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace zm
} // namespace scream

#endif
