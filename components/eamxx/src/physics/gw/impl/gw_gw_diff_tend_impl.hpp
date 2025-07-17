#ifndef GW_GW_DIFF_TEND_IMPL_HPP
#define GW_GW_DIFF_TEND_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_diff_tend. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_diff_tend(
// Inputs
const Int& ncol,
const Int& pver,
const Int& kbot,
const Int& ktop,
const uview_1d<const Spack>& q,
const Spack& dt,
// Outputs
const uview_1d<Spack>& dq)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
