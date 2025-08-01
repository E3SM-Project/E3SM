#ifndef GW_GWD_PROJECT_TAU_IMPL_HPP
#define GW_GWD_PROJECT_TAU_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gwd_project_tau. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gwd_project_tau(
// Inputs
const Int& pver,
const Int& pgwv,
const Int& ncol,
const uview_1d<const Int>& tend_level,
const uview_1d<const Spack>& tau,
const uview_1d<const Spack>& ubi,
const uview_1d<const Spack>& c,
const uview_1d<const Spack>& xv,
const uview_1d<const Spack>& yv,
// Outputs
const uview_1d<Spack>& taucd)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
