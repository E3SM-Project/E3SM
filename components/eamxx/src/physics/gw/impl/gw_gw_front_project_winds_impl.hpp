#ifndef GW_GW_FRONT_PROJECT_WINDS_IMPL_HPP
#define GW_GW_FRONT_PROJECT_WINDS_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_front_project_winds. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_front_project_winds(
// Inputs
const Int& pver,
const Int& ncol,
const Int& kbot,
const uview_1d<const Spack>& u,
const uview_1d<const Spack>& v,
// Outputs
const uview_1d<Spack>& xv,
const uview_1d<Spack>& yv,
const uview_1d<Spack>& ubm,
const uview_1d<Spack>& ubi)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
