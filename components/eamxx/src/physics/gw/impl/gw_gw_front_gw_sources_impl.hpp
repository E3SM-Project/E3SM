#ifndef GW_GW_FRONT_GW_SOURCES_IMPL_HPP
#define GW_GW_FRONT_GW_SOURCES_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_front_gw_sources. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_front_gw_sources(
// Inputs
const Int& pver,
const Int& pgwv,
const Int& ncol,
const Int& kbot,
const uview_1d<const Spack>& frontgf,
// Outputs
const uview_1d<Spack>& tau)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
