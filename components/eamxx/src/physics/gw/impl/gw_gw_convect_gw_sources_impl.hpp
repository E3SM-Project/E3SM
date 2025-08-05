#ifndef GW_GW_CONVECT_GW_SOURCES_IMPL_HPP
#define GW_GW_CONVECT_GW_SOURCES_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_convect_gw_sources. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_convect_gw_sources(
// Inputs
const Int& pver,
const Int& pgwv,
const Int& ncol,
const uview_1d<const Spack>& lat,
const Spack& hdepth_min,
const uview_1d<const Spack>& hdepth,
const uview_1d<const Int>& mini,
const uview_1d<const Int>& maxi,
const uview_1d<const Spack>& netdt,
const uview_1d<const Spack>& uh,
const uview_1d<const Int>& storm_speed,
const uview_1d<const Spack>& maxq0,
const uview_1d<const Spack>& umin,
const uview_1d<const Spack>& umax,
// Outputs
const uview_1d<Spack>& tau)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
