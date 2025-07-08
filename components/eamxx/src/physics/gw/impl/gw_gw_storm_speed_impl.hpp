#ifndef GW_GW_STORM_SPEED_IMPL_HPP
#define GW_GW_STORM_SPEED_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_storm_speed. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_storm_speed(
// Inputs
const Int& pver,
const Int& ncol,
const Spack& storm_speed_min,
const uview_1d<const Spack>& ubm,
const uview_1d<const Int>& mini,
const uview_1d<const Int>& maxi,
// Outputs
const uview_1d<Int>& storm_speed,
const uview_1d<Spack>& uh,
const uview_1d<Spack>& umin,
const uview_1d<Spack>& umax)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
