#ifndef GW_GW_HEATING_DEPTH_IMPL_HPP
#define GW_GW_HEATING_DEPTH_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_heating_depth. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_heating_depth(
// Inputs
const Int& pver,
const Int& ncol,
const Spack& maxq0_conversion_factor,
const Spack& hdepth_scaling_factor,
const bool& use_gw_convect_old,
const uview_1d<const Spack>& zm,
const uview_1d<const Spack>& netdt,
// Outputs
const uview_1d<Int>& mini,
const uview_1d<Int>& maxi,
const uview_1d<Spack>& hdepth,
const uview_1d<Spack>& maxq0_out,
const uview_1d<Spack>& maxq0)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
