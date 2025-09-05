#ifndef GW_GW_BERES_SRC_IMPL_HPP
#define GW_GW_BERES_SRC_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_beres_src. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_beres_src(
// Inputs
const Int& pver,
const Int& pgwv,
const Int& ncol,
const uview_1d<const Spack>& lat,
const uview_1d<const Spack>& u,
const uview_1d<const Spack>& v,
const uview_1d<const Spack>& netdt,
const uview_1d<const Spack>& zm,
// Outputs
const uview_1d<Int>& src_level,
const uview_1d<Int>& tend_level,
const uview_1d<Spack>& tau,
const uview_1d<Spack>& ubm,
const uview_1d<Spack>& ubi,
const uview_1d<Spack>& xv,
const uview_1d<Spack>& yv,
const uview_1d<Spack>& c,
const uview_1d<Spack>& hdepth,
const uview_1d<Spack>& maxq0_out,
// Inputs
const Spack& maxq0_conversion_factor,
const Spack& hdepth_scaling_factor,
const Spack& hdepth_min,
const Spack& storm_speed_min,
const bool& use_gw_convect_old)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
