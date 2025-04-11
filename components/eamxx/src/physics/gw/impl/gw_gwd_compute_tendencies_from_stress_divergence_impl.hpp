#ifndef GW_GWD_COMPUTE_TENDENCIES_FROM_STRESS_DIVERGENCE_IMPL_HPP
#define GW_GWD_COMPUTE_TENDENCIES_FROM_STRESS_DIVERGENCE_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gwd_compute_tendencies_from_stress_divergence. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
KOKKOS_FUNCTION
void Functions<S,D>::gwd_compute_tendencies_from_stress_divergence(const Int& pver, const Int& -pgwv:pgwv, const Int& 0:pver, const Int& -ngwv:ngwv, const Int& ncol, const Int& ngwv, const bool& do_taper, const Spack& dt, const Spack& effgw, const uview_1d<const Int>& tend_level, const uview_1d<const Spack>& lat, const uview_1d<const Spack>& dpm, const uview_1d<const Spack>& rdpm, const uview_1d<const Spack>& c, const uview_1d<const Spack>& ubm, const uview_1d<const Spack>& t, const uview_1d<const Spack>& nm, const uview_1d<const Spack>& xv, const uview_1d<const Spack>& yv, const uview_1d<Spack>& tau, const uview_1d<Spack>& gwut, const uview_1d<Spack>& utgw, const uview_1d<Spack>& vtgw)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
