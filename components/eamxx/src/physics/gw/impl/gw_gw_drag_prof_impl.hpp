#ifndef GW_GW_DRAG_PROF_IMPL_HPP
#define GW_GW_DRAG_PROF_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_drag_prof. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_drag_prof(
// Inputs
const Int& pver,
const Int& pgwv,
const Int& ncol,
const Int& ngwv,
const uview_1d<const Int>& src_level,
const uview_1d<const Int>& tend_level,
const bool& do_taper,
const Spack& dt,
const uview_1d<const Spack>& lat,
const uview_1d<const Spack>& t,
const uview_1d<const Spack>& ti,
const uview_1d<const Spack>& pmid,
const uview_1d<const Spack>& pint,
const uview_1d<const Spack>& dpm,
const uview_1d<const Spack>& rdpm,
const uview_1d<const Spack>& piln,
const uview_1d<const Spack>& rhoi,
const uview_1d<const Spack>& nm,
const uview_1d<const Spack>& ni,
const uview_1d<const Spack>& ubm,
const uview_1d<const Spack>& ubi,
const uview_1d<const Spack>& xv,
const uview_1d<const Spack>& yv,
const Spack& effgw,
const uview_1d<const Spack>& c,
const uview_1d<const Spack>& kvtt,
const uview_1d<const Spack>& q,
const uview_1d<const Spack>& dse,
// Inputs/Outputs
const uview_1d<Spack>& tau,
// Outputs
const uview_1d<Spack>& utgw,
const uview_1d<Spack>& vtgw,
const uview_1d<Spack>& ttgw,
const uview_1d<Spack>& qtgw,
const uview_1d<Spack>& taucd,
const uview_1d<Spack>& egwdffi,
const uview_1d<Spack>& gwut,
const uview_1d<Spack>& dttdf,
const uview_1d<Spack>& dttke)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
