#ifndef GW_GWD_PRECALC_RHOI_IMPL_HPP
#define GW_GWD_PRECALC_RHOI_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gwd_precalc_rhoi. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gwd_precalc_rhoi(
// Inputs
const Int& pver,
const Int& pgwv,
const Int& ncol,
const Int& ngwv,
const Spack& dt,
const uview_1d<const Int>& tend_level,
const uview_1d<const Spack>& pmid,
const uview_1d<const Spack>& pint,
const uview_1d<const Spack>& t,
const uview_1d<const Spack>& gwut,
const uview_1d<const Spack>& ubm,
const uview_1d<const Spack>& nm,
const uview_1d<const Spack>& rdpm,
const uview_1d<const Spack>& c,
const uview_1d<const Spack>& q,
const uview_1d<const Spack>& dse,
// Outputs
const uview_1d<Spack>& egwdffi,
const uview_1d<Spack>& qtgw,
const uview_1d<Spack>& dttdf,
const uview_1d<Spack>& dttke,
const uview_1d<Spack>& ttgw)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
