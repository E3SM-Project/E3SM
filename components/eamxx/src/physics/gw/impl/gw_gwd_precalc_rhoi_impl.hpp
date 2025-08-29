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
  const Real& dt,
  const uview_1d<const Int>& tend_level,
  const uview_1d<const Real>& pmid,
  const uview_1d<const Real>& pint,
  const uview_1d<const Real>& t,
  const uview_1d<const Real>& gwut,
  const uview_1d<const Real>& ubm,
  const uview_1d<const Real>& nm,
  const uview_1d<const Real>& rdpm,
  const uview_1d<const Real>& c,
  const uview_1d<const Real>& q,
  const uview_1d<const Real>& dse,
  // Outputs
  const uview_1d<Real>& egwdffi,
  const uview_1d<Real>& qtgw,
  const uview_1d<Real>& dttdf,
  const uview_1d<Real>& dttke,
  const uview_1d<Real>& ttgw)
{

}

} // namespace gw
} // namespace scream

#endif
