#ifndef GW_VD_LU_DECOMP_IMPL_HPP
#define GW_VD_LU_DECOMP_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw vd_lu_decomp. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::vd_lu_decomp(
  // Inputs
  const Int& pver,
  const uview_1d<const Real>& ksrf,
  const uview_1d<const Real>& kv,
  const uview_1d<const Real>& tmpi,
  const uview_1d<const Real>& rpdel,
  const Real& ztodt,
  const Real& gravit,
  const uview_1d<const Real>& cc_top,
  const Int& ntop,
  const Int& nbot,
  const uview_1d<const Real>& cpairv,
  // Outputs
  const uview_1d<Real>& decomp_ca,
  const uview_1d<Real>& decomp_cc,
  const uview_1d<Real>& decomp_dnom,
  const uview_1d<Real>& decomp_ze)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
