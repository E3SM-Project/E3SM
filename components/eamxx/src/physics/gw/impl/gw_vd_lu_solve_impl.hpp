#ifndef GW_VD_LU_SOLVE_IMPL_HPP
#define GW_VD_LU_SOLVE_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw vd_lu_solve. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::vd_lu_solve(
  // Inputs
  const Int& pver,
  const uview_1d<const Real>& decomp_ca,
  const uview_1d<const Real>& decomp_cc,
  const uview_1d<const Real>& decomp_dnom,
  const uview_1d<const Real>& decomp_ze,
  const Int& ntop,
  const Int& nbot,
  const uview_1d<const Real>& cd_top,
  // Inputs/Outputs
  const uview_1d<Real>& q)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace gw
} // namespace scream

#endif
