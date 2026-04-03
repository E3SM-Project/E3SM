#ifndef ZM_ZM_CALC_OUTPUT_TEND_IMPL_HPP
#define ZM_ZM_CALC_OUTPUT_TEND_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_calc_output_tend. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_calc_output_tend(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& msg, // number of levels to ignore at model top
  const Int& jt, // level index of updraft top
  const Int& mx, // level index of updraft base
  const Real& dsubcld, // sub-cloud layer thickness
  const uview_1d<const Real>& p_del, // pressure thickness
  const uview_1d<const Real>& s_int, // ambient interface dry energy
  const uview_1d<const Real>& q_int, // ambient interface specific humidity
  const uview_1d<const Real>& s_upd, // updraft dry energy
  const uview_1d<const Real>& q_upd, // updraft specific humidity
  const uview_1d<const Real>& mflx_up, // updraft mass flux
  const uview_1d<const Real>& detr_up, // updraft detrainment
  const uview_1d<const Real>& mflx_dn, // downdraft mass flux
  const uview_1d<const Real>& s_dnd, // downdraft dry energy
  const uview_1d<const Real>& q_dnd, // downdraft specific humidity
  const uview_1d<const Real>& ql, // cloud liquid water
  const uview_1d<const Real>& evp, // evaporation
  const uview_1d<const Real>& cu, // updraft condensation
  // Outputs
  const uview_1d<Real>& dsdt, // output tendency for dry energy
  const uview_1d<Real>& dqdt, // output tendency for specific humidity
  const uview_1d<Real>& dl) // output tendency for cloud liquid water
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace zm
} // namespace scream

#endif
