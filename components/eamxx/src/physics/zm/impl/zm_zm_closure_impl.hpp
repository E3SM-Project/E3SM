#ifndef ZM_ZM_CLOSURE_IMPL_HPP
#define ZM_ZM_CLOSURE_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_closure. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_closure(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& msg, // number of levels to ignore at model top
  const Real& cape_threshold_in, // CAPE threshold for "cloud work function" (i.e. A)
  const Int& lcl, // index of lcl
  const Int& lel, // index of launch leve
  const Int& jt, // top of updraft
  const Int& mx, // base of updraft
  const Real& dsubcld, // thickness of subcloud layer
  const uview_1d<const Real>& z_mid, // altitude (m)
  const uview_1d<const Real>& z_int, // height of interface levels
  const uview_1d<const Real>& p_mid, // ambient pressure (mb)
  const uview_1d<const Real>& p_del, // pressure thickness of layers
  const uview_1d<const Real>& t_mid, // ambient temperature
  const uview_1d<const Real>& s_mid, // ambient dry energy (normalized)
  const uview_1d<const Real>& q_mid, // ambient specific humidity
  const uview_1d<const Real>& qs, // ambient saturation specific humidity
  const uview_1d<const Real>& ql, // ambient liquid water mixing ratio
  const uview_1d<const Real>& s_int, // env. normalized dry energy at intrfcs
  const uview_1d<const Real>& q_int, // environment specific humidity at interfaces
  const Real& t_pcl_lcl, // parcel temperature at LCL
  const uview_1d<const Real>& t_pcl, // parcel temperature
  const uview_1d<const Real>& q_pcl_sat, // parcel specific humidity
  const uview_1d<const Real>& s_upd, // updraft dry energy (normalized)
  const uview_1d<const Real>& q_upd, // updraft specific humidity
  const uview_1d<const Real>& mflx_net, // net convective mass flux
  const uview_1d<const Real>& detr_up, // detrainment from updraft
  const uview_1d<const Real>& mflx_up, // updraft mass flux
  const uview_1d<const Real>& mflx_dn, // dndraft mass flux
  const uview_1d<const Real>& q_dnd, // dndraft specific humidity
  const uview_1d<const Real>& s_dnd, // dndraft dry energy
  const Real& cape, // convective available potential energy
  // Outputs
  Real& cld_base_mass_flux) // cloud base mass flux
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace zm
} // namespace scream

#endif
