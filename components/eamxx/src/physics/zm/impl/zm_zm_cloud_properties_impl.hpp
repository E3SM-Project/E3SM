#ifndef ZM_ZM_CLOUD_PROPERTIES_IMPL_HPP
#define ZM_ZM_CLOUD_PROPERTIES_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_cloud_properties. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_cloud_properties(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& msg, // number of levels to ignore at model top
  const Int& limcnv, // convection limiting level
  const uview_1d<const Real>& p_mid, // env pressure at mid-point
  const uview_1d<const Real>& z_mid, // env altitude at mid-point
  const uview_1d<const Real>& z_int, // env altitude at interface
  const uview_1d<const Real>& t_mid, // env temperature
  const uview_1d<const Real>& s_mid, // env dry energy of env [K] (normalized)
  const uview_1d<const Real>& s_int, // interface values of dry stat energy
  const uview_1d<const Real>& q_mid, // env specific humidity
  const Real& landfrac, // Land fraction
  const Real& tpert_g, // PBL temperature perturbation
  const Int& jb, // updraft base level
  const Int& lel, // updraft parcel launch level
  // Outputs
  Int& jt, // updraft plume top
  Int& jlcl, // updraft lifting cond level
  Int& j0, // level where detrainment begins (starting at h_env_min)
  Int& jd, // level of downdraft
  const uview_1d<Real>& mflx_up, // updraft mass flux
  const uview_1d<Real>& entr_up, // entrainment rate of updraft
  const uview_1d<Real>& detr_up, // detrainement rate of updraft
  const uview_1d<Real>& mflx_dn, // downdraft mass flux
  const uview_1d<Real>& entr_dn, // downdraft entrainment rate
  const uview_1d<Real>& mflx_net, // net mass flux
  const uview_1d<Real>& s_upd, // updraft dry energy [K] (normalized)
  const uview_1d<Real>& q_upd, // updraft specific humidity [kg/kg]
  const uview_1d<Real>& ql, // updraft liq water
  const uview_1d<Real>& s_dnd, // dndraft dry energy [K] (normalized)
  const uview_1d<Real>& q_dnd, // dndraft specific humidity [kg/kg]
  const uview_1d<Real>& qst, // env saturation mixing ratio
  const uview_1d<Real>& cu, // condensation rate
  const uview_1d<Real>& evp, // evaporation rate
  const uview_1d<Real>& pflx, // precipitation flux thru layer
  const uview_1d<Real>& rprd) // rate of production of precip at that layer
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace zm
} // namespace scream

#endif
