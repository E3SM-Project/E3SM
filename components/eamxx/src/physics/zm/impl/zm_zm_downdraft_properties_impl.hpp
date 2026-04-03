#ifndef ZM_ZM_DOWNDRAFT_PROPERTIES_IMPL_HPP
#define ZM_ZM_DOWNDRAFT_PROPERTIES_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_downdraft_properties. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_downdraft_properties(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& msg, // number of levels to ignore at model top
  const Int& jb, // updraft base level
  // Inputs/Outputs
  Int& jt, // updraft top level
  // Inputs
  const Int& j0, // level where updraft begins detraining
  // Inputs/Outputs
  Int& jd, // level of downdraft
  // Inputs
  const uview_1d<const Real>& z_int, // env altitude at interface
  const uview_1d<const Real>& dz, // layer thickness
  const uview_1d<const Real>& s_mid, // env dry energy of env [K] (normalized)
  const uview_1d<const Real>& q_mid, // env specific humidity
  const uview_1d<const Real>& h_env, // ambient env moist stat energy
  const uview_1d<const Real>& lambda, // fractional entrainment
  const Real& lambda_max, // fractional entrainment max
  const uview_1d<const Real>& qsthat, // interface interpolated qst
  const uview_1d<const Real>& hsthat, // interface interpolated hst
  const uview_1d<const Real>& gamhat, // interface interpolated gamma
  const uview_1d<const Real>& rprd, // rate of production of precip at that layer
  const uview_1d<const Real>& mflx_up, // updraft mass flux
  // Inputs/Outputs
  const uview_1d<Real>& mflx_dn, // downdraft mass flux
  const uview_1d<Real>& entr_dn, // downdraft entrainment rate
  const uview_1d<Real>& s_dnd, // dndraft dry energy [K] (normalized)
  const uview_1d<Real>& q_dnd, // dndraft specific humidity [kg/kg]
  const uview_1d<Real>& h_dnd, // dndraft moist energy
  const uview_1d<Real>& q_dnd_sat, // dndraft saturation specific humdity
  const uview_1d<Real>& evp, // evaporation rate
  Real& totevp) // total evap   for dndraft proportionality factor - see eq (4.106)
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace zm
} // namespace scream

#endif
