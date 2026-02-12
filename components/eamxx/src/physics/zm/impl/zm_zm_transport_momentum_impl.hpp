#ifndef ZM_ZM_TRANSPORT_MOMENTUM_IMPL_HPP
#define ZM_ZM_TRANSPORT_MOMENTUM_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_transport_momentum. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_transport_momentum(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const Int& pver, // number of mid-point levels
  const Int& pverp, // number of interface levels
  const uview_2d<const Real>& wind_in, // input Momentum array
  const Int& nwind, // number of tracers to transport
  const uview_1d<const Real>& mu, // mass flux up
  const uview_1d<const Real>& md, // mass flux down
  const uview_1d<const Real>& du, // mass detraining from updraft
  const uview_1d<const Real>& eu, // mass entraining from updraft
  const uview_1d<const Real>& ed, // mass entraining from downdraft
  const uview_1d<const Real>& dp, // gathered pressure delta between interfaces
  const Int& jt, // index of cloud top for each column
  const Int& mx, // index of cloud top for each column
  const Int& ideep, // gathering array
  const Int& il1g, // gathered min ncol index
  const Int& il2g, // gathered max ncol index
  const Real& dt, // time step in seconds : 2*delta_t
  // Outputs
  const uview_2d<Real>& wind_tend, // output momentum tendency
  const uview_2d<Real>& pguall, // apparent force from  updraft PG
  const uview_2d<Real>& pgdall, // apparent force from  downdraft PG
  const uview_2d<Real>& icwu, // in-cloud winds in updraft
  const uview_2d<Real>& icwd, // in-cloud winds in downdraft
  const uview_1d<Real>& seten) // dry energy tendency)
{
  
}

} // namespace zm
} // namespace scream

#endif
