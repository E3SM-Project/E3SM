#ifndef ZM_ZM_CONV_EVAP_IMPL_HPP
#define ZM_ZM_CONV_EVAP_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm zm_conv_evap. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::zm_conv_evap(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Real& time_step, // model time step                         [s]
  const uview_1d<const Real>& p_mid, // midpoint pressure                       [Pa]
  const uview_1d<const Real>& p_del, // layer thickness                         [Pa]
  const uview_1d<const Real>& t_mid, // temperature                             [K]
  const uview_1d<const Real>& q_mid, // water vapor                             [kg/kg]
  const uview_1d<const Real>& prdprec, // precipitation production                [kg/kg/s]
  const uview_1d<const Real>& cldfrc, // cloud fraction
  // Inputs/Outputs
  const uview_1d<Real>& tend_s, // heating rate                            [J/kg/s]
  const uview_1d<Real>& tend_q, // water vapor tendency                    [kg/kg/s]
  // Outputs
  const uview_1d<Real>& tend_s_snwprd, // Heating rate of snow production         [J/kg/s]
  const uview_1d<Real>& tend_s_snwevmlt, // Heating rate of snow evap/melt          [J/kg/s]
  // Inputs/Outputs
  Real& prec, // Convective-scale prec rate              [m/s]
  // Outputs
  Real& snow, // Convective-scale snow rate              [m/s]
  const uview_1d<Real>& ntprprd, // net precip production in layer          [kg/kg/s]
  const uview_1d<Real>& ntsnprd, // net snow production in layer            [kg/kg/s]
  const uview_1d<Real>& flxprec, // Convective flux of prec at interfaces   [kg/m2/s]
  const uview_1d<Real>& flxsnow) // Convective flux of snow at interfaces   [kg/m2/s])
{
  // TODO
  // Note, argument types may need tweaking. Generator is not always able to tell what needs to be packed
}

} // namespace zm
} // namespace scream

#endif
