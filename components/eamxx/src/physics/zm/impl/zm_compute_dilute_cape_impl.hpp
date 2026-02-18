#ifndef ZM_COMPUTE_DILUTE_CAPE_IMPL_HPP
#define ZM_COMPUTE_DILUTE_CAPE_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm compute_dilute_cape. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_dilute_cape(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const Int& ncol, // number of atmospheric columns (actual)
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& num_cin, // num of negative buoyancy regions that are allowed before the conv. top and CAPE calc are completed
  const Int& num_msg, // index of highest level convection is allowed
  const uview_1d<const Real>& sp_humidity_in, // specific humidity
  const uview_1d<const Real>& temperature_in, // temperature
  const uview_1d<const Real>& zmid, // altitude/height at mid-levels
  const uview_1d<const Real>& pmid, // pressure at mid-levels
  const uview_1d<const Real>& pint, // pressure at interfaces
  const Int& pblt, // index of pbl top used as upper limit index of max MSE search
  const Real& tpert, // perturbation temperature by pbl processes
  const bool& calc_msemax_klev, // true for normal procedure, otherwise use prev_msemax_klev from 1st call
  const Int& prev_msemax_klev, // values of msemax_klev from previous call for dcape closure
  const bool& use_input_tq_mx, // if .true., use input values of prev_msemax_klev, q_mx, t_mx in the CAPE calculation
  // Inputs/Outputs
  const uview_1d<Real>& parcel_qsat, // parcel saturation mixing ratio
  Int& msemax_klev, // index of max MSE at parcel launch level
  Int& lcl_klev, // index of lifting condensation level (i.e. cloud bottom)
  Int& eql_klev, // index of equilibrium level (i.e. cloud top)
  Real& cape, // convective available potential energy
  Real& q_mx, // specified sp humidity to apply at level of max MSE if use_input_tq_mx=.true.
  Real& t_mx, // specified temperature to apply at level of max MSE if use_input_tq_mx=.true.)
  // Outputs
  const uview_1d<Real>& parcel_temp, // parcel temperature
  Real& lcl_temperature) // lifting condensation level (LCL) temperature
{
}

} // namespace zm
} // namespace scream

#endif
