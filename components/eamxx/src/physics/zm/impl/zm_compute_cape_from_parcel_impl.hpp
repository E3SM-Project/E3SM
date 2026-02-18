#ifndef ZM_COMPUTE_CAPE_FROM_PARCEL_IMPL_HPP
#define ZM_COMPUTE_CAPE_FROM_PARCEL_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm compute_cape_from_parcel. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_cape_from_parcel(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& pverp, // number of interface vertical levels
  const Int& num_cin, // num of negative buoyancy regions that are allowed before the conv. top and CAPE calc are completed
  const Int& num_msg, // number of missing moisture levels at the top of model
  const uview_1d<const Real>& temperature, // temperature
  const uview_1d<const Real>& tv, // virtual temperature
  const uview_1d<const Real>& zmid, // height/altitude at mid-levels
  const uview_1d<const Real>& sp_humidity, // specific humidity
  const uview_1d<const Real>& pint, // pressure at interfaces
  const Int& msemax_klev, // index of max MSE at parcel launch level
  const Real& lcl_pmid, // lifting condensation level (LCL) pressure
  const Int& lcl_klev, // lifting condensation level (LCL) index
  // Inputs/Outputs
  const uview_1d<Real>& parcel_qsat, // parcel saturation mixing ratio
  const uview_1d<Real>& parcel_temp, // parcel temperature
  const uview_1d<Real>& parcel_vtemp, // parcel virtual temperature
  Int& eql_klev, // index of equilibrium level (i.e. cloud top)
  Real& cape) // convective available potential energy
{
}

} // namespace zm
} // namespace scream

#endif
