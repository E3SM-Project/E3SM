#ifndef ZM_COMPUTE_DILUTE_PARCEL_IMPL_HPP
#define ZM_COMPUTE_DILUTE_PARCEL_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm compute_dilute_parcel. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::compute_dilute_parcel(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& num_msg, // number of missing moisture levels at the top of model
  const Int& klaunch, // index of parcel launch level based on max MSE
  const uview_1d<const Real>& pmid, // ambient env pressure at cell center
  const uview_1d<const Real>& temperature, // ambient env temperature at cell center
  const uview_1d<const Real>& sp_humidity, // ambient env specific humidity at cell center
  const Real& tpert, // PBL temperature perturbation
  const Int& pblt, // index of pbl depth
  // Inputs/Outputs
  const uview_1d<Real>& parcel_temp, // Parcel temperature
  const uview_1d<Real>& parcel_vtemp, // Parcel virtual temperature
  const uview_1d<Real>& parcel_qsat, // Parcel water vapour (sat value above lcl)
  Real& lcl_pmid, // lifting condensation level (LCL) pressure
  Real& lcl_temperature, // lifting condensation level (LCL) temperature
  Int& lcl_klev) // lifting condensation level (LCL) vertical index
{
}

} // namespace zm
} // namespace scream

#endif
