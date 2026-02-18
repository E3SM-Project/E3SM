#ifndef ZM_FIND_MSE_MAX_IMPL_HPP
#define ZM_FIND_MSE_MAX_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace zm {

/*
 * Implementation of zm find_mse_max. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::find_mse_max(
  // Inputs
  const MemberType& team,
  const Int& pver, // number of mid-point vertical levels
  const Int& num_msg, // number of missing moisture levels at the top of model
  const Int& msemax_top_k, // upper limit index of max MSE search
  const bool& pergro_active, // flag for perturbation growth test (pergro)
  const uview_1d<const Real>& temperature, // environement temperature
  const uview_1d<const Real>& zmid, // height/altitude at mid-levels
  const uview_1d<const Real>& sp_humidity, // specific humidity
  // Inputs/Outputs
  Int& msemax_klev, // index of max MSE at parcel launch level
  Real& mse_max_val) // value of max MSE at parcel launch level)
{
}

} // namespace zm
} // namespace scream

#endif
