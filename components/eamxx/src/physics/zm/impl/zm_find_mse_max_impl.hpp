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
  const ZmRuntimeOpt& runtime_opt,
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
  //----------------------------------------------------------------------------
  // Purpose: find level of max moist static energy for parcel initialization
  //----------------------------------------------------------------------------
  // Local variables
  constexpr Real pergro_rhd_threshold = -1.e-4; // MSE difference threshold for perturbation growth test

  //----------------------------------------------------------------------------
  // initialize values
  mse_max_val = 0.0;
  msemax_klev = 0;
  const Int bot_layer = pver - 1 - runtime_opt.mx_bot_lyr_adj; // set lower limit to search for launch level with max MSE

  //----------------------------------------------------------------------------
  // Use parallel_reduce to find max moist static energy
  // The need for current mse_max_val means this loop has to be serial
  Kokkos::single(Kokkos::PerTeam(team), [&] {
    for (Int k = bot_layer; k >= num_msg; --k) {
      // calculate moist static energy
      const Real mse_env = PC::Cpair.value * temperature(k) +
                           PC::gravit.value * zmid(k) +
                           PC::LatVap.value * sp_humidity(k);
      // Reset max moist static energy level when relative difference exceeds threshold
      if (pergro_active) {
        const Real pergro_rhd = (mse_env - mse_max_val) / (mse_env + mse_max_val);
        if (k >= msemax_top_k && pergro_rhd > pergro_rhd_threshold) {
          mse_max_val = mse_env;
          msemax_klev = k;
        }
      }
      else {
        // find level and value of max moist static energy
        if (k >= msemax_top_k && mse_env > mse_max_val) {
          mse_max_val = mse_env;
          msemax_klev = k;
        }
      }
    }
  });
  team.team_barrier();
}

} // namespace zm
} // namespace scream

#endif
