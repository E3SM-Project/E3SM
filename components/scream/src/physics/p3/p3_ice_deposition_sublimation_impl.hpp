#ifndef P3_ICE_DEPOSITION_SUBLIMATION_IMPL_HPP
#define P3_ICE_DEPOSITION_SUBLIMATION_IMPL_HPP

#include "physics/p3/p3_functions.hpp"
#include "physics/share/physics_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_deposition_sublimation(
  const Spack& qi_incld, const Spack& ni_incld, const Spack& T_atm,   const Spack& qv_sat_l,
  const Spack& qv_sat_i,         const Spack& epsi,        const Spack& abi, const Spack& qv,
  Spack& qv2qi_vapdep_tend, Spack& qi2qv_sublim_tend, Spack& ni_sublim_tend, Spack& qc2qi_berg_tend,
  const Smask& context)
{
  constexpr Scalar QSMALL   = C::QSMALL;
  constexpr Scalar t_zerodegc = C::t_zerodegc;
  const auto oabi           = 1 / abi;

  const auto qi_incld_not_small = qi_incld >= QSMALL && context;
  const auto qi_incld_small     = qi_incld < QSMALL  && context;

  //Compute deposition/sublimation
  qv2qi_vapdep_tend.set(qi_incld_not_small,epsi * oabi * (qv - qv_sat_i));

  //Split into deposition or sublimation.
  //if "t" is greater than 0 degree celcius and qv2qi_vapdep_tend is positive
  const auto t_gt_t_zerodegc_pos_qv2qi_vapdep_tend = (T_atm < t_zerodegc && qv2qi_vapdep_tend > 0);

  qi2qv_sublim_tend.set(qi_incld_not_small && t_gt_t_zerodegc_pos_qv2qi_vapdep_tend, 0);

  //make qi2qv_sublim_tend positive for consistency with other evap/sub processes
  qi2qv_sublim_tend.set(qi_incld_not_small && !t_gt_t_zerodegc_pos_qv2qi_vapdep_tend, -min(qv2qi_vapdep_tend,0));
  qv2qi_vapdep_tend.set(qi_incld_not_small && !t_gt_t_zerodegc_pos_qv2qi_vapdep_tend, 0);

  //sublimation occurs @ any T. Not so for berg.
  const auto t_lt_t_zerodegc = T_atm < t_zerodegc;

  //Compute bergeron rate assuming cloud for whole step.
  qc2qi_berg_tend.set(qi_incld_not_small && t_lt_t_zerodegc, max(epsi*oabi*(qv_sat_l - qv_sat_i), 0));
  qc2qi_berg_tend.set(qi_incld_not_small && !t_lt_t_zerodegc, 0);

  if (qi_incld_not_small.any()) {
    ni_sublim_tend.set(qi_incld_not_small, qi2qv_sublim_tend*(ni_incld/qi_incld));
  }

  //if qi_incld is small (i.e. !qi_incld_not_small is true)
  qc2qi_berg_tend.set(qi_incld_small, 0);
  qv2qi_vapdep_tend.set (qi_incld_small, 0);
  qi2qv_sublim_tend.set (qi_incld_small, 0);
  ni_sublim_tend.set (qi_incld_small, 0);
}

} // namespace p3
} // namespace scream

#endif // P3_ICE_DEPOSITION_SUBLIMATION_IMPL_HPP
