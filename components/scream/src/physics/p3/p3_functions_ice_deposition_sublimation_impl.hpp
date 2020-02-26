#ifndef P3_FUNCTIONS_ICE_DEPOSITION_SUBLIMATION_IMPL.HPP
#define P3_FUNCTIONS_ICE_DEPOSITION_SUBLIMATION_IMPL.HPP

#include "p3_functions.hpp"
#include "p3_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_deposition_sublimation(const Spack& qitot_incld, const Spack& nitot_incld, const Spack& t,   const Spack& qvs,
			     const Spack& qvi,         const Spack& epsi,        const Spack& abi, const Spack& qv,
			     Spack& qidep, Spack& qisub, Spack& nisub, Spack& qiberg)
{
  constexpr Scalar QSMALL   = C::QSMALL;
  constexpr Scalar ZERODEGC = C::ZERODEGC;
  const auto oabi           = sp(1.0)/abi;

  const auto qitot_incld_not_small = qitot_incld >= QSMALL;

  //Compute deposition/sublimation
  qidep.set(qitot_incld_not_small,epsi * oabi * (qv - qvi));

  //Split into deposition or sublimation.
  //if "t" is greater than 0 degree celcius and qidep is positive
  const auto t_gt_zerodegc_pos_qidep = (t < ZERODEGC && qidep > 0);

  qisub.set(qitot_incld_not_small && t_gt_zerodegc_pos_qidep, 0);

  //make qisub positive for consistency with other evap/sub processes
  qisub.set(qitot_incld_not_small && !t_gt_zerodegc_pos_qidep, -pack::min(qidep,0));
  qidep.set(qitot_incld_not_small && !t_gt_zerodegc_pos_qidep, 0);

  //sublimation occurs @ any T. Not so for berg.
  const auto t_lt_zerodegc = t < ZERODEGC;

  //Compute bergeron rate assuming cloud for whole step.
  qiberg.set(qitot_incld_not_small && t_lt_zerodegc, pack::max(epsi*oabi*(qvs - qvi), 0));
  qiberg.set(qitot_incld_not_small && !t_lt_zerodegc, 0);

  nisub.set(qitot_incld_not_small, qisub*(nitot_incld/qitot_incld));

  //if qitot_incld is small (i.e. !qitot_incld_not_small is true)
  qiberg.set(!qitot_incld_not_small, 0);
  qidep.set(!qitot_incld_not_small, 0);
  qisub.set(!qitot_incld_not_small, 0);
  nisub.set(!qitot_incld_not_small, 0);
}


} // namespace p3
} // namespace scream

#endif
