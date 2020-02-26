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

  const auto oabi            = sp(1.0)/abi;
  const auto qitot_incld_not_small = qitot_incld >= QSMALL;

  if( qitot_incld_not_small.any()){
    //Compute deposition/sublimation
    qidep = epsi * oabi * (qv - qvi);
    //Split into deposition or sublimation.

    const auto t_qidep_cond = (t < ZERODEGC && qidep > 0);
    qisub.set(t_qidep_cond,0);

    //make qisub positive for consistency with other evap/sub processes
    qisub.set(!t_qidep_cond, -pack::min(qidep,0));
    qidep.set(!t_qidep_cond,0);

    //sublimation occurs @ any T. Not so for berg.
    const auto t_lt_ZERODEGC = t < ZERODEGC ;
    if (t_lt_ZERODEGC.any()){
      //Compute bergeron rate assuming cloud for whole step.
      qiberg = pack::max(epsi*oabi*(qvs - qvi), 0);
    }
    else{ //T>frz
      qiberg=0;
    } //T<frz
    nisub = qisub*(nitot_incld/qitot_incld);
  }
  else{
    qiberg = 0;
    qidep  = 0;
    qisub  = 0;
    nisub  = 0;
  }
}


} // namespace p3
} // namespace scream

#endif
