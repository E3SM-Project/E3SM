#ifndef P3_FUNCTIONS_GET_TIME_SPACE_PHYS_VARIABLES_IMPL.HPP
#define P3_FUNCTIONS_GET_TIME_SPACE_PHYS_VARIABLES_IMPL.HPP

#include "p3_functions.hpp"
#include "p3_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::get_time_space_phys_variables(const Spack& t, const Spack& pres, const Spack& rho, const Spack& xxlv, const Spack& xxls,
				const Spack& qvs, const Spack& qvi, Spack& mu, Spack& dv, Spack& sc, Spack& dqsdt,
				Spack& dqsidt, Spack& ab, Spack& abi, Spack& kap, Spack& eii)
{
  //time/space varying physical variables
  mu = sp(1.496e-6) * pow(t,sp(1.5))/(t+120);
  dv = sp(8.794e-5) * pow(t,sp(1.81))/pres;
  sc = mu/(rho*dv);

  constexpr Scalar RV     = C::RV;
  constexpr Scalar INV_CP = C::INV_CP;
  constexpr Scalar tval1  = 253.15;
  constexpr Scalar tval2  = 268.15;

  const auto dum = 1/(RV*square(t));
  dqsdt  = xxlv*qvs*dum;
  dqsidt = xxls*qvi*dum;
  ab     = 1+dqsdt*xxlv*INV_CP;
  abi    = 1+dqsidt*xxls*INV_CP;
  kap    = sp(1.414e+3)*mu;

  //very simple temperature dependent aggregation efficiency
  const auto t_lt_tval1 = t < tval1;
  const auto t_lt_tval2 = t < tval2;

  eii.set(t_lt_tval1,sp(0.1));
  eii.set(!t_lt_tval1 && t_lt_tval2,sp(0.1)+(t-sp(253.15))/15*sp(0.9));
  eii.set(!t_lt_tval1 && !t_lt_tval2,1);
}


} // namespace p3
} // namespace scream

#endif
