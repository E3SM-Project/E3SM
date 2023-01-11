#ifndef P3_GET_TIME_SPACE_PHYS_VARIABLES_IMPL_HPP
#define P3_GET_TIME_SPACE_PHYS_VARIABLES_IMPL_HPP

#include "p3_functions.hpp"
#include "physics/share/physics_constants.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::get_time_space_phys_variables(
  const Spack& T_atm, const Spack& pres, const Spack& rho, const Spack& latent_heat_vapor, const Spack& latent_heat_sublim,
  const Spack& qv_sat_l, const Spack& qv_sat_i, Spack& mu, Spack& dv, Spack& sc, Spack& dqsdt,
  Spack& dqsidt, Spack& ab, Spack& abi, Spack& kap, Spack& eii,
  const Smask& context)
{
  //time/space varying physical variables
  mu.set(context, sp(1.496e-6) * pow(T_atm,sp(1.5))/(T_atm+120));
  dv.set(context, sp(8.794e-5) * pow(T_atm,sp(1.81))/pres);
  sc.set(context, mu/(rho*dv));

  constexpr Scalar RV     = C::RV;
  constexpr Scalar INV_CP = C::INV_CP;
  constexpr Scalar tval1  = 253.15;
  constexpr Scalar tval2  = 273.15;
  constexpr Scalar dtval  = 20; //this is tval2-tval1, but specifying here as int to be BFB with F90.

  const auto dum = 1/(RV*square(T_atm));
  dqsdt.set(context, latent_heat_vapor*qv_sat_l*dum);
  dqsidt.set(context, latent_heat_sublim*qv_sat_i*dum);
  ab.set(context, 1+dqsdt*latent_heat_vapor*INV_CP);
  abi.set(context, 1+dqsidt*latent_heat_sublim*INV_CP);
  kap.set(context, sp(1.414e+3)*mu);

  //very simple temperature dependent aggregation efficiency
  const auto t_lt_tval1 = T_atm < tval1;
  const auto t_lt_tval2 = T_atm < tval2;

  eii.set(t_lt_tval1 && context,sp(0.001));
  eii.set(!t_lt_tval1 && t_lt_tval2 && context, sp(0.001)+( T_atm-sp(tval1) )*( sp(0.3) - sp(0.001) )/dtval);
  eii.set(!t_lt_tval1 && !t_lt_tval2 && context, sp(0.3) );
}

} // namespace p3
} // namespace scream

#endif
