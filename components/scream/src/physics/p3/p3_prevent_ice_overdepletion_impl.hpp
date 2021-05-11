#ifndef P3_PREVENT_ICE_OVERDEPLETION_IMPL_HPP
#define P3_PREVENT_ICE_OVERDEPLETION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "physics_functions.hpp" // also for ETI when not GPU
#include "physics_saturation_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 ice overdepletion prevention.
 * Clients should NOT #include this file, but include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::prevent_ice_overdepletion(
  const Spack& pres, const Spack& T_atm, const Spack& qv, const Spack& latent_heat_sublim, const Scalar& inv_dt,
  Spack& qv2qi_vapdep_tend, Spack& qi2qv_sublim_tend, const Smask& range_mask,
  const Smask& context)
{
  using physics = scream::physics::Functions<Scalar, Device>;

  constexpr Scalar cp = C::CP;
  constexpr Scalar rv = C::RH2O;

  const auto dumqv_sat_i = physics::qv_sat(T_atm,pres,true, range_mask);
  const auto qdep_satadj = (qv-dumqv_sat_i) /
    (1 + square(latent_heat_sublim) * dumqv_sat_i / (cp * rv * square(T_atm))) * inv_dt;
  qv2qi_vapdep_tend.set(context, qv2qi_vapdep_tend * min(1, max(0,  qdep_satadj) / max(qv2qi_vapdep_tend, sp(1.e-20))));
  qi2qv_sublim_tend.set(context, qi2qv_sublim_tend * min(1, max(0, -qdep_satadj) / max(qi2qv_sublim_tend, sp(1.e-20))));
}

} // namespace p3
} // namespace scream

#endif
