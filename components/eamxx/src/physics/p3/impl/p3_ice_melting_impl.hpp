#ifndef P3_ICE_MELTING_IMPL_HPP
#define P3_ICE_MELTING_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs
#include "physics/share/physics_saturation_impl.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_melting(
  const Spack& rho, const Spack& T_atm, const Spack& pres, const Spack& rhofaci,
  const Spack& table_val_qi2qr_melting, const Spack& table_val_qi2qr_vent_melt, const Spack& latent_heat_vapor, const Spack& latent_heat_fusion,
  const Spack& dv, const Spack& sc, const Spack& mu, const Spack& kap,
  const Spack& qv, const Spack& qi_incld, const Spack& ni_incld,
  Spack& qi2qr_melt_tend, Spack& ni2nr_melt_tend, const Smask& context)
{
  // Notes Left over from WRF Version:
  // need to add back accelerated melting due to collection of ice mass by rain (pracsw1)
  // note 'f1pr' values are normalized, so we need to multiply by N
  // currently enhanced melting from collision is neglected
  // include RH dependence

  using physics = scream::physics::Functions<Scalar, Device>;

  const auto Pi     = C::Pi;
  const auto QSMALL = C::QSMALL;
  const auto Tmelt  = C::Tmelt;

  //Find cells above freezing AND which have ice
  const auto has_melt_qi = (qi_incld >= QSMALL ) && (T_atm > Tmelt) && context;

  if (has_melt_qi.any()) {
    //    Note that qsat0 should be with respect to liquid. Confirmed F90 code did this.
    const auto qsat0 = physics::qv_sat(Spack(Tmelt), pres, false, context, physics::MurphyKoop, "p3::ice_melting"); //"false" here means NOT saturation w/ respect to ice.

    qi2qr_melt_tend.set(has_melt_qi, ( (table_val_qi2qr_melting+table_val_qi2qr_vent_melt*cbrt(sc)*sqrt(rhofaci*rho/mu))
			     *((T_atm-Tmelt)*kap-rho*latent_heat_vapor*dv*(qsat0-qv))
			     * 2 * Pi /latent_heat_fusion)*ni_incld );

    //make sure qi2qr_melt_tend is always negative
    qi2qr_melt_tend = max(qi2qr_melt_tend, 0);

    //Reduce nj in proportion to decrease in qi mass. Prev line makes sure it always has the right sign.
    ni2nr_melt_tend.set(has_melt_qi, qi2qr_melt_tend*(ni_incld/qi_incld) );
  }
}

} // namespace p3
} // namespace scream

#endif
