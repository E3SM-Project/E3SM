#ifndef P3_ICE_CLDLIQ_WET_GROWTH_IMPL_HPP
#define P3_ICE_CLDLIQ_WET_GROWTH_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp" // also for ETI not on GPU
#include "physics/share/physics_saturation_impl.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_cldliq_wet_growth(
  const Spack& rho, const Spack& temp, const Spack& pres, const Spack& rhofaci, const Spack& table_val_qi2qr_melting,
  const Spack& table_val_qi2qr_vent_melt, const Spack& latent_heat_vapor, const Spack& latent_heat_fusion, const Spack& dv,
  const Spack& kap, const Spack& mu, const Spack& sc, const Spack& qv, const Spack& qc_incld,
  const Spack& qi_incld, const Spack& ni_incld, const Spack& qr_incld,
  Smask& log_wetgrowth, Spack& qr2qi_collect_tend, Spack& qc2qi_collect_tend, Spack& qc_growth_rate, Spack& nr_ice_shed_tend, Spack& qc2qr_ice_shed_tend, const Smask& context)
{
  using physics = scream::physics::Functions<Scalar, Device>;

  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar tmelt  = C::Tmelt;
  constexpr Scalar twopi  = C::Pi*2;
  constexpr Scalar zero   = C::ZERO;
  constexpr Scalar one    = C::ONE;
  constexpr Scalar cpw    = C::CpLiq;

  const auto t_is_negative = temp < tmelt;
  const auto qi_incld_ge_small = qi_incld >= qsmall;
  const auto qc_qr_incld_ge_small = (qc_incld + qr_incld) >= sp(1.0e-6);
  const auto qc2qi_collect_tend_qr2qi_collect_tend_ge_small = (qc2qi_collect_tend + qr2qi_collect_tend) >= sp(1.0e-10);

  const auto any_if     = qi_incld_ge_small && qc_qr_incld_ge_small && t_is_negative && context;
  const auto any_if_col = any_if && qc2qi_collect_tend_qr2qi_collect_tend_ge_small && context;

  const Spack zerodeg{tmelt};

  Spack qsat0{0.};
  Spack dum{0.};
  Spack dum1{0.};

  if (any_if.any()) {
    qsat0 = physics::qv_sat( zerodeg,pres, false, context, physics::MurphyKoop, "p3::ice_cldliq_wet_growth" );

    qc_growth_rate.set(any_if,
               ((table_val_qi2qr_melting+table_val_qi2qr_vent_melt*cbrt(sc)*sqrt(rhofaci*rho/mu))*
                twopi*(rho*latent_heat_vapor*dv*(qsat0-qv)-(temp-tmelt)*kap)/
                (latent_heat_fusion+cpw*(temp-tmelt)))*ni_incld);

    qc_growth_rate.set(any_if,
               max(qc_growth_rate, zero));

    dum = max(zero, (qc2qi_collect_tend+qr2qi_collect_tend)-qc_growth_rate);

    auto const dum_ge_small = dum >= sp(1.0e-10) && context;

    if (dum_ge_small.any()) {
      nr_ice_shed_tend.set(any_if && dum_ge_small,
                 nr_ice_shed_tend+dum*sp(1.923e+6));

      dum1 = one/(qc2qi_collect_tend+qr2qi_collect_tend);

      qc2qr_ice_shed_tend.set(any_if_col && dum_ge_small,
                qc2qr_ice_shed_tend+dum*qc2qi_collect_tend*dum1);

      qc2qi_collect_tend.set(any_if_col && dum_ge_small,
			     max(0,qc2qi_collect_tend-dum*qc2qi_collect_tend*dum1));

      qr2qi_collect_tend.set(any_if_col && dum_ge_small,
			     max(0,qr2qi_collect_tend-dum*qr2qi_collect_tend*dum1));
    }

    log_wetgrowth = any_if && dum_ge_small;
  }
}

} // namespace p3
} // namespace scream

#endif // P3_ICE_CLDLIQ_WET_GROWTH_IMPL_HPP
