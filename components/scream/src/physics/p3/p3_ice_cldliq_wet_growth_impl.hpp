#ifndef P3_ICE_CLDLIQ_WET_GROWTH_IMPL_HPP
#define P3_ICE_CLDLIQ_WET_GROWTH_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "physics_functions.hpp" // also for ETI not on GPU
#include "physics_saturation_impl.hpp"

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_cldliq_wet_growth(
  const Spack& rho, const Spack& temp, const Spack& pres, const Spack& rhofaci, const Spack& f1pr05,
  const Spack& f1pr14, const Spack& xxlv, const Spack& xlf, const Spack& dv,
  const Spack& kap, const Spack& mu, const Spack& sc, const Spack& qv, const Spack& qc_incld,
  const Spack& qitot_incld, const Spack& nitot_incld, const Spack& qr_incld,
  Smask& log_wetgrowth, Spack& qrcol, Spack& qccol, Spack& qwgrth, Spack& nrshdr, Spack& qcshd,
  const Smask& context)
{
  using physics = scream::physics::Functions<Scalar, Device>;

  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar tmelt  = C::Tmelt;
  constexpr Scalar twopi  = C::Pi*2;
  constexpr Scalar zero   = C::ZERO;
  constexpr Scalar one    = C::ONE;
  constexpr Scalar cpw    = C::CpLiq;

  const auto t_is_negative = temp < tmelt;
  const auto qitot_incld_ge_small = qitot_incld >= qsmall;
  const auto qc_qr_incld_ge_small = (qc_incld + qr_incld) >= sp(1.0e-6);
  const auto qccol_qrcol_ge_small = (qccol + qrcol) >= sp(1.0e-10);

  const auto any_if     = qitot_incld_ge_small && qc_qr_incld_ge_small && t_is_negative && context;
  const auto any_if_col = any_if && qccol_qrcol_ge_small && context;

  const Spack zerodeg{tmelt};

  Spack qsat0{0.};
  Spack dum{0.};
  Spack dum1{0.};

  if (any_if.any()) {
    qsat0 = physics::qv_sat( zerodeg,pres,0 );

    qwgrth.set(any_if,
               ((f1pr05+f1pr14*pack::cbrt(sc)*sqrt(rhofaci*rho/mu))*
                twopi*(rho*xxlv*dv*(qsat0-qv)-(temp-tmelt)*kap)/
                (xlf+cpw*(temp-tmelt)))*nitot_incld);

    qwgrth.set(any_if,
               pack::max(qwgrth, zero));

    dum = pack::max(zero, (qccol+qrcol)-qwgrth);

    auto const dum_ge_small = dum >= sp(1.0e-10) && context;

    if (dum_ge_small.any()) {
      nrshdr.set(any_if && dum_ge_small,
                 nrshdr+dum*sp(1.923e+6));

      dum1 = one/(qccol+qrcol);

      qcshd.set(any_if_col && dum_ge_small,
                qcshd+dum*qccol*dum1);

      qccol.set(any_if_col && dum_ge_small,
                qccol-dum*qccol*dum1);

      qrcol.set(any_if_col && dum_ge_small,
                qrcol-dum*qrcol*dum1);
    }

    log_wetgrowth = any_if && dum_ge_small;
  }
}

} // namespace p3
} // namespace scream

#endif
