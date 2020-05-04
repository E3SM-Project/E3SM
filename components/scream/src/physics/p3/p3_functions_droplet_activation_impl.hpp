#ifndef P3_FUNCTIONS_DROPLET_ACTIVATION_IMPL_HPP
#define P3_FUNCTIONS_DROPLET_ACTIVATION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "p3_functions_math_impl.hpp"
#include <iomanip>      // std::setprecision


namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
:: droplet_activation(const Spack& temp, const Spack& pres, const Spack& qv, const Spack& qc,
                      const Spack& inv_rho,const Spack& sup, const Spack& xxlv, const Spack& npccn,
                      const bool& log_predictNc, const Scalar& odt,
                      Spack& qcnuc, Spack& ncnuc)
{
  constexpr Scalar cons7  = C::CONS7;
  constexpr Scalar nccnst = C::NCCNST;
  constexpr Scalar rv     = C::RH2O;
  constexpr Scalar inv_cp = C::INV_CP;

  const auto sup_gt_small   = sup > sp(1.e-6);
  const auto sup_if_log     = sup_gt_small && log_predictNc;
  const auto sup_if_not_log = sup_gt_small && !log_predictNc;

  Spack dum{0}, dumqvs{0}, dqsdt{0}, ab{0};

  //.................................................................
  // droplet activation

  ncnuc.set(sup_if_log, npccn);
  qcnuc.set(sup_if_log, ncnuc*cons7);

  // for specified Nc, make sure droplets are present if conditions are supersaturated
  // this is not applied at the first time step, since saturation adjustment is applied at the first step

  dum    = nccnst*inv_rho*cons7-qc;
  dum    = max(0, dum);
  dumqvs = qv_sat(temp, pres, false);
  dqsdt  = xxlv*dumqvs/(rv*temp*temp);
  ab     = 1 + dqsdt*xxlv*inv_cp;
  dum    = min(dum,(qv-dumqvs)/ab);  // limit overdepletion of supersaturation
  qcnuc.set(sup_if_not_log, dum*odt);
}

} // namespace p3
} // namespace scream

#endif
