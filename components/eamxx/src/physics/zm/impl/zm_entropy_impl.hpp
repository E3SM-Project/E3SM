#ifndef ZM_ENTROPY_IMPL_HPP
#define ZM_ENTROPY_IMPL_HPP

#include "zm_functions.hpp" // for ETI only but harmless for GPU

#include <ekat_math_utils.hpp>

namespace scream {
namespace zm {

/*
 * Implementation of zm entropy. Clients should NOT
 * #include this file, but include zm_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
Real Functions<S,D>::entropy(
  // Inputs
  const MemberType& team,
  const Real& tk,
  const Real& p,
  const Real& qtot)
{
  //----------------------------------------------------------------------------
  // Purpose: function to calculate entropy following:
  //
  //    Raymond, D. J., and A. M. Blyth, 1992: Extension of the Stochastic Mixing
  //       Model to Cumulonimbus Clouds. J. Atmos. Sci., 49, 1968–1983
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Local variables
  Real qv;    // water vapor mixing ratio
  Real qst;   // saturated vapor mixing ratio
  Real e;     // water vapor pressure
  Real est;   // saturated vapor pressure
  Real L;     // latent heat of vaporization
  //----------------------------------------------------------------------------

  // Calculate latent heat of vaporization - note T is converted to centigrade
  L = PC::LatVap.value - (PC::CpLiq.value - ZMC::cpwv)*(tk - PC::Tmelt.value);

  // Use saturation mixing ratio to partition qtot into vapor part only
  qsat_hPa(tk, p, est, qst);
  qv = ekat::impl::min(qtot,qst);
  e = qv*p / (PC::ep_2.value + qv);

  std::cout << "JGF tk=" << tk << ", p=" << p << ", qtot=" << qtot << ", L=" << L << ", est=" << est << ", qst=" << qst << ", qv=" << qv << ", e=" << e << std::endl;

  // calculate entropy per unit mass of dry air - Eq. 1
  return (
    (PC::Cpair.value + qtot*PC::CpLiq.value)*std::log(tk/PC::Tmelt.value)
    - PC::Rair.value*std::log( (p-e)/ZMC::pref)
    + L*qv/tk - qv*PC::RH2O.value*std::log(qv/qst));
}

} // namespace zm
} // namespace scream

#endif
