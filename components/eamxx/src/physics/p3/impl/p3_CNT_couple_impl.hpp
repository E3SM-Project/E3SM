#ifndef P3_CNT_COUPLE_IMPL_HPP
#define P3_CNT_COUPLE_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::CNT_couple(
  const Spack& frzimm, const Spack& frzcnt,
  const Spack& frzdep, const Spack& rho,
  const Spack& qc_incld, const Spack& nc_incld,
  const int Iflag,
  Spack& ncheti_cnt, Spack& qcheti_cnt,
  Spack& nicnt, Spack& qicnt,
  Spack& ninuc_cnt, Spack& qinuc_cnt)
{
  constexpr Scalar pi       = C::Pi;
  constexpr Scalar rho_h2o  = C::RHO_H2O;
  constexpr Scalar qsmall   = 1.0e-18; // BAD_CONSTANT!
  constexpr Scalar piov3    = pi/3.0;
  constexpr Scalar  mi0     = 4.0*piov3*900.0*1.0e-18; // BAD_CONSTANT!
 
  const Spack Zero(0.0);
  // minimum mass of new crystal due to freezing of cloud droplets done
  // externally (kg)
 
  const Scalar mi0l_min = (4.0/3.0)*pi*rho_h2o*(4.0e-6)*(4.0e-6)*(4.0e-6);
  Spack mi0l = qc_incld/ekat::max(nc_incld,1.0e6/rho);
  mi0l = ekat::max(mi0l_min, mi0l);
  
  const auto mask = qc_incld > qsmall;
  switch (Iflag) {
    case 1:  // cloud droplet immersion freezing
      ncheti_cnt.set(mask, frzimm*1.0e6/rho /* frzimm input is in [#/cm3] */ , Zero);
      qcheti_cnt.set(mask, ncheti_cnt*mi0l, Zero);
      break;
    case 2:  // deposition freezing / contact freezing
      nicnt.set(mask, frzcnt*1.0e6/rho, Zero);
      qicnt.set(mask, nicnt*mi0l, Zero);
      ninuc_cnt.set(mask, frzdep*1.0e6/rho, Zero);
      qinuc_cnt.set(mask, ninuc_cnt*mi0, Zero);
      break;
    default:
      break;
  }
}
} // namespace p3
} // namespace scream

#endif
