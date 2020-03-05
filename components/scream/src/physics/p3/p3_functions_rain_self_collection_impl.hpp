#ifndef P3_FUNCTIONS_RAIN_SELF_COLLECTION_IMPL_HPP
#define P3_FUNCTIONS_RAIN_SELF_COLLECTION_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_self_collection(const Spack& rho, const Spack& qr_incld, const Spack& nr_incld, Spack& nrslf)
{

  // ------------------------------------------------------
  // self-collection and breakup of rain
  // (breakup following modified Verlinde and Cotton scheme)

  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar rhow   = C::RhoH2O;
  constexpr Scalar pi     = C::Pi;

  const auto qr_incld_not_small = qr_incld >= qsmall; 
  
  if (qr_incld_not_small.any()){
    const Real dum1 = sp(280.e-6); 
    const auto dum2 = pack::cbrt((qr_incld)/(pi*rhow*nr_incld)); 
    const auto dum2_lt_dum1 = dum2 < dum1;
    const auto dum3 =sp(5.78)*nr_incld*qr_incld*rho; 
    if (dum2_lt_dum1.any()){
      nrslf.set(qr_incld_not_small && dum2_lt_dum1, dum3);//sp(5.78)*nr_incld*qr_incld*rho);
    }
    const auto dum2_ge_dum1 = dum2 >= dum1;  
    if (dum2_ge_dum1.any()){
      const auto dum = 2-pack::exp(2300*(dum2-dum1));
      nrslf.set(qr_incld_not_small && dum2_ge_dum1, dum*dum3);
    }
  }
}

} // namespace p3
} // namespace scream

#endif
