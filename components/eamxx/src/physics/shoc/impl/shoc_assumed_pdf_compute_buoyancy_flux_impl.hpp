#ifndef SHOC_SHOC_ASSUMED_PDF_COMPUTE_BUOYANCY_FLUX_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_COMPUTE_BUOYANCY_FLUX_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_compute_buoyancy_flux. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_compute_buoyancy_flux(
  const Spack& wthlsec,
  const Spack& wqwsec,
  const Spack& pval,
  const Spack& wqls,
  Spack&       wthv_sec)
{
  const Scalar basepres = C::P0.value;
  const Scalar rair     = C::Rair.value;
  const Scalar rv       = C::RV.value;
  const Scalar cp       = C::CP.value;
  const Scalar lcond    = C::LatVap.value;
  const Scalar basetemp = C::basetemp;
  const Scalar epsterm  = rair/rv;

  wthv_sec = wthlsec + ((1 - epsterm)/epsterm)*basetemp*wqwsec
             + ((lcond/cp)*ekat::pow(basepres/pval, (rair/cp))
             - (1/epsterm)*basetemp)*wqls;
}

} // namespace shoc
} // namespace scream

#endif
