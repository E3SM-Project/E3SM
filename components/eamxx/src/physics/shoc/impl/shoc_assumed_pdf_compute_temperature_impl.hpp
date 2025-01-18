#ifndef SHOC_SHOC_ASSUMED_PDF_COMPUTE_TEMPERATURE_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_COMPUTE_TEMPERATURE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_compute_temperature. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_compute_temperature(
  const Spack& thl1,
  const Spack& pval,
  Spack&       Tl1)
{
  constexpr Scalar basepres = C::P0;
  constexpr Scalar rair = C::Rair;
  constexpr Scalar cp = C::CP;
  Tl1 = thl1/(ekat::pow(basepres/pval,(rair/cp)));
}

} // namespace shoc
} // namespace scream

#endif
