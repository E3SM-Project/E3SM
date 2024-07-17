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
KOKKOS_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_compute_temperature(
  const Spack& thl1,
  const Spack& pval,
  Spack&       Tl1)
{
  Tl1 = thl1/(ekat::pow(C::P0/pval,(C::Rair/C::CP)));
}

} // namespace shoc
} // namespace scream

#endif
