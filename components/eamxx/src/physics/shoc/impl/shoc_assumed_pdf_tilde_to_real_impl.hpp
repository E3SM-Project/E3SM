#ifndef SHOC_SHOC_ASSUMED_PDF_TILDE_TO_REAL_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_TILDE_TO_REAL_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_tilde_to_real. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 *
 * Convert tilde variables to "real" variables
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_tilde_to_real(
  const Spack& w_first,
  const Spack& sqrtw2,
  Spack&       w1)
{
  w1 *= sqrtw2;
  w1 += w_first;
}

} // namespace shoc
} // namespace scream

#endif
