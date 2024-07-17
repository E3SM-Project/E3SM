#ifndef SHOC_SHOC_ASSUMED_PDF_COMPUTE_SGS_LIQUID_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_COMPUTE_SGS_LIQUID_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_compute_sgs_liquid. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_compute_sgs_liquid(
  const Spack& a,
  const Spack& ql1,
  const Spack& ql2,
  Spack&       shoc_ql)
{
  shoc_ql = ekat::max(0, a*ql1 + (1 - a)*ql2);
}

} // namespace shoc
} // namespace scream

#endif
