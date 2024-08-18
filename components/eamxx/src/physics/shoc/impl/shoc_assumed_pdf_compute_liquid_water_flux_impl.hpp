#ifndef SHOC_SHOC_ASSUMED_PDF_COMPUTE_LIQUID_WATER_FLUX_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_COMPUTE_LIQUID_WATER_FLUX_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_compute_liquid_water_flux. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_compute_liquid_water_flux(
  const Spack& a,
  const Spack& w1_1,
  const Spack& w_first,
  const Spack& ql1,
  const Spack& w1_2,
  const Spack& ql2,
  Spack&       wqls)
{
  wqls = a*((w1_1 - w_first)*ql1) + (1 - a)*((w1_2 - w_first)*ql2);
}

} // namespace shoc
} // namespace scream

#endif
