#ifndef SHOC_SHOC_ASSUMED_PDF_COMPUTE_CLOUD_LIQUID_VARIANCE_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_COMPUTE_CLOUD_LIQUID_VARIANCE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_compute_cloud_liquid_variance. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_compute_cloud_liquid_variance(
  const Spack& a,
  const Spack& s1,
  const Spack& ql1,
  const Spack& C1,
  const Spack& std_s1,
  const Spack& s2,
  const Spack& ql2,
  const Spack& C2,
  const Spack& std_s2,
  const Spack& shoc_ql,
  Spack&       shoc_ql2)
{
  shoc_ql2 = ekat::max(0, a*(s1*ql1 + C1*ekat::square(std_s1))
                          + (1 - a)*(s2*ql2 + C2*ekat::square(std_s2))
                          - ekat::square(shoc_ql));
}

} // namespace shoc
} // namespace scream

#endif
