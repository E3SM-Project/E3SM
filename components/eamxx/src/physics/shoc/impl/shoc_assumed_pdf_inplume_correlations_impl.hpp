#ifndef SHOC_SHOC_ASSUMED_PDF_INPLUME_CORRELATIONS_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_INPLUME_CORRELATIONS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_inplume_correlations. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 *
 * Find within-plume correlation
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_inplume_correlations(
  const Spack& sqrtqw2_1,
  const Spack& sqrtthl2_1,
  const Spack& a,
  const Spack& sqrtqw2_2,
  const Spack& sqrtthl2_2,
  const Spack& qwthlsec,
  const Spack& qw1_1,
  const Spack& qw_first,
  const Spack& thl1_1,
  const Spack& thl_first,
  const Spack& qw1_2,
  const Spack& thl1_2,
  Spack&       r_qwthl_1)
{
  r_qwthl_1 = 0;

  const Spack testvar = a*sqrtqw2_1*sqrtthl2_1 + (1 - a)*sqrtqw2_2*sqrtthl2_2;
  const auto testvar_ne_zero = testvar != 0;
  if (testvar_ne_zero.any()) {
    r_qwthl_1.set(testvar_ne_zero,
                  ekat::max(-1,
                            ekat::min(1, (qwthlsec - a*(qw1_1 - qw_first)*(thl1_1 - thl_first)
                                          - (1 - a)*(qw1_2 - qw_first)*(thl1_2 - thl_first))/testvar)));
  }
}

} // namespace shoc
} // namespace scream

#endif
