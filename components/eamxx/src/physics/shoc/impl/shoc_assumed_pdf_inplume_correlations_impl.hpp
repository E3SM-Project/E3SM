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
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_inplume_correlations(
  const Pack& sqrtqw2_1,
  const Pack& sqrtthl2_1,
  const Pack& a,
  const Pack& sqrtqw2_2,
  const Pack& sqrtthl2_2,
  const Pack& qwthlsec,
  const Pack& qw1_1,
  const Pack& qw_first,
  const Pack& thl1_1,
  const Pack& thl_first,
  const Pack& qw1_2,
  const Pack& thl1_2,
  Pack&       r_qwthl_1)
{
  r_qwthl_1 = 0;

  const Pack testvar = a*sqrtqw2_1*sqrtthl2_1 + (1 - a)*sqrtqw2_2*sqrtthl2_2;
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
