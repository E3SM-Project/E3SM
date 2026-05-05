#ifndef SHOC_SHOC_ASSUMED_PDF_COMPUTE_QS_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_COMPUTE_QS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU
#include "physics_functions.hpp"

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_compute_qs. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_compute_qs(
  const Pack& Tl1_1,
  const Pack& Tl1_2,
  const Pack& pval,
  const Mask& active_entries,
  Pack&       qs1,
  Pack&       beta1,
  Pack&       qs2,
  Pack&       beta2)
{
  const Scalar rair  = C::Rair.value;
  const Scalar rv    = C::RV.value;
  const Scalar cp    = C::CP.value;
  const Scalar lcond = C::LatVap.value;

  // Compute MurphyKoop_svp
  const int liquid = 0;
  const Pack esval1_1 = scream::physics::Functions<S,D>::MurphyKoop_svp(Tl1_1,liquid,active_entries,"shoc::shoc_assumed_pdf (Tl1_1)");
  const Pack esval1_2 = scream::physics::Functions<S,D>::MurphyKoop_svp(Tl1_2,liquid,active_entries,"shoc::shoc_assumed_pdf (Tl1_2)");
  const Pack lstarn(lcond);

  qs1 = sp(0.622)*esval1_1/ekat::max(esval1_1, pval - esval1_1);
  beta1 = (rair/rv)*(lstarn/(rair*Tl1_1))*(lstarn/(cp*Tl1_1));

  // Only compute qs2 and beta2 if the two plumes are not equal
  const Mask condition = (Tl1_1 != Tl1_2);
  qs2 = qs1;
  beta2 = beta1;

  qs2.set(condition, sp(0.622)*esval1_2/ekat::max(esval1_2, pval - esval1_2));
  beta2.set(condition, (rair/rv)*(lstarn/(rair*Tl1_2))*(lstarn/(cp*Tl1_2)));
}

} // namespace shoc
} // namespace scream

#endif
