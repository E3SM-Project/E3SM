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
  const Spack& Tl1_1,
  const Spack& Tl1_2,
  const Spack& pval,
  const Smask& active_entries,
  Spack&       qs1,
  Spack&       beta1,
  Spack&       qs2,
  Spack&       beta2)
{
  const Scalar rair = C::Rair;
  const Scalar rv = C::RV;
  const Scalar cp = C::CP;
  const Scalar lcond = C::LatVap;

  // Compute MurphyKoop_svp
  const int liquid = 0;
  const Spack esval1_1 = scream::physics::Functions<S,D>::MurphyKoop_svp(Tl1_1,liquid,active_entries,"shoc::shoc_assumed_pdf (Tl1_1)");
  const Spack esval1_2 = scream::physics::Functions<S,D>::MurphyKoop_svp(Tl1_2,liquid,active_entries,"shoc::shoc_assumed_pdf (Tl1_2)");
  const Spack lstarn(lcond);

  qs1 = sp(0.622)*esval1_1/ekat::max(esval1_1, pval - esval1_1);
  beta1 = (rair/rv)*(lstarn/(rair*Tl1_1))*(lstarn/(cp*Tl1_1));

  // Only compute qs2 and beta2 if the two plumes are not equal
  const Smask condition = (Tl1_1 != Tl1_2);
  qs2 = qs1;
  beta2 = beta1;

  qs2.set(condition, sp(0.622)*esval1_2/ekat::max(esval1_2, pval - esval1_2));
  beta2.set(condition, (rair/rv)*(lstarn/(rair*Tl1_2))*(lstarn/(cp*Tl1_2)));
}

} // namespace shoc
} // namespace scream

#endif
