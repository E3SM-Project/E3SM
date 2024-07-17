#ifndef SHOC_SHOC_ASSUMED_PDF_COMPUTE_S_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_COMPUTE_S_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_compute_s. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_compute_s(
  const Spack& qw1,
  const Spack& qs,
  const Spack& beta,
  const Spack& pval,
  const Spack& thl2,
  const Spack& qw2,
  const Spack& sqrtthl2,
  const Spack& sqrtqw2,
  const Spack& r_qwthl,
  Spack&       s,
  Spack&       std_s,
  Spack&       qn,
  Spack&       C)
{
  const Scalar rair = C::Rair;
  const Scalar basepres = C::P0;
  const Scalar cp = C::CP;
  const Scalar lcond = C::LatVap;
  const Scalar pi = C::Pi;

  const Scalar sqrt2(std::sqrt(Scalar(2.0))), sqrt2pi(std::sqrt(2*pi));

  const Spack cthl=((1 + beta*qw1)/ekat::square(1 + beta*qs))*(cp/lcond)*
                    beta*qs*ekat::pow(pval/basepres, (rair/cp));
  const Spack cqt = 1/(1 + beta*qs);

  std_s = ekat::sqrt(ekat::max(0,
                               ekat::square(cthl)*thl2
                               + ekat::square(cqt)*qw2 - 2*cthl*sqrtthl2*cqt*sqrtqw2*r_qwthl));
  const auto std_s_not_small = std_s > std::sqrt(std::numeric_limits<Scalar>::min()) * 100;
  s = qw1-qs*((1 + beta*qw1)/(1 + beta*qs));
  if (std_s_not_small.any()) {
    C.set(std_s_not_small, sp(0.5)*(1 + ekat::erf(s/(sqrt2*std_s))));
  }
  C.set(!std_s_not_small && s > 0, 1);
  const auto std_s_C_not_small = std_s_not_small && C != 0;
  if (std_s_C_not_small.any()) {
    qn.set(std_s_C_not_small, s*C+(std_s/sqrt2pi)*ekat::exp(-sp(0.5)*ekat::square(s/std_s)));
  }
  qn.set(!std_s_not_small && s > 0, s);

  // Checking to prevent empty clouds
  const auto qn_le_zero = qn <= 0;
  C.set(qn_le_zero,0);
  qn.set(qn_le_zero,0);
}

} // namespace shoc
} // namespace scream

#endif
