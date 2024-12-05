#ifndef SHOC_SHOC_ASSUMED_PDF_QW_PARAMETERS_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_QW_PARAMETERS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_qw_parameters. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 *
 * Find parameters for total water mixing ratio
 */

template<typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_qw_parameters(
  const Spack& wqwsec,
  const Spack& sqrtw2,
  const Spack& Skew_w,
  const Spack& sqrtqt,
  const Spack& qwsec,
  const Spack& w1_2,
  const Spack& w1_1,
  const Spack& qw_first,
  const Spack& a,
  const Scalar rt_tol,
  const Scalar w_thresh,
  Spack&       qw1_1,
  Spack&       qw1_2,
  Spack&       qw2_1,
  Spack&       qw2_2,
  Spack&       sqrtqw2_1,
  Spack&       sqrtqw2_2)
{
  qw1_1 = qw_first;
  qw1_2 = qw_first;
  qw2_1 = 0;
  qw2_2 = 0;
  sqrtqw2_1 = 0;
  sqrtqw2_2 = 0;

  const Smask condition = qwsec > (rt_tol*rt_tol) && ekat::abs(w1_2 - w1_1) > w_thresh;

  const Spack corrtest2 = ekat::max(-1, ekat::min(1, wqwsec/(sqrtw2*sqrtqt)));
  const Spack tmp_val_1(-corrtest2/w1_1), tmp_val_2(-corrtest2/w1_2);

  const auto tsign = ekat::abs(tmp_val_1 - tmp_val_2);
  Spack Skew_qw(0);
  Skew_qw.set(tsign>sp(0.4), sp(1.2)*Skew_w);
  Skew_qw.set(tsign>sp(0.2) && tsign<=sp(0.4), (((sp(1.2)*Skew_w)/sp(0.2))*(tsign-sp(0.2))));

  if (condition.any()) {
    qw2_1.set(condition,
            ekat::min(100,
                      ekat::max(0, (3*tmp_val_1*(1 - a*ekat::square(tmp_val_2) - (1 - a)*ekat::square(tmp_val_1))
                                    - (Skew_qw - a*ekat::cube(tmp_val_2) - (1 - a)*ekat::cube(tmp_val_1)))
                                    /(3*a*(tmp_val_1 - tmp_val_2))))*qwsec);
    qw2_2.set(condition,
            ekat::min(100,
                      ekat::max(0, (-3*tmp_val_2*(1 - a*ekat::square(tmp_val_2) - (1 - a)*ekat::square(tmp_val_1))
                                    + (Skew_qw - a*ekat::cube(tmp_val_2) - (1 - a)*ekat::cube(tmp_val_1)))
                                    /(3*(1 - a)*(tmp_val_1 - tmp_val_2))))*qwsec);
  }

  qw1_1.set(condition, tmp_val_2*sqrtqt+qw_first);
  qw1_2.set(condition, tmp_val_1*sqrtqt+qw_first);

  sqrtqw2_1.set(condition, ekat::sqrt(qw2_1));
  sqrtqw2_2.set(condition, ekat::sqrt(qw2_2));
}

} // namespace shoc
} // namespace scream

#endif
