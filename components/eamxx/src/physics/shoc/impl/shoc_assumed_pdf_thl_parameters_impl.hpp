#ifndef SHOC_SHOC_ASSUMED_PDF_THL_PARAMETERS_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_THL_PARAMETERS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc_assumed_pdf_thl_parameters. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 *
 * Find parameters for thetal
 */

template<typename S, typename D>
KOKKOS_INLINE_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_thl_parameters(
  const Pack& wthlsec,
  const Pack& sqrtw2,
  const Pack& sqrtthl,
  const Pack& thlsec,
  const Pack& thl_first,
  const Pack& w1_1,
  const Pack& w1_2,
  const Pack& Skew_w,
  const Pack& a,
  const Scalar thl_tol,
  const Scalar w_thresh,
  Pack&       thl1_1,
  Pack&       thl1_2,
  Pack&       thl2_1,
  Pack&       thl2_2,
  Pack&       sqrtthl2_1,
  Pack&       sqrtthl2_2)
{
  thl1_1 = thl_first;
  thl1_2 = thl_first;
  thl2_1 = 0;
  thl2_2 = 0;
  sqrtthl2_1 = 0;
  sqrtthl2_2 = 0;

  const Mask condition =  thlsec > (thl_tol*thl_tol) && ekat::abs(w1_2 - w1_1) > w_thresh;

  const Pack corrtest1 = ekat::max(-1, ekat::min(1, wthlsec/(sqrtw2*sqrtthl)));
  const Pack tmp_val_1(-corrtest1/w1_1), tmp_val_2(-corrtest1/w1_2);

  Pack Skew_thl(0);
  if (SC::dothetal_skew == true) {
    const auto tsign = ekat::abs(tmp_val_1 - tmp_val_2);
    Skew_thl.set(tsign>sp(0.4), sp(1.2)*Skew_w);
    Skew_thl.set(tsign>sp(0.2) && tsign<=sp(0.4), (((sp(1.2)*Skew_w)/sp(0.2))*(tsign-sp(0.2))));
  }

  if (condition.any()) {
    thl2_1.set(condition,
                ekat::min(100,
                          ekat::max(0, (3*tmp_val_1*(1 - a*ekat::square(tmp_val_2) - (1-a)*ekat::square(tmp_val_1))
                                        - (Skew_thl - a*ekat::cube(tmp_val_2) - (1 - a)*ekat::cube(tmp_val_1)))
                                    /(3*a*(tmp_val_1 - tmp_val_2))))*thlsec);
    thl2_2.set(condition,
                ekat::min(100,
                        ekat::max(0, (-3*tmp_val_2*(1 - a*ekat::square(tmp_val_2)
                                      - (1 - a)*ekat::square(tmp_val_1))
                                      + (Skew_thl - a*ekat::cube(tmp_val_2) - (1 - a)*ekat::cube(tmp_val_1)))
                                      /(3*(1 - a)*(tmp_val_1 - tmp_val_2))))*thlsec);

    thl1_1.set(condition, tmp_val_2*sqrtthl+thl_first);
    thl1_2.set(condition, tmp_val_1*sqrtthl+thl_first);

    sqrtthl2_1.set(condition, ekat::sqrt(thl2_1));
    sqrtthl2_2.set(condition, ekat::sqrt(thl2_2));
  }
}

} // namespace shoc
} // namespace scream

#endif
