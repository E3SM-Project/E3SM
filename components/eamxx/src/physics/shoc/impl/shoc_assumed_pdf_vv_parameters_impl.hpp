#ifndef SHOC_SHOC_ASSUMED_PDF_VV_PARAMETERS_IMPL_HPP
#define SHOC_SHOC_ASSUMED_PDF_VV_PARAMETERS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

#include <iomanip>

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_assumed_pdf_vv_parameters. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 *
 * Find parameters for vertical velocity
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_assumed_pdf_vv_parameters(
  const Spack& w_first,
  const Spack& w_sec,
  const Spack& w3var,
  const Scalar w_tol_sqd,
  Spack&       Skew_w,
  Spack&       w1_1,
  Spack&       w1_2,
  Spack&       w2_1,
  Spack&       w2_2,
  Spack&       a)
{
  Skew_w = 0;
  w1_1 = w_first;
  w1_2 = w_first;
  w2_1 = 0;
  w2_2 = 0;
  a = 0.5;

  const Smask condition = w_sec > w_tol_sqd;

  const Scalar tmp_val(0.4);
  const Scalar one_m_tmp_val(1 - tmp_val);
  const Scalar sqrtw2t(std::sqrt(1-tmp_val));

  Skew_w.set(condition, w3var/ekat::sqrt(ekat::cube(w_sec)));
  a.set(condition,
        ekat::max(sp(0.01),
                  ekat::min(sp(0.99),
                            sp(0.5)*(1 - Skew_w*ekat::sqrt(1/(4*(one_m_tmp_val*one_m_tmp_val*one_m_tmp_val)
                                                              + ekat::square(Skew_w)))))));

  w1_1.set(condition, ekat::sqrt((1 - a)/a)*sqrtw2t);
  w1_2.set(condition, -1*ekat::sqrt(a/(1 - a))*sqrtw2t);
  w2_1.set(condition, tmp_val*w_sec);
  w2_2.set(condition, tmp_val*w_sec);
}

} // namespace shoc
} // namespace scream

#endif
