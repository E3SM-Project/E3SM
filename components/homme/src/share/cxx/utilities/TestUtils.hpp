/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_TEST_UTILS_HPP
#define HOMMEXX_TEST_UTILS_HPP

#include "Types.hpp"
#include <functional>

namespace Homme {

template <typename rngAlg, typename PDF>
void genRandArray(int *const x, int length, rngAlg &engine, PDF &&pdf) {
  for (int i = 0; i < length; ++i) {
    x[i] = pdf(engine);
  }
}

template <typename rngAlg, typename PDF>
void genRandArray(Real *const x, int length, rngAlg &engine, PDF &&pdf) {
  for (int i = 0; i < length; ++i) {
    x[i] = pdf(engine);
  }
}

template <typename rngAlg, typename PDF>
void genRandArray(Scalar *const x, int length, rngAlg &engine, PDF &&pdf) {
  for (int i = 0; i < length; ++i) {
    for (int j = 0; j < VECTOR_SIZE; ++j) {
      x[i][j] = pdf(engine);
    }
  }
}

template <typename ViewType, typename rngAlg, typename PDF>
typename std::enable_if<Kokkos::is_view<ViewType>::value, void>::type
genRandArray(ViewType view, rngAlg &engine, PDF &&pdf,
             std::function<bool(typename ViewType::HostMirror)> constraint) {
  typename ViewType::HostMirror mirror = Kokkos::create_mirror_view(view);
  do {
    genRandArray(mirror.data(), view.size(), engine, pdf);
  } while (constraint(mirror) == false);
  Kokkos::deep_copy(view, mirror);
}

template <typename ViewType, typename rngAlg, typename PDF>
typename std::enable_if<Kokkos::is_view<ViewType>::value, void>::type
genRandArray(ViewType view, rngAlg &engine, PDF &&pdf) {
  genRandArray(view, engine, pdf,
               [](typename ViewType::HostMirror) { return true; });
}

template <typename FPType>
Real compare_answers(FPType target, FPType computed,
                     FPType relative_coeff = 1.0) {
  Real denom = 1.0;
  if (relative_coeff > 0.0 && target != 0.0) {
    denom = relative_coeff * std::fabs(target);
  }

  return std::fabs(target - computed) / denom;
}

} // namespace Homme

#endif // HOMMEXX_TEST_UTILS_HPP
