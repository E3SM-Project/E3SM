#ifndef SCREAM_TEST_UTILS_HPP
#define SCREAM_TEST_UTILS_HPP

#include <share/scream_types.hpp>
#include <share/scream_pack.hpp>
#include <share/field/field.hpp>

namespace scream
{

namespace util
{

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

template <typename rngAlg, typename PDF, typename ScalarType, int N>
void genRandArray(pack::Pack<ScalarType,N> *const x, int length, rngAlg &engine, PDF &&pdf) {
  for (int i = 0; i < length; ++i) {
    for (int j = 0; j < N; ++j) {
      x[i][j] = pdf(engine);
    }
  }
}

template <typename ViewType, typename rngAlg, typename PDF>
typename std::enable_if<Kokkos::is_view<ViewType>::value, void>::type
genRandArray(ViewType view, rngAlg &engine, PDF &&pdf) {
  typename ViewType::HostMirror mirror = Kokkos::create_mirror_view(view);
  genRandArray(mirror.data(), view.size(), engine, pdf);
  Kokkos::deep_copy(view, mirror);
}

template <typename FieldType, typename rngAlg, typename PDF>
typename std::enable_if<is_scream_field<FieldType>::value, void>::type
genRandArray(FieldType& field, rngAlg &engine, PDF &&pdf) {
  error::runtime_check(field.is_allocated(), "Error! Field must be allocated before it can be filled.\n");
  genRandArray(field.get_view(), engine, pdf);
}

} // namespace util

} // namespace scream

#endif // SCREAM_TEST_UTILS_HPP
