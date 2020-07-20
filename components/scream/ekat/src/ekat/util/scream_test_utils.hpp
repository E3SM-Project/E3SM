#ifndef EKAT_TEST_UTILS_HPP
#define EKAT_TEST_UTILS_HPP

#include <ekat/scream_types.hpp>
#include <ekat/scream_pack.hpp>

#include <catch2/catch.hpp>

namespace scream {
namespace util {

#ifdef EKAT_STRICT_FP
constexpr bool StrictFP = true;
#else
constexpr bool StrictFP = false;
#endif

struct TestSession {
  static TestSession& get () {
    static TestSession s;
    return s;
  }

  std::map<std::string,std::string> params;
private:
  TestSession() = default;
};

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

// Do an == check between a scalar result and a packed-C++ result.
// Expect BFB except when C++ pksize > 1 and fp model is not strict
template <bool StrictFP, int Packsize, typename Scalar>
void catch2_req_pk_sensitive(const Scalar lhs, const Scalar rhs)
{
  if (StrictFP) {
    REQUIRE(lhs == rhs);
  } else {
    if (Packsize > 1) {
      REQUIRE(lhs == Approx(rhs));
    } else {
      REQUIRE(lhs == rhs);
    }
  }
}

} // namespace util
} // namespace scream

#endif // EKAT_TEST_UTILS_HPP
