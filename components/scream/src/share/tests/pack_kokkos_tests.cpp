#include "catch2/catch.hpp"

#include "share/scream_pack_kokkos.hpp"
#include "share/scream_config.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"

namespace {

TEST_CASE("index", "scream::pack") {
}

template <typename View, int rank>
using OnlyRank = typename std::enable_if<View::Rank == rank>::type;

TEST_CASE("scalarize", "scream::pack") {
  using scream::pack::Pack;
  using scream::pack::scalarize;

  typedef Kokkos::View<Pack<double, 16>*> Array1;
  typedef Kokkos::View<Pack<double, 32>**> Array2;

  {
    const Array1 a1("a1", 10);
    const auto a2 = scalarize(a1);
    static_assert(decltype(a2)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a2.extent_int(0) == 160);
  }  

  {
    const Array2 a1("a1", 10, 4);
    const auto a2 = scalarize(a1);
    static_assert(decltype(a2)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a2.extent_int(0) == 10);
    REQUIRE(a2.extent_int(1) == 128);
  }
}

template <typename View>
OnlyRank<View, 1> fill (const View& a) {
  const auto m = Kokkos::create_mirror_view(a);
  for (int i = 0; i < m.extent_int(0); ++i)
    m(i) = i;
  Kokkos::deep_copy(a, m);
}

template <typename View>
OnlyRank<View, 2> fill (const View& a) {
  const auto m = Kokkos::create_mirror_view(a);
  for (int i = 0; i < m.extent_int(0); ++i)
    for (int j = 0; j < m.extent_int(1); ++j)
      a(i,j) = j*m.extent_int(0) + i;
  Kokkos::deep_copy(a, m);
}

template <typename VA, typename VB>
OnlyRank<VA, 1> compare (const VA& a, const VB& b) {
  const auto ma = Kokkos::create_mirror_view(a);
  const auto mb = Kokkos::create_mirror_view(b);
  for (int i = 0; i < ma.extent_int(0); ++i)
    REQUIRE(ma(i) == mb(i));
}

template <typename VA, typename VB>
OnlyRank<VA, 2> compare (const VA& a, const VB& b) {
  const auto ma = Kokkos::create_mirror_view(a);
  const auto mb = Kokkos::create_mirror_view(b);
  for (int i = 0; i < ma.extent_int(0); ++i)
    for (int j = 0; j < ma.extent_int(1); ++j)
      REQUIRE(ma(i,j) == mb(i,j));
}

template <int repack_size, typename Src, typename Dst>
OnlyRank<Src, 1> repack_test (const Src& a_src, const Dst& a) {
  static_assert(Dst::traits::memory_traits::Unmanaged, "Um");
  static_assert(Dst::value_type::n == repack_size, "Pack::n");
  REQUIRE(a.extent_int(0) == (a_src(0).n/repack_size)*a_src.extent_int(0));
  compare(scalarize(a_src), scalarize(a));
}

template <int repack_size, typename Src, typename Dst>
OnlyRank<Src, 2> repack_test (const Src& a_src, const Dst& a) {
  static_assert(Dst::traits::memory_traits::Unmanaged, "Um");
  static_assert(Dst::value_type::n == repack_size, "Pack::n");
  REQUIRE(a.extent_int(0) == a_src.extent_int(0));
  REQUIRE(a.extent_int(1) == (a_src(0,0).n/repack_size)*a_src.extent_int(1));
  compare(scalarize(a_src), scalarize(a));
}

TEST_CASE("repack", "scream::pack") {
  using scream::pack::Pack;
  using scream::pack::repack;

  typedef Kokkos::View<Pack<double, 16>*> Array1;
  typedef Kokkos::View<Pack<double, 32>**> Array2;

  {
    const Array1 a1("a1", 10);
    fill(a1);

    const auto a2 = repack<8>(a1);
    repack_test<8>(a1, a2);

    const auto a3 = repack<4>(a2);
    repack_test<4>(a2, a3);

    const auto a4 = repack<2>(a3);
    repack_test<2>(a3, a4);

    const auto a5 = repack<2>(a1);
    repack_test<2>(a1, a5);

    const auto a6 = repack<16>(a1);
    repack_test<16>(a1, a6);
  }

  {
    const Array2 a1("a1", 10, 4);
    fill(a1);

    const auto a2 = repack<8>(a1);
    repack_test<8>(a1, a2);

    const auto a3 = repack<4>(a2);
    repack_test<4>(a2, a3);

    const auto a4 = repack<2>(a3);
    repack_test<2>(a3, a4);

    const auto a5 = repack<2>(a1);
    repack_test<2>(a1, a5);

    const auto a6 = repack<32>(a1);
    repack_test<32>(a1, a6);
  }
}

} // empty namespace
