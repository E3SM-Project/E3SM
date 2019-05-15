#include <catch2/catch.hpp>

#include "share/scream_pack_kokkos.hpp"
#include "share/scream_config.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/util/scream_kokkos_utils.hpp"

namespace {

template <typename View, int rank>
using OnlyRank = typename std::enable_if<View::Rank == rank>::type;

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
      m(i,j) = j*m.extent_int(0) + i;
  Kokkos::deep_copy(a, m);
}

template <typename VA, typename VB>
OnlyRank<VA, 1> compare (const VA& a, const VB& b) {
  const auto ma = Kokkos::create_mirror_view(a);
  const auto mb = Kokkos::create_mirror_view(b);
  Kokkos::deep_copy(ma, a);
  Kokkos::deep_copy(mb, b);
  for (int i = 0; i < ma.extent_int(0); ++i)
    REQUIRE(ma(i) == mb(i));
}

template <typename VA, typename VB>
OnlyRank<VA, 2> compare (const VA& a, const VB& b) {
  const auto ma = Kokkos::create_mirror_view(a);
  const auto mb = Kokkos::create_mirror_view(b);
  Kokkos::deep_copy(ma, a);
  Kokkos::deep_copy(mb, b);
  for (int i = 0; i < ma.extent_int(0); ++i)
    for (int j = 0; j < ma.extent_int(1); ++j)
      REQUIRE(ma(i,j) == mb(i,j));
}

TEST_CASE("index", "scream::pack") {
  using scream::pack::scalarize;
  using scream::pack::Pack;

  {
    static constexpr int pack_size = 16;
    using IdxPack = Pack<int, pack_size>;
    Kokkos::View<double*> data("data", 100);
    fill(data);
    IdxPack idx;
    for (int i = 0; i < pack_size; ++i) idx[i] = 2*i;
    int nerr = 0;
    Kokkos::parallel_reduce(
      1, KOKKOS_LAMBDA (const int /* unused */, int& nerr) {
        const auto data_idx = index(data, idx);
        for (int i = 0; i < pack_size; ++i)
          if (data_idx[i] != idx[i])
            ++nerr;
      },
      nerr);
    REQUIRE(nerr == 0);
  }

  {
    static constexpr int pack_size = 8;
    using IdxPack = Pack<int, pack_size>;
    Kokkos::View<double**> data("data", 19, 24);
    fill(data);
    IdxPack i0, i1;
    for (int i = 0; i < pack_size; ++i) i0[i] = 2*i;
    for (int i = 0; i < pack_size; ++i) i1[i] = 3*i;
    int nerr = 0;
    Kokkos::parallel_reduce(
      1, KOKKOS_LAMBDA (const int /* unused */, int& nerr) {
        const auto data_idx = index(data, i0, i1);
        for (int i = 0; i < pack_size; ++i)
          if (data_idx[i] != data(i0[i], i1[i]))
            ++nerr;
      },
      nerr);
    REQUIRE(nerr == 0);
  }
}

TEST_CASE("scalarize", "scream::pack") {
  using scream::pack::Pack;
  using scream::pack::scalarize;

  typedef Kokkos::View<Pack<double, 16>*> Array1;
  typedef Kokkos::View<Pack<double, 32>**> Array2;
  typedef Kokkos::View<Pack<double, 8>***> Array3;
  typedef Kokkos::View<Pack<double, 24>****> Array4;

  {
    const Array1 a1("a1", 10);
    const auto a2 = scalarize(a1);
    typedef decltype(a2) VT;
    static_assert(VT::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a2.extent_int(0) == 160);
  }

  {
    const Array2 a1("a2", 10, 4);
    const auto a2 = scalarize(a1);
    typedef decltype(a2) VT;
    static_assert(VT::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a2.extent_int(0) == 10);
    REQUIRE(a2.extent_int(1) == 128);
  }

  {
    const Array3 a1("a3", 3, 2, 4);
    const auto a2 = scalarize(a1);
    typedef decltype(a2) VT;
    static_assert(VT::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a2.extent_int(0) == 3);
    REQUIRE(a2.extent_int(1) == 2);
    REQUIRE(a2.extent_int(2) == 32);
  }

  {
    const Array4 a1("a4", 3, 2, 4, 2);
    const auto a2 = scalarize(a1);
    typedef decltype(a2) VT;
    static_assert(VT::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a2.extent_int(0) == 3);
    REQUIRE(a2.extent_int(1) == 2);
    REQUIRE(a2.extent_int(2) == 4);
    REQUIRE(a2.extent_int(3) == 48);
  }
}

template <int repack_size, typename Src, typename Dst>
OnlyRank<Src, 1> repack_test (const Src& a_src, const Dst& a) {
  static_assert(Dst::traits::memory_traits::Unmanaged, "Um");
  static_assert(Dst::value_type::n == repack_size, "Pack::n");
  REQUIRE(a.extent_int(0) == (Src::value_type::n/repack_size)*a_src.extent_int(0));
  compare(scalarize(a_src), scalarize(a));
}

template <int repack_size, typename Src, typename Dst>
OnlyRank<Src, 2> repack_test (const Src& a_src, const Dst& a) {
  static_assert(Dst::traits::memory_traits::Unmanaged, "Um");
  static_assert(Dst::value_type::n == repack_size, "Pack::n");
  REQUIRE(a.extent_int(0) == a_src.extent_int(0));
  REQUIRE(a.extent_int(1) == (Src::value_type::n/repack_size)*a_src.extent_int(1));
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

TEST_CASE("kokkos_packs", "scream::pack") {
  using namespace scream;
  using namespace scream::pack;

  using TestBigPack = Pack<Real, 16>;

  using ExeSpace = typename KokkosTypes<DefaultDevice>::ExeSpace;
  using MemberType = typename KokkosTypes<DefaultDevice>::MemberType;

  int nerr = 0;
  const int num_bigs = 17;

  typename KokkosTypes<DefaultDevice>::template view_1d<TestBigPack> test_k_array("test_k_array", num_bigs);
  Kokkos::parallel_reduce("unittest_pack",
                          util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, 1),
                          KOKKOS_LAMBDA(const MemberType& team, int& total_errs) {

    int nerrs_local = 0;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_bigs), [&] (int i) {
      test_k_array(i) = i;
    });

    auto small = repack<4>(test_k_array);
    if (small.extent(0) != 4 * num_bigs) ++nerrs_local;

    team.team_barrier();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_bigs*4), [&] (int i) {
      for (int p = 0; p < 4; ++p) {
        if (small(i)[p] != i / 4) ++nerrs_local;
      }
    });

    auto big = repack<16>(small);
    if (big.extent(0) != num_bigs) ++nerrs_local;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_bigs*4), [&] (int i) {
      for (int p = 0; p < 4; ++p) {
        small(i)[p] = p * i;
      }
    });

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_bigs*4), [&] (int i) {
      auto mask = small(i) >= (2 * i);
      for (int p = 0; p < 4; ++p) {
        if (i == 0) {
          if (!mask[p]) ++nerrs_local;
        }
        else {
          if (mask[p] != (p >= 2)) ++nerrs_local;
        }
      }
    });

    total_errs += nerrs_local;
  }, nerr);

  // NOTE: catch2 documentation says that its assertion macros are not
  // thread safe, so we have to put them outside of kokkos kernels.
  REQUIRE(nerr == 0);
}

} // namespace
