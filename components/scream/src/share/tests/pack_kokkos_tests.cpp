#include "catch2/catch.hpp"

#include "share/scream_pack_kokkos.hpp"
#include "share/scream_config.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"

namespace {

TEST_CASE("index", "scream::pack") {
}

TEST_CASE("scalarize", "scream::pack") {
}

TEST_CASE("repack", "scream::pack") {
  using scream::pack::Pack;
  using scream::pack::repack;

  typedef Kokkos::View<Pack<double, 16>*> Array1;
  typedef Kokkos::View<Pack<double, 16>**> Array2;

  {
    const Array1 a1("a1", 10);
    const auto a2 = repack<8>(a1);
    static_assert(decltype(a2)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a2.extent_int(0) == 2*a1.extent_int(1));
  }

  {
    const Array2 a1("a1", 10, 4);
    const auto a2 = repack<8>(a1);
    static_assert(decltype(a2)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a2.extent_int(0) ==   a1.extent_int(0));
    REQUIRE(a2.extent_int(1) == 2*a1.extent_int(1));
  }
}

} // empty namespace
