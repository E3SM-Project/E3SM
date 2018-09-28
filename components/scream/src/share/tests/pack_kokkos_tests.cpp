#include "catch2/catch.hpp"

#include "share/scream_pack_kokkos.hpp"
#include "share/scream_config.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"

namespace {

TEST_CASE("index", "scream::pack") {
}

TEST_CASE("scalarize", "scream::pack") {
  using scream::pack::Pack;
  using scream::pack::scalarize;

  typedef Kokkos::View<Pack<double, 16>*> Array1;
  typedef Kokkos::View<Pack<double, 16>**> Array2;

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
    REQUIRE(a2.extent_int(1) == 64);
  }  
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
    REQUIRE(a2.extent_int(0) == 2*a1.extent_int(0));
    const auto a3 = repack<4>(a2);
    static_assert(decltype(a3)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a3.extent_int(0) == 4*a1.extent_int(0));    
    const auto a4 = repack<2>(a3);
    static_assert(decltype(a4)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a4.extent_int(0) == 8*a1.extent_int(0));    
    const auto a5 = repack<2>(a1);
    static_assert(decltype(a5)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a4.extent_int(0) == 8*a1.extent_int(0));    
    const auto a6 = repack<16>(a1);
    static_assert(decltype(a6)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a6.extent_int(0) == a1.extent_int(0));    
  }

  {
    const Array2 a1("a1", 10, 4);
    const auto a2 = repack<8>(a1);
    static_assert(decltype(a2)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a2.extent_int(0) ==   a1.extent_int(0));
    REQUIRE(a2.extent_int(1) == 2*a1.extent_int(1));
    const auto a3 = repack<4>(a2);
    static_assert(decltype(a3)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a3.extent_int(0) ==   a1.extent_int(0));
    REQUIRE(a3.extent_int(1) == 4*a1.extent_int(1));
    const auto a4 = repack<2>(a3);
    static_assert(decltype(a4)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a4.extent_int(0) ==   a1.extent_int(0));
    REQUIRE(a4.extent_int(1) == 8*a1.extent_int(1));
    const auto a5 = repack<2>(a1);
    static_assert(decltype(a5)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a5.extent_int(0) ==   a1.extent_int(0));
    REQUIRE(a5.extent_int(1) == 8*a1.extent_int(1));
    const auto a6 = repack<16>(a1);
    static_assert(decltype(a6)::traits::memory_traits::Unmanaged, "Um");
    REQUIRE(a6.extent_int(0) ==   a1.extent_int(0));
    REQUIRE(a6.extent_int(1) ==   a1.extent_int(1));
  }
}

} // empty namespace
