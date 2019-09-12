#include <catch2/catch.hpp>

#include "share/util/scream_utils.hpp"
#include "share/scream_pack.hpp"
#include "share/scream_kokkos_meta.hpp"

#include <cmath>

namespace {

TEST_CASE("precision", "util") {
  CHECK_FALSE(scream::util::is_single_precision<double>::value);
  CHECK(scream::util::is_single_precision<float>::value);
}

TEST_CASE("ipow_scalar", "util") {
  for (int i = 0; i < 5; ++i) {
    REQUIRE(scream::util::ipow_scalar(2, i) == std::pow(2, i));
  }
  double base = 1.167692662218503654003143310546875E+06;
  double invbase = 1/base;
  REQUIRE(scream::util::ipow_scalar(base, 3) == base*base*base);
  REQUIRE(scream::util::ipow_scalar(base, -3) == scream::util::ipow_scalar(invbase, 3));

  // This is case that caused problems
  REQUIRE(scream::util::ipow_scalar(base, 3) != std::pow(base, 3));
}

// This is just a compilation test.
TEST_CASE("Unmanaged", "scream::ko") {
  using scream::ko::Unmanaged;

  {
    typedef Kokkos::View<double*> V;
    V v("v", 10);
    typedef Unmanaged<V> VUm;
    VUm v_um(v);
    static_assert( ! V::traits::memory_traits::Unmanaged, "Um");
    static_assert(VUm::traits::memory_traits::Unmanaged, "Um");
  }

  {
    typedef Kokkos::View<scream::pack::Pack<double, SCREAM_PACK_SIZE>***,
                         Kokkos::LayoutLeft,
                         Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::RandomAccess> >
      V;
    V v("v", 2, 3, 4);
    typedef Unmanaged<V> VUm;
    static_assert(VUm::traits::memory_traits::RandomAccess, "Um");
    static_assert(VUm::traits::memory_traits::Unmanaged, "Um");
    VUm v_um(v);
    typedef Unmanaged<VUm> VUmUm;
    static_assert(VUmUm::traits::memory_traits::RandomAccess, "Um");
    static_assert(VUmUm::traits::memory_traits::Unmanaged, "Um");
    static_assert( ! VUmUm::traits::memory_traits::Atomic, "Um");
    static_assert( ! VUmUm::traits::memory_traits::Aligned, "Um");
    static_assert( ! VUmUm::traits::memory_traits::Restrict, "Um");
    VUmUm v_umum(v);
  }

  {
    typedef Kokkos::View<scream::pack::Pack<int, SCREAM_PACK_SIZE>[10],
                         Kokkos::HostSpace,
                         Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::Aligned | Kokkos::Restrict> >
      V;
    static_assert( ! V::traits::memory_traits::Unmanaged, "Um");
    V v("v");
    typedef Unmanaged<V>::const_type CVUm;
    static_assert(CVUm::traits::memory_traits::Atomic, "Um");
    static_assert(CVUm::traits::memory_traits::Aligned, "Um");
    static_assert(CVUm::traits::memory_traits::Restrict, "Um");
    static_assert(CVUm::traits::memory_traits::Unmanaged, "Um");

    using Kokkos::Impl::ViewMapping;
    static_assert(ViewMapping<CVUm::traits, V::traits, void>::is_assignable,
                  "CVUm <- V");
    static_assert( ! ViewMapping<V::traits, CVUm::traits, void>::is_assignable,
                  "V </- CVUm");
    static_assert(ViewMapping<CVUm::traits, Unmanaged<V>::traits, void>::is_assignable,
                  "CVUm <- VUm");
    static_assert( ! ViewMapping<Unmanaged<V>::traits, CVUm::traits, void>::is_assignable,
                  "VUm </- CVUm");
    CVUm cv_um(v);
  }
}

} // empty namespace
