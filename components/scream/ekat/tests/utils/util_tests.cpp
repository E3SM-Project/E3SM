#include <catch2/catch.hpp>

#include "ekat/util/scream_utils.hpp"
#include "ekat/util/string_utils.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/scream_kokkos_meta.hpp"

namespace {

TEST_CASE("precision", "util") {
  CHECK_FALSE(scream::util::is_single_precision<double>::value);
  CHECK(scream::util::is_single_precision<float>::value);
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
    typedef Kokkos::View<scream::pack::Pack<double, EKAT_PACK_SIZE>***,
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
    typedef Kokkos::View<scream::pack::Pack<int, EKAT_PACK_SIZE>[10],
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

TEST_CASE("string","string") {
  using namespace scream;

  util::CaseInsensitiveString cis1 = "field_1";
  util::CaseInsensitiveString cis2 = "fIeLd_1";
  util::CaseInsensitiveString cis3 = "field_2";
  util::CaseInsensitiveString cis4 = "feld_1";

  REQUIRE (cis1==cis2);
  REQUIRE (cis1!=cis3);
  REQUIRE (cis4<=cis1);
  REQUIRE (cis4<cis1);
}

} // empty namespace
