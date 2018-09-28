#include "catch2/catch.hpp"

#include "share/util/scream_utils.hpp"
#include "share/scream_pack.hpp"
#include "share/scream_kokkos_meta.hpp"

namespace {

TEST_CASE("precision", "util") {
  int nerr = 0;
  if (scream::util::is_single_precision<double>::value) nerr++;
  if ( ! scream::util::is_single_precision<float>::value) nerr++;
  REQUIRE(nerr == 0);
}

// This is just a compilation test.
TEST_CASE("Unmanaged", "scream::ko") {
  using scream::ko::Unmanaged;

  {
    typedef Kokkos::View<double*> V;
    V v("v", 10);
    typedef Unmanaged<V>::type VUm;
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
    typedef Unmanaged<V>::type VUm;
    static_assert(VUm::traits::memory_traits::RandomAccess, "Um");
    static_assert(VUm::traits::memory_traits::Unmanaged, "Um");
    VUm v_um(v);
    typedef Unmanaged<VUm>::type VUmUm;
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
    typedef Unmanaged<V>::type VUm;
    static_assert(VUm::traits::memory_traits::Atomic, "Um");
    static_assert(VUm::traits::memory_traits::Aligned, "Um");
    static_assert(VUm::traits::memory_traits::Restrict, "Um");
    static_assert(VUm::traits::memory_traits::Unmanaged, "Um");
    VUm v_um(v);
  }
}

} // empty namespace
