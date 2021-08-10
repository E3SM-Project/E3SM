#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

#include <catch2/catch.hpp>

#include <iostream>

namespace scream {

#ifdef SCREAM_FORCE_RUN_FAIL
TEST_CASE("force_fail", "[fake_infra_test]")
{
  REQUIRE(false); // force this test to fail
}
#endif

#ifdef SCREAM_FORCE_RUN_FPE_FAIL
TEST_CASE("force_fpe", "[fake_infra_test]")
{
  float foo = 42.0;
  float bar = foo / 0.0;
  std::cout << bar << std::endl;

  REQUIRE(true);
}
#endif

#ifdef SCREAM_FORCE_RUN_KOKKOS_OOB
TEST_CASE("force_kokkos_oob", "[fake_infra_test]")
{
  using namespace ekat;
  using KT = KokkosTypes<HostDevice>;

  using view_2d = typename KT::template view_2d<double>;

  view_2d v2d("", 10, 10);
  auto subv = ekat::subview(v2d, 5);

  std::cout << subv.extent(0) << std::endl;
  std::cout << subv(11) << std::endl;

  REQUIRE(true);
}
#endif

TEST_CASE("pass", "[fake_infra_test]")
{
  REQUIRE(true);
}

} // empty namespace
