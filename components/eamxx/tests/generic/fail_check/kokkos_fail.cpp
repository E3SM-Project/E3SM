#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

#include <catch2/catch.hpp>

#include <iostream>

namespace scream {

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

} // empty namespace
