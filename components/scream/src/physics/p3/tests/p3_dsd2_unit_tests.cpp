#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

/*
 * Unit-tests for p3 ice table functions.
 */

template <typename D>
struct UnitWrap::UnitTest<D>::TestDsd2 {

  static void run_cloud_bfb()
  {
  }

  static void run_cloud_phys()
  {
  }

  static void run_rain_bfb()
  {
  }

  static void run_rain_phys()
  {
  }
};

}
}
}

namespace {

TEST_CASE("p3_cloud_dsd2", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDsd2;

  TD::run_cloud_phys();
  TD::run_cloud_bfb();
}

TEST_CASE("p3_rain_dsd2", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDsd2;

  TD::run_rain_phys();
  TD::run_rain_bfb();
}

}
