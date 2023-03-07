#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestW3diagThirdMoms {

  static void run_property()
  {

    // Tests for the SHOC function:
    //   w3_diag_third_shoc_moment

    // TEST ONE
    // Do series of tests to make sure output is as expected

    // aa0 term
    constexpr static Real aa0_test1 = -0.2;
    // aa1 term
    constexpr static Real aa1_test1 = 5.65;
    // x0 term
    constexpr static Real x0_test1 = -4.31;
    // y0 term
    constexpr static Real x1_test1 = 41.05;
    // f5 term
    constexpr static Real f5_test1 = 4;

    // Initialize data structure for bridging to F90
    W3DiagThirdShocMomentData SDS;

    // Load up the data
    SDS.aa0 = aa0_test1;
    SDS.aa1 = aa1_test1;
    SDS.x0 = x0_test1;
    SDS.x1 = x1_test1;
    SDS.f5 = f5_test1;

    // Call the fortran implementation
    w3_diag_third_shoc_moment(SDS);

    // Verify result is negative
    REQUIRE(SDS.w3 < 0);

    // TEST TWO
    // Modify parameters to decrease w3
    // decrease this term
    constexpr static Real aa1_test2 = 2.65;

    SDS.aa1 = aa1_test2;

    // Call the fortran implementation
    w3_diag_third_shoc_moment(SDS);

    // Verify result has decreased
    REQUIRE(SDS.w3 < SDS.aa1);

    // TEST THREE
    // Modify parameters to get positive result
    // x0 term
    constexpr static Real x0_test3 = -4.31;
    // y0 term
    constexpr static Real x1_test3 = -41.05;
    // f5 term
    constexpr static Real f5_test3 = -4;

    SDS.x0 = x0_test3;
    SDS.x1 = x1_test3;
    SDS.f5 = f5_test3;

    // Call the fortran implementation
    w3_diag_third_shoc_moment(SDS);

    // Verify the result is positive
    REQUIRE(SDS.w3 > 0);

  }

  static void run_bfb()
  {
    // TODO
  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace{

TEST_CASE("shoc_w3_diag_third_moms_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestW3diagThirdMoms;

  TestStruct::run_property();
}

TEST_CASE("shoc_w3_diag_third_moms_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestW3diagThirdMoms;

  TestStruct::run_bfb();
}

} // namespace
