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
struct UnitWrap::UnitTest<D>::TestShocPdfTildetoReal {

  static void run_property()
  {
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_tilde_to_real

    // TEST ONE
    // If variance of vertical velocity is zero then
    //  verify that vertical velocity is equal grid mean

    // Define the grid mean vertical velocity [m/s]
    static constexpr Real w_first = 1;
    // Define the standard deviation of vertical velocity [m/s]
    static constexpr Real sqrtw2 = 0;
    // Define the normalized input of vertical velocity
    Real w1 = 0.1;

    // Initialize data structure for bridging to F90
    ShocAssumedPdfTildeToRealData SDS;

    // Fill the test data
    SDS.w_first = w_first;
    SDS.sqrtw2 = sqrtw2;
    SDS.w1 = w1;

    // Call the fortran implementation
    shoc_assumed_pdf_tilde_to_real(SDS);

    // Check the test, verify that vertical velocity is equal
    //  to the grid mean value

    REQUIRE(SDS.w1 == SDS.w_first);

    // TEST TWO
    // Given a series of tests with increasing values of standard
    //  deviation, verify that the gaussian value of w1 also
    //  is gradually increasing

    // Use value from test one
    SDS.w_first = w_first;
    // Initialize with value from test one
    SDS.sqrtw2 = sqrtw2;
    // Define an increment to incrase standard deviation [m/s]
    static constexpr Real incr = 0.1;
    // Define number of tests we want to do
    static constexpr Int num_tests = 10;

    // Require at least two tests
    REQUIRE(num_tests >= 2);

    // Initialize to a large neg value
    Real w1_previous = -999;

    for (Int s = 0; s < num_tests; ++s){
      // must initialize w1 at every test,
      //  since it is input/output
      SDS.w1 = w1;
      // increase value of sqrtw2 at every test
      SDS.sqrtw2 = SDS.sqrtw2 + incr;

      REQUIRE(SDS.sqrtw2 >= 0.0);

      // Call the fortran implementation
      shoc_assumed_pdf_tilde_to_real(SDS);

      // Make sure test value is greater than
      //  previous iteration
      REQUIRE(SDS.w1 > w1_previous);

      // Save the result of this sample
      w1_previous = SDS.w1;

    }
  }

  static void run_bfb()
  {
    // TODO
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_pdf_tildetoreal_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfTildetoReal;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_tildetoreal_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfTildetoReal;

  TestStruct::run_bfb();
}

} // namespace
