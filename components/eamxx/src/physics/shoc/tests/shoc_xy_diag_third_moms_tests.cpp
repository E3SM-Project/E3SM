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
struct UnitWrap::UnitTest<D>::TestXYdiagThirdMoms {

  static void run_property()
  {

    // Tests for the SHOC function:
    //   x_y_terms_diag_third_shoc_moment

    // TEST ONE
    // Call the function twice, with two different sets of terms.
    //  In the second test the fterms should be larger.  Verify
    //  that the x0 and y0 terms are the same, but that the
    //  x1 and y1 terms have increased.

    // buoyancy term (isotropy squared * brunt vaisalla frequency)
    constexpr static Real buoy_sgs2 = -100;
    // f0 term
    constexpr static Real f0_test1a = 13000;
    // f3 term
    constexpr static Real f1_test1a = 8500;
    // f4 term
    constexpr static Real f2_test1a = 30;

    // Initialize data structure for bridging to F90
    XYTermsDiagThirdShocMomentData SDS;

    // Load up the data
    SDS.buoy_sgs2 = buoy_sgs2;
    SDS.f0 = f0_test1a;
    SDS.f1 = f1_test1a;
    SDS.f2 = f2_test1a;

    // Call the fortran implementation
    x_y_terms_diag_third_shoc_moment(SDS);

    // Save test results
    Real x0_test1a = SDS.x0;
    Real y0_test1a = SDS.y0;
    Real x1_test1a = SDS.x1;
    Real y1_test1a = SDS.y1;

    // Now load up data for second part of test
    // Feed in LARGER values
    SDS.f0 = 1.2*f0_test1a;
    SDS.f1 = 1.2*f1_test1a;
    SDS.f2 = 1.2*f2_test1a;

    // Call the fortran implementation
    x_y_terms_diag_third_shoc_moment(SDS);

    // Now check the result

    // x0 and y0 terms should NOT have changed
    REQUIRE(SDS.x0 == x0_test1a);
    REQUIRE(SDS.y0 == y0_test1a);

    // x1 and y1 terms should have increased
    REQUIRE(SDS.x1 > x1_test1a);
    REQUIRE(SDS.y1 > y1_test1a);

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

TEST_CASE("shoc_xy_diag_third_moms_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestXYdiagThirdMoms;

  TestStruct::run_property();
}

TEST_CASE("shoc_xy_diag_third_moms_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestXYdiagThirdMoms;

  TestStruct::run_bfb();
}

} // namespace
