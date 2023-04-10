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
struct UnitWrap::UnitTest<D>::TestAAdiagThirdMoms {

  static void run_property()
  {

    // Tests for the SHOC function:
    //   aa_terms_diag_third_shoc_moment

    // Run test two times.  One where 0 parameters stay the same
    //  and the other where 1 parameters vary.  Verify that the aa0
    //  term remains constant between the two tests.

    // omega0 term
    constexpr static Real omega0 = 7.54e-2;
    // omega1 term
    constexpr static Real omega1 = 5.37e-3;
    // omega2 term
    constexpr static Real omega2 = 0.54;
    // x0 term
    constexpr static Real x0 = -4.31;
    // y0 term
    constexpr static Real y0 = 22.72;
    // x1 term
    constexpr static Real x1_test1a = 41.05;
    // y1 term
    constexpr static Real y1_test1a = 375.69;

    // Initialize data structure for bridging to F90
    AaTermsDiagThirdShocMomentData SDS;

    // Load up the data
    SDS.omega0 = omega0;
    SDS.x0 = x0;
    SDS.y0 = y0;
    SDS.omega1 = omega1;
    SDS.x1 = x1_test1a;
    SDS.y1 = y1_test1a;
    SDS.omega2 = omega2;

    // Call the fortran implementation
    aa_terms_diag_third_shoc_moment(SDS);

    // Save the output
    Real aa0_test1a = SDS.aa0;
    Real aa1_test1a = SDS.aa1;

    // Load up data where only the 1 terms have varies
    SDS.x1 = 0.3*x1_test1a;
    SDS.y1 = 0.3*y1_test1a;

    // Call the fortran implementation
    aa_terms_diag_third_shoc_moment(SDS);

    // Check the result
    REQUIRE(SDS.aa0 == aa0_test1a);
    REQUIRE(SDS.aa1 < aa1_test1a);

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

TEST_CASE("shoc_aa_diag_third_moms_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestAAdiagThirdMoms;

  TestStruct::run_property();
}

TEST_CASE("shoc_aa_diag_third_moms_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestAAdiagThirdMoms;

  TestStruct::run_bfb();

}

} // namespace
