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
struct UnitWrap::UnitTest<D>::TestOmegadiagThirdMoms {

  static void run_property()
  {

    // Tests for the SHOC function:
    //   omega_terms_diag_third_shoc_moment

    // TEST ONE
    // Call the function twice, with two different sets of terms.
    //  In the second test the fterms should be larger.  Verify
    //  that the omega0 and omega1 terms are the same, but that the
    //  omega2 term has increased.

    // buoyancy term (isotropy squared * brunt vaisalla frequency)
    constexpr static Real buoy_sgs2 = -100;
    // f3 term
    constexpr static Real f3_test1a = 14;
    // f4 term
    constexpr static Real f4_test1a = 5;

    // Initialize data structure for bridging to F90
    OmegaTermsDiagThirdShocMomentData SDS;

    // Load up the data
    SDS.buoy_sgs2 = buoy_sgs2;
    SDS.f3 = f3_test1a;
    SDS.f4 = f4_test1a;

    // Call the fortran implementation
    omega_terms_diag_third_shoc_moment(SDS);

    // Save test results
    Real omega0_test1a = SDS.omega0;
    Real omega1_test1a = SDS.omega1;
    Real omega2_test1a = SDS.omega2;

    // Now load up data for second part of test
    // Feed in LARGER values
    SDS.f3 = 2*f3_test1a;
    SDS.f4 = 2*f4_test1a;

    // Call the fortran implementation
    omega_terms_diag_third_shoc_moment(SDS);

    // Now check the result

    // omega0 and omega1 should NOT have changed
    REQUIRE(SDS.omega0 == omega0_test1a);
    REQUIRE(SDS.omega1 == omega1_test1a);

    // omega2 should have increased
    REQUIRE(SDS.omega2 > omega2_test1a);

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

TEST_CASE("shoc_omega_diag_third_moms_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestOmegadiagThirdMoms;

  TestStruct::run_property();
}

TEST_CASE("shoc_omega_diag_third_moms_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestOmegadiagThirdMoms;

  TestStruct::run_bfb();
}

} // namespace
