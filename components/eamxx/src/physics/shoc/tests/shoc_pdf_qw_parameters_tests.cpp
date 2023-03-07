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
struct UnitWrap::UnitTest<D>::TestShocQwParameters {

  static void run_property()
  {
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_qw_parameters

    // TEST ONE
    // run the test two times given identical inputs, except
    // one test where wqwsec is larger.  Verify that differece in
    // gaussians is also larger

    // Define the flux of total water [kg/kg m/s] small value
    static constexpr Real wqwsec_small = 1e-4;
    // Define the flux of total water [kg/kg m/s] large value
    static constexpr Real wqwsec_large = 5e-4;
    // Define the standard deviation of vertical velocity [m/s]
    static constexpr Real sqrtw2_test1 = 0.5;
    // Define the variance of total water [kg^2/kg^2]
    static constexpr Real qwsec_test1 = 1e-4;
    // Define grid mean total water [kg/kg]
    static constexpr Real qw_first_test1 = 1.5e-2;
    // Define first gaussian vertical velocity [m/s]
    static constexpr Real w1_1_test1 = 1.5;
    // Define second gaussian vertical velocity [m/s]
    static constexpr Real w1_2_test1 = -0.1;
    // Define vertical velocity skewness [-]
    static constexpr Real Skew_w_test1 = 1;
    // Define fraction of first gaussian
    static constexpr Real a_test1 = 0.2;

    // conversion factor from kg/kg to g/kg
    static constexpr Real qvconv=1000;

    // Define reasonable bounds checking for output
    static constexpr Real qw_bound_low = 1e-3; // [kg/kg]
    static constexpr Real qw_bound_high = 5e-2; // [kg/kg]
    static constexpr Real qw2_bound_high = 1e-2; // [kg^2/kg^2]

    // Be sure this value is actually larger
    REQUIRE(wqwsec_large > wqwsec_small);

    // Initialize data structure for bridging to F90
    ShocAssumedPdfQwParametersData SDS;

    // Fill the test data
    SDS.wqwsec = wqwsec_small;
    SDS.sqrtw2 = sqrtw2_test1;
    SDS.qwsec = qwsec_test1;
    SDS.sqrtqt = sqrt(qwsec_test1);
    SDS.qw_first = qw_first_test1;
    SDS.w1_1 = w1_1_test1;
    SDS.w1_2 = w1_2_test1;
    SDS.skew_w = Skew_w_test1;
    SDS.a = a_test1;

    // Verify input is physical
    REQUIRE(SDS.sqrtw2 >= 0);
    REQUIRE(SDS.sqrtqt >= 0);
    REQUIRE(SDS.qw_first > 0);

    // Make sure vertical velocity data is consistent
    REQUIRE(SDS.w1_1 > 0);
    REQUIRE(SDS.w1_2 < 0);
    if (SDS.skew_w > 0){
      REQUIRE(abs(SDS.w1_1) > abs(SDS.w1_2));
      REQUIRE(SDS.a < 0.5);
    }
    else if (SDS.skew_w < 0){
      REQUIRE(abs(SDS.w1_1) < abs(SDS.w1_2));
      REQUIRE(SDS.a > 0.5);
    }
    else if (SDS.skew_w == 0){
      REQUIRE(abs(SDS.w1_1) == abs(SDS.w1_2));
      REQUIRE(SDS.a == 0);
    }

    // Call Fortran implementation
    shoc_assumed_pdf_qw_parameters(SDS);

    // Save absolute difference between the two gaussian moistures
    Real qwgaus_diff_result1 = abs(qvconv*SDS.qw1_2 - qvconv*SDS.qw1_1);

    // Now laod up value for the large wqwsec test
    SDS.wqwsec = wqwsec_large;

    // Call Fortran implementation
    shoc_assumed_pdf_qw_parameters(SDS);

    // Save absolute difference between the two gaussian temps
    Real qwgaus_diff_result2 = abs(qvconv*SDS.qw1_2 - qvconv*SDS.qw1_1);

    // Now check the result
    REQUIRE(qwgaus_diff_result2 > qwgaus_diff_result1);

    // Make sure output is within reasonable bounds
    REQUIRE( (SDS.qw1_1 > qw_bound_low && SDS.qw1_2 > qw_bound_low) );
    REQUIRE( (SDS.qw1_1 < qw_bound_high && SDS.qw1_2 < qw_bound_high) );
    REQUIRE( (SDS.qw2_1 > 0 && SDS.qw2_2 > 0) );
    REQUIRE( (SDS.qw2_2 < qw2_bound_high && SDS.qw2_2 < qw2_bound_high) );
    REQUIRE( (SDS.sqrtqw2_1 > 0 && SDS.sqrtqw2_2 > 0) );
    REQUIRE( (SDS.sqrtqw2_1 < std::sqrt(qw2_bound_high) && SDS.sqrtqw2_2 < std::sqrt(qw2_bound_high)) );

    // TEST TWO
    // Run the test two times given idential inputs, except one test
    //  where qwsec is larger.  Verify that qw2_1 and qw2_2 are larger

    // Define the variance of total water [kg^2/kg^2] small
    static constexpr Real qwsec_small = 1e-4;
    // Define the variance of total water [kg^2/kg^2] large
    static constexpr Real qwsec_large = 2e-4;

    REQUIRE(qwsec_large > qwsec_small);

    // Fill in test data
    SDS.qwsec = qwsec_small;
    SDS.sqrtqt = sqrt(qwsec_small);

    // Call the Fortran implementaiton
    shoc_assumed_pdf_qw_parameters(SDS);

    // Save the result to compare with next test
    Real qw2_1_result1 = SDS.qw2_1;
    Real qw2_2_result1 = SDS.qw2_2;

    // Now load up input, with larger variances
    SDS.qwsec = qwsec_large;
    SDS.sqrtqt = sqrt(qwsec_large);

    // Call the fortran implementation
    shoc_assumed_pdf_qw_parameters(SDS);

    // Now check the result
    REQUIRE(SDS.qw2_1 > qw2_1_result1);
    REQUIRE(SDS.qw2_2 > qw2_2_result1);
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

TEST_CASE("shoc_pdf_qw_parameters_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocQwParameters;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_qw_parameters_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocQwParameters;

  TestStruct::run_bfb();
}

} // namespace
