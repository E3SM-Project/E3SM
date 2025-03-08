#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/eamxx_types.hpp"

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
struct UnitWrap::UnitTest<D>::TestShocVVParameters {

  static void run_property()
  {
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_vv_parameters

    // TEST ONE
    // Given a symmetric distribution, verify that gaussian fractions
    //  are equal and that the two gaussians have identical vertical
    //  velocities (but opposite sign).

    // Define the grid mean vertical velocity [m/s]
    static constexpr Real w_first_sym = 1;
    // Define the standard deviation of vertical velocity [m2/s2]
    static constexpr Real w_sec_sym = 1;
    // Define the third moment of vertical velocity [m3/s3]
    static constexpr Real w3var_sym = 0;

    // Initialize data structure for bridging to F90
    ShocAssumedPdfVvParametersData SDS;

    // Fill the test data
    SDS.w_first = w_first_sym;
    SDS.w_sec = w_sec_sym;
    SDS.w3var = w3var_sym;
    SDS.w_tol_sqd = 0;

    // Verify input is physical
    REQUIRE(SDS.w_sec >= 0);

    // Call the fortran implementation
    shoc_assumed_pdf_vv_parameters(SDS);

    // Check the results

    // Verify that gaussian fractions are equal
    REQUIRE(SDS.a == 0.5);
    // Verify that guassian vertical velocities are
    //   identical but opposite sign
    REQUIRE(SDS.w1_1 == -1.0*SDS.w1_2);

    // TEST TWO
    // Given highly positive skewed distribution, verify that the gaussian
    //  fractions and vertical velocities are as expected

    // Define the grid mean vertical velocity [m/s]
    static constexpr Real w_first_skew = 0;
    // Define the standard deviation of vertical velocity [m2/s2]
    static constexpr Real w_sec_skew = 1;
    // Define the third moment of vertical velocity [m3/s3]
    static constexpr Real w3var_skew = 5;

    // Fill the test data
    SDS.w_first = w_first_skew;
    SDS.w_sec = w_sec_skew;
    SDS.w3var = w3var_skew;

    // Verify input is physicsl
    REQUIRE(SDS.w_sec >= 0);
    // For this test we want w_first to be zero exactly
    REQUIRE(SDS.w_first == 0);
    // For this test we want w3 to be positive
    REQUIRE(SDS.w3var > 0);

    // Call the fortran implementation
    shoc_assumed_pdf_vv_parameters(SDS);

    // Verify that first gaussian is less than half
    REQUIRE(SDS.a < 0.5);

    // Verify that first gaussian absolute value is larger
    //  than second gaussian
    REQUIRE(abs(SDS.w1_1) > abs(SDS.w1_2));

    // TEST THREE
    // Given highly negative skewed distribution, verify that the gaussian
    //  fractions and vertical velocities are as expected

    // Define the grid mean vertical velocity [m/s]
    static constexpr Real w_first_neg = 0;
    // Define the standard deviation of vertical velocity [m2/s2]
    static constexpr Real w_sec_neg = 1;
    // Define the third moment of vertical velocity [m3/s3]
    static constexpr Real w3var_neg = -5;

    // Fill the test data
    SDS.w_first = w_first_neg;
    SDS.w_sec = w_sec_neg;
    SDS.w3var = w3var_neg;

    // Verify input is physicsl
    REQUIRE(SDS.w_sec >= 0);
    // For this test we want w_first to be zero exactly
    REQUIRE(SDS.w_first == 0);
    // For this test we want w3 to be negative
    REQUIRE(SDS.w3var < 0);

    // Call the fortran implementation
    shoc_assumed_pdf_vv_parameters(SDS);

    // Verify that first gaussian is greater than half
    REQUIRE(SDS.a > 0.5);

    // Verify that second gaussian absolute value is larger
    //  than first gaussian
    REQUIRE(abs(SDS.w1_1) < abs(SDS.w1_2));
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

TEST_CASE("shoc_pdf_vv_parameters_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocVVParameters;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_vv_parameters_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocVVParameters;

  TestStruct::run_bfb();
}

} // namespace
