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
struct UnitWrap::UnitTest<D>::TestShocPdfCompCldVar {

  static void run_property()
  {
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_compute_cloud_liquid_variance

    // TEST ONE
    // standard deviation test.  Given two tests with identical inputs,
    //  but one where s1 and s2 are larger, verify that SGS cloud water
    //  variance is larger

    // Define gaussian fraction [-]
    static constexpr Real a = 0.3;
    // Define s1
    static constexpr Real s1 = 0.03;
    // Define cloud liquid water for gaussian one [kg/kg]
    static constexpr Real ql1 = 0.01;
    // define cloud fraction for gaussian one [-]
    static constexpr Real C1 = 0.2;
    // Define s2
    static constexpr Real s2 = 0.003;
    // Define cloud liquid water for gaussian two [kg/kg]
    static constexpr Real ql2 = 0.001;
    // define cloud fraction for gaussian one [-]
    static constexpr Real C2 = 0.1;
    // define grid mean shoc_ql [kg/kg]
    static constexpr Real shoc_ql = 0.005;

    // Now define standard deviations of the gaussians.
    //    Here we define two sets for the different tests
    static constexpr Real std_s1_small = 0.01;
    static constexpr Real std_s2_small = 0.001;

    // Define larger quantities for the second part of the test
    static constexpr Real std_s1_large = 0.02;
    static constexpr Real std_s2_large = 0.002;

    // Verify these are in fact larger
    REQUIRE(std_s1_large > std_s1_small);
    REQUIRE(std_s2_large > std_s2_small);

    // Initialize data structure for bridging to F90
    ShocAssumedPdfComputeCloudLiquidVarianceData SDS;

    // load the data for the first part of the test
    SDS.a = a;
    SDS.s1 = s1;
    SDS.ql1 = ql1;
    SDS.c1 = C1;
    SDS.s2 = s2;
    SDS.ql2 = ql2;
    SDS.c2 = C2;
    SDS.shoc_ql = shoc_ql;

    SDS.std_s1 = std_s1_small;
    SDS.std_s2 = std_s2_small;

    // Verify inputs
    REQUIRE(SDS.a > 0);
    REQUIRE(SDS.a < 1);
    REQUIRE(SDS.s1 >= 0);
    REQUIRE(SDS.c1 >= 0);
    REQUIRE(SDS.ql1 >= 0);
    REQUIRE(SDS.s2 >= 0);
    REQUIRE(SDS.c2 >= 0);
    REQUIRE(SDS.ql2 >= 0);
    REQUIRE(SDS.shoc_ql >= 0);
    REQUIRE(SDS.std_s1 >= 0);
    REQUIRE(SDS.std_s2 >= 0);

    // Call the fortran implementation
    shoc_assumed_pdf_compute_cloud_liquid_variance(SDS);

    // Save the result for the first part of the test
    Real shoc_ql2_small = SDS.shoc_ql2;

    // Now load the data needed for the second part of the test
    SDS.std_s1 = std_s1_large;
    SDS.std_s2 = std_s2_large;

    // Call the fortran implementation
    shoc_assumed_pdf_compute_cloud_liquid_variance(SDS);

    // Verify that the liquid water variance is larger when the standard
    //   deviations of s are larger

    REQUIRE(SDS.shoc_ql2 > shoc_ql2_small);

    // TEST TWO
    // Grid mean ql test.  Given identical inputs with the exception of
    //  that in one set of inputs the grid mean of cloud water is larger, verify
    //  that shoc_ql2 is small in this case

    // For inputs we will recycle the inputs from the last test used.
    //  The exception being the grid mean cloud water

    static constexpr Real shoc_ql_small = 0.005;
    static constexpr Real shoc_ql_large = 0.006;

    // Verify that values are expected
    REQUIRE(shoc_ql_large > shoc_ql_small);

    // Load input data
    SDS.shoc_ql = shoc_ql_small;

    // Call the fortran implementation
    shoc_assumed_pdf_compute_cloud_liquid_variance(SDS);

    // save the output for this test
    shoc_ql2_small = SDS.shoc_ql2;

    // Load data for second part of this test
    SDS.shoc_ql = shoc_ql_large;

    // Call the fortran implementation
    shoc_assumed_pdf_compute_cloud_liquid_variance(SDS);

    // Check the result.  The cloud water variance should be SMALLER
    //  if the grid mean is larger for this case
    REQUIRE(SDS.shoc_ql2 < shoc_ql2_small);
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

TEST_CASE("shoc_pdf_compute_cldvar_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfCompCldVar;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_compute_cldvar_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfCompCldVar;

  TestStruct::run_bfb();
}

} // namespace
