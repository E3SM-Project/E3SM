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
struct UnitWrap::UnitTest<D>::TestShocPdfComputeQs {

  static void run_property()
  {
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_compute_qs

    // TEST ONE
    // Zero Skewness tests.  Given inputs where Temperatures for
    //  the two gaussians are the same, verify that outputs are the same

    // Define temperature of Gaussian 1 [K]
    static constexpr Real Tl1_1_eq = 290;
    // Define temperature of Gaussian 2 [K]
    static constexpr Real Tl1_2_eq = Tl1_1_eq;
    // Define pressure value [Pa]
    static constexpr Real pval_eq = 85000;

    // Initialize data structure for bridging to F90
    ShocAssumedPdfComputeQsData SDS;

    // Fill in input data
    SDS.tl1_1 = Tl1_1_eq;
    SDS.tl1_2 = Tl1_2_eq;
    SDS.pval = pval_eq;

    // Check the inputs
    REQUIRE(SDS.tl1_1 == SDS.tl1_2);
    REQUIRE(SDS.pval > 0);

    // Call the fortran implementation
    shoc_assumed_pdf_compute_qs(SDS);

    // Check the result
    REQUIRE(SDS.beta1 == SDS.beta2);
    REQUIRE(SDS.qs1 == SDS.qs2);

    // TEST TWO
    // Pressure test.  Use the data from the last test and modify
    //  the pressure level to be lower and verify that qs1 is lower.

    static constexpr Real pval_high = 100000;

    // Save the result from last test
    Real qs1_test1 = SDS.qs1;

    // Load new input
    SDS.pval = pval_high;

    // Verify new pressure is greater than last test
    REQUIRE(SDS.pval > pval_eq);

    // Call the fortran implementation
    shoc_assumed_pdf_compute_qs(SDS);

    // Check the result
    REQUIRE(SDS.qs1 < qs1_test1);

    // TEST THREE
    // Skewed test.  using pressre input from last test, but using
    //  a skewed temperature distribution.  Verify that the gaussian
    //  with lower input temperature also has lower qs value

    // Define temperature of Gaussian 1 [K]
    static constexpr Real Tl1_1_skew = 290;
    // Define temperature of Gaussian 2 [K]
    static constexpr Real Tl1_2_skew = 300;

    // Load new input
    SDS.tl1_1 = Tl1_1_skew;
    SDS.tl1_2 = Tl1_2_skew;

    // Verify Gaussians are not equal
    REQUIRE(SDS.tl1_1 != SDS.tl1_2);

    // Call the fortran implementation
    shoc_assumed_pdf_compute_qs(SDS);

    // Check the result
    if (SDS.tl1_1 < SDS.tl1_2){
      REQUIRE(SDS.qs1 < SDS.qs2);
      REQUIRE(SDS.beta1 > SDS.beta2);
    }
    else{
      REQUIRE(SDS.qs1 > SDS.qs2);
      REQUIRE(SDS.beta1 < SDS.beta2);
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

TEST_CASE("shoc_pdf_compute_qs_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeQs;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_compute_qs_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeQs;

  TestStruct::run_bfb();
}

} // namespace
