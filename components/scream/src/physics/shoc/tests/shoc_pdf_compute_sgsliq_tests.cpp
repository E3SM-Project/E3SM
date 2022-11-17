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
struct UnitWrap::UnitTest<D>::TestShocPdfComputeSgsLiq {

  static void run_property()
  {
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_compute_sgs_liquid

    // TEST ONE
    // Symmetric test.  For conditions of zero skewness, verify that
    //  that the SGS cloud condensate is equal to that of the two
    //  Gaussians.

    // Define the gaussian fraction
    static constexpr Real a = 0.5;
    // Define the liquid water mixing ratio, gaussian 1 [kg/kg]
    static constexpr Real ql1 = 0.02;
    // Define the liquid water mixing ratio, gaussian 2 [kg/kg]
    // for this test, set equal to the first
    static constexpr Real ql2 = ql1;

    // Initialize data structure for bridging to F90
    ShocAssumedPdfComputeSgsLiquidData SDS;

    // Fill in the data
    SDS.a = a;
    SDS.ql1 = ql1;
    SDS.ql2 = ql2;

    // Check the inputs
    // For this test should be symmetric
    REQUIRE(SDS.a == 0.5);
    REQUIRE(SDS.ql1 == SDS.ql2);

    // Call the fortran implementation
    shoc_assumed_pdf_compute_sgs_liquid(SDS);

    // Check the result
    REQUIRE(SDS.shoc_ql == SDS.ql1);

    // TEST TWO
    // Negative test.  Assume that this routine is wrongly fed negative values
    //  of ql1 and ql2.  Verify that the result is zero.

    // Define the gaussian fraction
    static constexpr Real a_neg = 0.2;
    // Define the liquid water mixing ratio, gaussian 1 [kg/kg]
    static constexpr Real ql1_neg = -0.02;
    // Define the liquid water mixing ratio, gaussian 2 [kg/kg]
    static constexpr Real ql2_neg =-0.01;

    // Fill in the data
    SDS.a = a_neg;
    SDS.ql1 = ql1_neg;
    SDS.ql2 = ql2_neg;

    // Check the inputs
    REQUIRE(SDS.a > 0);
    REQUIRE(SDS.a < 1);

    // For this test we want mixing ratios to be negative
    REQUIRE(SDS.ql1 < 0);
    REQUIRE(SDS.ql2 < 0);

    // Call the fortran implementation
    shoc_assumed_pdf_compute_sgs_liquid(SDS);

    // Check the result.  Make sure mixing ratio is zero
    REQUIRE(SDS.shoc_ql == 0);

    // TEST THREE
    // Positively skewed test. For test with positive skewness, verify that
    //  the resulting SGS cloud water lies between the two input values, given
    //  a dry downdraft.

    // Define the liquid water mixing ratio, gaussian 1 [kg/kg]
    static constexpr Real ql1_skew = 0.02;
    // Define the liquid water mixing ratio, gaussian 2 [kg/kg]
    static constexpr Real ql2_skew =0;

    // Fill in the data
    SDS.a = a_neg;
    SDS.ql1 = ql1_skew;
    SDS.ql2 = ql2_skew;

    // Check the data
    REQUIRE (SDS.a < 0.5);
    REQUIRE (SDS.ql1 > 0);
    REQUIRE (SDS.ql2 == 0);

    // Call the fortran implementation
    shoc_assumed_pdf_compute_sgs_liquid(SDS);

    // Check the result.  Make sure SDS value is between the lower and
    //  upper bound of the Gaussians
    REQUIRE(SDS.shoc_ql > SDS.ql2);
    REQUIRE(SDS.shoc_ql < SDS.ql1);
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

TEST_CASE("shoc_pdf_compute_sgsliq_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeSgsLiq;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_compute_sgsliq_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeSgsLiq;

  TestStruct::run_bfb();
}

} // namespace
