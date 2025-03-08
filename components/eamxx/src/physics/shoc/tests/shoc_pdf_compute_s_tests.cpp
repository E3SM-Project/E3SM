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
struct UnitWrap::UnitTest<D>::TestShocPdfComputeS {

  static void run_property()
  {
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_compute_s

    // TEST ONE
    // Saturation test.  Given inputs that should produce saturation,
    //  verify this condition is met.

    // define Gaussian total water [kg/kg]
    static constexpr Real qw1 = 0.020;
    // define qs for gaussian 1 [kg/kg]
    static constexpr Real qs1 = 0.018;
    // define beta for gaussian 1
    static constexpr Real beta = 0.17;
    // define pressure [Pa]
    static constexpr Real pval = 85000;
    // define temperature variance for gaussian 1 [K^2]
    static constexpr Real thl2 = 2;
    // define total water variance for gaussian 1 [kg^2/kg^2]
    static constexpr Real qw2 = 1e-5;
    // define correlation between total water and temperature
    static constexpr Real r_qwthl = 0.5;

    // Initialize data structure for bridging to F90
    ShocAssumedPdfComputeSData SDS;

    // fill in data
    SDS.qw1 = qw1;
    SDS.qs1 = qs1;
    SDS.beta = beta;
    SDS.pval = pval;
    SDS.thl2 = thl2;
    SDS.qw2 = qw2;
    SDS.sqrtthl2 = sqrt(thl2);
    SDS.sqrtqw2 = sqrt(qw2);
    SDS.r_qwthl = r_qwthl;

    // Call the fortran implementation
    shoc_assumed_pdf_compute_s(SDS);

    // Verify saturation has occurred
    REQUIRE(SDS.s > 0);
    REQUIRE(SDS.std_s > 0);
    REQUIRE(SDS.qn > 0);
    REQUIRE(SDS.c > 0);

    // TEST TWO
    // Variance test.  Given the same inputs as the previous test, but
    //  with the variances doubled, verify that s is the same as the previous
    //  test (because this does not depend on variance), but that std_s and
    //  qn have INCREASES but C has DECREASED.

    // Save output from previous test

    Real s_test1 = SDS.s;
    Real std_s_test1 = SDS.std_s;
    Real qn_test1 = SDS.qn;
    Real C_test1 = SDS.c;

    // Fill in new values for variance, double the values
    //  uses in the previous test

    // define temperature variance for gaussian 1 [K^2]
    static constexpr Real thl2_doub = 2*thl2;
    // define total water variance for gaussian 1 [kg^2/kg^2]
    static constexpr Real qw2_doub = 2*qw2;

    REQUIRE(thl2_doub > thl2);
    REQUIRE(qw2_doub > qw2);

    SDS.thl2 = thl2_doub;
    SDS.qw2 = qw2_doub;
    SDS.sqrtthl2 = sqrt(thl2_doub);
    SDS.sqrtqw2 = sqrt(qw2_doub);

    // Call the fortran implementation
    shoc_assumed_pdf_compute_s(SDS);

    // check the result
    REQUIRE(SDS.s == s_test1);
    REQUIRE(SDS.std_s > std_s_test1);
    REQUIRE(SDS.qn > qn_test1);
    REQUIRE(SDS.c < C_test1);

    //TEST THREE
    //Unsaturated test.  Given input with unsaturated conditions, verify
    // that the output is as expected.
    static constexpr Real qw1_unsat = 0.016;
    // define qs for gaussian 1 [kg/kg]
    static constexpr Real qs1_unsat = 0.018;
    // define temperature variance for gaussian 1 [K^2]
    static constexpr Real thl2_unsat = 0;
    // define total water variance for gaussian 1 [kg^2/kg^2]
    static constexpr Real qw2_unsat = 0;

    // load data
    SDS.qw1 = qw1_unsat;
    SDS.qs1 = qs1_unsat;
    SDS.thl2 = thl2_unsat;
    SDS.qw2 = qw2_unsat;
    SDS.sqrtthl2 = sqrt(thl2_unsat);
    SDS.sqrtqw2 = sqrt(qw2_unsat);

    // Check data
    REQUIRE(SDS.qw1 < SDS.qs1);

    // Call the fortran implementation
    shoc_assumed_pdf_compute_s(SDS);

    // Verify output
    REQUIRE(SDS.s < 0);
    REQUIRE(SDS.c == 0);
    REQUIRE(SDS.qn == 0);
    REQUIRE(SDS.std_s == 0);
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

TEST_CASE("shoc_pdf_compute_s_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeS;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_compute_s_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeS;

  TestStruct::run_bfb();
}

} // namespace
