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
struct UnitWrap::UnitTest<D>::TestShocInPlumeCorr {

  static void run_property()
  {
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_inplume_correlations

    // Test One
    // Zero correlation test.  Test conditions that should always
    //  yield a zero correlation

    // standard deviation moisture gaussian 1 [kg/kg]
    static constexpr Real sqrtqw2_1_zero = 0;
    // standard deviation moisture gaussian 2 [kg/kg]
    static constexpr Real sqrtqw2_2_zero = 0;
    // standard deviation temperature gaussian 1 [K]
    static constexpr Real sqrtthl2_1_zero = 0;
    // standard deviation moisture gaussian 2 [K]
    static constexpr Real sqrtthl2_2_zero = 0;
    // covariance of temperature and moisture [kg/kg K]
    static constexpr Real qwthlsec_zero = 1.e-5;
    // moisture grid mean
    static constexpr Real qw_first_zero = 0.02;
    // moisture first gaussian [kg/kg]
    static constexpr Real qw1_1_zero = 0.020;
    // moisture second gaussian [kg/kg]
    static constexpr Real qw1_2_zero = 0.020;
    // Temperatur grid mean
    static constexpr Real thl_first_zero = 290;
    // Temperature first gaussian [K]
    static constexpr Real thl1_1_zero = 292;
    // Temperature second gaussian [K]
    static constexpr Real thl1_2_zero = 289;
    // Define gaussian fraction [-]
    static constexpr Real a_zero = 0.2;

    // Initialize data structure for bridging to F90
    ShocAssumedPdfInplumeCorrelationsData SDS;

    // Fill in test data
    SDS.sqrtqw2_1 = sqrtqw2_1_zero;
    SDS.sqrtqw2_2 = sqrtqw2_2_zero;
    SDS.sqrtthl2_1 = sqrtthl2_1_zero;
    SDS.sqrtthl2_2 = sqrtthl2_2_zero;
    SDS.qwthlsec = qwthlsec_zero;
    SDS.qw_first = qw_first_zero;
    SDS.qw1_1 = qw1_1_zero;
    SDS.qw1_2 = qw1_2_zero;
    SDS.thl_first = thl_first_zero;
    SDS.thl1_1 = thl1_1_zero;
    SDS.thl1_2 = thl1_2_zero;
    SDS.a = a_zero;

    // Call Fortran implementation
    shoc_assumed_pdf_inplume_correlations(SDS);

    // Check the result
    // Verify that correlation is zero
    REQUIRE(SDS.r_qwthl_1 == 0.0);

    // Test Two
    // Test conditions that should always yield a positive result
    // Note: some output is recycled

    // standard deviation moisture gaussian 1 [kg/kg]
    static constexpr Real sqrtqw2_1_pos = 3e-4;
    // standard deviation moisture gaussian 2 [kg/kg]
    static constexpr Real sqrtqw2_2_pos = 1e-4;
    // standard deviation temperature gaussian 1 [K]
    static constexpr Real sqrtthl2_1_pos = 0.7;
    // standard deviation moisture gaussian 2 [K]
    static constexpr Real sqrtthl2_2_pos = 0.2;
    // covariance of temperature and moisture [kg/kg K]
    static constexpr Real qwthlsec_pos = 1.e-5;

    // Fill in test data
    SDS.sqrtqw2_1 = sqrtqw2_1_pos;
    SDS.sqrtqw2_2 = sqrtqw2_2_pos;
    SDS.sqrtthl2_1 = sqrtthl2_1_pos;
    SDS.sqrtthl2_2 = sqrtthl2_2_pos;
    SDS.qwthlsec = qwthlsec_pos;

    // Call Fortran implementation
    shoc_assumed_pdf_inplume_correlations(SDS);

    // Check the result
    // Verify correlation is positive and not greater than one
    REQUIRE(SDS.r_qwthl_1 > 0);
    REQUIRE(SDS.r_qwthl_1 <= 1);

    // Test Three
    // Test conditions that should always yield a negative result

    // covariance of temperature and moisture [kg/kg K]
    static constexpr Real qwthlsec_neg = -1.e-5;

    // Fill in test data
    SDS.qwthlsec = qwthlsec_neg;

    // Call Fortran implementation
    shoc_assumed_pdf_inplume_correlations(SDS);

    // Check the result
    // Verify correlation is negative and not less than negative one
    REQUIRE(SDS.r_qwthl_1 < 0);
    REQUIRE(SDS.r_qwthl_1 >= -1);
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

TEST_CASE("shoc_pdf_inplume_corr_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocInPlumeCorr;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_inplume_corr_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocInPlumeCorr;

  TestStruct::run_bfb();
}

} // namespace
