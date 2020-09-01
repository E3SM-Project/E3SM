#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
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

    // standard deviation moisture gaussian 1 [kg/kg]
    static constexpr Real sqrtqw2_1 = 3e-3;
    // standard deviation moisture gaussian 2 [kg/kg]
    static constexpr Real sqrtqw2_2 = 1e-3;    
    // standard deviation temperature gaussian 1 [K]
    static constexpr Real sqrtthl2_1 = 0.7;
    // standard deviation moisture gaussian 2 [K]
    static constexpr Real sqrtthl2_2 = 0.2;
    // covariance of temperature and moisture [kg/kg K]
    static constexpr Real qwthlsec = 1.e-3;
    // moisture grid mean
    static constexpr Real qw_first = 0.02;
    // moisture first gaussian [kg/kg]
    static constexpr Real qw1_1 = 0.02;
    // moisture second gaussian [kg/kg]
    static constexpr Real qw1_2 = 0.02;
    // Temperatur grid mean
    static constexpr Real thl_first = 290;
    // Temperature first gaussian [K]
    static constexpr Real thl1_1 = 292;
    // Temperature second gaussian [K]
    static constexpr Real thl1_2 = 290;
    // Define gaussian fraction [-]
    static constexpr Real a = 0.2;
    
    // Initialize data structure for bridging to F90
    SHOCPDFinplumeData SDS;
    
    // Fill in test data
    SDS.sqrtqw2_1 = sqrtqw2_1;
    SDS.sqrtqw2_2 = sqrtqw2_2;
    SDS.sqrtthl2_1 = sqrtthl2_1;
    SDS.sqrtthl2_2 = sqrtthl2_2;
    SDS.qwthlsec = qwthlsec;
    SDS.qw_first = qw_first;
    SDS.qw1_1 = qw1_1;
    SDS.qw1_2 = qw1_2;
    SDS.thl_first = thl_first;
    SDS.thl1_1 = thl1_1;
    SDS.thl1_2 = thl1_2;
    SDS.a = a;
    
    // Verify input is expected and physical
    // For this test we want all temperatures and moisture equal
//    REQUIRE(SDS.thl_first == SDS.thl1_1);
//    REQUIRE(SDS.thl_first == SDS.thl1_2);
//    REQUIRE(SDS.qw_first == SDS.qw1_1);
//    REQUIRE(SDS.qw_first == SDS.qw1_2);    
    
    // Call Fortran implementation 
    shoc_assumed_pdf_inplume_correlations(SDS); 
    
    // Check the result
    // Verify that correlation is zero
//    REQUIRE(SDS.r_qwthl_1 == 0.0);      
      
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

TEST_CASE("shoc_pdf_inplume_corr_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocInPlumeCorr;

}

} // namespace
