#include "catch2/catch.hpp"

//#include "share/scream_types.hpp"
#include <algorithm>
#include <array>
#include <random>
#include <thread>

#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/scream_arch.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/util/scream_utils.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "shoc_unit_tests_common.hpp"

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
    static constexpr Real wqwsec_small = 0.0001;
    // Define the flux of total water [kg/kg m/s] small value
    static constexpr Real wqwsec_large = 0.0005;    
    // Define the standard deviation of vertical velocity [m/s]
    static constexpr Real sqrtw2_test1 = 0.5;
    // Define the variance of total water [kg^2/kg^2] small
    static constexpr Real qwsec_test1 = 1e-5;
    // Define grid mean total water [kg/kg]
    static constexpr Real qw_first_test1 = 0.015;
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
    
    // Be sure this value is actually larger
    REQUIRE(wqwsec_large > wqwsec_small);
    
    // Initialize data structure for bridging to F90
    SHOCPDFqwparamData SDS;
    
    // Fill the test data
    SDS.wqwsec = wqwsec_small;
    SDS.sqrtw2 = sqrtw2_test1;
    SDS.qwsec = qwsec_test1;
    SDS.sqrtqt = sqrt(qwsec_test1);
    SDS.qw_first = qw_first_test1;
    SDS.w1_1 = w1_1_test1;
    SDS.w1_2 = w1_2_test1;
    SDS.Skew_w = Skew_w_test1;
    SDS.a = a_test1;
    
    // Verify input is physical
    REQUIRE(SDS.sqrtw2 >= 0);
    REQUIRE(SDS.sqrtqt >= 0);
    REQUIRE(SDS.qw_first > 0);
    
    // Make sure vertical velocity data is consistent
    REQUIRE(SDS.w1_1 > 0);
    REQUIRE(SDS.w1_2 < 0);
    if (SDS.Skew_w > 0){
      REQUIRE(abs(SDS.w1_1) > abs(SDS.w1_2));
      REQUIRE(SDS.a < 0.5);
    }
    else if (SDS.Skew_w < 0){
      REQUIRE(abs(SDS.w1_1) < abs(SDS.w1_2));
      REQUIRE(SDS.a > 0.5);    
    }
    else if (SDS.Skew_w == 0){
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

TEST_CASE("shoc_pdf_qw_parameters_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocQwParameters;

}

} // namespace
