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
struct UnitWrap::UnitTest<D>::TestShocThlParameters {

  static void run_property()
  {
  
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_thl_parameters

    // TEST ONE
    // run the test two times given identical inputs, except 
    // one test where thlsec is larger.  Verify that thl2_1 and 
    // th2_2 are larger

    // Define the flux of thetal [K m/s]
    static constexpr Real wthlsec_test1 = 0.02;
    // Define the variance of vertical velocity [m/s]
    static constexpr Real sqrtw2_test1 = 0.5;
    // Define the variance of thetal [K] small value
    static constexpr Real thlsec_small = 1;
    // Define the standard deviation of thetal [K] large value
    static constexpr Real thlsec_large = 1.5;
    // Define grid mean thetal [K]
    static constexpr Real thl_first_test1 = 290;
    // Define first gaussian vertical velocity [m/s]
    static constexpr Real w1_1_test1 = 1.5;
    // Define second gaussian vertical velocity [m/s]
    static constexpr Real w1_2_test1 = -0.1;
    // Define vertical velocity skewness [-]
    static constexpr Real Skew_w_test1 = 3;
    // Define fraction of first gaussian
    static constexpr Real a_test1 = 0.2;
    // Define logical 
    static constexpr bool dothetal_skew = false; 
    
    // Be sure this value is actually larger
    REQUIRE(thlsec_large > thlsec_small);
    
    // Initialize data structure for bridging to F90
    SHOCPDFthlparamData SDS;
    
    // Fill the test data
    SDS.wthlsec = wthlsec_test1;
    SDS.sqrtw2 = sqrtw2_test1;
    SDS.thlsec = thlsec_small;
    SDS.sqrtthl = sqrt(thlsec_small);
    SDS.thl_first = thl_first_test1;
    SDS.w1_1 = w1_1_test1;
    SDS.w1_2 = w1_2_test1;
    SDS.Skew_w = Skew_w_test1;
    SDS.a = a_test1;
    SDS.dothetal_skew = dothetal_skew;
    
    // Verify input is physical
    REQUIRE(SDS.sqrtw2 >= 0);
    REQUIRE(SDS.sqrtthl >= 0);
    REQUIRE(SDS.thl_first > 0);
    
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

    // Call the fortran implementation
    shoc_assumed_pdf_thl_parameters(SDS);
    
    // Save some results to compare with next test
    // save the gaussian temperature variances
    Real thl2_1_result1 = SDS.thl2_1;
    Real thl2_2_result1 = SDS.thl2_2;
    
    // Now load up input data with a larger temperature standard deviation
    SDS.thlsec = thlsec_large;
    SDS.sqrtthl = sqrt(thlsec_large); 

    // Call the fortran implementation again 
    shoc_assumed_pdf_thl_parameters(SDS);    
    
    // Save some results to compare with the old test
    // save the gaussian temperature variances
    Real thl2_1_result2 = SDS.thl2_1;
    Real thl2_2_result2 = SDS.thl2_2;    
    
    // Now check the result
    REQUIRE(thl2_1_result2 > thl2_1_result1);
    REQUIRE(thl2_2_result2 > thl2_2_result1);
    
    // TEST TWO
    // For a case with identical input except wthlsec, verify that 
    // the case with higher wthlsec has larger |thl1_2-thl1_1|
    
    // Define the flux of thetal [K m/s], small value
    static constexpr Real wthlsec_small = 0.02;
    // Define the large value
    static constexpr Real wthlsec_large = 0.05;
    
    REQUIRE(wthlsec_large > wthlsec_small);
    
    //load the value for the small test, all other inputs
    //  will be recycled from last test
    SDS.wthlsec = wthlsec_small;
    
    // Call Fortran implementation 
    shoc_assumed_pdf_thl_parameters(SDS); 
    
    // Save absolute difference between the two gaussian temps
    Real thlgaus_diff_result1 = abs(SDS.thl1_2 - SDS.thl1_1);
    
    // Now laod up value for the large wthlsec test
    SDS.wthlsec = wthlsec_large;

    // Call Fortran implementation 
    shoc_assumed_pdf_thl_parameters(SDS);
    
    // Save absolute difference between the two gaussian temps
    Real thlgaus_diff_result2 = abs(SDS.thl1_2 - SDS.thl1_1); 
    
    // Now check the result
    REQUIRE(thlgaus_diff_result2 > thlgaus_diff_result1);       
      
  }
  
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_pdf_thl_parameters_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocThlParameters;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_thl_parameters_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocThlParameters;

}

} // namespace
