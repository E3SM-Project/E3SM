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
struct UnitWrap::UnitTest<D>::TestShocPdfComputeTemp {

  static void run_property()
  {
  
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_compute_temperature

    // TEST ONE
    // Keep liquid water potential temperature constant and
    //  decrease pressure by 10000 Pa each iteration.  Verify
    //  that the liquid water temperature behaves as expected 
    //  (i.e. temperature decreases with decreasing pressure and
    //  value is as expected relative to base pressure).  

    // Input liquid water potential temperature [K]
    static constexpr Real thl1 = 305;
    // Input basepressure [Pa]
    static constexpr Real basepres = 100000;
    // Input value of pval [Pa]
    Real pval = 110000;
    
    // Decrease base pres by a certain amount each test
    static constexpr Real presincr = -10000;
    
    Real Tl1_save;
    
    // Initialize data structure for bridging to F90
    SHOCPDFcomptempData SDS;
    
    // Fill in data
    SDS.thl1 = thl1;
    SDS.basepres = basepres;
    SDS.pval = pval;
    
    Int num_tests = SDS.pval/abs(presincr);
    
    REQUIRE(num_tests > 1);
    REQUIRE(presincr < 0);
    // Make sure our starting pressure is greater than 
    //  basepres just so we test a range
    REQUIRE(SDS.pval > SDS.basepres);
    
    for (Int s = 0; s < num_tests; ++s){
      
      // make sure base pres is greater than zero
      REQUIRE(basepres > 0);    

      // Call the fortran implementation
      shoc_assumed_pdf_compute_temperature(SDS);   
      
      // Check the result   
      // If pressure is greater than basepressure then
      //  make sure that temperature is greater than thetal
      if (SDS.pval > SDS.basepres){
        REQUIRE(SDS.Tl1 > SDS.thl1);
      } 
      // otherwise temperature should be less than thetal
      else if(SDS.pval < SDS.basepres){
        REQUIRE(SDS.Tl1 < SDS.thl1);
      }
      // otherwise if they are equal the temperatures
      //  should be equal
      else
      {
        REQUIRE(SDS.Tl1 == SDS.thl1);
      } 
      
      // Make sure temperature are decreasing
      if (s > 0){
        REQUIRE(SDS.Tl1 < Tl1_save);
      }
      
      // Save result to make sure that temperatures
      //  are decreasing as pressure decreases
      Tl1_save = SDS.Tl1; 
      
      // Decrease pressure value
      SDS.pval = SDS.pval+presincr;
     
    }

  }
  
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_pdf_computetemp_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeTemp;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_computetemp_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeTemp;

}

} // namespace
