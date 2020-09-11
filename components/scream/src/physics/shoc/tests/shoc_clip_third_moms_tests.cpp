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
struct UnitWrap::UnitTest<D>::TestClipThirdMoms {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlevi    = 5;

    // Tests for the SHOC function:
    //   clipping_diag_third_shoc_moments
  
    //  Test to be sure that very high values of w3
    //    are reduced but still the same sign   
  
    // Define the second moment of vertical velocity
    static constexpr Real w_sec_zi[nlevi] = {0.3, 0.4, 0.5, 0.4, 0.1}; 
    // Define the third moment of vertical velocity
    static constexpr Real w3_in[nlevi] = {0.2, 99999, -99999, 0.2, 0.05}; 
    
    // Define a local logical
    bool w3_large;

    // Initialize data structure for bridging to F90
    SHOCClipthirdmomsData SDS(shcol, nlevi);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol() > 0);
    REQUIRE(SDS.nlevi() > 1);
    
    // load up the input data
    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.w_sec_zi[offset] = w_sec_zi[n];
	SDS.w3[offset] = w3_in[n];
      }
    }
    
    // Verify the data
    // For this particular test wen want to be sure that
    //   the input has relatively small values of w2 (let's say
    //   all smaller than 1, which is reasonable) and let's 
    //   be sure that w3 has some very large unreasonable values
    for(Int s = 0; s < shcol; ++s) {
      // Initialize conditional to make sure there are
      //   large values of w3 in input
      w3_large = false;
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        REQUIRE(SDS.w_sec_zi[offset] <= 1);
	if (abs(SDS.w3[offset]) > 1000){
	  w3_large = true;
	}
        SDS.w_sec_zi[offset] = w_sec_zi[n];
	SDS.w3[offset] = w3_in[n];
      }
      REQUIRE(w3_large == true);
    }
    
    // Call the fortran implementation    
    clipping_diag_third_shoc_moments(SDS);
    
    // Check the result
    // For large values of w3, verify that the result has been 
    //  reduced and is of the same sign
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
	
	if (abs(w3_in[n]) > 1000){
	  REQUIRE(abs(SDS.w3[offset]) < abs(w3_in[n]));
	  if (w3_in[n] < 0){
	    REQUIRE(SDS.w3[offset] < 0);
	  }
	} 
	
      }
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

namespace{

TEST_CASE("shoc_clip_third_moms_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestClipThirdMoms;
  
  TestStruct::run_property();
}

TEST_CASE("shoc_clip_third_moms_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestClipThirdMoms;

  TestStruct::run_bfb();
}

} // namespace
