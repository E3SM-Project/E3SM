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
struct UnitWrap::UnitTest<D>::TestShocCheckTke {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    
    // Tests for SHOC subroutine 
    //  check_tke
    
    // Test to make sure that given negative values of TKE 
    // that these values are corrected
    
    static constexpr Real tke_input[nlev] = {-0.3, 0.4, -100, -1.5, 0.4};
    
    // Intialize data structure for bridging to F90
    SHOCCheckTkeData SDS(shcol, nlev);
    
    // Load the data
    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < SDS.nlev; ++n) {
        const auto offset = n + s * SDS.nlev;
	
	SDS.tke[offset] = tke_input[n];
      }
    }
    
    // call the fortran implementation 
    check_tke(SDS);
    
    // Check the result against the input values
    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < SDS.nlev; ++n) {
        const auto offset = n + s * SDS.nlev;
	
	// if input TKE was less than zero, verify it was adjusted
	if (tke_input[n] < 0){
	  REQUIRE(SDS.tke[offset] > 0);
	}
	// Else make sure TKE remains untouched
	else{
	  REQUIRE(SDS.tke[offset] == tke_input[n]);
	}	
      }
    } 

  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_check_tke_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocCheckTke;

  TestStruct::run_property();
}

TEST_CASE("shoc_check_tke_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocCheckTke;

}

} // namespace
