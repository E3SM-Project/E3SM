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
struct UnitWrap::UnitTest<D>::TestShocLinearInt {

  static void run_property()
  {
    static constexpr Real gravit  = scream::physics::Constants<Real>::gravit;
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;
    
    // TEST ONE 
    // Test interpolation goign from nlev to nlevi grid

    // Define the interface height grid [m]
    static constexpr Real zi_grid[nlevi] = {12500., 7500., 3000., 750., 250.0, 0.};
    // Define the liquid water potential temperature [K]
    //  on the midpoint grid
    static constexpr Real thetal_zt[nlev] = {320.0, 310.0, 300.0, 300.0, 306.0};
    // Define minimum threshold for this experiment
    static constexpr Real mintresh = 0.0;

    // Initialzie data structure for bridgeing to F90
    SHOCLinearintData SDS(shcol, nlev, nlevi, minthresh);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.nlevi - SDS.nlev == 1);
    REQUIRE(SDS.shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;
        
	// For zt grid heights compute as midpoint
	//  between interface heights
	SDS.x1[offset] = 0.5*(zi_grid[n]+zi_grid[n+1]);
	SDS.y1[offset] = thetal_zt[n];
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < SDS.nlevi; ++n) {
	const auto offset   = n + s * SDS.nlevi;
	SDS.x2[offset] = zi_grid[n];
      }
    }

    // Check that the inputs make sense

    // Check that zt decreases upward
    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < SDS.nlev - 1; ++n) {
	const auto offset = n + s * SDS.nlev;
	REQUIRE(SDS.x1[offset + 1] - SDS.x1[offset] < 0.0);
      }

      // Check that zi decreases upward
      for(Int n = 0; n < SDS.nlevi - 1; ++n) {
	const auto offset = n + s * SDS.nlevi;
	REQUIRE(SDS.x2[offset + 1] - SDS.x2[offset] < 0.0);
      }
    }

    // Call the fortran implementation
    linear_interp(SDS);
    
    // First check that all output temperatures are greater than zero

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;
	REQUIRE(SDS.y2[offset] > 0.0);
	// check boundary points
	// First upper boundary
	if (n == 0){
	  const auto uppergradient = SDS.y1[offset] - SDS.y1[offset+1];
	  if (uppergradient > 0.0){
	    // if upper gradient is positive then make sure that 
	    //  the temperature of the highest nlevi layer is greater
	    //  than the tempature of the highest nlev layer
	    REQUIRE(SDS.y2[offset] > SDS.y1[offset]);
	  }
	  if (uppergradient < 0.0){
	    REQUIRE(SDS.y2[offset] < SDS.y1[offset]);
	  }
	  if (uppergradient == 0.0){
	    REQUIRE(SDS.y2[offset] == SDS.y2[offset]);
	  }
	}
	
        // Now the lower boundary
	else if (n == nlev){
	  const auto lowergradient = SDS.y1[offset-1] - SDS.y1[offset];
	  if (lowergradient > 0.0){
	    // if lower gradient is positive then make sure that 
	    //  the temperature of the lowest nlevi layer is lower
	    //  than the tempature of the lowest nlev layer
	    REQUIRE(SDS.y2[offset+1] < SDS.y1[offset]);
	  }
	  if (lowergradient < 0.0){
	    REQUIRE(SDS.y2[offset+1] > SDS.y1[offset]);
	  }
	  if (uppergradient == 0.0){
	    REQUIRE(SDS.y2[offset+1] == SDS.y2[offset]);
	  }
	}
	
	// Now make sure all points are bounded as expected
	else{
	
	}
		
      }
    }

  }
  
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_linear_interp_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocLinearInt;

  TestStruct::run_property();
}

TEST_CASE("shoc_linear_interp_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocLinearInt;

}

} // namespace
