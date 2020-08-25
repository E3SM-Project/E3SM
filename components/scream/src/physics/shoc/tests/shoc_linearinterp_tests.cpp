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
struct UnitWrap::UnitTest<D>::TestShocLinearInt {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int km1     = 5;
    static constexpr auto km2   = km1 + 1;
    
    // TEST ONE 
    // Test interpolation going from midpoint grid to interface grid.  Note that
    //  in this case the nlev grid is denoted by variable km1 and the nlevi 
    //  grid is represented by variable km2.  This is because the linear
    //  interp routine can interpolate from midpoint to interface grid or
    //  vice versa so notation should be flexible.   
    
    // Note that in this test we are testing TWO profils.  The first column
    //  will be loaded up with potential temperature values while the second
    //  column will be loaded up with meridional wind values.  

    // Define the interface height grid [m]
    static constexpr Real zi_grid[km2] = {12500., 7500., 3000., 750., 250.0, 0.};
    // Define the liquid water potential temperature [K]
    //  on the midpoint grid
    static constexpr Real thetal_zt[km1] = {320.0, 310.0, 300.0, 300.0, 306.0};
    // Define meridional wind on midpoint grid
    static constexpr Real u_wind_zt[km1] = {-2.0, -0.5, 0.0, 5.0, 2.0};
    // Define minimum threshold for this experiment, since we are
    //  dealing with negative winds in a column then set to a large
    //  negative value since we do not want to clip.  
    static constexpr Real minthresh = -99999.0;

    // Initialzie data structure for bridgeing to F90
    SHOCLinearintData SDS(shcol, km1, km2, minthresh);

    // For this test we need exactly two columns
    REQUIRE(SDS.shcol == 2);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < km1; ++n) {
	const auto offset = n + s * km1;
        
	// For zt grid heights compute as midpoint
	//  between interface heights
	SDS.x1[offset] = 0.5*(zi_grid[n]+zi_grid[n+1]);
	if (s == 0){
	  SDS.y1[offset] = thetal_zt[n];
	}
	else{
	  SDS.y1[offset] = u_wind_zt[n];
	}
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < km2; ++n) {
	const auto offset   = n + s * km2;
	SDS.x2[offset] = zi_grid[n];
      }
    }

    // Check that the inputs make sense

    // Check that zt decreases upward
    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < km1 - 1; ++n) {
	const auto offset = n + s * km1;
	REQUIRE(SDS.x1[offset + 1] - SDS.x1[offset] < 0.0);
      }

      // Check that zi decreases upward
      for(Int n = 0; n < km2 - 1; ++n) {
	const auto offset = n + s * km2;
	REQUIRE(SDS.x2[offset + 1] - SDS.x2[offset] < 0.0);
      }
    }

    // Call the fortran implementation
    linear_interp(SDS);
    
    // First check that all output temperatures are greater than zero

    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < km1; ++n) {
	const auto offset = n + s * km1;
	const auto offseti = n + s * km2;
	
	// check boundary points
	// First upper boundary
	if (n == 0){
	  const auto uppergradient = SDS.y1[offset] - SDS.y1[offset+1];
	  if (uppergradient > 0.0){
	    // if upper gradient is positive then make sure that 
	    //  the temperature of the highest nlevi layer is greater
	    //  than the tempature of the highest nlev layer
	    REQUIRE(SDS.y2[offseti] > SDS.y1[offset]);
	  }
	  if (uppergradient < 0.0){
	    REQUIRE(SDS.y2[offseti] < SDS.y1[offset]);
	  }
	  if (uppergradient == 0.0){
	    REQUIRE(SDS.y2[offseti] == SDS.y2[offseti]);
	  }
	}
	
        // Now the lower boundary
	else if (n == km1-1){
	  const auto lowergradient = SDS.y1[offset-1] - SDS.y1[offset];
	  if (lowergradient > 0.0){
	    // if lower gradient is positive then make sure that 
	    //  the temperature of the lowest nlevi layer is lower
	    //  than the tempature of the lowest nlev layer
	    REQUIRE(SDS.y2[offseti+1] < SDS.y1[offset]);
	  }
	  if (lowergradient < 0.0){
	    REQUIRE(SDS.y2[offseti+1] > SDS.y1[offset]);
	  }
	  if (lowergradient == 0.0){
	    REQUIRE(SDS.y2[offseti+1] == SDS.y2[offset]);
	  }
	}
	
	// Now make sure all points are bounded as expected
	else{
	  const auto gradient = SDS.y1[offset-1] - SDS.y1[offset];
	  if (gradient == 0.0){
	    REQUIRE(SDS.y2[offseti] == SDS.y1[offset]);
	  } 
	  else if (gradient > 0.0){
	    REQUIRE(SDS.y2[offseti] < SDS.y1[offset-1]);
	    REQUIRE(SDS.y2[offseti] > SDS.y1[offset]);
	  }
	  else {
	    REQUIRE(SDS.y2[offseti] > SDS.y1[offset-1]);
	    REQUIRE(SDS.y2[offseti] < SDS.y1[offset]);	    
	  }
	}
		
      }
    }

//  TEST TWO
//  Now we test going from interface grid to mid point grid

    // Initialzie data structure for bridgeing to F90
    // NOTE that km2 and km1 grid must be switched here.
    // Must initialize a new data structure since km1 and km2 are swapped.
    SHOCLinearintData SDS2(shcol, km2, km1, minthresh); 

    // Fill in test data on zi_grid.
    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < km2; ++n) {
	const auto offset = n + s * km2;
        
	// Load up stuff on interface grid.  Here we
	//  are going to use information from the last test
	//  to initialize our grid
	SDS2.x1[offset] = SDS.x2[offset];
	SDS2.y1[offset] = SDS.y2[offset];
      }

      // Load up stuff on midpoint grid, use output from 
      //  the last test here.  
      for(Int n = 0; n < km1; ++n) {
	const auto offset   = n + s * km1;
	SDS2.x2[offset] = SDS.x1[offset];
      }
    }
    
    linear_interp(SDS2);
    
    // Check the result, make sure output is bounded correctly
   
    for(Int s = 0; s < SDS2.shcol; ++s) {
      for(Int n = 0; n < km1; ++n) {
	const auto offset = n + s * km1;
	const auto offseti = n + s * km2;

	// compute gradient
	const auto gradient = SDS2.y1[offseti] - SDS2.y1[offseti+1];

        if (gradient == 0.0){
	  REQUIRE(SDS2.y2[offset] == SDS2.y1[offseti]);
	} 
	else if (gradient > 0.0){
	  REQUIRE(SDS2.y2[offset] > SDS2.y1[offseti+1]);
	  REQUIRE(SDS2.y2[offset] < SDS2.y1[offseti]);
	}
	else {
	  REQUIRE(SDS2.y2[offset] < SDS2.y1[offseti+1]);
	  REQUIRE(SDS2.y2[offset] > SDS2.y1[offseti]);	    
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
