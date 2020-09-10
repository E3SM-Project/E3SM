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
struct UnitWrap::UnitTest<D>::TestShocLength {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr Int nlevi    = nlev+1;

    // Tests for the upper level SHOC function:
    //   shoc_length
    
    // TEST ONE
    // Standard input and TKE test.  
    // With reasonable inputs, verify outputs follow.  
    // In addition, feed columns higher values of TKE and the length
    // scale for these columns should be gradually larger, given all other
    // inputs are equal.  

    // PBL height [m]
    static constexpr Real pblh = 1000.0;
    // Define the host grid box size x-direction [m]
    static constexpr Real host_dx = 3000; 
    // Defin the host grid box size y-direction [m]
    static constexpr Real host_dy = 5000; 
    // Define the heights on the zt grid [m]
    static constexpr Real zi_grid[nlevi] = {9000, 5000, 1500, 900, 500, 0};
    // Virtual potential temperature on thermo grid [K]
    static constexpr Real thv[nlev] = {315, 310, 305, 300, 295};    
    // Buoyancy flux [K m/s]
    static constexpr Real wthv_sec[nlev] = {0.02, 0.01, 0.04, 0.02, 0.05};   
    // Turbulent kinetc energy [m2/s2]
    static constexpr Real tke[nlev] = {0.1, 0.15, 0.2, 0.25, 0.3};
    
    // compute geometric grid mesh
    const auto grid_mesh = sqrt(host_dx*host_dy);
    
    // Grid stuff to compute based on zi_grid
    Real zt_grid[nlev];
    Real dz_zt[nlev];
    Real dz_zi[nlevi];
    // Compute heights on midpoint grid 
    for(Int n = 0; n < nlev; ++n) {
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
      dz_zt[n] = zi_grid[n] - zi_grid[n+1];
      if (n == 0){
        dz_zi[n] = 0;
      }
      else{
        dz_zi[n] = zt_grid[n-1] - zt_grid[n];
      }
    }
    // set upper condition for dz_zi
    dz_zi[nlevi-1] = zt_grid[nlev-1];

    // Initialize data structure for bridging to F90
    SHOCLengthData SDS(shcol, nlev, nlevi);
    
    // Load up input data
    for(Int s = 0; s < SDS.shcol; ++s) {    
      SDS.host_dx[s] = host_dx;
      SDS.host_dy[s] = host_dy;
      SDS.pblh[s] = pblh;
      // Fill in test data on zt_grid.     
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;
	
	// for subsequent columns, increase TKE
	SDS.tke[offset] = (s+1)*tke[n];
	
	SDS.zt_grid[offset] = zt_grid[n];
	SDS.thv[offset] = thv[n];
	SDS.wthv_sec[offset] = wthv_sec[n];
	SDS.dz_zt[offset] = dz_zt[n];	    
      }
      
      // Fill in test data on zi_grid
      for(Int n = 0; n < SDS.nlevi; ++n) {
	const auto offset = n + s * SDS.nlevi;
	
	SDS.zi_grid[offset] = zi_grid[n];
	SDS.dz_zi[offset] = dz_zi[n];
      }   
    }
    
    // Check input data
    // for this test we need more than one column
    REQUIRE(SDS.shcol > 1);
    REQUIRE(SDS.nlevi == SDS.nlev+1);    
    
    for(Int s = 0; s < SDS.shcol; ++s) {     
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;
	
	REQUIRE(SDS.zt_grid[offset] > 0);
	REQUIRE(SDS.tke[offset] > 0);
	REQUIRE(SDS.thv[offset] > 0);
	REQUIRE(SDS.dz_zt[offset] > 0);
	
	// Make sure TKE is larger in next column over
	if (s < SDS.shcol-1){
	  // get offset for "neighboring" column
	  const auto offsets = n + (s+1) * SDS.nlev;
	  REQUIRE(SDS.tke[offsets] > SDS.tke[offset]);
	}
      }
      
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;
	
	REQUIRE(SDS.zi_grid[offset] >= 0);
	REQUIRE(SDS.dz_zi[offset] >= 0);
      }
    }
    
    // Call the Fortran implementation 
    shoc_length(SDS);

    // Verify output
    for(Int s = 0; s < SDS.shcol; ++s) {   
      for(Int n = 0; n < SDS.nlev; ++n) { 
        const auto offset = n + s * SDS.nlev;
        // Require mixing length is greater than zero and is
        //  less than geometric grid mesh length + 1 m
        REQUIRE(SDS.shoc_mix[offset] > 0.0); 
        REQUIRE(SDS.shoc_mix[offset] < 1.0+grid_mesh);
	
	// Be sure brunt vaisalla frequency is reasonable
	REQUIRE(SDS.brunt[offset] < 1);
	REQUIRE(SDS.brunt[offset] > -1);
	
	// Make sure length scale is larger when TKE is larger
	if (s < SDS.shcol-1){
	  // get offset for "neighboring" column
	  const auto offsets = n + (s+1) * SDS.nlev;
	  if (SDS.tke[offsets] > SDS.tke[offset]){
	    REQUIRE(SDS.shoc_mix[offsets] > SDS.shoc_mix[offset]);
	  }
	  else{
	    REQUIRE(SDS.shoc_mix[offsets] < SDS.shoc_mix[offset]);
	  }
	}	
      }
    }
    
    // TEST TWO
    // Small grid mesh test.  Given a very small grid mesh, verify that
    //  the length scale is confined to this value.  Input from first
    //  test will be recycled into this one.
     
    // Define the host grid box size x-direction [m]
    static constexpr Real host_dx_small = 3;
    // Defin the host grid box size y-direction [m]
    static constexpr Real host_dy_small = 5;
    
    // compute geometric grid mesh
    const auto grid_mesh_small = sqrt(host_dx_small*host_dy_small);
    
    // Load new input
    for(Int s = 0; s < SDS.shcol; ++s) {
      SDS.host_dx[s] = host_dx_small;
      SDS.host_dy[s] = host_dy_small;
    }
    
    // call fortran implentation
    shoc_length(SDS);
    
    // Verify output
    for(Int s = 0; s < SDS.shcol; ++s) {   
      for(Int n = 0; n < SDS.nlev; ++n) { 
        const auto offset = n + s * SDS.nlev;
        // Require mixing length is greater than zero and is
        //  less than geometric grid mesh length + 1 m
        REQUIRE(SDS.shoc_mix[offset] > 0.0); 
        REQUIRE(SDS.shoc_mix[offset] < 1.0+grid_mesh_small);	
      }
    }

  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocLength;

  TestStruct::run_property();
}

TEST_CASE("shoc_length_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocLength;

}

} // namespace
