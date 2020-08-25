#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_arch.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestLInfShocLength {

  static void run_property()
  {
    static constexpr Int shcol    = 3;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC function:
    //   compute_l_inf_shoc_length   

    // Multi-column test, where the input heights for zt_grid
    //  are increased uniformly by 100 m per column to verify
    //  that l_inf always gets larger per column.        

    // Grid difference centered on thermo grid [m]
    static constexpr Real dz_zt[nlev] = {100.0, 100.0, 100.0, 100.0, 100.0};
    // Grid height centered on thermo grid [m]
    static constexpr Real zt_grid[nlev] = {500.0, 400.0, 300.0, 200.0, 100.0};
    // Turbulent kinetic energy [m2/s2]
    static constexpr Real tke[nlev] = {0.1, 0.15, 0.2, 0.25, 0.3};

    // Initialize data structure for bridgeing to F90
    SHOCInflengthData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    //  At least two columns are needed for this test!
    REQUIRE(SDS.shcol > 1);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;

	SDS.dz_zt[offset] = dz_zt[n];
	// Testing identical columns but one with larger zt heights.
	//  first column set as "base" column, and the others
	//  to a larger value of TKE uniformly.
	SDS.zt_grid[offset] = (s*100.0)+zt_grid[n];  
	SDS.tke[offset] = tke[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;
	// Be sure that relevant variables are greater than zero
	REQUIRE(SDS.dz_zt[offset] > 0.0);
	REQUIRE(SDS.tke[offset] > 0.0); 
	REQUIRE(SDS.zt_grid[offset] > 0.0);
	if (n < nlev-1){
          // check that zt increases upward
          REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0.0);       
	}      
	if (s < shcol-1){
          // Verify that zt_grid is offset larger column by column
          const auto offsets = n + (s+1) * SDS.nlev;
	  REQUIRE(SDS.zt_grid[offset] < SDS.zt_grid[offsets]);
	}
      } 
    }

    // Call the fortran implementation
    compute_l_inf_shoc_length(SDS);

    // Check the results
    // Make sure that conv_vel is negative
    for(Int s = 0; s < SDS.shcol; ++s) {   
      REQUIRE(SDS.l_inf[s] > 0.0);
    } 

    // Make sure that l_inf is getting larger
    //  per column 
    for (Int s = 0; s < SDS.shcol-1; ++s){
      REQUIRE(SDS.l_inf[s] < SDS.l_inf[s+1]);
    }
  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_l_inf_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestLInfShocLength;

  TestStruct::run_property();
}

TEST_CASE("shoc_l_inf_length_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestLInfShocLength;

}

} // namespace
