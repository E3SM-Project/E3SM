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
struct UnitWrap::UnitTest<D>::TestShocUpdateDse {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the SHOC function
    //     update_host_dse

    // FIRST TEST and SECOND TEST
    //  1) Verify that DSE increases with height and 2) verify that 
    //   points that contain cloud condensate result in a larger
    //   dse if all other inputs are equal.  

    // Liquid water potential temperature [K]
    static constexpr Real thlm[nlev] = {350.0, 325.0, 315.0, 310.0, 300.0};
    // Exner function [-]
    static constexpr Real exner[nlev] = {0.1, 0.3, 0.5, 0.7, 1.0};
    // Cloud condensate [g/kg] (converted to kg/kg when input)
    static constexpr Real shoc_ql[nlev] = {0.005, 0.08, 0.04, 0.03, 0.001};
    // Mid-point heights [m]
    static constexpr Real zt_grid[nlev] = {9000.0, 6000.0, 4000.0, 2000.0, 1000.0};
    // Surface geopotential
    static constexpr Real phis = 100.0;

    // Initialzie data structure for bridgeing to F90
    SHOCEnergydseData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    // for this test we need exactly two columns
    REQUIRE(SDS.shcol == 2);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < SDS.shcol; ++s) {
      SDS.phis[s] = phis;
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;

	SDS.thlm[offset] = thlm[n];
	SDS.zt_grid[offset] = zt_grid[n];
	SDS.exner[offset] = exner[n];

	// convert to [kg/kg]
	// Force the first column of cloud liquid 
	//   to be zero!
	SDS.shoc_ql[offset] = s*shoc_ql[n]/1000.0;

      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < SDS.shcol; ++s) {
      REQUIRE(SDS.phis[s] >= 0.0);
      for (Int n = 0; n < SDS.nlev; ++n){
	const auto offset = n + s * SDS.nlev;

	REQUIRE(SDS.thlm[offset] > 0.0);
	REQUIRE(SDS.shoc_ql[offset] >= 0.0);
	REQUIRE(SDS.exner[offset] > 0.0);
	REQUIRE(SDS.zt_grid[offset] >= 0.0);

	// make sure the two columns are different and
	//  as expected for the relevant variables
	if (s == 0){
          const auto offsets = n + (s+1) * SDS.nlev;

          REQUIRE(SDS.shoc_ql[offset] == 0.0);
          REQUIRE(SDS.shoc_ql[offsets] > SDS.shoc_ql[offset]);
	}    

	// Check that heights and thlm increase upward
	if (n < nlev-1){
          REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0.0);
	  REQUIRE(SDS.thlm[offset + 1] - SDS.thlm[offset] < 0.0);
	}  

      }
    }

    // Call the fortran implementation
    update_host_dse(SDS);

    // Check test
    for(Int s = 0; s < SDS.shcol; ++s) {
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;
	// Verify dse is reasonable
	REQUIRE(SDS.host_dse[offset] > 0.0);

	// Verify that dse increases with height upward
	if (n < nlev-1){
          REQUIRE(SDS.host_dse[offset + 1] - SDS.host_dse[offset] < 0.0);
	} 

	// Verify that dse is greater in points with condensate loading
	if (s == 0){
          const auto offsets = n + (s+1) * SDS.nlev;
	  REQUIRE(SDS.host_dse[offsets] > SDS.host_dse[offset]);
	}       
      }
    }
  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_energy_host_dse_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocUpdateDse;

  TestStruct::run_property();
}

TEST_CASE("shoc_energy_host_dse_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocUpdateDse;

}

} // namespace
