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
struct UnitWrap::UnitTest<D>::TestShocEnergyDseFixer {

  static void run_property()
  {
    static constexpr Int shcol    = 3;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the SHOC function
    //     shoc_energy_dse_fixer
    
    // TEST ONE
    // For columns that are identical EXCEPT for the shoctop indicee, 
    //  verify that given a positive value of energy imbalance that 
    //  columns with a higher SHOC top were subject to more energy removal.

    // Host model dry static energy [K]
    static constexpr Real host_dse_input[nlev] = {350.0, 325.0, 315.0, 310.0, 300.0};

    // Energy disbalance.  For this test we assume all columns have
    //  the same disbalance magnitude.  
    static constexpr Real se_dis = 0.1;
    
    // level indicee of SHOC top layer
    static constexpr Int shoctop[shcol] = {5, 3, 1}; 


    // Initialzie data structure for bridgeing to F90
    SHOCEnergydsefixerData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    // for this test we need exactly three columns
    REQUIRE(SDS.shcol == 3);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < SDS.shcol; ++s) {
      SDS.shoctop[s] = shoctop[s];
      SDS.se_dis[s] = se_dis[s];
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;

	SDS.host_dse[offset] = host_dse_input[n];

      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < SDS.shcol; ++s) {
      // For this test we WANT se_dis > 0
      REQUIRE(SDS.se_dis[s] > 0.0);
      REQUIRE(SDS.shoc_top >= 1);
      REQUIRE(SDS.shoc_top <= nlev);
      for (Int n = 0; n < SDS.nlev; ++n){
	const auto offset = n + s * SDS.nlev;

	REQUIRE(SDS.host_dse[offset] > 0.0);

      }
    }

    // Call the fortran implementation
    shoc_energy_dse_fixer(SDS);

  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_energy_dse_fixer_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyDseFixer;

  TestStruct::run_property();
}

TEST_CASE("shoc_energy_dse_fixer_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyDseFixer;

}

} // namespace
