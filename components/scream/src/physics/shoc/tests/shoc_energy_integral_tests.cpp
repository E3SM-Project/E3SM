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
struct UnitWrap::UnitTest<D>::TestShocEnergyInt {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC function
    //     shoc_energy_integrals

    // FIRST TEST and SECOND TEST
    //   FIRST: kinetic energy test.  Given at least two columns with
    //   identical inputs but with increasing winds in each column
    //   verify that the column with stronger winds results
    //   in a higher kinetic energy test.
    //   SECOND: vapor test.  verify that columns with no vapor
    //   or moisture result in zero integral outputs.  If the
    //   column is positive then verify that wl_int >= wv_int

    // Define host model dry static energy [J kg-1]
    static constexpr Real host_dse[nlev] = {350e3, 325e3, 315e3, 310e3, 300e3};
    // Defin the pressure difference [hPa] (converted to Pa later)
    static constexpr Real pdel[nlev]={100.0, 75.0, 50.0, 25.0, 10.0};
    // Define zonal wind on nlev grid [m/s]
    static constexpr Real u_wind[nlev] = {-4.0, -4.0, -3.0, -1.0, -2.0};
    // Define meridional wind on nlev grid [m/s]
    static constexpr Real v_wind[nlev] = {1.0, 2.0, 3.0, 4.0, 5.0};
    // Define the total water mixing ratio [g/kg] (converted to kg/kg later)
    static constexpr Real rtm[nlev] = {7.0, 9.0, 11.0, 15.0, 20.0};
    // Define cloud mixing ratio [g/kg] (converted to kg/kg later)
    static constexpr Real rcm[nlev] = {0.001, 0.01, 0.04, 0.002, 0.001};

    // Initialzie data structure for bridgeing to F90
    SHOCEnergyintData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    // for this test we need exactly two columns
    REQUIRE( (SDS.shcol() == shcol && SDS.nlev() == nlev) );
    REQUIRE(shcol == 2);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	// Add one degree K in the second column
	SDS.host_dse[offset] = s+host_dse[n];

	// convert to [kg/kg]
	// Force the first column of cloud liquid
	//   to be zero!
	SDS.rcm[offset] = s*rcm[n]/1000.0;
	SDS.rtm[offset] = rtm[n]/1000.0;

	// convert to Pa
	SDS.pdel[offset] = pdel[n]*100.0;

	// Increase winds with increasing columns
	SDS.u_wind[offset] = (1.0+s)*u_wind[n];
	SDS.v_wind[offset] = (1.0+s)*v_wind[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
	const auto offset = n + s * nlev;

	REQUIRE(SDS.host_dse[offset] > 0.0);
	REQUIRE(SDS.rcm[offset] >= 0.0);
	REQUIRE(SDS.rtm[offset] > 0.0);

	// make sure the two columns are different and
	//  as expected for the relevant variables
	if (s == 0){
          const auto offsets = n + (s+1) * nlev;

          REQUIRE(abs(SDS.u_wind[offsets]) > abs(SDS.u_wind[offset]));
          REQUIRE(abs(SDS.v_wind[offsets]) > abs(SDS.v_wind[offset]));
          REQUIRE(SDS.rcm[offset] == 0.0);
          REQUIRE(SDS.rcm[offsets] > SDS.rcm[offset]);
	}
      }
    }

    // Call the fortran implementation
    shoc_energy_integrals(SDS);

    // Check test
    for(Int s = 0; s < shcol; ++s) {
      // Verify relevant integrals are reasonable
      REQUIRE(SDS.se_int[s] > 0.0);
      REQUIRE(SDS.ke_int[s] > 0.0);
      REQUIRE(SDS.wv_int[s] > 0.0);
      REQUIRE(SDS.wl_int[s] >= 0.0);
      // Do column comparison tests
      if (s == 0){
	REQUIRE(SDS.ke_int[s+1] > SDS.ke_int[s]);
	REQUIRE(SDS.wl_int[s+1] > SDS.wl_int[s]);
	REQUIRE(SDS.wl_int[s] == 0.0);
	REQUIRE(SDS.se_int[s+1] > SDS.se_int[s]);
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

namespace {

TEST_CASE("shoc_energy_integrals_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyInt;

  TestStruct::run_property();
}

TEST_CASE("shoc_energy_integrals_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyInt;

  TestStruct::run_bfb();
}

} // namespace
