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
struct UnitWrap::UnitTest<D>::TestCompShocConvVel {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static const auto approx_zero = Approx(0.0).margin(1e-16);

    // Tests for the SHOC function:
    //   compute_conv_vel_shoc_length

    // Note that the output conv_vel from this subroutine
    //   represents the integrated buoyancy flux and not technically
    //   the convective velocity scale, which is the cubed root
    //   of conv_vel.

    // FIRST TEST
    // Verify that given negatively buoyant conditions that the
    //   result is also negative.

    // PBL height [m]
    // for this test set to be higher than highest zt_grid value
    static constexpr Real pblh = 1000.0;
    // Grid difference centered on thermo grid [m]
    static constexpr Real dz_zt[nlev] = {100.0, 100.0, 100.0, 100.0, 100.0};
    // Grid height centered on thermo grid [m]
    static constexpr Real zt_grid[nlev] = {500.0, 400.0, 300.0, 200.0, 100.0};
    // Virtual potential temperature on interface grid [K]
    static constexpr Real thv[nlev] = {310.0, 305.0, 300.0, 300.0, 295.0};
    // Buoyancy flux [K m/s]
    static constexpr Real wthv_sec[nlev] = {-0.02, -0.01, -0.04, -0.02, -0.05};

    // Initialzie data structure for bridgeing to F90
    SHOCConvvelData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol() == shcol && SDS.nlev() == nlev) );
    REQUIRE(shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      SDS.pblh[s] = pblh;
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	SDS.dz_zt[offset] = dz_zt[n];
	SDS.zt_grid[offset] = zt_grid[n];
	SDS.thv[offset] = thv[n];
	SDS.wthv_sec[offset] = wthv_sec[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;
	// Be sure that relevant variables are greater than zero
	REQUIRE(SDS.dz_zt[offset] > 0.0);
	REQUIRE(SDS.thv[offset] > 0.0);
	// For this negative buoyancy test, verify that all parcels
	//  have negative buoyancy flux.
	REQUIRE(SDS.wthv_sec[offset] < 0.0);
	if (n < nlev-1){
          // check that zt increases upward
          REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0.0);
	}
      }
    }

    // Call the fortran implementation
    compute_conv_vel_shoc_length(SDS);

    // Check the results
    // Make sure that conv_vel is negative
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.conv_vel[s] < 0.0);
    }

    // SECOND TEST
    // Symmetrical input test.  Here we feed input values of wthv_sec
    //  that are symmetrical around zero.  This is to verify that the
    //  result of conv_vel is also zero.  Here we need to assume
    //  a well mixed layer and constant profile of dz_zt

    static constexpr Real dz_zt_sym[nlev] = {100.0, 100.0, 100.0, 100.0, 100.0};
    // Virtual potential temperature on interface grid [K]
    static constexpr Real thv_sym[nlev] = {300.0, 300.0, 300.0, 300.0, 300.0};
    // Buoyancy flux [K m/s]
    static constexpr Real wthv_sec_sym[nlev] = {0.04, 0.02, 0.00, -0.02, -0.04};

    // Fill in new test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	SDS.dz_zt[offset] = dz_zt_sym[n];
	SDS.thv[offset] = thv_sym[n];
	SDS.wthv_sec[offset] = wthv_sec_sym[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      Real wthv_sum = 0.0;
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;
	// Be sure that relevant variables are greater than zero
	REQUIRE(SDS.dz_zt[offset] > 0.0);
	REQUIRE(SDS.thv[offset] > 0.0);
	wthv_sum += SDS.wthv_sec[offset];
      }
      // Make sure inputs of buoyancy flux sum to zero
      REQUIRE(wthv_sum == approx_zero);
    }

    // Call the fortran implementation
    compute_conv_vel_shoc_length(SDS);

    // Check the results
    // Make sure that conv_vel is zero
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.conv_vel[s] == approx_zero);
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

TEST_CASE("shoc_conv_vel_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompShocConvVel;

  TestStruct::run_property();
}

TEST_CASE("shoc_conv_vel_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompShocConvVel;

  TestStruct::run_bfb();
}

} // namespace
