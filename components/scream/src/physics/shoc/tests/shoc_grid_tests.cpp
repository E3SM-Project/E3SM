#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"

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
struct UnitWrap::UnitTest<D>::TestShocGrid {

  static void run_property()
  {
    static constexpr Real gravit  = scream::physics::Constants<Real>::gravit;
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Define the midpoint height grid [m]
    static constexpr Real zt_pts[nlev] = {10000., 5000., 1000., 500., 100.};
    // Define the interface height grid [m]
    static constexpr Real zi_pts[nlevi] = {12500., 7500., 3000., 750., 250.0, 0.};
    // Define the air density [kg/m3]
    static constexpr Real density_zt[nlev] = {0.4, 0.6, 0.8, 1.0, 1.2};

    // Initialzie data structure for bridgeing to F90
    SHOCGridData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol() == shcol && SDS.nlev() == nlev && SDS.nlevi()) );
    REQUIRE(nlevi - nlev == 1);
    REQUIRE(shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	SDS.zt_grid[offset] = zt_pts[n];
	SDS.pdel[offset]    = density_zt[n] * gravit * (zi_pts[n]-zi_pts[n+1]);
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < nlevi; ++n) {
	const auto offset   = n + s * nlevi;
	SDS.zi_grid[offset] = zi_pts[n];
      }
    }

    // Check that the inputs make sense

    // Check that zt decreases upward
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev - 1; ++n) {
	const auto offset = n + s * nlev;
	REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0.0);
      }

      // Check that zi decreases upward
      for(Int n = 0; n < nlevi - 1; ++n) {
	const auto offset = n + s * nlevi;
	REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0.0);
      }
    }

    // Call the fortran implementation
    shoc_grid(SDS);

    // First check that dz is correct
    for(Int s = 0; s < shcol; ++s) {
      Real zt_sum = 0;
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;
	REQUIRE(SDS.dz_zt[offset] > 0);
	REQUIRE(SDS.dz_zt[offset] == zi_pts[n] - zi_pts[n+1]);
	zt_sum += SDS.dz_zt[offset];
      }
      // Check that the sum of dz_zt is equal to the largest zi
      REQUIRE(zt_sum == SDS.zi_grid[0]);
    }

    for(Int s = 0; s < shcol; ++s) {
      const auto s_offset = s * nlevi;
      REQUIRE(SDS.dz_zi[s_offset] == 0);
      REQUIRE(SDS.dz_zi[s_offset + nlevi - 1] == zt_pts[nlev-1]);

      Real zi_sum = 0;
      for(Int n = 1; n < nlevi - 1; ++n) {
	const auto offset = n + s * nlevi;
	REQUIRE(SDS.dz_zi[offset] > 0.0);
	REQUIRE(SDS.dz_zi[offset] == zt_pts[n-1] - zt_pts[n]);
	zi_sum += SDS.dz_zi[offset];
      }
      // Check that the sum of dz_zi is equal to the largest zt
      zi_sum += SDS.dz_zi[nlevi - 1];
      REQUIRE(zi_sum == SDS.zt_grid[0]);
    }

    // Now check density
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	// check that the density is consistent with the hydrostatic approximation
	REQUIRE(abs(SDS.rho_zt[offset] - density_zt[n]) <= std::numeric_limits<Real>::epsilon());

	// check that the density has physically realistic values
	REQUIRE(SDS.rho_zt[offset] <= 2.0);
	REQUIRE(SDS.rho_zt[offset] > 0.0);
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

TEST_CASE("shoc_grid_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocGrid;

  TestStruct::run_property();
}

TEST_CASE("shoc_grid_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocGrid;

  TestStruct::run_bfb();
}

} // namespace
