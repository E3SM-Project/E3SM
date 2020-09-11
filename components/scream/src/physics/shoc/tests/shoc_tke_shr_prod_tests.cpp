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
struct UnitWrap::UnitTest<D>::TestShocShearProd {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the subroutine compute_shr_prod in the SHOC
    //   TKE module.

    // FIRST TEST
    //  For first tests input a sheared profile for both wind
    //  components, one with zonal winds increasing with height
    //  at a constant rate per GRID box and another with meridional
    //  winds decreasing at a constant rate per GRID box. The
    //  grid prescribed will be a stretched grid. Here we want to
    //  validate that the shear term DECREASES with height.
    //  NOTE: sterm boundaries will be returned as ZERO.

    // Define height thickness on nlevi grid [m]
    //   NOTE: First indicee is zero because it is never used
    //   Do a stretched grid
    static constexpr Real dz_zi[nlevi] = {0.0, 500., 200., 100., 50., 10.};
    // Define zonal wind on nlev grid [m/s]
    static constexpr Real u_wind_shr[nlev] = {2.0, 1.0, 0.0, -1.0, -2.0};
    // Define meridional wind on nlev grid [m/s]
    static constexpr Real v_wind_shr[nlev] = {1.0, 2.0, 3.0, 4.0, 5.0};

    // Initialzie data structure for bridgeing to F90
    SHOCTkeshearData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol() == shcol && SDS.nlev() == nlev && SDS.nlevi() == nlevi) );
    REQUIRE(nlevi - nlev == 1);
    REQUIRE(shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	SDS.u_wind[offset] = u_wind_shr[n];
	SDS.v_wind[offset] = v_wind_shr[n];
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < nlevi; ++n) {
	const auto offset   = n + s * nlevi;
	SDS.dz_zi[offset] = dz_zi[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlevi; ++n){
	const auto offset = n + s * nlevi;
	// Make sure top level dz_zi value is zero
	if (n == 0){
          REQUIRE(SDS.dz_zi[offset] == 0.0);
	}
	// Otherwise, should be greater than zero
	else{
          REQUIRE(SDS.dz_zi[offset] > 0.0);
	}
      }
    }

    // Call the fortran implementation
    compute_shr_prod(SDS);

    // Check test
    for(Int s = 0; s < shcol; ++s) {
      // First check that sterm is ALWAYS greater than
      //  zero for non boundary points, but exactly zero
      //  for boundary points.
      for(Int n = 0; n < nlevi; ++n) {
	const auto offset = n + s * nlevi;
	if (n == 0 || n == nlevi-1){
          // Boundary point check
          REQUIRE(SDS.sterm[offset] == 0.0);
	}
	else{
          REQUIRE(SDS.sterm[offset] > 0.0);
	}
      }
      // Now validate that shear term is ALWAYS
      //  decreasing with height for these inputs, keeping
      //  in mind to exclude boundary points, which should be zero
      for(Int n = 1; n < nlevi-2; ++n){
	const auto offset = n + s * nlevi;
	REQUIRE(SDS.sterm[offset]-SDS.sterm[offset+1] < 0.0);
      }
    }

    // SECOND TEST
    // For second test we input wind profiles that are
    // constant with height to validate that shear production
    // term is zero everywhere.

    // Define zonal wind on nlev grid [m/s]
    static constexpr Real u_wind_cons[nlev] = {10.0, 10.0, 10.0, 10.0, 10.0};
    // Define meridional wind on nlev grid [m/s]
    static constexpr Real v_wind_cons[nlev] = {-5.0, -5.0, -5.0, -5.0, -5.0};

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;

	SDS.u_wind[offset] = u_wind_cons[n];
	SDS.v_wind[offset] = v_wind_cons[n];
      }
    }

    // Call the fortran implementation
    compute_shr_prod(SDS);

    // Check test
    // Verify that shear term is zero everywhere
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
	const auto offset = n + s * nlevi;
	REQUIRE(SDS.sterm[offset] == 0.0);
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

TEST_CASE("shoc_tke_shr_prod_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocShearProd;

  TestStruct::run_property();
}

TEST_CASE("shoc_tke_shr_prod_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocShearProd;

  TestStruct::run_bfb();
}

} // namespace
