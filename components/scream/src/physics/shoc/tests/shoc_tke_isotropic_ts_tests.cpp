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
struct UnitWrap::UnitTest<D>::TestShocIsotropicTs {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 1;

    // Tests for the subroutine isotropic_ts
    //   in the SHOC TKE module.

    // For this routine all inputs are on midpoint grid and
    //  there are no vertical derivatives, therefore we will
    //  consider one vertical level per test.  Each column will
    //  be loaded up with a different test.

    // FIRST TEST
    // Stability test.  Verify that two points with the same inputs
    //  but one with brunt > 0 and the other with brunt < 0 that the
    //  timescale with brunt > 0 is less.

    // Integrated brunt vaisalla
    static constexpr Real brunt_int_st = 0.1;
    // TKE [m2/s2]
    static constexpr Real tke_st = 0.4;
    // Dissipation rate [m2/s3]
    static constexpr Real diss_st = 0.1;
    // Brunt Vaisalla frequency [/s]
    static constexpr Real brunt_st[shcol] = {-0.004, 0.004};

    // Initialzie data structure for bridgeing to F90
    SHOCIsotropicData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < SDS.shcol; ++s) {
      // Column only input
      SDS.brunt_int[s] = brunt_int_st;
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;

	SDS.tke[offset] = tke_st;
	SDS.a_diss[offset] = diss_st;
	SDS.brunt[offset] = brunt_st[s];
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < SDS.shcol; ++s) {
      for (Int n = 0; n < SDS.nlev; ++n){
	const auto offset = n + s * SDS.nlev;
	// Should be greater than zero
	REQUIRE(SDS.tke[offset] > 0.0);
	REQUIRE(SDS.a_diss[offset] > 0.0);
      }
    }

    // Call the fortran implementation
    isotropic_ts(SDS);

    // Check to make sure that column with positive
    //  brunt vaisalla frequency is smaller
    for(Int s = 0; s < SDS.shcol-1; ++s) {
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;
	// Get value corresponding to next column
	const auto offsets = n + (s+1) * SDS.nlev;
	if(SDS.brunt[offset] < 0.0 & SDS.brunt[offsets] > 0.0){
          REQUIRE(SDS.isotropy[offset] > SDS.isotropy[offsets]);
	}
	else{
          REQUIRE(SDS.isotropy[offset] < SDS.isotropy[offsets]);
	}
      }
    }

    // SECOND TEST
    // Dissipation test.  Verify that two points with the same inputs
    //  but one with smaller dissipation rate, that that point has
    //  a higher isotropy timescale.

    // Integrated brunt vaisalla
    static constexpr Real brunt_int_diss = 0.1;
    // TKE [m2/s2]
    static constexpr Real tke_diss = 0.4;
    // Dissipation rate [m2/s3]
    static constexpr Real diss_diss[shcol] = {0.1, 0.2};
    // Brunt Vaisalla frequency [/s]
    static constexpr Real brunt_diss = 0.004;

    // Fill in test data on zt_grid.
    for(Int s = 0; s < SDS.shcol; ++s) {
      SDS.brunt_int[s] = brunt_int_diss;
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;

	SDS.tke[offset] = tke_diss;
	SDS.a_diss[offset] = diss_diss[s];
	SDS.brunt[offset] = brunt_diss;
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < SDS.shcol; ++s) {
      for (Int n = 0; n < SDS.nlev; ++n){
	const auto offset = n + s * SDS.nlev;
	// Should be greater than zero
	REQUIRE(SDS.tke[offset] > 0.0);
	REQUIRE(SDS.a_diss[offset] > 0.0);
      }
    }

    // Call the fortran implementation
    isotropic_ts(SDS);

    // Check to make sure that column with positive
    //  brunt vaisalla frequency is smaller
    for(Int s = 0; s < SDS.shcol-1; ++s) {
      for(Int n = 0; n < SDS.nlev; ++n) {
	const auto offset = n + s * SDS.nlev;
	// Get value corresponding to next column
	const auto offsets = n + (s+1) * SDS.nlev;
	if(SDS.a_diss[offset] < SDS.a_diss[offsets]){
          REQUIRE(SDS.isotropy[offset] > SDS.isotropy[offsets]);
	}
	else{
          REQUIRE(SDS.isotropy[offset] < SDS.isotropy[offsets]);
	}
      }
    }
  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_tke_isotropic_ts_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocIsotropicTs;

  TestStruct::run_property();
}

TEST_CASE("shoc_tke_isotropic_ts_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocIsotropicTs;

}

} // namespace
