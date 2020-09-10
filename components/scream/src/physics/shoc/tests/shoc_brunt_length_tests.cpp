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
struct UnitWrap::UnitTest<D>::TestCompBruntShocLength {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the SHOC function:
    //   compute_brunt_shoc_length

    // Test for the Brunt Vaissalla frequency.
    // Note that input temperature profiles are selected
    //  deliberately so that it includes a well mixed layer,
    //  an unstable layer, and a conditionally unstable layer
    //  to test a range of conditions.

    // Grid difference centered on thermo grid [m]
    static constexpr Real dz_zt[nlev] = {100.0, 75.0, 50.0, 25.0, 10.0};
    // Virtual potential temperature on interface grid [K]
    static constexpr Real thv_zi[nlevi] = {310.0, 305.0, 300.0, 300.0, 295.0, 305.0};

    // Initialize data structure for bridgeing to F90
    SHOCBruntlengthData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol() == shcol && SDS.nlev() == nlev && SDS.nlevi() == nlevi) );
    REQUIRE(nlevi - nlev == 1);
    REQUIRE(shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.dz_zt[offset] = dz_zt[n];
        // For theta_v on thermo grid just take the vertical average
        //  of thv_zi for this simple test.  Just used as a reference
        //  in this subroutine.
        SDS.thv[offset] = 0.5*(thv_zi[n]+thv_zi[n+1]);
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset   = n + s * nlevi;
        SDS.thv_zi[offset] = thv_zi[n];
      }
    }

    // Check that the inputs make sense

    // Be sure that relevant variables are greater than zero
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev - 1; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.dz_zt[offset] > 0.0);
        REQUIRE(SDS.thv[offset] > 0.0);
      }

      for(Int n = 0; n < nlevi - 1; ++n) {
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.thv_zi[offset] > 0.0);
      }
    }

    // Call the fortran implementation
    compute_brunt_shoc_length(SDS);

    // Check the results
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // Validate that brunt vaisalla frequency
        //  is the correct sign given atmospheric conditions

        // well mixed layer
        if (thv_zi[n] - thv_zi[n+1] == 0.0){
          REQUIRE(SDS.brunt[offset] == 0.0);
        }
        // unstable layer
        if (thv_zi[n] - thv_zi[n+1] < 0.0){
          REQUIRE(SDS.brunt[offset] < 0.0);
        }
        // stable layer
        if (thv_zi[n] - thv_zi[n+1] > 0.0){
          REQUIRE(SDS.brunt[offset] > 0.0);
        }

        // Validate that values fall within some
        //  reasonable bounds for this variable.
        REQUIRE(SDS.brunt[offset < 1.0]);
        REQUIRE(SDS.brunt[offset > -1.0]);
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

namespace{

TEST_CASE("shoc_brunt_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompBruntShocLength;

  TestStruct::run_property();
}

TEST_CASE("shoc_brunt_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompBruntShocLength;

  TestStruct::run_bfb();
}

} // namespace
