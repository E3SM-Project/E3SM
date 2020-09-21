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
struct UnitWrap::UnitTest<D>::TestShocDiagObklen {

  static void run_property()
  {
    static constexpr Int shcol    = 5;

    // Tests for the SHOC function:
    //   shoc_diag_obklen

    // TEST ONE
    // Load up input arrays with a variety of conditions to test
    //  a plethora of regimes and conditions.

    //  For this function we test several one layer columns

    // Surface sensible heat flux (K m/s)
    static constexpr Real wthl_sfc[shcol] = {-0.03, 0.02, 0, 0.05, -0.05};
    // Surface latent heat flux (kg/kg m/s)
    static constexpr Real wqw_sfc[shcol] = {0, 1e-4, -2e-4, 0.007, -0.007};
    // Surface momentum flux (u-direction) [m2/s2]
    static constexpr Real uw_sfc[shcol] = {-0.02, 0.05, 0.1, 0, -0.9};
    // Surface momentum flux (v-direction) [m2/s2]
    static constexpr Real vw_sfc[shcol] = {-0.02, 0.05, 0.1, 0, -0.9};
    // Surface potential temperature [K]
    static constexpr Real thl_sfc[shcol] = {300, 250, 310, 270, 305};
    // Surface cloud liquid water [kg/kg]
    static constexpr Real cldliq_sfc[shcol] = {0, 1e-4, 2.5e-5, 0, 3.5e-3};
    // Surface water vapor [kg/kg]
    static constexpr Real qv_sfc[shcol] = {1e-2, 2e-2, 1.5e-2, 9e-3, 5e-3};

    // Initialzie data structure for bridgeing to F90
    SHOCObklenData SDS(shcol);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol() == shcol);
    REQUIRE(shcol > 0);

    // Fill in test data, all one dimensional
    for(Int s = 0; s < shcol; ++s) {
      SDS.wthl_sfc[s] = wthl_sfc[s];
      SDS.wqw_sfc[s] = wqw_sfc[s];
      SDS.uw_sfc[s] = uw_sfc[s];
      SDS.vw_sfc[s] = vw_sfc[s];
      SDS.thl_sfc[s] = thl_sfc[s];
      SDS.cldliq_sfc[s] = cldliq_sfc[s];
      SDS.qv_sfc[s] = qv_sfc[s];
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.thl_sfc[s] > 150);
      REQUIRE( (SDS.cldliq_sfc[s] >= 0 && SDS.cldliq_sfc[s] < 0.05) );
      REQUIRE( (SDS.qv_sfc[s] > 0 && SDS.qv_sfc[s] < 0.1) );
    }

    // Call the fortran implementation
    shoc_diag_obklen(SDS);

    // Check the result

    for (Int s = 0; s < shcol; ++s){
      // Verify that ustar is ALWAYS positive and greater than zero
      REQUIRE(SDS.ustar[s] > 0);
      // if surface moisture flux is zero then surface kinematic
      //  flux should be equal to surface flux
      if (SDS.wqw_sfc[s] == 0){
        REQUIRE(SDS.kbfs[s] == SDS.wthl_sfc[s]);
      }
      // Verify that Obukhov length is opposite sign of surface
      //  kinematic buoyancy flux
      if (SDS.kbfs[s] > 0){
        REQUIRE(SDS.obklen[s] < 0);
      }
      else{
        REQUIRE(SDS.obklen[s] > 0);
      }

      // Verify output falls within some reasonable bounds
      REQUIRE( (SDS.ustar[s] > -10 && SDS.ustar[s] < 10) );
      REQUIRE( (SDS.kbfs[s] > -10 && SDS.kbfs[s] < 10) );
      REQUIRE( (SDS.obklen[s] > -10 && SDS.obklen[s] < 1000));
    }

    // TEST TWO
    // Surface flux and dimensionless obuhov length relationship test
    // Feel columns all constant values, except for surface
    //  fluxes.  Be sure that as the sum of the surface fluxes
    //  INCREASES that the Obukhov length scale DECREASES

    // Surface sensible heat flux (K m/s)
    static constexpr Real wthl_sfc_test2[shcol] = {-3e-2, -2e-2, 1e-2, 5e-2, 7e-2};
    // Surface momentum flux (u-direction) [m2/s2]
    static constexpr Real uw_sfc_test2 = 0.03;
    // Surface momentum flux (v-direction) [m2/s2]
    static constexpr Real vw_sfc_test2 = -0.02;
    // Surface potential temperature [K]
    static constexpr Real thl_sfc_test2 = 300;
    // Surface cloud liquid water [kg/kg]
    static constexpr Real cldliq_sfc_test2 = 2.5e-5;
    // Surface water vapor [kg/kg]
    static constexpr Real qv_sfc_test2 = 2e-2;

    // Fill in test data, all one dimensional
    for(Int s = 0; s < shcol; ++s) {
      SDS.wthl_sfc[s] = wthl_sfc_test2[s];
      // Set surface vapor flux proportional to heat flux
      SDS.wqw_sfc[s] = wthl_sfc_test2[s]/100;
      SDS.uw_sfc[s] = uw_sfc_test2;
      SDS.vw_sfc[s] = vw_sfc_test2;
      SDS.thl_sfc[s] = thl_sfc_test2;
      SDS.cldliq_sfc[s] = cldliq_sfc_test2;
      SDS.qv_sfc[s] = qv_sfc_test2;
    }

    // Verify that sum of surface fluxes increases with columns
    for(Int s = 0; s < shcol-1; ++s) {
      REQUIRE(SDS.wthl_sfc[s+1]+SDS.wqw_sfc[s+1] >
              SDS.wthl_sfc[s]+SDS.wqw_sfc[s]);
    }

    // Call the fortran implementation
    shoc_diag_obklen(SDS);

    // Verify that DIMENSIONLESS obukhov length decreases as columns
    //   increases due to the increasing surface fluxes

    // Define the lowest model layer grid height [m]
    Real zt_low = 50;
    for(Int s = 0; s < shcol-1; ++s) {
      REQUIRE(zt_low/SDS.obklen[s+1] < zt_low/SDS.obklen[s]);
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

TEST_CASE("shoc_diag_obklen_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocDiagObklen;

  TestStruct::run_property();
}

TEST_CASE("shoc_diagobklen_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocDiagObklen;

  TestStruct::run_bfb();
}

} // namespace
