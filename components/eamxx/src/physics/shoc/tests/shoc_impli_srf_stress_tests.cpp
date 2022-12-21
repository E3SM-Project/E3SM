#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
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
struct UnitWrap::UnitTest<D>::TestImpSfcStress {

  static void run_property()
  {
    static constexpr Int shcol = 5;

    // Tests for the SHOC subroutine
    //   impli_srf_stress_term

    // TEST ONE
    // Feed in several columns worth of data and make sure
    //  the output is consistent.

    // Surface density on the zi grid [kg/m3]
    static constexpr Real rho_zi_sfc[shcol] = {1.2, 1.0, 0.9, 1.1, 1.15};
    // Surface moment flux, zonal direction [m3/s3]
    static constexpr Real uw_sfc[shcol] = {0.03, -0.03, 0.1, 0, -0.1};
    // Surface moment flux, meridional direction [m3/s3]
    static constexpr Real vw_sfc[shcol] = {-0.01, -0.01, 0.3, 0, -0.3};
    // Surface wind, zonal direction [m/s]
    static constexpr Real u_wind_sfc[shcol] = {5, -5, 0, 2, -10};
    // Surface wind, meridional direction [m/s]
    static constexpr Real v_wind_sfc[shcol] = {-10, 2, 20, 0, 1};

    // Initialize data structure for bridging to F90
    ImpliSrfStressTermData SDS(shcol);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol == shcol);
    REQUIRE(shcol > 1);

    // Fill in test data, column only
    for(Int s = 0; s < shcol; ++s) {
      SDS.rho_zi_sfc[s] = rho_zi_sfc[s];
      SDS.uw_sfc[s] = uw_sfc[s];
      SDS.vw_sfc[s] = vw_sfc[s];
      SDS.u_wind_sfc[s] = u_wind_sfc[s];
      SDS.v_wind_sfc[s] = v_wind_sfc[s];
    }

    // Call the fortran implementation
    impli_srf_stress_term(SDS);

    // Verify that output is reasonable
    for(Int s = 0; s < shcol; ++s) {
      // term should be greater than zero and less than one given
      //  reasonable input values
      REQUIRE( (SDS.ksrf[s] > 0 && SDS.ksrf[s] < 1) );
    }

    // TEST TWO
    // Given inputs that are identical but the absolute value of
    //   the surface fluxes are INCREASING, verify ksrf value is larger.
    //   Can recycle input from surface fluxes from last test.

    // Fill in test data, column only
    for(Int s = 0; s < shcol; ++s) {
      SDS.rho_zi_sfc[s] = 1.2; // density [kg/m3]
      SDS.u_wind_sfc[s] = 5; // zonal wind [m/s]
      SDS.v_wind_sfc[s] = -10; // meridional wind [m/s]
    }

    // Call the fortran implementation
    impli_srf_stress_term(SDS);

    Real stress1, stress2;
    // Verify that output is as expected and reasonable
    for(Int s = 0; s < shcol; ++s) {
      // term should be greater than zero and less than one given
      //  reasonable input values
      REQUIRE( (SDS.ksrf[s] > 0 && SDS.ksrf[s] < 1) );
      if (s < shcol-1){
        stress1 = uw_sfc[s]*uw_sfc[s] + uw_sfc[s]*uw_sfc[s];
        stress2 = uw_sfc[s+1]*uw_sfc[s+1] + uw_sfc[s+1]*uw_sfc[s+1];
        if (stress1 > stress2){
          REQUIRE(SDS.ksrf[s] > SDS.ksrf[s+1]);
        }
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

TEST_CASE("shoc_imp_sfc_stress_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestImpSfcStress;

  TestStruct::run_property();
}

TEST_CASE("shoc_imp_sfc_stress_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestImpSfcStress;

  TestStruct::run_bfb();
}

} // namespace
