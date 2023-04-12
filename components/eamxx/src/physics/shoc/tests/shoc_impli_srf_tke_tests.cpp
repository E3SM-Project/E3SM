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
struct UnitWrap::UnitTest<D>::TestImpTkeSfcStress {

  static void run_property()
  {
    static constexpr Int shcol = 5;

    // Tests for the SHOC subroutine
    //   impli_srf_stress_term

    // TEST ONE
    // Feed in several columns worth of data and make sure
    //  the output is consistent.  Make sure that columns with higher
    //  surface stress result in greater surface tke flux.

    // Surface moment flux, zonal direction [m3/s3]
    static constexpr Real uw_sfc[shcol] = {0.03, -0.03, 0.1, 0, -0.1};
    // Surface moment flux, meridional direction [m3/s3]
    static constexpr Real vw_sfc[shcol] = {-0.01, -0.01, 0.3, 0, -0.3};

    // Initialize data structure for bridging to F90
    TkeSrfFluxTermData SDS(shcol);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol == shcol);
    REQUIRE(shcol > 1);

    // Fill in test data, column only
    for(Int s = 0; s < shcol; ++s) {
      SDS.uw_sfc[s] = uw_sfc[s];
      SDS.vw_sfc[s] = vw_sfc[s];
    }

    // Call the fortran implementation
    tke_srf_flux_term(SDS);

    Real stress1, stress2;
    // Verify that output is as expected and reasonable
    for(Int s = 0; s < shcol; ++s) {
      // term should be greater than zero and less than one given
      //  reasonable input values
      REQUIRE(SDS.wtke_sfc[s] > 0);
      if (s < shcol-1){
        stress1 = uw_sfc[s]*uw_sfc[s] + uw_sfc[s]*uw_sfc[s];
        stress2 = uw_sfc[s+1]*uw_sfc[s+1] + uw_sfc[s+1]*uw_sfc[s+1];
        if (stress1 > stress2){
          REQUIRE(SDS.wtke_sfc[s] > SDS.wtke_sfc[s+1]);
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

TEST_CASE("shoc_imp_tkesfc_stress_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestImpTkeSfcStress;

  TestStruct::run_property();
}

TEST_CASE("shoc_imp_tkesfc_stress_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestImpTkeSfcStress;

  TestStruct::run_bfb();
}

} // namespace
