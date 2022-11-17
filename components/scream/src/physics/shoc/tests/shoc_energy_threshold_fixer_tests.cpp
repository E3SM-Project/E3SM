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
struct UnitWrap::UnitTest<D>::TestShocEnergyThreshFixer {

  static void run_property()
  {
    static constexpr Real mintke = scream::shoc::Constants<Real>::mintke;
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the SHOC function
    //     shoc_energy_threshold_fixer

    // TEST ONE
    // Set up a reasonable profile verify results are as expected

    // Host model TKE  [m2/s2]
    Real tke_input[nlev] = {mintke, mintke, 0.01, 0.4, 0.5};
    //  Pressure at interface [Pa]
    Real pint[nlevi] = {500e2, 600e2, 700e2, 800e2, 900e2, 1000e2};

    // Integrated total energy after SHOC.
    static constexpr Real te_a = 100;
    // Integrated total energy before SHOC
    static constexpr Real te_b = 110;

    // convert pressure to Pa
    for(Int n = 0; n < nlevi; ++n) {
      pint[n] = pint[n];
    }

    // Initialize data structure for bridging to F90
    ShocEnergyThresholdFixerData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(SDS.shcol > 1);
    REQUIRE(nlev+1 == nlevi);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      SDS.te_a[s] = te_a;
      SDS.te_b[s] = te_b;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.tke[offset] = tke_input[n];
      }

      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.pint[offset] = pint[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;

        REQUIRE(SDS.tke[offset] >= mintke);
      }
    }

    // Call the fortran implementation
    shoc_energy_threshold_fixer(SDS);

    // Verify the result
    for(Int s = 0; s < shcol; ++s) {
      // Make sure value of shoctop is within reasonable range
      REQUIRE(SDS.shoctop[s] < nlev);
      REQUIRE(SDS.shoctop[s] > 1);

      // Verify that shoctop represents what we want it to
      // Make sure that thickness that bounds shoctop is positive
      const auto offset_stopi = (SDS.shoctop[s]-1) + s * nlevi;
      const auto offset_bot = (nlevi-1) + s * nlevi;
      REQUIRE(SDS.pint[offset_bot] - SDS.pint[offset_stopi] > 0.0);

      if (SDS.shoctop[s] < nlev){
        const auto offset_stop = (SDS.shoctop[s]-1) + s * nlev;
        REQUIRE(SDS.tke[offset_stop] == mintke);
        REQUIRE(SDS.tke[offset_stop+1] > mintke);
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

TEST_CASE("shoc_energy_threshold_fixer_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyThreshFixer;

  TestStruct::run_property();
}

TEST_CASE("shoc_energy_threshold_fixer_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyThreshFixer;

  TestStruct::run_bfb();
}

} // namespace
