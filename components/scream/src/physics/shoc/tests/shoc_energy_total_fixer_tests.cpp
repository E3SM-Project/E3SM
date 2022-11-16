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
struct UnitWrap::UnitTest<D>::TestShocTotEnergyFixer {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the SHOC function
    //     shoc_energy_total_fixer

    // FIRST TEST
    // Surface flux test.  Have two columns, one with zero surface fluxes
    //   and the other with positive surface fluxes.  Verify the column
    //   with surface fluxes has greater total energy "before".

    // Timestep [s]
    static constexpr Real dtime = 300;
    // Number of macmic steps
    static constexpr Int nadv = 2;
    // Air density [km/m3]
    static constexpr Real rho_zt[nlev] = {0.4, 0.6, 0.7, 0.9, 1.0};
    // Interface heights [m]
    static constexpr Real zi_grid[nlevi] = {11000, 7500, 5000, 3000, 1500, 0};
    // Define integrated static energy, kinetic energy, water vapor,
    //  and liquid water respectively
    static constexpr Real se = 200;
    static constexpr Real ke = 150;
    static constexpr Real wv = 0.5;
    static constexpr Real wl = 0.1;
    // Define surface sensible heat flux [K m/s]
    static constexpr Real wthl_sfc = 0.5;
    // Define surface total water flux [kg/kg m/s]
    static constexpr Real wqw_sfc = 0.01;

    // Initialize data structure for bridging to F90
    ShocEnergyTotalFixerData SDS(shcol, nlev, nlevi, dtime, nadv);

    // Test that the inputs are reasonable.
    // for this test we need exactly two columns
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi && SDS.dtime == dtime && SDS.nadv == nadv) );
    REQUIRE(shcol == 2);
    REQUIRE(nlevi == nlev+1);

    for(Int s = 0; s < shcol; ++s) {
      // Set before and after integrals equal
      SDS.se_a[s] = se;
      SDS.se_b[s] = se;
      SDS.ke_a[s] = ke;
      SDS.ke_b[s] = ke;
      SDS.wv_a[s] = wv;
      SDS.wv_b[s] = wv;
      SDS.wl_a[s] = wl;
      SDS.wl_b[s] = wl;

      // Make first column be zero for the surface fluxes
      SDS.wthl_sfc[s] = s*wthl_sfc;
      SDS.wqw_sfc[s] = s*wqw_sfc;

      // Fill in test data on zt_grid.
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // For zt grid, set as midpoint of zi grid
        SDS.zt_grid[offset] = 0.5*(zi_grid[n]+zi_grid[n+1]);
        SDS.rho_zt[offset] = rho_zt[n];
      }
      // Fill in test data on zi_grid.
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.zi_grid[offset] = zi_grid[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;

        REQUIRE(SDS.zt_grid[offset] >= 0);
        REQUIRE(SDS.rho_zt[offset] > 0);

        // Check that heights increase upward
        if (n > nlev-1){
          REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
        }
      }
      for (Int n = 0; n < nlevi; ++n){
        const auto offset = n + s * nlevi;

        REQUIRE(SDS.zi_grid[offset] >= 0);

        // Check that heights increase upward
        if (n > nlevi-1){
          REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0);
        }
      }
    }

    // Call the fortran implementation
    shoc_energy_total_fixer(SDS);

    // Check test

    // For first column verify that total energies are the same
    REQUIRE(SDS.te_a[0] == SDS.te_b[0]);

    // Verify that second column "before" energy is greater than
    //  the first column, since here we have active surface fluxes
    REQUIRE(SDS.te_b[1] > SDS.te_b[0]);
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

TEST_CASE("shoc_energy_total_fixer_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocTotEnergyFixer;

  TestStruct::run_property();
}

TEST_CASE("shoc_energy_total_fixer_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocTotEnergyFixer;

  TestStruct::run_bfb();
}

} // namespace
