#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "physics/share/physics_constants.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_setup_random_test.hpp"

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
    static constexpr Real zt_pts[nlev] = {10000, 5000, 1000, 500, 100};
    // Define the interface height grid [m]
    static constexpr Real zi_pts[nlevi] = {12500, 7500, 3000, 750, 250, 0};
    // Define the air density [kg/m3]
    static constexpr Real density_zt[nlev] = {0.4, 0.6, 0.8, 1.0, 1.2};

    // Initialize data structure for bridging to F90
    ShocGridData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
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
        REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
      }

      // Check that zi decreases upward
      for(Int n = 0; n < nlevi - 1; ++n) {
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0);
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
        REQUIRE(SDS.rho_zt[offset] <= 2);
        REQUIRE(SDS.rho_zt[offset] > 0);
      }
    }
  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ShocGridData f90_data[] = {
      ShocGridData(10, 71, 72),
      ShocGridData(10, 12, 13),
      ShocGridData(7,  16, 17),
      ShocGridData(2, 7, 8),
    };

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ShocGridData cxx_data[] = {
      ShocGridData(f90_data[0]),
      ShocGridData(f90_data[1]),
      ShocGridData(f90_data[2]),
      ShocGridData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      shoc_grid(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      shoc_grid_f(d.shcol, d.nlev, d.nlevi, d.zt_grid, d.zi_grid, d.pdel, d.dz_zt, d.dz_zi, d.rho_zt);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(ShocGridData);
      for (Int i = 0; i < num_runs; ++i) {
        ShocGridData& d_f90 = f90_data[i];
        ShocGridData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.dz_zt); ++k) {
          REQUIRE(d_f90.dz_zt[k] == d_cxx.dz_zt[k]);
          REQUIRE(d_f90.rho_zt[k] == d_cxx.rho_zt[k]);
        }
        for (Int k = 0; k < d_f90.total(d_f90.dz_zi); ++k) {
          REQUIRE(d_f90.dz_zi[k] == d_cxx.dz_zi[k]);
        }
      }
    }
  } // run_bfb
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
