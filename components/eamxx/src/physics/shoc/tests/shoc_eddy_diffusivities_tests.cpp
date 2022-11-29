#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
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
struct UnitWrap::UnitTest<D>::TestShocEddyDiff {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 1;

    // Tests for the subroutine eddy_diffusivities
    //   in the SHOC TKE module.

    // For this routine all inputs are on midpoint grid and
    //  there are no vertical derivatives, therefore we will
    //  consider one vertical level per test.  Each column will
    //  be loaded up with a different test.

    // FIRST TEST
    // Boundary layer regime test.  Input that have identical values
    //  except for the Monin Obukhov length, in this case will be positive
    //  and negative.  Test to make sure that the resulting diffusivites
    //  are DIFFERENT

    // Monin Obukov length [m]
    static constexpr Real obklen_reg[shcol] = {-1, 1};
    // PBL depth [m]
    static constexpr Real pblh = 1000;
    // zt_grid [m]
    static constexpr Real zt_grid = 200;
    // SHOC Mixing length [m]
    static constexpr Real shoc_mix_reg = 1000;
    // Shear term [s-2]
    static constexpr Real sterm_zt_reg = 0.1;
    // Return to isotropy timescale [s]
    static constexpr Real isotropy_reg = 500;
    // Turbulent kinetic energy [m2/s2]
    static constexpr Real tke_reg = 0.4;

    // Initialize data structure for bridging to F90
    EddyDiffusivitiesData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    REQUIRE(shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      // Column only input
      SDS.pblh[s] = pblh;
      SDS.obklen[s] = obklen_reg[s];
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.tke[offset] = tke_reg;
        SDS.zt_grid[offset] = zt_grid;
        SDS.shoc_mix[offset] = shoc_mix_reg;
        SDS.sterm_zt[offset] = sterm_zt_reg;
        SDS.isotropy[offset] = isotropy_reg;
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] > 0);
      // Make sure point we are testing is within PBL
      REQUIRE(SDS.zt_grid[s] < SDS.pblh[s]);
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Should be greater than zero
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.zt_grid[offset] > 0);
        REQUIRE(SDS.shoc_mix[offset] > 0);
        REQUIRE(SDS.isotropy[offset] > 0);
        REQUIRE(SDS.sterm_zt[offset] > 0);
      }
    }

    // Call the fortran implementation
    eddy_diffusivities(SDS);

    // Check to make sure the answers in the columns are different
    for(Int s = 0; s < shcol-1; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Get value corresponding to next column
        const auto offsets = n + (s+1) * nlev;
        REQUIRE(SDS.tk[offset] != SDS.tk[offsets]);
        REQUIRE(SDS.tkh[offset] != SDS.tkh[offsets]);
      }
    }

    // SECOND TEST
    // Stable boundary layer test.  Given Obukhov length,
    //   verify that each regime behaves as expected when the relevant
    //   inputs are modified.  Column with larger mixing length and
    //   shear term should produce larger diffusivity values.

    // Monin Obukov length [m]
    static constexpr Real obklen_stab[shcol] = {1, 1};
    // SHOC Mixing length [m]
    static constexpr Real shoc_mix_stab[shcol] = {100, 150};
    // Shear term [s-2]
    static constexpr Real sterm_zt_stab[shcol] = {0.1, 0.2};
    // Return to isotropy timescale [s]
    static constexpr Real isotropy_stab = 500;
    // Turbulent kinetic energy [m2/s2]
    static constexpr Real tke_stab = 0.4;

    // Verify that input mixing length and shear term
    //   are increasing in each column for this
    //   test to be valid
    for(Int s = 0; s < shcol-1; ++s) {
      REQUIRE(shoc_mix_stab[s+1] > shoc_mix_stab[s]);
      REQUIRE(sterm_zt_stab[s+1] > sterm_zt_stab[s]);
    }

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      // Column only input
      SDS.obklen[s] = obklen_stab[s];
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.tke[offset] = tke_stab;
        SDS.shoc_mix[offset] = shoc_mix_stab[s];
        SDS.sterm_zt[offset] = sterm_zt_stab[s];
        SDS.isotropy[offset] = isotropy_stab;
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      // Make sure we are testing stable boundary layer
      REQUIRE(SDS.obklen[s] > 0);
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Should be greater than zero
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.shoc_mix[offset] > 0);
        REQUIRE(SDS.isotropy[offset] > 0);
        REQUIRE(SDS.sterm_zt[offset] > 0);
      }
    }

    // Call the fortran implementation
    eddy_diffusivities(SDS);

    // Check to make sure the answers in the columns are larger
    //   when the length scale and shear term are larger
    for(Int s = 0; s < shcol-1; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Get value corresponding to next column
        const auto offsets = n + (s+1) * nlev;
        if (SDS.shoc_mix[offset] < SDS.shoc_mix[offsets] &
            SDS.sterm_zt[offset] < SDS.sterm_zt[offsets]){
          REQUIRE(SDS.tk[offset] < SDS.tkh[offsets]);
          REQUIRE(SDS.tkh[offset] < SDS.tkh[offsets]);
        }
      }
    }

    // THIRD TEST
    // Unstable boundary layer test.  Given Obukhov length,
    //   verify that each regime behaves as expected when the relevant
    //   inputs are modified.

    // Monin Obukov length [m]
    static constexpr Real obklen_ustab[shcol] = {-1, -1};
    // SHOC Mixing length [m]
    static constexpr Real shoc_mix_ustab = 500;
    // Shear term [s-2]
    static constexpr Real sterm_zt_ustab = 0.1;
    // Return to isotropy timescale [s]
    static constexpr Real isotropy_ustab[shcol] = {500, 550};
    // Turbulent kinetic energy [m2/s2]
    static constexpr Real tke_ustab[shcol] = {0.4, 0.5};

    // Verify that input isotropy and tke
    //   are increasing in each column for this
    //   test to be valid
    for(Int s = 0; s < shcol-1; ++s) {
      REQUIRE(isotropy_ustab[s+1] > isotropy_ustab[s]);
      REQUIRE(tke_ustab[s+1] > tke_ustab[s]);
    }

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      // Column only input
      SDS.obklen[s] = obklen_ustab[s];
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.tke[offset] = tke_ustab[s];
        SDS.shoc_mix[offset] = shoc_mix_ustab;
        SDS.sterm_zt[offset] = sterm_zt_ustab;
        SDS.isotropy[offset] = isotropy_ustab[s];
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      // Make sure we are testing unstable boundary layer
      REQUIRE(SDS.obklen[s] < 0);
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Should be greater than zero
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.shoc_mix[offset] > 0);
        REQUIRE(SDS.isotropy[offset] > 0);
        REQUIRE(SDS.sterm_zt[offset] > 0);
      }
    }

    // Call the fortran implementation
    eddy_diffusivities(SDS);

    // Check to make sure the diffusivities are smaller
    //  in the columns where isotropy and tke are smaller
    for(Int s = 0; s < shcol-1; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Get value corresponding to next column
        const auto offsets = n + (s+1) * nlev;
        if (SDS.tke[offset] < SDS.tke[offsets] &
            SDS.isotropy[offset] < SDS.isotropy[offsets]){
          REQUIRE(SDS.tk[offset] < SDS.tk[offsets]);
          REQUIRE(SDS.tkh[offset] < SDS.tkh[offsets]);
        }
      }
    }
  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    EddyDiffusivitiesData f90_data[] = {
      EddyDiffusivitiesData(10, 71),
      EddyDiffusivitiesData(10, 12),
      EddyDiffusivitiesData(7,  16),
      EddyDiffusivitiesData(2, 7),
    };

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    EddyDiffusivitiesData cxx_data[] = {
      EddyDiffusivitiesData(f90_data[0]),
      EddyDiffusivitiesData(f90_data[1]),
      EddyDiffusivitiesData(f90_data[2]),
      EddyDiffusivitiesData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      eddy_diffusivities(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      eddy_diffusivities_f(d.nlev, d.shcol, d.obklen, d.pblh, d.zt_grid, d.shoc_mix, d.sterm_zt, d.isotropy, d.tke, d.tkh, d.tk);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(EddyDiffusivitiesData);
      for (Int i = 0; i < num_runs; ++i) {
        EddyDiffusivitiesData& d_f90 = f90_data[i];
        EddyDiffusivitiesData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.tkh); ++k) {
          REQUIRE(d_f90.tkh[k] == d_cxx.tkh[k]);
          REQUIRE(d_f90.tk[k] == d_cxx.tk[k]);
        }
      }
    }
  } // run_bfb
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_tke_eddy_diffusivities_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEddyDiff;

  TestStruct::run_property();
}

TEST_CASE("shoc_tke_eddy_diffusivities_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEddyDiff;

  TestStruct::run_bfb();
}

} // namespace
