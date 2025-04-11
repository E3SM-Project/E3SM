#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/eamxx_types.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

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
struct UnitWrap::UnitTest<D>::TestShocEddyDiff : public UnitWrap::UnitTest<D>::Base {

  void run_property()
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
    // Boundary layer regime test.  Input surface temperature with a value that is clearly
    //  "running away" (i.e. a very low value).  This test verifies that the diffusivities
    //  returned are different given that all other conditions are the same.  This is to very
    //  that SHOC can repair a runaway surface temperature.

    // Surface temperature [K]
    static constexpr Real tabs[shcol] = {100, 250};
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
      SDS.tabs[s] = tabs[s];
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

    // Call the C++ implementation
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
    // Runaway stable boundary layer test.  Given Obukhov length,
    //   verify that each regime behaves as expected when the relevant
    //   inputs are modified.  Column with larger mixing length and
    //   shear term should produce larger diffusivity values.

    // Surface temperature [K]
    static constexpr Real tabs_stab[shcol] = {100, 100};
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
      SDS.tabs[s] = tabs_stab[s];
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
      // Make sure we are testing a runaway stable boundary layer
      REQUIRE(SDS.tabs[s] < 180);
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Should be greater than zero
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.shoc_mix[offset] > 0);
        REQUIRE(SDS.isotropy[offset] > 0);
        REQUIRE(SDS.sterm_zt[offset] > 0);
      }
    }

    // Call the C++ implementation
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
    // Default boundary layer test.
    //   verify that each regime behaves as expected when the relevant
    //   inputs are modified for the default definitions of diffusivities.

    // Surface temperature [K]
    static constexpr Real tabs_ustab[shcol] = {300, 300};
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
      SDS.tabs[s] = tabs_ustab[s];
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
      // Make sure we are testing default boundary layer diffusivities
      REQUIRE(SDS.tabs[s] > 220);
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Should be greater than zero
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.shoc_mix[offset] > 0);
        REQUIRE(SDS.isotropy[offset] > 0);
        REQUIRE(SDS.sterm_zt[offset] > 0);
      }
    }

    // Call the C++ implementation
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


  void run_bfb()
  {
    auto engine = Base::get_engine();

    EddyDiffusivitiesData baseline_data[] = {
      EddyDiffusivitiesData(10, 71),
      EddyDiffusivitiesData(10, 12),
      EddyDiffusivitiesData(7,  16),
      EddyDiffusivitiesData(2, 7),
    };

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    EddyDiffusivitiesData cxx_data[] = {
      EddyDiffusivitiesData(baseline_data[0]),
      EddyDiffusivitiesData(baseline_data[1]),
      EddyDiffusivitiesData(baseline_data[2]),
      EddyDiffusivitiesData(baseline_data[3]),
    };

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      eddy_diffusivities(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(EddyDiffusivitiesData);
      for (Int i = 0; i < num_runs; ++i) {
        EddyDiffusivitiesData& d_baseline = baseline_data[i];
        EddyDiffusivitiesData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.tkh); ++k) {
          REQUIRE(d_baseline.tkh[k] == d_cxx.tkh[k]);
          REQUIRE(d_baseline.tk[k] == d_cxx.tk[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (auto& d : cxx_data) {
        d.write(Base::m_fid);
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

  TestStruct().run_property();
}

TEST_CASE("shoc_tke_eddy_diffusivities_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEddyDiff;

  TestStruct().run_bfb();
}

} // namespace
