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
struct UnitWrap::UnitTest<D>::TestShocDiagThird : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr Int nlevi    = nlev+1;

    // Tests for the SHOC function:
    //   diag_third_shoc_moments

    // Upper level test for w3

    // Convective boundary layer test with extreme gradients
    // Most physics of the test have been applied in midlayer function.
    //  Here we give the profile very small dz and hence extreme gardients
    //  to be sure that w3 is clipped appropriately.

    // Then we will increase w_sec and tke for a second test to verify that
    //  the magnitude of w3 has increased.

    // Define vertical velocity second moment [m2/s2]
    static constexpr Real w_sec[nlev] = {0.2, 0.3, 0.5, 0.3, 0.1};
    // Define potential temperature second moment [K2]
    static constexpr Real thl_sec[nlevi] = {0.5, 0.9, 1.2, 0.8, 0.4, 0.3};
    // Define vertical flux of temperature [K m/s]
    static constexpr Real wthl_sec[nlevi] = {0.003, -0.03, -0.04, -0.01, 0.01, 0.03};
    // Define the heights on the zi grid [m]
    static constexpr Real zi_grid[nlevi] = {900, 500, 150, 90, 50, 0};
    // Define the return to isotropy timescale [s]
    static constexpr Real isotropy[nlev] = {1000, 2000, 4000, 1000, 500};
    // Define the brunt vaisalla frequency
    static constexpr Real brunt[nlev] = {4e-5, 3e-5, 2e-5, 2e-5, -1e-5};
    // Define the potential temperature on zi grid [K]
    static constexpr Real thetal[nlev] = {330, 320, 310, 300, 305};

    // Define TKE [m2/s2], compute from w_sec
    Real tke[nlev];

    // Grid stuff to compute based on zi_grid
    Real zt_grid[nlev];
    Real dz_zt[nlev];
    Real dz_zi[nlevi];
    // Compute heights on midpoint grid
    for(Int n = 0; n < nlev; ++n) {
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
      tke[n] = 1.5*w_sec[n];
      dz_zt[n] = zi_grid[n] - zi_grid[n+1];
      if (n == 0){
        dz_zi[n] = 0;
      }
      else{
        dz_zi[n] = zt_grid[n-1] - zt_grid[n];
      }
    }
    // set upper condition for dz_zi
    dz_zi[nlevi-1] = zt_grid[nlev-1];

    // Initialize data structure for bridging to F90
    DiagThirdShocMomentsData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable.
    // For this test shcol MUST be at least 2
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(SDS.nlevi == SDS.nlev+1);

    // Load up the new data
    for(Int s = 0; s < shcol; ++s) {
      // Fill in test data on zt_grid.
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.w_sec[offset] = w_sec[n];
        SDS.dz_zt[offset] = dz_zt[n];
        SDS.zt_grid[offset] = zt_grid[n];
        SDS.tke[offset] = tke[n];

        SDS.isotropy[offset] = isotropy[n];
        SDS.brunt[offset] = brunt[n];
        SDS.thetal[offset] = thetal[n];
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.dz_zi[offset] = dz_zi[n];
        SDS.zi_grid[offset] = zi_grid[n];
        SDS.thl_sec[offset] = thl_sec[n];
        SDS.wthl_sec[offset] = wthl_sec[n];

      }
    }

    // Check that the inputs make sense
    // Load up the new data
    for(Int s = 0; s < shcol; ++s) {
      // Fill in test data on zt_grid.
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.w_sec[offset] >= 0);
        // for this test make sure w_sec <= 1
        REQUIRE(SDS.w_sec[offset] <= 1);

        REQUIRE(SDS.dz_zt[offset] > 0);
        REQUIRE(SDS.zt_grid[offset] > 0);
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.isotropy[offset] >= 0);
        REQUIRE(SDS.thetal[offset] >= 0);
      }

      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        REQUIRE(SDS.dz_zi[offset] >= 0);
        REQUIRE(SDS.zi_grid[offset] >= 0);
        REQUIRE(SDS.thl_sec[offset] >= 0);
      }
    }

    // Call the C++ implementation
    diag_third_shoc_moments(SDS);

    // Check to make sure there is at least one
    //  positive w3 value for convective boundary layer
    bool is_skew;
    // Verify that boundary points are zero
    for(Int s = 0; s < shcol; ++s) {
      is_skew = false; // Initialize
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        // For this test make sure w3 has been clipped.
        //  Given input w2, this should be less than 1
        REQUIRE(abs(SDS.w3[offset]) < 1);

        if (SDS.w3[offset] > 0){
          is_skew = true;
        }
      }
      REQUIRE(is_skew == true);
    }

    // Now save the result from first column
    Real w3_test1[nlevi*shcol];
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        w3_test1[offset] = SDS.w3[offset];
      }
    }

    // Load up new and increased TKE values
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlevi;

        SDS.w_sec[offset] = 10*w_sec[n];
        // update new TKE value
        SDS.tke[offset] = 1.5*SDS.w_sec[offset];
      }
    }

    // Call the C++ implementation
    diag_third_shoc_moments(SDS);

    // Verify that new result is greater or equal in magnitude
    //  that the result from test one
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
        if (n != 0 && n != nlevi-1){
          REQUIRE(std::abs(SDS.w3[offset]) >= std::abs(w3_test1[offset]));
        }
      }
    }

  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    DiagThirdShocMomentsData SDS_baseline[] = {
      //               shcol, nlev, nlevi
      DiagThirdShocMomentsData(10, 71, 72),
      DiagThirdShocMomentsData(10, 12, 13),
      DiagThirdShocMomentsData(7,  16, 17),
      DiagThirdShocMomentsData(2, 7, 8),
    };

    // Generate random input data
    for (auto& d : SDS_baseline) {
      d.randomize(engine, {{d.thetal, {300, 301}}});
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    DiagThirdShocMomentsData SDS_cxx[] = {
      DiagThirdShocMomentsData(SDS_baseline[0]),
      DiagThirdShocMomentsData(SDS_baseline[1]),
      DiagThirdShocMomentsData(SDS_baseline[2]),
      DiagThirdShocMomentsData(SDS_baseline[3]),
    };

    static constexpr Int num_runs = sizeof(SDS_baseline) / sizeof(DiagThirdShocMomentsData);

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : SDS_baseline) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      diag_third_shoc_moments(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        DiagThirdShocMomentsData& d_baseline = SDS_baseline[i];
        DiagThirdShocMomentsData& d_cxx = SDS_cxx[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.w3); ++k) {
          REQUIRE(d_baseline.w3[k] == d_cxx.w3[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        SDS_cxx[i].write(Base::m_fid);
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_diag_third_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocDiagThird;

  TestStruct().run_property();
}

TEST_CASE("shoc_diag_third_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocDiagThird;

  TestStruct().run_bfb();
}

} // namespace
