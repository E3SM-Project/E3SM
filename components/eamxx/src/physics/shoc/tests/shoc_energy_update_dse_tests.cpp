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
struct UnitWrap::UnitTest<D>::TestShocUpdateDse : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC function
    //     update_host_dse

    // FIRST TEST and SECOND TEST
    //  1) Verify that DSE increases with height and 2) verify that
    //   points that contain cloud condensate result in a larger
    //   dse if all other inputs are equal.

    // Liquid water potential temperature [K]
    static constexpr Real thlm[nlev] = {350, 325, 315, 310, 300};
    // Inverse of the Exner function [-]
    static constexpr Real inv_exner[nlev] = {1.22, 1.15, 1.10, 1.05, 1.0};
    // Cloud condensate [kg/kg]
    static constexpr Real shoc_ql[nlev] = {5e-6, 8e-5, 4e-4, 3e-4, 1e-6};
    // Mid-point heights [m]
    static constexpr Real zt_grid[nlev] = {9000, 6000, 4000, 2000, 1000};
    // Surface geopotential
    static constexpr Real phis = 100;

    // Set lower/upper bounds to check output DSE
    static constexpr Real dse_upper = 5e5; // [J/kg]
    static constexpr Real dse_lower = 1e5; // [J/kg]

    // Initialize data structure for bridging to F90
    UpdateHostDseData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    // for this test we need exactly two columns
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    REQUIRE(shcol == 2);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      SDS.phis[s] = phis;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.thlm[offset] = thlm[n];
        SDS.zt_grid[offset] = zt_grid[n];
        SDS.inv_exner[offset] = inv_exner[n];

        // Force the first column of cloud liquid
        //   to be zero!
        SDS.shoc_ql[offset] = s*shoc_ql[n];

      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.phis[s] >= 0.0);
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;

        REQUIRE(SDS.thlm[offset] > 0);
        REQUIRE(SDS.shoc_ql[offset] >= 0);
        REQUIRE(SDS.inv_exner[offset] > 0);
        REQUIRE(SDS.zt_grid[offset] >= 0);

        // make sure the two columns are different and
        //  as expected for the relevant variables
        if (s == 0){
          const auto offsets = n + (s+1) * nlev;

          REQUIRE(SDS.shoc_ql[offset] == 0);
          REQUIRE(SDS.shoc_ql[offsets] > SDS.shoc_ql[offset]);
        }

        // Check that heights, thlm, and inverse of exner increase upward
        if (n < nlev-1){
          REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
          REQUIRE(SDS.thlm[offset + 1] - SDS.thlm[offset] < 0);
          REQUIRE(SDS.inv_exner[offset + 1] - SDS.inv_exner[offset] < 0);
        }

      }
    }

    // Call the C++ implementation
    update_host_dse(SDS);

    // Check test
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Verify dse is reasonable
        REQUIRE(SDS.host_dse[offset] > 0);

        // Verify that dse increases with height upward
        if (n < nlev-1){
          REQUIRE(SDS.host_dse[offset + 1] - SDS.host_dse[offset] < 0);
        }

        // Verify that DSE falls within some reasonable bounds
        REQUIRE(SDS.host_dse[offset] > dse_lower);
        REQUIRE(SDS.host_dse[offset] < dse_upper);

        // Verify that dse is greater in points with condensate loading
        if (s == 0){
          const auto offsets = n + (s+1) * nlev;
          REQUIRE(SDS.host_dse[offsets] > SDS.host_dse[offset]);
        }
      }
    }
  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    UpdateHostDseData SDS_baseline[] = {
      //               shcol, nlev
      UpdateHostDseData(10, 71),
      UpdateHostDseData(10, 12),
      UpdateHostDseData(7,  16),
      UpdateHostDseData(2, 7),
    };

    // Generate random input data
    for (auto& d : SDS_baseline) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    UpdateHostDseData SDS_cxx[] = {
      UpdateHostDseData(SDS_baseline[0]),
      UpdateHostDseData(SDS_baseline[1]),
      UpdateHostDseData(SDS_baseline[2]),
      UpdateHostDseData(SDS_baseline[3]),
    };

    static constexpr Int num_runs = sizeof(SDS_baseline) / sizeof(UpdateHostDseData);

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : SDS_baseline) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      update_host_dse(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        UpdateHostDseData& d_baseline = SDS_baseline[i];
        UpdateHostDseData& d_cxx = SDS_cxx[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.host_dse); ++k) {
          REQUIRE(d_baseline.host_dse[k] == d_cxx.host_dse[k]);
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

TEST_CASE("shoc_energy_host_dse_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocUpdateDse;

  TestStruct().run_property();
}

TEST_CASE("shoc_energy_host_dse_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocUpdateDse;

  TestStruct().run_bfb();
}

} // namespace
