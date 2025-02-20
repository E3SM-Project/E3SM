#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/share/physics_constants.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
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
struct UnitWrap::UnitTest<D>::TestShocDiagObklen : public UnitWrap::UnitTest<D>::Base {

  void run_property()
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

    // Define bounds for reasonable output
    static constexpr Real obklen_lower = -10;
    static constexpr Real obklen_upper = 1000;
    static constexpr Real ustar_bound = 10;
    static constexpr Real kbfs_bound = 10;

    // Initialize data structure for bridging to F90
    ShocDiagObklenData SDS(shcol);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol == shcol);
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

    // Call the C++ implementation
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
      REQUIRE(std::abs(SDS.ustar[s]) < ustar_bound);
      REQUIRE(std::abs(SDS.kbfs[s]) < kbfs_bound);
      REQUIRE( (SDS.obklen[s] > obklen_lower && SDS.obklen[s] < obklen_upper));
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

    // Call the C++ implementation
    shoc_diag_obklen(SDS);

    // Verify that DIMENSIONLESS obukhov length decreases as columns
    //   increases due to the increasing surface fluxes

    // Define the lowest model layer grid height [m]
    Real zt_low = 50;
    for(Int s = 0; s < shcol-1; ++s) {
      REQUIRE(zt_low/SDS.obklen[s+1] < zt_low/SDS.obklen[s]);
    }

  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    ShocDiagObklenData SDS_baseline[] = {
      //             shcol
      ShocDiagObklenData(12),
      ShocDiagObklenData(10),
      ShocDiagObklenData(7),
      ShocDiagObklenData(2)
    };

    // Generate random input data
    for (auto& d : SDS_baseline) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    ShocDiagObklenData SDS_cxx[] = {
      ShocDiagObklenData(SDS_baseline[0]),
      ShocDiagObklenData(SDS_baseline[1]),
      ShocDiagObklenData(SDS_baseline[2]),
      ShocDiagObklenData(SDS_baseline[3])
    };

    static constexpr Int num_runs = sizeof(SDS_baseline) / sizeof(ShocDiagObklenData);

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : SDS_baseline) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      shoc_diag_obklen(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ShocDiagObklenData& d_baseline = SDS_baseline[i];
        ShocDiagObklenData& d_cxx = SDS_cxx[i];
        for (Int s = 0; s < d_baseline.shcol; ++s) {
          REQUIRE(d_baseline.ustar[s] == d_cxx.ustar[s]);
          REQUIRE(d_baseline.kbfs[s] == d_cxx.kbfs[s]);
          REQUIRE(d_baseline.obklen[s] == d_cxx.obklen[s]);
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

TEST_CASE("shoc_diag_obklen_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocDiagObklen;

  TestStruct().run_property();
}

TEST_CASE("shoc_diag_obklen_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocDiagObklen;

  TestStruct().run_bfb();
}

} // namespace
