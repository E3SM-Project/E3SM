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
struct UnitWrap::UnitTest<D>::TestShocEnergyInt : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC function
    //     shoc_energy_integrals

    // FIRST TEST and SECOND TEST
    //   FIRST: kinetic energy test.  Given at least two columns with
    //   identical inputs but with increasing winds in each column
    //   verify that the column with stronger winds results
    //   in a higher kinetic energy test.
    //   SECOND: vapor test.  verify that columns with no vapor
    //   or moisture result in zero integral outputs.  If the
    //   column is positive then verify that wl_int >= wv_int

    // Define host model dry static energy [J kg-1]
    static constexpr Real host_dse[nlev] = {350e3, 325e3, 315e3, 310e3, 300e3};
    // Defin the pressure difference [Pa]
    static constexpr Real pdel[nlev]={100e2, 75e2, 50e2, 25e2, 10e2};
    // Define zonal wind on nlev grid [m/s]
    static constexpr Real u_wind[nlev] = {-4, -4, -3, -1, -2};
    // Define meridional wind on nlev grid [m/s]
    static constexpr Real v_wind[nlev] = {1, 2, 3, 4, 5};
    // Define the total water mixing ratio [kg/kg]
    static constexpr Real rtm[nlev] = {7e-3, 9e-3, 11e-3, 15e-3, 20e-3};
    // Define cloud mixing ratio [g/kg] (converted to kg/kg later)
    static constexpr Real rcm[nlev] = {1e-6, 1e-5, 4e-5, 2e-6, 1e-6};

    // Initialize data structure for bridging to F90
    ShocEnergyIntegralsData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    // for this test we need exactly two columns
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    REQUIRE(shcol == 2);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // Add one degree K in the second column
        SDS.host_dse[offset] = s+host_dse[n];

        // Force the first column of cloud liquid
        //   to be zero!
        SDS.rcm[offset] = s*rcm[n];
        SDS.rtm[offset] = rtm[n];

        // convert to Pa
        SDS.pdel[offset] = pdel[n];

        // Increase winds with increasing columns
        SDS.u_wind[offset] = (1.0+s)*u_wind[n];
        SDS.v_wind[offset] = (1.0+s)*v_wind[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;

        REQUIRE(SDS.host_dse[offset] > 0);
        REQUIRE(SDS.rcm[offset] >= 0);
        REQUIRE(SDS.rtm[offset] > 0);

        // make sure the two columns are different and
        //  as expected for the relevant variables
        if (s == 0){
          const auto offsets = n + (s+1) * nlev;

          REQUIRE(abs(SDS.u_wind[offsets]) > abs(SDS.u_wind[offset]));
          REQUIRE(abs(SDS.v_wind[offsets]) > abs(SDS.v_wind[offset]));
          REQUIRE(SDS.rcm[offset] == 0);
          REQUIRE(SDS.rcm[offsets] > SDS.rcm[offset]);
        }
      }
    }

    // Call the C++ implementation
    shoc_energy_integrals(SDS);

    // Check test
    for(Int s = 0; s < shcol; ++s) {
      // Verify relevant integrals are reasonable
      REQUIRE(SDS.se_int[s] > 0);
      REQUIRE(SDS.ke_int[s] > 0);
      REQUIRE(SDS.wv_int[s] > 0);
      REQUIRE(SDS.wl_int[s] >= 0);
      // Do column comparison tests
      if (s == 0){
        REQUIRE(SDS.ke_int[s+1] > SDS.ke_int[s]);
        REQUIRE(SDS.wl_int[s+1] > SDS.wl_int[s]);
        REQUIRE(SDS.wl_int[s] == 0);
        REQUIRE(SDS.se_int[s+1] > SDS.se_int[s]);
      }
    }
  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    ShocEnergyIntegralsData SDS_baseline[] = {
      //               shcol, nlev
      ShocEnergyIntegralsData(10, 71),
      ShocEnergyIntegralsData(10, 12),
      ShocEnergyIntegralsData(7,  16),
      ShocEnergyIntegralsData(2, 7),
    };

    // Generate random input data
    for (auto& d : SDS_baseline) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    ShocEnergyIntegralsData SDS_cxx[] = {
      ShocEnergyIntegralsData(SDS_baseline[0]),
      ShocEnergyIntegralsData(SDS_baseline[1]),
      ShocEnergyIntegralsData(SDS_baseline[2]),
      ShocEnergyIntegralsData(SDS_baseline[3]),
    };

    static constexpr Int num_runs = sizeof(SDS_baseline) / sizeof(ShocEnergyIntegralsData);

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : SDS_baseline) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      shoc_energy_integrals(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ShocEnergyIntegralsData& d_baseline = SDS_baseline[i];
        ShocEnergyIntegralsData& d_cxx = SDS_cxx[i];
        for (Int c = 0; c < d_baseline.shcol; ++c) {
          REQUIRE(d_baseline.se_int[c] == d_cxx.se_int[c]);
          REQUIRE(d_baseline.ke_int[c] == d_cxx.ke_int[c]);
          REQUIRE(d_baseline.wv_int[c] == d_cxx.wv_int[c]);
          REQUIRE(d_baseline.wl_int[c] == d_cxx.wl_int[c]);
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

TEST_CASE("shoc_energy_integrals_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyInt;

  TestStruct().run_property();
}

TEST_CASE("shoc_energy_integrals_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyInt;

  TestStruct().run_bfb();
}

} // namespace
