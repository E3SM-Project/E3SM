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
struct UnitWrap::UnitTest<D>::TestImpCompTmpi : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlevi    = 6;

    // Tests for the SHOC subroutine
    //   compute_tmpi

    // TEST
    // Load up two columns, one with smaller dz.  The case with
    //  smaller dz should return a HIGHER value of tmpi

    // Define height thickness on nlevi grid [m]
    //   NOTE: First indicee is zero because it is never used
    //   Do a stretched grid
    static constexpr Real dz_zi[nlevi] = {0, 500, 200, 100, 50, 10};
    // Define density on zi grid [kg/m3]
    static constexpr Real rho_zi[nlevi] = {0.5, 0.6, 0.8, 0.9, 1.0, 1.2};
    // Define timestep [s]
    static constexpr Real dtime = 300;

    // Initialize data structure for bridging to F90
    ComputeTmpiData SDS(shcol, nlevi, dtime);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlevi == nlevi) );
    REQUIRE(shcol == 2);
    // Need exactly two columns for this test
    REQUIRE(SDS.dtime > 0);

    // Fill in test data on zi_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        // Feed second column SMALLER dz values
        SDS.dz_zi[offset] = dz_zi[n]/(1+s);
        SDS.rho_zi[offset] = rho_zi[n];
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlevi; ++n){
        const auto offset = n + s * nlevi;
        const auto offsets = n + (1+s)*nlevi;
        // make sure that density is in reasonable bounds
        REQUIRE( (SDS.rho_zi[offset] > 0 && SDS.rho_zi[offset] < 1.5) );
        // Make sure top level dz_zi value is zero
        if (n == 0){
          REQUIRE(SDS.dz_zi[offset] == 0);
        }
        // Otherwise, should be greater than zero
        else{
          REQUIRE(SDS.dz_zi[offset] > 0);
          // Verify that the second column has smaller dz values
          if (s < shcol-1){
            REQUIRE(SDS.dz_zi[offset] > SDS.dz_zi[offsets]);
          }
        }
      }
    }

    // Call the C++ implementation
    compute_tmpi(SDS);

    // Verify result
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlevi; ++n){
        const auto offset = n + s * nlevi;
        const auto offsets = n + (1+s)*nlevi;

        // Verify that highest model point is ALWAYS zero
        if (n == 0){
          REQUIRE(SDS.tmpi[offset] == 0);
        }
        else{
          // Verify bounds are reasonable
          REQUIRE(SDS.tmpi[offset] > 0);
          // By dimensional analysis, this value should really
          //  not be greater than an order of magnitude higher than
          //  the timestep, with reasonable inputs
          REQUIRE(SDS.tmpi[offset] < dtime*10);
          // Verify that tmpi values are larger when dz is smaller
          if (s < shcol-1){
            REQUIRE(SDS.tmpi[offset] < SDS.tmpi[offsets]);
          }
        }
      }
    }

  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    ComputeTmpiData baseline_data[] = {
      //          shcol, nlevi, dtime
      ComputeTmpiData(10, 72, 1),
      ComputeTmpiData(10, 13, 10),
      ComputeTmpiData(7,  17, 125),
      ComputeTmpiData(2,   8, 300)
    };

    // Generate random input data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    ComputeTmpiData cxx_data[] = {
      ComputeTmpiData(baseline_data[0]),
      ComputeTmpiData(baseline_data[1]),
      ComputeTmpiData(baseline_data[2]),
      ComputeTmpiData(baseline_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ComputeTmpiData);

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      compute_tmpi(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ComputeTmpiData& d_baseline = baseline_data[i];
        ComputeTmpiData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.tmpi); ++k) {
          REQUIRE(d_baseline.tmpi[k] == d_cxx.tmpi[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        cxx_data[i].write(Base::m_fid);
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_imp_comp_tmpi_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestImpCompTmpi;

  TestStruct().run_property();
}

TEST_CASE("shoc_imp_comp_tmpi_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestImpCompTmpi;

  TestStruct().run_bfb();
}

} // namespace
