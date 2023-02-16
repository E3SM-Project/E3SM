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
struct UnitWrap::UnitTest<D>::TestShocIntColStab {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;

    // Tests for the subroutine integ_column_stability
    //   in the SHOC TKE module.

    // FIRST TEST
    // Symmetric test.  Verify that given a profile
    //  of symetric inputs for brunt vaisalla and dz profile
    //  that is uniform with height that the column integrated
    //  value is 0.0.

    // Define height thickness on nlev grid [m]
    //   Do a uniform grid for symetric test
    static constexpr Real dz_zt[nlev] = {30, 30, 30, 30, 30};
    // Define Pressure [Pa]
    Real pres[nlev] = {850e2, 900e2, 925e2, 950e2, 1000e2};
    // Brunt Vaisalla frequency [/s]
    static constexpr Real brunt_sym[nlev] = {-0.5, -0.25, 0.0, 0.25, 0.5};

    // Initialize data structure for bridging to F90
    IntegColumnStabilityData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    REQUIRE(shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.dz_zt[offset] = dz_zt[n];
        SDS.pres[offset] = pres[n];
        SDS.brunt[offset] = brunt_sym[n];
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Should be greater than zero
        REQUIRE(SDS.dz_zt[offset] > 0);
        // Make sure all pressure levels are in the
        //  lower troposphere for this test
        REQUIRE(SDS.pres[offset] > 80000);
      }
    }

    // Call the fortran implementation
    integ_column_stability(SDS);

    // Check test
    //  Verify that output is zero
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.brunt_int[s] == 0);
    }

    // SECOND TEST
    // For set of inputs where brunt is negative at all
    //   points, then verify that output is negative

    // Brunt Vaisalla frequency [/s]
    static constexpr Real brunt_neg[nlev] = {-0.3, -0.4, -0.1, -10.0, -0.5};

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        SDS.brunt[offset] = brunt_neg[n];
      }
    }

    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // All points should be less than zero
        REQUIRE(SDS.brunt[offset] < 0);
      }
    }

    // Call the fortran implementation
    integ_column_stability(SDS);

    // Check test
    //  Verify that output is negative
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.brunt_int[s] < 0);
    }
  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    //declare data for the f90 function call
    IntegColumnStabilityData f90_data[] = {
      IntegColumnStabilityData(10, 71),
      IntegColumnStabilityData(10, 12),
      IntegColumnStabilityData(7,  16),
      IntegColumnStabilityData(2,  7 ),
    };

    //Generate random data
    for (auto &d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data (if any) is in original state
    IntegColumnStabilityData cxx_data[] = {
      IntegColumnStabilityData(f90_data[0]),
      IntegColumnStabilityData(f90_data[1]),
      IntegColumnStabilityData(f90_data[2]),
      IntegColumnStabilityData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto &d : f90_data) {
      // expects data in C layout
      integ_column_stability(d);
    }

    // Get data from cxx
    for (auto &d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      integ_column_stability_f(d.nlev, d.shcol, d.dz_zt, d.pres, d.brunt, d.brunt_int);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(IntegColumnStabilityData);
      for (Int i = 0; i < num_runs; ++i) {
        IntegColumnStabilityData& d_f90 = f90_data[i];
        IntegColumnStabilityData& d_cxx = cxx_data[i];
        for (Int c = 0; c < d_f90.shcol; ++c) {
          REQUIRE(d_f90.brunt_int[c] == d_cxx.brunt_int[c]);
        }
      }
    }
  } //run_bfb
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_tke_column_stab_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocIntColStab;

  TestStruct::run_property();
}

TEST_CASE("shoc_tke_column_stab_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocIntColStab;

  TestStruct::run_bfb();
}

} // namespace
