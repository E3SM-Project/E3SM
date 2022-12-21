#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "ekat/util/ekat_arch.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestLInfShocLength {

  static void run_property()
  {
    static constexpr Int shcol    = 3;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC function:
    //   compute_l_inf_shoc_length

    // Multi-column test, where the input heights for zt_grid
    //  are increased uniformly by 100 m per column to verify
    //  that l_inf always gets larger per column.

    // Grid difference centered on thermo grid [m]
    static constexpr Real dz_zt[nlev] = {100, 100, 100, 100, 100};
    // Grid height centered on thermo grid [m]
    static constexpr Real zt_grid[nlev] = {500, 400, 300, 200, 100};
    // Turbulent kinetic energy [m2/s2]
    static constexpr Real tke[nlev] = {0.1, 0.15, 0.2, 0.25, 0.3};

    // Set upper bound for checking output [m]
    static constexpr Real l_inf_bound = 500;

    // Initialize data structure for bridgeing to F90
    ComputeLInfShocLengthData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    //  At least two columns are needed for this test!
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    REQUIRE(shcol > 1);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.dz_zt[offset] = dz_zt[n];
        // Testing identical columns but one with larger zt heights.
        //  first column set as "base" column, and the others
        //  to a larger value of TKE uniformly.
        SDS.zt_grid[offset] = (s*100.0)+zt_grid[n];
        SDS.tke[offset] = tke[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Be sure that relevant variables are greater than zero
        REQUIRE(SDS.dz_zt[offset] > 0);
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.zt_grid[offset] > 0);
        if (n < nlev-1){
          // check that zt increases upward
          REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
        }
        if (s < shcol-1){
          // Verify that zt_grid is offset larger column by column
          const auto offsets = n + (s+1) * nlev;
          REQUIRE(SDS.zt_grid[offset] < SDS.zt_grid[offsets]);
        }
      }
    }

    // Call the fortran implementation
    compute_l_inf_shoc_length(SDS);

    // Check the results
    // Make sure result is bounded correctly
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.l_inf[s] > 0);
      REQUIRE(SDS.l_inf[s] < l_inf_bound);
    }

    // Make sure that l_inf is getting larger
    //  per column
    for (Int s = 0; s < shcol-1; ++s){
      REQUIRE(SDS.l_inf[s] < SDS.l_inf[s+1]);
    }
  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ComputeLInfShocLengthData SDS_f90[] = {
      //               shcol, nlev
      ComputeLInfShocLengthData(10, 71),
      ComputeLInfShocLengthData(10, 12),
      ComputeLInfShocLengthData(7,  16),
      ComputeLInfShocLengthData(2, 7),
    };

    // Generate random input data
    for (auto& d : SDS_f90) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ComputeLInfShocLengthData SDS_cxx[] = {
      ComputeLInfShocLengthData(SDS_f90[0]),
      ComputeLInfShocLengthData(SDS_f90[1]),
      ComputeLInfShocLengthData(SDS_f90[2]),
      ComputeLInfShocLengthData(SDS_f90[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      compute_l_inf_shoc_length(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      d.transpose<ekat::TransposeDirection::c2f>();
      // expects data in fortran layout
      compute_l_inf_shoc_length_f(d.nlev,d.shcol,d.zt_grid,d.dz_zt,d.tke,d.l_inf);
      d.transpose<ekat::TransposeDirection::f2c>();
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(SDS_f90) / sizeof(ComputeLInfShocLengthData);
      for (Int i = 0; i < num_runs; ++i) {
        ComputeLInfShocLengthData& d_f90 = SDS_f90[i];
        ComputeLInfShocLengthData& d_cxx = SDS_cxx[i];
        for (Int c = 0; c < d_f90.shcol; ++c) {
          REQUIRE(d_f90.l_inf[c] == d_cxx.l_inf[c]);
        }
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_l_inf_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestLInfShocLength;

  TestStruct::run_property();
}

TEST_CASE("shoc_l_inf_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestLInfShocLength;

  TestStruct::run_bfb();
}

} // namespace
