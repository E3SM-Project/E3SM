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
struct UnitWrap::UnitTest<D>::TestCompShocMixLength {

  static void run_property()
  {
    static constexpr Int shcol    = 3;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC function:
    //   compute_shoc_mix_shoc_length

    // Multi column test will verify that 1) mixing length increases
    // with height given brunt vaisalla frequency and
    // TKE are constant with height and 2) Columns with larger
    // TKE values produce a larger length scale value

    // Define TKE [m2/s2] that will be used for each column
    static constexpr Real tke_cons = 0.1;
    // Define the brunt vasailla frequency [s-1]
    static constexpr Real brunt_cons = 0.001;
    // Define the assymptoic length [m]
    static constexpr Real l_inf = 100;
    // Define the heights on the zt grid [m]
    static constexpr Real zt_grid[nlev] = {5000, 3000, 2000, 1000, 500};

    // Initialize data structure for bridging to F90
    ComputeShocMixShocLengthData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    // For this test shcol MUST be at least 2
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    REQUIRE(shcol > 1);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      SDS.l_inf[s] = l_inf;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // for the subsequent columns, increase
        //  the amount of TKE
        SDS.tke[offset] = (1.0+s)*tke_cons;
        SDS.brunt[offset] = brunt_cons;
        SDS.zt_grid[offset] = zt_grid[n];
      }
    }

    // Check that the inputs make sense

    // Be sure that relevant variables are greater than zero
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.l_inf[s] > 0);
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.zt_grid[offset] > 0);
        if (s < shcol-1){
          // Verify that tke is larger column by column
          const auto offsets = n + (s+1) * nlev;
          REQUIRE(SDS.tke[offset] < SDS.tke[offsets]);
        }
      }

      // Check that zt increases upward
      for(Int n = 0; n < nlev - 1; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
      }
    }

    // Call the fortran implementation
    compute_shoc_mix_shoc_length(SDS);

    // Check the results
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Validate shoc_mix greater than zero everywhere
        REQUIRE(SDS.shoc_mix[offset] > 0);
        if (s < shcol-1){
          // Verify that mixing length increases column by column
          const auto offsets = n + (s+1) * nlev;
          REQUIRE(SDS.shoc_mix[offset] < SDS.shoc_mix[offsets]);
        }
      }

      // Check that mixing length increases upward
      for(Int n = 0; n < nlev - 1; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.shoc_mix[offset + 1] - SDS.shoc_mix[offset] < 0);
      }
    }
  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ComputeShocMixShocLengthData SDS_f90[] = {
      //               shcol, nlev
      ComputeShocMixShocLengthData(10, 71),
      ComputeShocMixShocLengthData(10, 12),
      ComputeShocMixShocLengthData(7,  16),
      ComputeShocMixShocLengthData(2, 7)
    };

    // Generate random input data
    for (auto& d : SDS_f90) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ComputeShocMixShocLengthData SDS_cxx[] = {
      ComputeShocMixShocLengthData(SDS_f90[0]),
      ComputeShocMixShocLengthData(SDS_f90[1]),
      ComputeShocMixShocLengthData(SDS_f90[2]),
      ComputeShocMixShocLengthData(SDS_f90[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      compute_shoc_mix_shoc_length(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      d.transpose<ekat::TransposeDirection::c2f>();
      // expects data in fortran layout
      compute_shoc_mix_shoc_length_f(d.nlev, d.shcol,
                                     d.tke, d.brunt,
                                     d.zt_grid,
                                     d.l_inf, d.shoc_mix);
      d.transpose<ekat::TransposeDirection::f2c>();
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(SDS_f90) / sizeof(ComputeShocMixShocLengthData);
      for (Int i = 0; i < num_runs; ++i) {
        ComputeShocMixShocLengthData& d_f90 = SDS_f90[i];
        ComputeShocMixShocLengthData& d_cxx = SDS_cxx[i];
        for (Int k = 0; k < d_f90.total(d_f90.shoc_mix); ++k) {
          REQUIRE(d_f90.shoc_mix[k] == d_cxx.shoc_mix[k]);
        }
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_mix_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompShocMixLength;

  TestStruct::run_property();
}

TEST_CASE("shoc_mix_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompShocMixLength;

  TestStruct::run_bfb();
}

} // namespace
