#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"

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
struct UnitWrap::UnitTest<D>::TestCompShocConvVel {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static const auto approx_zero = Approx(0.0).margin(1e-16);

    // Tests for the SHOC function:
    //   compute_conv_vel_shoc_length

    // Note that the output conv_vel from this subroutine
    //   represents the integrated buoyancy flux and not technically
    //   the convective velocity scale, which is the cubed root
    //   of conv_vel.

    // FIRST TEST
    // Buoyancy test.  Verify that given positively buoyant conditions
    //   in one column and negatively buoyant conditions in another
    //   that the result is as expected (positive and negative,
    //   respectively).

    // PBL height [m]
    // for this test set to be higher than highest zt_grid value
    static constexpr Real pblh = 1000;
    // Grid difference centered on thermo grid [m]
    static constexpr Real dz_zt[nlev] = {100, 100, 100, 100, 100};
    // Grid height centered on thermo grid [m]
    static constexpr Real zt_grid[nlev] = {500, 400, 300, 200, 100};
    // Virtual potential temperature on interface grid [K]
    static constexpr Real thv[nlev] = {310, 305, 300, 300, 295};
    // Buoyancy flux [K m/s]
    static constexpr Real wthv_sec[nlev] = {0.02, 0.01, 0.04, 0.02, 0.05};

    // Bounds to check result [m/s]
    static constexpr Real conv_vel_bound = 10;

    // Initialize data structure for bridging to F90
    SHOCConvvelData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol() == shcol && SDS.nlev() == nlev) );
    // For this test we want exactly two columns
    REQUIRE(shcol == 2);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      SDS.pblh[s] = pblh;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.dz_zt[offset] = dz_zt[n];
        SDS.zt_grid[offset] = zt_grid[n];
        SDS.thv[offset] = thv[n];
        SDS.wthv_sec[offset] = wthv_sec[n];
        // for the second column, give negative values
        if (s == 1){
          SDS.wthv_sec[offset] = -1*SDS.wthv_sec[offset];
        }
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Be sure that relevant variables are greater than zero
        REQUIRE(SDS.dz_zt[offset] > 0);
        REQUIRE(SDS.thv[offset] > 0);
        // Make sure sign of buoyancy flux is as expected
        if (s == 0){
          REQUIRE(SDS.wthv_sec[offset] > 0);
        }
        else{
          REQUIRE(SDS.wthv_sec[offset] < 0);
        }
        if (n < nlev-1){
          // check that zt increases upward
          REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
        }
      }
    }

    // Call the fortran implementation
    compute_conv_vel_shoc_length(SDS);

    // Check the results
    // Make sure that conv_vel is positive in first
    //  column and negative in the second
    REQUIRE(SDS.conv_vel[0] > 0);
    REQUIRE(SDS.conv_vel[1] < 0);

    // Make sure result falls within reasonable bounds
    for (Int s = 0; s < shcol; ++s){
      REQUIRE(SDS.conv_vel[s] < std::abs(conv_vel_bound));
    }

    // SECOND TEST
    // Symmetrical input test.  Here we feed input values of wthv_sec
    //  that are symmetrical around zero.  This is to verify that the
    //  result of conv_vel is also zero.  Here we need to assume
    //  a well mixed layer and constant profile of dz_zt

    static constexpr Real dz_zt_sym[nlev] = {100, 100, 100, 100, 100};
    // Virtual potential temperature on interface grid [K]
    static constexpr Real thv_sym[nlev] = {300, 300, 300, 300, 300};
    // Buoyancy flux [K m/s]
    static constexpr Real wthv_sec_sym[nlev] = {0.04, 0.02, 0.00, -0.02, -0.04};

    // Fill in new test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.dz_zt[offset] = dz_zt_sym[n];
        SDS.thv[offset] = thv_sym[n];
        SDS.wthv_sec[offset] = wthv_sec_sym[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      Real wthv_sum = 0;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Be sure that relevant variables are greater than zero
        REQUIRE(SDS.dz_zt[offset] > 0);
        REQUIRE(SDS.thv[offset] > 0);
        wthv_sum += SDS.wthv_sec[offset];
      }
      // Make sure inputs of buoyancy flux sum to zero
      REQUIRE(wthv_sum == approx_zero);
    }

    // Call the fortran implementation
    compute_conv_vel_shoc_length(SDS);

    // Check the results
    // Make sure that conv_vel is zero
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.conv_vel[s] == approx_zero);
    }
  }

  static void run_bfb()
  {
    SHOCConvvelData SDS_f90[] = {
      //           shcol, nlev
      SHOCConvvelData(10, 71),
      SHOCConvvelData(10, 12),
      SHOCConvvelData(7,  16),
      SHOCConvvelData(2, 7)
    };

    static constexpr Int num_runs = sizeof(SDS_f90) / sizeof(SHOCConvvelData);

    // Generate random input data
    for (auto& d : SDS_f90) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    SHOCConvvelData SDS_cxx[] = {
      SHOCConvvelData(SDS_f90[0]),
      SHOCConvvelData(SDS_f90[1]),
      SHOCConvvelData(SDS_f90[2]),
      SHOCConvvelData(SDS_f90[3])
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      compute_conv_vel_shoc_length(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      d.transpose<ekat::TransposeDirection::c2f>();
      // expects data in fortran layout
      compute_conv_vel_shoc_length_f(d.nlev(),d.shcol(),d.pblh,d.zt_grid,
                                     d.dz_zt,d.thv,d.wthv_sec,d.conv_vel);
      d.transpose<ekat::TransposeDirection::f2c>();
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      SHOCConvvelData& d_f90 = SDS_f90[i];
      SHOCConvvelData& d_cxx = SDS_cxx[i];
      for (Int c = 0; c < d_f90.dim1; ++c) {
        REQUIRE(d_f90.conv_vel[c] == d_cxx.conv_vel[c]);
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_conv_vel_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompShocConvVel;

  TestStruct::run_property();
}

TEST_CASE("shoc_conv_vel_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompShocConvVel;

  TestStruct::run_bfb();
}

} // namespace
