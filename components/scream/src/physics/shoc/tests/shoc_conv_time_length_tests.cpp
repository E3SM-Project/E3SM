#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
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
struct UnitWrap::UnitTest<D>::TestCompShocConvTime {

  static void run_property()
  {
    static constexpr Int shcol    = 5;

    // Tests for the SHOC function:
    //   compute_conv_time_shoc_length

    // FIRST TEST
    //  Verify that no matter the input of conv_vel that
    //  time scale is always return as positive.
    //  The main concern is that if input conv_vel is negative that
    //  we should avoid taking the cube root of this value in the
    //  subroutine call.

    //  For this function we test several one layer columns

    // PBL height [m]
    static constexpr Real pblh[shcol] = {1000, 400, 10, 500, 300};
    // Integrated convective velocity [m3/s3]
    static constexpr Real conv_vel[shcol] = {10, -3.5, 0.1, -100, -0.4};
    // Upper bound of expected result [s]
    static constexpr Real tscale_upper = 10000;

    // Initialize data structure for bridging to F90
    SHOCConvtimeData SDS(shcol);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol() == shcol);
    REQUIRE(shcol > 0);

    // Fill in test data, all one dimensional
    for(Int s = 0; s < shcol; ++s) {
      SDS.pblh[s] = pblh[s];
      SDS.conv_vel[s] = conv_vel[s];
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] > 0.0);
    }

    // Call the fortran implementation
    compute_conv_time_shoc_length(SDS);

    // Check the results
    // Make sure that timescale is positive always
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.tscale[s] > 0.0);
      REQUIRE(SDS.tscale[s] < tscale_upper);
    }

    // SECOND TEST
    // Constant PBL depth test
    // Given a set of inputs where PBL depth is constant but
    //  convective velocity scale changes, verify that the
    //  turn over time scale is inversely proportional to
    //  convective scale.  All values of convective scale
    //  must be positive for this test.

    // PBL height [m]
    static constexpr Real pblh_cons1[shcol] = {500.0, 500.0, 500.0, 500.0, 500.0};
    // Integrated convective velocity [m3/s3]
    static constexpr Real conv_vel_cons1[shcol] = {10.0, 3.5, 0.1, 100.0, 0.4};

    // Fill in test data, all one dimensional
    for(Int s = 0; s < shcol; ++s) {
      SDS.pblh[s] = pblh_cons1[s];
      SDS.conv_vel[s] = conv_vel_cons1[s];
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] > 0.0);
      // For this test all values of conv_vel must be positive
      REQUIRE(SDS.conv_vel[s] > 0.0);
    }

    // Call the fortran implementation
    compute_conv_time_shoc_length(SDS);

    // Verify Results
    // Make sure that if conv_vel is smaller than in a neighboring
    //  column, that the tscale is larger
    for (Int s = 0; s < shcol-1; ++s){
      if (SDS.tscale[s] > SDS.tscale[s+1]){
        REQUIRE(SDS.conv_vel[s] < SDS.conv_vel[s+1]);
      }
    }

    // Make sure bounds fall in reasonable range
    for (Int s = 0; s < shcol; ++s){
      REQUIRE(SDS.tscale[s] < tscale_upper);
      REQUIRE(SDS.tscale[s] > 0);
    }

    // THIRD TEST
    // Constant conv_vel test
    // Given a set of inputs where conv_vel is constant but
    //  boundary layer depth changes, verify that the
    //  turn over time scale is proportional to
    //  PBL height.  All values of convective scale
    //  must be positive for this test.

    // PBL height [m]
    static constexpr Real pblh_cons2[shcol] = {500.0, 100.0, 10.0, 5000.0, 750.0};
    // Integrated convective velocity [m3/s3]
    static constexpr Real conv_vel_cons2[shcol] = {0.5, 0.5, 0.5, 0.5, 0.5};

    // Fill in test data, all one dimensional
    for(Int s = 0; s < shcol; ++s) {
      SDS.pblh[s] = pblh_cons2[s];
      SDS.conv_vel[s] = conv_vel_cons2[s];
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] > 0.0);
      // For this test all values of conv_vel must be positive
      REQUIRE(SDS.conv_vel[s] > 0.0);
    }

    // Call the fortran implementation
    compute_conv_time_shoc_length(SDS);

    // Verify Results
    // Make sure that if tscale is larger than in a
    //  neighboring column, that the PBL depth is larger
    for (Int s = 0; s < shcol-1; ++s){
      if (SDS.tscale[s] > SDS.tscale[s+1]){
        REQUIRE(SDS.pblh[s] > SDS.pblh[s+1]);
      }
    }

    // Make sure bounds fall in reasonable range
    for (Int s = 0; s < shcol; ++s){
      REQUIRE(SDS.tscale[s] < tscale_upper);
      REQUIRE(SDS.tscale[s] > 0);
    }

  }

  static void run_bfb()
  {
    // TODO
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_conv_time_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompShocConvTime;

  TestStruct::run_property();
}

TEST_CASE("shoc_conv_time_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompShocConvTime;

  TestStruct::run_bfb();
}

} // namespace
