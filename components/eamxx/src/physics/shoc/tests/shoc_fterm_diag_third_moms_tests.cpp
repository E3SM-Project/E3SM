#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
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
struct UnitWrap::UnitTest<D>::TestFtermdiagThirdMoms {

  static void run_property()
  {

    // Tests for the SHOC function:
    //   f0_to_f5_diag_third_shoc_moment

    // TEST ONE
    // Zero test.  Given no gradients, verify that relevant
    //  terms are zero.

    // 1/grid spacing [m-1]
    constexpr static Real thedz = 0.1;
    // 1/grid spacing for two grids [m-1]
    constexpr static Real thedz2 = 0.05;
    // bet2 term (ggr/thetal)
    constexpr static Real bet2 = 0.0327;
    // return to isotropy timescale [s]
    constexpr static Real iso = 1000;
    // liquid water flux [K m/s]
    constexpr static Real wthl_sec_zero = 0.01;
    // thetal variance [K^2]
    constexpr static Real thl_sec_zero = 2;
    // vertical velocity variance [m2/s2]
    constexpr static Real w_sec_zero = 0.4;
    // TKE [m2/s2]
    constexpr static Real tke_zero = 0.5;

    // Initialize data structure for bridging to F90
    F0ToF5DiagThirdShocMomentData SDS;

    // Fill in data
    SDS.thedz = thedz;
    SDS.thedz2 = thedz2;
    SDS.bet2 = bet2;
    SDS.iso = iso;
    SDS.isosqrd = iso*iso;
    // for the following moments, feed each level
    //  the same value for this test
    SDS.wthl_sec = wthl_sec_zero;
    SDS.wthl_sec_kc = wthl_sec_zero;
    SDS.wthl_sec_kb = wthl_sec_zero;
    SDS.thl_sec_kc = thl_sec_zero;
    SDS.thl_sec_kb = thl_sec_zero;
    SDS.w_sec = w_sec_zero;
    SDS.w_sec_kc = w_sec_zero;
    SDS.w_sec_zi = w_sec_zero;
    SDS.tke = tke_zero;
    SDS.tke_kc = tke_zero;

    // Be sure inputs are as we expect
    REQUIRE(SDS.thedz > 0);
    REQUIRE(SDS.thedz2 > 0);
    REQUIRE(SDS.wthl_sec_kc == SDS.wthl_sec_kb);
    REQUIRE(SDS.thl_sec_kc == SDS.thl_sec_kb);
    REQUIRE(SDS.w_sec_kc == SDS.w_sec);
    REQUIRE(SDS.tke_kc == SDS.tke);

    // Call the fortran implementation
    f0_to_f5_diag_third_shoc_moment(SDS);

    // Check result, make sure all outputs are zero
    REQUIRE(SDS.f0 == 0);
    REQUIRE(SDS.f1 == 0);
    REQUIRE(SDS.f2 == 0);
    REQUIRE(SDS.f3 == 0);
    REQUIRE(SDS.f4 == 0);
    REQUIRE(SDS.f5 == 0);

    // TEST TWO
    // Positive gradient test.  Feed the function values of the second
    //  moments with positive gradients.  All fterms should have positive values

    // liquid water flux [K m/s]
    constexpr static Real wthl_sec = 0.01;
    // liquid water flux [K m/s] above
    constexpr static Real wthl_sec_kc = 0.02;
    // liquid water flux [K m/s] below
    constexpr static Real wthl_sec_kb = 0;
    // thetal variance [K^2] above
    constexpr static Real thl_sec_kc = 2.5;
    // thetal variance [K^2]
    constexpr static Real thl_sec_kb = 1.7;
    // vertical velocity variance [m2/s2]
    constexpr static Real w_sec = 0.4;
    // vertical velocity variance [m2/s2] above
    constexpr static Real w_sec_kc = 0.5;
    // TKE [m2/s2]
    constexpr static Real tke = 0.5;
    // TKE [m2/s2] above
    constexpr static Real tke_kc = 0.55;

    // Feed in data
    SDS.wthl_sec = wthl_sec;
    SDS.wthl_sec_kc = wthl_sec_kc;
    SDS.wthl_sec_kb = wthl_sec_kb;
    SDS.thl_sec_kc = thl_sec_kc;
    SDS.thl_sec_kb = thl_sec_kb;
    SDS.w_sec = w_sec;
    SDS.w_sec_kc = w_sec_kc;
    SDS.w_sec_zi = w_sec;
    SDS.tke = tke;
    SDS.tke_kc = tke_kc;

    // Verify input is what we want for this test
    REQUIRE(wthl_sec > 0);
    REQUIRE(wthl_sec_kc > wthl_sec_kb);
    REQUIRE(thl_sec_kc > thl_sec_kb);
    REQUIRE(w_sec_kc > w_sec);
    REQUIRE(tke_kc > tke);

    // Call the fortran implementation
    f0_to_f5_diag_third_shoc_moment(SDS);

    // Check result, make sure all outputs are greater than zero
    REQUIRE(SDS.f0 > 0);
    REQUIRE(SDS.f1 > 0);
    REQUIRE(SDS.f2 > 0);
    REQUIRE(SDS.f3 > 0);
    REQUIRE(SDS.f4 > 0);
    REQUIRE(SDS.f5 > 0);

  }

  static void run_bfb()
  {
    // TODO
  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace{

TEST_CASE("shoc_fterm_diag_third_moms_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestFtermdiagThirdMoms;

  TestStruct::run_property();
}

TEST_CASE("shoc_fterm_diag_third_moms_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestFtermdiagThirdMoms;

  TestStruct::run_bfb();
}

} // namespace
