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
struct UnitWrap::UnitTest<D>::TestShocEnergyDseFixer {

  static void run_property()
  {
    static constexpr Int shcol    = 6;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC function
    //     shoc_energy_dse_fixer

    // TEST
    // For columns that are identical EXCEPT for the shoctop indicee,
    //  verify that given a positive value of energy imbalance that
    //  columns with a higher SHOC top were subject to more energy removal.

    // Host model dry static energy [J kg-1]
    static constexpr Real host_dse_input[nlev] = {350e3, 325e3, 315e3, 310e3, 300e3};

    // Energy disbalance.  For this test we assume all columns have
    //  the same disbalance magnitude.
    static constexpr Real se_dis = 0.1;

    // level indicee of SHOC top layer
    static constexpr Int shoctop[shcol] = {5, 3, 1, 2, 4, 4};

    // Initialize data structure for bridging to F90
    ShocEnergyDseFixerData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    // for this test we need exactly six columns
    REQUIRE( (SDS.shcol == 6 && SDS.nlev == nlev) );

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      SDS.shoctop[s] = shoctop[s];
      SDS.se_dis[s] = se_dis;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.host_dse[offset] = host_dse_input[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      // For this test we WANT se_dis > 0
      REQUIRE(SDS.se_dis[s] > 0.0);
      REQUIRE(SDS.shoctop[s] >= 1);
      REQUIRE(SDS.shoctop[s] <= nlev);
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;

        REQUIRE(SDS.host_dse[offset] > 0.0);
      }
    }

    // Call the fortran implementation
    shoc_energy_dse_fixer(SDS);

    // Check the results
    Real temp_sum[shcol];
    for(Int s = 0; s < shcol; ++s) {
      temp_sum[s] = 0.0;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        temp_sum[s] += SDS.host_dse[offset];
      }
    }

    // Verify that as shoctop values get lower that the
    // summation of temperatures also gets lower.  This is proportionally
    // to the amount of energy we expect to be removed from a column.
    for (Int s = 0; s < shcol-1; ++s) {
      if (shoctop[s] < shoctop[s+1]){
        REQUIRE(temp_sum[s] < temp_sum[s+1]);
      }
      else if (shoctop[s] > shoctop[s+1]){
        REQUIRE(temp_sum[s] > temp_sum[s+1]);
      }
      else{
        REQUIRE(temp_sum[s] == temp_sum[s+1]);
      }
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

TEST_CASE("shoc_energy_dse_fixer_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyDseFixer;

  TestStruct::run_property();
}

TEST_CASE("shoc_energy_dse_fixer_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyDseFixer;

  TestStruct::run_bfb();
}

} // namespace
