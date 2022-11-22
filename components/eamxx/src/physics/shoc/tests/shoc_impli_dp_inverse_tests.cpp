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
struct UnitWrap::UnitTest<D>::TestImpDpInverse {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC subroutine
    //   dp_inverse

    // TEST
    // Load up two columns, one with smaller dz.  The case with
    //  smaller dz should return a HIGHER value of rdp_zt

    // Define height thickness on nlev grid [m]
    static constexpr Real dz_zt[nlev] = {500, 200, 100, 50, 10};
    // Define density on zt grid [kg/m3]
    static constexpr Real rho_zt[nlev] = {0.6, 0.8, 0.9, 1.0, 1.2};

    // Initialize data structure for bridging to F90
    DpInverseData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    REQUIRE(shcol == 2);

    // Fill in test data on zi_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // Feed second column SMALLER dz values
        SDS.dz_zt[offset] = dz_zt[n]/(1+s);
        SDS.rho_zt[offset] = rho_zt[n];
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        const auto offsets = n + (1+s)*nlev;
        // make sure that density is in reasonable bounds
        REQUIRE( (SDS.rho_zt[offset] > 0 && SDS.rho_zt[offset] < 1.5) );
        REQUIRE(SDS.dz_zt[offset] > 0);
        // Verify that the second column has smaller dz values
        if (s < shcol-1){
          REQUIRE(SDS.dz_zt[offset] > SDS.dz_zt[offsets]);
        }
      }
    }

    // Call the fortran implementation
    dp_inverse(SDS);

    // Verify result
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        const auto offsets = n + (1+s)*nlev;

        // Verify bounds are reasonable
        REQUIRE(SDS.rdp_zt[offset] > 0);
        // By dimensional analysis the value of rdp_zt
        //  with reasonable inputs should be less than 1
        REQUIRE(SDS.rdp_zt[offset] < 1);
        // Verify that rdp_zt values are larger when dz is smaller
        if (s < shcol-1){
          REQUIRE(SDS.rdp_zt[offset] < SDS.rdp_zt[offsets]);
        }
      }
    }

  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    DpInverseData f90_data[] = {
      //            shcol, nlev
      DpInverseData(10, 71),
      DpInverseData(10, 12),
      DpInverseData(7,  16),
      DpInverseData(2,   7)
    };

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    DpInverseData cxx_data[] = {
      DpInverseData(f90_data[0]),
      DpInverseData(f90_data[1]),
      DpInverseData(f90_data[2]),
      DpInverseData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      dp_inverse(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      dp_inverse_f(d.nlev, d.shcol, d.rho_zt, d.dz_zt, d.rdp_zt);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(DpInverseData);
      for (Int i = 0; i < num_runs; ++i) {
        DpInverseData& d_f90 = f90_data[i];
        DpInverseData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.rdp_zt); ++k) {
          REQUIRE(d_f90.rdp_zt[k] == d_cxx.rdp_zt[k]);
        }
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_imp_dp_inverse_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestImpDpInverse;

  TestStruct::run_property();
}

TEST_CASE("shoc_imp_dp_inverse_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestImpDpInverse;

  TestStruct::run_bfb();
}

} // namespace
