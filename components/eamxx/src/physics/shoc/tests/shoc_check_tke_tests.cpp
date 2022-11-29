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
struct UnitWrap::UnitTest<D>::TestShocCheckTke {

  static void run_property()
  {
    static constexpr Real mintke = scream::shoc::Constants<Real>::mintke;
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;

    // Tests for SHOC subroutine
    //  check_tke

    // Test to make sure that given negative values of TKE
    // that these values are corrected

    static constexpr Real tke_input[nlev] = {-0.3, 0.4, -100, -1.5, 0.4};

    // Intialize data structure for bridging to F90
    CheckTkeData SDS(shcol, nlev);

    // Load the data
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.tke[offset] = tke_input[n];
      }
    }

    // Check some input
    REQUIRE((SDS.shcol > 0 && SDS.nlev > 0));

    // call the fortran implementation
    check_tke(SDS);

    // Check the result against the input values
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // if input TKE was less than zero, verify it was adjusted
        if (tke_input[n] < 0){
          REQUIRE(SDS.tke[offset] >= mintke);
        }
        // Else make sure TKE remains untouched
        else{
          REQUIRE(SDS.tke[offset] == tke_input[n]);
        }
      }
    }

  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    CheckTkeData SDS_f90[] = {
      //               shcol, nlev
      CheckTkeData(10, 71),
      CheckTkeData(10, 12),
      CheckTkeData(7,  16),
      CheckTkeData(2,   7),
    };

    // Generate random input data
    for (auto& d : SDS_f90) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    CheckTkeData SDS_cxx[] = {
      CheckTkeData(SDS_f90[0]),
      CheckTkeData(SDS_f90[1]),
      CheckTkeData(SDS_f90[2]),
      CheckTkeData(SDS_f90[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      check_tke(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      d.transpose<ekat::TransposeDirection::c2f>();
      // expects data in fortran layout
      check_tke_f(d.nlev, d.shcol, d.tke);
      d.transpose<ekat::TransposeDirection::f2c>();
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(SDS_f90) / sizeof(CheckTkeData);
      for (Int i = 0; i < num_runs; ++i) {
        CheckTkeData& d_f90 = SDS_f90[i];
        CheckTkeData& d_cxx = SDS_cxx[i];
        for (Int k = 0; k < d_f90.total(d_f90.tke); ++k) {
          REQUIRE(d_f90.tke[k]    == d_cxx.tke[k]);
        }
      }
    }
  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_check_tke_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocCheckTke;

  TestStruct::run_property();
}

TEST_CASE("shoc_check_tke_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocCheckTke;

  TestStruct::run_bfb();

}

} // namespace
