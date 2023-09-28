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
struct UnitWrap::UnitTest<D>::TestClipThirdMoms {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlevi    = 5;

    // Tests for the SHOC function:
    //   clipping_diag_third_shoc_moments

    //  Test to be sure that very high values of w3
    //    are reduced but still the same sign

    // Define the second moment of vertical velocity
    static constexpr Real w_sec_zi[nlevi] = {0.3, 0.4, 0.5, 0.4, 0.1};
    // Define the third moment of vertical velocity
    static constexpr Real w3_in[nlevi] = {0.2, 99999, -99999, 0.2, 0.05};

    // Define a local logical
    bool w3_large;

    // Initialize data structure for bridging to F90
    ClippingDiagThirdShocMomentsData SDS(shcol, nlevi);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol > 0);
    REQUIRE(SDS.nlevi > 1);

    // load up the input data
    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.w_sec_zi[offset] = w_sec_zi[n];
        SDS.w3[offset] = w3_in[n];
      }
    }

    // Verify the data
    // For this particular test wen want to be sure that
    //   the input has relatively small values of w2 (let's say
    //   all smaller than 1, which is reasonable) and let's
    //   be sure that w3 has some very large unreasonable values
    for(Int s = 0; s < shcol; ++s) {
      // Initialize conditional to make sure there are
      //   large values of w3 in input
      w3_large = false;
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        REQUIRE(SDS.w_sec_zi[offset] <= 1);
        if (abs(SDS.w3[offset]) > 1000){
          w3_large = true;
        }
        SDS.w_sec_zi[offset] = w_sec_zi[n];
        SDS.w3[offset] = w3_in[n];
      }
      REQUIRE(w3_large == true);
    }

    // Call the fortran implementation
    clipping_diag_third_shoc_moments(SDS);

    // Check the result
    // For large values of w3, verify that the result has been reduced
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        if (abs(w3_in[n]) > 1000){
          REQUIRE(abs(SDS.w3[offset]) < abs(w3_in[n]));
        }

      }
    }

  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ClippingDiagThirdShocMomentsData SDS_f90[] = {
      //               shcol, nlevi
      ClippingDiagThirdShocMomentsData(10, 72),
      ClippingDiagThirdShocMomentsData(10, 13),
      ClippingDiagThirdShocMomentsData(7, 17),
      ClippingDiagThirdShocMomentsData(2, 8),
    };

    // Generate random input data
    for (auto& d : SDS_f90) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ClippingDiagThirdShocMomentsData SDS_cxx[] = {
      ClippingDiagThirdShocMomentsData(SDS_f90[0]),
      ClippingDiagThirdShocMomentsData(SDS_f90[1]),
      ClippingDiagThirdShocMomentsData(SDS_f90[2]),
      ClippingDiagThirdShocMomentsData(SDS_f90[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      clipping_diag_third_shoc_moments(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      d.transpose<ekat::TransposeDirection::c2f>();
      // expects data in fortran layout
      clipping_diag_third_shoc_moments_f(d.nlevi,d.shcol,d.w_sec_zi,d.w3);
      d.transpose<ekat::TransposeDirection::f2c>();
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(SDS_f90) / sizeof(ClippingDiagThirdShocMomentsData);
      for (Int i = 0; i < num_runs; ++i) {
        ClippingDiagThirdShocMomentsData& d_f90 = SDS_f90[i];
        ClippingDiagThirdShocMomentsData& d_cxx = SDS_cxx[i];
        for (Int k = 0; k < d_f90.total(d_f90.w3); ++k) {
          REQUIRE(d_f90.w3[k] == d_cxx.w3[k]);
        }
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace{

TEST_CASE("shoc_clip_third_moms_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestClipThirdMoms;

  TestStruct::run_property();
}

TEST_CASE("shoc_clip_third_moms_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestClipThirdMoms;

  TestStruct::run_bfb();
}

} // namespace
