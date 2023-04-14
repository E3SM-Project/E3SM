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
struct UnitWrap::UnitTest<D>::TestCheckShocLength {

  static void run_property()
  {
    static constexpr Real maxlen = scream::shoc::Constants<Real>::maxlen;
    static constexpr Real minlen = scream::shoc::Constants<Real>::minlen;
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC function:
    //   check_length_scale_shoc_length

    // Test this simple function.  Feed the routine values of
    //  mixing length that are zero and much larger than the
    //  host grid box size to make sure that the function
    //  corrects these erroneous values.

    // Define the host grid box size x-direction [m]
    static constexpr Real host_dx = 3000;
    // Defin the host grid box size y-direction [m]
    static constexpr Real host_dy = 5000;
    // Define mixing length [m]
    //   In profile include values of zero and very large values
    Real shoc_mix[nlev] = {50000, 4000, 2000, 0, 500};

    // Initialize data structure for bridging to F90
    CheckLengthScaleShocLengthData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    REQUIRE(shcol > 0);

    // compute geometric grid mesh
    const auto grid_mesh = sqrt(host_dx*host_dy);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      SDS.host_dx[s] = host_dx;
      SDS.host_dy[s] = host_dy;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.shoc_mix[offset] = shoc_mix[n];
      }
    }

    // Check that the inputs make sense

    // Be sure that relevant variables are greater than zero
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.host_dx[s] > 0);
      REQUIRE(SDS.host_dy[s] > 0);
    }

    // Call the fortran implementation
    check_length_scale_shoc_length(SDS);

    // Check the results
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Require mixing length is greater than zero and is
        //  less than geometric grid mesh length + 1 m
        REQUIRE(SDS.shoc_mix[offset] >= minlen);
        REQUIRE(SDS.shoc_mix[offset] <= maxlen);
        REQUIRE(SDS.shoc_mix[offset] < 1.0+grid_mesh);
      }
    }
  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    CheckLengthScaleShocLengthData SDS_f90[] = {
      //             shcol, nlev
      CheckLengthScaleShocLengthData(10, 71),
      CheckLengthScaleShocLengthData(10, 12),
      CheckLengthScaleShocLengthData(7,  16),
      CheckLengthScaleShocLengthData(2, 7),
    };

    // Generate random input data
    for (auto& d : SDS_f90) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    CheckLengthScaleShocLengthData SDS_cxx[] = {
      CheckLengthScaleShocLengthData(SDS_f90[0]),
      CheckLengthScaleShocLengthData(SDS_f90[1]),
      CheckLengthScaleShocLengthData(SDS_f90[2]),
      CheckLengthScaleShocLengthData(SDS_f90[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      check_length_scale_shoc_length(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      d.transpose<ekat::TransposeDirection::c2f>();
      // expects data in fortran layout
      check_length_scale_shoc_length_f(d.nlev,d.shcol,d.host_dx,d.host_dy,d.shoc_mix);
      d.transpose<ekat::TransposeDirection::f2c>();
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(SDS_f90) / sizeof(CheckLengthScaleShocLengthData);
      for (Int i = 0; i < num_runs; ++i) {
        CheckLengthScaleShocLengthData& d_f90 = SDS_f90[i];
        CheckLengthScaleShocLengthData& d_cxx = SDS_cxx[i];
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

namespace{

TEST_CASE("shoc_check_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCheckShocLength;

  TestStruct::run_property();
}

TEST_CASE("shoc_check_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCheckShocLength;

  TestStruct::run_bfb();
}

} // namespace
