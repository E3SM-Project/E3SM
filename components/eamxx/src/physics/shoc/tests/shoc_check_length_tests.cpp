#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/eamxx_types.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

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
struct UnitWrap::UnitTest<D>::TestCheckShocLength : public UnitWrap::UnitTest<D>::Base {

  void run_property()
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

    // Call the C++ implementation.
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

  void run_bfb()
  {
    auto engine = Base::get_engine();

    CheckLengthScaleShocLengthData SDS_baseline[] = {
      //             shcol, nlev
      CheckLengthScaleShocLengthData(10, 71),
      CheckLengthScaleShocLengthData(10, 12),
      CheckLengthScaleShocLengthData(7,  16),
      CheckLengthScaleShocLengthData(2, 7),
    };

    // Generate random input data
    for (auto& d : SDS_baseline) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    CheckLengthScaleShocLengthData SDS_cxx[] = {
      CheckLengthScaleShocLengthData(SDS_baseline[0]),
      CheckLengthScaleShocLengthData(SDS_baseline[1]),
      CheckLengthScaleShocLengthData(SDS_baseline[2]),
      CheckLengthScaleShocLengthData(SDS_baseline[3]),
    };

    static constexpr Int num_runs = sizeof(SDS_baseline) / sizeof(CheckLengthScaleShocLengthData);

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : SDS_baseline) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      check_length_scale_shoc_length(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        CheckLengthScaleShocLengthData& d_baseline = SDS_baseline[i];
        CheckLengthScaleShocLengthData& d_cxx = SDS_cxx[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.shoc_mix); ++k) {
          REQUIRE(d_baseline.shoc_mix[k] == d_cxx.shoc_mix[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        SDS_cxx[i].write(Base::m_fid);
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

  TestStruct().run_property();
}

TEST_CASE("shoc_check_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCheckShocLength;

  TestStruct().run_bfb();
}

} // namespace
