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
struct UnitWrap::UnitTest<D>::TestCompBruntShocLength : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the SHOC function:
    //   compute_brunt_shoc_length

    // Test for the Brunt Vaissalla frequency.
    // Note that input temperature profiles are selected
    //  deliberately so that it includes a well mixed layer,
    //  an unstable layer, and a conditionally unstable layer
    //  to test a range of conditions.

    // Grid difference centered on thermo grid [m]
    static constexpr Real dz_zt[nlev] = {100, 75, 50, 25, 10};
    // Virtual potential temperature on interface grid [K]
    static constexpr Real thv_zi[nlevi] = {310, 305, 300, 300, 295, 305};

    // Define reasonable bound for output
    static constexpr Real brunt_bound = 1;

    // Initialize data structure for bridging to F90
    ComputeBruntShocLengthData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(nlevi - nlev == 1);
    REQUIRE(shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.dz_zt[offset] = dz_zt[n];
        // For theta_v on thermo grid just take the vertical average
        //  of thv_zi for this simple test.  Just used as a reference
        //  in this subroutine.
        SDS.thv[offset] = 0.5*(thv_zi[n]+thv_zi[n+1]);
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset   = n + s * nlevi;
        SDS.thv_zi[offset] = thv_zi[n];
      }
    }

    // Check that the inputs make sense

    // Be sure that relevant variables are greater than zero
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev - 1; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.dz_zt[offset] > 0);
        REQUIRE(SDS.thv[offset] > 0);
      }

      for(Int n = 0; n < nlevi - 1; ++n) {
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.thv_zi[offset] > 0);
      }
    }

    // Call the C++ implementation.
    compute_brunt_shoc_length(SDS);

    // Check the results
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // Validate that brunt vaisalla frequency
        //  is the correct sign given atmospheric conditions

        // well mixed layer
        if (thv_zi[n] - thv_zi[n+1] == 0){
          REQUIRE(SDS.brunt[offset] == 0);
        }
        // unstable layer
        if (thv_zi[n] - thv_zi[n+1] < 0){
          REQUIRE(SDS.brunt[offset] < 0);
        }
        // stable layer
        if (thv_zi[n] - thv_zi[n+1] > 0){
          REQUIRE(SDS.brunt[offset] > 0);
        }

        // Validate that values fall within some
        //  reasonable bounds for this variable.
        REQUIRE(std::abs(SDS.brunt[offset]) < brunt_bound);
      }
    }
  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    ComputeBruntShocLengthData SDS_baseline[] = {
      //               shcol, nlev, nlevi
      ComputeBruntShocLengthData(10, 71, 72),
      ComputeBruntShocLengthData(10, 12, 13),
      ComputeBruntShocLengthData(7,  16, 17),
      ComputeBruntShocLengthData(2, 7, 8),
    };

    // Generate random input data
    for (auto& d : SDS_baseline) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    ComputeBruntShocLengthData SDS_cxx[] = {
      ComputeBruntShocLengthData(SDS_baseline[0]),
      ComputeBruntShocLengthData(SDS_baseline[1]),
      ComputeBruntShocLengthData(SDS_baseline[2]),
      ComputeBruntShocLengthData(SDS_baseline[3]),
    };

    static constexpr Int num_runs = sizeof(SDS_baseline) / sizeof(ComputeBruntShocLengthData);

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : SDS_baseline) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      compute_brunt_shoc_length(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ComputeBruntShocLengthData& d_baseline = SDS_baseline[i];
        ComputeBruntShocLengthData& d_cxx = SDS_cxx[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.brunt); ++k) {
          REQUIRE(d_baseline.brunt[k] == d_cxx.brunt[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        SDS_cxx[i].write(Base::m_ofile);
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace{

TEST_CASE("shoc_brunt_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompBruntShocLength;

  TestStruct().run_property();
}

TEST_CASE("shoc_brunt_length_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCompBruntShocLength;

  TestStruct().run_bfb();
}

} // namespace
