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
struct UnitWrap::UnitTest<D>::TestShocCheckTke : public UnitWrap::UnitTest<D>::Base {

  void run_property()
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

    // Call the C++ implementation.
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

  void run_bfb()
  {
    auto engine = Base::get_engine();

    CheckTkeData SDS_baseline[] = {
      //               shcol, nlev
      CheckTkeData(10, 71),
      CheckTkeData(10, 12),
      CheckTkeData(7,  16),
      CheckTkeData(2,   7),
    };

    // Generate random input data
    for (auto& d : SDS_baseline) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    CheckTkeData SDS_cxx[] = {
      CheckTkeData(SDS_baseline[0]),
      CheckTkeData(SDS_baseline[1]),
      CheckTkeData(SDS_baseline[2]),
      CheckTkeData(SDS_baseline[3]),
    };

    static constexpr Int num_runs = sizeof(SDS_baseline) / sizeof(CheckTkeData);

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : SDS_baseline) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      check_tke(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        CheckTkeData& d_baseline = SDS_baseline[i];
        CheckTkeData& d_cxx = SDS_cxx[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.tke); ++k) {
          REQUIRE(d_baseline.tke[k]    == d_cxx.tke[k]);
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

namespace {

TEST_CASE("shoc_check_tke_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocCheckTke;

  TestStruct().run_property();
}

TEST_CASE("shoc_check_tke_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocCheckTke;

  TestStruct().run_bfb();

}

} // namespace
