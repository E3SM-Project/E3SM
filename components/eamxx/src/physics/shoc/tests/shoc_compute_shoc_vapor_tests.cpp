#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestComputeShocVapor : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC subroutine:
    //  compute_shoc_vapor

    // Test
    // Given profile of a variety of conditions,
    //  verify that the output is as expected

    // Total water mixing ratio [kg/kg]
    static constexpr Real qw[nlev] = {1e-2, 1.2e-2, 1.5e-2, 1.7e-2, 2.0e-2};
    // Liquid water mixing ratio [kg/kg]
    static constexpr Real ql[nlev] = {0, 0, 1.5e-4, 2e-3, 0};

    // Initialize data structure for bridging to F90
    ComputeShocVaporData SDS(shcol, nlev);

    // Load input data
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.qw[offset] = qw[n];
        SDS.ql[offset] = ql[n];

      }
    }

    // Check that inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // total water should be greater than zero
        REQUIRE(SDS.qw[offset] > 0);
        // cloud water greater than or equal to zero
        REQUIRE(SDS.ql[offset] >= 0);
        // total water should be greater than cloud water
        REQUIRE(SDS.qw[offset] > SDS.ql[offset]);

      }
    }

    // Call the C++ implementation
    compute_shoc_vapor(SDS);

    // Verify the result
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // vapor should be greater than zero
        REQUIRE(SDS.qv[offset] > 0);
        // if cloud present vapor should be less than total water
        if (SDS.ql[offset] > 0){
          REQUIRE(SDS.qv[offset] < SDS.qw[offset]);
        }
        // else they should be equal
        else{
          REQUIRE(SDS.qv[offset] == SDS.qw[offset]);
        }

      }
    }

  } // run_property

  void run_bfb()
  {
    auto engine = Base::get_engine();

    ComputeShocVaporData baseline_data[] = {
      //              shcol, nlev
      ComputeShocVaporData(10, 71),
      ComputeShocVaporData(10, 12),
      ComputeShocVaporData(7,  16),
      ComputeShocVaporData(2,   7),
    };

    // Generate random input data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    ComputeShocVaporData cxx_data[] = {
      ComputeShocVaporData(baseline_data[0]),
      ComputeShocVaporData(baseline_data[1]),
      ComputeShocVaporData(baseline_data[2]),
      ComputeShocVaporData(baseline_data[3]),
    };

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      compute_shoc_vapor(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ComputeShocVaporData);
      for (Int i = 0; i < num_runs; ++i) {
        ComputeShocVaporData& d_baseline = baseline_data[i];
        ComputeShocVaporData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.qv); ++k) {
          REQUIRE(d_baseline.qv[k] == d_cxx.qv[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (auto& d : cxx_data) {
        d.write(Base::m_fid);
      }
    }
  } // run_bfb
};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("compute_shoc_vapor_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeShocVapor;

  TestStruct().run_property();
}

TEST_CASE("compute_shoc_vapor_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeShocVapor;

  TestStruct().run_bfb();
}

} // empty namespace
