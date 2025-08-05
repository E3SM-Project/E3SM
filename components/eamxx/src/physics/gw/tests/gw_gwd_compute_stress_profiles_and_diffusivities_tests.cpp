#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "physics/gw/gw_functions.hpp"
#include "physics/gw/tests/infra/gw_test_data.hpp"

#include "gw_unit_tests_common.hpp"

#include <ekat_pack.hpp>

namespace scream {
namespace gw {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestGwdComputeStressProfilesAndDiffusivities : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up inputs
    GwdComputeStressProfilesAndDiffusivitiesData baseline_data[] = {
      //                                        ncol
      GwdComputeStressProfilesAndDiffusivitiesData(2, init_data[0]),
      GwdComputeStressProfilesAndDiffusivitiesData(3, init_data[1]),
      GwdComputeStressProfilesAndDiffusivitiesData(4, init_data[2]),
      GwdComputeStressProfilesAndDiffusivitiesData(5, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwdComputeStressProfilesAndDiffusivitiesData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      // ni must be very small or else we risk a FPE due to a huge exp
      d.randomize(engine, { {d.ni, {1.E-08, 2.E-08}} });
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwdComputeStressProfilesAndDiffusivitiesData test_data[] = {
      GwdComputeStressProfilesAndDiffusivitiesData(baseline_data[0]),
      GwdComputeStressProfilesAndDiffusivitiesData(baseline_data[1]),
      GwdComputeStressProfilesAndDiffusivitiesData(baseline_data[2]),
      GwdComputeStressProfilesAndDiffusivitiesData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gwd_compute_stress_profiles_and_diffusivities(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwdComputeStressProfilesAndDiffusivitiesData& d_baseline = baseline_data[i];
        GwdComputeStressProfilesAndDiffusivitiesData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.tau); ++k) {
          REQUIRE(d_baseline.total(d_baseline.tau) == d_test.total(d_test.tau));
          REQUIRE(d_baseline.tau[k] == d_test.tau[k]);
        }

      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        test_data[i].write(Base::m_ofile);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace gw
} // namespace scream

namespace {

TEST_CASE("gwd_compute_stress_profiles_and_diffusivities_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwdComputeStressProfilesAndDiffusivities;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
