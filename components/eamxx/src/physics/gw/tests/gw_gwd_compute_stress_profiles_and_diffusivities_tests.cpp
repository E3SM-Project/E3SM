#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/gw/gw_functions.hpp"
#include "physics/gw/tests/infra/gw_test_data.hpp"

#include "gw_unit_tests_common.hpp"

namespace scream {
namespace gw {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestGwdComputeStressProfilesAndDiffusivities : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    GwInit init_data[] = {
          // pver, pgwv,   dc, orog_only, molec_diff, tau_0_ubc, nbot_molec, ktop, kbotbg, fcrit2, kwv
      GwInit(  72,   20, 0.75,     false,      false,     false,         16,   60,     16,    .67, 6.28e-5),
      GwInit(  72,   20, 0.75,     true ,      false,     true ,         16,   60,     16,    .67, 6.28e-5),
      GwInit(  72,   20, 0.75,     false,      true ,     true ,         16,   60,     16,    .67, 6.28e-5),
      GwInit(  72,   20, 0.75,     true ,      true ,     false,         16,   60,     16,    .67, 6.28e-5),
    };

    for (auto& d : init_data) {
      d.randomize(engine);
    }

    // Set up inputs
    GwdComputeStressProfilesAndDiffusivitiesData baseline_data[] = {
      GwdComputeStressProfilesAndDiffusivitiesData(2, 10, init_data[0]),
      GwdComputeStressProfilesAndDiffusivitiesData(3, 11, init_data[1]),
      GwdComputeStressProfilesAndDiffusivitiesData(4, 12, init_data[2]),
      GwdComputeStressProfilesAndDiffusivitiesData(5, 13, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwdComputeStressProfilesAndDiffusivitiesData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
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
        d.read(Base::m_fid);
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
        test_data[i].write(Base::m_fid);
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
