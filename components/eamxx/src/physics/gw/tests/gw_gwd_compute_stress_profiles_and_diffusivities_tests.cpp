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
    for (Int i = 0; i < num_runs; ++i) {
      auto& d = baseline_data[i];
      // ni must be very small or else we risk a FPE due to a huge exp
      d.randomize(engine, { {d.ni, {1.E-06, 2.E-06}}, {d.src_level, {init_data[i].ktop+1, init_data[i].kbotbg-1}}, {d.ubi, {2.E-04, 3.E-04}}, {d.c, {1.E-04, 2.E-04}} });
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
      if (this->m_baseline_action == GENERATE) {
        gwd_compute_stress_profiles_and_diffusivities_f(d);
      }
      else {
        gwd_compute_stress_profiles_and_diffusivities(d);
      }
    }

    const auto margin = std::numeric_limits<Real>::epsilon() *
      (ekat::is_single_precision<Real>::value ? 1000 : 1);

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwdComputeStressProfilesAndDiffusivitiesData& d_baseline = baseline_data[i];
        GwdComputeStressProfilesAndDiffusivitiesData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.tau) == d_test.total(d_test.tau));
        for (Int k = 0; k < d_baseline.total(d_baseline.tau); ++k) {
          // We must add a tolerance here since we are doing the operations
          // in a different order in order to improve the amount of parallel
          // computation. This tol can be removed once we are no longer using
          // fortran to generate baselines.
          REQUIRE(d_baseline.tau[k] == Approx(d_test.tau[k]).margin(margin));
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
