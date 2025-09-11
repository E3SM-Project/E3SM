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
struct UnitWrap::UnitTest<D>::TestGwdComputeTendenciesFromStressDivergence : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up inputs
    GwdComputeTendenciesFromStressDivergenceData baseline_data[] = {
      //                                        ncol, do_taper,   dt, effgw, init
      GwdComputeTendenciesFromStressDivergenceData(2,    false,  0.4,   0.3, init_data[0]),
      GwdComputeTendenciesFromStressDivergenceData(3,    false,  0.4,   0.3, init_data[1]),
      GwdComputeTendenciesFromStressDivergenceData(4,    true ,  0.4,   0.3, init_data[2]),
      GwdComputeTendenciesFromStressDivergenceData(5,    true ,  0.4,   0.3, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwdComputeTendenciesFromStressDivergenceData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (Int i = 0; i < num_runs; ++i) {
      GwdComputeTendenciesFromStressDivergenceData& d = baseline_data[i];
      d.randomize(engine,  { {d.tend_level, {init_data[i].ktop+1, init_data[i].kbotbg-1}} });
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwdComputeTendenciesFromStressDivergenceData test_data[] = {
      GwdComputeTendenciesFromStressDivergenceData(baseline_data[0]),
      GwdComputeTendenciesFromStressDivergenceData(baseline_data[1]),
      GwdComputeTendenciesFromStressDivergenceData(baseline_data[2]),
      GwdComputeTendenciesFromStressDivergenceData(baseline_data[3]),
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
        gwd_compute_tendencies_from_stress_divergence_f(d);
      }
      else {
        gwd_compute_tendencies_from_stress_divergence(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwdComputeTendenciesFromStressDivergenceData& d_baseline = baseline_data[i];
        GwdComputeTendenciesFromStressDivergenceData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.tau) == d_test.total(d_test.tau));
        for (Int k = 0; k < d_baseline.total(d_baseline.tau); ++k) {
          REQUIRE(d_baseline.tau[k] == d_test.tau[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.gwut) == d_test.total(d_test.gwut));
        for (Int k = 0; k < d_baseline.total(d_baseline.gwut); ++k) {
          REQUIRE(d_baseline.gwut[k] == d_test.gwut[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.utgw));
        REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.vtgw));
        for (Int k = 0; k < d_baseline.total(d_baseline.utgw); ++k) {
          REQUIRE(d_baseline.utgw[k] == d_test.utgw[k]);
          REQUIRE(d_baseline.vtgw[k] == d_test.vtgw[k]);
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

TEST_CASE("gwd_compute_tendencies_from_stress_divergence_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwdComputeTendenciesFromStressDivergence;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
