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
struct UnitWrap::UnitTest<D>::TestGwdComputeTendenciesFromStressDivergence : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    GwdComputeTendenciesFromStressDivergenceData baseline_data[] = {
      // ncol, pver, pgwv, ngwv, do_taper, dt, effgw
      GwdComputeTendenciesFromStressDivergenceData(72),
      GwdComputeTendenciesFromStressDivergenceData(63),
      GwdComputeTendenciesFromStressDivergenceData(47),
      GwdComputeTendenciesFromStressDivergenceData(31),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwdComputeTendenciesFromStressDivergenceData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before fortran calls so that
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
        d.read(Base::m_fid);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gwd_compute_tendencies_from_stress_divergence(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwdComputeTendenciesFromStressDivergenceData& d_baseline = baseline_data[i];
        GwdComputeTendenciesFromStressDivergenceData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.tau); ++k) {
          REQUIRE(d_baseline.total(d_baseline.tau) == d_test.total(d_test.tau));
          REQUIRE(d_baseline.tau[k] == d_test.tau[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.gwut); ++k) {
          REQUIRE(d_baseline.total(d_baseline.gwut) == d_test.total(d_test.gwut));
          REQUIRE(d_baseline.gwut[k] == d_test.gwut[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.utgw); ++k) {
          REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.utgw));
          REQUIRE(d_baseline.utgw[k] == d_test.utgw[k]);
          REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.vtgw));
          REQUIRE(d_baseline.vtgw[k] == d_test.vtgw[k]);
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

TEST_CASE("gwd_compute_tendencies_from_stress_divergence_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwdComputeTendenciesFromStressDivergence;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
