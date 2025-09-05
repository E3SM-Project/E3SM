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
struct UnitWrap::UnitTest<D>::TestGwdProjectTau : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up inputs
    GwdProjectTauData baseline_data[] = {
      //             ncol
      GwdProjectTauData(2, init_data[0]),
      GwdProjectTauData(3, init_data[1]),
      GwdProjectTauData(4, init_data[2]),
      GwdProjectTauData(5, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwdProjectTauData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwdProjectTauData test_data[] = {
      GwdProjectTauData(baseline_data[0]),
      GwdProjectTauData(baseline_data[1]),
      GwdProjectTauData(baseline_data[2]),
      GwdProjectTauData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gwd_project_tau(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwdProjectTauData& d_baseline = baseline_data[i];
        GwdProjectTauData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.taucd); ++k) {
          REQUIRE(d_baseline.total(d_baseline.taucd) == d_test.total(d_test.taucd));
          REQUIRE(d_baseline.taucd[k] == d_test.taucd[k]);
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

TEST_CASE("gwd_project_tau_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwdProjectTau;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
