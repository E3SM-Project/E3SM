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
struct UnitWrap::UnitTest<D>::TestGwDiffTend : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up inputs
    GwDiffTendData baseline_data[] = {
      //          ncol, kbot, ktop,  dt
      GwDiffTendData(5,   59,   20, 0.1, init_data[0]),
      GwDiffTendData(6,   58,   19, 0.2, init_data[1]),
      GwDiffTendData(7,   57,   18, 0.3, init_data[2]),
      GwDiffTendData(8,   56,   17, 0.4, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwDiffTendData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwDiffTendData test_data[] = {
      GwDiffTendData(baseline_data[0]),
      GwDiffTendData(baseline_data[1]),
      GwDiffTendData(baseline_data[2]),
      GwDiffTendData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gw_diff_tend(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwDiffTendData& d_baseline = baseline_data[i];
        GwDiffTendData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.dq); ++k) {
          REQUIRE(d_baseline.total(d_baseline.dq) == d_test.total(d_test.dq));
          REQUIRE(d_baseline.dq[k] == d_test.dq[k]);
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

TEST_CASE("gw_diff_tend_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwDiffTend;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
