#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

#include <ekat_pack.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestFindMseMax : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    FindMseMaxData baseline_data[] = {
      //             pcols, ncol, pver, num_msg, pergro_active
      FindMseMaxData(    4,    4,   72 ,      30, false),
      FindMseMaxData(    4,    4,   72,       30, true),
      FindMseMaxData(    4,    4,   128,      30, false),
      FindMseMaxData(    4,    4,   128,      30, true),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(FindMseMaxData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    FindMseMaxData test_data[] = {
      FindMseMaxData(baseline_data[0]),
      FindMseMaxData(baseline_data[1]),
      FindMseMaxData(baseline_data[2]),
      FindMseMaxData(baseline_data[3]),
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
        find_mse_max_f(d);
      }
      else {
        find_mse_max(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        FindMseMaxData& d_baseline = baseline_data[i];
        FindMseMaxData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.mse_max_val) == d_test.total(d_test.mse_max_val));
        REQUIRE(d_baseline.total(d_baseline.mse_max_val) == d_test.total(d_test.msemax_klev));
        for (Int k = 0; k < d_baseline.total(d_baseline.mse_max_val); ++k) {
          REQUIRE(d_baseline.mse_max_val[k] == d_test.mse_max_val[k]);
          REQUIRE(d_baseline.msemax_klev[k] == d_test.msemax_klev[k]);
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
} // namespace zm
} // namespace scream

namespace {

TEST_CASE("find_mse_max_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestFindMseMax;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
