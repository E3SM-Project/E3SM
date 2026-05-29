#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZmStateUpdate : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    //                            pcols, ncol, pver, pverp, dt
    ZmStateUpdateData baseline_data[] = {
      ZmStateUpdateData(4,    4,  72,   73,    1800.0),
      ZmStateUpdateData(4,    4,  128,  129,   1800.0),
      ZmStateUpdateData(4,    4,  72,   73,    3600.0),
      ZmStateUpdateData(4,    4,  128,  129,   3600.0),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmStateUpdateData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmStateUpdateData test_data[] = {
      ZmStateUpdateData(baseline_data[0]),
      ZmStateUpdateData(baseline_data[1]),
      ZmStateUpdateData(baseline_data[2]),
      ZmStateUpdateData(baseline_data[3]),
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
        zm_state_update_f(d);
      }
      else {
        zm_state_update(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmStateUpdateData& d_baseline = baseline_data[i];
        ZmStateUpdateData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.zm) == d_test.total(d_test.zm));
        for (Int k = 0; k < d_baseline.total(d_baseline.zm); ++k) {
          REQUIRE(d_baseline.zm[k] == d_test.zm[k]);
          REQUIRE(d_baseline.t[k]  == d_test.t[k]);
          REQUIRE(d_baseline.qv[k] == d_test.qv[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.zi) == d_test.total(d_test.zi));
        for (Int k = 0; k < d_baseline.total(d_baseline.zi); ++k) {
          REQUIRE(d_baseline.zi[k] == d_test.zi[k]);
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

TEST_CASE("zm_state_update_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmStateUpdate;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
