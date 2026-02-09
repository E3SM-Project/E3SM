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
struct UnitWrap::UnitTest<D>::TestZmTransportTracer : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ZmTransportTracerData baseline_data[] = {
      //                    pcols, ncol, pver, ncnst, il1g, il2g, dt
      ZmTransportTracerData(    4,    4,   72,     6,    1,    2, 0.5),
      ZmTransportTracerData(    4,    4,   72,     6,    2,    3, 1.5),
      ZmTransportTracerData(    4,    4,   72,     6,    4,    2, 2.5),
      ZmTransportTracerData(    4,    4,  128,     6,    4,    2, 2.5),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmTransportTracerData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmTransportTracerData test_data[] = {
      ZmTransportTracerData(baseline_data[0]),
      ZmTransportTracerData(baseline_data[1]),
      ZmTransportTracerData(baseline_data[2]),
      ZmTransportTracerData(baseline_data[3]),
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
        zm_transport_tracer_f(d);
      }
      else {
        zm_transport_tracer(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmTransportTracerData& d_baseline = baseline_data[i];
        ZmTransportTracerData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.dqdt) == d_test.total(d_test.dqdt));
        for (Int k = 0; k < d_baseline.total(d_baseline.dqdt); ++k) {
          REQUIRE(d_baseline.dqdt[k] == d_test.dqdt[k]);
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

TEST_CASE("zm_transport_tracer_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmTransportTracer;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
