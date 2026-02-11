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
struct UnitWrap::UnitTest<D>::TestZmTransportMomentum : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ZmTransportMomentumData baseline_data[] = {
      //                      pcols, ncol, pver, pverp, nwind, il1g, il2g, dt
      ZmTransportMomentumData(    4,    4,   72,    73,     4,    0,    3, 0.5),
      ZmTransportMomentumData(    4,    4,   72,    73,     4,    0,    3, 0.5),
      ZmTransportMomentumData(    4,    4,   72,    73,     4,    0,    3, 0.5),
      ZmTransportMomentumData(    4,    4,  128,   129,     4,    0,    3, 0.5),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmTransportMomentumData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmTransportMomentumData test_data[] = {
      ZmTransportMomentumData(baseline_data[0]),
      ZmTransportMomentumData(baseline_data[1]),
      ZmTransportMomentumData(baseline_data[2]),
      ZmTransportMomentumData(baseline_data[3]),
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
        zm_transport_momentum_f(d);
      }
      else {
        zm_transport_momentum(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmTransportMomentumData& d_baseline = baseline_data[i];
        ZmTransportMomentumData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.wind_tend) == d_test.total(d_test.wind_tend));
        REQUIRE(d_baseline.total(d_baseline.wind_tend) == d_test.total(d_test.pguall));
        REQUIRE(d_baseline.total(d_baseline.wind_tend) == d_test.total(d_test.pgdall));
        REQUIRE(d_baseline.total(d_baseline.wind_tend) == d_test.total(d_test.icwu));
        REQUIRE(d_baseline.total(d_baseline.wind_tend) == d_test.total(d_test.icwd));
        for (Int k = 0; k < d_baseline.total(d_baseline.wind_tend); ++k) {
          REQUIRE(d_baseline.wind_tend[k] == d_test.wind_tend[k]);
          REQUIRE(d_baseline.pguall[k] == d_test.pguall[k]);
          REQUIRE(d_baseline.pgdall[k] == d_test.pgdall[k]);
          REQUIRE(d_baseline.icwu[k] == d_test.icwu[k]);
          REQUIRE(d_baseline.icwd[k] == d_test.icwd[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.seten) == d_test.total(d_test.seten));
        for (Int k = 0; k < d_baseline.total(d_baseline.seten); ++k) {
          REQUIRE(d_baseline.seten[k] == d_test.seten[k]);
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

TEST_CASE("zm_transport_momentum_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmTransportMomentum;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
