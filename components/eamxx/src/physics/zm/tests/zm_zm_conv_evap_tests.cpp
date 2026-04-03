#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZmConvEvap : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ZmConvEvapData baseline_data[] = {
      //             pcols, ncol, pver, pverp, time_step
      ZmConvEvapData(    4,    4,   72,    73, 2.0),
      ZmConvEvapData(    4,    4,   72,    73, 3.0),
      ZmConvEvapData(    4,    4,  128,    73, 4.0),
      ZmConvEvapData(    4,    4,  128,    73, 5.0),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmConvEvapData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmConvEvapData test_data[] = {
      ZmConvEvapData(baseline_data[0]),
      ZmConvEvapData(baseline_data[1]),
      ZmConvEvapData(baseline_data[2]),
      ZmConvEvapData(baseline_data[3]),
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
        zm_conv_evap_f(d);
      }
      else {
        zm_conv_evap(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmConvEvapData& d_baseline = baseline_data[i];
        ZmConvEvapData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.tend_s) == d_test.total(d_test.tend_s));
        REQUIRE(d_baseline.total(d_baseline.tend_s) == d_test.total(d_test.tend_q));
        REQUIRE(d_baseline.total(d_baseline.tend_s) == d_test.total(d_test.tend_s_snwprd));
        REQUIRE(d_baseline.total(d_baseline.tend_s) == d_test.total(d_test.tend_s_snwevmlt));
        REQUIRE(d_baseline.total(d_baseline.tend_s) == d_test.total(d_test.ntprprd));
        REQUIRE(d_baseline.total(d_baseline.tend_s) == d_test.total(d_test.ntsnprd));
        for (Int k = 0; k < d_baseline.total(d_baseline.tend_s); ++k) {
          REQUIRE(d_baseline.tend_s[k] == d_test.tend_s[k]);
          REQUIRE(d_baseline.tend_q[k] == d_test.tend_q[k]);
          REQUIRE(d_baseline.tend_s_snwprd[k] == d_test.tend_s_snwprd[k]);
          REQUIRE(d_baseline.tend_s_snwevmlt[k] == d_test.tend_s_snwevmlt[k]);
          REQUIRE(d_baseline.ntprprd[k] == d_test.ntprprd[k]);
          REQUIRE(d_baseline.ntsnprd[k] == d_test.ntsnprd[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.prec));
        REQUIRE(d_baseline.total(d_baseline.prec) == d_test.total(d_test.snow));
        for (Int k = 0; k < d_baseline.total(d_baseline.prec); ++k) {
          REQUIRE(d_baseline.prec[k] == d_test.prec[k]);
          REQUIRE(d_baseline.snow[k] == d_test.snow[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.flxprec) == d_test.total(d_test.flxprec));
        REQUIRE(d_baseline.total(d_baseline.flxprec) == d_test.total(d_test.flxsnow));
        for (Int k = 0; k < d_baseline.total(d_baseline.flxprec); ++k) {
          REQUIRE(d_baseline.flxprec[k] == d_test.flxprec[k]);
          REQUIRE(d_baseline.flxsnow[k] == d_test.flxsnow[k]);
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

TEST_CASE("zm_conv_evap_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmConvEvap;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
