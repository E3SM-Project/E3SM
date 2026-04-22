#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZmCalcOutputTend : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ZmCalcOutputTendData baseline_data[] = {
      //                   pcols, ncol, pver, pverp, msg
      ZmCalcOutputTendData(    4,    4,   72,    73,  10),
      ZmCalcOutputTendData(    4,    4,   72,    73,  20),
      ZmCalcOutputTendData(    4,    4,  128,   129,  30),
      ZmCalcOutputTendData(    4,    4,  128,   129,  40),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmCalcOutputTendData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmCalcOutputTendData test_data[] = {
      ZmCalcOutputTendData(baseline_data[0]),
      ZmCalcOutputTendData(baseline_data[1]),
      ZmCalcOutputTendData(baseline_data[2]),
      ZmCalcOutputTendData(baseline_data[3]),
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
        zm_calc_output_tend_f(d);
      }
      else {
        zm_calc_output_tend(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmCalcOutputTendData& d_baseline = baseline_data[i];
        ZmCalcOutputTendData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.dsdt) == d_test.total(d_test.dsdt));
        REQUIRE(d_baseline.total(d_baseline.dsdt) == d_test.total(d_test.dqdt));
        REQUIRE(d_baseline.total(d_baseline.dsdt) == d_test.total(d_test.dl));
        for (Int k = 0; k < d_baseline.total(d_baseline.dsdt); ++k) {
          REQUIRE(d_baseline.dsdt[k] == d_test.dsdt[k]);
          REQUIRE(d_baseline.dqdt[k] == d_test.dqdt[k]);
          REQUIRE(d_baseline.dl[k] == d_test.dl[k]);
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

TEST_CASE("zm_calc_output_tend_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmCalcOutputTend;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
