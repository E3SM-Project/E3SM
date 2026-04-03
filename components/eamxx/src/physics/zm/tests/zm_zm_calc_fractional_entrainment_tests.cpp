#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZmCalcFractionalEntrainment : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ZmCalcFractionalEntrainmentData baseline_data[] = {
      //                              pcols, ncol, pver, pverp, msg
      ZmCalcFractionalEntrainmentData(    4,    4,   72,    73, 10),
      ZmCalcFractionalEntrainmentData(    4,    4,   72,    73, 20),
      ZmCalcFractionalEntrainmentData(    4,    4,  128,   129, 30),
      ZmCalcFractionalEntrainmentData(    4,    4,  128,   129, 40),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmCalcFractionalEntrainmentData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmCalcFractionalEntrainmentData test_data[] = {
      ZmCalcFractionalEntrainmentData(baseline_data[0]),
      ZmCalcFractionalEntrainmentData(baseline_data[1]),
      ZmCalcFractionalEntrainmentData(baseline_data[2]),
      ZmCalcFractionalEntrainmentData(baseline_data[3]),
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
        zm_calc_fractional_entrainment_f(d);
      }
      else {
        zm_calc_fractional_entrainment(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmCalcFractionalEntrainmentData& d_baseline = baseline_data[i];
        ZmCalcFractionalEntrainmentData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.h_env_min) == d_test.total(d_test.h_env_min));
        REQUIRE(d_baseline.total(d_baseline.h_env_min) == d_test.total(d_test.lambda_max));
        REQUIRE(d_baseline.total(d_baseline.h_env_min) == d_test.total(d_test.j0));
        for (Int k = 0; k < d_baseline.total(d_baseline.h_env_min); ++k) {
          REQUIRE(d_baseline.h_env_min[k] == d_test.h_env_min[k]);
          REQUIRE(d_baseline.lambda_max[k] == d_test.lambda_max[k]);
          REQUIRE(d_baseline.j0[k] == d_test.j0[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.lambda) == d_test.total(d_test.lambda));
        for (Int k = 0; k < d_baseline.total(d_baseline.lambda); ++k) {
          REQUIRE(d_baseline.lambda[k] == d_test.lambda[k]);
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

TEST_CASE("zm_calc_fractional_entrainment_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmCalcFractionalEntrainment;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
