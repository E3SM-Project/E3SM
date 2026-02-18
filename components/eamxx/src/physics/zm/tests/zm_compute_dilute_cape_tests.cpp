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
struct UnitWrap::UnitTest<D>::TestComputeDiluteCape : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ComputeDiluteCapeData baseline_data[] = {
      //                    pcols, ncol, pver, pverp, num_cin, num_msg, calc_msemax_klev, use_input_tq_mx
      ComputeDiluteCapeData(    4,    4,   72,    73,      10,      30,            false, false),
      ComputeDiluteCapeData(    4,    4,   72,    73,      10,      30,            false, true),
      ComputeDiluteCapeData(    4,    4,   72,    73,      10,      30,            true,  false),
      ComputeDiluteCapeData(    4,    4,   72,    73,      10,      30,            true,  true),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ComputeDiluteCapeData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ComputeDiluteCapeData test_data[] = {
      ComputeDiluteCapeData(baseline_data[0]),
      ComputeDiluteCapeData(baseline_data[1]),
      ComputeDiluteCapeData(baseline_data[2]),
      ComputeDiluteCapeData(baseline_data[3]),
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
        compute_dilute_cape_f(d);
      }
      else {
        compute_dilute_cape(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ComputeDiluteCapeData& d_baseline = baseline_data[i];
        ComputeDiluteCapeData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.parcel_temp) == d_test.total(d_test.parcel_temp));
        REQUIRE(d_baseline.total(d_baseline.parcel_temp) == d_test.total(d_test.parcel_qsat));
        for (Int k = 0; k < d_baseline.total(d_baseline.parcel_temp); ++k) {
          REQUIRE(d_baseline.parcel_temp[k] == d_test.parcel_temp[k]);
          REQUIRE(d_baseline.parcel_qsat[k] == d_test.parcel_qsat[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.lcl_temperature) == d_test.total(d_test.lcl_temperature));
        REQUIRE(d_baseline.total(d_baseline.lcl_temperature) == d_test.total(d_test.cape));
        REQUIRE(d_baseline.total(d_baseline.lcl_temperature) == d_test.total(d_test.q_mx));
        REQUIRE(d_baseline.total(d_baseline.lcl_temperature) == d_test.total(d_test.t_mx));
        REQUIRE(d_baseline.total(d_baseline.lcl_temperature) == d_test.total(d_test.msemax_klev));
        REQUIRE(d_baseline.total(d_baseline.lcl_temperature) == d_test.total(d_test.lcl_klev));
        REQUIRE(d_baseline.total(d_baseline.lcl_temperature) == d_test.total(d_test.eql_klev));
        for (Int k = 0; k < d_baseline.total(d_baseline.lcl_temperature); ++k) {
          REQUIRE(d_baseline.lcl_temperature[k] == d_test.lcl_temperature[k]);
          REQUIRE(d_baseline.cape[k] == d_test.cape[k]);
          REQUIRE(d_baseline.q_mx[k] == d_test.q_mx[k]);
          REQUIRE(d_baseline.t_mx[k] == d_test.t_mx[k]);
          REQUIRE(d_baseline.msemax_klev[k] == d_test.msemax_klev[k]);
          REQUIRE(d_baseline.lcl_klev[k] == d_test.lcl_klev[k]);
          REQUIRE(d_baseline.eql_klev[k] == d_test.eql_klev[k]);
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

TEST_CASE("compute_dilute_cape_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeDiluteCape;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
