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
struct UnitWrap::UnitTest<D>::TestComputeCapeFromParcel : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ComputeCapeFromParcelData baseline_data[] = {
      //                        pcols, ncol, pver, pverp, num_cin, num_msg
      ComputeCapeFromParcelData(    4,    4,   72,    73,      10, 30),
      ComputeCapeFromParcelData(    4,    4,   72,    73,      10, 40),
      ComputeCapeFromParcelData(    4,    4,   128,  129,      10, 70),
      ComputeCapeFromParcelData(    4,    4,   128,  129,      10, 80),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ComputeCapeFromParcelData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ComputeCapeFromParcelData test_data[] = {
      ComputeCapeFromParcelData(baseline_data[0]),
      ComputeCapeFromParcelData(baseline_data[1]),
      ComputeCapeFromParcelData(baseline_data[2]),
      ComputeCapeFromParcelData(baseline_data[3]),
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
        compute_cape_from_parcel_f(d);
      }
      else {
        compute_cape_from_parcel(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ComputeCapeFromParcelData& d_baseline = baseline_data[i];
        ComputeCapeFromParcelData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.parcel_qsat) == d_test.total(d_test.parcel_qsat));
        REQUIRE(d_baseline.total(d_baseline.parcel_qsat) == d_test.total(d_test.parcel_temp));
        REQUIRE(d_baseline.total(d_baseline.parcel_qsat) == d_test.total(d_test.parcel_vtemp));
        for (Int k = 0; k < d_baseline.total(d_baseline.parcel_qsat); ++k) {
          REQUIRE(d_baseline.parcel_qsat[k] == d_test.parcel_qsat[k]);
          REQUIRE(d_baseline.parcel_temp[k] == d_test.parcel_temp[k]);
          REQUIRE(d_baseline.parcel_vtemp[k] == d_test.parcel_vtemp[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.cape) == d_test.total(d_test.cape));
        REQUIRE(d_baseline.total(d_baseline.cape) == d_test.total(d_test.eql_klev));
        for (Int k = 0; k < d_baseline.total(d_baseline.cape); ++k) {
          REQUIRE(d_baseline.cape[k] == d_test.cape[k]);
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

TEST_CASE("compute_cape_from_parcel_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeCapeFromParcel;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
