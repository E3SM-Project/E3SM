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
struct UnitWrap::UnitTest<D>::TestComputeDiluteParcel : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ComputeDiluteParcelData baseline_data[] = {
      //                      pcols, ncol, pver, num_msg
      ComputeDiluteParcelData(    4,    4,   72,     30),
      ComputeDiluteParcelData(    4,    4,   72,     20),
      ComputeDiluteParcelData(    4,    4,  128,     40),
      ComputeDiluteParcelData(    4,    4,  128,     80),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ComputeDiluteParcelData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ComputeDiluteParcelData test_data[] = {
      ComputeDiluteParcelData(baseline_data[0]),
      ComputeDiluteParcelData(baseline_data[1]),
      ComputeDiluteParcelData(baseline_data[2]),
      ComputeDiluteParcelData(baseline_data[3]),
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
        compute_dilute_parcel_f(d);
      }
      else {
        compute_dilute_parcel(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ComputeDiluteParcelData& d_baseline = baseline_data[i];
        ComputeDiluteParcelData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.parcel_temp) == d_test.total(d_test.parcel_temp));
        REQUIRE(d_baseline.total(d_baseline.parcel_temp) == d_test.total(d_test.parcel_vtemp));
        REQUIRE(d_baseline.total(d_baseline.parcel_temp) == d_test.total(d_test.parcel_qsat));
        for (Int k = 0; k < d_baseline.total(d_baseline.parcel_temp); ++k) {
          REQUIRE(d_baseline.parcel_temp[k] == d_test.parcel_temp[k]);
          REQUIRE(d_baseline.parcel_vtemp[k] == d_test.parcel_vtemp[k]);
          REQUIRE(d_baseline.parcel_qsat[k] == d_test.parcel_qsat[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.lcl_pmid) == d_test.total(d_test.lcl_pmid));
        REQUIRE(d_baseline.total(d_baseline.lcl_pmid) == d_test.total(d_test.lcl_temperature));
        REQUIRE(d_baseline.total(d_baseline.lcl_pmid) == d_test.total(d_test.lcl_klev));
        for (Int k = 0; k < d_baseline.total(d_baseline.lcl_pmid); ++k) {
          REQUIRE(d_baseline.lcl_pmid[k] == d_test.lcl_pmid[k]);
          REQUIRE(d_baseline.lcl_temperature[k] == d_test.lcl_temperature[k]);
          REQUIRE(d_baseline.lcl_klev[k] == d_test.lcl_klev[k]);
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

TEST_CASE("compute_dilute_parcel_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeDiluteParcel;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
