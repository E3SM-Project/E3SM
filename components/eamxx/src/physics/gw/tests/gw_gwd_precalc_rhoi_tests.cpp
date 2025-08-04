#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "physics/gw/gw_functions.hpp"
#include "physics/gw/tests/infra/gw_test_data.hpp"

#include "gw_unit_tests_common.hpp"

#include <ekat_pack.hpp>

namespace scream {
namespace gw {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestGwdPrecalcRhoi : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up inputs
    GwdPrecalcRhoiData baseline_data[] = {
      //             pcnst, ncol,  dt
      GwdPrecalcRhoiData(5,    2,  .4, init_data[0]),
      GwdPrecalcRhoiData(6,    3,  .8, init_data[1]),
      GwdPrecalcRhoiData(7,    4, 1.4, init_data[2]),
      GwdPrecalcRhoiData(8,    5, 2.4, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwdPrecalcRhoiData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwdPrecalcRhoiData test_data[] = {
      GwdPrecalcRhoiData(baseline_data[0]),
      GwdPrecalcRhoiData(baseline_data[1]),
      GwdPrecalcRhoiData(baseline_data[2]),
      GwdPrecalcRhoiData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gwd_precalc_rhoi(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwdPrecalcRhoiData& d_baseline = baseline_data[i];
        GwdPrecalcRhoiData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.egwdffi); ++k) {
          REQUIRE(d_baseline.total(d_baseline.egwdffi) == d_test.total(d_test.egwdffi));
          REQUIRE(d_baseline.egwdffi[k] == d_test.egwdffi[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.qtgw); ++k) {
          REQUIRE(d_baseline.total(d_baseline.qtgw) == d_test.total(d_test.qtgw));
          REQUIRE(d_baseline.qtgw[k] == d_test.qtgw[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.dttdf); ++k) {
          REQUIRE(d_baseline.total(d_baseline.dttdf) == d_test.total(d_test.dttdf));
          REQUIRE(d_baseline.dttdf[k] == d_test.dttdf[k]);
          REQUIRE(d_baseline.total(d_baseline.dttdf) == d_test.total(d_test.dttke));
          REQUIRE(d_baseline.dttke[k] == d_test.dttke[k]);
          REQUIRE(d_baseline.total(d_baseline.dttdf) == d_test.total(d_test.ttgw));
          REQUIRE(d_baseline.ttgw[k] == d_test.ttgw[k]);
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
} // namespace gw
} // namespace scream

namespace {

TEST_CASE("gwd_precalc_rhoi_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwdPrecalcRhoi;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
