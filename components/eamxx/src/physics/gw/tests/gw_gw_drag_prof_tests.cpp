#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/gw/gw_functions.hpp"
#include "physics/gw/tests/infra/gw_test_data.hpp"

#include "gw_unit_tests_common.hpp"

namespace scream {
namespace gw {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestGwDragProf : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    GwInit init_data[] = {
          // pver, pgwv,   dc, orog_only, molec_diff, tau_0_ubc, nbot_molec, ktop, kbotbg, fcrit2, kwv
      GwInit(  72,   20, 0.75,     false,      false,     false,         16,   60,     16,    .67, 6.28e-5),
      GwInit(  72,   20, 0.75,     true ,      false,     true ,         16,   60,     16,    .67, 6.28e-5),
      GwInit(  72,   20, 0.75,     false,      true ,     true ,         16,   60,     16,    .67, 6.28e-5),
      GwInit(  72,   20, 0.75,     true ,      true ,     false,         16,   60,     16,    .67, 6.28e-5),
    };

    for (auto& d : init_data) {
      d.randomize(engine);
    }

    // Set up inputs
    GwDragProfData baseline_data[] = {
      GwDragProfData(5, 2, 10, true, .4, 1.8, init_data[0]),
      GwDragProfData(6, 3, 11, false, .8, 2.4, init_data[1]),
      GwDragProfData(7, 4, 12, true, 1.4, 3.4, init_data[2]),
      GwDragProfData(8, 5, 13, false, 2.4, 4.4, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwDragProfData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwDragProfData test_data[] = {
      GwDragProfData(baseline_data[0]),
      GwDragProfData(baseline_data[1]),
      GwDragProfData(baseline_data[2]),
      GwDragProfData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_fid);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gw_drag_prof(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwDragProfData& d_baseline = baseline_data[i];
        GwDragProfData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.tau); ++k) {
          REQUIRE(d_baseline.total(d_baseline.tau) == d_test.total(d_test.tau));
          REQUIRE(d_baseline.tau[k] == d_test.tau[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.utgw); ++k) {
          REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.utgw));
          REQUIRE(d_baseline.utgw[k] == d_test.utgw[k]);
          REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.vtgw));
          REQUIRE(d_baseline.vtgw[k] == d_test.vtgw[k]);
          REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.ttgw));
          REQUIRE(d_baseline.ttgw[k] == d_test.ttgw[k]);
          REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.dttdf));
          REQUIRE(d_baseline.dttdf[k] == d_test.dttdf[k]);
          REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.dttke));
          REQUIRE(d_baseline.dttke[k] == d_test.dttke[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.qtgw); ++k) {
          REQUIRE(d_baseline.total(d_baseline.qtgw) == d_test.total(d_test.qtgw));
          REQUIRE(d_baseline.qtgw[k] == d_test.qtgw[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.taucd); ++k) {
          REQUIRE(d_baseline.total(d_baseline.taucd) == d_test.total(d_test.taucd));
          REQUIRE(d_baseline.taucd[k] == d_test.taucd[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.egwdffi); ++k) {
          REQUIRE(d_baseline.total(d_baseline.egwdffi) == d_test.total(d_test.egwdffi));
          REQUIRE(d_baseline.egwdffi[k] == d_test.egwdffi[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.gwut); ++k) {
          REQUIRE(d_baseline.total(d_baseline.gwut) == d_test.total(d_test.gwut));
          REQUIRE(d_baseline.gwut[k] == d_test.gwut[k]);
        }

      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        test_data[i].write(Base::m_fid);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace gw
} // namespace scream

namespace {

TEST_CASE("gw_drag_prof_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwDragProf;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
