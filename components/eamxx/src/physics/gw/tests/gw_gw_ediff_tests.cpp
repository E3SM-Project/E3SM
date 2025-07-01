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
struct UnitWrap::UnitTest<D>::TestGwEdiff : public UnitWrap::UnitTest<D>::Base {

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
    GwEdiffData baseline_data[] = {
      // ncol, ngwv, kbot, ktop, dt
      GwEdiffData(5, 2, 20, 59, 0.1, init_data[0]),
      GwEdiffData(6, 3, 21, 58, 0.2, init_data[1]),
      GwEdiffData(7, 4, 22, 57, 0.3, init_data[2]),
      GwEdiffData(8, 5, 23, 56, 0.4, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwEdiffData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwEdiffData test_data[] = {
      GwEdiffData(baseline_data[0]),
      GwEdiffData(baseline_data[1]),
      GwEdiffData(baseline_data[2]),
      GwEdiffData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gw_ediff(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwEdiffData& d_baseline = baseline_data[i];
        GwEdiffData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.decomp_ca); ++k) {
          REQUIRE(d_baseline.total(d_baseline.decomp_ca) == d_test.total(d_test.decomp_ca));
          REQUIRE(d_baseline.decomp_ca[k] == d_test.decomp_ca[k]);
          REQUIRE(d_baseline.total(d_baseline.decomp_cc) == d_test.total(d_test.decomp_cc));
          REQUIRE(d_baseline.decomp_cc[k] == d_test.decomp_cc[k]);
          REQUIRE(d_baseline.total(d_baseline.decomp_dnom) == d_test.total(d_test.decomp_dnom));
          REQUIRE(d_baseline.decomp_dnom[k] == d_test.decomp_dnom[k]);
          REQUIRE(d_baseline.total(d_baseline.decomp_ze) == d_test.total(d_test.decomp_ze));
          REQUIRE(d_baseline.decomp_ze[k] == d_test.decomp_ze[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.egwdffi); ++k) {
          REQUIRE(d_baseline.total(d_baseline.egwdffi) == d_test.total(d_test.egwdffi));
          REQUIRE(d_baseline.egwdffi[k] == d_test.egwdffi[k]);
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

TEST_CASE("gw_ediff_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwEdiff;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
