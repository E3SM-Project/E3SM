#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
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
    auto init_data = get_common_init_data(engine);

    // Set up inputs
    GwEdiffData baseline_data[] = {
      //       ncol, kbot, ktop,  dt
      GwEdiffData(5,   59,   20, 0.1, init_data[0]),
      GwEdiffData(6,   58,   19, 0.2, init_data[1]),
      GwEdiffData(7,   57,   18, 0.3, init_data[2]),
      GwEdiffData(8,   56,   17, 0.4, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwEdiffData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (Int i = 0; i < num_runs; ++i) {
      auto& d = baseline_data[i];
      d.randomize(engine, { {d.tend_level, {d.ktop+1, d.kbot-1}} });
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
      if (this->m_baseline_action == GENERATE) {
        gw_ediff_f(d);
      }
      else {
        gw_ediff(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwEdiffData& d_baseline = baseline_data[i];
        GwEdiffData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.decomp_ca) == d_test.total(d_test.decomp_ca));
        REQUIRE(d_baseline.total(d_baseline.decomp_cc) == d_test.total(d_test.decomp_cc));
        REQUIRE(d_baseline.total(d_baseline.decomp_dnom) == d_test.total(d_test.decomp_dnom));
        REQUIRE(d_baseline.total(d_baseline.decomp_ze) == d_test.total(d_test.decomp_ze));
        for (Int k = 0; k < d_baseline.total(d_baseline.decomp_ca); ++k) {
          REQUIRE(d_baseline.decomp_ca[k] == d_test.decomp_ca[k]);
          REQUIRE(d_baseline.decomp_cc[k] == d_test.decomp_cc[k]);
          REQUIRE(d_baseline.decomp_dnom[k] == d_test.decomp_dnom[k]);
          REQUIRE(d_baseline.decomp_ze[k] == d_test.decomp_ze[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.egwdffi) == d_test.total(d_test.egwdffi));
        for (Int k = 0; k < d_baseline.total(d_baseline.egwdffi); ++k) {
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
