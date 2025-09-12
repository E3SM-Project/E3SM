#include "catch2/catch.hpp"

#include "physics/gw/gw_functions.hpp"
#include "physics/gw/tests/infra/gw_test_data.hpp"
#include "share/core/eamxx_types.hpp"

#include "gw_unit_tests_common.hpp"

#include <ekat_pack.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {
namespace gw {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestVdLuDecomp : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up inputs
    VdLuDecompData baseline_data[] = {
      // ncol, ntop, nbot, ztodt
      VdLuDecompData(2, 20, 50, 0.5, init_data[0]),
      VdLuDecompData(3, 19, 51, 0.5, init_data[1]),
      VdLuDecompData(4, 18, 52, 0.5, init_data[2]),
      VdLuDecompData(5, 17, 53, 0.5, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(VdLuDecompData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    VdLuDecompData test_data[] = {
      VdLuDecompData(baseline_data[0]),
      VdLuDecompData(baseline_data[1]),
      VdLuDecompData(baseline_data[2]),
      VdLuDecompData(baseline_data[3]),
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
        vd_lu_decomp_f(d);
      }
      else {
        vd_lu_decomp(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        VdLuDecompData& d_baseline = baseline_data[i];
        VdLuDecompData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.decomp_ca) == d_test.total(d_test.decomp_ca));
        REQUIRE(d_baseline.total(d_baseline.decomp_ca) == d_test.total(d_test.decomp_cc));
        REQUIRE(d_baseline.total(d_baseline.decomp_ca) == d_test.total(d_test.decomp_dnom));
        REQUIRE(d_baseline.total(d_baseline.decomp_ca) == d_test.total(d_test.decomp_ze));
        for (Int k = 0; k < d_baseline.total(d_baseline.decomp_ca); ++k) {
          REQUIRE(d_baseline.decomp_ca[k] == d_test.decomp_ca[k]);
          REQUIRE(d_baseline.decomp_cc[k] == d_test.decomp_cc[k]);
          REQUIRE(d_baseline.decomp_dnom[k] == d_test.decomp_dnom[k]);
          REQUIRE(d_baseline.decomp_ze[k] == d_test.decomp_ze[k]);
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

TEST_CASE("vd_lu_decomp_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestVdLuDecomp;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
