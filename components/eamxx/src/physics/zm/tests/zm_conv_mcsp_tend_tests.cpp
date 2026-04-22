#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZmConvMcspTend : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ZmConvMcspTendData baseline_data[] = {
      //             pcols, ncol, pver, pverp, ztodt
      ZmConvMcspTendData(4,    4,   72,    73, 2.0),
      ZmConvMcspTendData(4,    4,   72,    73, 3.0),
      ZmConvMcspTendData(4,    4,   128,  129, 4.0),
      ZmConvMcspTendData(4,    4,   128,  129, 5.0),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmConvMcspTendData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmConvMcspTendData test_data[] = {
      ZmConvMcspTendData(baseline_data[0]),
      ZmConvMcspTendData(baseline_data[1]),
      ZmConvMcspTendData(baseline_data[2]),
      ZmConvMcspTendData(baseline_data[3]),
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
        zm_conv_mcsp_tend_f(d);
      }
      else {
        zm_conv_mcsp_tend(d);
      }
    }

    // zm_conv_mcsp_test does a few sum reductions and we can't guarantee
    // order of operations consistency with fortran, so we need approx
    const auto margin = std::numeric_limits<Real>::epsilon() *
      (ekat::is_single_precision<Real>::value ? 1000 : 1);

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmConvMcspTendData& d_baseline = baseline_data[i];
        ZmConvMcspTendData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.ptend_s) == d_test.total(d_test.ptend_s));
        REQUIRE(d_baseline.total(d_baseline.ptend_s) == d_test.total(d_test.ptend_q));
        REQUIRE(d_baseline.total(d_baseline.ptend_s) == d_test.total(d_test.ptend_u));
        REQUIRE(d_baseline.total(d_baseline.ptend_s) == d_test.total(d_test.ptend_v));
        REQUIRE(d_baseline.total(d_baseline.ptend_s) == d_test.total(d_test.mcsp_dt_out));
        REQUIRE(d_baseline.total(d_baseline.ptend_s) == d_test.total(d_test.mcsp_dq_out));
        REQUIRE(d_baseline.total(d_baseline.ptend_s) == d_test.total(d_test.mcsp_du_out));
        REQUIRE(d_baseline.total(d_baseline.ptend_s) == d_test.total(d_test.mcsp_dv_out));
        for (Int k = 0; k < d_baseline.total(d_baseline.ptend_s); ++k) {
          REQUIRE(d_baseline.ptend_s[k] == Approx(d_test.ptend_s[k]).margin(margin));
          REQUIRE(d_baseline.ptend_q[k] == d_test.ptend_q[k]);
          REQUIRE(d_baseline.ptend_u[k] == d_test.ptend_u[k]);
          REQUIRE(d_baseline.ptend_v[k] == d_test.ptend_v[k]);
          REQUIRE(d_baseline.mcsp_dt_out[k] == Approx(d_test.mcsp_dt_out[k]).margin(margin));
          REQUIRE(d_baseline.mcsp_dq_out[k] == d_test.mcsp_dq_out[k]);
          REQUIRE(d_baseline.mcsp_du_out[k] == d_test.mcsp_du_out[k]);
          REQUIRE(d_baseline.mcsp_dv_out[k] == d_test.mcsp_dv_out[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.mcsp_freq) == d_test.total(d_test.mcsp_freq));
        REQUIRE(d_baseline.total(d_baseline.mcsp_freq) == d_test.total(d_test.mcsp_shear));
        REQUIRE(d_baseline.total(d_baseline.mcsp_freq) == d_test.total(d_test.zm_depth));
        for (Int k = 0; k < d_baseline.total(d_baseline.mcsp_freq); ++k) {
          REQUIRE(d_baseline.mcsp_freq[k] == d_test.mcsp_freq[k]);
          REQUIRE(d_baseline.mcsp_shear[k] == d_test.mcsp_shear[k]);
          REQUIRE(d_baseline.zm_depth[k] == d_test.zm_depth[k]);
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

TEST_CASE("zm_conv_mcsp_tend_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmConvMcspTend;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
