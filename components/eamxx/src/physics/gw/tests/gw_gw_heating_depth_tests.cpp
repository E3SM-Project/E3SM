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
struct UnitWrap::UnitTest<D>::TestGwHeatingDepth : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up convect init data
    GwConvectInitData front_init_data[] = {
                  //    maxh,        maxuh, plev_src_wnd, init
      GwConvectInitData(  10,           20,          .10, init_data[0]),
      GwConvectInitData(  11,           21,          .11, init_data[1]),
      GwConvectInitData(  12,           22,          .12, init_data[2]),
      GwConvectInitData(  13,           23,          .13, init_data[3]),
    };

    for (auto& d : front_init_data) {
      d.randomize(engine);
    }

    // Set up inputs
    GwHeatingDepthData baseline_data[] = {
      //               ncol, maxq0_conversion_factor, hdepth_scaling_factor, use_gw_convect_old
      GwHeatingDepthData(2,                    0.42,                  0.43,         false, front_init_data[0]),
      GwHeatingDepthData(3,                    0.43,                  0.44,         false, front_init_data[1]),
      GwHeatingDepthData(4,                    0.44,                  0.45,         true,  front_init_data[2]),
      GwHeatingDepthData(5,                    0.45,                  0.46,         true,  front_init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwHeatingDepthData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine, { {d.zm, {5000., 21000.}},  {d.netdt, {-1.432438750782E-04, 1.432438750782E-04}} });
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwHeatingDepthData test_data[] = {
      // TODO
      GwHeatingDepthData(baseline_data[0]),
      GwHeatingDepthData(baseline_data[1]),
      GwHeatingDepthData(baseline_data[2]),
      GwHeatingDepthData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gw_heating_depth(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwHeatingDepthData& d_baseline = baseline_data[i];
        GwHeatingDepthData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.hdepth); ++k) {
          REQUIRE(d_baseline.total(d_baseline.hdepth) == d_test.total(d_test.hdepth));
          REQUIRE(d_baseline.hdepth[k] == d_test.hdepth[k]);
          REQUIRE(d_baseline.total(d_baseline.hdepth) == d_test.total(d_test.maxq0_out));
          REQUIRE(d_baseline.maxq0_out[k] == d_test.maxq0_out[k]);
          REQUIRE(d_baseline.total(d_baseline.hdepth) == d_test.total(d_test.maxq0));
          REQUIRE(d_baseline.maxq0[k] == d_test.maxq0[k]);
          REQUIRE(d_baseline.total(d_baseline.hdepth) == d_test.total(d_test.mini));
          REQUIRE(d_baseline.mini[k] == d_test.mini[k]);
          REQUIRE(d_baseline.total(d_baseline.hdepth) == d_test.total(d_test.maxi));
          REQUIRE(d_baseline.maxi[k] == d_test.maxi[k]);
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

TEST_CASE("gw_heating_depth_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwHeatingDepth;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
