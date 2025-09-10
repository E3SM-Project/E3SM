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
struct UnitWrap::UnitTest<D>::TestGwStormSpeed : public UnitWrap::UnitTest<D>::Base {

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
    GwStormSpeedData baseline_data[] = {
      //             ncol, storm_speed_min
      GwStormSpeedData(10,              .5, front_init_data[0]),
      GwStormSpeedData(11,             1.5, front_init_data[1]),
      GwStormSpeedData(12,             2.5, front_init_data[2]),
      GwStormSpeedData(13,             3.5, front_init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwStormSpeedData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine, { {d.mini, {10, 20}}, {d.maxi, {40, 50}}, {d.ubm, {5., 20.}} });
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwStormSpeedData test_data[] = {
      GwStormSpeedData(baseline_data[0]),
      GwStormSpeedData(baseline_data[1]),
      GwStormSpeedData(baseline_data[2]),
      GwStormSpeedData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gw_storm_speed(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwStormSpeedData& d_baseline = baseline_data[i];
        GwStormSpeedData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.uh); ++k) {
          REQUIRE(d_baseline.total(d_baseline.uh) == d_test.total(d_test.uh));
          REQUIRE(d_baseline.uh[k] == d_test.uh[k]);
          REQUIRE(d_baseline.total(d_baseline.uh) == d_test.total(d_test.umin));
          REQUIRE(d_baseline.umin[k] == d_test.umin[k]);
          REQUIRE(d_baseline.total(d_baseline.uh) == d_test.total(d_test.umax));
          REQUIRE(d_baseline.umax[k] == d_test.umax[k]);
          REQUIRE(d_baseline.total(d_baseline.uh) == d_test.total(d_test.storm_speed));
          REQUIRE(d_baseline.storm_speed[k] == d_test.storm_speed[k]);
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

TEST_CASE("gw_storm_speed_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwStormSpeed;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
