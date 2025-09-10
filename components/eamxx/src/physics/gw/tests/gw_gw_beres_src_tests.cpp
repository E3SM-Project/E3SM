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
struct UnitWrap::UnitTest<D>::TestGwBeresSrc : public UnitWrap::UnitTest<D>::Base {

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
    GwBeresSrcData baseline_data[] = {
      //           ncol, maxq0_conversion_factor, hdepth_scaling_factor, hdepth_min, storm_speed_min, use_gw_convect_old
      GwBeresSrcData(10,                     0.1,                   0.2,      2000.,              1., false,  front_init_data[0]),
      GwBeresSrcData(11,                     0.2,                   0.3,      3000.,              2., false,  front_init_data[1]),
      GwBeresSrcData(12,                     0.3,                   0.4,      4000.,              3., true,   front_init_data[2]),
      GwBeresSrcData(13,                     0.4,                   0.5,      5000.,              4., true,   front_init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwBeresSrcData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine, {
          {d.lat, {-1., 1.}},
          {d.netdt, {-1.432438750782E-04, 1.432438750782E-04}},
          {d.u, {-6.355356031720E+01, 2.735168920349E+01}},
          {d.v, {-5.154620900385E+00, 2.358902991534E+00}},
          {d.zm, {5.444482895112E+01, 6.240534656097E+04}},
        });
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwBeresSrcData test_data[] = {
      GwBeresSrcData(baseline_data[0]),
      GwBeresSrcData(baseline_data[1]),
      GwBeresSrcData(baseline_data[2]),
      GwBeresSrcData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gw_beres_src(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwBeresSrcData& d_baseline = baseline_data[i];
        GwBeresSrcData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.tau); ++k) {
          REQUIRE(d_baseline.total(d_baseline.tau) == d_test.total(d_test.tau));
          REQUIRE(d_baseline.tau[k] == d_test.tau[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.ubm); ++k) {
          REQUIRE(d_baseline.total(d_baseline.ubm) == d_test.total(d_test.ubm));
          REQUIRE(d_baseline.ubm[k] == d_test.ubm[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.ubi); ++k) {
          REQUIRE(d_baseline.total(d_baseline.ubi) == d_test.total(d_test.ubi));
          REQUIRE(d_baseline.ubi[k] == d_test.ubi[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.xv); ++k) {
          REQUIRE(d_baseline.total(d_baseline.xv) == d_test.total(d_test.xv));
          REQUIRE(d_baseline.xv[k] == d_test.xv[k]);
          REQUIRE(d_baseline.total(d_baseline.xv) == d_test.total(d_test.yv));
          REQUIRE(d_baseline.yv[k] == d_test.yv[k]);
          REQUIRE(d_baseline.total(d_baseline.xv) == d_test.total(d_test.hdepth));
          REQUIRE(d_baseline.hdepth[k] == d_test.hdepth[k]);
          REQUIRE(d_baseline.total(d_baseline.xv) == d_test.total(d_test.maxq0_out));
          REQUIRE(d_baseline.maxq0_out[k] == d_test.maxq0_out[k]);
          REQUIRE(d_baseline.total(d_baseline.xv) == d_test.total(d_test.src_level));
          REQUIRE(d_baseline.src_level[k] == d_test.src_level[k]);
          REQUIRE(d_baseline.total(d_baseline.xv) == d_test.total(d_test.tend_level));
          REQUIRE(d_baseline.tend_level[k] == d_test.tend_level[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.c); ++k) {
          REQUIRE(d_baseline.total(d_baseline.c) == d_test.total(d_test.c));
          REQUIRE(d_baseline.c[k] == d_test.c[k]);
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

TEST_CASE("gw_beres_src_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwBeresSrc;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
