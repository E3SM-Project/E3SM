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
struct UnitWrap::UnitTest<D>::TestGwCmSrc : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up init data
    GwFrontInitData front_init_data[] = {
                // taubgnd,  frontgfc_in, kfront_in, init
      GwFrontInitData(  .1,           .4,        10, init_data[0]),
      GwFrontInitData(  .2,           .5,        11, init_data[1]),
      GwFrontInitData(  .3,           .6,        12, init_data[2]),
      GwFrontInitData(  .4,           .7,        13, init_data[3]),
    };

    for (auto& d : front_init_data) {
      d.randomize(engine);
    }

    // Set up inputs
    GwCmSrcData baseline_data[] = {
      //       ncol, kbot
      GwCmSrcData(2,   53, front_init_data[0]),
      GwCmSrcData(3,   54, front_init_data[1]),
      GwCmSrcData(4,   55, front_init_data[2]),
      GwCmSrcData(5,   56, front_init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwCmSrcData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwCmSrcData test_data[] = {
      GwCmSrcData(baseline_data[0]),
      GwCmSrcData(baseline_data[1]),
      GwCmSrcData(baseline_data[2]),
      GwCmSrcData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gw_cm_src(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwCmSrcData& d_baseline = baseline_data[i];
        GwCmSrcData& d_test = test_data[i];
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

TEST_CASE("gw_cm_src_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwCmSrc;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
