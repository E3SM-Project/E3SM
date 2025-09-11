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
struct UnitWrap::UnitTest<D>::TestGwOroSrc : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up inputs
    GwOroSrcData baseline_data[] = {
      //        ncol
      GwOroSrcData(2, init_data[0]),
      GwOroSrcData(3, init_data[1]),
      GwOroSrcData(4, init_data[2]),
      GwOroSrcData(5, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwOroSrcData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwOroSrcData test_data[] = {
      GwOroSrcData(baseline_data[0]),
      GwOroSrcData(baseline_data[1]),
      GwOroSrcData(baseline_data[2]),
      GwOroSrcData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gw_oro_src(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwOroSrcData& d_baseline = baseline_data[i];
        GwOroSrcData& d_test = test_data[i];
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

TEST_CASE("gw_oro_src_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwOroSrc;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
