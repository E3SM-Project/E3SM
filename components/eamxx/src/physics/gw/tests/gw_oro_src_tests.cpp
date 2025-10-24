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
      d.randomize(engine, { {d.sgh, {2, 7}} });
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
      if (this->m_baseline_action == GENERATE) {
        gw_oro_src_f(d);
      }
      else {
        gw_oro_src(d);
      }
    }

    // We need a tolerance since the order of operations is different from f90.
    // This tol can be removed once we are no longer using
    // fortran to generate baselines.
    const auto margin = std::numeric_limits<Real>::epsilon() *
      (ekat::is_single_precision<Real>::value ? 1000 : 1);


    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwOroSrcData& d_baseline = baseline_data[i];
        GwOroSrcData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.tau) == d_test.total(d_test.tau));
        for (Int k = 0; k < d_baseline.total(d_baseline.tau); ++k) {
          REQUIRE(d_baseline.tau[k] == Approx(d_test.tau[k]).margin(margin));
        }
        REQUIRE(d_baseline.total(d_baseline.ubm) == d_test.total(d_test.ubm));
        for (Int k = 0; k < d_baseline.total(d_baseline.ubm); ++k) {
          REQUIRE(d_baseline.ubm[k] == Approx(d_test.ubm[k]).margin(margin));
        }
        REQUIRE(d_baseline.total(d_baseline.ubi) == d_test.total(d_test.ubi));
        for (Int k = 0; k < d_baseline.total(d_baseline.ubi); ++k) {
          REQUIRE(d_baseline.ubi[k] == Approx(d_test.ubi[k]).margin(margin));
        }
        REQUIRE(d_baseline.total(d_baseline.xv) == d_test.total(d_test.xv));
        REQUIRE(d_baseline.total(d_baseline.xv) == d_test.total(d_test.yv));
        REQUIRE(d_baseline.total(d_baseline.xv) == d_test.total(d_test.src_level));
        REQUIRE(d_baseline.total(d_baseline.xv) == d_test.total(d_test.tend_level));
        for (Int k = 0; k < d_baseline.total(d_baseline.xv); ++k) {
          REQUIRE(d_baseline.xv[k] == Approx(d_test.xv[k]).margin(margin));
          REQUIRE(d_baseline.yv[k] == Approx(d_test.yv[k]).margin(margin));
          REQUIRE(d_baseline.src_level[k] == d_test.src_level[k]);
          REQUIRE(d_baseline.tend_level[k] == d_test.tend_level[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.c) == d_test.total(d_test.c));
        for (Int k = 0; k < d_baseline.total(d_baseline.c); ++k) {
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
