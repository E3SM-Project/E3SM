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
struct UnitWrap::UnitTest<D>::TestGwConvectGwSources : public UnitWrap::UnitTest<D>::Base {

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
    GwConvectGwSourcesData baseline_data[] = {
      //                   ncol, hdepth_min
      GwConvectGwSourcesData(10,     2000., front_init_data[0]),
      GwConvectGwSourcesData(11,     3000., front_init_data[1]),
      GwConvectGwSourcesData(12,     4000., front_init_data[2]),
      GwConvectGwSourcesData(13,     5000., front_init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwConvectGwSourcesData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwConvectGwSourcesData test_data[] = {
      GwConvectGwSourcesData(baseline_data[0]),
      GwConvectGwSourcesData(baseline_data[1]),
      GwConvectGwSourcesData(baseline_data[2]),
      GwConvectGwSourcesData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) {
      gw_convect_gw_sources(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwConvectGwSourcesData& d_baseline = baseline_data[i];
        GwConvectGwSourcesData& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.tau); ++k) {
          REQUIRE(d_baseline.total(d_baseline.tau) == d_test.total(d_test.tau));
          REQUIRE(d_baseline.tau[k] == d_test.tau[k]);
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

TEST_CASE("gw_convect_gw_sources_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwConvectGwSources;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
