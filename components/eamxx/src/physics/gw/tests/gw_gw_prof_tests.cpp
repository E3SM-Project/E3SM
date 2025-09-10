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
struct UnitWrap::UnitTest<D>::TestGwProf : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up inputs
    GwProfData baseline_data[] = {
      //      ncol, cpair
      GwProfData(2,  .4, init_data[0]),
      GwProfData(3,  .8, init_data[1]),
      GwProfData(4, 1.4, init_data[2]),
      GwProfData(5, 2.4, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(GwProfData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    GwProfData test_data[] = {
      GwProfData(baseline_data[0]),
      GwProfData(baseline_data[1]),
      GwProfData(baseline_data[2]),
      GwProfData(baseline_data[3]),
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
        gw_prof_f(d);
      }
      else {
        gw_prof(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        GwProfData& d_baseline = baseline_data[i];
        GwProfData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.rhoi) == d_test.total(d_test.rhoi));
        REQUIRE(d_baseline.total(d_baseline.rhoi) == d_test.total(d_test.ti));
        REQUIRE(d_baseline.total(d_baseline.rhoi) == d_test.total(d_test.ni));
        for (Int k = 0; k < d_baseline.total(d_baseline.rhoi); ++k) {
          REQUIRE(d_baseline.rhoi[k] == d_test.rhoi[k]);
          REQUIRE(d_baseline.ti[k] == d_test.ti[k]);
          REQUIRE(d_baseline.ni[k] == d_test.ni[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.nm) == d_test.total(d_test.nm));
        for (Int k = 0; k < d_baseline.total(d_baseline.nm); ++k) {
          REQUIRE(d_baseline.nm[k] == d_test.nm[k]);
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

TEST_CASE("gw_prof_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwProf;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
