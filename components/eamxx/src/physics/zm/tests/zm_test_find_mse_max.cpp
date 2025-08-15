#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"
#include "physics/zm/tests/infra/zm_test_data_functions.hpp"

#include "zm_unit_tests_common.hpp"

#include <ekat_pack.hpp>

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::Test_zm_find_mse_max : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up baseline data
    zm_data_find_mse_max baseline_data[] = {
      //                    pcols, ncol, pver, num_msg, pergro_active
      zm_data_find_mse_max( 8,     1,     8,   0,       false ),
      zm_data_find_mse_max( 8,     2,    16,   1,       false ),
      zm_data_find_mse_max( 8,     8,    32,   2,       false ),
      // zm_data_find_mse_max( 8,     8,    64,   3,       true,  ),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(zm_data_find_mse_max);

    // Generate input data - baseline
    for (auto& d : baseline_data) {
      for (auto i = decltype(d.ncol){0}; i < d.ncol; ++i) {
        d.msemax_top_k[i] = 0;
      }
      zm_test_data_generate_profile( engine, d.pver, d.ncol, d.zmid, d.temperature, d.sp_humidity );
    }

    // Create copies of data for use by test
    // (needs to happen before read calls so that inout data is in original state)
    zm_data_find_mse_max test_data[] = {
      zm_data_find_mse_max( baseline_data[0] ),
      zm_data_find_mse_max( baseline_data[1] ),
      zm_data_find_mse_max( baseline_data[2] ),
      // zm_data_find_mse_max(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from test
    for (auto& d : test_data) { zm_find_mse_max(d); }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        zm_data_find_mse_max& d_baseline = baseline_data[i];
        zm_data_find_mse_max& d_test = test_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.msemax_klev); ++k) {
          REQUIRE(d_baseline.msemax_klev[k] == d_test.msemax_klev[k]);
        }
        for (Int k = 0; k < d_baseline.total(d_baseline.mse_max_val); ++k) {
          REQUIRE(d_baseline.mse_max_val[k] == d_test.mse_max_val[k]);
        }
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        test_data[i].write(Base::m_ofile);
      }
    }

  } // zm_find_mse_max

};

} // namespace unit_test
} // namespace zm
} // namespace scream

namespace {

TEST_CASE("zm_find_mse_max", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::Test_zm_find_mse_max;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
