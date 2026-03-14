#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZmConvMcspCalculateShear : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ZmConvMcspCalculateShearData baseline_data[] = {
      //                       pcols, ncol, pver
      ZmConvMcspCalculateShearData(4,    4,  72),
      ZmConvMcspCalculateShearData(4,    4,  128),
      ZmConvMcspCalculateShearData(4,    4,  72),
      ZmConvMcspCalculateShearData(4,    4,  128),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmConvMcspCalculateShearData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmConvMcspCalculateShearData test_data[] = {
      ZmConvMcspCalculateShearData(baseline_data[0]),
      ZmConvMcspCalculateShearData(baseline_data[1]),
      ZmConvMcspCalculateShearData(baseline_data[2]),
      ZmConvMcspCalculateShearData(baseline_data[3]),
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
        zm_conv_mcsp_calculate_shear_f(d);
      }
      else {
        zm_conv_mcsp_calculate_shear(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmConvMcspCalculateShearData& d_baseline = baseline_data[i];
        ZmConvMcspCalculateShearData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.mcsp_shear) == d_test.total(d_test.mcsp_shear));
        for (Int k = 0; k < d_baseline.total(d_baseline.mcsp_shear); ++k) {
          REQUIRE(d_baseline.mcsp_shear[k] == d_test.mcsp_shear[k]);
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

TEST_CASE("zm_conv_mcsp_calculate_shear_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmConvMcspCalculateShear;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
