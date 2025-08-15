#include "catch2/catch.hpp"

#include "tms_unit_tests_common.hpp"
#include "tms_functions.hpp"
#include "tms_test_data.hpp"

#include "share/util/eamxx_setup_random_test.hpp"
#include "share/eamxx_types.hpp"

namespace scream {
namespace tms {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestComputeTMS : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    // Should property tests be created?
  } // run_property

  void run_bfb()
  {
    auto engine = Base::get_engine();

    ComputeTMSData baseline_data[] = {
      //             ncols, nlevs
      ComputeTMSData(12,    72),
      ComputeTMSData(8,     12),
      ComputeTMSData(7,     16),
      ComputeTMSData(2,     7)
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ComputeTMSData);

    // Generate random input data
    for (auto& d : baseline_data) {
      d.randomize(engine, { {d.sgh, {0.5, 1.5}} });
    }

    // Create copies of data for use by cxx. Needs to happen before read calls so that
    // inout data is in original state
    ComputeTMSData cxx_data[] = {
      ComputeTMSData(baseline_data[0]),
      ComputeTMSData(baseline_data[1]),
      ComputeTMSData(baseline_data[2]),
      ComputeTMSData(baseline_data[3])
    };

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      compute_tms(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (int r = 0; r<num_runs; ++r) {
        ComputeTMSData& d_baseline = baseline_data[r];
        ComputeTMSData& d_cxx = cxx_data[r];
        REQUIRE(d_baseline.total(d_baseline.ksrf) == d_cxx.total(d_cxx.ksrf));
        REQUIRE(d_baseline.total(d_baseline.taux) == d_cxx.total(d_cxx.taux));
        REQUIRE(d_baseline.total(d_baseline.tauy) == d_cxx.total(d_cxx.tauy));

        for (int i=0; i<d_baseline.total(d_baseline.ksrf); ++i) {
          REQUIRE(d_baseline.ksrf[i] == d_cxx.ksrf[i]);
          REQUIRE(d_baseline.taux[i] == d_cxx.taux[i]);
          REQUIRE(d_baseline.tauy[i] == d_cxx.tauy[i]);
        }
      }
    }
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        cxx_data[i].write(Base::m_ofile);
      }
    }
  } // run_bfb
};

} // namespace unit_test
} // namespace tms
} // namespace scream

namespace {

TEST_CASE("compute_tms_bfb", "tms")
{
  using TestStruct = scream::tms::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeTMS;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
