#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

#include <ekat_pack.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestIentropy : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    IentropyData baseline_data[] = {
      //              s,    p,  qt, tfg, t (output, set to zero), qst (output, set to zero)
      IentropyData(3000, 1000, .01, 270, 0                      , 0),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(IentropyData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    IentropyData test_data[] = {
      IentropyData(baseline_data[0]),
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
        ientropy_f(d);
      }
      else {
        ientropy(d);
      }
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        IentropyData& d_baseline = baseline_data[i];
        IentropyData& d_test = test_data[i];
        REQUIRE(d_baseline.t == d_test.t);
        REQUIRE(d_baseline.qst == d_test.qst);
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

TEST_CASE("ientropy_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIentropy;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
