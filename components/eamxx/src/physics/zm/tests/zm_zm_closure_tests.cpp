#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZmClosure : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ZmClosureData baseline_data[] = {
      //            pcols, ncol, pver, pverp, msg, cape_threshold_in
      ZmClosureData(    1,    1,   72,    73,  10,    5.0),
      ZmClosureData(    1,    1,   72,    73,  10,    6.0),
      ZmClosureData(    1,    1,  128,   129,  10,    7.0),
      ZmClosureData(    1,    1,  128,   129,  10,    8.0),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmClosureData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmClosureData test_data[] = {
      ZmClosureData(baseline_data[0]),
      ZmClosureData(baseline_data[1]),
      ZmClosureData(baseline_data[2]),
      ZmClosureData(baseline_data[3]),
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
        zm_closure_f(d);
      }
      else {
        zm_closure(d);
      }
    }

    // TODO - Remove?
    const auto margin = std::numeric_limits<Real>::epsilon() *
      (ekat::is_single_precision<Real>::value ? 1000 : 1);

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmClosureData& d_baseline = baseline_data[i];
        ZmClosureData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.cld_base_mass_flux) == d_test.total(d_test.cld_base_mass_flux));
        for (Int k = 0; k < d_baseline.total(d_baseline.cld_base_mass_flux); ++k) {
          REQUIRE(d_baseline.cld_base_mass_flux[k] == Approx(d_test.cld_base_mass_flux[k]).margin(margin));
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

TEST_CASE("zm_closure_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmClosure;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
