#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

namespace scream {
namespace zm {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestZmDowndraftProperties : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up inputs
    ZmDowndraftPropertiesData baseline_data[] = {
      //                        pcols, ncol, pver, pverp, msg
      ZmDowndraftPropertiesData(    4,    4,   72,    73, 10),
      ZmDowndraftPropertiesData(    4,    4,   72,    73, 10),
      ZmDowndraftPropertiesData(    4,    4,  128,   129, 10),
      ZmDowndraftPropertiesData(    4,    4,  128,   129, 10),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ZmDowndraftPropertiesData);

    // Generate random input data
    // Alternatively, you can use the baseline_data constructors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    ZmDowndraftPropertiesData test_data[] = {
      ZmDowndraftPropertiesData(baseline_data[0]),
      ZmDowndraftPropertiesData(baseline_data[1]),
      ZmDowndraftPropertiesData(baseline_data[2]),
      ZmDowndraftPropertiesData(baseline_data[3]),
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
        zm_downdraft_properties_f(d);
      }
      else {
        zm_downdraft_properties(d);
      }
    }

    // TODO - Remove?
    const auto margin = std::numeric_limits<Real>::epsilon() *
      (ekat::is_single_precision<Real>::value ? 1000 : 1);

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmDowndraftPropertiesData& d_baseline = baseline_data[i];
        ZmDowndraftPropertiesData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.mflx_dn) == d_test.total(d_test.mflx_dn));
        REQUIRE(d_baseline.total(d_baseline.mflx_dn) == d_test.total(d_test.entr_dn));
        REQUIRE(d_baseline.total(d_baseline.mflx_dn) == d_test.total(d_test.s_dnd));
        REQUIRE(d_baseline.total(d_baseline.mflx_dn) == d_test.total(d_test.q_dnd));
        REQUIRE(d_baseline.total(d_baseline.mflx_dn) == d_test.total(d_test.h_dnd));
        REQUIRE(d_baseline.total(d_baseline.mflx_dn) == d_test.total(d_test.q_dnd_sat));
        REQUIRE(d_baseline.total(d_baseline.mflx_dn) == d_test.total(d_test.evp));
        for (Int k = 0; k < d_baseline.total(d_baseline.mflx_dn); ++k) {
          REQUIRE(d_baseline.mflx_dn[k] == d_test.mflx_dn[k]);
          REQUIRE(d_baseline.entr_dn[k] == d_test.entr_dn[k]);
          REQUIRE(d_baseline.s_dnd[k] == d_test.s_dnd[k]);
          REQUIRE(d_baseline.q_dnd[k] == d_test.q_dnd[k]);
          REQUIRE(d_baseline.h_dnd[k] == d_test.h_dnd[k]);
          REQUIRE(d_baseline.q_dnd_sat[k] == d_test.q_dnd_sat[k]);
          REQUIRE(d_baseline.evp[k] == Approx(d_test.evp[k]).margin(margin));
        }
        REQUIRE(d_baseline.total(d_baseline.totevp) == d_test.total(d_test.totevp));
        REQUIRE(d_baseline.total(d_baseline.totevp) == d_test.total(d_test.jt));
        REQUIRE(d_baseline.total(d_baseline.totevp) == d_test.total(d_test.jd));
        for (Int k = 0; k < d_baseline.total(d_baseline.totevp); ++k) {
          REQUIRE(d_baseline.totevp[k] == Approx(d_test.totevp[k]).margin(margin));
          REQUIRE(d_baseline.jt[k] == d_test.jt[k]);
          REQUIRE(d_baseline.jd[k] == d_test.jd[k]);
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

TEST_CASE("zm_downdraft_properties_bfb", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmDowndraftProperties;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
