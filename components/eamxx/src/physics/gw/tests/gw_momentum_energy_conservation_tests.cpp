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
struct UnitWrap::UnitTest<D>::TestMomentumEnergyConservation : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    // Set up init data
    auto init_data = get_common_init_data(engine);

    // Set up inputs
    MomentumEnergyConservationData baseline_data[] = {
      //                          ncol, dt
      MomentumEnergyConservationData(2, .4, init_data[0]),
      MomentumEnergyConservationData(3, .8, init_data[1]),
      MomentumEnergyConservationData(4, 1.4, init_data[2]),
      MomentumEnergyConservationData(5, 2.4, init_data[3]),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(MomentumEnergyConservationData);

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (Int i = 0; i < num_runs; ++i) {
      auto& d = baseline_data[i];
      d.randomize(engine, { {d.tend_level, {init_data[i].ktop+1, init_data[i].kbotbg-1}} });
    }

    // Create copies of data for use by test. Needs to happen before read calls so that
    // inout data is in original state
    MomentumEnergyConservationData test_data[] = {
      MomentumEnergyConservationData(baseline_data[0]),
      MomentumEnergyConservationData(baseline_data[1]),
      MomentumEnergyConservationData(baseline_data[2]),
      MomentumEnergyConservationData(baseline_data[3]),
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
        momentum_energy_conservation_f(d);
      }
      else {
        momentum_energy_conservation(d);
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
        MomentumEnergyConservationData& d_baseline = baseline_data[i];
        MomentumEnergyConservationData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.dudt) == d_test.total(d_test.dudt));
        REQUIRE(d_baseline.total(d_baseline.dudt) == d_test.total(d_test.dsdt));
        REQUIRE(d_baseline.total(d_baseline.dudt) == d_test.total(d_test.dsdt));
        for (Int k = 0; k < d_baseline.total(d_baseline.dudt); ++k) {
          REQUIRE(d_baseline.dudt[k] == Approx(d_test.dudt[k]).margin(margin));
          REQUIRE(d_baseline.dvdt[k] == Approx(d_test.dvdt[k]).margin(margin));
          REQUIRE(d_baseline.dsdt[k] == Approx(d_test.dsdt[k]).margin(margin));
        }
        REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.utgw));
        REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.vtgw));
        REQUIRE(d_baseline.total(d_baseline.utgw) == d_test.total(d_test.ttgw));
        for (Int k = 0; k < d_baseline.total(d_baseline.utgw); ++k) {
          REQUIRE(d_baseline.utgw[k] == Approx(d_test.utgw[k]).margin(margin));
          REQUIRE(d_baseline.vtgw[k] == Approx(d_test.vtgw[k]).margin(margin));
          REQUIRE(d_baseline.ttgw[k] == Approx(d_test.ttgw[k]).margin(margin));
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

TEST_CASE("momentum_energy_conservation_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestMomentumEnergyConservation;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
