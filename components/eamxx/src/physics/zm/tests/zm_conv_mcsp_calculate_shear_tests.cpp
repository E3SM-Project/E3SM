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
      zm_conv_mcsp_calculate_shear(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ZmConvMcspCalculateShearData& d_baseline = baseline_data[i];
        ZmConvMcspCalculateShearData& d_test = test_data[i];
        REQUIRE(d_baseline.total(d_baseline.shear_u) == d_test.total(d_test.shear_u));
        for (Int k = 0; k < d_baseline.total(d_baseline.shear_u); ++k) {
          REQUIRE(d_baseline.shear_u[k] == d_test.shear_u[k]);
        }
        REQUIRE(d_baseline.total(d_baseline.shear_v) == d_test.total(d_test.shear_v));
        for (Int k = 0; k < d_baseline.total(d_baseline.shear_v); ++k) {
          REQUIRE(d_baseline.shear_v[k] == d_test.shear_v[k]);
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

TEST_CASE("zm_shear_scaling", "[zm]")
{
  using ZMF = scream::zm::Functions<scream::Real, scream::DefaultDevice>;
  using scream::Real;

  ZMF::ZmRuntimeOpt opt;
  opt.mcsp_use_full_shear    = true;
  opt.use_shear_scaling_mcsp = true;
  opt.use_shear_scaling_evap = true;
  opt.mcsp_shear_min         = 4;
  opt.mcsp_shear_optimal_min = 12;
  opt.mcsp_shear_optimal_max = 30;
  opt.mcsp_shear_max         = 55;
  opt.evap_shear_min         = 5;
  opt.evap_shear_max         = 25;
  opt.evap_efficiency_min    = 1.0e-6;
  opt.evap_efficiency_max    = 3.0e-6;

  // shear metric
  REQUIRE(opt.calc_shear_metric(3, 4) == Approx(5));
  opt.mcsp_use_full_shear = false;
  REQUIRE(opt.calc_shear_metric(-3, 4) == -3);
  opt.mcsp_use_full_shear = true;
  // undefined shear (surface above the storm reference level) maps to the fill value
  REQUIRE(opt.calc_shear_metric(ZMF::ZMC::MCSP_shear_undefined, ZMF::ZMC::MCSP_shear_undefined)
          == scream::constants::fill_value<Real>);
  REQUIRE(opt.calc_mcsp_shear_scaling(scream::constants::fill_value<Real>) == 0);
  REQUIRE(opt.calc_evap_efficiency(scream::constants::fill_value<Real>) == Approx(3.0e-6));

  // trapezoid shape - zero outside, one on the plateau, linear ramps
  REQUIRE(opt.calc_mcsp_shear_scaling(0)    == 0);
  REQUIRE(opt.calc_mcsp_shear_scaling(4)    == 0);
  REQUIRE(opt.calc_mcsp_shear_scaling(8)    == Approx(0.5));
  REQUIRE(opt.calc_mcsp_shear_scaling(12)   == 1);
  REQUIRE(opt.calc_mcsp_shear_scaling(20)   == 1);
  REQUIRE(opt.calc_mcsp_shear_scaling(30)   == 1);
  REQUIRE(opt.calc_mcsp_shear_scaling(42.5) == Approx(0.5));
  REQUIRE(opt.calc_mcsp_shear_scaling(55)   == 0);
  REQUIRE(opt.calc_mcsp_shear_scaling(100)  == 0);
  REQUIRE(opt.calc_mcsp_shear_scaling(-999) == 0); // sentinel for surface above the storm reference level
  REQUIRE(opt.calc_mcsp_shear_scaling(-8)   == Approx(0.5)); // signed zonal shear uses the magnitude

  // evaporation scaling is a monotonic ramp - zero below, one above, linear between
  REQUIRE(opt.calc_evap_shear_scaling(0)    == 0);
  REQUIRE(opt.calc_evap_shear_scaling(5)    == 0);
  REQUIRE(opt.calc_evap_shear_scaling(15)   == Approx(0.5));
  REQUIRE(opt.calc_evap_shear_scaling(25)   == 1);
  REQUIRE(opt.calc_evap_shear_scaling(100)  == 1);
  REQUIRE(opt.calc_evap_shear_scaling(-15)  == Approx(0.5)); // signed zonal shear uses the magnitude

  // evaporation efficiency interpolates between min and max
  REQUIRE(opt.calc_evap_efficiency(0)   == Approx(1.0e-6));
  REQUIRE(opt.calc_evap_efficiency(15)  == Approx(2.0e-6));
  REQUIRE(opt.calc_evap_efficiency(25)  == Approx(3.0e-6));
  REQUIRE(opt.calc_evap_efficiency(100) == Approx(3.0e-6));

  // without shear scaling the max efficiency is always used (legacy ke)
  opt.use_shear_scaling_evap = false;
  REQUIRE(opt.calc_evap_efficiency(0)  == 3.0e-6);
  REQUIRE(opt.calc_evap_efficiency(20) == 3.0e-6);
}

} // empty namespace
