#include "catch2/catch.hpp"

#include "share/core/eamxx_types.hpp"
#include "physics/zm/zm_functions.hpp"
#include "physics/zm/tests/infra/zm_test_data.hpp"

#include "zm_unit_tests_common.hpp"

#include <cmath>
#include <limits>

namespace scream {
namespace zm {
namespace unit_test {

/*
 * Tests for the optional humidity and temperature scaling of the ZM
 * autoconversion coefficient (c0). These are not baseline tests: the legacy
 * (scaling disabled) path is covered by the zm_cloud_properties_bfb and
 * zm_conv_main_bfb baselines, while the tests here exercise the enabled paths.
 */
template <typename D>
struct UnitWrap::UnitTest<D>::TestZmC0Scaling : public UnitWrap::UnitTest<D>::Base {

  using ZmRuntimeOpt = typename Functions::ZmRuntimeOpt;
  using ZMC          = typename Functions::ZMC;

  static ZmRuntimeOpt make_opts(const bool rh_on, const bool t_on)
  {
    ZmRuntimeOpt opts;
    opts.c0_lnd              = ZMC::c0;
    opts.c0_ocn              = ZMC::c0;
    opts.use_rh_scaling_c0   = rh_on;
    opts.use_temp_scaling_c0 = t_on;
    opts.c0_rh_min           = ZMC::c0_rh_min_default;
    opts.c0_sf_lo            = ZMC::c0_sf_lo_default;
    opts.c0_sf_hi            = ZMC::c0_sf_hi_default;
    opts.c0_sf_p_bot         = ZMC::c0_sf_p_bot_default;
    opts.c0_sf_p_top         = ZMC::c0_sf_p_top_default;
    opts.c0_t_lo             = ZMC::c0_t_lo_default;
    opts.c0_t_hi             = ZMC::c0_t_hi_default;
    opts.c0_t_min_frac       = ZMC::c0_t_min_frac_default;
    return opts;
  }

  //----------------------------------------------------------------------------
  // factor functions: disabled paths are exact, ramps hit their end points
  void run_factors()
  {
    const Real tol = std::numeric_limits<Real>::epsilon() * 10;

    // disabled: both factors exactly one and c0 exactly the land/ocean blend
    {
      const auto opts = make_opts(false, false);
      for (Real sf : {Real(-1), Real(0), Real(0.3), Real(0.7), Real(1), Real(2)}) {
        REQUIRE(opts.calc_c0_rh_factor(sf) == 1);
      }
      for (Real t : {Real(180), Real(253.15), Real(260), Real(273.15), Real(300)}) {
        REQUIRE(opts.calc_c0_temp_factor(t) == 1);
      }
      const Real c0_base = opts.calc_c0_base(0.3);
      REQUIRE(opts.calc_c0(c0_base, Real(0.5), Real(200)) == c0_base);
    }

    // humidity ramp
    {
      const auto opts = make_opts(true, false);
      const Real lo = opts.c0_sf_lo, hi = opts.c0_sf_hi, fmin = opts.c0_rh_min;
      REQUIRE(opts.calc_c0_rh_factor(lo - 0.2) == fmin);
      REQUIRE(opts.calc_c0_rh_factor(lo)       == Approx(fmin).margin(tol));
      REQUIRE(opts.calc_c0_rh_factor(hi)       == Approx(1).margin(tol));
      REQUIRE(opts.calc_c0_rh_factor(hi + 0.1) == 1);
      REQUIRE(opts.calc_c0_rh_factor(Real(0.5)*(lo+hi)) == Approx(Real(0.5)*(1+fmin)).margin(tol));
      // monotone non-decreasing and bounded by [c0_rh_min,1]
      Real prev = 0;
      for (Int n = 0; n <= 100; ++n) {
        const Real f = opts.calc_c0_rh_factor(Real(n)/100);
        REQUIRE(f >= fmin);
        REQUIRE(f <= 1);
        REQUIRE(f >= prev);
        prev = f;
      }
      // temperature has no effect when only the humidity scaling is enabled
      REQUIRE(opts.calc_c0_temp_factor(Real(200)) == 1);
    }

    // temperature ramp
    {
      const auto opts = make_opts(false, true);
      const Real lo = opts.c0_t_lo, hi = opts.c0_t_hi, fmin = opts.c0_t_min_frac;
      REQUIRE(opts.calc_c0_temp_factor(lo - 20) == fmin);
      REQUIRE(opts.calc_c0_temp_factor(lo)      == fmin);
      REQUIRE(opts.calc_c0_temp_factor(hi)      == 1);
      REQUIRE(opts.calc_c0_temp_factor(hi + 20) == 1);
      REQUIRE(opts.calc_c0_temp_factor(Real(0.5)*(lo+hi)) == Approx(Real(0.5)*(1+fmin)).margin(tol));
      // the humidity factor is ignored when only the temperature scaling is enabled
      const Real c0_base = opts.calc_c0_base(0);
      REQUIRE(opts.calc_c0(c0_base, Real(0.3), hi + 1) == c0_base);
    }

    // both enabled: c0 is the product of the base value and both factors
    {
      const auto opts = make_opts(true, true);
      const Real c0_base = opts.calc_c0_base(0.6);
      const Real rh_fac  = opts.calc_c0_rh_factor(Real(0.6));
      for (Real t : {Real(240), Real(263.15), Real(290)}) {
        const Real expect = c0_base * rh_fac * opts.calc_c0_temp_factor(t);
        REQUIRE(opts.calc_c0(c0_base, rh_fac, t) == Approx(expect).margin(tol*c0_base));
      }
    }
  } // run_factors

  //----------------------------------------------------------------------------
  // Fill d with a deterministic, conditionally unstable tropical sounding so
  // that zm_cloud_properties produces an active updraft with liquid water and
  // rain production. Uniform layers of thickness dz from the surface to z_top,
  // a 6.5 K/km lapse rate up to an isothermal stratosphere, and relative
  // humidity decreasing from 90% near the surface to 30% aloft. Columns differ
  // in surface temperature and land fraction. k=0 is the model top.
  static void set_convecting_sounding(ZmCloudPropertiesData& d)
  {
    const Real grav  = 9.80616, cpair = 1004.64, rdair = 287.042;
    const Real z_top = 20000, dz = z_top / d.pver;
    for (Int i = 0; i < d.pcols; ++i) {
      const Real t_sfc = 300 + Real(0.5)*i;
      d.landfrac[i] = d.pcols > 1 ? Real(i) / (d.pcols - 1) : Real(0);
      d.tpert_g[i]  = 0.5;
      for (Int k = 0; k < d.pverp; ++k) {
        d.z_int[i*d.pverp + k] = (d.pverp - 1 - k) * dz;
      }
      // integrate the hypsometric equation upward from the surface
      Real p_below = 1000; // [mb]
      Real z_below = 0;
      for (Int k = d.pver - 1; k >= 0; --k) {
        const Int  n  = i*d.pver + k;
        const Real z  = Real(0.5) * (d.z_int[i*d.pverp + k] + d.z_int[i*d.pverp + k + 1]);
        const Real t  = std::max(t_sfc - Real(6.5e-3)*z, t_sfc - Real(6.5e-3)*15000);
        const Real p  = p_below * std::exp(-grav * (z - z_below) / (rdair * t));
        const Real rh = z < 2000 ? Real(0.9) : (z < 10000 ? Real(0.9) - Real(0.6)*(z - 2000)/8000 : Real(0.3));
        // saturation specific humidity from Bolton (1980)
        const Real es = Real(6.112) * std::exp(Real(17.67)*(t - Real(273.15)) / (t - Real(29.65)));
        const Real qs = Real(0.622) * es / std::max(p - Real(0.378)*es, es);
        d.z_mid[n] = z;
        d.t_mid[n] = t;
        d.p_mid[n] = p;
        d.q_mid[n] = rh * qs;
        d.s_mid[n] = t + grav / cpair * z;
        p_below = p;
        z_below = z;
      }
      // interface dry static energy: mean of the adjacent layers
      for (Int k = 0; k < d.pver; ++k) {
        const Int n = i*d.pver + k;
        d.s_int[n] = k == 0 ? d.s_mid[n] : Real(0.5) * (d.s_mid[n-1] + d.s_mid[n]);
      }
      // updraft base in the second layer above the surface, launch level near 13 km
      d.jb[i]  = d.pver - 2;
      d.lel[i] = Int(d.pver * (z_top - 13000) / z_top);
    }
  }

  //----------------------------------------------------------------------------
  // zm_cloud_properties with the scaling enabled
  void run_cloud_properties()
  {
    auto engine = Base::get_engine();

    // The deterministic sounding must exercise the rain production and liquid
    // budget checks with non-neutral scaling, the random data only has to pass
    // them where it happens to convect
    struct DataSet {
      ZmCloudPropertiesData d;
      bool require_active;
    };
    std::vector<DataSet> base_data = {
      //                     pcols, ncol, pver, pverp, msg, limcnv
      {ZmCloudPropertiesData(    4,    4,   40,    41,   2,   3), true},
      {ZmCloudPropertiesData(    4,    4,   72,    73,  10,  40), false},
      {ZmCloudPropertiesData(    4,    4,  128,   129,  30,  60), false},
    };
    for (auto& ds : base_data) {
      if (ds.require_active) {
        set_convecting_sounding(ds.d);
      } else {
        ds.d.randomize(engine);
      }
    }

    struct Case {
      ZmC0ScalingTestOpts c0;
      bool rh_on, t_on;
    };

    const Real rh_fac = 0.4;
    const Real t_frac = 0.3;
    std::vector<Case> cases(3);
    cases[0].c0.use_rh_scaling_c0   = true;  cases[0].c0.c0_rh_fac = rh_fac;  // humidity only
    cases[1].c0.use_temp_scaling_c0 = true;  cases[1].c0.c0_t_min_frac = t_frac; // temperature only
    cases[2].c0.use_rh_scaling_c0   = true;  cases[2].c0.c0_rh_fac = rh_fac;  // both
    cases[2].c0.use_temp_scaling_c0 = true;  cases[2].c0.c0_t_min_frac = t_frac;
    for (auto& c : cases) {
      c.rh_on = c.c0.use_rh_scaling_c0;
      c.t_on  = c.c0.use_temp_scaling_c0;
    }

    for (const auto& ds : base_data) {
      const auto& d0 = ds.d;
      const Int pcols = d0.pcols, pver = d0.pver, pverp = d0.pverp;

      // reference run with the scaling disabled
      ZmCloudPropertiesData d_off(d0);
      zm_cloud_properties(d_off);

      // enabled with neutral factors must reproduce the disabled path exactly
      {
        ZmCloudPropertiesData d_neutral(d0);
        ZmC0ScalingTestOpts neutral;
        neutral.use_rh_scaling_c0   = true;
        neutral.use_temp_scaling_c0 = true;
        neutral.c0_rh_fac           = 1;
        neutral.c0_t_min_frac       = 1;
        zm_cloud_properties(d_neutral, neutral);
        for (Int n = 0; n < d_off.total(d_off.ql); ++n) {
          REQUIRE(d_neutral.ql[n]   == d_off.ql[n]);
          REQUIRE(d_neutral.rprd[n] == d_off.rprd[n]);
          REQUIRE(d_neutral.evp[n]  == d_off.evp[n]);
        }
      }

      for (const auto& c : cases) {
        ZmCloudPropertiesData d(d0);
        zm_cloud_properties(d, c.c0);

        ZmRuntimeOpt opts = make_opts(c.rh_on, c.t_on);
        if (c.c0.c0_t_min_frac >= 0) opts.c0_t_min_frac = c.c0.c0_t_min_frac;

        // tolerance for host recomputation of device results
        const Real rtol = ekat::is_single_precision<Real>::value ? 1e-4 : 1e-10;

        Int n_active = 0, n_scaled = 0;
        for (Int i = 0; i < pcols; ++i) {
          // updraft properties upstream of the liquid budget must not depend on c0
          REQUIRE(d.jt[i] == d_off.jt[i]);
          const Real c0_base = opts.calc_c0_base(d.landfrac[i]);
          for (Int k = 0; k < pver; ++k) {
            const Int n = i*pver + k;
            REQUIRE(d.mflx_up[n] == d_off.mflx_up[n]);
            REQUIRE(d.cu[n]      == d_off.cu[n]);
            REQUIRE(d.ql[n] >= 0);

            // zm_cloud_properties returns rprd net of downdraft evaporation,
            // so add evp back to recover the gross autoconversion rate
            const Real rprd_auto = d.rprd[n] + d.evp[n];
            if (rprd_auto > 0 && d.mflx_up[n] > 0 && k+1 < pver) {
              ++n_active;
              const Real c0_k = opts.calc_c0(c0_base, c.c0.c0_rh_fac, d.t_mid[n]);
              if (c0_k < c0_base * (1 - rtol)) ++n_scaled;
              // rain production uses the scaled c0
              const Real rprd_expect = c0_k * d.mflx_up[n] * d.ql[n];
              REQUIRE(std::abs(rprd_auto - rprd_expect) <= rtol * (std::abs(rprd_expect) + std::abs(d.evp[n])));
              // and the liquid budget closes with the same c0:
              // flux in + condensation = flux out + detrainment + rain production
              const Real dz   = d.z_int[i*pverp + k] - d.z_int[i*pverp + k + 1];
              const Real t_in = d.mflx_up[n+1]*d.ql[n+1];
              const Real t_cu = dz*d.cu[n];
              const Real t_out = d.mflx_up[n]*d.ql[n];
              const Real t_det = dz*d.detr_up[n]*d.ql[n+1];
              const Real t_rp  = dz*rprd_auto;
              const Real scale = std::abs(t_in) + std::abs(t_cu) + std::abs(t_out)
                               + std::abs(t_det) + std::abs(t_rp) + dz*std::abs(d.evp[n]);
              REQUIRE(std::abs((t_in + t_cu) - (t_out + t_det + t_rp)) <= rtol * scale);
            }
          }
        }
        if (ds.require_active) {
          // fail rather than pass silently if the enabled path was not exercised
          REQUIRE(n_active > 0);
          REQUIRE(n_scaled > 0);
        }
      }
    }
  } // run_cloud_properties

};

} // namespace unit_test
} // namespace zm
} // namespace scream

namespace {

TEST_CASE("zm_c0_scaling_factors", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmC0Scaling;

  TestStruct t;
  t.run_factors();
}

TEST_CASE("zm_cloud_properties_c0_scaling", "[zm]")
{
  using TestStruct = scream::zm::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestZmC0Scaling;

  TestStruct t;
  t.run_cloud_properties();
}

} // empty namespace
