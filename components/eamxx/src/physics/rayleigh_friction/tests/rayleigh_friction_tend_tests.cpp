#include "catch2/catch.hpp"

#include "rayleigh_friction_unit_tests_common.hpp"
#include "rayleigh_friction_functions.hpp"
#include "rayleigh_friction_test_data.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/core/eamxx_setup_random_test.hpp"
#include "share/core/eamxx_types.hpp"

#include <cmath>
#include <limits>

namespace scream {
namespace rayleigh_friction {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestRayleighFrictionTend : public UnitWrap::UnitTest<D>::Base {

  // Randomize with physically reasonable ranges
  template <typename Engine>
  static void randomize(Engine& engine, RayleighFrictionTendData& d, const Real otau_max)
  {
    d.randomize(engine, { {d.otau,   {0, otau_max}},
                          {d.u_wind, {-100, 100}},
                          {d.v_wind, {-100, 100}},
                          {d.t_mid,  {150, 320}} });
  }

  void run_property()
  {
    using C = physics::Constants<Real>;
    const Real cpair = C::Cpair.value;
    const Real eps   = std::numeric_limits<Real>::epsilon();

    auto engine = Base::get_engine();

    // Zero decay rate must leave the state untouched
    {
      RayleighFrictionTendData d(5, 72, 300);
      randomize(engine, d, 0);
      for (int k=0; k<d.nlevs; ++k) d.otau[k] = 0;
      RayleighFrictionTendData d0(d);

      rayleigh_friction_tend(d);

      for (int n=0; n<d.total(d.u_wind); ++n) {
        REQUIRE(d.u_wind[n] == d0.u_wind[n]);
        REQUIRE(d.v_wind[n] == d0.v_wind[n]);
        REQUIRE(d.t_mid[n]  == d0.t_mid[n]);
      }
    }

    // ncols, nlevs, dt, max decay rate [1/s]
    // The last case is a very strong damping limit (otau*dt >> 1)
    struct Case { int ncols; int nlevs; Real dt; Real otau_max; };
    const Case cases[] = {
      {10, 72,  1800, 1/(86400*Real(5))},
      {7,  128, 300,  1/(86400*Real(0.1))},
      {3,  16,  100,  1},
    };

    for (const auto& c : cases) {
      RayleighFrictionTendData d(c.ncols, c.nlevs, c.dt);
      randomize(engine, d, c.otau_max);
      RayleighFrictionTendData d0(d);

      rayleigh_friction_tend(d);

      for (int i=0; i<d.ncols; ++i) {
        for (int k=0; k<d.nlevs; ++k) {
          const int n = i*d.nlevs + k;
          const Real u0 = d0.u_wind[n], v0 = d0.v_wind[n], t0 = d0.t_mid[n];
          const Real u1 = d.u_wind[n],  v1 = d.v_wind[n],  t1 = d.t_mid[n];

          // Euler backward damping: u_new = u_old / (1 + otau*dt)
          const Real c2 = 1 / (1 + d.otau[k]*d.dt);
          REQUIRE(std::abs(u1 - c2*u0) <= 10*eps*std::abs(u0));
          REQUIRE(std::abs(v1 - c2*v0) <= 10*eps*std::abs(v0));

          // Winds are damped without changing direction
          REQUIRE(std::abs(u1) <= std::abs(u0));
          REQUIRE(std::abs(v1) <= std::abs(v0));
          REQUIRE(u1*u0 >= 0);
          REQUIRE(v1*v0 >= 0);

          // Friction heats the air
          REQUIRE(t1 >= t0);

          // Total energy (kinetic + enthalpy) is conserved
          const Real ke0 = (u0*u0 + v0*v0) / 2;
          const Real ke1 = (u1*u1 + v1*v1) / 2;
          const Real energy_err = (ke1 + cpair*t1) - (ke0 + cpair*t0);
          REQUIRE(std::abs(energy_err) <= 10*eps*(ke0 + cpair*t0));
        }
      }
    }
  } // run_property

  void run_bfb()
  {
    auto engine = Base::get_engine();

    RayleighFrictionTendData baseline_data[] = {
      //                       ncols, nlevs, dt
      RayleighFrictionTendData(12,    72,    1800),
      RayleighFrictionTendData(8,     128,   300),
      RayleighFrictionTendData(7,     16,    100),
      RayleighFrictionTendData(2,     7,     600),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(RayleighFrictionTendData);

    // Generate random input data
    for (auto& d : baseline_data) {
      randomize(engine, d, 1/(86400*Real(0.1)));
    }

    // Create copies of data for use by cxx. Needs to happen before read calls so that
    // inout data is in original state
    RayleighFrictionTendData cxx_data[] = {
      RayleighFrictionTendData(baseline_data[0]),
      RayleighFrictionTendData(baseline_data[1]),
      RayleighFrictionTendData(baseline_data[2]),
      RayleighFrictionTendData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      rayleigh_friction_tend(d);
    }

    // Verify BFB results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (int r = 0; r<num_runs; ++r) {
        RayleighFrictionTendData& d_baseline = baseline_data[r];
        RayleighFrictionTendData& d_cxx = cxx_data[r];
        REQUIRE(d_baseline.total(d_baseline.u_wind) == d_cxx.total(d_cxx.u_wind));
        for (int n=0; n<d_baseline.total(d_baseline.u_wind); ++n) {
          REQUIRE(d_baseline.u_wind[n] == d_cxx.u_wind[n]);
          REQUIRE(d_baseline.v_wind[n] == d_cxx.v_wind[n]);
          REQUIRE(d_baseline.t_mid[n]  == d_cxx.t_mid[n]);
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
} // namespace rayleigh_friction
} // namespace scream

namespace {

TEST_CASE("rayleigh_friction_tend_property", "rayleigh_friction")
{
  using TestStruct = scream::rayleigh_friction::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestRayleighFrictionTend;

  TestStruct t;
  t.run_property();
}

TEST_CASE("rayleigh_friction_tend_bfb", "rayleigh_friction")
{
  using TestStruct = scream::rayleigh_friction::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestRayleighFrictionTend;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
