#include "catch2/catch.hpp"

#include "rayleigh_friction_unit_tests_common.hpp"
#include "rayleigh_friction_functions.hpp"
#include "rayleigh_friction_test_data.hpp"

#include "share/core/eamxx_types.hpp"

#include <cmath>
#include <limits>

namespace scream {
namespace rayleigh_friction {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestRayleighFrictionInit : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    const Real eps = std::numeric_limits<Real>::epsilon();

    // raytau0 = 0 means no Rayleigh friction, so otau must be identically zero
    {
      RayleighFrictionInitData d(72, 2, 0, 0);
      rayleigh_friction_init(d);
      for (int k=0; k<d.nlevs; ++k) {
        REQUIRE(d.otau[k] == 0);
      }
    }

    // nlevs, rayk0, raykrange, raytau0
    RayleighFrictionInitData test_data[] = {
      RayleighFrictionInitData(72,  2,  0,   5),  // EAM default
      RayleighFrictionInitData(128, 2,  0,   5),
      RayleighFrictionInitData(72,  5,  0,   2),
      RayleighFrictionInitData(72,  10, 3.5, 1),
      RayleighFrictionInitData(16,  1,  2,   0.5),
    };

    for (auto& d : test_data) {
      rayleigh_friction_init(d);

      const Real otau0  = 1 / (86400*d.raytau0);
      const Real krange = (d.raykrange == 0) ? (d.rayk0 - 1) / Real(2) : d.raykrange;

      for (int k=0; k<d.nlevs; ++k) {
        // Profile must be bounded by [0, otau0]
        REQUIRE(d.otau[k] >= 0);
        REQUIRE(d.otau[k] <= otau0);

        // Profile must be monotonically non-increasing away from the model top
        if (k > 0) {
          REQUIRE(d.otau[k] <= d.otau[k-1]);
        }

        // Check against the analytic profile (k is 0-based, rayk0 is 1-based)
        const Real x = (d.rayk0 - (k+1)) / krange;
        const Real expected = otau0 * (1 + std::tanh(x)) / 2;
        REQUIRE(std::abs(d.otau[k] - expected) <= 10*eps*otau0);
      }

      // Decay rate at the center level is exactly half of the top value
      REQUIRE(std::abs(d.otau[d.rayk0-1] - otau0/2) <= 10*eps*otau0);

      // With the default width, x = 2 at the model top
      if (d.raykrange == 0) {
        REQUIRE(std::abs(d.otau[0] - otau0*(1 + std::tanh(Real(2)))/2) <= 10*eps*otau0);
      }
    }
  } // run_property

  void run_bfb()
  {
    // Init has no random inputs, but the baseline file stores the seed
    auto engine = Base::get_engine();
    (void) engine;

    RayleighFrictionInitData baseline_data[] = {
      //                       nlevs, rayk0, raykrange, raytau0
      RayleighFrictionInitData(72,    2,     0,         5),
      RayleighFrictionInitData(128,   2,     0,         5),
      RayleighFrictionInitData(72,    10,    3.5,       1),
      RayleighFrictionInitData(16,    1,     2,         0.5),
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(RayleighFrictionInitData);

    // Create copies of data for use by cxx. Needs to happen before read calls so that
    // inout data is in original state
    RayleighFrictionInitData cxx_data[] = {
      RayleighFrictionInitData(baseline_data[0]),
      RayleighFrictionInitData(baseline_data[1]),
      RayleighFrictionInitData(baseline_data[2]),
      RayleighFrictionInitData(baseline_data[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_ifile);
      }
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      rayleigh_friction_init(d);
    }

    // Verify BFB results
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (int r = 0; r<num_runs; ++r) {
        RayleighFrictionInitData& d_baseline = baseline_data[r];
        RayleighFrictionInitData& d_cxx = cxx_data[r];
        REQUIRE(d_baseline.total(d_baseline.otau) == d_cxx.total(d_cxx.otau));
        for (int k=0; k<d_baseline.total(d_baseline.otau); ++k) {
          REQUIRE(d_baseline.otau[k] == d_cxx.otau[k]);
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

TEST_CASE("rayleigh_friction_init_property", "rayleigh_friction")
{
  using TestStruct = scream::rayleigh_friction::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestRayleighFrictionInit;

  TestStruct t;
  t.run_property();
}

TEST_CASE("rayleigh_friction_init_bfb", "rayleigh_friction")
{
  using TestStruct = scream::rayleigh_friction::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestRayleighFrictionInit;

  TestStruct t;
  t.run_bfb();
}

} // empty namespace
