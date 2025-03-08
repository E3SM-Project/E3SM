#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/eamxx_types.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestShocTke : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Real mintke = scream::shoc::Constants<Real>::mintke;
    static constexpr Real maxtke = scream::shoc::Constants<Real>::maxtke;
    static constexpr Real maxiso = scream::shoc::Constants<Real>::maxiso;
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the subroutine shoc_tke.  This is the top
    //  level function for the TKE module in SHOC.  Thus since
    //  all associated functions here contain property tests, the
    //  purpose here is mostly to expose whether errors have been
    //  made in the input/output fields in these routines

    // TEST ONE
    // For this test we will perform a TKE growth test.
    //  Will provide inputs that should always encourage TKE to
    //  grow at each level.

    // Define height thickness on nlevi grid [m]
    //   NOTE: First indicee is zero because it is never used
    //   Do a stretched grid

    // Timestep [s]
    Real dtime = 300;
    // Buoyancy flux [K m/s]
    Real wthv_sec[nlev] = {0.05, 0.04, 0.03, 0.02, 0.03};
    // Length Scale [m]
    Real shoc_mix[nlev] = {1000, 750, 500, 400, 300};
    // Define zonal wind on nlev grid [m/s]
    Real u_wind[nlev] = {2, 1, 0, -1, -2};
    // Define meridional wind on nlev grid [m/s]
    Real v_wind[nlev] = {1, 2, 3, 4, 5};
    // Define surface temperature [K] (value irrelevant, just make sure it's physical)
    Real tabs[nlev] ={300, 300, 300, 300, 300};
    // Define thickness on the interface grid [m]
    Real dz_zi[nlevi] = {0, 100, 100, 100, 100, 50};
    // Define thickness on the thermo grid [m]
    Real dz_zt[nlev] = {100, 100, 100, 100, 100};
    // Define the midpoint height grid [m]
    Real zt_grid[nlev] = {450, 350, 250, 150, 50};
    // Define the interface height grid [m]
    Real zi_grid[nlevi] = {500, 400, 300, 200, 100, 0};
    // Define PBL height [m]
    Real pblh = 1000;
    // Define pressure [Pa]
    Real pres[nlev] = {85000, 87500, 90000, 95000, 100000};

    // Define TKE initial [m2/s2]
    Real tke_init[nlev] = {0.01, 0.01, 0.01, 0.01, 0.01};
    // Define initial eddy diffusivities
    Real tkh[nlev] = {1, 1, 1, 1, 1};
    Real tk[nlev];

    // Set Tk equal to tkh
    for(Int n = 0; n < nlev; ++n){
      tk[n] = tkh[n];
    }

    // Initialize data structure for bridging to F90
    ShocTkeData SDS(shcol, nlev, nlevi, dtime);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.nlevi - SDS.nlev == 1);
    REQUIRE(SDS.shcol > 0);

    SDS.dtime = dtime;
    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      SDS.pblh[s] = pblh;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.wthv_sec[offset] = wthv_sec[n];
        SDS.shoc_mix[offset] = shoc_mix[n];
        SDS.u_wind[offset] = u_wind[n];
        SDS.v_wind[offset] = v_wind[n];
        SDS.dz_zt[offset] = dz_zt[n];
        SDS.zt_grid[offset] = zt_grid[n];
        SDS.pres[offset] = pres[n];
        SDS.tke[offset] = tke_init[n];
        SDS.tkh[offset] = tkh[n];
        SDS.tk[offset] = tk[n];
        SDS.tabs[offset] = tabs[n];
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset   = n + s * nlevi;
        SDS.dz_zi[offset] = dz_zi[n];
        SDS.zi_grid[offset] = zi_grid[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      // nlevi loop
      for (Int n = 0; n < nlevi; ++n){
        const auto offset = n + s * nlevi;
        // Make sure top level dz_zi value is zero
        if (n == 0){
          REQUIRE(SDS.dz_zi[offset] == 0);
        }
        // Otherwise, should be greater than zero
        else{
          REQUIRE(SDS.dz_zi[offset] > 0);
        }
        // Check that zi increases updward
        if (n < nlevi-1){
          REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0);
        }
      }
      // nlev loop
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Check that zt increases upward
        if (n < nlev-1){
          REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
        }
        REQUIRE(SDS.dz_zt[offset] > 0);
        REQUIRE(SDS.shoc_mix[offset] > 0);
        REQUIRE(SDS.tke[offset] >= mintke);
        REQUIRE(SDS.tkh[offset] > 0);
        REQUIRE(SDS.tk[offset] > 0);
        REQUIRE(SDS.pres[offset] > 0);
      }
    }

    // Call the C++ implementation
    shoc_tke(SDS);

    // Check test
    // Make sure that TKE has increased everwhere relative
    //   to initial value.  Also make sure that outputs fall
    //   within some reasonable bounds.

    // Make array to save the result of TKE
    Real tke_test1[nlev*shcol];

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.tke[offset] > tke_init[n]);
        REQUIRE(SDS.tke[offset] >= mintke);
        REQUIRE(SDS.tke[offset] <= maxtke);
        REQUIRE(SDS.tkh[offset] > 0);
        REQUIRE(SDS.tk[offset] > 0);
        REQUIRE(SDS.isotropy[offset] >= 0);
        REQUIRE(SDS.isotropy[offset] <= maxiso);
        tke_test1[offset] = SDS.tke[offset];
      }
    }

    // TEST TWO
    // Decay test.  Now starting with the TKE from TEST ONE in
    // its spun up state, feed inputs that should always make
    // TKE decay.

    // New inputs below, all others are recycled

    // characteristics of new input are negative buoyancy flux,
    //  small length scale, and a uniform wind profile.

    // Buoyancy flux [K m/s]
    Real wthv_sec_decay[nlev] = {-0.05, -0.04, -0.03, -0.02, -0.03};
    // Length Scale [m]
    Real shoc_mix_decay[nlev] = {100, 75, 50, 40, 30};
    // Define zonal wind on nlev grid [m/s]
    Real u_wind_decay[nlev] = {1, 1, 1, 1, 1};
    // Define meridional wind on nlev grid [m/s]
    Real v_wind_decay[nlev] = {-2, -2, -2, -2, -2};

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.wthv_sec[offset] = wthv_sec_decay[n];
        SDS.shoc_mix[offset] = shoc_mix_decay[n];
        SDS.u_wind[offset] = u_wind_decay[n];
        SDS.v_wind[offset] = v_wind_decay[n];

      }
    }

    for(Int s = 0; s < shcol; ++s) {
      // nlev loop
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // be sure wind has no gradient
        if (n < nlev-1){
          REQUIRE(SDS.u_wind[offset + 1] - SDS.u_wind[offset] == 0);
          REQUIRE(SDS.v_wind[offset + 1] - SDS.v_wind[offset] == 0);
        }
        // Be sure that buoyancy flux is less than zero
        REQUIRE(SDS.wthv_sec[offset] < 0);
        // Be sure length scale is reasonably small
        REQUIRE(SDS.shoc_mix[offset] <= 100);
      }
    }

    // Call the C++ implementation
    shoc_tke(SDS);

    // Check the result

    // Verify ALL outputs are reasonable and that TKE has decayed
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.tke[offset] < tke_test1[offset]);
        REQUIRE(SDS.tke[offset] >= mintke);
        REQUIRE(SDS.tke[offset] <= maxtke);
        REQUIRE(SDS.tkh[offset] > 0);
        REQUIRE(SDS.tk[offset] > 0);
        REQUIRE(SDS.isotropy[offset] >= 0);
        REQUIRE(SDS.isotropy[offset] <= maxiso);
      }
    }
  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    ShocTkeData baseline_data[] = {
      ShocTkeData(10, 71, 72, 300),
      ShocTkeData(10, 12, 13, 100),
      ShocTkeData(7,  16, 17, 50),
      ShocTkeData(2, 7, 8, 5),
    };

    // Generate random input data
    // Alternatively, you can use the baseline_data construtors/initializer lists to hardcode data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    ShocTkeData cxx_data[] = {
      ShocTkeData(baseline_data[0]),
      ShocTkeData(baseline_data[1]),
      ShocTkeData(baseline_data[2]),
      ShocTkeData(baseline_data[3]),
    };

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      shoc_tke(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ShocTkeData);
      for (Int i = 0; i < num_runs; ++i) {
        ShocTkeData& d_baseline = baseline_data[i];
        ShocTkeData& d_cxx = cxx_data[i];
        REQUIRE(d_baseline.total(d_baseline.tke) == d_cxx.total(d_cxx.tke));
        REQUIRE(d_baseline.total(d_baseline.tke) == d_cxx.total(d_cxx.tk));
        REQUIRE(d_baseline.total(d_baseline.tke) == d_cxx.total(d_cxx.tkh));
        REQUIRE(d_baseline.total(d_baseline.tke) == d_cxx.total(d_cxx.isotropy));
        for (Int k = 0; k < d_baseline.total(d_baseline.tke); ++k) {
          REQUIRE(d_baseline.tke[k] == d_cxx.tke[k]);
          REQUIRE(d_baseline.tk[k] == d_cxx.tk[k]);
          REQUIRE(d_baseline.tkh[k] == d_cxx.tkh[k]);
          REQUIRE(d_baseline.isotropy[k] == d_cxx.isotropy[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (auto& d : cxx_data) {
        d.write(Base::m_fid);
      }
    }
  } // run_bfb
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_tke_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocTke;

  TestStruct().run_property();
}

TEST_CASE("shoc_tke_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocTke;

  TestStruct().run_bfb();

}

} // namespace
