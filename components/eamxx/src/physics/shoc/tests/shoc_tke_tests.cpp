#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_setup_random_test.hpp"

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
struct UnitWrap::UnitTest<D>::TestShocTke {

  static void run_property()
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
    // Define Obukov length [m]
    Real obklen = -1;
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
      SDS.obklen[s] = obklen;
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

    // Call the fortran implementation
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

    // Call the fortran implementation
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

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ShocTkeData f90_data[] = {
      ShocTkeData(10, 71, 72, 300),
      ShocTkeData(10, 12, 13, 100),
      ShocTkeData(7,  16, 17, 50),
      ShocTkeData(2, 7, 8, 5),
    };

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ShocTkeData cxx_data[] = {
      ShocTkeData(f90_data[0]),
      ShocTkeData(f90_data[1]),
      ShocTkeData(f90_data[2]),
      ShocTkeData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      shoc_tke(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      shoc_tke_f(d.shcol, d.nlev, d.nlevi, d.dtime, d.wthv_sec, d.shoc_mix, d.dz_zi, d.dz_zt,
                 d.pres, d.u_wind, d.v_wind, d.brunt, d.obklen, d.zt_grid, d.zi_grid, d.pblh,
                 d.tke, d.tk, d.tkh, d.isotropy);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(ShocTkeData);
      for (Int i = 0; i < num_runs; ++i) {
        ShocTkeData& d_f90 = f90_data[i];
        ShocTkeData& d_cxx = cxx_data[i];
        REQUIRE(d_f90.total(d_f90.tke) == d_cxx.total(d_cxx.tke));
        REQUIRE(d_f90.total(d_f90.tke) == d_cxx.total(d_cxx.tk));
        REQUIRE(d_f90.total(d_f90.tke) == d_cxx.total(d_cxx.tkh));
        REQUIRE(d_f90.total(d_f90.tke) == d_cxx.total(d_cxx.isotropy));
        for (Int k = 0; k < d_f90.total(d_f90.tke); ++k) {
          REQUIRE(d_f90.tke[k] == d_cxx.tke[k]);
          REQUIRE(d_f90.tk[k] == d_cxx.tk[k]);
          REQUIRE(d_f90.tkh[k] == d_cxx.tkh[k]);
          REQUIRE(d_f90.isotropy[k] == d_cxx.isotropy[k]);
        }
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

  TestStruct::run_property();
}

TEST_CASE("shoc_tke_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocTke;

  TestStruct::run_bfb();

}

} // namespace
