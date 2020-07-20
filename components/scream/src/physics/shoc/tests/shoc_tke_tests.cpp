#include "catch2/catch.hpp"

//#include "share/scream_types.hpp"
#include <algorithm>
#include <array>
#include <random>
#include <thread>

#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/scream_arch.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/util/scream_utils.hpp"
#include "physics/common/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

TEST_CASE("shoc_tke", "shoc") {
  constexpr Real gravit  = scream::physics::Constants<Real>::gravit;
  constexpr Int shcol    = 2;
  constexpr Int nlev     = 5;
  constexpr auto nlevi   = nlev + 1;
  
  // Tests for the subroutine shoc_tke.  This is the top
  //  level function for the TKE module in SHOC.  Thus since
  //  all associated functions here contain property tests, the
  //  purpose here is mostly to expose whether errors have been 
  //  made in the input/output fields in these routines  

  // For this test we will perform a TKE growth test.  
  //  Will provide inputs that should encourage TKE to grow
  //  at each level.  

  // Define height thickness on nlevi grid [m]
  //   NOTE: First indicee is zero because it is never used
  //   Do a stretched grid

  // Timestep [s]
  Real dtime = 300.0;
  // Buoyancy flux [K m/s]
  Real wthv_sec[nlev] = {0.5, 0.4, 0.3, 0.2, 0.3};
  // Length Scale [m]
  Real shoc_mix[nlev] = {1000.0, 750.0, 500.0, 400.0, 300.0};
  // Define zonal wind on nlev grid [m/s]
  Real u_wind[nlev] = {2.0, 1.0, 0.0, -1.0, -2.0};
  // Define meridional wind on nlev grid [m/s]
  Real v_wind[nlev] = {1.0, 2.0, 3.0, 4.0, 5.0};
  // Define Obukov length [m] 
  Real obukov = -1.0; 
  // Define thickness on the interface grid [m]
  Real dz_zi[nlevi] = {0.0, 100.0, 100.0, 100.0, 100.0, 50.0};
  // Define thickness on the thermo grid [m]
  Real dz_zt[nlev] = {100.0, 100.0, 100.0, 100.0, 100.0};
  // Define the midpoint height grid [m]
  Real zt_grid[nlev] = {450.0, 350.0, 250.0, 150.0, 50.0};
  // Define the interface height grid [m]
  Real zi_grid[nlevi] = {500.0, 400.0, 300.0, 200.0, 100.0, 0.0};
  // Define PBL height
  Real pblh = 1000.0; 
  // Define pressure [hPa]
  Real pres[nlev] = {850.0, 875.0, 900.0, 950.0, 1000.0};
  // Define brunt Vaisalla frequency
  Real brunt[nlev] = {0.1, 0.1, 0.1, 0.1, 0.1};
  
  // Define TKE initial [m2/s2]
  Real tke[nlev] = {0.01, 0.01, 0.01, 0.01, 0.01};
  // Define initial eddy diffusivities
  Real tkh[nlev] = {1.0, 1.0, 1.0, 1.0, 1.0};
  Real tk[nlev];
  
  // Convert pressure to [Pa]
  pres[:]=pres[:]*100.0;
  tk[:]=tkh[:];

  // Initialzie data structure for bridgeing to F90
  SHOCTKEData SDS(shcol, nlev, nlevi);

  // Test that the inputs are reasonable.
  REQUIRE(SDS.nlevi - SDS.nlev == 1);
  REQUIRE(SDS.shcol > 0);

  SDS.dtime = dtime;
  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    SDS.pblh[s] = pblh;
    SDS.obukov[s] = obukov;
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      SDS.wthv_sec[offset] = wthv_sec[n];
      SDS.shoc_mix[offset] = shoc_mix[n];
      SDS.u_wind[offset] = u_wind[n];
      SDS.v_wind[offset] = v_wind[n];
      SDS.dz_zt[offset] = dz_zt[n];
      SDS.zt_grid[offset] = zt_grid[n];
      SDS.pres[offset] = pres[n];
      SDS.tke[offset] = tke[n];
      SDS.tkh[offset] = tkh[n];
      SDS.tk[offset] = tk[n];
    }

    // Fill in test data on zi_grid.
    for(Int n = 0; n < SDS.nlevi; ++n) {
      const auto offset   = n + s * SDS.nlevi;
      SDS.dz_zi[offset] = dz_zi[n];
      SDS.zi_grid[offset] = zi_grid[n];
    }
  }

  // Check that the inputs make sense

  for(Int s = 0; s < SDS.shcol; ++s) {
    // nlevi loop
    for (Int n = 0; n < SDS.nlevi; ++n){
      const auto offset = n + s * SDS.nlevi;
      // Make sure top level dz_zi value is zero
      if (n == 0){
        REQUIRE(SDS.dz_zi[offset] == 0.0);
      } 
      // Otherwise, should be greater than zero
      else{
        REQUIRE(SDS.dz_zi[offset] > 0.0);
      }
      // Check that zi increases updward
      REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0.0);
    }
    // nlev loop 
    for (Int n = 0; n < SDS.nlev; ++n){
      const auto offset = n + s * SDS.nlev;
      // Check that zt increases upward
      REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0.0);
      REQUIRE(SDS.dz_zt[offset] > 0.0);
      REQUIRE(SDS.shoc_mix[offset] > 0.0);
      REQUIRE(SDS.tke[offset] > 0.0);
      REQUIRE(SDS.tkh[offset] > 0.0);
      REQUIRE(SDS.tk[offset] > 0.0);
      REQUIRE(SDS.pres[offset] > 0.0);
    }
  }

  // Call the fortran implementation
  shoc_tke(nlev, SDS);

  // Check test
  // Make sure that TKE has increased everwhere relative 
  //   to initial value.  Also make sure that outputs fall 
  //   within some reasonable bounds.
  
  for(Int s = 0; s < shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlevi;
      REQUIRE(SDS.tke[offset] > tke[n]);
      REQUIRE(SDS.tke[offset] > 0.0 &
              SDS.tke[offset] <= 50.0);
      REQUIRE(SDS.tkh[offset] > 0.0);
      REQUIRE(SDS.tk[offset] > 0.0);
      REQUIRE(SDS.isotropy[offset] > 0.0); 
    }
  }
  
}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
