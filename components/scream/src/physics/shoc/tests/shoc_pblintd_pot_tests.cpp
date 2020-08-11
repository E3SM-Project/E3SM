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
#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

TEST_CASE("shoc_tke_shr_prod", "shoc") {
  constexpr Int shcol    = 2;
  constexpr Int nlev     = 1;
  
  // Tests for the SHOC function:
  //   pblintd_init_pot

  // FIRST TEST
  //  Dry atmosphere test.  Verify that in a dry atmosphere
  //  with no water vapor or liquid water loading that the output
  //  thv (virtual potential temperature) is equal to the input
  //  thl (liquid water potential temperature).  

  // Define the liquid water potential temperature [K]
  Real thl_dry[shcol] = {300.0, 295.0};
  // Define the water vapor [kg/kg]
  Real qv = 0.0;
  // Define the liquid water mixing ratio [kg/kg]
  Real ql = 0.0;

  // Initialize data structure for bridging to F90
  SHOCPblpotData SDS(shcol, nlev);

  // Test that the inputs are reasonable.
  REQUIRE(SDS.shcol > 0);

  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      SDS.thl[offset] = thl_dry[n];
      SDS.ql[offset] = ql;
      SDS.q[offset] = qv;
    }
  }

  // Check that the inputs are expected
  for(Int s = 0; s < SDS.shcol; ++s) {
    for (Int n = 0; n < SDS.nlev; ++n){
      const auto offset = n + s * SDS.nlevi;
      // Make sure top level dz_zi value is zero
      REQUIRE(SDS.thl[offset] > 0.0);
      // For this test require ql and q be zero
      REQUIRE(SDS.q[offset] == 0.0);
      REQUIRE(SDS.ql[offset] == 0.0);
    }
  }

  // Call the fortran implementation
  pblintd_init_pot(nlev, SDS);

  // Check the result.  
  // Verify that virtual potential temperature is exactly
  //   equal to liquid water potential temperature 
  for(Int s = 0; s < SDS.shcol; ++s) {
    for (Int n = 0; n < SDS.nlev; ++n){
      const auto offset = n + s * SDS.nlevi;
      // Make sure top level dz_zi value is zero
      REQUIRE(SDS.thl[offset] == SDS.thv[offset]);
    }
  }  
    
  // SECOND TEST 
  // For two parcels with identical inputs, but one with condensate
  //  loading, verify that the parcel with condensate has a lower
  //  virtual potential temperature than the other.

  // Define the liquid water potential temperature [K]
  Real thl_parcel[shcol] = {290.0, 290.0};
  // Define the water vapor [g/kg]
  Real qv_parcel = 10.0;
  // Define the liquid water mixing ratio [g/kg]
  Real ql_parcel[shcol] = {0.4, 0.0};  

  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      SDS.thl[offset] = thl_parcel[n];
      // convert the following to kg/kg
      SDS.ql[offset] = ql_parcel[n]/1000.0;
      SDS.q[offset] = qv_parcel/1000.0;
    }
  }  
  
  // Check that the inputs are expected
  for(Int s = 0; s < SDS.shcol; ++s) {
    for (Int n = 0; n < SDS.nlev; ++n){
      const auto offset = n + s * SDS.nlevi;
      // Make sure top level dz_zi value is zero
      REQUIRE(SDS.thl[offset] > 0.0);
      // Make sure unit conversion resulted in 
      //   reasonable magnitudes
      REQUIRE(SDS.q[offset] < 0.1);
      REQUIRE(SDS.ql[offset] < 1.0);
    }
  }
  
  // Call the fortran implementation
  pblintd_init_pot(nlev, SDS);
  
  // Check test
  // Verify that column with condensate loading 
  //   results in a lower virtual potential temperature
  for(Int s = 0; s < SDS.shcol-1; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      // Get value corresponding to next column
      const auto offsets = n + (s+1) * SDS.nlev;
      if(SDS.ql[offset] > SDS.ql[offsets]){
        REQUIRE(SDS.thv[offset] < SDS.thv[offsets]);
      }
      else{
        REQUIRE(SDS.thv[offset] > SDS.thv[offsets]);
      }
    }
  }
 
}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
