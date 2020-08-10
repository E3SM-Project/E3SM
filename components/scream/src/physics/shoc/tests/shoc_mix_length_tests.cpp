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

TEST_CASE("shoc_mix_length", "shoc") {
  constexpr Int shcol    = 3;
  constexpr Int nlev     = 5;

  // Tests for the SHOC function:
  // compute_shoc_mix_shoc_length
  
  // Multi column test will verify that 1) mixing length increases 
  // with height given brunt vaisalla frequency and 
  // TKE are constant with height and 2) Columns will increase with
  // TKE values to ensure that mixing length is also increasing
  
  // Define TKE [m2/s2] that will be used for each column
  Real tke_cons = 0.1;  
  // Define the brunt vasailla frequency [s-1]
  Real brunt_cons = 0.001;
  // Define the assymptoic length [m]
  Real l_inf = 100.0;
  // Define the overturning timescale [s]
  Real tscale = 300.0; 
  // Define the heights on the zt grid [m]
  Real zt_grid[nlev] = {5000.0, 3000.0, 2000.0, 1000.0, 500.0};

  // Initialize data structure for bridgeing to F90
  SHOCMixlengthData SDS(shcol, nlev);

  // Test that the inputs are reasonable.
  // For this test shcol MUST be at least 2
  REQUIRE(SDS.shcol > 1);

  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    SDS.l_inf[s] = l_inf;
    SDS.tscale[s] = tscale;
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      // for the second column, double the 
      //  amount of TKE
      SDS.tke[offset] = (1.0+s)*tke_cons;
      SDS.brunt[offset] = brunt_cons;
      SDS.zt_grid[offset] = zt_grid[n];
    }
  }

  // Check that the inputs make sense

  // Be sure that relevant variables are greater than zero
  for(Int s = 0; s < SDS.shcol; ++s) {
    REQUIRE(SDS.l_inf[s] > 0.0);
    REQUIRE(SDS.tscale[s] > 0.0);
    for(Int n = 0; n < SDS.nlev - 1; ++n) {
      const auto offset = n + s * SDS.nlev;
      REQUIRE(SDS.tke[offset] > 0.0);
      REQUIRE(SDS.zt_grid[offset] > 0.0);
      if (s < shcol-1){
        // Verify that tke is larger column by column
        const auto offsets = n + (s+1) * SDS.nlev;
        REQUIRE(SDS.tke[offset] < SDS.tke[offsets]);
      }
    }
    
    // Check that zt increases upward
    for(Int n = 0; n < SDS.nlev - 1; ++n) {
      const auto offset = n + s * SDS.nlev;
      REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0.0);
    }        
  }

  // Call the fortran implementation
  compute_shoc_mix_shoc_length(nlev, SDS);

  // Check the results
  for(Int s = 0; s < SDS.shcol; ++s) {   
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;   
      // Validate shoc_mix greater than zero everywhere
      REQUIRE(SDS.shoc_mix[offset] > 0.0);  
      if (s < shcol-1){
        // Verify that mixing length increases column by column
        const auto offsets = n + (s+1) * SDS.nlev;
        REQUIRE(SDS.shoc_mix[offset] < SDS.shoc_mix[offsets]);
      }   
    }
    
    // Check that mixing length increases upward
    for(Int n = 0; n < SDS.nlev - 1; ++n) {
      const auto offset = n + s * SDS.nlev;
      REQUIRE(SDS.shoc_mix[offset + 1] - SDS.shoc_mix[offset] < 0.0);
    }
    
  }

}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
