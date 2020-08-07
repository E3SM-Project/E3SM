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

TEST_CASE("shoc_l_inf_length", "shoc") {
  constexpr Int shcol    = 2;
  constexpr Int nlev     = 5;

  // Tests for the compute_l_inf_shoc_length scale, which is the
  //  asymptotic length scale parameter   
  
  // Test two columns with identical profiles, EXCEPT 
  //  one which has larger TKE values.  Verify that the column 
  //  with the larger TKE values produces a larger asymptotic 
  //  length scale.       

  // Grid difference centered on thermo grid [m]
  Real dz_zt[nlev] = {100.0, 100.0, 100.0, 100.0, 100.0};
  // Grid height centered on thermo grid [m]
  Real zt_grid[nlev] = {500.0, 400.0, 300.0, 200.0, 100.0};
  // Turbulent kinetic energy [m2/s2]
  Real tke[nlev] = {0.1, 0.15, 0.2, 0.25, 0.3};

  // Initialzie data structure for bridgeing to F90
  SHOCInflengthData SDS(shcol, nlev);

  // Test that the inputs are reasonable.
  //  At least two columns are needed for this test!
  REQUIRE(SDS.shcol > 1);

  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      SDS.dz_zt[offset] = dz_zt[n];
      SDS.zt_grid[offset] = zt_grid[n];
      // Testing identical columns but one with larger TKE
      //  set first column as "base" column, and the others
      //  to a larger value of TKE uniformly.  
      if (s == 0){
        SDS.tke[offset] = tke[n];
      }
      else{
        SDS.tke[offset] = 2.0*tke[n];
      }
    }
  }

  // Check that the inputs make sense

  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      // Be sure that relevant variables are greater than zero
      REQUIRE(SDS.dz_zt[offset] > 0.0);
      REQUIRE(SDS.tke[offset] > 0.0); 
      REQUIRE(SDS.zt_grid[offset] > 0.0);
      if (n < nlev){
        // check that zt increases upward
        REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0.0);       
      }      
    } 
  }

  // Call the fortran implementation
  compute_l_inf_shoc_length(nlev, SDS);

  // Check the results
  // Make sure that conv_vel is negative
  for(Int s = 0; s < SDS.shcol; ++s) {   
    REQUIRE(SDS.l_inf[s] > 0.0);
  } 
  
  // Make sure that l_inf is greater than l_inf
  //  in the first column 
  for (Int s = 1; s < SDS.shcol; ++s){
    REQUIRE(SDS.l_inf[s] > SDS.l_inf[0]);
  }

}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
