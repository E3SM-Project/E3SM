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

TEST_CASE("shoc_vertflux", "shoc") {
  constexpr Int shcol    = 2;
  constexpr Int nlev     = 4;
  constexpr auto nlevi   = nlev + 1;
   
  //NOTE: This routine does not compute the vertical fluxes
  // for boundary points.  Input Grid values that are never used
  // are set to ZERO so as not to confuse the test data
  // that is actually used (and is also an implicit test to be
  // sure that operations are not done at the boundary).
  // Output values at boundary set to something other than 
  // zero to check that they are not modified.
    
  // Define delta z on the nlevi grid [m]
  Real dz_zi[nlevi] = {0, 100.0, 50.0, 20.0, 0.0};
  // Eddy coefficients on the nlevi grid [m2s-1]
  Real tkh_zi[nlevi] = {0.0, 10.0, 10.0, 10.0, 0.0};
  // Invar (in this example we use potential temperature) [K]
  //   this variable should be on the nlev grid
  Real invar[nlev] = {320.0, 315.0, 315.0, 316.0};
  // Initialize vertflux on the nlevi grid, set boundaries
  //   to large values to make sure they are not modified
  Real vertflux[nlevi] = {100., 0., 0., 0., 100.};
  
  // NOTE: the input profile invar contains data to test conditions
  //   of 1) and unstable boundary layer, 2) well mixed
  //   layer and 3) conditionally stable layer.  

  // Initialzie data structure for bridgeing to F90
  SHOCVertfluxData SDS(shcol, nlev, nlevi);

  // Test that the inputs are reasonable.
  REQUIRE(SDS.nlevi - SDS.nlev == 1);
  REQUIRE(SDS.shcol > 0);
  
  // Fill in test data 
  for(Int s = 0; s < SDS.shcol; ++s) {
    // First on the nlev grid
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      
      SDS.invar[offset] = invar[n]; 
    }
    
    // Now for data on the nlevi grid
    for(Int n = 0; n < SDS.nlevi; ++n) {
      const auto offset = n + s * SDS.nlevi;
      
      SDS.tkh_zi[offset] = tkh_zi[n];
      SDS.dz_zi[offset] = dz_zi[n];
      SDS.vertflux[offset] = vertflux[n];
    }
  }

  // Check that the inputs make sense
  
  // Check to make sure that dz_zi are tkh_zi
  //  (outside of the boundaries) are greater than zero 
  for(Int s = 0; s < SDS.shcol; ++s) {
    // do NOT check boundaries!
    for(Int n = 1; n < SDS.nlevi; ++n) {
      const auto offset = n + s * SDS.nlev;
      REQUIRE(SDS.dz_zi[offset] > 0.0);
      REQUIRE(SDS.tkh_zi[offset] > 0.0);
    }
  }
  
  // Call the fortran implementation
  calc_shoc_vertflux(nlev, SDS);
  
  // Check the results
  for(Int s = 0; s < SDS.shcol; ++s) {   
    for(Int n = 0; n < SDS.nlevi; ++n) {
      const auto offset = n + s * SDS.nlevi;

      // validate that the boundary points 
      //   have NOT been modified      
      if (n == 0 || n == nlevi){
        REQUIRE(SDS.vertflux[offset] == 100.0);
      }
      else{
      
      // Validate Downgradient assumption for 
      //   various possible layers
      
        // conditionally stable layer
        if ((invar[n-1] - invar[n]) > 0.0){
          REQUIRE(SDS.vertflux[offset] < 0.0);
        } 
        // well mixed layer
        if ((invar[n-1] - invar[n]) == 0.0){
          REQUIRE(SDS.vertflux[offset] == 0.0);
        }
        // unstable layer
        if ((invar[n-1] - invar[n]) < 0.0){
          REQUIRE(SDS.vertflux[offset] > 0.0);
        }
      }
    } 
  }

}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
