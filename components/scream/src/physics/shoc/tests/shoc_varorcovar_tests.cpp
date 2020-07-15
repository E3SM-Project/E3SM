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

TEST_CASE("shoc_varorcovar", "shoc") {
  constexpr Int shcol    = 2;
  constexpr Int nlev     = 4;
  constexpr auto nlevi   = nlev + 1;
   
  //NOTE: This routine does not compute the (co)variance
  // for boundary points.  Input Grid values that are never used
  // are set to ZERO so as not to confuse the test data
  // that is actually used (and is also an implicit test to be
  // sure that operations are not done at the boundary).
  // Output values at boundary set to something other than 
  // zero to check that they are not modified.
  
  //IN: tunefac, isotropy_zi, tkh_zi, dz_zi, 
  //IN: invar1, invar2
  //OUT: varorcovar
    
  // Define delta z on the nlevi grid [m]
  Real dz_zi[nlevi] = {0.0, 100.0, 50.0, 20.0, 0.0};
  // Eddy coefficients on the nlevi grid [m2s-1]
  Real tkh_zi[nlevi] = {0.0, 10.0, 10.0, 10.0, 0.0};
  // Return to isotropy timescale [s]
  Real isotropy_zi[nlevi] = {0.0, 100.0, 100.0, 100.0, 0.0};
  // Invar1 (in this example we use potential temperature) [K]
  //   this variable should be on the nlev grid
  Real invar1[nlev] = {320.0, 315.0, 315.0, 316.0};
  // Invar2 this variable should be on the nlev grid
  Real invar2[nlev];
  // Initialize (co)variance on the nlevi grid, set boundaries
  //   to large values to make sure they are not modified
  Real varorcovar[nlevi] = {100., 0., 0., 0., 100.};
  
  // First test we set invar2 equal to invar 1 to test
  //  the variance calculation.
  invar2[:] = invar1[:];  

  // Initialzie data structure for bridgeing to F90
  SHOCVarorcovarData SDS(shcol, nlev, nlevi);

  // Test that the inputs are reasonable.
  REQUIRE(SDS.nlevi - SDS.nlev == 1);
  REQUIRE(SDS.shcol > 0);
  
  // Fill in test data 
  for(Int s = 0; s < SDS.shcol; ++s) {
    // First on the nlev grid
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      
      SDS.invar1[offset] = invar1[n]; 
      SDS.invar2[offset] = invar2[n];
    }
    
    // Now for data on the nlevi grid
    for(Int n = 0; n < SDS.nlevi; ++n) {
      const auto offset = n + s * SDS.nlevi;
      
      SDS.tkh_zi[offset] = tkh_zi[n];
      SDS.isotropy_zi[offset] = isotropy_zi[n];
      SDS.dz_zi[offset] = dz_zi[n];
      SDS.varorcovar[offset] = varorcovar[n];
    }
  }

  // Check that the inputs make sense
  
  // Check to make sure that dz_zi, tkh_zi, and isotropy_zi
  //  (outside of the boundaries) are greater than zero 
  for(Int s = 0; s < SDS.shcol; ++s) {
    // do NOT check boundaries!
    for(Int n = 1; n < SDS.nlevi; ++n) {
      const auto offset = n + s * SDS.nlev;
      REQUIRE(SDS.dz_zi[offset] > 0.0);
      REQUIRE(SDS.tkh_zi[offset] > 0.0);
      REQUIRE(SDS.isotropy_zi[offset] > 0.0);
    }
  }
  
  // Call the fortran implementation for variance
  calc_shoc_varorcovar(nlev, SDS);
  
  // Check the results
  for(Int s = 0; s < SDS.shcol; ++s) {   
    for(Int n = 0; n < SDS.nlevi; ++n) {
      const auto offset = n + s * SDS.nlevi;

      // validate that the boundary points 
      //   have NOT been modified      
      if (n == 0 || n == nlevi){
        REQUIRE(SDS.varorcovar[offset] == 100.0);
      }
      else{
      
      // Validate that all values are greater to 
      //   or equal to zero 
      
        REQUIRE(SDS.varorcovar[offset] >= 0.0);
      
        // well mixed layer test
        if ((invar1[n-1] - invar1[n]) == 0.0){
          REQUIRE(SDS.varorcovar[offset] == 0.0);
        }
      }
    } 
  }

  // Now we check the covariance calculation
  // Change input values of invar2 to be total water [g/kg]
  Real invar2_2[nlev] = {20.0, 15.0, 10.0, 5.0};
  
  // convert to [kg/kg]
  invar2_2[:] = invar2_2[:]/1000.0;

  // Update invar2 to be total water
  for(Int s = 0; s < SDS.shcol; ++s) {
    // First on the nlev grid
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
       
      SDS.invar2[offset] = invar2_2[n];
    }   
  }
  
  // Call the fortran implementation for covariance
  calc_shoc_varorcovar(nlev, SDS);
  
  // Check the results
  for(Int s = 0; s < SDS.shcol; ++s) {   
    for(Int n = 0; n < SDS.nlevi; ++n) {
      const auto offset = n + s * SDS.nlevi;

      // validate that the boundary points 
      //   have NOT been modified      
      if (n == 0 || n == nlevi){
        REQUIRE(SDS.varorcovar[offset] == 100.);
      }
      else{  
      
        // well mixed layer test
        if ((invar1[n-1] - invar1[n]) == 0.0 || \
	    (invar2[n-1] - invar2[n]) == 0.0){
          REQUIRE(SDS.varorcovar[offset] == 0.0);
        }
	
	// validate values are negative if potential
	//  temperature decreases with height
        if ((invar1[n-1] - invar1[n]) < 0.0 && \
	    (invar2[n-1] - invar2[n]) > 0.0){
          REQUIRE(SDS.varorcovar[offset] < 0.0);
        }	
	
      }
    } 
  }
  
  // Next check that vertical differencing is done correctly
  //  Assume that isotropy and tkh are constant with height
  //  Assign values of potential temperature that are increasing
  //  at a constant rate per GRID box.  Assign dz that INCREASES
  //  with height and check that varorcovar DECREASES with height
  
  Real invar1_2[nlev] = {320.0, 315.0, 310.0, 305.0};
  
  // Update invar1 and invar2 to be identical
  for(Int s = 0; s < SDS.shcol; ++s) {
    // First on the nlev grid
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
       
      SDS.invar1[offset] = invar1_2[n];
      SDS.invar2[offset] = invar1_2[n];
    }   
  }  

  // Call the fortran implementation for variance
  calc_shoc_varorcovar(nlev, SDS);
  
  // Check the results
  for(Int s = 0; s < SDS.shcol; ++s) {   
    for(Int n = 1; n < SDS.nlev-1; ++n) {
      const auto offset = n + s * SDS.nlevi;  

      // Validate that values of varorcovar
      //  are decreasing with height
      REQUIRE(SDS.varorcovar[offset]-\
              SDS.varorcovar[offset+1] < 0.0);        
    
    }
  }
  
}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
