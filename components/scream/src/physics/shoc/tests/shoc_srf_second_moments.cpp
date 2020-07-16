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

TEST_CASE("shoc_srf_second_moments", "shoc") {
  constexpr Int shcol    = 5;
   
  //NOTE: This routine computes the surface properties
  // needed for the lower boundary condition for the
  // second order moments.  Returned values are 
  // surface friction velocity and surface convective
  // velocity.
 
  // Here all inputs are shcol only, so we will test
  //   out a variety of conditions with this code
 
  // Surface sensible heat flux [K m/s] 
  Real wthl_sfc[shcol] = {0.01, -0.5, -0.3, 0.02};
  // Surface momentum flux (u-direction) [m2/s2]
  Real uw_sfc[shcol] = {0.0, -0.1, 0.2, -0.2};
  // Surface momentum flux (v-direction) [m2/s2]
  Real vw_sfc[shcol] = {0.0, -0.2, 0.1, 0.1};

  // Initialzie data structure for bridgeing to F90
  SHOCsrfsecondmomsData SDS(shcol);

  // Test that the inputs are reasonable.
  REQUIRE(SDS.shcol > 0);
  
  // Fill in test data 
  for(Int s = 0; s < SDS.shcol; ++s) {      
    SDS.wthl_sfc[s] = wthl_sfc; 
    SDS.uw_sfc[s] = uw_sfc;
    SDS.vw_sfc[s] = vw_sfc;
  }
  
  // Call the fortran implementation for variance
  diag_second_moments_srf(nlev, SDS);
  
  // Check the results
  
  for(Int s = 0; s < SDS.shcol; ++s) {   
    // if surface momentum fluxes are zero make
    //  sure that friction velocity is also zero
    if((SDS.uw_sfc[s] == 0.0) && 
       (SDS.uw_sfc[s] == 0.0)){
      REQUIRE(SDS.ustar2[s] == 0.0);         
    } 
    else{
      // if either momentum fluxes are not zero then
      //   make sure friction velocity is positive
      REQUIRE(SDS.ustar2[s] > 0.0);
    }
    
    // if surface heat flux is less than zero 
    //  then verify that convective scale is zero
    if(SDS.wthl_sfc[s] < 0.0){
      REQUIRE(SDS.wstar[s] == 0.0);
    }
    else{
      // otherwise, make sure it is always positive
      REQUIRE(SDS.wstar[s] > 0.0);
    }
  }

}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
