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

TEST_CASE("diag_second_moments", "shoc") {
  constexpr Int shcol    = 2;
  constexpr Int nlev     = 5;
  constexpr auto nlevi   = nlev + 1;
  constexpr int num_tracer = 2;
   
  //Tests for the subroutine diag_second_moments.  
  //  This routine calls the routines calc_shoc_vertflux
  //  and calc_shoc_varorcovar.  These routines test the physics
  //  of these routines, thus this test is provided mainly
  //  to ensure input/output correctness and not necessarily
  //  to test the physical formulation.  

  // Define liquid water potential temperature [K]
  Real thetal[nlev] = {310.0, 308.0, 305.0, 300.0, 301.0};
  // Define total water mixint ratio [g/kg]
  Real qw[nlev] = {16.0, 17.0, 18.0, 19.0, 20.0};
  // Define tracer [no units]
  Real tracer[nlev] = {100.0, 40.0, 60.0, 700.0, 3.};
  // Define zonal wind on nlev grid [m/s]
  Real u_wind[nlev] = {2.0, 1.0, 0.0, -1.0, -2.0};
  // Define meridional wind on nlev grid [m/s]
  Real v_wind[nlev] = {1.0, 2.0, 3.0, 4.0, 5.0};
  // Length Scale [m]
  Real shoc_mix[nlev] = {1000.0, 750.0, 500.0, 400.0, 300.0};    
  // Define the midpoint height grid [m]
  Real zt_grid[nlev] = {450.0, 350.0, 250.0, 150.0, 50.0};
  // Define the interface height grid [m]
  Real zi_grid[nlevi] = {500.0, 400.0, 300.0, 200.0, 100.0, 0.0};
  // Define thickness on the interface grid [m]
  Real dz_zi[nlevi] = {0.0, 100.0, 100.0, 100.0, 100.0, 50.0};  
  // Define TKE initial [m2/s2]
  Real tke[nlev] = {0.1, 0.2, 0.4, 0.02, 0.3};
  // Define initial eddy diffusivities
  Real tkh[nlev] = {1.0, 1.0, 1.0, 1.0, 1.0};
  // Define the return to isotropy timescale [s]
  Real isotropy[nlev] = {100.0, 200.0, 300.0, 200.0, 100.0};
  Real tk[nlev]; 
  
  // Convert total water to [kg/kg]
  qw[:]=qw[:]/1000.0;
  
  // Set eddy diffusivity of momentum equal to that of heat
  tk[:]=tkh[:]; 
  
  // NOTE: the input profile invar contains data to test conditions
  //   of 1) and unstable boundary layer, 2) well mixed
  //   layer and 3) conditionally stable layer.  

  // Initialzie data structure for bridgeing to F90
  SHOCdiagSecondMomsSDS(shcol, nlev, nlevi,numtracer);

  // Test that the inputs are reasonable.
  REQUIRE(SDS.nlevi - SDS.nlev == 1);
  REQUIRE(SDS.shcol > 0);
  
  // Fill in test data 
  for(Int s = 0; s < SDS.shcol; ++s) {
    // First on the nlev grid
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      
      SDS.thetal[offset] = thetal[n];
      SDS.qw[offset] = qw[n];
      SDS.tracer[offset] = tracer[n];
      SDS.u_wind[offset] = u_wind[n];
      SDS.v_wind[offset] = v_wind[n];
      SDS.shoc_mix[offset] = shoc_mix[n];
      SDS.zt_grid[offset] = zt_grid[n];
      SDS.tke[offset] = tke[n];
      SDS.tkh[offset] = tkh[n];
      SDS.tk[offset] = tk[n];
      SDS.isotropy[offset] = isotropy[n]; 
      
      // Tracer input
      for (Int t = 0; t < SDS.num_tracer; ++t){
        // unsure about this
        const auto offset_tracer = t + n + s * SDS.nlev;
	SDS.tracer[offset_tracer] = tracer[n];
      }
    }
    
    // Now for data on the nlevi grid
    for(Int n = 0; n < SDS.nlevi; ++n) {
      const auto offset = n + s * SDS.nlevi;
      
      SDS.zi_grid[offset] = zi_grid[n];
      SDS.dz_zi[offset] = dz_zi[n];
      
      // Initialize the second moments
      SDS.thl_sec[offset] = 0.0;
      SDS.qw_sec[offset] = 0.0;
      SDS.qwthl_sec[offset] = 0.0;
      SDS.wthl_sec[offset] = 0.0;
      SDS.wqw_sec[offset] = 0.0;
      SDS.uw_sec[offset] = 0.0;
      SDS.vw_sec[offset] = 0.0;
      SDS.wtke_sec[offset] = 0.0;

      // Tracer input
      for (Int t = 0; t < SDS.num_tracer; ++t){
        // unsure about this
        const auto offset_tracer = t + n + s * SDS.nlevi;
	SDS.tracer[offset_tracer] = 0.0;
      }      
      
    }
  }

  // Check that the inputs make sense
  
  // Check to make sure that dz_zi are tkh_zi
  //  (outside of the boundaries) are greater than zero 
  for(Int s = 0; s < SDS.shcol; ++s) {
    // do NOT check boundaries!
    for(Int n = 1; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      
      REQUIRE(SDS.tkh[offset] > 0.0);
      REQUIRE(SDS.tk[offset] > 0.0);
      REQUIRE(SDS.tke[offset] > 0.0);
      REQUIRE(SDS.shoc_mix[offset] > 0.0);
    }
    for(Int n = 1; n < SDS.nlevi; ++n) {
      const auto offset = n + s * SDS.nlev;
      REQUIRE(SDS.dz_zi[offset] > 0.0);
      REQUIRE(SDS.tkh[offset] > 0.0);
    }    
    
  }
  
  // Call the fortran implementation
  diag_second_moments(nlev, SDS);

}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
