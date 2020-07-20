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

TEST_CASE("shoc_tke_isotropic_ts", "shoc") {
  constexpr Int shcol    = 2;
  constexpr Int nlev     = 1;
  
  // Tests for the subroutine isotropic_ts
  //   in the SHOC TKE module. 

  // For this routine all inputs are on midpoint grid and
  //  there are no vertical derivatives, therefore we will
  //  consider one vertical level per test.  Each column will
  //  be loaded up with a different test.
  
  // FIRST TEST
  // Stability test.  Verify that two points with the same inputs
  //  but one with brunt > 0 and the other with brunt < 0 that the
  //  timescale with brunt > 0 is less.   
  
  // Integrated brunt vaisalla 
  Real brunt_int_st[shcol] = {0.1, 0.1};
  // TKE [m2/s2]
  Real tke_st[shcol] = {0.4, 0.4};
  // Dissipation rate [m2/s3]
  Real diss_st[shcol] = {0.1, 0.1};
  // Brunt Vaisalla frequency [/s]
  Real brunt_st[shcol] = {0.004, -0.004}; 

  // Initialzie data structure for bridgeing to F90
  SHOCIsoTsData SDS(shcol, nlev);

  // Test that the inputs are reasonable.
  REQUIRE(SDS.shcol > 0);

  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    SDS.brunt_int[s] = brunt_int_st[s];
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      SDS.tke[offset] = tke_st[n];
      SDS.diss[offset] = diss_st[n];
    }
  }

  // Check that the inputs make sense
  for(Int s = 0; s < SDS.shcol; ++s) {
    for (Int n = 0; n < SDS.nlev; ++n){
      const auto offset = n + s * SDS.nlevi;
      // Should be greater than zero 
      REQUIRE(SDS.tke[offset] > 0.0);
      REQUIRE(SDS.diss[offset] > 0.0);
    }
  }

  // Call the fortran implementation
  isotropic_ts(nlev, SDS);

  // Check to make sure that column with positive 
  //  brunt vaisalla frequency is smaller
  for(Int s = 0; s < SDS.shcol-1; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      // Get value corresponding to next column
      const auto offsets = n = (s+1) * SDS.nlev;     
      if(SDS.brunt[offset] < 0.0 & SDS.brunt[offsets] > 0.0){
        REQUIRE(SDS.isotropy[offset] < SDS.isotropy[offsets]);
      }
    }
  }  

  // SECOND TEST
  // Dissipation test.  Verify that two points with the same inputs
  //  but one with smaller dissipation rate, that that point has 
  //  a higher isotropy timescale   
  
  // Integrated brunt vaisalla 
  Real brunt_int_diss[shcol] = {0.1, 0.1};
  // TKE [m2/s2]
  Real tke_diss[shcol] = {0.4, 0.4};
  // Dissipation rate [m2/s3]
  Real diss_diss[shcol] = {0.1, 0.2};
  // Brunt Vaisalla frequency [/s]
  Real brunt_diss[shcol] = {0.004, 0.004}; 
  
  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    SDS.brunt_int[s] = brunt_int_diss[s];
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      SDS.tke[offset] = tke_diss[n];
      SDS.diss[offset] = diss_diss[n];
    }
  }

  // Check that the inputs make sense
  for(Int s = 0; s < SDS.shcol; ++s) {
    for (Int n = 0; n < SDS.nlev; ++n){
      const auto offset = n + s * SDS.nlevi;
      // Should be greater than zero 
      REQUIRE(SDS.tke[offset] > 0.0);
      REQUIRE(SDS.diss[offset] > 0.0);
    }
  }   
  
  // Call the fortran implementation
  isotropic_ts(nlev, SDS);   
  
  // Check to make sure that column with positive 
  //  brunt vaisalla frequency is smaller
  for(Int s = 0; s < SDS.shcol-1; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      // Get value corresponding to next column
      const auto offsets = n = (s+1) * SDS.nlev;     
      if(SDS.diss[offset] < SDS.diss[offsets]){
        REQUIRE(SDS.isotropy[offset] > SDS.isotropy[offsets]);
      }
    }
  }  
  
}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
