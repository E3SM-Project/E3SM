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

TEST_CASE("shoc_tke_adv_sgs_tke", "shoc") {
  constexpr Int shcol    = 2;
  constexpr Int nlev     = 1;
  
  // Tests for the subroutine adv_sgs_tke 
  //   in the SHOC TKE module.     
  
  // For this routine all inputs are on midpoint grid and
  //  there are no vertical derivatives, therefore we will
  //  consider one vertical level per test. Each column will
  //  be loaded up with a different test (growth, decay)  

  // FIRST TEST
  // TKE growth/decay test.  Given the correct conditions, verify
  //  that TKE grows/decays over a timestep.  For this we choose a
  //  large/small mixing length, positive buoyancy flux, positive
  //  shear term, and large eddy diffusivity.  
  
  // Values in the FIRST column represents the growth test
  // Values in the SECOND column represents the decay test 

  // Define timestep [s]
  Real dtime = 20.0;
  // SHOC mixing length [m]
  Real shoc_mix_gr[shcol] = {20000.0, 20.};
  // Buoyancy flux [K m/s]
  Real wthv_sec_gr[shcol] = {0.5, -0.5};
  // Shear production term [s-2]
  Real sterm_gr[shcol] = {0.5, 0.0};
  // TKE initial value
  Real tke_init_gr[shcol] = {0.004, 0.4};

  // Initialzie data structure for bridgeing to F90
  SHOCAdvTKEData SDS(shcol, nlev);

  // Test that the inputs are reasonable.
  REQUIRE(SDS.shcol > 0);

  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      SDS.dtime = 20.0;
      SDS.shoc_mix[offset] = shoc_mix_gr[s];
      SDS.wthv_sec[offset] = wthv_sec_gr[s];
      SDS.sterm[offset] = sterm_gr[s];
      SDS.tke[offset] = tke_init_gr[s];
    }
  }

  // Check that the inputs make sense
  for(Int s = 0; s < SDS.shcol; ++s) {
    for (Int n = 0; n < SDS.nlev; ++n){
      const auto offset = n + s * SDS.nlevi;
      // time step, mixing length, TKE values,
      // shear terms should all be greater than zero 
      REQUIRE(SDS.dtime > 0.0);
      REQUIRE(SDS.shoc_mix[offset] > 0.0);
      REQUIRE(SDS.tke[offset] > 0.0);
      REQUIRE(SDS.sterm[offset] > 0.0);
    }
  }

  // Call the fortran implementation
  adv_sgs_tke(nlev, SDS);
  
  // Check to make sure that there has been
  //  TKE growth
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      
      if (s == 0){
      REQUIRE(SDS.tke[offset] > tke_init_gr[s]);
      }
      else{
        REQUIRE(SDS.tke[offset] < tke_init_gr[s]);
      }
    }
  }
  
  // SECOND TEST
  // TKE Dissipation test.  Given input values that are identical
  //  in two columns, verify that the dissipation rate is higher
  //  when the length scale is lower.  

  // SHOC mixing length [m]
  Real shoc_mix_diss[shcol] = {1000.0, 500.};
  // Buoyancy flux [K m/s]
  Real wthv_sec_diss[shcol] = {0.1, 0.1};
  // Shear production term [s-2]
  Real sterm_diss[shcol] = {0.01, 0.01};
  // TKE initial value
  Real tke_init_diss[shcol] = {0.1, 0.1};      

  // Fill in test data on zt_grid.
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;

      SDS.shoc_mix[offset] = shoc_mix_diss[s];
      SDS.wthv_sec[offset] = wthv_sec_diss[s];
      SDS.sterm[offset] = sterm_diss[s];
      SDS.tke[offset] = tke_init_diss[s];
    }
  }

  // Check that the inputs make sense
  for(Int s = 0; s < SDS.shcol; ++s) {
    for (Int n = 0; n < SDS.nlev; ++n){
      const auto offset = n + s * SDS.nlevi;
      // time step, mixing length, TKE values,
      // shear terms should all be greater than zero 
      REQUIRE(SDS.dtime > 0.0);
      REQUIRE(SDS.shoc_mix[offset] > 0.0);
      REQUIRE(SDS.tke[offset] > 0.0);
      REQUIRE(SDS.sterm[offset] > 0.0);
    }
  }
  
  // Call the fortran implementation
  adv_sgs_tke(nlev, SDS);  
  
  // Check to make sure that the column with 
  //  the smallest length scale has larger 
  //  dissipation rate
  Real mix_small = 40000.0; // Initialize to something large
  Real diss_large = 0.0; // Initialize to something small
  for(Int s = 0; s < SDS.shcol; ++s) {
    for(Int n = 0; n < SDS.nlev; ++n) {
      const auto offset = n + s * SDS.nlev;
      
      if (SDS.shoc_mix[offset] < mix_small){
        REQUIRE(SDS.diss[offset] > diss_large);
	mix_small = SDS.shoc_mix[offset];
	diss_large = SDS.diss[offset];
      }
    }
  }  
  
}

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream
