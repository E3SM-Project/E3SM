#include "catch2/catch.hpp"

//#include "share/scream_types.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/util/scream_utils.hpp"
#include "ekat/scream_kokkos.hpp"
#include "ekat/scream_pack.hpp"
#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/util/scream_arch.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"

#include "physics/common/physics_constants.hpp"

#include "shoc_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace shoc {
namespace unit_test {

TEST_CASE("shoc_grid", "shoc"){

  constexpr Real gravit = scream::physics::Constants<Real>::gravit;

  const Int shcol = 1; 
  const Int nlev = 128; 
  const auto nlevi = nlev + 1;

  std::array<Real, nlev * shcol> pdel, zt_grid, dz_zt, rho_zt;
  std::array<Real, nlevi * shcol> zi_grid, dz_zi; 

  SHOCGridData SDS(shcol, nlev, nlevi);

  //std::cout << "PRINT in test \t" << SDS.shcol << "\n";

  //SDS.shcol = shcol; 
  //SDS.nlev = nlev; 
  //SDS.nlevi = nlevi;

  //Copy pointers from arrays into SDS
  //SDS.pdel = pdel.data();
  //SDS.zt_grid = zt_grid.data();
  //SDS.dz_zt = dz_zt.data();
  //SDS.rho_zt = rho_zt.data();
  //SDS.zi_grid = zi_grid.data(); 
  //SDS.dz_zi = dz_zi.data();


  //const Real dz = 50.0;
  // Fill in test data on zt_grid
  //for(Int n=nlev-1; n >= 0; --n){
  //  for(Int s=0; s < shcol; ++s){

  //    const auto p = s + n*shcol;
 
  //    SDS.zt_grid[p] = dz*n + dz/2;   
  //    SDS.pdel[p] = -gravit*dz; 
  //  }
 // }

  // Fill in test data on zi_grid
  //for(Int n=nlevi-1; n >= 0; --n){
  //  for(Int s=0; s < shcol; ++s){
  //     const auto p = s + n*shcol;
  //     SDS.zi_grid[p] = 0.0;
  //    }
  //}

  // Test that the inputs are resonable.
  //REQUIRE(SDS.nlevi - SDS.nlev == 1); 
  //REQUIRE(SDS.shcol>0);

  //Call the fortran implementation
  shoc_grid(nlev, SDS);



}


}//namespace unit_test
}//namespace shoc
}//namespace scream

