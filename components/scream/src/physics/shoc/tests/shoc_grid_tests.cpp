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
  const Real density = 1.0; 






  std::array<Real, nlev * shcol> pdel, zt_grid, dz_zt, rho_zt;
  std::array<Real, nlevi * shcol> zi_grid, dz_zi; 

  SHOCGridData SDS(shcol, nlev, nlevi);

  // Test that the inputs are resonable.
  REQUIRE(SDS.nlevi - SDS.nlev == 1); 
  REQUIRE(SDS.shcol>0);

  const Real dz = 50.0;
  // Fill in test data on zt_grid
  for(Int n=0; n < SDS.nlev; ++n){
    for(Int s=0; s < SDS.shcol; ++s){
      const auto nft = SDS.nlev - 1 - n;
      
      SDS.zt_grid[s, n] = nft * dz + dz/2; 
      SDS.pdel[s,n] = density*gravit*dz; 
    }
  }

  // Fill in test data on zi_grid
  for(Int n=0; n < SDS.nlevi; ++n){
    for(Int s=0; s < SDS.shcol; ++s){
       const auto nft = SDS.nlevi - 1 - n;
       
       SDS.zi_grid[s,n] =  dz * nft;
      }
  }

  // Here we sould check that the inputs make sense whe

 // Check that zt decreases upward
 int nerror = 0; 
  for(Int n=0; n < SDS.nlev; ++n){
    {
      for(Int s=0; s < SDS.shcol; ++s){
        if (SDS.zt_grid[s,n+1] - SDS.zt_grid[s,n] > 0.0){
          nerror++;
        }
      }
    }
  }  
  REQUIRE(nerror == 0);
 
 //Check that zi decreases upward
  for(Int n=0; n < SDS.nlevi-1; ++n){
    {
      for(Int s=0; s < SDS.shcol; ++s){
        if (SDS.zi_grid[s,n+1] - SDS.zi_grid[s,n] > 0.0){
          nerror++;
        }
      }
    }
  }  
 REQUIRE(nerror == 0);

  //Call the fortran implementation
  shoc_grid(nlev, SDS);

  //First check that dz is correctå
  for(Int n=0; n < SDS.nlev; ++n){
    for(Int s=0; s < shcol; ++s){
      if (SDS.dz_zt[s,n] != dz){
        nerror++;
      }
    }
  } 
  REQUIRE(nerror == 0);
  
  for(Int s=0; s < shcol; ++s){
    REQUIRE(SDS.dz_zi[s,0] == 0);
    REQUIRE(SDS.dz_zi[s,SDS.nlevi-1] == dz/2.0);
  }


  for(Int n=1; n < SDS.nlevi-1; ++n){
    for(Int s=0; s < shcol; ++s){
      if (SDS.dz_zi[s,n] != dz){
        nerror++;
      }
    }
  }
  REQUIRE(nerror == 0);
  
  //Now check density 
  for(Int n=0; n < SDS.nlev; ++n){
    for(Int s=0; s < shcol; ++s){
      if(SDS.rho_zt[s,n] != density){
        nerror++;
      }
    }
  }
  REQUIRE(nerror == 0);

}


}//namespace unit_test
}//namespace shoc
}//namespace scream

