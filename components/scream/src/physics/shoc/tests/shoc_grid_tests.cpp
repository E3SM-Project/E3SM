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

  const Int shcol = 2; 
  const Int nlev = 128; 
  const auto nlevi = nlev + 1;
  const Real density = 1.0; 


  SHOCGridData SDS(shcol, nlev, nlevi);

  // Test that the inputs are resonable.
  REQUIRE(SDS.nlevi - SDS.nlev == 1); 
  REQUIRE(SDS.shcol>0);

  const Real dz = 50.0;
  // Fill in test data on zt_grid
  for(Int s=0; s < SDS.shcol; ++s){
    for(Int n=0; n < SDS.nlev; ++n){
      const auto nft = SDS.nlev - 1 - n;
      const auto offset = n + s*SDS.nlev; 
      
      SDS.zt_grid[offset] = nft * dz + dz/2; 
      SDS.pdel[offset] = density*gravit*dz; 
    }
  }

  // Fill in test data on zi_grid
  for(Int s=0; s < SDS.shcol; ++s){
    for(Int n=0; n < SDS.nlevi; ++n){
       const auto nft = SDS.nlevi - 1 - n;
       const auto offset = n + s*SDS.nlevi;

       SDS.zi_grid[offset] =  dz * nft;
      }
  }

  // Here we sould check that the inputs make sense whe

 // Check that zt decreases upward
 int nerror = 0; 
 for(Int s=0; s < SDS.shcol; ++s){
   for(Int n=0; n < SDS.nlev-1; ++n){    
     const auto offset = n + s*SDS.nlev;
     if (SDS.zt_grid[offset+1] - SDS.zt_grid[offset] > 0.0){
          nerror++;
          }
      }
  }  
  REQUIRE(nerror == 0);
 
 //Check that zi decreases upward
 for(Int s=0; s < SDS.shcol; ++s){
   for(Int n=0; n < SDS.nlevi-1; ++n){
     const auto offset = n + s*SDS.nlevi;
     if (SDS.zi_grid[offset+1] - SDS.zi_grid[offset] > 0.0){
          nerror++;
      }
    }
  }  
  REQUIRE(nerror == 0);

  //Call the fortran implementation
  shoc_grid(nlev, SDS);

  //First check that dz is correct
  for(Int s=0; s < shcol; ++s){
    for(Int n=0; n < SDS.nlev; ++n){
     const auto offset = n + s*SDS.nlev;
      if (SDS.dz_zt[offset] != dz){
        nerror++;
      }
    }
  } 
  REQUIRE(nerror == 0);
  
  for(Int s=0; s < shcol; ++s){
    const auto offset = s*SDS.nlevi;
    REQUIRE(SDS.dz_zi[offset] == 0);
    REQUIRE(SDS.dz_zi[offset + SDS.nlevi-1] == dz/2.0);
  }

  for(Int s=0; s < shcol; ++s){
    for(Int n=1; n < SDS.nlevi-1; ++n){
      const auto offset = n + s*SDS.nlevi;
      if (SDS.dz_zi[offset] != dz){
        nerror++;
      }
    }
  }
  REQUIRE(nerror == 0);
  
  //Now check density 
  for(Int s=0; s < shcol; ++s){
    for(Int n=0; n < SDS.nlev; ++n){
      const auto offset = n + s*SDS.nlev;
      if(SDS.rho_zt[offset] != density){
        nerror++;
      }
    }
  }
  REQUIRE(nerror == 0);

}

}//namespace unit_test
}//namespace shoc
}//namespace scream

