#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestPblintdInitPot {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 1;  
  
    // Tests for the SHOC function:
    //   pblintd_init_pot
    
    // FIRST TEST
    //  Dry atmosphere test.  Verify that in a dry atmosphere
    //  with no water vapor or liquid water loading that the output
    //  thv (virtual potential temperature) is equal to the input
    //  thl (liquid water potential temperature).  

    // Define the liquid water potential temperature [K]
    Real thl_dry[shcol] = {300, 295};
    // Define the water vapor [kg/kg]
    Real qv = 0.0;
    // Define the liquid water mixing ratio [kg/kg]
    Real ql = 0.0;

    // Initialize data structure for bridging to F90
    SHOCPblintdInitPotData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol() == shcol);
    REQUIRE(SDS.nlev() == nlev);
    // for this test require two columns
    REQUIRE(shcol == 2);
    
    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.thl[offset] = thl_dry[s];
        SDS.ql[offset] = ql;
        SDS.q[offset] = qv;
      }
    }   
    
     // Check that the inputs are expected
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        REQUIRE(SDS.thl[offset] > 0);
        // For this test require ql and q be zero
        REQUIRE(SDS.q[offset] == 0);
        REQUIRE(SDS.ql[offset] == 0);
      }
    }
    
    // call the fortran implementation
    shoc_pblintd_init_pot(SDS);
    
    // Check the result.  
    // Verify that virtual potential temperature is idential 
    //  to potential temperatuer
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        REQUIRE(SDS.thl[offset] == SDS.thv[offset]);
      }
    }
        
  }

  static void run_bfb()
  {
    SHOCPblintdInitPotData pblintd_init_pot_data_f90[] = {
      //                     shcol, nlev
      SHOCPblintdInitPotData(36,  72),
      SHOCPblintdInitPotData(72,  72),
      SHOCPblintdInitPotData(128, 72),
      SHOCPblintdInitPotData(256, 72),
    };

    static constexpr Int num_runs = sizeof(pblintd_init_pot_data_f90) / sizeof(SHOCPblintdInitPotData);

    for (auto& d : pblintd_init_pot_data_f90) {
      d.randomize();
    }

    SHOCPblintdInitPotData pblintd_init_pot_data_cxx[] = {
      SHOCPblintdInitPotData(pblintd_init_pot_data_f90[0]),
      SHOCPblintdInitPotData(pblintd_init_pot_data_f90[1]),
      SHOCPblintdInitPotData(pblintd_init_pot_data_f90[2]),
      SHOCPblintdInitPotData(pblintd_init_pot_data_f90[3]),
    };

    for (auto& d : pblintd_init_pot_data_f90) {
      // expects data in C layout
      shoc_pblintd_init_pot(d);
    }

    for (auto& d : pblintd_init_pot_data_cxx) {
      shoc_pblintd_init_pot_f(d.shcol(), d.nlev(), d.thl, d.ql, d.q, d.thv);
    }

    for (Int i = 0; i < num_runs; ++i) {
      Int shcol = pblintd_init_pot_data_cxx[i].shcol();
      Int nlev  = pblintd_init_pot_data_cxx[i].nlev();
      for (Int j = 0; j < shcol; ++j ) {
        for (Int k = 0; k < nlev; ++k) {
          REQUIRE(pblintd_init_pot_data_f90[i].thv[j*k] == pblintd_init_pot_data_cxx[i].thv[j*k]);
        }
      }
    }
  }

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("shoc_pblintd_init_pot_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdInitPot;

  TestStruct::run_property();
}

TEST_CASE("shoc_pblintd_init_pot_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdInitPot;

  TestStruct::run_bfb();
}

}  // namespace
