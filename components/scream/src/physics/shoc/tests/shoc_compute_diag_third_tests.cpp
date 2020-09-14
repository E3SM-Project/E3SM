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
struct UnitWrap::UnitTest<D>::TestShocCompDiagThird {

  static void run_property()
  {
    static constexpr Int shcol    = 1;
    static constexpr Int nlev     = 5;
    static constexpr Int nlevi    = nlev+1;

    // Tests for the SHOC function:
    //   compute_diag_third_shoc_moment

    // Define vertical velocity second moment [m2/s2]
    static constexpr Real w_sec_zi[nlevi] = {0.2, 0.3, 0.5, 0.4, 0.3, 0.1};
    // Define potential temperature second moment [K2]
    static constexpr Real thl_sec[nlevi] = {0.5, 0.9, 1.2, 0.8, 0.4, 0.3};
    // Define second moment total water second moment [kg^2/kg^2]
    static constexpr Real qw_sec[nlevi] = {0.1e-6, 0.3e-6, 1.4e-6, 0.5e-6, 0.4e-6};
    // Define covarance of thetal and qw [K kg/kg]
    static constexpr Real qwthl_sec[nlevi] = {1e-3, 3e-3, 1.4e-3, 5e-3, 4e-3};
    // Define vertical flux of temperature [K m/s]
    static constexpr Real wthl_sec[nlevi] = {0.003, -0.03, -0.04, -0.01, 0.01, 0.03};
    // Define the heights on the zi grid [m]
//    static constexpr Real zi_grid[nlevi] = {9000, 5000, 1500, 900, 500, 0};
    static constexpr Real zi_grid[nlevi] = {900, 500, 150, 90, 50, 0};
    // Define the return to isotropy timescale [s]
    static constexpr Real isotropy_zi[nlevi] = {2000, 3000, 5000, 2000, 1000, 500};
    // Define the brunt vaisalla frequency
    static constexpr Real brunt_zi[nlevi] = {4e-5, 3e-5, 2e-5, 2e-5, -1e-5};
    // Define the potential temperature on zi grid [K]
    static constexpr Real thetal_zi[nlevi] = {330, 325, 320, 310, 300, 301};
    // Define the buoyancy flux on zi grid [K m/s]
    static constexpr Real wthv_sec_zi[nlevi] = {0.002, 0.03, 0.04, 0.01, 0.02, 0.04};
    // Define the length scale on the zi grid [m]
    static constexpr Real shoc_mix_zi[nlevi] = {600, 750, 1500, 1200, 750, 20};
    
    // Define TKE [m2/s2], compute from w_sec
    Real tke[nlev];

    // Grid stuff to compute based on zi_grid
    Real zt_grid[nlev];
    Real dz_zt[nlev];
    Real dz_zi[nlevi];
    Real w_sec[nlev];
    // Compute heights on midpoint grid
    for(Int n = 0; n < nlev; ++n) {
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
      w_sec[n] = 0.5*(w_sec_zi[n]+w_sec_zi[n+1]);
      tke[n] = 1.5*w_sec[n];
      dz_zt[n] = zi_grid[n] - zi_grid[n+1];
      if (n == 0){
        dz_zi[n] = 0;
      }
      else{
        dz_zi[n] = zt_grid[n-1] - zt_grid[n];
      }
    }
    // set upper condition for dz_zi
    dz_zi[nlevi-1] = zt_grid[nlev-1];

    // Initialize data structure for bridging to F90
    SHOCCompThirdMomData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable.
    // For this test shcol MUST be at least 2
    REQUIRE( (SDS.shcol() == shcol && SDS.nlev() == nlev && SDS.nlevi() == nlevi) );
    REQUIRE(SDS.nlevi() == SDS.nlev()+1);
    
    // Load up the new data
    for(Int s = 0; s < shcol; ++s) { 
      // Fill in test data on zt_grid.     
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
	
	SDS.w_sec[offset] = w_sec[n];
	SDS.dz_zt[offset] = dz_zt[n];
	SDS.zt_grid[offset] = zt_grid[n];
	SDS.tke[offset] = tke[n];
	
      }

      // Fill in test data on zi_grid.     
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
	
	SDS.dz_zi[offset] = dz_zi[n];
	SDS.zi_grid[offset] = zi_grid[n];
	SDS.thl_sec[offset] = thl_sec[n];
	SDS.qw_sec[offset] = qw_sec[n];
	SDS.qwthl_sec[offset] = qwthl_sec[n];
	SDS.wthl_sec[offset] = wthl_sec[n];
	
	SDS.w_sec_zi[offset] = w_sec_zi[n];
	SDS.isotropy_zi[offset] = isotropy_zi[n];
	SDS.brunt_zi[offset] = brunt_zi[n];
	SDS.thetal_zi[offset] = thetal_zi[n];
	SDS.wthv_sec_zi[offset] = wthv_sec_zi[n];
	SDS.shoc_mix_zi[offset] = shoc_mix_zi[n];
	
      }
    }      

    // Check that the inputs make sense
    // Load up the new data
    for(Int s = 0; s < shcol; ++s) { 
      // Fill in test data on zt_grid.    
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
	
        REQUIRE(SDS.w_sec[offset] >= 0);
	REQUIRE(SDS.dz_zt[offset] > 0);
	REQUIRE(SDS.zt_grid[offset] > 0);
//	REQUIRE(SDS.tke[offset] > 0);  
      }
      
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
	
	REQUIRE(SDS.dz_zi[offset] >= 0);
	REQUIRE(SDS.zi_grid[offset] >= 0);
	REQUIRE(SDS.thl_sec[offset] >= 0);
	REQUIRE(SDS.qw_sec[offset] >= 0);
	REQUIRE(SDS.w_sec_zi[offset] >= 0);
	REQUIRE(SDS.isotropy_zi[offset] >= 0);
	REQUIRE(SDS.thetal_zi[offset] >= 0);
	REQUIRE(SDS.shoc_mix_zi[offset] >= 0);
	
      }    
      
    }

    // Call the fortran implementation
    compute_diag_third_shoc_moment(SDS);

  }

  static void run_bfb()
  {
    // TODO
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_comp_diag_third_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocCompDiagThird;

  TestStruct::run_property();
}

TEST_CASE("shoc_comp_diag_third_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocCompDiagThird;

  TestStruct::run_bfb();
}

} // namespace
