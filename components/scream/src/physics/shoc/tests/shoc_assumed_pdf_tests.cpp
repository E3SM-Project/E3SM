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
struct UnitWrap::UnitTest<D>::TestShocAssumedPdf {

  static void run_property()
  {
    static constexpr Int shcol    = 1;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the top level subroutine
    //   shoc_assumed_pdf

    // Define input data
    // Liquid water potential temperature [K]
    static constexpr Real thetal[nlev] = {320, 315, 310, 305, 300};
    // Total water mixing ratio [kg/kg]
    static constexpr Real qw[nlev] = {0.010, 0.011, 0.012, 0.014, 0.020};
    // Pressure [Pa]
    static constexpr Real pres[nlev] = {85000, 87500, 90000, 95000, 100000};
    // Define the heights on the zt grid [m]
    static constexpr Real zi_grid[nlevi] = {2000, 1500, 1000, 500, 150, 0};    
    // Variance of thetal [K^2]
//    static constexpr Real thl_sec[nlevi] = {0, 0, 0, 0, 0, 0};
    // Variance of total water mixing ratio [kg^2/kg^2]
//    static constexpr Real 

    Real zt_grid[nlev];
    // Compute heights on midpoint grid
    for(Int n = 0; n < nlev; ++n) {
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
    }

    // Initialize data structure for bridging to F90
    SHOCAssumedpdfData SDS(shcol, nlev, nlevi);
    
    // Load input data
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
	const auto offset = n + s * nlev;
       
        SDS.thetal[offset] = thetal[n];
        SDS.qw[offset] = qw[n];
        SDS.pres[offset] = pres[n];
        SDS.zt_grid[offset] = zt_grid[n];
        SDS.w_field[offset] = 0;
        
      }
      
      for(Int n = 0; n < nlevi; ++n) {
	const auto offset = n + s * nlevi;
       
        SDS.thl_sec[offset] = 0;
        SDS.qw_sec[offset] = 0;
        SDS.wthl_sec[offset] = 0;
        SDS.wqw_sec[offset] = 0;
        SDS.qwthl_sec[offset] = 0;
        SDS.w_sec[offset] = 0.0004;
        SDS.w3[offset] = 0;
        SDS.zi_grid[offset] = zi_grid[n];
        
      }          
    }    

    // Test that the inputs are reasonable.
    REQUIRE(SDS.nlevi() - SDS.nlev() == 1);
    REQUIRE(SDS.shcol() > 0);

    // Call the fortran implementation
    shoc_assumed_pdf(SDS);
       
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

TEST_CASE("shoc_assumed_pdf_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocAssumedPdf;

  TestStruct::run_property();
}

TEST_CASE("shoc_assumed_pdf_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocAssumedPdf;

  TestStruct::run_bfb();

}

} // namespace
