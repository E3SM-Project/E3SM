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
#include "physics/share/physics_constants.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestShocPdfCompBuoyFlux {

  static void run_property()
  {
  
    static constexpr Real epsterm  = scream::physics::Constants<Real>::ep_2;
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_compute_buoyancy_flux

    // TEST ONE
    // Dry atmostphere test.  If the atmosphere has no moisture, then 
    //  verify that wthvsec = wthlsec
    
    // Input data
    // Define the heat flux [K m/s]
    static constexpr Real wthlsec = 0.02;
    // Define the latent heat flux [kg/kg m/s]
    static constexpr Real wqwsec_dry = 0;
    // Define the liquid water flux [kg/kg m/s]
    static constexpr Real wqls_dry = 0;
    // Define the input pressure [Pa]
    static constexpr Real pval = 85000;
    
    // Initialize data structure for bridging to F90
    SHOCPDFcompbuoyfluxData SDS;
    
    // Load input data
    SDS.wthlsec = wthlsec;
    SDS.wqwsec = wqwsec_dry;
    SDS.wqls = wqls_dry;
    SDS.pval = pval;
    SDS.epsterm = epsterm;
    
    // Call the fortran implementation
    shoc_assumed_pdf_compute_buoyancy_flux(SDS);
    
    // Check the result
    // Verify that buoyancy flux is equal to heat flux
    REQUIRE(SDS.wthv_sec == SDS.wthlsec);
    
    // TEST TWO
    // Positive in cloud buoyancy test
    // Given two tests and inputs that are identical but one test
    //  with positive liquid water flux and the other with zero liquid
    //  water flux, verify that buoyancy flux is greater when cloud is
    //  present and in updrafts.  
    
    // Some inputs will be recycled from previous test
    
    // Define the latent heat flux [kg/kg m/s]
    static constexpr Real wqwsec_cloud = 0.004;
    
    // Load the data
    SDS.wqwsec = wqwsec_cloud;
    
    // For this test verify that liquid water flux is zero
    REQUIRE(SDS.wqls == 0);
    
    // Call the fortran implementation
    shoc_assumed_pdf_compute_buoyancy_flux(SDS);
    
    // Save the result
    Real wthv_sec_nocloud = SDS.wthv_sec;
    
    // Define the liquid water flux [kg/kg m/s]
    static constexpr Real wqls_cloud = 0.0004;    
    
    SDS.wqls = wqls_cloud;
    
    // For this test be sure that liquid water flux is positive
    REQUIRE(SDS.wqls > 0);
    
    // Call the fortran implementation
    shoc_assumed_pdf_compute_buoyancy_flux(SDS);
    
    // Make sure that this buoyancy flux is greater than the buoyancy
    //   flux with no incloud updraft 
    REQUIRE(SDS.wthv_sec > wthv_sec_nocloud);
    
    // Also do some consistency checks.  Such as buoyancy flux
    //  should be greater (in this case) than the heat flux and
    //  the liquid water flux
    REQUIRE(SDS.wthv_sec > SDS.wthlsec);
    REQUIRE(SDS.wthv_sec > SDS.wqls);

  }
  
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_pdf_compute_buoyflux_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfCompBuoyFlux;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_compute_buoyflux_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfCompBuoyFlux;

}

} // namespace
