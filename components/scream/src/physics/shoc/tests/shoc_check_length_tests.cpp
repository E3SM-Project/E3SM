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
struct UnitWrap::UnitTest<D>::TestCheckShocLength {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC function:
    //   check_length_scale_shoc_length
  
    // Test this simple function.  Feed the routine values of
    //  mixing length that are zero and much larger than the
    //  host grid box size to make sure that the function
    //  corrects these erroneous values.  
  
    // Define the host grid box size x-direction [m]
    static constexpr Real host_dx = 3000.0; 
    // Defin the host grid box size y-direction [m]
    static constexpr Real host_dy = 5000.0; 
    // Define mixing length [m]
    //   In profile include values of zero and very large values
    Real shoc_mix[nlev] = {50000.0, 4000.0, 2000.0, 0.0, 500.0};

    // Initialize data structure for bridgeing to F90
    SHOCMixcheckData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol > 0);
  
    // compute geometric grid mesh
    const auto grid_mesh = sqrt(host_dx*host_dy);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < SDS.shcol; ++s) {
      SDS.host_dx[s] = host_dx;
      SDS.host_dy[s] = host_dy;
      for(Int n = 0; n < SDS.nlev; ++n) {
        const auto offset = n + s * SDS.nlev;

        SDS.shoc_mix[offset] = shoc_mix[n];
      }
    }

    // Check that the inputs make sense

    // Be sure that relevant variables are greater than zero
    for(Int s = 0; s < SDS.shcol; ++s) {
      REQUIRE(SDS.host_dx[s] > 0.0);
      REQUIRE(SDS.host_dy[s] > 0.0);
    }

    // Call the fortran implementation
    check_length_scale_shoc_length(SDS);

    // Check the results
    for(Int s = 0; s < SDS.shcol; ++s) {   
      for(Int n = 0; n < SDS.nlev; ++n) {
        const auto offset = n + s * SDS.nlev;   
        // Require mixing length is greater than zero and is
        //  less than geometric grid mesh length + 1 m
        REQUIRE(SDS.shoc_mix[offset] > 0.0); 
        REQUIRE(SDS.shoc_mix[offset] < 1.0+grid_mesh); 
      }   
    }
  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace{

TEST_CASE("shoc_check_length_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCheckShocLength;
  
  TestStruct::run_property();
}

TEST_CASE("shoc_check_length_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCheckShocLength;

}

} // namespace
