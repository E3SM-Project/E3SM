#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
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
struct UnitWrap::UnitTest<D>::TestFtermInputThirdMoms {

  static void run_property()
  {

    // Tests for the SHOC function:
    //   fterms_input_for_diag_third_shoc_moment

    // TEST
    // Given inputs, verify that output is reasonable

    // grid spacing on interface grid [m]
    constexpr static Real dz_zi = 100;
    // grid spacing on midpoint grid [m]
    constexpr static Real dz_zt = 80;
    // grid spacing on adjacent midpoint grid [m]
    constexpr static Real dz_zt_kc = 120;
    // Return to isotropic timescale [s]
    constexpr static Real isotropy_zi = 1000;
    // Brunt vaisalla frequency [s]
    constexpr static Real brunt_zi = -0.05;
    // Potential temperature on interface grid [K]
    constexpr static Real thetal_zi = 300;

    // Initialize data structure for bridging to F90
    FtermsInputForDiagThirdShocMomentData SDS;

    SDS.dz_zi = dz_zi;
    SDS.dz_zt = dz_zt;
    SDS.dz_zt_kc = dz_zt_kc;
    SDS.isotropy_zi = isotropy_zi;
    SDS.brunt_zi = brunt_zi;
    SDS.thetal_zi = thetal_zi;

    // Check that input is physical
    REQUIRE(SDS.dz_zi > 0);
    REQUIRE(SDS.dz_zt > 0);
    REQUIRE(SDS.dz_zt_kc > 0);
    REQUIRE(SDS.isotropy_zi > 0);
    REQUIRE(SDS.thetal_zi > 0);

    // Call the fortran implementation
    fterms_input_for_diag_third_shoc_moment(SDS);

    // Verify the result

    // Check that thedz2 is smaller than thedz.
    REQUIRE(SDS.thedz2 < SDS.thedz);

    // Check that bet2 is smaller than thetal
    REQUIRE(SDS.bet2 < SDS.thetal_zi);

    // Be sure that iso and isosqrd relationships hold
    if (SDS.isotropy_zi > 1){
      REQUIRE(SDS.isosqrd > SDS.iso);
    }
    else{
      REQUIRE(SDS.isosqrd < SDS.iso);
    }

  }

  static void run_bfb()
  {
    // TODO
  }

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace{

TEST_CASE("shoc_fterm_input_third_moms_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestFtermInputThirdMoms;
  
  TestStruct::run_property();
}

TEST_CASE("shoc_fterm_input_third_moms_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestFtermInputThirdMoms;

  TestStruct::run_bfb();
}

} // namespace
