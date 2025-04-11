#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/eamxx_types.hpp"

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
struct UnitWrap::UnitTest<D>::TestShocPdfComputeTemp {

  static void run_property()
  {
    // Property tests for the SHOC function
    //  shoc_assumed_pdf_compute_temperature

    // TEST ONE
    // Keep liquid water potential temperature constant and
    //  decrease pressure by 10000 Pa each iteration.  Verify
    //  that the liquid water temperature behaves as expected
    //  (i.e. temperature decreases with decreasing pressure and
    //  value is as expected relative to base pressure).

    // Input liquid water potential temperature [K]
    static constexpr Real thl1 = 305;
    // Input basepressure [Pa]
    static constexpr Real basepres = C::P0;
    // Input value of pval [Pa]
    Real pval = 110000;

    // Decrease base pres by a certain amount each test
    static constexpr Real presincr = -10000;

    // Define reasonable bounds for checking
    static constexpr Real Tl_lower_bound = 150;
    static constexpr Real Tl_upper_bound = 350;

    Real Tl1_save;

    // Initialize data structure for bridging to F90
    ShocAssumedPdfComputeTemperatureData SDS;

    // Fill in data
    SDS.thl1 = thl1;
    SDS.pval = pval;

    Int num_tests = SDS.pval/abs(presincr);

    REQUIRE(num_tests > 1);
    REQUIRE(presincr < 0);
    // Make sure our starting pressure is greater than
    //  basepres just so we test a range
    REQUIRE(SDS.pval > basepres);

    for (Int s = 0; s < num_tests; ++s){

      // make sure base pres is greater than zero
      REQUIRE(basepres > 0);

      // Call the fortran implementation
      shoc_assumed_pdf_compute_temperature(SDS);

      // Check the result
      // Make sure temperature falls within reasonable bound
      REQUIRE(SDS.tl1 < Tl_upper_bound);
      REQUIRE(SDS.tl1 > Tl_lower_bound);

      // If pressure is greater than basepressure then
      //  make sure that temperature is greater than thetal
      if (SDS.pval > basepres){
        REQUIRE(SDS.tl1 > SDS.thl1);
      }
      // otherwise temperature should be less than thetal
      else if(SDS.pval < basepres){
        REQUIRE(SDS.tl1 < SDS.thl1);
      }
      // otherwise if they are equal the temperatures
      //  should be equal
      else
      {
        REQUIRE(SDS.tl1 == SDS.thl1);
      }

      // Make sure temperature are decreasing
      if (s > 0){
        REQUIRE(SDS.tl1 < Tl1_save);
      }

      // Save result to make sure that temperatures
      //  are decreasing as pressure decreases
      Tl1_save = SDS.tl1;

      // Decrease pressure value
      SDS.pval = SDS.pval+presincr;

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

namespace {

TEST_CASE("shoc_pdf_computetemp_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeTemp;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_computetemp_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfComputeTemp;

  TestStruct::run_bfb();
}

} // namespace
