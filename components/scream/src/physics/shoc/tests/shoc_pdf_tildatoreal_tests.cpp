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
struct UnitWrap::UnitTest<D>::TestShocPdfTildatoReal {

  static void run_property()
  {

    // Define the midpoint height grid [m]
    static constexpr Real w_first = 1;
    static constexpr Real sqrtw2 = 0.5;
    Real w1 = 0.1;
    
    SHOCPDFtildaData SDS;
    
    SDS.w_first = w_first;
    SDS.sqrtw2 = sqrtw2;
    SDS.w1 = w1;

    shoc_assumed_pdf_tilda_to_real(SDS);

    // Check the test here

  }
  
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_pdf_tildatoreal_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfTildatoReal;

  TestStruct::run_property();
}

TEST_CASE("shoc_pdf_tildatoreal_b4b", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocPdfTildatoReal;

}

} // namespace
