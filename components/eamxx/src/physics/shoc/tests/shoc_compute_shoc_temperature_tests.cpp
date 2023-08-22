#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_setup_random_test.hpp"

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
struct UnitWrap::UnitTest<D>::TestComputeShocTemp {

  static void run_property()
  {
    // TODO: Add property test?
  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ComputeShocTempData f90_data[] = {
      //            shcol, nlev
      ComputeShocTempData(10, 71),
      ComputeShocTempData(10, 12),
      ComputeShocTempData(7,  16),
      ComputeShocTempData(2,   7)
    };

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ComputeShocTempData cxx_data[] = {
      ComputeShocTempData(f90_data[0]),
      ComputeShocTempData(f90_data[1]),
      ComputeShocTempData(f90_data[2]),
      ComputeShocTempData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      compute_shoc_temperature(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      compute_shoc_temperature_f(d.shcol, d.nlev, d.thetal, d.ql, d.inv_exner, d.tabs);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(ComputeShocTempData);
      for (Int i = 0; i < num_runs; ++i) {
        ComputeShocTempData& d_f90 = f90_data[i];
        ComputeShocTempData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.tabs); ++k) {
          REQUIRE(d_f90.tabs[k] == d_cxx.tabs[k]);
        }
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

// TEST_CASE("shoc_imp_dp_inverse_property", "shoc")
// {
//   using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeShocTemp;

//   TestStruct::run_property();
// }

TEST_CASE("shoc_compute_shoc_temperature_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeShocTemp;

  TestStruct::run_bfb();
}

} // namespace
