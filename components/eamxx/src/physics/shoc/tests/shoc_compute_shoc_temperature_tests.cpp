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
    static constexpr Int shcol    = 1;
    static constexpr Int nlev     = 3;

    // Test One
    // Given Exner value = 1 and cloud liquid = 0 everywhere
    //  verify that all points of absolute temperature (tabs)
    //  are exactly equal to liquid water potential temperature (thetal)

    // Inverse Exner value [-]
    Real inv_exner_first[nlev] = {1, 1, 1};
    // Liquid water potential temperature [K]
    Real thetal_first[nlev] = {300, 290, 280};
    // Liquid water mixing ratio [kg/kg]
    Real ql_first[nlev] = {0, 0, 0};

    ComputeShocTempData SDS(shcol, nlev);

    REQUIRE(SDS.shcol > 0);
    REQUIRE(SDS.nlev > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.inv_exner[offset] = inv_exner_first[n];
        SDS.thetal[offset] = thetal_first[n];
        SDS.ql[offset] = ql_first[n];
      }
    }

    // Check that inputs are as expected
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.ql[offset] == 0);
        REQUIRE(SDS.inv_exner[offset] == 1);
        REQUIRE(SDS.thetal[offset] > 0);
      }
    }

    // Call the fortran implementation
    compute_shoc_temperature(SDS);

    // Require that absolute temperature is equal to thetal
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.tabs[offset] == SDS.thetal[offset]);
      }
    }

    // Test Two
    // Given profiles all with cloud liquid water greater than zero,
    //  AND inverse exner functions equal to 1, ensure that tabs is greater that
    //  absolute temperature is greater than (or equal to) liquid water potential temperature

    // Inverse Exner value [-]
    Real inv_exner_sec[nlev] = {1, 1, 1};
    // Liquid water potential temperature [K]
    Real thetal_sec[nlev] = {300, 290, 280};
    // Liquid water mixing ratio [kg/kg]
    Real ql_sec[nlev] = {3e-5, 1e-10, 1e-3};

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.inv_exner[offset] = inv_exner_sec[n];
        SDS.thetal[offset] = thetal_sec[n];
        SDS.ql[offset] = ql_sec[n];
      }
    }

    // Check that inputs are as expected
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.ql[offset] > 0);
        REQUIRE(SDS.inv_exner[offset] == 1);
        REQUIRE(SDS.thetal[offset] > 0);
      }
    }

    // Call the fortran implementation
    compute_shoc_temperature(SDS);

    // Require that absolute temperature is greather than thetal
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.tabs[offset] >= SDS.thetal[offset]);
      }
    }

    // Test Three
    // Given "realistic" atmospheric profiles with thetal increasing
    //  with height, as well as inv_exner function; verify that
    //  absolute temperature is decreasing with height.  Give
    //  all levels the same amount of cloud liquid water

    // Inverse Exner value [-]
    Real inv_exner_third[nlev] = {1.1, 1.5, 2};
    // Liquid water potential temperature [K]
    Real thetal_third[nlev] = {300, 350, 400};
    // Liquid water mixing ratio [kg/kg]
    Real ql_third[nlev] = {1e-5, 1e-5, 1e-5};

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.inv_exner[offset] = inv_exner_third[n];
        SDS.thetal[offset] = thetal_third[n];
        SDS.ql[offset] = ql_third[n];
      }
    }

    // Check that inputs are as expected
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.ql[offset] > 0);
        REQUIRE(SDS.inv_exner[offset] > 1);
        REQUIRE(SDS.thetal[offset] > 0);
      }
    }

    // Check that inputs are changing with height as expected
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 1; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        const auto offsetl = (n-1) + s * nlev;

        // Verify inverse exner and thetal are increasing with height
        REQUIRE(SDS.inv_exner[offset] > SDS.inv_exner[offsetl]);
        REQUIRE(SDS.thetal[offset] > SDS.thetal[offsetl]);
      }
    }

    // Call the fortran implementation
    compute_shoc_temperature(SDS);

    // Require that absolute temperature be less than thetal
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.tabs[offset] < SDS.thetal[offset]);
      }
    }

    // Check that tabs is decreasing with height as expected
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 1; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        const auto offsetl = (n-1) + s * nlev;

        REQUIRE(SDS.tabs[offset] < SDS.tabs[offsetl]);
      }
    }

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

TEST_CASE("shoc_compute_shoc_temperature_property", "shoc")
 {
   using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeShocTemp;

   TestStruct::run_property();
 }

TEST_CASE("shoc_compute_shoc_temperature_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeShocTemp;

  TestStruct::run_bfb();
}

} // namespace
