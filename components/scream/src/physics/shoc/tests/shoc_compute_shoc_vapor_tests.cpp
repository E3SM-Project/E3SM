#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestComputeShocVapor {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;

    // Tests for the SHOC subroutine:
    //  compute_shoc_vapor

    // Test
    // Given profile of a variety of conditions,
    //  verify that the output is as expected

    // Total water mixing ratio [kg/kg]
    static constexpr Real qw[nlev] = {1e-2, 1.2e-2, 1.5e-2, 1.7e-2, 2.0e-2};
    // Liquid water mixing ratio [kg/kg]
    static constexpr Real ql[nlev] = {0, 0, 1.5e-4, 2e-3, 0};

    // Initialize data structure for bridging to F90
    ComputeShocVaporData SDS(shcol, nlev);

    // Load input data
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.qw[offset] = qw[n];
        SDS.ql[offset] = ql[n];

      }
    }

    // Check that inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // total water should be greater than zero
        REQUIRE(SDS.qw[offset] > 0);
        // cloud water greater than or equal to zero
        REQUIRE(SDS.ql[offset] >= 0);
        // total water should be greater than cloud water
        REQUIRE(SDS.qw[offset] > SDS.ql[offset]);

      }
    }

    // Call the fortran implementation
    compute_shoc_vapor(SDS);

    // Verify the result
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // vapor should be greater than zero
        REQUIRE(SDS.qv[offset] > 0);
        // if cloud present vapor should be less than total water
        if (SDS.ql[offset] > 0){
          REQUIRE(SDS.qv[offset] < SDS.qw[offset]);
        }
        // else they should be equal
        else{
          REQUIRE(SDS.qv[offset] == SDS.qw[offset]);
        }

      }
    }

  } // run_property

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ComputeShocVaporData f90_data[] = {
      //              shcol, nlev
      ComputeShocVaporData(10, 71),
      ComputeShocVaporData(10, 12),
      ComputeShocVaporData(7,  16),
      ComputeShocVaporData(2,   7),
    };

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ComputeShocVaporData cxx_data[] = {
      ComputeShocVaporData(f90_data[0]),
      ComputeShocVaporData(f90_data[1]),
      ComputeShocVaporData(f90_data[2]),
      ComputeShocVaporData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      compute_shoc_vapor(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      compute_shoc_vapor_f(d.shcol, d.nlev, d.qw, d.ql, d.qv);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(ComputeShocVaporData);
      for (Int i = 0; i < num_runs; ++i) {
        ComputeShocVaporData& d_f90 = f90_data[i];
        ComputeShocVaporData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.qv); ++k) {
          REQUIRE(d_f90.qv[k] == d_cxx.qv[k]);
        }
      }
    }
  } // run_bfb
};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("compute_shoc_vapor_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeShocVapor;

  TestStruct::run_property();
}

TEST_CASE("compute_shoc_vapor_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeShocVapor;

  TestStruct::run_bfb();
}

} // empty namespace
