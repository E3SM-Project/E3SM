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
struct UnitWrap::UnitTest<D>::TestPblintdCheckPblh {

  static void run_property()
  {
    static constexpr auto ustar_min = scream::shoc::Constants<Scalar>::ustar_min;
    static constexpr Int shcol = 5;
    static constexpr Int nlev = 4;
    static constexpr Int nlevi = nlev+1;

    // Tests for the subroutine pblintd_check_pblh

    // TEST ONE
    // Simple function and a simple test to be sure that "undefined"
    //   input values of PBLH are modified by this routine.

    // Define mid point height [m]
    static constexpr Real z[nlev] = {1500, 1000, 500, 20};
    // Define surface friction velocity [m/s]
    static constexpr Real ustar[shcol] = {0.1, 4.0, 0.9, 2.0, ustar_min};
    // Define check array
    static constexpr bool check[shcol] = {true, false, true, false, true};
    // Define PBL depth [m]
    // Make some values purposefully "undefined"
    static constexpr Real pblh[shcol] = {-20, 500, 1000, -20, 700};

    // Initialize data structure for bridging to F90
    PblintdCheckPblhData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(shcol > 0);
    REQUIRE(nlevi == nlev+1);

    // Fill in test data
    for(Int s = 0; s < shcol; ++s) {
      SDS.ustar[s] = ustar[s];
      SDS.check[s] = check[s];
      SDS.pblh[s] = pblh[s];
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        SDS.z[offset] = z[n];
      }
    }

    // check to make sure the input data makes sense
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.ustar[s] >= ustar_min);
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.z[offset] > 0);
      }
    }

    // Call the fortran implementation
    pblintd_check_pblh(SDS);

    // Check the result
    // Check that PBL height is greater than zero.  This is an
    //  important check to determine if the undefined values were modified.
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] > 0);
    }
  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    PblintdCheckPblhData f90_data[] = {
      PblintdCheckPblhData(36,  72, 73),
      PblintdCheckPblhData(72,  72, 73),
      PblintdCheckPblhData(128, 72, 73),
      PblintdCheckPblhData(256, 72, 73),
    };

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine, { {d.check, {1, 1}} });
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    PblintdCheckPblhData cxx_data[] = {
      PblintdCheckPblhData(f90_data[0]),
      PblintdCheckPblhData(f90_data[1]),
      PblintdCheckPblhData(f90_data[2]),
      PblintdCheckPblhData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      pblintd_check_pblh(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      pblintd_check_pblh_f(d.shcol, d.nlev, d.nlevi, d.nlev, d.z, d.ustar, d.check, d.pblh);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(PblintdCheckPblhData);
      for (Int i = 0; i < num_runs; ++i) {
        PblintdCheckPblhData& d_f90 = f90_data[i];
        PblintdCheckPblhData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.pblh); ++k) {
          REQUIRE(d_f90.total(d_f90.pblh) == d_cxx.total(d_cxx.pblh));
          REQUIRE(d_f90.pblh[k] == d_cxx.pblh[k]);
        }
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("pblintd_check_pblh_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdCheckPblh;

  TestStruct::run_property();
}

TEST_CASE("pblintd_check_pblh_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdCheckPblh;

  TestStruct::run_bfb();
}

} // empty namespace
