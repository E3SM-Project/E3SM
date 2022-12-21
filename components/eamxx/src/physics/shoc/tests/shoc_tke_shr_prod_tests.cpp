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
struct UnitWrap::UnitTest<D>::TestShocShearProd {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the subroutine compute_shr_prod in the SHOC
    //   TKE module.

    // FIRST TEST
    //  For first tests input a sheared profile for both wind
    //  components, one with zonal winds increasing with height
    //  at a constant rate per GRID box and another with meridional
    //  winds decreasing at a constant rate per GRID box. The
    //  grid prescribed will be a stretched grid. Here we want to
    //  validate that the shear term DECREASES with height.
    //  NOTE: sterm boundaries will be returned as ZERO.

    // Define height thickness on nlevi grid [m]
    //   NOTE: First indicee is zero because it is never used
    //   Do a stretched grid
    static constexpr Real dz_zi[nlevi] = {0, 500, 200, 100, 50, 10};
    // Define zonal wind on nlev grid [m/s]
    static constexpr Real u_wind_shr[nlev] = {2, 1, 0, -1, -2};
    // Define meridional wind on nlev grid [m/s]
    static constexpr Real v_wind_shr[nlev] = {1, 2, 3, 4, 5};

    // Initialize data structure for bridging to F90
    ComputeShrProdData SDS(shcol, nlev, nlevi);

    // Define upper limit for reasonable bounds checking
    static constexpr Real sterm_upper_bound = 1e-2;

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(nlevi - nlev == 1);
    REQUIRE(shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.u_wind[offset] = u_wind_shr[n];
        SDS.v_wind[offset] = v_wind_shr[n];
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset   = n + s * nlevi;
        SDS.dz_zi[offset] = dz_zi[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlevi; ++n){
        const auto offset = n + s * nlevi;
        // Make sure top level dz_zi value is zero
        if (n == 0){
          REQUIRE(SDS.dz_zi[offset] == 0);
        }
        // Otherwise, should be greater than zero
        else{
          REQUIRE(SDS.dz_zi[offset] > 0);
        }
      }
    }

    // Call the fortran implementation
    compute_shr_prod(SDS);

    // Check test
    for(Int s = 0; s < shcol; ++s) {
      // First check that sterm is ALWAYS greater than
      //  zero for non boundary points, but exactly zero
      //  for boundary points.
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
        // Make sure output falls within reasonable bound
        REQUIRE(SDS.sterm[offset] < sterm_upper_bound);
        if (n == 0 || n == nlevi-1){
          // Boundary point check
          REQUIRE(SDS.sterm[offset] == 0);
        }
        else{
          REQUIRE(SDS.sterm[offset] > 0);
        }
      }
      // Now validate that shear term is ALWAYS
      //  decreasing with height for these inputs, keeping
      //  in mind to exclude boundary points, which should be zero
      for(Int n = 1; n < nlevi-2; ++n){
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.sterm[offset]-SDS.sterm[offset+1] < 0);
      }
    }

    // SECOND TEST
    // For second test we input wind profiles that are
    // constant with height to validate that shear production
    // term is zero everywhere.

    // Define zonal wind on nlev grid [m/s]
    static constexpr Real u_wind_cons = 10;
    // Define meridional wind on nlev grid [m/s]
    static constexpr Real v_wind_cons = -5;

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.u_wind[offset] = u_wind_cons;
        SDS.v_wind[offset] = v_wind_cons;
      }
    }

    // Call the fortran implementation
    compute_shr_prod(SDS);

    // Check test
    // Verify that shear term is zero everywhere
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.sterm[offset] == 0);
      }
    }
  }

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ComputeShrProdData f90_data[] = {
      //            shcol, nlev
      ComputeShrProdData(10, 71, 72),
      ComputeShrProdData(10, 12, 13),
      ComputeShrProdData(7,  16, 17),
      ComputeShrProdData(2,   7, 8)
    };

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ComputeShrProdData cxx_data[] = {
      ComputeShrProdData(f90_data[0]),
      ComputeShrProdData(f90_data[1]),
      ComputeShrProdData(f90_data[2]),
      ComputeShrProdData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      compute_shr_prod(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      compute_shr_prod_f(d.nlevi, d.nlev, d.shcol, d.dz_zi, d.u_wind, d.v_wind, d.sterm);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(ComputeShrProdData);
      for (Int i = 0; i < num_runs; ++i) {
        ComputeShrProdData& d_f90 = f90_data[i];
        ComputeShrProdData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.sterm); ++k) {
          REQUIRE(d_f90.sterm[k] == d_cxx.sterm[k]);
        }
      }
    }
  } //run_bfb
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_tke_shr_prod_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocShearProd;

  TestStruct::run_property();
}

TEST_CASE("shoc_tke_shr_prod_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocShearProd;

  TestStruct::run_bfb();
}

} // namespace
