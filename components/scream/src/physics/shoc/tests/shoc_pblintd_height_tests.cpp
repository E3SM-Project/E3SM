#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestPblintdHeight {

  static void run_property()
  {
    const auto ustar_min = scream::shoc::Constants<Scalar>::ustar_min;
    static const auto approx_zero = Approx(0.0).margin(1e-16);
    static constexpr Int shcol = 4;
    static constexpr Int nlev = 9;

    // Tests for the subroutine pblintd_height
    // Perform a series of tests to ensure that subroutine returns values
    //  as expected

    // TEST ONE
    // Given input where the value of ustar is increasing, verify that
    //  the PBL depth also increases.  Input temperature should be stable

    // Define the midpoint heights [m]
    static constexpr Real z[nlev] = {5000, 4000, 3000, 2000, 1000, 750, 500, 250, 100};
    // Define the u wind [m/s]
    static constexpr Real u[nlev] = {10, 10, 10, 9, 9, 8, 7, 6, 5};
    // Define the v wind [m/s]
    static constexpr Real v[nlev] = {-10, -10, -10, -9, -9, -8, -7, -6, -5};
    // Define the virtual potential temperature [K]
    static constexpr Real thv[nlev] = {320, 315, 314, 312, 311, 310, 302, 302, 300};
    // Define the surface friction velocity [m4/s3]
    Real ustar = ustar_min;

    // Initialize dtata structure for bridging to F90
    PblintdHeightData SDS(shcol,nlev);

    // Require that the inputs are reasonable
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    // This test requires more than one column
    REQUIRE(shcol > 1);

    // Fill in test data on the zt_grid
    for(Int s = 0; s < shcol; ++s) {
      SDS.ustar[s] = ustar + s; // Add to ustar with each column
      SDS.check[s] = true;
      SDS.thv_ref[s] = thv[nlev-1];
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        SDS.z[offset] = z[n];
        SDS.u[offset] = u[n];
        SDS.v[offset] = v[n];
        SDS.thv[offset] = thv[n];
        // Should be initialized to zero everywhere
        SDS.rino[offset] = 0;
      }
    }

    // check to make sure the input data makes sense
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.ustar[s] >= ustar_min);
      REQUIRE(SDS.check[s] == true);
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.z[offset] > 0);
        REQUIRE(std::abs(SDS.u[offset] < 100));
        REQUIRE(std::abs(SDS.v[offset] < 100));
        REQUIRE(SDS.thv[offset] < 1000);
        REQUIRE(SDS.thv[offset] > 150);
      }
    }

    // check that as column increases then does ustar
    for(Int s = 0; s < shcol-1; ++s) {
      REQUIRE(SDS.ustar[s+1] > SDS.ustar[s]);
    }

    // Call the fortran implementation
    pblintd_height(SDS);

    // Check the result
    for(Int s = 0; s < shcol; ++s) {
      // Verify PBL height falls within reasonable bounds
      REQUIRE(SDS.pblh[s] > z[nlev-1]); // should be higher than lowest grid
      REQUIRE(SDS.pblh[s] < z[0]); // and lower than highest grid
      REQUIRE(SDS.check[s] == false);
    }

    for(Int s = 0; s < shcol-1; ++s) {
      // Verify that as ustar increases, PBLH increases
      REQUIRE(SDS.pblh[s+1] > SDS.pblh[s]);
    }

    // Save result for checking with next test
    Real pblh_save[shcol];
    for(Int s = 0; s < shcol-1; ++s) {
      pblh_save[s] = SDS.pblh[s];
    }

    // TEST TWO
    // Now use the same inputs as the last test BUT modify the surface
    //  temperature value so that it is unstable.  The PBLH should be
    //  larger everywhere.

    // Replace lowest input temperature to a value that is unstable
    for(Int s = 0; s < shcol; ++s) {
      SDS.check[s] = true; // reinitialize
      const auto offset = nlev-1 + s * nlev;
      SDS.thv[offset] = SDS.thv[offset-1] + 2;
    }

    // Call the fortran implementation
    pblintd_height(SDS);

    // Check the result
    for(Int s = 0; s < shcol; ++s) {
      // Verify PBL height falls within reasonable bounds
      REQUIRE(SDS.pblh[s] > z[nlev-1]); // should be higher than lowest grid
      REQUIRE(SDS.pblh[s] < z[0]); // and lower than highest grid
      REQUIRE(SDS.check[s] == false);
    }

    // Verify that the PBLh with unstable atmosphere is LARGER everywhere
    //  than the tests with conditionally stable atmosphere
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] > pblh_save[s]);
    }

    // TEST THREE
    // Given profiles where all variables are homogenous, verify that
    //  PBLH everywhere is returned as zero (in real code the subroutine
    //  will be called again to ammeliorate this, tested in upper level function
    //  for the PBLH).

    // Fill in test data on the zt_grid, all other data recycled
    for(Int s = 0; s < shcol; ++s) {
      SDS.check[s] = true;
      SDS.thv_ref[s] = 300;
      SDS.pblh[s] = 0;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Make all input, homogenous
        SDS.u[offset] = 5;
        SDS.v[offset] = -2;
        SDS.thv[offset] = 300;
        // Should be initialized to zero everywhere
        SDS.rino[offset] = 0;
      }
    }

    // Call the fortran implementation
    pblintd_height(SDS);

    // Check that PBLH is zero (not modified) everywhere
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] == approx_zero);
    }

  } // run_property

  static void run_bfb()
  {
    PblintdHeightData f90_data[] = {
      PblintdHeightData(10, 72),
      PblintdHeightData(10, 12),
      PblintdHeightData(7, 16),
      PblintdHeightData(2, 7),
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(PblintdHeightData);

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    PblintdHeightData cxx_data[] = {
      PblintdHeightData(f90_data[0]),
      PblintdHeightData(f90_data[1]),
      PblintdHeightData(f90_data[2]),
      PblintdHeightData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      pblintd_height(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      pblintd_height_f(d.shcol, d.nlev, d.z, d.u, d.v, d.ustar, d.thv, d.thv_ref, d.pblh, d.rino, d.check);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      PblintdHeightData& d_f90 = f90_data[i];
      PblintdHeightData& d_cxx = cxx_data[i];
      for (Int k = 0; k < d_f90.total(d_f90.pblh); ++k) {
        REQUIRE(d_f90.total(d_f90.pblh) == d_cxx.total(d_cxx.pblh));
        REQUIRE(d_f90.pblh[k] == d_cxx.pblh[k]);
        REQUIRE(d_f90.total(d_f90.pblh) == d_cxx.total(d_cxx.check));
        REQUIRE(d_f90.check[k] == d_cxx.check[k]);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("pblintd_height_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdHeight;

  TestStruct::run_property();
}

TEST_CASE("pblintd_height_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdHeight;

  TestStruct::run_bfb();
}

} // empty namespace
