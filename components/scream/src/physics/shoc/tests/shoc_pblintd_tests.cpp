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
struct UnitWrap::UnitTest<D>::TestPblintd {

  static void run_property()
  {
    const auto ustar_min = scream::shoc::Constants<Scalar>::ustar_min;
    static constexpr Int shcol = 5;
    static constexpr Int nlev = 5;
    static constexpr Int nlevi = nlev+1;

    // Tests for the subroutine pblintd

    // TEST
    // This is the top level routine to diagnose PBL height.  The individual routines
    //  have been throughly tests.  Thus, give subroutine realsitic inputs and be sure
    //  that the output PBL height is physical and expected.

    // Define the heights on the zi grid [m]
    static constexpr Real zi_grid[nlevi] = {3000, 2000, 1500, 1000, 500, 0};
    // Define the liquid water potential temperature [K]
    static constexpr Real thetal[nlev] = {320, 315, 310, 305, 303};
    // Define cloud water mixing ratio [kg/kg]
    static constexpr Real ql[nlev] = {0, 0, 1e-5, 2e-5, 0};
    // Define water vapor [kg/kg]
    static constexpr Real q[nlev] = {5e-3, 6e-3, 7e-3, 8e-3, 10e-3};
    // Define the zonal wind [m/s]
    static constexpr Real u_wind[nlev] = {4, 4, 2, 0, -1};
    // define the meridional wind [m/s]
    static constexpr Real v_wind[nlev] = {-2, -2, 1, 3, 0};
    // Define surface friction velocity [m/s]
    static constexpr Real ustar[shcol] = {0.1, 4.0, 0.9, 2.0, ustar_min};
    // Define obklen [m]
    static constexpr Real obklen[shcol] = {-10, 20, -100, 150, 10};
    // Surface kinematic buoyancy flux
    static constexpr Real kbfs[shcol] = {0.03, -0.03, 0.1, -0.1, -0.01};

    // First compute variables related to height
    Real zt_grid[nlev];
    for(Int n = 0; n < nlev; ++n) {
      // height on the midpoint grid
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
    }

    // Initialize data structure for bridging to F90
    PblintdData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(shcol > 0);
    REQUIRE(nlevi == nlev+1);

    // Fill in test data on the zt_grid
    for(Int s = 0; s < shcol; ++s) {
      SDS.ustar[s] = ustar[s];
      SDS.kbfs[s] = kbfs[s];
      SDS.obklen[s] = obklen[s];
      SDS.pblh[s] = 0; // Initialize to zero
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        SDS.z[offset] = zt_grid[n];
        SDS.thl[offset] = thetal[n];
        SDS.u[offset] = u_wind[n];
        SDS.v[offset] = v_wind[n];
        SDS.q[offset] = q[n];
        SDS.ql[offset] = ql[n];
      }
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
        SDS.zi[offset] = zi_grid[n];
      }
    }

    // check that inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev - 1; ++n) {
        const auto offset = n + s * nlev;
        // Check that zt increases upward
        REQUIRE(SDS.z[offset + 1] - SDS.z[offset] < 0);
      }

      // Check that zi increases upward
      for(Int n = 0; n < nlevi - 1; ++n) {
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.zi[offset + 1] - SDS.zi[offset] < 0);
      }
    }

    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.ustar[s] >= ustar_min);
      REQUIRE(std::abs(SDS.kbfs[s]) < 1);
      // Make sure the sign of the Monin Obukov length is
      //  consistent with the sign of kinematic buoyancy flux
      if (SDS.kbfs[s] < 0){
        REQUIRE(SDS.obklen[s] > 0);
      }
      else{
        REQUIRE(SDS.obklen[s] < 0);
      }

      for(Int n = 0; n < nlev - 1; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE( (SDS.thl[offset] > 150 && SDS.thl[offset] < 400) );
        REQUIRE( (SDS.ql[offset] >= 0 && SDS.ql[offset] < 1) );
        REQUIRE( (SDS.q[offset] >= 0 && SDS.q[offset] < 1) );
        REQUIRE( (SDS.cldn[offset] >= 0 && SDS.cldn[offset] < 1) );
      }
    }

    // Call the fortran implementation
    pblintd(SDS);

    // Make sure PBL height is reasonable
    // Should be larger than second lowest zi level and lower
    //   (or equal to) the highest zi level
    for(Int s = 0; s < shcol; ++s) {
      const auto lowlev = (nlevi-1)+ s * nlevi;
      const auto highlev = s * nlevi;
      REQUIRE(SDS.pblh[s] > SDS.zi[lowlev]);
      REQUIRE(SDS.pblh[s] < SDS.zi[highlev]);
    }

  }

  static void run_bfb()
  {
    PblintdData f90_data[] = {
      PblintdData(10, 71, 72),
      PblintdData(10, 12, 13),
      PblintdData(7,  16, 17),
      PblintdData(2, 7, 8),
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(PblintdData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    PblintdData cxx_data[] = {
      PblintdData(f90_data[0]),
      PblintdData(f90_data[1]),
      PblintdData(f90_data[2]),
      PblintdData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      pblintd(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      pblintd_f(d.shcol, d.nlev, d.nlevi, d.z, d.zi, d.thl, d.ql, d.q, d.u, d.v, d.ustar, d.obklen, d.kbfs, d.cldn, d.pblh);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      PblintdData& d_f90 = f90_data[i];
      PblintdData& d_cxx = cxx_data[i];
      for (Int k = 0; k < d_f90.total(d_f90.pblh); ++k) {
        REQUIRE(d_f90.pblh[k] == d_cxx.pblh[k]);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("pblintd_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintd;

  TestStruct::run_property();
}

TEST_CASE("pblintd_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintd;

  TestStruct::run_bfb();
}

} // empty namespace
