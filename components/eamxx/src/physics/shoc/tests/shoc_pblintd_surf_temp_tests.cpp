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
struct UnitWrap::UnitTest<D>::TestPblintdSurfTemp {

  static void run_property()
  {
    static constexpr auto ustar_min = scream::shoc::Constants<Scalar>::ustar_min;
    static constexpr Int shcol = 4;
    static constexpr Int nlev = 4;
    static constexpr Int nlevi = nlev+1;

    // Tests for the subroutine pblintd_surf_temp

    // Define mid point height [m]
    static constexpr Real z[nlev] = {1500, 1000, 500, 20};
    // Define virtual potential temperature [K]
    static constexpr Real thv[nlev]= {305, 302, 302, 300};
    // Define obklen [m]
    static constexpr Real obklen[shcol] = {-10, 20, -100, 150};
    // Define surf friction velocity [m4/s4]
    static constexpr Real ustar[shcol] = {ustar_min, 1, 5, 10};
    // Surface kinematic buoyancy flux
    static constexpr Real kbfs[shcol] = {0.03, -0.03, 0.1, -0.1};
    // Input check parameter
    static constexpr bool check[shcol] = {true, true, false, false};
    // Input bulk richardson number
    static constexpr Real rino = 0.1;

    // Initialize data structure for bridging to F90
    PblintdSurfTempData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(shcol > 0);
    REQUIRE(nlevi == nlev+1);

    // Fill in test data on the zt_grid
    for(Int s = 0; s < shcol; ++s) {
      SDS.ustar[s] = ustar[s];
      SDS.check[s] = check[s];
      SDS.kbfs[s] = kbfs[s];
      SDS.obklen[s] = obklen[s];
      SDS.pblh[s] = 0; // Initialize to zero
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        SDS.z[offset] = z[n];
        SDS.thv[offset] = thv[n];
        SDS.rino[offset] = rino;
      }
    }

    // check to make sure the input data makes sense
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
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.z[offset] > 0);
        REQUIRE(SDS.thv[offset] < 350);
        REQUIRE(SDS.thv[offset] > 290);
      }
    }

    // Call the fortran implementation
    pblintd_surf_temp(SDS);

    // Check the result
    for(Int s = 0; s < shcol; ++s) {
      // If OUTPUT check was performed, verify that surface diagnosed
      //  temperature is appropriate
      if (SDS.check[s] == true){
        REQUIRE(SDS.tlv[s] < 350);
        REQUIRE(SDS.tlv[s] > 290);

       // determine lowest level index for this column
       const auto low_lev = (nlev-1) + s * nlev;
       REQUIRE(SDS.rino[low_lev] == 0);
      }

      // If INPUT check was performed, verify that PBLH depth was
      //  set to the expected value, which should be the highest level
      //  of the data provided
      if (check[s] == true){
        REQUIRE(SDS.pblh[s] == SDS.z[0]);
      }
    }

  } // run_property

  static void run_bfb()
  {
    auto engine = setup_random_test();

    PblintdSurfTempData f90_data[] = {
      PblintdSurfTempData(6, 7, 8),
      PblintdSurfTempData(64, 72, 73),
      PblintdSurfTempData(128, 72, 73),
      PblintdSurfTempData(256, 72, 73),
    };

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine, { {d.obklen, {100., 200.}} });
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    PblintdSurfTempData cxx_data[] = {
      PblintdSurfTempData(f90_data[0]),
      PblintdSurfTempData(f90_data[1]),
      PblintdSurfTempData(f90_data[2]),
      PblintdSurfTempData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      pblintd_surf_temp(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      pblintd_surf_temp_f(d.shcol, d.nlev, d.nlevi, d.z, d.ustar, d.obklen, d.kbfs, d.thv, d.tlv, d.pblh, d.check, d.rino);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(PblintdSurfTempData);
      for (Int i = 0; i < num_runs; ++i) {
        PblintdSurfTempData& d_f90 = f90_data[i];
        PblintdSurfTempData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.tlv); ++k) {
          REQUIRE(d_f90.total(d_f90.tlv) == d_cxx.total(d_cxx.tlv));
          REQUIRE(d_f90.tlv[k] == d_cxx.tlv[k]);
          REQUIRE(d_f90.total(d_f90.tlv) == d_cxx.total(d_cxx.pblh));
          REQUIRE(d_f90.pblh[k] == d_cxx.pblh[k]);
          REQUIRE(d_f90.total(d_f90.tlv) == d_cxx.total(d_cxx.check));
          REQUIRE(d_f90.check[k] == d_cxx.check[k]);
        }
        for (Int k = 0; k < d_f90.total(d_f90.rino); ++k) {
          REQUIRE(d_f90.total(d_f90.rino) == d_cxx.total(d_cxx.rino));
          REQUIRE(d_f90.rino[k] == d_cxx.rino[k]);
        }
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("pblintd_surf_temp_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdSurfTemp;

  TestStruct::run_property();
}

TEST_CASE("pblintd_surf_temp_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdSurfTemp;

  TestStruct::run_bfb();
}

} // empty namespace
