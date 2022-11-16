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
struct UnitWrap::UnitTest<D>::TestDiagSecondMomentsLbycond {

  static void run_property()
  {
    // Property tests for the SHOC function
    //  diag_second_moments_lbycond

    static constexpr Int shcol = 5;

    // TEST
    // Feed function several columns worth of varying input
    //  and verify output is physical and as expected.

    // Surface heat flux [K m/s]
    static constexpr Real wthl_sfc[shcol] = {2e-2, 1e-1, 0, -4e-3, -2e-2};
    // Surface water flux [kg/kg m/s]
    static constexpr Real wqw_sfc[shcol] = {5e-3, 1e-4, 0, -1e-4, -5e-3};
    // Surface zonal wind flux [m2/s2]
    static constexpr Real uw_sfc[shcol] = {2e-2, -1e-2, 0, -5e-3, -4e-2};
    // Surface meridional wind flux [m2/s2]
    static constexpr Real vw_sfc[shcol] = {1e-2, 5e-2, 0, -1e-2, -1e-4};
    // Surface convective velocity scale [m/s]
    static constexpr Real wstar[shcol] = {1, 0.5, 0, 0, 0};
    // Surface friction velocity [m4/s4]
    Real ustar2[shcol];

    // establish reasonable bound for output checking
    static constexpr Real wtke_srf_ubound = 0.1; // [m3/s3]

    // Compute other inputs and check some inputs
    for (Int s = 0; s < shcol; ++s){
      // compute surface friction velocity
      ustar2[s] = std::sqrt(uw_sfc[s]*uw_sfc[s] + vw_sfc[s]*vw_sfc[s]);
      // Be sure wstar input is consistent with wthl_sfc
      if (wthl_sfc[s] <= 0){
        REQUIRE(wstar[s] == 0);
      }
    }

    // Initialize data structure for bridging to F90
    DiagSecondMomentsLbycondData SDS(shcol);

    // Load up input data to structure
    for (Int s = 0; s < shcol; ++s){
      SDS.wthl_sfc[s] = wthl_sfc[s];
      SDS.wqw_sfc[s] = wqw_sfc[s];
      SDS.uw_sfc[s] = uw_sfc[s];
      SDS.vw_sfc[s] = vw_sfc[s];
      SDS.wstar[s] = wstar[s];
      SDS.ustar2[s] = ustar2[s];
    }

    // Test that the inputs are reasonable
    REQUIRE(SDS.shcol == shcol);
    REQUIRE(shcol > 0);

    // Verify input is as expected
    for (Int s = 0; s < shcol; ++s){
      REQUIRE(SDS.ustar2[s] >= 0);
    }

    // Call the fortran implementation
    diag_second_moments_lbycond(SDS);

    // Verify output is as expected
    for (Int s = 0; s < shcol; ++s){
      // Some outputs should be equal to inputs
      REQUIRE(SDS.wthl_sec[s] == wthl_sfc[s]);
      REQUIRE(SDS.wqw_sec[s] == wqw_sfc[s]);
      REQUIRE(SDS.uw_sec[s] == uw_sfc[s]);
      REQUIRE(SDS.vw_sec[s] == vw_sfc[s]);

      // Verify bounds of wtke_sec make sense
      REQUIRE(SDS.wtke_sec[s] >= 0);
      REQUIRE(SDS.wtke_sec[s] < wtke_srf_ubound);

      // Verify second moments make sense
      if (SDS.wthl_sfc[s] == 0){
        REQUIRE(SDS.thl_sec[s] == 0);
      }
      else{
        REQUIRE(SDS.thl_sec[s] > 0);
      }

      if (SDS.wqw_sfc[s] == 0){
        REQUIRE(SDS.qw_sec[s] == 0);
      }
      else{
        REQUIRE(SDS.qw_sec[s] > 0);
      }

      if ( (SDS.wqw_sfc[s] == 0 && SDS.wthl_sfc[s] == 0) ){
        REQUIRE(SDS.qwthl_sec[s] == 0);
      }
      else{
        REQUIRE(SDS.qwthl_sec[s] > 0);
      }
    }

  } // run_property

  static void run_bfb()
  {
    auto engine = setup_random_test();

    DiagSecondMomentsLbycondData f90_data[] = {
      DiagSecondMomentsLbycondData(120),
      DiagSecondMomentsLbycondData(120),
      DiagSecondMomentsLbycondData(120),
      DiagSecondMomentsLbycondData(120),
    };

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    DiagSecondMomentsLbycondData cxx_data[] = {
      DiagSecondMomentsLbycondData(f90_data[0]),
      DiagSecondMomentsLbycondData(f90_data[1]),
      DiagSecondMomentsLbycondData(f90_data[2]),
      DiagSecondMomentsLbycondData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      diag_second_moments_lbycond(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      diag_second_moments_lbycond_f(d.shcol, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.ustar2, d.wstar,
         d.wthl_sec, d.wqw_sec, d.uw_sec, d.vw_sec, d.wtke_sec, d.thl_sec, d.qw_sec, d.qwthl_sec);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(DiagSecondMomentsLbycondData);
      for (Int i = 0; i < num_runs; ++i) {
        DiagSecondMomentsLbycondData& d_f90 = f90_data[i];
        DiagSecondMomentsLbycondData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.shcol; ++k) {
          REQUIRE(d_f90.wthl_sec[k] == d_cxx.wthl_sec[k]);
          REQUIRE(d_f90.wqw_sec[k] == d_cxx.wqw_sec[k]);
          REQUIRE(d_f90.uw_sec[k] == d_cxx.uw_sec[k]);
          REQUIRE(d_f90.vw_sec[k] == d_cxx.vw_sec[k]);
          REQUIRE(d_f90.wtke_sec[k] == d_cxx.wtke_sec[k]);
          REQUIRE(d_f90.thl_sec[k] == d_cxx.thl_sec[k]);
          REQUIRE(d_f90.qw_sec[k] == d_cxx.qw_sec[k]);
          REQUIRE(d_f90.qwthl_sec[k] == d_cxx.qwthl_sec[k]);
        }
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("diag_second_moments_lbycond_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDiagSecondMomentsLbycond;

  TestStruct::run_property();
}

TEST_CASE("diag_second_moments_lbycond_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDiagSecondMomentsLbycond;

  TestStruct::run_bfb();
}

} // empty namespace
