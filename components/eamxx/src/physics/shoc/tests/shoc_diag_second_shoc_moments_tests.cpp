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
struct UnitWrap::UnitTest<D>::TestDiagSecondShocMoments {

  static void run_property()
  {

    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Property test for the upper level function:
    //  diag_second_moments

    // TEST
    //  Feed routine profiles with well mixed, unstable, and stable layers
    //  Verify output is as physically expected and that output falls within
    //  a reasonable range.

    // Define the liquid water potential temperature [K]
    static constexpr Real thetal[nlev] = {310, 307, 302, 302, 303};
    // Define the total water mixing ratio [kg/kg]
    static constexpr Real qw[nlev] = {1e-2, 1.2e-2, 1.5e-2, 1.5e-2, 1.4e-2};
    // Define the zonal wind [m/s]
    static constexpr Real u_wind[nlev] = {4, 4, 2, 0, -1};
    // define the meridional wind [m/s]
    static constexpr Real v_wind[nlev] = {-2, -2, 1, 3, 0};
    // Define the TKE [m2/s2]
    static constexpr Real tke[nlev] = {0.2, 0.3, 0.5, 0.4, 0.1};
    // Define the heights on the zi grid [m]
    static constexpr Real zi_grid[nlevi] = {900, 500, 150, 90, 50, 0};
    // Define the return to isotropy timescale [s]
    static constexpr Real isotropy[nlev] = {1000, 2000, 4000, 1000, 500};
    // Define the eddy vicosity for heat and momentum [m2/s]
    static constexpr Real tkh[nlev] = {3, 10, 50, 30, 20};
    // Define the mixing length [m]
    static constexpr Real shoc_mix[nlev] = {300, 500, 2000, 1500, 300};

    // Surface heat flux [K m/s]
    static constexpr Real wthl_sfc[shcol] = {2e-2, -2e-2};
    // Surface water flux [kg/kg m/s]
    static constexpr Real wqw_sfc[shcol] = {5e-3, -5e-3};
    // Surface zonal wind flux [m2/s2]
    static constexpr Real uw_sfc[shcol] = {2e-2, -4e-2};
    // Surface meridional wind flux [m2/s2]
    static constexpr Real vw_sfc[shcol] = {1e-2, -1e-4};

    // establish reasonable bounds for checking result
    static constexpr Real thl_sec_ubound = 1000; // [K^s]
    static constexpr Real qwthl_sec_bound = 1e-1; // [K kg/kg]
    static constexpr Real wthl_sec_bound = 10; // [K m/s]
    static constexpr Real wqw_sec_bound = 1e-3; // [kg/kg m/s]
    static constexpr Real uw_sec_bound = 10; // [m2/s2]
    static constexpr Real vw_sec_bound = 10; // [m2/s2]
    static constexpr Real wtke_sec_bound = 1; // [m3/s3]
    static constexpr Real wtke_srf_ubound = 0.1; // [m3/s3]

    // Compute height information based on zi_grid
    Real zt_grid[nlev];
    Real dz_zi[nlevi];
    // Compute heights on midpoint grid
    for(Int n = 0; n < nlev; ++n) {
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
      if (n == 0){
        dz_zi[n] = 0;
      }
      else{
        dz_zi[n] = zt_grid[n-1] - zt_grid[n];
      }
    }
    // set upper condition for dz_zi
    dz_zi[nlevi-1] = zt_grid[nlev-1];

    DiagSecondShocMomentsData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(SDS.nlevi == SDS.nlev+1);

    // Load up the new data
    for(Int s = 0; s < shcol; ++s) {
      SDS.wthl_sfc[s] = wthl_sfc[s];
      SDS.wqw_sfc[s] = wqw_sfc[s];
      SDS.uw_sfc[s] = uw_sfc[s];
      SDS.vw_sfc[s] = vw_sfc[s];
      // Fill in test data on zt_grid.
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.zt_grid[offset] = zt_grid[n];
        SDS.tke[offset] = tke[n];
        SDS.thetal[offset] = thetal[n];
        SDS.qw[offset] = qw[n];
        SDS.u_wind[offset] = u_wind[n];
        SDS.v_wind[offset] = v_wind[n];
        SDS.isotropy[offset] = isotropy[n];
        SDS.shoc_mix[offset] = shoc_mix[n];
        // Intentionally feed tk and tkh the same value
        SDS.tk[offset] = tkh[n];
        SDS.tkh[offset] = tkh[n];
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.dz_zi[offset] = dz_zi[n];
        SDS.zi_grid[offset] = zi_grid[n];
      }

    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      // nlevi loop
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
        // Check that zi increases updward
        if (n < nlevi-1){
          REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0);
        }
      }
      // nlev loop
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Check that zt increases upward
        if (n < nlev-1){
          REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
        }
        REQUIRE(SDS.shoc_mix[offset] > 0);
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.tkh[offset] > 0);
        REQUIRE(SDS.tk[offset] > 0);
        REQUIRE(SDS.qw[offset] > 0);
        REQUIRE(SDS.thetal[offset] > 0);
        REQUIRE(SDS.isotropy[offset] > 0);
      }
    }

    // Call the fortran implementation
    diag_second_shoc_moments(SDS);

    // Verify output makes sense
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Verify TKE is in reasonable range
        REQUIRE(SDS.w_sec[offset] > 0);
        REQUIRE(SDS.w_sec[offset] < SDS.tke[offset]);
      }
      for (Int n = 0; n < nlevi; ++n){
        const auto offset = n + s * nlevi;

        // Validate upper boundary points
        if (n == 0){
          SDS.thl_sec[offset]=0;
          SDS.qw_sec[offset]=0;
          SDS.qwthl_sec[offset]=0;
          SDS.wthl_sec[offset]=0;
          SDS.wqw_sec[offset]=0;
          SDS.uw_sec[offset]=0;
          SDS.vw_sec[offset]=0;
          SDS.wtke_sec[offset]=0;
        }
        // Validate lower boundary points
        else if (n == nlevi-1){
          REQUIRE(SDS.thl_sec[offset] >= 0);
          REQUIRE(SDS.qw_sec[offset] >= 0);
          REQUIRE(SDS.qwthl_sec[offset] >= 0);

          // Verify bounds of wtke_sec make sense
          REQUIRE(SDS.wtke_sec[offset] >= 0);
          REQUIRE(SDS.wtke_sec[offset] < wtke_srf_ubound);

          // Some outputs should be equal to inputs
          REQUIRE(SDS.wthl_sec[offset] == wthl_sfc[s]);
          REQUIRE(SDS.wqw_sec[offset] == wqw_sfc[s]);
          REQUIRE(SDS.uw_sec[offset] == uw_sfc[s]);
          REQUIRE(SDS.vw_sec[offset] == vw_sfc[s]);
        }
        else{
          // Validate Downgradient assumption for
          //   various possible layers.  Note that these physics
          //   have been tested in lower level functions, but
          //   as a double check, will do here for ONE variable.
          // conditionally stable layer
          if ((thetal[n-1] - thetal[n]) > 0){
            REQUIRE(SDS.wthl_sec[offset] < 0);
          }
          // well mixed layer
          if ((thetal[n-1] - thetal[n]) == 0){
            REQUIRE(SDS.wthl_sec[offset] == 0);
          }
          // unstable layer
          if ((thetal[n-1] - thetal[n]) < 0){
            REQUIRE(SDS.wthl_sec[offset] > 0);
          }

          // Validate the covariance formulation.
          // well mixed layer test
          if ((thetal[n-1] - thetal[n]) == 0 ||
              (qw[n-1] - qw[n]) == 0){
            REQUIRE(SDS.qwthl_sec[offset] == 0);
          }

          // validate values are NEGATIVE if potential
          //  temperature INCREASES with height and total water
          //  DECREASES with height
          if ((thetal[n-1] - thetal[n]) > 0 &&
              (qw[n-1] - qw[n]) < 0){
            REQUIRE(SDS.qwthl_sec[offset] < 0);
          }

          // validate values are POSITIVE if both
          //   potential temperature and total water
          //   DECREASE with height
          if ((thetal[n-1] - thetal[n]) < 0 &&
              (qw[n-1] - qw[n]) < 0){
            REQUIRE(SDS.qwthl_sec[offset] > 0);
          }

          // Now check that all output falls in reasonable bounds
          REQUIRE(SDS.thl_sec[offset] >= 0);
          REQUIRE(SDS.thl_sec[offset] < thl_sec_ubound);

          REQUIRE(SDS.qw_sec[offset] >= 0);
          REQUIRE(SDS.qw_sec[offset] < thl_sec_ubound);

          REQUIRE(std::abs(SDS.qwthl_sec[offset]) < qwthl_sec_bound);
          REQUIRE(std::abs(SDS.wthl_sec[offset]) < wthl_sec_bound);
          REQUIRE(std::abs(SDS.qwthl_sec[offset]) < qwthl_sec_bound);
          REQUIRE(std::abs(SDS.wqw_sec[offset]) < wqw_sec_bound);
          REQUIRE(std::abs(SDS.uw_sec[offset]) < uw_sec_bound);
          REQUIRE(std::abs(SDS.vw_sec[offset]) < vw_sec_bound);
          REQUIRE(std::abs(SDS.wtke_sec[offset]) < wtke_sec_bound);
        }
      }
    }

  } // run_property

  static void run_bfb()
  {
    auto engine = setup_random_test();

    DiagSecondShocMomentsData f90_data[] = {
      DiagSecondShocMomentsData(36,  72, 73),
      DiagSecondShocMomentsData(72,  72, 73),
      DiagSecondShocMomentsData(128, 72, 73),
      DiagSecondShocMomentsData(256, 72, 73),
    };

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    DiagSecondShocMomentsData cxx_data[] = {
      DiagSecondShocMomentsData(f90_data[0]),
      DiagSecondShocMomentsData(f90_data[1]),
      DiagSecondShocMomentsData(f90_data[2]),
      DiagSecondShocMomentsData(f90_data[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      diag_second_shoc_moments(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      diag_second_shoc_moments_f(d.shcol, d.nlev, d.nlevi, d.thetal, d.qw, d.u_wind, d.v_wind, d.tke, d.isotropy,
                d.tkh, d.tk, d.dz_zi, d.zt_grid, d.zi_grid, d.shoc_mix, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.thl_sec,
                d.qw_sec, d.wthl_sec, d.wqw_sec, d.qwthl_sec, d.uw_sec, d.vw_sec, d.wtke_sec, d.w_sec);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(DiagSecondShocMomentsData);
      for (Int i = 0; i < num_runs; ++i) {
        DiagSecondShocMomentsData& d_f90 = f90_data[i];
        DiagSecondShocMomentsData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.w_sec); ++k) {
          REQUIRE(d_f90.w_sec[k] == d_cxx.w_sec[k]);
        }
        for (Int k = 0; k < d_f90.total(d_f90.thl_sec); ++k) {
          REQUIRE(d_f90.thl_sec[k] == d_cxx.thl_sec[k]);
          REQUIRE(d_f90.qw_sec[k] == d_cxx.qw_sec[k]);
          REQUIRE(d_f90.wthl_sec[k] == d_cxx.wthl_sec[k]);
          REQUIRE(d_f90.wqw_sec[k] == d_cxx.wqw_sec[k]);
          REQUIRE(d_f90.qwthl_sec[k] == d_cxx.qwthl_sec[k]);
          REQUIRE(d_f90.uw_sec[k] == d_cxx.uw_sec[k]);
          REQUIRE(d_f90.vw_sec[k] == d_cxx.vw_sec[k]);
          REQUIRE(d_f90.wtke_sec[k] == d_cxx.wtke_sec[k]);
        }
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("diag_second_shoc_moments_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDiagSecondShocMoments;

  TestStruct::run_property();
}

TEST_CASE("diag_second_shoc_moments_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestDiagSecondShocMoments;

  TestStruct::run_bfb();
}

} // empty namespace
