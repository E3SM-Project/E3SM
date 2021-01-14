#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"

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
struct UnitWrap::UnitTest<D>::TestShocAssumedPdf {

  static void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the top level subroutine
    //   shoc_assumed_pdf

    // Tests will start simple, and gradually add complexity to test
    //  the physics.
    // NOTE: for this test we want exactly two columns.

    // TEST ONE
    // No SGS variability test.  Given inputs where there is a saturated
    //  profile but NO SGS variability in the scalar fluxes or variances
    //  (i.e. all second and third moment terms are zero), then verify that
    //  that cloud fraction is either 1 or 0 and that the SGS variability
    //  outputs are also zero everywhere.

    // Define input data

    // Note that the moisture and height profiles below represent that
    //  of the BOMEX case, but the temperatures are much colder, to encourage
    //  there to be points with ample cloud produced for this test.

    // Liquid water potential temperature [K]
    static constexpr Real thetal[nlev] = {303, 300, 298, 298, 300};
    // Total water mixing ratio [kg/kg]
    static constexpr Real qw[nlev] = {0.003, 0.004, 0.011, 0.016, 0.017};
    // Pressure [Pa]
    static constexpr Real pres[nlev] = {70000, 80000, 85000, 90000, 100000};
    // Define the heights on the zt grid [m]
    static constexpr Real zi_grid[nlevi] = {3000, 2000, 1500, 1000, 500, 0};

    // All variances will be given zero or minimum threshold inputs

    // Define some reasonable bounds for output
    static constexpr Real wqls_bound = 0.1;
    static constexpr Real wthv_sec_bound = 10;
    static constexpr Real shoc_ql2_bound = 0.1;
    static constexpr Real shoc_ql_bound = 0.1;

    Real zt_grid[nlev];
    // Compute heights on midpoint grid
    for(Int n = 0; n < nlev; ++n) {
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
    }

    // Initialize data structure for bridging to F90
    ShocAssumedPdfData SDS(shcol, nlev, nlevi);

    // Load input data
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.thetal[offset] = thetal[n];
        SDS.qw[offset] = qw[n];
        SDS.pres[offset] = pres[n];
        SDS.zt_grid[offset] = zt_grid[n];
        SDS.w_field[offset] = 0;
        SDS.w_sec[offset] = 0.004;
      }

      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.thl_sec[offset] = 0;
        SDS.qw_sec[offset] = 0;
        SDS.wthl_sec[offset] = 0;
        SDS.wqw_sec[offset] = 0;
        SDS.qwthl_sec[offset] = 0;
        SDS.w3[offset] = 0;
        SDS.zi_grid[offset] = zi_grid[n];
      }
    }

    // Check that the inputs make sense

    // Load input data
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.qw[offset] > 0);
        REQUIRE(SDS.thetal[offset] > 0);
        REQUIRE(SDS.pres[offset] > 0);
        REQUIRE(SDS.zt_grid[offset] > 0);
      }

      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.zi_grid[offset] >= 0);
      }
    }

    // Check that zt increase updward and pres decrease upward
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev - 1; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
        REQUIRE(SDS.pres[offset + 1] - SDS.pres[offset] > 0);
      }

      // Check that zi increase upward
      for(Int n = 0; n < nlevi - 1; ++n) {
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0);
      }

    }

    // Test that the inputs are reasonable.
    REQUIRE(SDS.nlevi - SDS.nlev == 1);
    // For this test we want exactly two columns
    REQUIRE(SDS.shcol == 2);

    // Call the fortran implementation
    shoc_assumed_pdf(SDS);

    // Verify the result
    // Make sure cloud fraction is either 1 or 0 and all
    //  SGS terms are zero.

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE( (SDS.shoc_cldfrac[offset] == 0  || SDS.shoc_cldfrac[offset] == 1) );
        REQUIRE(SDS.wqls[offset] == 0);
        REQUIRE(SDS.wthv_sec[offset] == 0);
        REQUIRE(SDS.shoc_ql2[offset] == 0);
        REQUIRE(SDS.shoc_ql[offset] >= 0);

      }
    }

    // TEST TWO
    // Add in Scalar fluxes.  This should give us a nonzero
    //  buoyancy flux everywhere but other SGS terms should remain zero

    // We will assume turbulence with a uniform profile

    // Flux of liquid water [K m/s]
    static constexpr Real wthl_sec = -0.03;
    // Flux of total water [m/s kg/kg]
    static constexpr Real wqw_sec = 0.0002;

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.wthl_sec[offset] = wthl_sec;
        SDS.wqw_sec[offset] = wqw_sec;
      }
    }

    // Call the fortran implementation
    shoc_assumed_pdf(SDS);

    // Verify the result
    // Make sure cloud fraction is either 1 or 0 and all
    //  SGS terms are zero, EXCEPT wthv_sec.

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE( (SDS.shoc_cldfrac[offset] == 0  || SDS.shoc_cldfrac[offset] == 1) );
        REQUIRE(SDS.wqls[offset] == 0);
        REQUIRE(SDS.wthv_sec[offset] != 0);
        REQUIRE(std::abs(SDS.wthv_sec[offset] < wthv_sec_bound));
        REQUIRE(SDS.shoc_ql2[offset] == 0);
        REQUIRE(SDS.shoc_ql[offset] >= 0);
        REQUIRE(SDS.shoc_ql[offset] < shoc_ql_bound);
      }
    }

    // TEST THREE and FOUR
    // Add in Scalar variances, and POSITIVE vertical velocity skewness test.

    // Add strong scalar variances as such that will produce cloud at every level.

    // For the first column feed in zero vertical velocity skewness.
    // For the second column feed in large veriticle velocity skewss.
    // Verify that for points where cloud fraction was < 0.5 in the first column,
    //  that cloud fraction then vice versa for points with cloud fraction > 0.5.

    // Thetal variance [K^2]
    static constexpr Real thl_sec = 2;
    // total water variance [kg^2/kg^2]
    static constexpr Real qw_sec = 0.0002;
    // Temperature and total water covariance [K kg/kg]
    static constexpr Real qwthl_sec = 1e-5;
    // Vertical velocity variance [m2/s2]
    static constexpr Real w_sec = 0.2;

    // Load input data
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.w_sec[offset] = w_sec;
      }

      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.thl_sec[offset] = thl_sec;
        SDS.qw_sec[offset] = qw_sec;
        SDS.qwthl_sec[offset] = qwthl_sec;
        SDS.w3[offset] = s*1.0;
      }
    }

    // Call the fortran implementation
    shoc_assumed_pdf(SDS);

    // Check the result

    // With such a turbulence and scalar profile, this should have
    //  encouraged cloud everywhere.  Verify that this is true.

    // Also verify that output lies within some reasonable bounds.

    // Then verify vertical velocity skewness info

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
        const auto offsets = n + (s+1) * nlevi;
        if (s < shcol-1){
          // Verify input w3 is greater in subsequent columns
          REQUIRE(SDS.w3[offsets] > SDS.w3[offset]);
        }
      }

      // Verify output falls within reasonable bounds.  For this positive
      //  vertical velocity skewness test and with the give inputs there
      //  should be cloud everywhere and all flux terms should be positive.

      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        const auto offsets = n + (s+1) * nlev;

        REQUIRE( (SDS.shoc_cldfrac[offset] > 0  || SDS.shoc_cldfrac[offset] < 1) );
        REQUIRE(SDS.wqls[offset] > 0);
        REQUIRE(SDS.wthv_sec[offset] > 0);
        REQUIRE(SDS.shoc_ql2[offset] > 0);
        REQUIRE(SDS.shoc_ql[offset] > 0);

        REQUIRE(SDS.wqls[offset] < 0.1);
        REQUIRE(SDS.wthv_sec[offset] < wthv_sec_bound);
        REQUIRE(SDS.shoc_ql2[offset] < shoc_ql2_bound);
        REQUIRE(SDS.shoc_ql[offset] < shoc_ql_bound);

        // Now verify that the relationships in a strongly positive vertical
        //  velocity flux regime hold true, relative to a symmetric vertical
        //  velocity regime.
        if (s < shcol-1){

          if (SDS.shoc_cldfrac[offset] < 0.5){
            REQUIRE(SDS.shoc_cldfrac[offsets] < SDS.shoc_cldfrac[offset]);
          }
          else if (SDS.shoc_cldfrac[offset] > 0.5){
            REQUIRE(SDS.shoc_cldfrac[offsets] > SDS.shoc_cldfrac[offset]);
          }

          // In addition, in a positive skewness environment, the following
          //  should also be true

          // Grid mean liquid water decreased
          REQUIRE(SDS.shoc_ql[offsets] < SDS.shoc_ql[offset]);
          // liquid water flux increased
          REQUIRE(SDS.wqls[offsets] > SDS.wqls[offset]);
          // buoyancy flux increased
          REQUIRE(SDS.wthv_sec[offsets] > SDS.wthv_sec[offset]);
          // liquid water variance increased
          REQUIRE(SDS.shoc_ql2[offsets] > SDS.shoc_ql2[offset]);
        }
      }
    }

    // TEST FIVE
    // Negative vertical velocity skewness.

    // Using same input as the above test, feed one column with zero skeweness
    //  and another test with negative vertical velocity skewness and verify
    //  result is physical.

    // Load input data
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.w3[offset] = s*-1.0;
      }
    }

    // Call the fortran implementation
    shoc_assumed_pdf(SDS);

    // Check the result

    // Verify that output lies within some reasonable bounds.

    // Then verify vertical velocity skewness info

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
        const auto offsets = n + (s+1) * nlevi;
        if (s < shcol-1){
          // Verify input w3 is greater in subsequent columns
          REQUIRE(SDS.w3[offsets] < SDS.w3[offset]);
        }
      }

      // Verify output falls within reasonable bounds
      //   For this negative vertical velocity test some variables
      //   will be expected to be less than zero.

      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        const auto offsets = n + (s+1) * nlev;

        REQUIRE( (SDS.shoc_cldfrac[offset] >= 0  || SDS.shoc_cldfrac[offset] < 1) );
        REQUIRE(SDS.shoc_ql2[offset] > 0);
        REQUIRE(SDS.shoc_ql[offset] >= 0);

        REQUIRE(std::abs(SDS.wqls[offset] ) < wqls_bound);
        REQUIRE(std::abs(SDS.wthv_sec[offset]) < wthv_sec_bound);
        REQUIRE(SDS.shoc_ql2[offset] < shoc_ql2_bound);
        REQUIRE(SDS.shoc_ql[offset] < shoc_ql_bound);

        // Now verify that the relationships in a strongly negative vertical
        //  velocity flux regime hold true, relative to a symmetric vertical
        //  velocity regime.
        if (s < shcol-1){

          if (SDS.shoc_cldfrac[offset] < 0.5){
            REQUIRE(SDS.shoc_cldfrac[offsets] < SDS.shoc_cldfrac[offset]);
          }
          else if (SDS.shoc_cldfrac[offset] > 0.5){
            REQUIRE(SDS.shoc_cldfrac[offsets] > SDS.shoc_cldfrac[offset]);
          }

          // In addition, in a positive skewness environment, the following
          //  should also be true

          // Grid mean liquid water decreased
          REQUIRE(SDS.shoc_ql[offsets] < SDS.shoc_ql[offset]);
          // if cloud present, verify liquid water and buoyancy flux is negative
          if (SDS.shoc_ql[offsets] > 0){
            REQUIRE(SDS.wqls[offsets] < SDS.wqls[offset]);
            REQUIRE(SDS.wqls[offsets] < 0);

            REQUIRE(SDS.wthv_sec[offsets] < SDS.wthv_sec[offset]);
            REQUIRE(SDS.wthv_sec[offsets] < 0);
          }

          // liquid water variance increased
          REQUIRE(SDS.shoc_ql2[offsets] > SDS.shoc_ql2[offset]);
        }
      }
    }

  }

  static void run_bfb()
  {
    ShocAssumedPdfData SDS_f90[] = {
      //              shcol, nlev, nlevi
      ShocAssumedPdfData(10, 71, 72),
      ShocAssumedPdfData(10, 12, 13),
      ShocAssumedPdfData(7,  16, 17),
      ShocAssumedPdfData(2, 7, 8),
    };

  static constexpr Int num_runs = sizeof(SDS_f90) / sizeof(ShocAssumedPdfData);

    // Generate random input data
    for (auto& d : SDS_f90) {
      d.randomize({ {d.thetal, {500, 700}} });

      // Generate grid as decreasing set of points.
      // Allows interpolated values to stay withing
      // reasonable range, avoiding errors in
      // shoc_assumed_pdf.
      const Real upper = 10;
      const Real lower = 0;
      for (Int s = 0; s < d.shcol; ++s) {
        for (Int k=0; k<d.nlevi; ++k) {
          const auto zi_k = upper - k*(upper-lower)/(d.nlevi-1);
          d.zi_grid[k+s*d.nlevi] = zi_k;

          if (k!=d.nlevi-1) {
            const auto zi_kp1 = upper - (k+1)*(upper-lower)/(d.nlevi-1);
            d.zt_grid[k+s*d.nlev] = 0.5*(zi_k + zi_kp1);
          }
        }
      }
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ShocAssumedPdfData SDS_cxx[] = {
      ShocAssumedPdfData(SDS_f90[0]),
      ShocAssumedPdfData(SDS_f90[1]),
      ShocAssumedPdfData(SDS_f90[2]),
      ShocAssumedPdfData(SDS_f90[3]),
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : SDS_f90) {
      // expects data in C layout
      shoc_assumed_pdf(d);
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      d.transpose<ekat::TransposeDirection::c2f>();
      // expects data in fortran layout
      shoc_assumed_pdf_f(d.shcol, d.nlev, d.nlevi, d.thetal, d.qw, d.w_field,
                         d.thl_sec, d.qw_sec, d.wthl_sec, d.w_sec, d.wqw_sec,
                         d.qwthl_sec, d.w3, d.pres, d.zt_grid, d.zi_grid,
                         d.shoc_cldfrac, d.shoc_ql, d.wqls, d.wthv_sec, d.shoc_ql2);
      d.transpose<ekat::TransposeDirection::f2c>();
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      ShocAssumedPdfData& d_f90 = SDS_f90[i];
      ShocAssumedPdfData& d_cxx = SDS_cxx[i];
      for (Int k = 0; k < d_f90.total(d_f90.wqls); ++k) {
        REQUIRE(d_f90.shoc_cldfrac[k] == d_cxx.shoc_cldfrac[k]);
        REQUIRE(d_f90.shoc_ql[k] == d_cxx.shoc_ql[k]);
        REQUIRE(d_f90.wqls[k] == d_cxx.wqls[k]);
        REQUIRE(d_f90.wthv_sec[k] == d_cxx.wthv_sec[k]);
        REQUIRE(d_f90.shoc_ql2[k] == d_cxx.shoc_ql2[k]);
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_assumed_pdf_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocAssumedPdf;

  TestStruct::run_property();
}

TEST_CASE("shoc_assumed_pdf_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocAssumedPdf;

  TestStruct::run_bfb();
}

} // namespace
