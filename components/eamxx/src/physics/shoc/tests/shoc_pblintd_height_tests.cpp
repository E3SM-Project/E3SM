#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestPblintdHeight : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr auto ustar_min = scream::shoc::Constants<Scalar>::ustar_min;
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
    PblintdHeightData SDS(shcol,nlev,nlev);

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

    // Call the C++ implementation
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
    for(Int s = 0; s < shcol; ++s) {
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

    // Call the C++ implementation
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

    // Call the C++ implementation
    pblintd_height(SDS);

    // Check that PBLH is zero (not modified) everywhere
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] == approx_zero);
    }

  } // run_property

  void run_bfb()
  {
    auto engine = Base::get_engine();

    Int npbl_rand = rand()%72 + 1;

    PblintdHeightData baseline_data[] = {
      PblintdHeightData(10, 72, 1),
      PblintdHeightData(10, 72, 72),
      PblintdHeightData(10, 72, npbl_rand),
      PblintdHeightData(10, 12, 1),
      PblintdHeightData(7, 16, 1),
      PblintdHeightData(2, 7, 1),
    };

    // Generate random input data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    PblintdHeightData cxx_data[] = {
      PblintdHeightData(baseline_data[0]),
      PblintdHeightData(baseline_data[1]),
      PblintdHeightData(baseline_data[2]),
      PblintdHeightData(baseline_data[3]),
      PblintdHeightData(baseline_data[4]),
      PblintdHeightData(baseline_data[5]),
    };

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : baseline_data) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      pblintd_height(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(PblintdHeightData);
      for (Int i = 0; i < num_runs; ++i) {
        PblintdHeightData& d_baseline = baseline_data[i];
        PblintdHeightData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.pblh); ++k) {
          REQUIRE(d_baseline.total(d_baseline.pblh) == d_cxx.total(d_cxx.pblh));
          REQUIRE(d_baseline.pblh[k] == d_cxx.pblh[k]);
          REQUIRE(d_baseline.total(d_baseline.pblh) == d_cxx.total(d_cxx.check));
          REQUIRE(d_baseline.check[k] == d_cxx.check[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (auto& d : cxx_data) {
        d.write(Base::m_fid);
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

  TestStruct().run_property();
}

TEST_CASE("pblintd_height_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdHeight;

  TestStruct().run_bfb();
}

} // empty namespace
