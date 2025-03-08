#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/eamxx_types.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

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
struct UnitWrap::UnitTest<D>::TestShocIsotropicTs : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Real maxiso = scream::shoc::Constants<Real>::maxiso;
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 1;

    // Tests for the subroutine isotropic_ts
    //   in the SHOC TKE module.

    // For this routine all inputs are on midpoint grid and
    //  there are no vertical derivatives, therefore we will
    //  consider one vertical level per test.  Each column will
    //  be loaded up with a different test.

    // FIRST TEST
    // Stability test.  Verify that two points with the same inputs
    //  but one with brunt > 0 and the other with brunt < 0 that the
    //  timescale with brunt > 0 is less.

    // Integrated brunt vaisalla
    static constexpr Real brunt_int_st = 0.1;
    // TKE [m2/s2]
    static constexpr Real tke_st = 0.4;
    // Dissipation rate [m2/s3]
    static constexpr Real diss_st = 0.1;
    // Brunt Vaisalla frequency [/s]
    static constexpr Real brunt_st[shcol] = {-4e-3, 4e-3};

    // Initialzie data structure for bridgeing to F90
    IsotropicTsData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev) );
    REQUIRE(shcol > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      // Column only input
      SDS.brunt_int[s] = brunt_int_st;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.tke[offset] = tke_st;
        SDS.a_diss[offset] = diss_st;
        SDS.brunt[offset] = brunt_st[s];
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Should be greater than zero
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.a_diss[offset] > 0);
      }
    }

    // Call the C++ implementation
    isotropic_ts(SDS);

    // Check that output falls within reasonable bounds
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.isotropy[offset] <= maxiso);
        REQUIRE(SDS.isotropy[offset] >= 0);
      }
    }

    // Check to make sure that column with positive
    //  brunt vaisalla frequency is smaller
    for(Int s = 0; s < shcol-1; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Get value corresponding to next column
        const auto offsets = n + (s+1) * nlev;
        if(SDS.brunt[offset] < 0 & SDS.brunt[offsets] > 0){
          REQUIRE(SDS.isotropy[offset] > SDS.isotropy[offsets]);
        }
        else{
          REQUIRE(SDS.isotropy[offset] < SDS.isotropy[offsets]);
        }
      }
    }

    // SECOND TEST
    // Dissipation test.  Verify that two points with the same inputs
    //  but one with smaller dissipation rate, that that point has
    //  a higher isotropy timescale.

    // Integrated brunt vaisalla
    static constexpr Real brunt_int_diss = 0.1;
    // TKE [m2/s2]
    static constexpr Real tke_diss = 0.4;
    // Dissipation rate [m2/s3]
    static constexpr Real diss_diss[shcol] = {0.1, 0.2};
    // Brunt Vaisalla frequency [/s]
    static constexpr Real brunt_diss = 4e-3;

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      SDS.brunt_int[s] = brunt_int_diss;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.tke[offset] = tke_diss;
        SDS.a_diss[offset] = diss_diss[s];
        SDS.brunt[offset] = brunt_diss;
      }
    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Should be greater than zero
        REQUIRE(SDS.tke[offset] > 0);
        REQUIRE(SDS.a_diss[offset] > 0);
      }
    }

    // Call the C++ implementation
    isotropic_ts(SDS);

    // Check that output falls within reasonable bounds
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.isotropy[offset] <= maxiso);
        REQUIRE(SDS.isotropy[offset] >= 0);
      }
    }

    // Check to make sure that column with positive
    //  brunt vaisalla frequency is smaller
    for(Int s = 0; s < shcol-1; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Get value corresponding to next column
        const auto offsets = n + (s+1) * nlev;
        if(SDS.a_diss[offset] < SDS.a_diss[offsets]){
          REQUIRE(SDS.isotropy[offset] > SDS.isotropy[offsets]);
        }
        else{
          REQUIRE(SDS.isotropy[offset] < SDS.isotropy[offsets]);
        }
      }
    }
  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    IsotropicTsData baseline_data[] = {
      //            shcol, nlev
      IsotropicTsData(10, 71),
      IsotropicTsData(10, 12),
      IsotropicTsData(7,  16),
      IsotropicTsData(2,   7)
    };

    // Generate random input data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    IsotropicTsData cxx_data[] = {
      IsotropicTsData(baseline_data[0]),
      IsotropicTsData(baseline_data[1]),
      IsotropicTsData(baseline_data[2]),
      IsotropicTsData(baseline_data[3]),
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
      isotropic_ts(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(IsotropicTsData);
      for (Int i = 0; i < num_runs; ++i) {
        IsotropicTsData& d_baseline = baseline_data[i];
        IsotropicTsData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.isotropy); ++k) {
          REQUIRE(d_baseline.isotropy[k] == d_cxx.isotropy[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (auto& d : cxx_data) {
        d.write(Base::m_fid);
      }
    }
  }//run_bfb

};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_tke_isotropic_ts_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocIsotropicTs;

  TestStruct().run_property();
}

TEST_CASE("shoc_tke_isotropic_ts_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocIsotropicTs;

  TestStruct().run_bfb();
}

} // namespace
