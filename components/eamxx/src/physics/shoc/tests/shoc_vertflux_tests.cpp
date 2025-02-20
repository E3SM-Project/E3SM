#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "physics/share/physics_constants.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
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
struct UnitWrap::UnitTest<D>::TestCalcShocVertflux : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 4;
    static constexpr auto nlevi   = nlev + 1;

    // Propert tests for SHOC subroutine
    //   calc_shoc_vertflux

    //NOTE: This routine does not compute the vertical fluxes
    // for boundary points.  Input Grid values that are never used
    // are set to ZERO so as not to confuse the test data
    // that is actually used (and is also an implicit test to be
    // sure that operations are not done at the boundary).
    // Output values at boundary set to something other than
    // zero to check that they are not modified.

    // Define delta z on the nlevi grid [m]
    static constexpr Real dz_zi[nlevi] = {0, 100, 50, 20, 0};
    // Eddy coefficients on the nlevi grid [m2s-1]
    static constexpr Real tkh_zi[nlevi] = {0, 10, 10, 10, 0};
    // Invar (in this example we use potential temperature) [K]
    //   this variable should be on the nlev grid
    static constexpr Real invar[nlev] = {320, 315, 315, 316};
    // Initialize vertflux on the nlevi grid, set boundaries
    //   to large values to make sure they are not modified
    Real vertflux[nlevi] = {100, 0, 0, 0, 100};

    // NOTE: the input profile invar contains data to test conditions
    //   of 1) and unstable boundary layer, 2) well mixed
    //   layer and 3) conditionally stable layer.

    // Initialzie data structure for bridging to F90
    CalcShocVertfluxData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable.
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(nlevi - nlev == 1);
    REQUIRE(shcol > 0);

    // Fill in test data
    for(Int s = 0; s < shcol; ++s) {
      // First on the nlev grid
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.invar[offset] = invar[n];
      }

      // Now for data on the nlevi grid
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.tkh_zi[offset] = tkh_zi[n];
        SDS.dz_zi[offset] = dz_zi[n];
        SDS.vertflux[offset] = vertflux[n];
      }
    }

    // Check that the inputs make sense

    // Check to make sure that dz_zi are tkh_zi
    //  (outside of the boundaries) are greater than zero
    for(Int s = 0; s < shcol; ++s) {
      // do NOT check boundaries!
      for(Int n = 1; n < nlevi-1; ++n) {
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.dz_zi[offset] > 0);
        REQUIRE(SDS.tkh_zi[offset] > 0);
      }
    }

    // Call the C++ implementation
    calc_shoc_vertflux(SDS);

    // Check the results
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        // validate that the boundary points
        //   have NOT been modified
        if (n == 0 || n == nlevi-1){
          REQUIRE(SDS.vertflux[offset] == 100);
        }
        else{

          // Validate Downgradient assumption for
          //   various possible layers

          // conditionally stable layer
          if ((invar[n-1] - invar[n]) > 0){
            REQUIRE(SDS.vertflux[offset] < 0);
          }
          // well mixed layer
          if ((invar[n-1] - invar[n]) == 0){
            REQUIRE(SDS.vertflux[offset] == 0);
          }
          // unstable layer
          if ((invar[n-1] - invar[n]) < 0){
            REQUIRE(SDS.vertflux[offset] > 0);
          }
        }
      }
    }
  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    CalcShocVertfluxData SDS_baseline[] = {
      //               shcol, nlev, nlevi
      CalcShocVertfluxData(10, 71, 72),
      CalcShocVertfluxData(10, 12, 13),
      CalcShocVertfluxData(7,  16, 17),
      CalcShocVertfluxData(2, 7, 8),
    };

    // Generate random input data
    for (auto& d : SDS_baseline) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    CalcShocVertfluxData SDS_cxx[] = {
      CalcShocVertfluxData(SDS_baseline[0]),
      CalcShocVertfluxData(SDS_baseline[1]),
      CalcShocVertfluxData(SDS_baseline[2]),
      CalcShocVertfluxData(SDS_baseline[3]),
    };

    static constexpr Int num_runs = sizeof(SDS_baseline) / sizeof(CalcShocVertfluxData);

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : SDS_baseline) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      calc_shoc_vertflux(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        CalcShocVertfluxData& d_baseline = SDS_baseline[i];
        CalcShocVertfluxData& d_cxx = SDS_cxx[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.vertflux); ++k) {
          REQUIRE(d_baseline.vertflux[k] == d_cxx.vertflux[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        SDS_cxx[i].write(Base::m_fid);
      }
    }
  }

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("shoc_vertflux_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCalcShocVertflux;

  TestStruct().run_property();
}

TEST_CASE("shoc_vertflux_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestCalcShocVertflux;

  TestStruct().run_bfb();
}

} // namespace
