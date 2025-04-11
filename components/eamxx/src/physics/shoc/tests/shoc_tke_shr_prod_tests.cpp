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
struct UnitWrap::UnitTest<D>::TestShocShearProd : public UnitWrap::UnitTest<D>::Base {

  void run_property()
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

    // Call the C++ implementation
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

    // Call the C++ implementation
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

  void run_bfb()
  {
    auto engine = Base::get_engine();

    ComputeShrProdData baseline_data[] = {
      //            shcol, nlev
      ComputeShrProdData(10, 71, 72),
      ComputeShrProdData(10, 12, 13),
      ComputeShrProdData(7,  16, 17),
      ComputeShrProdData(2,   7, 8)
    };

    // Generate random input data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    ComputeShrProdData cxx_data[] = {
      ComputeShrProdData(baseline_data[0]),
      ComputeShrProdData(baseline_data[1]),
      ComputeShrProdData(baseline_data[2]),
      ComputeShrProdData(baseline_data[3]),
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
      compute_shr_prod(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ComputeShrProdData);
      for (Int i = 0; i < num_runs; ++i) {
        ComputeShrProdData& d_baseline = baseline_data[i];
        ComputeShrProdData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.sterm); ++k) {
          REQUIRE(d_baseline.sterm[k] == d_cxx.sterm[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (auto& d : cxx_data) {
        d.write(Base::m_fid);
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

  TestStruct().run_property();
}

TEST_CASE("shoc_tke_shr_prod_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocShearProd;

  TestStruct().run_bfb();
}

} // namespace
