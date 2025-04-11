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
struct UnitWrap::UnitTest<D>::TestComputeShocTemp : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Int shcol    = 1;
    static constexpr Int nlev     = 3;

    // Test One
    // Given Exner value = 1 and cloud liquid = 0 everywhere
    //  verify that all points of absolute temperature (tabs)
    //  are exactly equal to liquid water potential temperature (thetal)

    // Inverse Exner value [-]
    Real inv_exner_first[nlev] = {1, 1, 1};
    // Liquid water potential temperature [K]
    Real thetal_first[nlev] = {300, 290, 280};
    // Liquid water mixing ratio [kg/kg]
    Real ql_first[nlev] = {0, 0, 0};

    ComputeShocTempData SDS(shcol, nlev);

    REQUIRE(SDS.shcol > 0);
    REQUIRE(SDS.nlev > 0);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.inv_exner[offset] = inv_exner_first[n];
        SDS.thetal[offset] = thetal_first[n];
        SDS.ql[offset] = ql_first[n];
      }
    }

    // Check that inputs are as expected
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.ql[offset] == 0);
        REQUIRE(SDS.inv_exner[offset] == 1);
        REQUIRE(SDS.thetal[offset] > 0);
      }
    }

    // Call the C++ implementation
    compute_shoc_temperature(SDS);

    // Require that absolute temperature is equal to thetal
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.tabs[offset] == SDS.thetal[offset]);
      }
    }

    // Test Two
    // Given profiles all with cloud liquid water greater than zero,
    //  AND inverse exner functions equal to 1, ensure that tabs is greater that
    //  absolute temperature is greater than (or equal to) liquid water potential temperature

    // Inverse Exner value [-]
    Real inv_exner_sec[nlev] = {1, 1, 1};
    // Liquid water potential temperature [K]
    Real thetal_sec[nlev] = {300, 290, 280};
    // Liquid water mixing ratio [kg/kg]
    Real ql_sec[nlev] = {3e-5, 1e-10, 1e-3};

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.inv_exner[offset] = inv_exner_sec[n];
        SDS.thetal[offset] = thetal_sec[n];
        SDS.ql[offset] = ql_sec[n];
      }
    }

    // Check that inputs are as expected
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.ql[offset] > 0);
        REQUIRE(SDS.inv_exner[offset] == 1);
        REQUIRE(SDS.thetal[offset] > 0);
      }
    }

    // Call the C++ implementation
    compute_shoc_temperature(SDS);

    // Require that absolute temperature is greather than thetal
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.tabs[offset] >= SDS.thetal[offset]);
      }
    }

    // Test Three
    // Given "realistic" atmospheric profiles with thetal increasing
    //  with height, as well as inv_exner function; verify that
    //  absolute temperature is decreasing with height.  Give
    //  all levels the same amount of cloud liquid water

    // Inverse Exner value [-]
    Real inv_exner_third[nlev] = {1.1, 1.5, 2};
    // Liquid water potential temperature [K]
    Real thetal_third[nlev] = {300, 350, 400};
    // Liquid water mixing ratio [kg/kg]
    Real ql_third[nlev] = {1e-5, 1e-5, 1e-5};

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.inv_exner[offset] = inv_exner_third[n];
        SDS.thetal[offset] = thetal_third[n];
        SDS.ql[offset] = ql_third[n];
      }
    }

    // Check that inputs are as expected
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.ql[offset] > 0);
        REQUIRE(SDS.inv_exner[offset] > 1);
        REQUIRE(SDS.thetal[offset] > 0);
      }
    }

    // Check that inputs are changing with height as expected
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 1; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        const auto offsetl = (n-1) + s * nlev;

        // Verify inverse exner and thetal are increasing with height
        REQUIRE(SDS.inv_exner[offset] > SDS.inv_exner[offsetl]);
        REQUIRE(SDS.thetal[offset] > SDS.thetal[offsetl]);
      }
    }

    // Call the C++ implementation
    compute_shoc_temperature(SDS);

    // Require that absolute temperature be less than thetal
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        REQUIRE(SDS.tabs[offset] < SDS.thetal[offset]);
      }
    }

    // Check that tabs is decreasing with height as expected
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 1; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        const auto offsetl = (n-1) + s * nlev;

        REQUIRE(SDS.tabs[offset] < SDS.tabs[offsetl]);
      }
    }

  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    ComputeShocTempData baseline_data[] = {
      //            shcol, nlev
      ComputeShocTempData(10, 71),
      ComputeShocTempData(10, 12),
      ComputeShocTempData(7,  16),
      ComputeShocTempData(2,   7)
    };

    // Generate random input data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    ComputeShocTempData cxx_data[] = {
      ComputeShocTempData(baseline_data[0]),
      ComputeShocTempData(baseline_data[1]),
      ComputeShocTempData(baseline_data[2]),
      ComputeShocTempData(baseline_data[3]),
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
      compute_shoc_temperature(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(ComputeShocTempData);
      for (Int i = 0; i < num_runs; ++i) {
        ComputeShocTempData& d_baseline = baseline_data[i];
        ComputeShocTempData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.tabs); ++k) {
          REQUIRE(d_baseline.tabs[k] == d_cxx.tabs[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (auto& d : cxx_data) {
        d.write(Base::m_fid);
      }
    }
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_compute_shoc_temperature_property", "shoc")
 {
   using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeShocTemp;

   TestStruct().run_property();
 }

TEST_CASE("shoc_compute_shoc_temperature_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeShocTemp;

  TestStruct().run_bfb();
}

} // namespace
