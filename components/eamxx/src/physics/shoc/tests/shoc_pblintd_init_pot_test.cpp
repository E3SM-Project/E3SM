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
struct UnitWrap::UnitTest<D>::TestPblintdInitPot : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Int shcol    = 2;
    static constexpr Int nlev     = 1;

    // Tests for the SHOC function:
    //   pblintd_init_pot

    // FIRST TEST
    //  Dry atmosphere test.  Verify that in a dry atmosphere
    //  with no water vapor or liquid water loading that the output
    //  thv (virtual potential temperature) is equal to the input
    //  thl (liquid water potential temperature).

    // Define the liquid water potential temperature [K]
    Real thl_dry[shcol] = {300, 295};
    // Define the water vapor [kg/kg]
    Real qv = 0.0;
    // Define the liquid water mixing ratio [kg/kg]
    Real ql = 0.0;

    // Initialize data structure for bridging to F90
    PblintdInitPotData SDS(shcol, nlev);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol == shcol);
    REQUIRE(SDS.nlev == nlev);
    // for this test require two columns
    REQUIRE(shcol == 2);

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.thl[offset] = thl_dry[s];
        SDS.ql[offset] = ql;
        SDS.q[offset] = qv;
      }
    }

     // Check that the inputs are expected
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        REQUIRE(SDS.thl[offset] > 0);
        // For this test require ql and q be zero
        REQUIRE(SDS.q[offset] == 0);
        REQUIRE(SDS.ql[offset] == 0);
      }
    }

    // call the C++ implementation
    pblintd_init_pot(SDS);

    // Check the result.
    // Verify that virtual potential temperature is idential
    //  to potential temperature
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        REQUIRE(SDS.thl[offset] == SDS.thv[offset]);
      }
    }

    // SECOND TEST
    // For two parcels with identical inputs, but one with condensate
    //  loading, verify that the two parcels return different answers
    //  and fall within reasonable bounds

    // Define the liquid water potential temperature [K]
    Real thl_parcel[shcol] = {290, 290};
    // Define the water vapor [kg/kg]
    Real qv_parcel = 0.01;
    // Define the liquid water mixing ratio [kg/kg]
    Real ql_parcel[shcol] = {0, 4e-4};

    // Fill in test data on zt_grid.
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.thl[offset] = thl_parcel[s];
        // convert the following to kg/kg
        SDS.ql[offset] = ql_parcel[s];
        SDS.q[offset] = qv_parcel;
      }
    }

    // Check that the inputs are expected
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        // Make sure top level dz_zi value is zero
        REQUIRE(SDS.thl[offset] > 0.0);
        // Make sure input is within reasonable bounds
        REQUIRE(SDS.q[offset] < 0.1);
        REQUIRE(SDS.ql[offset] < 0.1);
      }
    }

    // Call the C++ implementation
    pblintd_init_pot(SDS);

    // Check test
    // Verify that column with condensate loading
    //   results in a lower virtual potential temperature
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Verify output falls within reasonable bounds
        // Parcels should be of greater magnitude than
        //  the input potential temperature, but within reason
        REQUIRE(SDS.thv[offset] > SDS.thl[offset]);
        REQUIRE(SDS.thv[offset] < SDS.thl[offset]+10);
        // Get value corresponding to next column
        if (s < shcol-1) {
          const auto offsets = n + (s+1) * nlev;
          REQUIRE(SDS.thv[offsets] != SDS.thv[offset]);
        }
      }
    }

  } // run_property

  void run_bfb()
  {
    auto engine = Base::get_engine();

    PblintdInitPotData pblintd_init_pot_data_baseline[] = {
      //                     shcol, nlev
      PblintdInitPotData(36,  72),
      PblintdInitPotData(72,  72),
      PblintdInitPotData(128, 72),
      PblintdInitPotData(256, 72),
    };

    for (auto& d : pblintd_init_pot_data_baseline) {
      d.randomize(engine);
    }

    PblintdInitPotData pblintd_init_pot_data_cxx[] = {
      PblintdInitPotData(pblintd_init_pot_data_baseline[0]),
      PblintdInitPotData(pblintd_init_pot_data_baseline[1]),
      PblintdInitPotData(pblintd_init_pot_data_baseline[2]),
      PblintdInitPotData(pblintd_init_pot_data_baseline[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : pblintd_init_pot_data_baseline) {
        d.read(Base::m_fid);
      }
    }

    for (auto& d : pblintd_init_pot_data_cxx) {
      pblintd_init_pot(d);
    }

    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(pblintd_init_pot_data_baseline) / sizeof(PblintdInitPotData);
      for (Int i = 0; i < num_runs; ++i) {
        Int shcol = pblintd_init_pot_data_cxx[i].shcol;
        Int nlev  = pblintd_init_pot_data_cxx[i].nlev;
        for (Int j = 0; j < shcol; ++j ) {
          for (Int k = 0; k < nlev; ++k) {
            REQUIRE(pblintd_init_pot_data_baseline[i].thv[j*k] == pblintd_init_pot_data_cxx[i].thv[j*k]);
          }
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (auto& d : pblintd_init_pot_data_cxx) {
        d.write(Base::m_fid);
      }
    }
  }

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("shoc_pblintd_init_pot_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdInitPot;

  TestStruct().run_property();
}

TEST_CASE("shoc_pblintd_init_pot_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdInitPot;

  TestStruct().run_bfb();
}

}  // namespace
