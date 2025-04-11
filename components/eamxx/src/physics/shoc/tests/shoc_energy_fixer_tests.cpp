#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "physics/share/physics_constants.hpp"
#include "shoc_constants.hpp"
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
struct UnitWrap::UnitTest<D>::TestShocEnergyFixer : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Real gravit  = scream::physics::Constants<Real>::gravit;
    static constexpr Real Cpair   = scream::physics::Constants<Real>::Cpair;
    static constexpr Real mintke  = scream::shoc::Constants<Real>::mintke;
    static constexpr Int shcol    = 3;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the SHOC function
    //     shoc_energy_fixer

    // TESTs ONE through THREE
    // No energy change and surface flux test.
    //  In one column, given inputs where the energy has not changed
    //  from timestep to timestep, verify that the host model dry static
    //  energy also has not changed.

    // In second column, give positive surface fluxes,
    //  all other information being the same.  Verify that energy was added
    //  in to the system and verify that energy was added into the system AND
    //  verify that energy was only adjust at the levels where TKE is greater
    //  than zero.

    // In third column, give negative surface fluxes and verify energy
    //  was removed from the system.

    // Timestep [s]
    static constexpr Real dtime = 300;
    // Number of macmic steps
    static constexpr Int nadv = 2;
    // Air density [km/m3]
    static constexpr Real rho_zt[nlev] = {0.4, 0.6, 0.7, 0.9, 1.0};
    // Interface heights [m]
    static constexpr Real zi_grid[nlevi] = {11000, 7500, 5000, 3000, 1500, 0};
    // Host model temperture [K]
    static constexpr Real host_temp_input[nlev] = {250, 275, 285, 290, 300};
    // Define TKE inputs
    static constexpr Real tke[nlev] = {mintke, mintke, 0.3, 0.4, 0.1};
   //  Pressure at interface [Pa]
    static constexpr Real pint[nlevi] = {50000, 60000, 70000, 80000, 90000, 100000};
    // Define integrated static energy, kinetic energy, water vapor,
    //  and liquid water respectively
    static constexpr Real se = 200;
    static constexpr Real ke = 150;
    static constexpr Real wv = 0.5;
    static constexpr Real wl = 0.1;
    // Define surface sensible heat flux [K m/s]
    static constexpr Real wthl_sfc[shcol] = {0, 0.5, -0.5};
    // Define surface total water flux [kg/kg m/s]
    static constexpr Real wqw_sfc[shcol] = {0, 0.001, -0.001};
    Real host_dse_input[nlev];
    Real zt_grid[nlev];

    // Initialize data structure for bridging to F90
    ShocEnergyFixerData SDS(shcol, nlev, nlevi, dtime, nadv);

    // Test that the inputs are reasonable.
    // for this test we need exactly two columns
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi && SDS.dtime == dtime && SDS.nadv == nadv) );
    // Want exactly three columns for this case
    REQUIRE(shcol == 3);
    REQUIRE(nlevi == nlev+1);

    // compute host model dry static energy
    for(Int n = 0; n < nlev; ++n) {
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
      host_dse_input[n] = Cpair*host_temp_input[n]+gravit*zt_grid[n];
    }

    for(Int s = 0; s < shcol; ++s) {
      // Set before and after integrals equal
      SDS.se_a[s] = se;
      SDS.se_b[s] = se;
      SDS.ke_a[s] = ke;
      SDS.ke_b[s] = ke;
      SDS.wv_a[s] = wv;
      SDS.wv_b[s] = wv;
      SDS.wl_a[s] = wl;
      SDS.wl_b[s] = wl;

      // Make first column be zero for the surface fluxes
      SDS.wthl_sfc[s] = wthl_sfc[s];
      SDS.wqw_sfc[s] = wqw_sfc[s];

      // Fill in test data on zt_grid.
      for(Int n = 0; n < nlev; ++n) {
      const auto offset = n + s * nlev;

        // For zt grid, set as midpoint of zi grid
        SDS.zt_grid[offset] = zt_grid[n];
        SDS.rho_zt[offset] = rho_zt[n];
        SDS.tke[offset] = tke[n];
        SDS.host_dse[offset] = host_dse_input[n];
      }

      // Fill in test data on zi_grid.
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.pint[offset] = pint[n];
        SDS.zi_grid[offset] = zi_grid[n];
      }
    }

    // Check that the inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;

        REQUIRE(SDS.zt_grid[offset] >= 0);
        REQUIRE(SDS.rho_zt[offset] > 0);
        REQUIRE(SDS.tke[offset] >= 0);

        // Check that heights increase upward
        if (n > nlev-1){
          REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
        }
      }
      for (Int n = 0; n < nlevi; ++n){
        const auto offset = n + s * nlevi;

        REQUIRE(SDS.zi_grid[offset] >= 0);

        // Check that heights increase upward
        if (n > nlevi-1){
          REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0);
        }
      }
    }

    // Call the C++ implementation
    shoc_energy_fixer(SDS);

    // Check test
    // Verify that the dry static energy has not changed if surface
    //  fluxes are zero, else verify host_dse has increased
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;

        // If TKE is at minimum threshold then make sure that
        //  host_dse has not been modified
        if (SDS.tke[offset] == mintke){
          REQUIRE(SDS.host_dse[offset] == host_dse_input[n]);
        }
        else{
          if (SDS.wthl_sfc[s] == 0){
            REQUIRE(SDS.wqw_sfc[s] == 0); // verify input
            REQUIRE(SDS.host_dse[offset] == host_dse_input[n]);
          }
          else if (SDS.wthl_sfc[s] < 0){
            REQUIRE(SDS.wqw_sfc[s] < 0); // verify input
            REQUIRE(SDS.host_dse[offset] < host_dse_input[n]);
          }
          else {
            REQUIRE(SDS.wqw_sfc[s] > 0);
            REQUIRE(SDS.host_dse[offset] > host_dse_input[n]);
          }
        }

      }
    }

    // TEST FOUR
    // Energy loss, gain test.
    // Now set surface fluxes to zero and set all *_a scalar arrays
    //  to be less/more than the *_b arrays.  This will signify that SHOC
    //  lost/gained energy during integration and thus host_dse should be
    //  INCREASED/DECREASED at all levels where TKE is sufficient

    // Define by how much each energy integral array will
    //  be perturbed
    static constexpr Real se_gainloss = 2;
    static constexpr Real ke_gainloss = 1;
    static constexpr Real wv_gainloss = 0.01;
    static constexpr Real wl_gainloss = 0.001;

    // For this test set to zero
    static constexpr Real wthl_sfc_gainloss = 0;
    static constexpr Real wqw_sfc_gainloss = 0;

    // Load up the data
    // This will alternate between a positive and negative
    Real gainloss_fac[shcol];
    // Initialize value
    gainloss_fac[0] = 1;
    for(Int s = 0; s < shcol; ++s) {
      SDS.wthl_sfc[s] = wthl_sfc_gainloss;
      SDS.wqw_sfc[s] = wqw_sfc_gainloss;

      REQUIRE( (gainloss_fac[s] == 1 || gainloss_fac[s] == -1) );

      SDS.se_a[s] = se+gainloss_fac[s]*se_gainloss;
      SDS.ke_a[s] = ke+gainloss_fac[s]*ke_gainloss;
      SDS.wv_a[s] = wv+gainloss_fac[s]*wv_gainloss;
      SDS.wl_a[s] = wl+gainloss_fac[s]*wl_gainloss;

      // Alternate between positive and negative loss,
      //  change sign for next column
      if (s < shcol-1){
        gainloss_fac[s+1] = gainloss_fac[s]*-1;
      }
    }

    // Call the C++ implementation
    shoc_energy_fixer(SDS);

    // Verify the result
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;

        // If TKE is at minimum threshold then make sure that
        //  host_dse has not been modified
        if (SDS.tke[offset] == mintke){
          REQUIRE(SDS.host_dse[offset] == host_dse_input[n]);
        }
        else{
          // If the system gained energy, make sure that host_dse
          //  has decreased to compensate for this
          if (gainloss_fac[s] > 0){
            REQUIRE(SDS.host_dse[offset] < host_dse_input[n]);
          }
          // Else, the vice versa of above
          else{
            REQUIRE(SDS.host_dse[offset] > host_dse_input[n]);
          }
        }

      }
    }

  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    ShocEnergyFixerData SDS_baseline[] = {
      //               shcol, nlev, nlevi, dtime, nadv
      ShocEnergyFixerData(10, 71, 72, 300, 2),
      ShocEnergyFixerData(10, 12, 13, 100, 10),
      ShocEnergyFixerData(7,  16, 17, 50, 1),
      ShocEnergyFixerData(2, 7, 8, 5, 5),
    };

    // Generate random input data
    for (auto& d : SDS_baseline) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    ShocEnergyFixerData SDS_cxx[] = {
      ShocEnergyFixerData(SDS_baseline[0]),
      ShocEnergyFixerData(SDS_baseline[1]),
      ShocEnergyFixerData(SDS_baseline[2]),
      ShocEnergyFixerData(SDS_baseline[3]),
    };

    static constexpr Int num_runs = sizeof(SDS_baseline) / sizeof(ShocEnergyFixerData);

    // Assume all data is in C layout

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : SDS_baseline) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (auto& d : SDS_cxx) {
      shoc_energy_fixer(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        ShocEnergyFixerData& d_baseline = SDS_baseline[i];
        ShocEnergyFixerData& d_cxx = SDS_cxx[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.host_dse); ++k) {
          REQUIRE(d_baseline.host_dse[k] == d_cxx.host_dse[k]);
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

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_energy_fixer_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyFixer;

  TestStruct().run_property();
}

TEST_CASE("shoc_energy_fixer_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyFixer;

  TestStruct().run_bfb();
}

} // namespace
