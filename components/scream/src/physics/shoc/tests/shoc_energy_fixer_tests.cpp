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
struct UnitWrap::UnitTest<D>::TestShocEnergyFixer {

  static void run_property()
  {
    static constexpr Real gravit  = scream::physics::Constants<Real>::gravit;
    static constexpr Real Cpair   = scream::physics::Constants<Real>::Cpair;
    static constexpr Int shcol    = 1;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;

    // Tests for the SHOC function
    //     shoc_energy_fixer
    
    // TEST ONE and TWO
    // No energy change and surface flux test.  
    //  In one column, given inputs where the energy has not changed
    //  from timestep to timestep, verify that the host model dry static
    //  energy also has not changed. 
    
    // In other, give positive surface fluxes, 
    //  all other information being the same.  Verify that energy was added 
    //  in to the system and verify that energy was added into the system AND
    //  verify that energy was only adjust at the levels where TKE is greater
    //  than zero.  

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
    static constexpr Real tke[nlev] = {0, 0, 0.3, 0.4, 0.1};
   //  Pressure at interface [Pa]
    static constexpr Real pint[nlevi] = {50000, 60000, 70000, 80000, 90000, 100000};
    // Define integrated static energy, kinetic energy, water vapor,
    //  and liquid water respectively
    static constexpr Real se = 200.0;
    static constexpr Real ke = 150.0;
    static constexpr Real wv = 0.5;
    static constexpr Real wl = 0.1;
    // Define surface sensible heat flux [K m/s]
    static constexpr Real wthl_sfc = 0.5;
    // Define surface total water flux [kg/kg m/s]
    static constexpr Real wqw_sfc = 0.001;
    Real host_dse_input[nlev];
    Real zt_grid[nlev];

    // Initialzie data structure for bridgeing to F90
    SHOCEnergyfixerData SDS(shcol, nlev, nlevi, dtime, nadv);

    // Test that the inputs are reasonable.
    // for this test we need exactly two columns
    REQUIRE( (SDS.shcol() == shcol && SDS.nlev() == nlev && SDS.nlevi() && SDS.dtime == dtime && SDS.nadv == nadv) );
    // Want exactly two columns for this case
    REQUIRE(shcol == 2);
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
      SDS.wthl_sfc[s] = s*wthl_sfc;
      SDS.wqw_sfc[s] = s*wqw_sfc;

      // Fill in test data on zt_grid.
      for(Int n = 0; n < nlev; ++n) {
      const auto offset = n + s * nlev;

        // For zt grid, set as midpoint of zi grid
        SDS.zt_grid[offset] = zt_grid[n];
        // For pdel, compute from pint
        SDS.pdel[offset] = pint[n+1]-pint[n];
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
        REQUIRE(SDS.pdel[offset] > 0);
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

    // Call the fortran implementation
    shoc_energy_fixer(SDS);

    // Check test
    // Verify that the dry static energy has not changed
    for(Int s = 0; s < shcol; ++s) {
      for (Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        
        REQUIRE(SDS.host_dse[offset] == host_dse_input[n]);
      }
    }
    
    // TEST TWO

    // For first column verify that total energies are the same
    // REQUIRE(SDS.te_a[0] == SDS.te_b[0]);

    // Verify that second column "before" energy is greater than
    //  the first column, since here we have active surface fluxes
    //  REQUIRE(SDS.te_b[1] > SDS.te_b[0]);
  }

  static void run_bfb()
  {
    // TODO
  }
};

}  // namespace unit_test
}  // namespace shoc
}  // namespace scream

namespace {

TEST_CASE("shoc_energy_fixer_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyFixer;

  TestStruct::run_property();
}

TEST_CASE("shoc_energy_fixer_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocEnergyFixer;

  TestStruct::run_bfb();
}

} // namespace
