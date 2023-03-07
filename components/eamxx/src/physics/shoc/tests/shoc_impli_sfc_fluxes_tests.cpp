#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
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
struct UnitWrap::UnitTest<D>::TestImpSfcFluxes {

  static void run_property()
  {
    static constexpr Int shcol = 5;
    static constexpr Int num_tracer = 10;

    // Tests for the SHOC subroutine
    //   sfc_fluxes

    // TEST
    // Feed in several columns worth of data and make sure
    //  the output is consistent.

    // Surface density on the zi grid [kg/m3]
    static constexpr Real rho_zi_sfc[shcol] = {1.2, 1.0, 0.9, 1.1, 1.15};
    // Rdp value on zt grid [ms^2/kg], same for all columns
    static constexpr Real rdp_zt_sfc = 8.5e-3;
    // heat flux at surface [K m/s]
    static constexpr Real wthl_sfc[shcol] = {0.03, -0.03, 0.1, 0, -0.1};
    // moisture flux at surface [kg/kg m/s]
    static constexpr Real wqw_sfc[shcol] = {2e-5, 1e-6, 0, -2e-5, 1e-4};
    // TKE flux at the surface [m3/s3]
    static constexpr Real wtke_sfc[shcol] = {4e-2, 1e-3, -2e-3, 0, -1e-3};

    // Supply input values
    // liquid water potential temperature [K]
    static constexpr Real thetal_in = 300;
    // total water mixing ratio [kg/kg]
    static constexpr Real qw_in = 0.015;
    // turbulent kinetic energy [m2/s2]
    static constexpr Real tke_in = 0.4;

    // time step [s]
    static constexpr Real dtime = 300;

    // Input for tracer (no units)
    Real tracer_in[num_tracer];

    // Feed tracer random data from 1 to 1000
    for(Int t = 0; t < num_tracer; ++t) {
      tracer_in[t] = rand()% 1000 + 1;
    }

    // Initialize data structure for bridging to F90
    SfcFluxesData SDS(shcol, num_tracer, dtime);

    // Test that the inputs are reasonable.
    REQUIRE(SDS.shcol == shcol);
    REQUIRE(SDS.num_tracer == num_tracer);
    REQUIRE(shcol > 1);

    // Fill in test data, column only
    for(Int s = 0; s < shcol; ++s) {
      SDS.rho_zi_sfc[s] = rho_zi_sfc[s];
      SDS.rdp_zt_sfc[s] = rdp_zt_sfc;
      SDS.wthl_sfc[s] = wthl_sfc[s];
      SDS.wqw_sfc[s] = wqw_sfc[s];
      SDS.wtke_sfc[s] = wtke_sfc[s];

      SDS.thetal[s] = thetal_in;
      SDS.qw[s] = qw_in;
      SDS.tke[s] = tke_in;

      for (Int t = 0; t < num_tracer; ++t){
        const auto offset = t + s * num_tracer;
        SDS.wtracer[offset] = tracer_in[t];
        // Feed tracer flux random data from -100 to 100
        //   note this is different for every point
        SDS.wtracer_sfc[offset] = rand()% 200 + (-100);
      }

    }

    // Check that the inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE( (SDS.thetal[s] > 150 && SDS.thetal[s] < 350) );
      REQUIRE( (SDS.qw[s] > 0.0001 && SDS.qw[s] < 0.05) );
      REQUIRE( (SDS.tke[s] > 0 && SDS.tke[s] < 10) );
      REQUIRE( (SDS.rdp_zt_sfc[s] > 0 && SDS.rdp_zt_sfc[s] < 1) );
      REQUIRE( (SDS.rho_zi_sfc[s] > 0 && SDS.rho_zi_sfc[s] < 2) );
      REQUIRE(std::abs(SDS.wthl_sfc[s]) < 1);
      REQUIRE(std::abs(SDS.wqw_sfc[s]) < 1e-3);
      REQUIRE(std::abs(SDS.wtke_sfc[s]) < 0.1);
      REQUIRE(SDS.dtime > 0);
    }

    // Call the fortran implementation
    sfc_fluxes(SDS);

    // Verify that output is reasonable
    for(Int s = 0; s < shcol; ++s) {

      // Verify output falls within reasonable bounds
      REQUIRE( (SDS.thetal[s] > 150 && SDS.thetal[s] < 350) );
      REQUIRE( (SDS.qw[s] > 0.0001 && SDS.qw[s] < 0.05) );
      REQUIRE( (SDS.tke[s] > 0 && SDS.tke[s] < 10) );

      // Based on surface flux input, make sure that
      //  temperature, moisture, and tke all have output
      //  that is expected with respect to the input

      // Check temperature
      if (wthl_sfc[s] > 0){
        REQUIRE(SDS.thetal[s] > thetal_in);
      }
      else if (wthl_sfc[s] < 0){
        REQUIRE(SDS.thetal[s] < thetal_in);
      }
      else{
        REQUIRE(SDS.thetal[s] == thetal_in);
      }

      // Check moisture
      if (wqw_sfc[s] > 0){
        REQUIRE(SDS.qw[s] > qw_in);
      }
      else if (wqw_sfc[s] < 0){
        REQUIRE(SDS.qw[s] < qw_in);
      }
      else{
        REQUIRE(SDS.qw[s] == qw_in);
      }

      // Check TKE
      if (wtke_sfc[s] > 0){
        REQUIRE(SDS.tke[s] > tke_in);
      }
      else if (wtke_sfc[s] < 0){
        REQUIRE(SDS.tke[s] < tke_in);
      }
      else{
        REQUIRE(SDS.tke[s] == tke_in);
      }

      // Check tracer
      for (Int t = 0; t < num_tracer; ++t){
        const auto offset = t + s * num_tracer;
        if (SDS.wtracer_sfc[offset] > 0){
          REQUIRE(SDS.wtracer[offset] > tracer_in[t]);
        }
        else if (SDS.wtracer_sfc[offset] < 0){
          REQUIRE(SDS.wtracer[offset] < tracer_in[t]);
        }
        else{
          REQUIRE(SDS.wtracer[offset] == tracer_in[t]);
        }
      }

    }

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

TEST_CASE("shoc_imp_sfc_fluxes_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestImpSfcFluxes;

  TestStruct::run_property();
}

TEST_CASE("shoc_imp_sfc_fluxes_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestImpSfcFluxes;

  TestStruct::run_bfb();
}

} // namespace
