#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestUpdatePrognosticsImplicit {

  static void run_property()
  {
    static constexpr Int shcol    = 5;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;
    static constexpr Int num_tracer = 10;
  
    // Tests for the subroutine update_prognostics_implicit

    // Timestep [s]
    static constexpr Real dtime = 300;
 
    // Define the heights on the zi grid [m]
    static constexpr Real zi_grid[nlevi] = {900, 500, 150, 90, 50, 0};
    // Define the eddy vicosity for heat and momentum [m2/s]
    static constexpr Real tkh[nlev] = {3, 10, 50, 30, 20};
    // Define air density on midpoint grid [kg/m3]
    static constexpr Real rho_zt[nlev] = {0.7, 0.8, 0.85, 0.9, 1.0};

    // heat flux at surface [K m/s], COLUMN ONLY variables
    static constexpr Real wthl_sfc[shcol] = {0.03, -0.03, 0.1, 0, -0.1};
    // moisture flux at surface [kg/kg m/s]
    static constexpr Real wqw_sfc[shcol] = {2e-5, 1e-6, 0, -2e-5, 1e-4};
    // TKE flux at the surface [m3/s3]
    static constexpr Real wtke_sfc[shcol] = {4e-2, 1e-3, -2e-3, 0, -1e-3};
    // Surface moment flux, zonal direction [m3/s3]
    static constexpr Real uw_sfc[shcol] = {0.03, -0.03, 0.1, 0, -0.1};
    // Surface moment flux, meridional direction [m3/s3]
    static constexpr Real vw_sfc[shcol] = {-0.01, -0.01, 0.3, 0, -0.3};
    
    // IN/OUT variables, PROFILES
    // Define the liquid water potential temperature [K]
    static constexpr Real thetal_in[nlev] = {310, 307, 302, 302, 303};
    // Define the total water mixing ratio [kg/kg]
    static constexpr Real qw_in[nlev] = {1e-2, 1.2e-2, 1.5e-2, 1.5e-2, 1.4e-2};
    // Define the zonal wind [m/s]
    static constexpr Real u_wind_in[nlev] = {4, 4, 2, 0, -1};
    // define the meridional wind [m/s]
    static constexpr Real v_wind_in[nlev] = {-2, -2, 1, 3, 0};
    // Define the TKE [m2/s2]
    static constexpr Real tke_in[nlev] = {0.2, 0.3, 0.5, 0.4, 0.1};

    // Input for tracer (no units)
    Real tracer_in[shcol,nlev,num_tracer];

    // Compute needed grid information from zi_grid
    // Grid stuff to compute based on zi_grid
    Real zt_grid[nlev];
    Real dz_zt[nlev];
    Real dz_zi[nlevi];
    Real w_sec[nlev];
    // Compute heights on midpoint grid
    for(Int n = 0; n < nlev; ++n) {
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
      dz_zt[n] = zi_grid[n] - zi_grid[n+1];
      if (n == 0){
        dz_zi[n] = 0;
      }
      else{
        dz_zi[n] = zt_grid[n-1] - zt_grid[n];
      }
    }
    // set upper condition for dz_zi
    dz_zi[nlevi-1] = zt_grid[nlev-1]; 
    
    // Initialize data structure for bridging to F90
    UpdatePrognosticsImplicitData SDS(shcol,nlev,nlevi,num_tracer,dtime);     

  } // run_property

  static void run_bfb()
  {
    UpdatePrognosticsImplicitData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(UpdatePrognosticsImplicitData);

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    UpdatePrognosticsImplicitData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      update_prognostics_implicit(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      update_prognostics_implicit_f(d.shcol, d.nlev, d.nlevi, d.num_tracer, d.dtime, d.dz_zt, d.dz_zi, d.rho_zt, d.zt_grid, d.zi_grid, d.tk, d.tkh, d.uw_sfc, d.vw_sfc, d.wthl_sfc, d.wqw_sfc, d.wtracer_sfc, d.thetal, d.qw, d.tracer, d.tke, d.u_wind, d.v_wind);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      UpdatePrognosticsImplicitData& d_f90 = f90_data[i];
      UpdatePrognosticsImplicitData& d_cxx = cxx_data[i];
      for (Int k = 0; k < d_f90.total(d_f90.thetal); ++k) {
        REQUIRE(d_f90.total(d_f90.thetal) == d_cxx.total(d_cxx.thetal));
        REQUIRE(d_f90.thetal[k] == d_cxx.thetal[k]);
        REQUIRE(d_f90.total(d_f90.thetal) == d_cxx.total(d_cxx.qw));
        REQUIRE(d_f90.qw[k] == d_cxx.qw[k]);
        REQUIRE(d_f90.total(d_f90.thetal) == d_cxx.total(d_cxx.tke));
        REQUIRE(d_f90.tke[k] == d_cxx.tke[k]);
        REQUIRE(d_f90.total(d_f90.thetal) == d_cxx.total(d_cxx.u_wind));
        REQUIRE(d_f90.u_wind[k] == d_cxx.u_wind[k]);
        REQUIRE(d_f90.total(d_f90.thetal) == d_cxx.total(d_cxx.v_wind));
        REQUIRE(d_f90.v_wind[k] == d_cxx.v_wind[k]);
      }
      for (Int k = 0; k < d_f90.total(d_f90.tracer); ++k) {
        REQUIRE(d_f90.total(d_f90.tracer) == d_cxx.total(d_cxx.tracer));
        REQUIRE(d_f90.tracer[k] == d_cxx.tracer[k]);
      }

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("update_prognostics_implicit_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestUpdatePrognosticsImplicit;

  TestStruct::run_property();
}

TEST_CASE("update_prognostics_implicit_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestUpdatePrognosticsImplicit;

  TestStruct::run_bfb();
}

} // empty namespace
