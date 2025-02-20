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
struct UnitWrap::UnitTest<D>::TestUpdatePrognosticsImplicit : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    static constexpr Int shcol    = 5;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;
    static constexpr Int num_tracer = 3;

    // Tests for the subroutine update_prognostics_implicit

    // TEST details
    //  Given a set of realistic model profiles verify that energy/water
    //    and tracers are conserved and that output falls within
    //    reasonable bounds.

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

    // establish reasonable bounds for checking input/output
    static constexpr Real thl_lbound = 200; // [K]
    static constexpr Real thl_ubound = 350; // [K]
    static constexpr Real qw_lbound = 1e-4; // [kg/kg]
    static constexpr Real qw_ubound = 5e-2; // [kg/kg]
    static constexpr Real tke_lbound = 0; // [m2/s2]
    static constexpr Real tke_ubound = 5; // [m2/s2]
    static constexpr Real rho_lbound = 0; // [kg/m3]
    static constexpr Real rho_ubound = 1.5; // [kg/m3]
    static constexpr Real wind_bounds = 5; // [m/s]

    static constexpr Real thresh_check = 1e-4;

    // Input for tracer (no units)
    Real tracer_in[shcol][nlev][num_tracer];

    // Define Integrals for energy/water conservation checking
    Real qw_int_b[shcol];
    Real qw_int_a[shcol];
    Real thl_int_b[shcol];
    Real thl_int_a[shcol];
    Real trc_int_b[shcol][num_tracer];
    Real trc_int_a[shcol][num_tracer];

    // Compute needed grid information from zi_grid
    // Grid stuff to compute based on zi_grid
    Real zt_grid[nlev];
    Real dz_zt[nlev];
    Real dz_zi[nlevi];
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

    // Load up tracer input array with random data
    //  ranging from values of 1 to 1000 (unitless)
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        for (Int t = 0; t < num_tracer; ++t){
          tracer_in[s][n][t] = (rand()%100 + 1)/1000.;
        }
      }
    }

    // Initialize data structure for bridging to F90
    UpdatePrognosticsImplicitData SDS(shcol,nlev,nlevi,num_tracer,dtime);
    // Note to validate this test we need to call the linear interp function
    //  to get the value of rho on zi grid at surface
    LinearInterpData SDSL(shcol,nlev,nlevi,0);

    // Test that the inputs are reasonable
    REQUIRE(SDS.shcol == shcol);
    REQUIRE(SDS.nlev == nlev);
    REQUIRE(SDS.nlevi == nlevi);
    REQUIRE(SDS.num_tracer == num_tracer);
    REQUIRE(shcol > 1);
    REQUIRE(nlev > 1);
    REQUIRE(nlevi == nlev+1);
    REQUIRE(num_tracer >= 1);

    // Fill in test data, first for column only input
    for(Int s = 0; s < shcol; ++s) {
      SDS.uw_sfc[s] = uw_sfc[s];
      SDS.vw_sfc[s] = vw_sfc[s];
      SDS.wthl_sfc[s] = wthl_sfc[s];
      SDS.wqw_sfc[s] = wqw_sfc[s];

      // Fill in tracer fluxes with random data from -10 to 10 (unitless)
      for (Int t = 0; t < num_tracer; ++t){
        const auto offset = t + s * num_tracer;
          SDS.wtracer_sfc[offset] = (rand()%20 + -10)/1000.;
      }

      // Fill in data on the nlev grid
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // TKH and TK get the same values on purpose
        SDS.tkh[offset] = tkh[n];
        SDS.tk[offset] = tkh[n];

        SDS.rho_zt[offset] = rho_zt[n];
        SDS.dz_zt[offset] = dz_zt[n];
        SDS.zt_grid[offset] = zt_grid[n];

        SDSL.y1[offset] = rho_zt[n];
        SDSL.x1[offset] = zt_grid[n];

        // Prognostic input/output variables
        SDS.thetal[offset] = thetal_in[n];
        SDS.qw[offset] = qw_in[n];
        SDS.u_wind[offset] = u_wind_in[n];
        SDS.v_wind[offset] = v_wind_in[n];
        SDS.tke[offset] = tke_in[n];

        for (Int t = 0; t < num_tracer; t++){
          const auto t_offset = t + offset * num_tracer;
          SDS.tracer[t_offset] = tracer_in[s][n][t];
        }
      }

      // Fill in data on the nlevi grid
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.zi_grid[offset] = zi_grid[n];
        SDS.dz_zi[offset] = dz_zi[n];

        SDSL.x2[offset] = zi_grid[n];
      }

    } // column loop

    // Check that inputs make sense
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(std::abs(SDS.wthl_sfc[s]) < 1);
      REQUIRE(std::abs(SDS.wqw_sfc[s]) < 1e-3);
      REQUIRE(std::abs(SDS.uw_sfc[s]) < 1);
      REQUIRE(std::abs(SDS.vw_sfc[s]) < 1);
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // Make sure inputs fall within reasonable bounds
        REQUIRE(SDS.dz_zt[offset] > 0);
        REQUIRE(SDS.zt_grid[offset] > 0);
        REQUIRE( (SDS.thetal[offset] > thl_lbound && SDS.thetal[offset] < thl_ubound) );
        REQUIRE( (SDS.qw[offset] > qw_lbound && SDS.qw[offset] < qw_ubound) );
        REQUIRE( (SDS.tke[offset] > tke_lbound && SDS.tke[offset] < tke_ubound) );
        REQUIRE( (SDS.rho_zt[offset] > rho_lbound && SDS.rho_zt[offset] < rho_ubound) );

        // While there is nothing unphysical with winds outside of these
        //  bounds, for this particular test we want to make sure the
        //  winds are modestly defined for checking later on.
        REQUIRE(std::abs(SDS.u_wind[offset]) < wind_bounds);
        REQUIRE(std::abs(SDS.v_wind[offset]) < wind_bounds);

        REQUIRE( (SDS.tkh[offset] > 0 && SDS.tkh[offset] > 0) );
        // For this case we assume diffusivity of heat equals the
        //  diffusivity of momentum.
        REQUIRE(SDS.tkh[offset] == SDS.tk[offset]);
      }

      // Make sure height arrays on zi grid are physical
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlevi;

        REQUIRE(SDS.zi_grid[offset] >= 0);
        REQUIRE(SDS.dz_zi[offset] >= 0);
      }

      // Check that zt increases in the upward direction
      for(Int n = 0; n < nlev - 1; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
      }

      // Check that zi increases in the upward direction
      for(Int n = 0; n < nlevi - 1; ++n) {
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0);
      }
    }

    // Compute integrals of input variables for energy check
    for (Int s = 0; s < shcol; ++s){
      qw_int_b[s] = 0; // Initialize
      thl_int_b[s] = 0; // Initialize
      for (Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        // Total water input integral
        qw_int_b[s] = qw_int_b[s] + SDS.qw[offset]*SDS.rho_zt[offset]
                                  *SDS.dz_zt[offset];
        // Potential temperature input integral
        thl_int_b[s] = thl_int_b[s] + SDS.thetal[offset]
                                  *SDS.rho_zt[offset]*SDS.dz_zt[offset];

        // Tracer input integrals
        for (Int t = 0; t < num_tracer; ++t){
          if (n == 0){
            trc_int_b[s][t] = 0;
          }
          const auto t_offset = t + offset * num_tracer;
          trc_int_b[s][t] = trc_int_b[s][t] + SDS.tracer[t_offset]*
                             SDS.rho_zt[offset]*SDS.dz_zt[offset];
        }
      }
    }

    // Call the C++ implementation
    update_prognostics_implicit(SDS);

    // Call linear interp to get rho value at surface for checking
    linear_interp(SDSL);

    // Check the result

    // First make sure that all output is within reasonable bounds and
    //  compute integrals of the outputs
    for(Int s = 0; s < shcol; ++s) {
      qw_int_a[s] = 0;
      thl_int_a[s] = 0;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // Make sure all output is reasonable
        REQUIRE( (SDS.thetal[offset] > thl_lbound && SDS.thetal[offset] < thl_ubound) );
        REQUIRE( (SDS.qw[offset] > qw_lbound && SDS.qw[offset] < qw_ubound) );
        REQUIRE( (SDS.tke[offset] > tke_lbound && SDS.tke[offset] < tke_ubound) );
        // Increase wind bounds by 2 m/s to allow for surface flux effects
        REQUIRE(std::abs(SDS.u_wind[offset] < wind_bounds+2));
        REQUIRE(std::abs(SDS.v_wind[offset] < wind_bounds+2));

        // Compute integrals of end result
        // Output total water integral
        qw_int_a[s] = qw_int_a[s] + SDS.rho_zt[offset]*SDS.dz_zt[offset]
                      *SDS.qw[offset];
        // Output potential temperature integral
        thl_int_a[s] = thl_int_a[s] + SDS.rho_zt[offset]*SDS.dz_zt[offset]
                      *SDS.thetal[offset];

        // Output tracer integral
        for (Int t = 0; t < num_tracer; ++t){
          if (n == 0){
            trc_int_a[s][t] = 0;
          }
          const auto t_offset = t + offset * num_tracer;
          trc_int_a[s][t] = trc_int_a[s][t] + SDS.tracer[t_offset]*
                             SDS.rho_zt[offset]*SDS.dz_zt[offset];
        }
      }

      // Get surface value of rho on zi grid
      Real rho_zi_srf, spurious;
      const auto offset_srf = (nlevi-1) + s * nlevi;
      rho_zi_srf = SDSL.y2[offset_srf];

      // Calculate the spurious source for total water
      spurious = (qw_int_a[s] - qw_int_b[s])/dtime
                   - rho_zi_srf*SDS.wqw_sfc[s];

      // Spurious source should be sufficiently small for water conservation
      REQUIRE(std::abs(spurious) < thresh_check);

      // Calculate the spurious source for thetal
      spurious = (thl_int_a[s] - thl_int_b[s])/dtime
                   - rho_zi_srf*SDS.wthl_sfc[s];

      // Spurious source should be sufficiently small for energy conservation
      REQUIRE(std::abs(spurious) < thresh_check);

      // Check that tracers were conserved during vertical transport
      for (Int t = 0; t < num_tracer; ++t){
        const auto t_offset = t + s * num_tracer;
        // Calculate spurious source
        spurious = (trc_int_a[s][t] - trc_int_b[s][t])/dtime
                   - rho_zi_srf*SDS.wtracer_sfc[t_offset];
        // Spurious source should be sufficiently small
        REQUIRE(std::abs(spurious) < thresh_check);
      }

    }

  } // run_property

  void run_bfb()
  {
    auto engine = Base::get_engine();

    UpdatePrognosticsImplicitData baseline_data[] = {
      UpdatePrognosticsImplicitData(10, 71, 72, 19, .5),
      UpdatePrognosticsImplicitData(10, 12, 13, 7, .25),
      UpdatePrognosticsImplicitData(7, 16, 17, 2, .1),
      UpdatePrognosticsImplicitData(2, 7, 8, 1, .1)
    };

    // Generate random input data
    for (auto& d : baseline_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    UpdatePrognosticsImplicitData cxx_data[] = {
      UpdatePrognosticsImplicitData(baseline_data[0]),
      UpdatePrognosticsImplicitData(baseline_data[1]),
      UpdatePrognosticsImplicitData(baseline_data[2]),
      UpdatePrognosticsImplicitData(baseline_data[3]),
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
      update_prognostics_implicit(d);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      static constexpr Int num_runs = sizeof(baseline_data) / sizeof(UpdatePrognosticsImplicitData);
      for (Int i = 0; i < num_runs; ++i) {
        UpdatePrognosticsImplicitData& d_baseline = baseline_data[i];
        UpdatePrognosticsImplicitData& d_cxx = cxx_data[i];

        REQUIRE(d_baseline.total(d_baseline.thetal) == d_cxx.total(d_cxx.thetal));
        REQUIRE(d_baseline.total(d_baseline.qw) == d_cxx.total(d_cxx.qw));
        REQUIRE(d_baseline.total(d_baseline.tke) == d_cxx.total(d_cxx.tke));
        REQUIRE(d_baseline.total(d_baseline.u_wind) == d_cxx.total(d_cxx.u_wind));
        REQUIRE(d_baseline.total(d_baseline.v_wind) == d_cxx.total(d_cxx.v_wind));
        for (Int k = 0; k < d_baseline.total(d_baseline.thetal); ++k) {
          REQUIRE(d_baseline.thetal[k] == d_cxx.thetal[k]);
          REQUIRE(d_baseline.qw[k] == d_cxx.qw[k]);
          REQUIRE(d_baseline.tke[k] == d_cxx.tke[k]);
          REQUIRE(d_baseline.u_wind[k] == d_cxx.u_wind[k]);
          REQUIRE(d_baseline.v_wind[k] == d_cxx.v_wind[k]);
        }

        REQUIRE(d_baseline.total(d_baseline.tracer) == d_cxx.total(d_cxx.tracer));
        for (Int k = 0; k < d_baseline.total(d_baseline.tracer); ++k) {
          REQUIRE(d_baseline.tracer[k] == d_cxx.tracer[k]);
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

TEST_CASE("update_prognostics_implicit_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestUpdatePrognosticsImplicit;

  TestStruct().run_property();
}

TEST_CASE("update_prognostics_implicit_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestUpdatePrognosticsImplicit;

  TestStruct().run_bfb();
}

} // empty namespace
