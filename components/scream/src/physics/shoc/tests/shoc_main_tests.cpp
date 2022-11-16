#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestShocMain {

  static void run_property()
  {
    static constexpr Real mintke = scream::shoc::Constants<Real>::mintke;
    static constexpr Real minlen = scream::shoc::Constants<Real>::minlen;
    static constexpr Real maxlen = scream::shoc::Constants<Real>::maxlen;
    static constexpr Real maxiso = scream::shoc::Constants<Real>::maxiso;
    static constexpr Real Cpair = scream::physics::Constants<Real>::Cpair;
    static constexpr Real gravit = scream::physics::Constants<Real>::gravit;
    static constexpr Real LatVap = scream::physics::Constants<Real>::LatVap;
    static constexpr Real Rair = scream::physics::Constants<Real>::Rair;
    static constexpr Real p0 = scream::physics::Constants<Real>::P0;

    static constexpr Int shcol    = 5;
    static constexpr Int nlev     = 5;
    static constexpr auto nlevi   = nlev + 1;
    static constexpr Int num_qtracers = 3;
    static constexpr Int nadv = 5;

    // Tests for the subroutine shoc_main
    //  This is the top level function of SHOC.  Note that all
    //  SHOC functions have property tests to throughly vet those
    //  functions.  Thus the purpose here is not to duplicate those
    //  tests but merly verify that SHOC returns values as expected
    //  to verify that these individual functions are stitched
    //  together properly.  Verifying that all outputs fall within the
    //  physically reasonable bounds given various surface forcing is
    //  sufficient to validate this.

    // Timestep [s]
    static constexpr Real dtime = 300;
    // host dx [m]
    static constexpr Real host_dx = 3000;
    // host dy [m]
    static constexpr Real host_dy = 3000;

    // Define PROFILE variables
    // Define the heights on the zi grid [m]
    static constexpr Real zi_grid[nlevi] = {3000, 2000, 1500, 1000, 500, 0};
    // Define pressures on the interface grid [Pa]
    static constexpr Real presi[nlevi] = {700e2, 800e2, 850e2, 900e2, 950e2, 1000e2};
    // Define temperature on the zt grid [K]
    static constexpr Real temp[nlev] = {290, 295, 297, 297, 300};
    // Define the large scale vertical velocity on zt grid [m/s]
    static constexpr Real w_field[nlev] = {5e-2, 4e-2, 3e-2, 2e-2, 1e-2};
    // Define the zonal wind [m/s]
    static constexpr Real u_wind[nlev] = {-4.61, -6, -8.75, -7.75, -7.75};
    // define the meridional wind [m/s]
    static constexpr Real v_wind[nlev] = {0, 0, 0, 0, 0};
    // Define the total water mixing ratio [kg/kg]
    static constexpr Real qw[nlev] = {1e-3, 4.2e-3, 10.7e-3, 16.3e-3, 17e-3};
    // Define the TKE [m2/s2]
    static constexpr Real tke[nlev] = {mintke, 0.1, 0.3, 0.2, 0.1};
    // Define the eddy vicosity for heat and momentum [m2/s]
    static constexpr Real tkh[nlev] = {3, 10, 50, 30, 20};
    // Buoyancy flux [K m/s]
    static constexpr Real wthv_sec[nlev] = {-0.02, 0.04, 0.03, -0.02, 0.03};
    // SHOC cloud liquid water [kg/kg]
    static constexpr Real shoc_ql[nlev] = {0, 0, 1e-3, 1e-4, 0};
    // SHOC cloud fraction [-]
    static constexpr Real shoc_cldfrac[nlev] = {0, 0, 0.8, 0.2, 0};

    // heat flux at surface [K m/s], COLUMN ONLY variables
    static constexpr Real wthl_sfc[shcol] = {0.03, -0.03, 0.1, 0, -0.1};
    // moisture flux at surface [kg/kg m/s]
    static constexpr Real wqw_sfc[shcol] = {2e-5, 1e-6, 0, -2e-5, 1e-4};
    // Surface moment flux, zonal direction [m3/s3]
    static constexpr Real uw_sfc[shcol] = {0.03, -0.03, 0.1, 0, -0.1};
    // Surface moment flux, meridional direction [m3/s3]
    static constexpr Real vw_sfc[shcol] = {-0.01, -0.01, 0.3, 0, -0.3};
    // Surfafce geopotential height
    static constexpr Real phis[shcol] = {100, 200, 0, 150, 500};

    // establish reasonable bounds for checking input/output
    static constexpr Real thl_lbound = 200; // [K]
    static constexpr Real thl_ubound = 350; // [K]
    static constexpr Real qw_lbound = 1e-4; // [kg/kg]
    static constexpr Real qw_ubound = 5e-2; // [kg/kg]
    static constexpr Real tke_lbound = 0; // [m2/s2]
    static constexpr Real tke_ubound = 5; // [m2/s2]
    static constexpr Real wind_bounds = 10; // [m/s]
    static constexpr Real dse_upper = 5e5; // [J/kg]
    static constexpr Real dse_lower = 1e5; // [J/kg]

    // Compute some inputs based on the above

    // Input for tracer (no units)
    Real tracer_in[shcol][nlev][num_qtracers];

    // First compute variables related to height
    Real pres[nlev], zt_grid[nlev], pdel[nlev];
    for(Int n = 0; n < nlev; ++n) {
      // height on the midpoint grid
      zt_grid[n] = 0.5*(zi_grid[n]+zi_grid[n+1]);
      // pressure on the midpoint grid
      pres[n] = 0.5*(presi[n]+presi[n+1]);
      // pressure thickness
      pdel[n] = presi[n+1] - presi[n];
    }

    // Compute variables related to temperature
    Real host_dse[shcol][nlev];
    Real thetal[nlev], thv[nlev], inv_exner[nlev];
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        // Compute the dry static energy
        host_dse[s][n] = Cpair*temp[n]+gravit*zt_grid[n]+phis[s];
      }
    }

    Real pot_temp, qv;
    for(Int n = 0; n < nlev; ++n) {
      // Compute potential temperature and water vapor
      pot_temp = temp[n]*std::pow(p0/pres[n],Rair/Cpair);
      qv = qw[n] - shoc_ql[n];
      // Liquid water potential temperature
      thetal[n] = pot_temp - (LatVap/Cpair)*shoc_ql[n];
      // Virtual potential temperature
      thv[n] = pot_temp * (1 + 0.61*qv - shoc_ql[n]);
      // Compute the inverse of the exner function
      const Real exner = std::pow(pres[n]/p0,Rair/Cpair);
      REQUIRE(exner > 0);
      inv_exner[n] = 1/exner;
    }

    // Load up tracer input array with random data
    //  ranging from values of 0 to 1e-1 (unitless)
    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev; ++n) {
        for (Int t = 0; t < num_qtracers; ++t){
          tracer_in[s][n][t] = (rand()%100 + 1)/1000.;
        }
      }
    }
    // Initialize data structure for bridging to F90
    ShocMainData SDS(shcol,nlev,nlevi,num_qtracers,dtime,nadv,0,0);

    // Test the inputs are good
    REQUIRE(SDS.shcol == shcol);
    REQUIRE(SDS.nlev == nlev);
    REQUIRE(SDS.nlevi == nlevi);
    REQUIRE(SDS.num_qtracers == num_qtracers);
    REQUIRE(SDS.dtime == dtime);
    REQUIRE(SDS.nadv == nadv);
    REQUIRE(shcol > 0);
    REQUIRE(nlev > 1);
    REQUIRE(nlevi == nlev+1);
    REQUIRE(num_qtracers >= 1);
    REQUIRE(dtime > 0);
    REQUIRE(nadv > 0);

    // Fill in test data, first for column only input
    for(Int s = 0; s < shcol; ++s) {
      SDS.uw_sfc[s] = uw_sfc[s];
      SDS.vw_sfc[s] = vw_sfc[s];
      SDS.wthl_sfc[s] = wthl_sfc[s];
      SDS.wqw_sfc[s] = wqw_sfc[s];
      SDS.phis[s] = phis[s];
      SDS.host_dx[s] = host_dx;
      SDS.host_dy[s] = host_dy;

      // Fill in tracer fluxes with random data from -1e-2 to 1e-2 (unitless)
      for (Int t = 0; t < num_qtracers; ++t){
        const auto offset = t + s * num_qtracers;
          SDS.wtracer_sfc[offset] = (rand()%20 + -10)/1000.;
      }

      // Fill in data on the nlev grid
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        SDS.zt_grid[offset] = zt_grid[n];
        SDS.pres[offset] = pres[n];
        SDS.pdel[offset] = pdel[n];
        SDS.thv[offset] = thv[n];
        SDS.w_field[offset] = w_field[n];
        SDS.inv_exner[offset] = inv_exner[n];
        SDS.shoc_ql[offset] = shoc_ql[n];
        SDS.shoc_cldfrac[offset] = shoc_cldfrac[n];
        SDS.qw[offset] = qw[n];
        SDS.thetal[offset] = thetal[n];
        SDS.u_wind[offset] = u_wind[n];
        SDS.v_wind[offset] = v_wind[n];
        SDS.tke[offset] = tke[n];
        SDS.wthv_sec[offset] = wthv_sec[n];
        SDS.host_dse[offset] = host_dse[s][n];

        // TKH and TK get the same values on purpose
        SDS.tkh[offset] = tkh[n];
        SDS.tk[offset] = tkh[n];

        for (Int t = 0; t < num_qtracers; t++){
          const auto t_offset = t + offset * num_qtracers;
          SDS.qtracers[t_offset] = tracer_in[s][n][t];
        }

      }

      // Fill in data on the nlevi grid
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        SDS.zi_grid[offset] = zi_grid[n];
        SDS.presi[offset] = presi[n];

      }
    }

    // Check that inputs make sense

    for(Int s = 0; s < shcol; ++s) {
      for(Int n = 0; n < nlev - 1; ++n) {
        const auto offset = n + s * nlev;
        // Check that zt increases upward
        REQUIRE(SDS.zt_grid[offset + 1] - SDS.zt_grid[offset] < 0);
        // Check that inverse of exner increases upward
        REQUIRE(SDS.inv_exner[offset + 1] - SDS.inv_exner[offset] < 0);
      }

      // Check that zi increases upward
      for(Int n = 0; n < nlevi - 1; ++n) {
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.zi_grid[offset + 1] - SDS.zi_grid[offset] < 0);
      }

      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;

        // Make sure inputs fall within reasonable bounds
        REQUIRE(SDS.zt_grid[offset] > 0);
        REQUIRE( (SDS.thetal[offset] > thl_lbound && SDS.thetal[offset] < thl_ubound) );
        REQUIRE( (SDS.qw[offset] > qw_lbound && SDS.qw[offset] < qw_ubound) );
        REQUIRE( (SDS.tke[offset] > tke_lbound && SDS.tke[offset] < tke_ubound) );

        // While there is nothing unphysical with winds outside of these
        //  bounds, for this particular test we want to make sure the
        //  winds are modestly defined for checking later on.
        REQUIRE(std::abs(SDS.u_wind[offset]) < wind_bounds);
        REQUIRE(std::abs(SDS.v_wind[offset]) < wind_bounds);
      }

    }

    // Call the fortran implementation
    shoc_main(SDS);

    // Make sure output falls within reasonable bounds
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE( (SDS.pblh[s] > 0 && SDS.pblh[s] <= zi_grid[0]) );
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        REQUIRE( (SDS.thetal[offset] > thl_lbound && SDS.thetal[offset] < thl_ubound) );
        REQUIRE( (SDS.qw[offset] > qw_lbound && SDS.qw[offset] < qw_ubound) );
        REQUIRE( (SDS.tke[offset] > tke_lbound && SDS.tke[offset] < tke_ubound) );
        // Increase wind bounds by 2 m/s to allow for surface flux effects
        REQUIRE(std::abs(SDS.u_wind[offset] < wind_bounds+2));
        REQUIRE(std::abs(SDS.v_wind[offset] < wind_bounds+2));

        REQUIRE( (SDS.shoc_mix[offset] >= minlen && SDS.shoc_mix[offset] <= maxlen) );
        REQUIRE( (SDS.isotropy[offset] >= 0 && SDS.isotropy[offset] < maxiso) );
        REQUIRE( (SDS.shoc_cldfrac[offset] >= 0 && SDS.shoc_cldfrac[offset] <= 1) );
        REQUIRE( (SDS.shoc_ql[offset] >= 0 && SDS.shoc_ql[offset] <= 1) );
        REQUIRE( (SDS.tk[offset] >= 0 && SDS.tk[offset] < 100) );
        REQUIRE( (SDS.tkh[offset] >= 0 && SDS.tkh[offset] < 100) );
        REQUIRE(std::abs(SDS.wthv_sec[offset]) < 1);
        REQUIRE(std::abs(SDS.brunt[offset] < 1));

        // Make sure there are no "empty" clouds
        if (SDS.shoc_cldfrac[offset] > 0){
          REQUIRE(SDS.shoc_ql[offset] > 0);
        }
        else if (SDS.shoc_cldfrac[offset] == 0){
          REQUIRE(SDS.shoc_ql[offset] == 0);
        }

        // Verify that dse increases with height upward
        if (n < nlev-1){
          REQUIRE(SDS.host_dse[offset + 1] - SDS.host_dse[offset] < 0);
        }

        // Verify that DSE falls within some reasonable bounds
        REQUIRE(SDS.host_dse[offset] > dse_lower);
        REQUIRE(SDS.host_dse[offset] < dse_upper);

        // Verify that w2 is less than tke
        REQUIRE(std::abs(SDS.w_sec[offset] < SDS.tke[offset]));

        // Verify tracer output is reasonable
        for (Int t = 0; t < num_qtracers; ++t){
          const auto t_offset = t + offset * num_qtracers;
          REQUIRE(SDS.qtracers[t_offset] < 1);
        }

      }

      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;

        // Make sure turbulence output is appropriate
        REQUIRE( (SDS.thl_sec[offset] >= 0 && SDS.thl_sec[offset] < 1e2 ) );
        REQUIRE( (SDS.qw_sec[offset] >= 0 && SDS.qw_sec[offset] < 1e-3) );
        REQUIRE(std::abs(SDS.wthl_sec[offset] < 1));
        REQUIRE(std::abs(SDS.wqw_sec[offset] < 1e-3));
        REQUIRE(std::abs(SDS.w3[offset] < 10));
        REQUIRE(std::abs(SDS.qwthl_sec[offset] < 1));
        REQUIRE(std::abs(SDS.uw_sec[offset] < 1));
        REQUIRE(std::abs(SDS.vw_sec[offset] < 1));
      }

    }

  } // run_property

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ShocMainData f90_data[] = {
      //           shcol, nlev, nlevi, num_qtracers, dtime, nadv, nbot_shoc, ntop_shoc(C++ indexing)
      ShocMainData(12,      72,    73,            5,   300,   15,        72, 0),
      ShocMainData(8,       12,    13,            3,   300,   10,         8, 3),
      ShocMainData(7,       16,    17,            3,   300,    1,        12, 0),
      ShocMainData(2,       7,      8,            2,   300,    5,         7, 4)
    };

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize(engine,
                  {
                    {d.presi, {700e2,1000e2}},
                    {d.tkh, {3,50}},
                    {d.tke, {0.1,0.3}},
                    {d.zi_grid, {0, 3000}},
                    {d.wthl_sfc, {0,1e-4}},
                    {d.wqw_sfc, {0,1e-6}},
                    {d.uw_sfc, {0,1e-2}},
                    {d.vw_sfc, {0,1e-4}},
                    {d.host_dx, {3000, 3000}},
                    {d.host_dy, {3000, 3000}},
                    {d.phis, {0, 500}}, // 500
                    {d.wthv_sec, {-0.02, 0.03}},
                    {d.qw, {1e-4, 5e-2}},
                    {d.u_wind, {-10, 0}},
                    {d.v_wind, {-10, 0}},
                    {d.shoc_ql, {0, 1e-3}},
                  });
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ShocMainData cxx_data[] = {
      ShocMainData(f90_data[0]),
      ShocMainData(f90_data[1]),
      ShocMainData(f90_data[2]),
      ShocMainData(f90_data[3])
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      shoc_main_with_init(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      const int npbl = shoc_init_f(d.nlev, d.pref_mid, d.nbot_shoc, d.ntop_shoc);

      shoc_main_f(d.shcol, d.nlev, d.nlevi, d.dtime, d.nadv, npbl, d.host_dx, d.host_dy,
                  d.thv, d.zt_grid, d.zi_grid, d.pres, d.presi, d.pdel, d.wthl_sfc,
                  d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.wtracer_sfc, d.num_qtracers,
                  d.w_field, d.inv_exner, d.phis, d.host_dse, d.tke, d.thetal, d.qw,
                  d.u_wind, d.v_wind, d.qtracers, d.wthv_sec, d.tkh, d.tk, d.shoc_ql,
                  d.shoc_cldfrac, d.pblh, d.shoc_mix, d.isotropy, d.w_sec, d.thl_sec,
                  d.qw_sec, d.qwthl_sec, d.wthl_sec, d.wqw_sec, d.wtke_sec, d.uw_sec,
                  d.vw_sec, d.w3, d.wqls_sec, d.brunt, d.shoc_ql2);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(ShocMainData);

      for (Int i = 0; i < num_runs; ++i) {
        ShocMainData& d_f90 = f90_data[i];
        ShocMainData& d_cxx = cxx_data[i];
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.host_dse));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.tke));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.thetal));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.qw));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.u_wind));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.v_wind));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.wthv_sec));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.tkh));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.tk));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.shoc_ql));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.shoc_cldfrac));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.shoc_mix));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.isotropy));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.w_sec));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.wqls_sec));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.brunt));
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.shoc_ql2));
        for (Int k = 0; k < d_f90.total(d_f90.host_dse); ++k) {
          REQUIRE(d_f90.host_dse[k] == d_cxx.host_dse[k]);
          REQUIRE(d_f90.tke[k] == d_cxx.tke[k]);
          REQUIRE(d_f90.thetal[k] == d_cxx.thetal[k]);
          REQUIRE(d_f90.qw[k] == d_cxx.qw[k]);
          REQUIRE(d_f90.u_wind[k] == d_cxx.u_wind[k]);
          REQUIRE(d_f90.v_wind[k] == d_cxx.v_wind[k]);
          REQUIRE(d_f90.wthv_sec[k] == d_cxx.wthv_sec[k]);
          REQUIRE(d_f90.tk[k] == d_cxx.tk[k]);
          REQUIRE(d_f90.shoc_ql[k] == d_cxx.shoc_ql[k]);
          REQUIRE(d_f90.shoc_cldfrac[k] == d_cxx.shoc_cldfrac[k]);
          REQUIRE(d_f90.shoc_mix[k] == d_cxx.shoc_mix[k]);
          REQUIRE(d_f90.isotropy[k] == d_cxx.isotropy[k]);
          REQUIRE(d_f90.w_sec[k] == d_cxx.w_sec[k]);
          REQUIRE(d_f90.wqls_sec[k] == d_cxx.wqls_sec[k]);
          REQUIRE(d_f90.brunt[k] == d_cxx.brunt[k]);
          REQUIRE(d_f90.shoc_ql2[k] == d_cxx.shoc_ql2[k]);
        }

        REQUIRE(d_f90.total(d_f90.qtracers) == d_cxx.total(d_cxx.qtracers));
        for (Int k = 0; k < d_f90.total(d_f90.qtracers); ++k) {
          REQUIRE(d_f90.qtracers[k] == d_cxx.qtracers[k]);
        }

        REQUIRE(d_f90.total(d_f90.pblh) == d_cxx.total(d_cxx.pblh));
        for (Int k = 0; k < d_f90.total(d_f90.pblh); ++k) {
          REQUIRE(d_f90.pblh[k] == d_cxx.pblh[k]);
        }

        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.thl_sec));
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.qw_sec));
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.qwthl_sec));
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.wthl_sec));
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.wqw_sec));
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.wtke_sec));
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.uw_sec));
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.vw_sec));
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.w3));
        for (Int k = 0; k < d_f90.total(d_f90.thl_sec); ++k) {
          REQUIRE(d_f90.thl_sec[k] == d_cxx.thl_sec[k]);
          REQUIRE(d_f90.qw_sec[k] == d_cxx.qw_sec[k]);
          REQUIRE(d_f90.qwthl_sec[k] == d_cxx.qwthl_sec[k]);
          REQUIRE(d_f90.wthl_sec[k] == d_cxx.wthl_sec[k]);
          REQUIRE(d_f90.wqw_sec[k] == d_cxx.wqw_sec[k]);
          REQUIRE(d_f90.wtke_sec[k] == d_cxx.wtke_sec[k]);
          REQUIRE(d_f90.uw_sec[k] == d_cxx.uw_sec[k]);
          REQUIRE(d_f90.vw_sec[k] == d_cxx.vw_sec[k]);
          REQUIRE(d_f90.w3[k] == d_cxx.w3[k]);
        }
      }
    }
  } // run_bfb
};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("shoc_main_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocMain;

  TestStruct::run_property();
}

TEST_CASE("shoc_main_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocMain;

  TestStruct::run_bfb();
}

} // empty namespace
