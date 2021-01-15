#include <iostream>
#include <cmath>
#include "catch2/catch.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "rrtmgp_test_utils.hpp"
#include "netcdf.h"
#include "mo_gas_concentrations.h"
#include "mo_fluxes.h"
#include "mo_cloud_optics.h"
#include "mo_garand_atmos_io.h"
#include "Intrinsics.h"
namespace {

    template <class T> double arrmin(T &arr) {
        double minval = arr.myData[0];
        for (int i = 0; i<arr.totElems(); i++) {
            if (arr.myData[i] < minval) {
                minval = arr.myData[i];
            }
        }
        return minval;
    }

    TEST_CASE("rrtmgp_init", "") {
        // Try to initialize RRTMGP with known gas names
        const int ngas = 8;
        const std::string gas_names[8] = {
            "h2o", "co2", "o3", "n2o",
            "co" , "ch4", "o2", "n2"
        };
        scream::rrtmgp::rrtmgp_initialize(ngas, gas_names);
        REQUIRE(1 == 1);

        // Sanity check to see if we were able to load the data correctly. 
        // Check integer part of reference pressure to avoid floating point 
        // differences.
        REQUIRE(int(scream::rrtmgp::k_dist_sw.press_ref(1)) == 109663);
        REQUIRE(int(scream::rrtmgp::k_dist_sw.press_ref(size(scream::rrtmgp::k_dist_sw.press_ref,1))) == 1);
    }

    TEST_CASE("rrtmgp_run", "") {
        // Get reference fluxes from input file; do this here so we can get ncol dimension
        std::string inputfile = "./data/rrtmgp-allsky.nc";
        real2d sw_flux_up_ref;
        real2d sw_flux_dn_ref;
        real2d sw_flux_dn_dir_ref;
        real2d lw_flux_up_ref;
        real2d lw_flux_dn_ref;
        rrtmgpTest::read_fluxes(inputfile, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

        // Get dimension sizes
        int ncol = sw_flux_up_ref.dimension[0];
        int nlev = sw_flux_up_ref.dimension[1];
        int nlay = nlev - 1;

        // Read in dummy Garand atmosphere; if this were an actual model simulation, 
        // these would be passed as inputs to the driver
        // NOTE: set ncol to size of col_flx dimension in the input file. This is so
        // that we can compare to the reference data provided in that file. Note that
        // this will copy the first column of the input data (the first profile) ncol
        // times. We will then fill some fraction of these columns with clouds for
        // the test problem.
        real2d p_lay;
        real2d t_lay;
        real2d p_lev;
        real2d t_lev;
        real2d col_dry;
        GasConcs gas_concs;
        real2d sfc_alb_dir;
        real2d sfc_alb_dif;
        real1d mu0;
        real2d lwp;
        real2d iwp;
        real2d rel;
        real2d rei;
        read_atmos(inputfile, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, ncol);

        // Check that data has the shape we expect
        REQUIRE(p_lay.dimension[0] == ncol);
        REQUIRE(p_lay.dimension[1] == nlay);
        REQUIRE(p_lev.dimension[1] == nlev);

        // Initialize absorption coefficients
        int ngas = gas_concs.get_num_gases();
        string1d gas_names_1d = gas_concs.get_gas_names();
        std::string gas_names[ngas];
        for (int igas = 0; igas < ngas; igas++) {
            gas_names[igas] = gas_names_1d(igas+1);
        }
        scream::rrtmgp::rrtmgp_initialize(ngas, gas_names);

        // Check that data was loaded properly; just sanity check a few values
        REQUIRE(int(p_lay(1,1)) == 100933);
        REQUIRE(int(p_lay(1,p_lay.dimension[1])) == 19);

        // Setup our dummy atmosphere based on the input data we read in
        rrtmgpTest::dummy_atmos(
                inputfile, ncol, p_lay, t_lay,
                sfc_alb_dir, sfc_alb_dif, mu0,
                lwp, iwp, rel, rei
            );

        // Setup flux outputs; In a real model run, the fluxes would be
        // input/outputs into the driver (persisting between calls), and
        // we would just have to setup the pointers to them in the
        // FluxesBroadband object
        real2d sw_flux_up ("sw_flux_up" ,ncol,nlay+1);
        real2d sw_flux_dn ("sw_flux_dn" ,ncol,nlay+1);
        real2d sw_flux_dn_dir("sw_flux_dn_dir",ncol,nlay+1);
        real2d lw_flux_up ("lw_flux_up" ,ncol,nlay+1);
        real2d lw_flux_dn ("lw_flux_dn" ,ncol,nlay+1);

        // Run RRTMGP code on dummy atmosphere; this might get ugly
        // Inputs should be atmosphere state, outputs should be fluxes
        // TODO: should absorption coefficients be an input, or should that be initialized
        // and kept in the scream::rrtmgp namespace?
        scream::rrtmgp::rrtmgp_main(
                p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, 
                sfc_alb_dir, sfc_alb_dif, mu0,
                lwp, iwp, rel, rei,
                sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
                lw_flux_up, lw_flux_dn);
 
        // Check values
        REQUIRE(rrtmgpTest::all_equals(sw_flux_up_ref    , sw_flux_up    ));
        REQUIRE(rrtmgpTest::all_equals(sw_flux_dn_ref    , sw_flux_dn    ));
        REQUIRE(rrtmgpTest::all_equals(sw_flux_dn_dir_ref, sw_flux_dn_dir));
        REQUIRE(rrtmgpTest::all_equals(lw_flux_up_ref    , lw_flux_up    ));
        REQUIRE(rrtmgpTest::all_equals(lw_flux_dn_ref    , lw_flux_dn    ));

        // ALTERNATIVELY: create a single or two-layer atmosphere to do a dummy calc

        // Clean up
        scream::rrtmgp::rrtmgp_finalize();
    }
   
} // empty namespace
