#include <iostream>
#include <cmath>
#include "catch2/catch.hpp"
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "rrtmgp_dummy_atmos.hpp"
#include "netcdf.h"
#include "mo_gas_concentrations.h"
#include "mo_fluxes.h"
#include "mo_cloud_optics.h"
#include "FortranIntrinsics.h"
namespace {

    // TODO: use YAKL intrinsics for this; this won't work on the GPU
    bool all_equals(real2d &arr1, real2d &arr2) {
        double tolerance = 0.01;
        /*
        real2d residual = arr1 - arr2;
        if (yakl::fortran::anyGT(residual, tolerance) || yakl::fortran::anyLT(residual, -tolerance)) {
            printf("max(arr1 - arr2) = %f\n", yakl::fortran::maxval(residual));
            return false;
        } else {
            return true;
        }
        */
        int nx = arr1.dimension[0];
        int ny = arr2.dimension[1];
        for (int i=1; i<nx+1; i++) {
            for (int j=1; j<ny+1; j++) {
                if (abs(arr1(i,j) - arr2(i,j)) > tolerance) {
                    printf("arr1 = %f, arr2 = %f\n", arr1(i,j), arr2(i,j));
                    return false;
                }
            }
        }
        return true;
    }

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
        scream::rrtmgp::rrtmgp_initialize();
        REQUIRE(1 == 1);

        // Sanity check to see if we were able to load the data correctly. 
        // Check integer part of reference pressure to avoid floating point 
        // differences.
        REQUIRE(int(scream::rrtmgp::k_dist_sw.press_ref(1)) == 109663);
        REQUIRE(int(scream::rrtmgp::k_dist_sw.press_ref(size(scream::rrtmgp::k_dist_sw.press_ref,1))) == 1);
    }

    TEST_CASE("rrtmgp_run", "") {
        // Initialize absorption coefficients
        scream::rrtmgp::rrtmgp_initialize();

        // Get reference fluxes from input file; do this here so we can get ncol dimension
        std::string inputfile = "./data/rrtmgp-allsky.nc";
        real2d sw_flux_up_ref;
        real2d sw_flux_dn_ref;
        real2d sw_flux_dn_dir_ref;
        real2d lw_flux_up_ref;
        real2d lw_flux_dn_ref;
        rrtmgp::read_fluxes(inputfile, sw_flux_up_ref, sw_flux_dn_ref, sw_flux_dn_dir_ref, lw_flux_up_ref, lw_flux_dn_ref );

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
        rrtmgp::dummy_atmos(
                inputfile, ncol, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, 
                sfc_alb_dir, sfc_alb_dif, mu0,
                lwp, iwp, rel, rei
            );

        // Check that data has the shape we expect
        REQUIRE(p_lay.dimension[0] == ncol);
        REQUIRE(p_lay.dimension[1] == nlay);
        REQUIRE(p_lev.dimension[1] == nlev);

        // Check that data was loaded properly; just sanity check a few values
        REQUIRE(int(p_lay(1,1)) == 100933);
        REQUIRE(int(p_lay(1,p_lay.dimension[1])) == 19);

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
        REQUIRE(all_equals(sw_flux_up_ref    , sw_flux_up    ));
        REQUIRE(all_equals(sw_flux_dn_ref    , sw_flux_dn    ));
        REQUIRE(all_equals(sw_flux_dn_dir_ref, sw_flux_dn_dir));
        REQUIRE(all_equals(lw_flux_up_ref    , lw_flux_up    ));
        REQUIRE(all_equals(lw_flux_dn_ref    , lw_flux_dn    ));

        // ALTERNATIVELY: create a single or two-layer atmosphere to do a dummy calc
    }
   
} // empty namespace
