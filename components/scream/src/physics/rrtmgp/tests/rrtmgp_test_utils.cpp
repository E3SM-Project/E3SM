#include <iostream>
#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "rrtmgp_test_utils.hpp"
#include "mo_gas_concentrations.h"
#include "mo_garand_atmos_io.h"
#include "mo_fluxes.h"
#include "mo_cloud_optics.h"
namespace rrtmgpTest {

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

    void dummy_atmos(
            std::string inputfile, 
            int ncol, real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev, GasConcs &gas_concs, real2d &col_dry, 
            real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0,
            real2d &lwp, real2d &iwp, real2d &rel, real2d &rei) {

        // Read example atmosphere profiles
        read_atmos(inputfile, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, ncol);

        // Setup boundary conditions, solar zenith angle, etc
        // NOTE: this stuff would come from the model in a real run
        int nbndsw = scream::rrtmgp::k_dist_sw.get_nband();
        sfc_alb_dir = real2d("sfc_alb_dir", nbndsw, ncol);
        sfc_alb_dif = real2d("sfc_alb_dif", nbndsw, ncol);

        // Ocean-ish values for surface albedos, just for example
        memset(sfc_alb_dir , 0.06_wp );
        memset(sfc_alb_dif , 0.06_wp );

        // Pick a solar zenith angle; this should come from the model
        mu0 = real1d("mu0", ncol);
        memset(mu0, 0.86_wp );

        // Get dummy cloud PHYSICAL properties. Note that this function call
        // needs the CloudOptics object only because it uses the min and max
        // valid values from the lookup tables for liquid and ice water path to
        // create a dummy atmosphere.
        dummy_clouds(scream::rrtmgp::cloud_optics_sw, p_lay, t_lay, lwp, iwp, rel, rei);
    }

    void dummy_clouds(
            CloudOptics &cloud_optics, real2d &p_lay, real2d &t_lay, 
            real2d &lwp, real2d &iwp, real2d &rel, real2d &rei) {

        // Problem sizes
        int ncol = t_lay.dimension[0];
        int nlay = t_lay.dimension[1];

        // Generate some fake liquid and ice water data. We pick values to be midway between
        // the min and max of the valid lookup table values for effective radii
        real rel_val = 0.5 * (cloud_optics.get_min_radius_liq() + cloud_optics.get_max_radius_liq());
        real rei_val = 0.5 * (cloud_optics.get_min_radius_ice() + cloud_optics.get_max_radius_ice());

        // Restrict clouds to troposphere (> 100 hPa = 100*100 Pa) and not very close to the ground (< 900 hPa), and
        // put them in 2/3 of the columns since that's roughly the total cloudiness of earth.
        // Set sane values for liquid and ice water path.
        rel = real2d("rel", ncol, nlay);
        rei = real2d("rei", ncol, nlay);
        lwp = real2d("lwp", ncol, nlay);
        iwp = real2d("iwp", ncol, nlay);
        real2d cloud_mask("cloud_mask", ncol, nlay);
        parallel_for( Bounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
            cloud_mask(icol,ilay) = p_lay(icol,ilay) > 100._wp * 100._wp && p_lay(icol,ilay) < 900._wp * 100._wp && mod(icol, 3) != 0;
            // Ice and liquid will overlap in a few layers
            lwp(icol,ilay) = merge(10._wp,  0._wp, cloud_mask(icol,ilay) && t_lay(icol,ilay) > 263._wp);
            iwp(icol,ilay) = merge(10._wp,  0._wp, cloud_mask(icol,ilay) && t_lay(icol,ilay) < 273._wp);
            rel(icol,ilay) = merge(rel_val, 0._wp, lwp(icol,ilay) > 0._wp);
            rei(icol,ilay) = merge(rei_val, 0._wp, iwp(icol,ilay) > 0._wp);
        });
    }

    // Function to read fluxes from input file so we can compare our answers against the reference
    void read_fluxes(
            std::string inputfile, 
            real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dn_dir,
            real2d &lw_flux_up, real2d &lw_flux_dn) {

        // Initialize netcdf reader
        yakl::SimpleNetCDF io;
        io.open(inputfile, yakl::NETCDF_MODE_READ);

        // Initialize arrays to hold fluxes
        int nlev = io.getDimSize("lev");
        int ncol = io.getDimSize("col_flx");
        sw_flux_up = real2d("sw_flux_up", ncol, nlev);
        sw_flux_dn = real2d("sw_flux_dn", ncol, nlev);
        sw_flux_dn_dir = real2d("sw_flux_dn_dir", ncol, nlev);
        lw_flux_up = real2d("lw_flux_up", ncol, nlev);
        lw_flux_dn = real2d("lw_flux_dn", ncol, nlev);

        // Read data
        io.read(sw_flux_up, "sw_flux_up_result");
        io.read(sw_flux_dn, "sw_flux_dn_result");
        io.read(sw_flux_dn_dir, "sw_flux_dir_result");
        io.read(lw_flux_up, "lw_flux_up_result");
        io.read(lw_flux_dn, "lw_flux_dn_result");
    }

}  // namespace rrtmgp
