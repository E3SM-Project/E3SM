#include "physics/rrtmgp/scream_rrtmgp_interface.hpp"
#include "physics/rrtmgp/rrtmgp_test_utils.hpp"
#include "physics/rrtmgp/mo_garand_atmos_io.h"
#include "physics/rrtmgp/simple_netcdf.hpp"

#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "cpp/rte/mo_fluxes.h"
#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"

#include <iostream>
#include <netcdf.h>

namespace rrtmgpTest {

    using yakl::fortran::parallel_for;
    using yakl::fortran::SimpleBounds;
    using yakl::intrinsics::mod;
    using yakl::intrinsics::merge;

    bool file_exists(const char *filename) {
        if (auto file = fopen(filename, "r")) {
            fclose(file);
            return true;
        } else {
            return false;
        }
    }

    // TODO: use YAKL intrinsics for this to avoid needing to make host copies
    bool all_close(real2d &arr1, real2d &arr2, double tolerance) {
        int nx = arr1.dimension[0];
        int ny = arr2.dimension[1];
        auto arr1_h = arr1.createHostCopy();
        auto arr2_h = arr2.createHostCopy();
        for (int i=1; i<nx+1; i++) {
            for (int j=1; j<ny+1; j++) {
                if (abs(arr1_h(i,j) - arr2_h(i,j)) > tolerance || std::isnan(arr1_h(i,j) - arr2_h(i,j))) {
                    printf("arr1 = %f, arr2 = %f at %i,%i\n", arr1_h(i,j), arr2_h(i,j), i, j);
                    return false;
                }
            }
        }
        return true;
    }


    void dummy_atmos(
            std::string inputfile, 
            int ncol, real2d &p_lay, real2d &t_lay,
            real1d &sfc_alb_dir_vis, real1d &sfc_alb_dir_nir,
            real1d &sfc_alb_dif_vis, real1d &sfc_alb_dif_nir,
            real1d &mu0,
            real2d &lwp, real2d &iwp, real2d &rel, real2d &rei, real2d &cld) {

        // Setup boundary conditions, solar zenith angle, etc
        // NOTE: this stuff would come from the model in a real run

        // Ocean-ish values for surface albedos, just for example
        memset(sfc_alb_dir_vis , 0.06_wp );
        memset(sfc_alb_dir_nir , 0.06_wp );
        memset(sfc_alb_dif_vis , 0.06_wp );
        memset(sfc_alb_dif_nir , 0.06_wp );


        // Pick a solar zenith angle; this should come from the model
        memset(mu0, 0.86_wp );

        // Get dummy cloud PHYSICAL properties. Note that this function call
        // needs the CloudOptics object only because it uses the min and max
        // valid values from the lookup tables for liquid and ice water path to
        // create a dummy atmosphere.
        dummy_clouds(scream::rrtmgp::cloud_optics_sw, p_lay, t_lay, lwp, iwp, rel, rei, cld);
    }

    void dummy_clouds(
            CloudOptics &cloud_optics, real2d &p_lay, real2d &t_lay, 
            real2d &lwp, real2d &iwp, real2d &rel, real2d &rei, real2d &cloud_mask) {

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
        // NOTE: these "sane" values are in g/m2!
        parallel_for( SimpleBounds<2>(nlay,ncol) , YAKL_LAMBDA (int ilay, int icol) {
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
            real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dir,
            real2d &lw_flux_up, real2d &lw_flux_dn) {

        // Initialize netcdf reader
        simple_netcdf::SimpleNetCDF io;
        io.open(inputfile, NC_NOWRITE);

        // Initialize arrays to hold fluxes
        int nlev = io.getDimSize("lev");
        int ncol = io.getDimSize("col_flx");
        sw_flux_up  = real2d("sw_flux_up" , ncol, nlev);
        sw_flux_dn  = real2d("sw_flux_dn" , ncol, nlev);
        sw_flux_dir = real2d("sw_flux_dir", ncol, nlev);
        lw_flux_up  = real2d("lw_flux_up" , ncol, nlev);
        lw_flux_dn  = real2d("lw_flux_dn" , ncol, nlev);

        // Read data
        io.read(sw_flux_up,  "sw_flux_up" );
        io.read(sw_flux_dn,  "sw_flux_dn" );
        io.read(sw_flux_dir, "sw_flux_dir");
        io.read(lw_flux_up,  "lw_flux_up" );
        io.read(lw_flux_dn,  "lw_flux_dn" );
    }

    void write_fluxes(
            std::string outputfile,
            real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dir,
            real2d &lw_flux_up, real2d &lw_flux_dn) {

        simple_netcdf::SimpleNetCDF io;
        io.create(outputfile);
        io.write(sw_flux_up , "sw_flux_up" , {"col_flx","lev"});
        io.write(sw_flux_dn , "sw_flux_dn" , {"col_flx","lev"});
        io.write(sw_flux_dir, "sw_flux_dir", {"col_flx","lev"});
        io.write(lw_flux_up , "lw_flux_up" , {"col_flx","lev"});
        io.write(lw_flux_dn , "lw_flux_dn" , {"col_flx","lev"});
        io.close();

    }

    // TODO: This function should instead take values, not file names,
    // because for the test we do not want to write to file
    int compare(std::string file1, std::string file2) {
        // Read data from baseline and test file
        real2d sw_flux_up_1;
        real2d sw_flux_dn_1;
        real2d sw_flux_dir_1;
        real2d lw_flux_up_1;
        real2d lw_flux_dn_1;
        read_fluxes(file1, sw_flux_up_1, sw_flux_dn_1, sw_flux_dir_1, lw_flux_up_1, lw_flux_dn_1);
        real2d sw_flux_up_2;
        real2d sw_flux_dn_2;
        real2d sw_flux_dir_2;
        real2d lw_flux_up_2;
        real2d lw_flux_dn_2;
        read_fluxes(file2, sw_flux_up_2, sw_flux_dn_2, sw_flux_dir_2, lw_flux_up_2, lw_flux_dn_2);

        // Check values
        int nerr = 0;
        if (!rrtmgpTest::all_close(sw_flux_up_1    , sw_flux_up_2 , 0.001)) nerr++;
        if (!rrtmgpTest::all_close(sw_flux_dn_1    , sw_flux_dn_2 , 0.001)) nerr++;
        if (!rrtmgpTest::all_close(sw_flux_dir_1   , sw_flux_dir_2, 0.001)) nerr++;
        if (!rrtmgpTest::all_close(lw_flux_up_1    , lw_flux_up_2 , 0.001)) nerr++;
        if (!rrtmgpTest::all_close(lw_flux_dn_1    , lw_flux_dn_2 , 0.001)) nerr++;
        return nerr;
    }

}  // namespace rrtmgp
