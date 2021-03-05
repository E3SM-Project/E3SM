#include "scream_rrtmgp_interface.hpp"
#include "mo_load_coefficients.h"
#include "mo_load_cloud_coefficients.h"
#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "cpp/rrtmgp/mo_gas_optics_rrtmgp.h"
#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"
#include "cpp/rte/mo_rte_sw.h"
#include "cpp/rte/mo_rte_lw.h"
#include "cpp/const.h"

namespace scream {
    namespace rrtmgp {

        OpticalProps2str get_cloud_optics_sw(CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist, real2d &p_lay, real2d &t_lay, real2d &lwp, real2d &iwp, real2d &rel, real2d &rei);
        OpticalProps1scl get_cloud_optics_lw(CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist, real2d &p_lay, real2d &t_lay, real2d &lwp, real2d &iwp, real2d &rel, real2d &rei);

        /*
         * Names of input files we will need.
         */
        std::string coefficients_file_sw = "./data/rrtmgp-data-sw-g224-2018-12-04.nc";
        std::string coefficients_file_lw = "./data/rrtmgp-data-lw-g256-2018-12-04.nc";
        std::string cloud_optics_file_sw = "./data/rrtmgp-cloud-optics-coeffs-sw.nc";
        std::string cloud_optics_file_lw = "./data/rrtmgp-cloud-optics-coeffs-lw.nc";

        /* 
         * Objects containing k-distribution information need to be initialized
         * once and then persist throughout the life of the program, so we
         * declare them here within the rrtmgp namespace.
         */
        GasOpticsRRTMGP k_dist_sw;
        GasOpticsRRTMGP k_dist_lw;

        /*
         * Objects containing cloud optical property look-up table information.
         * We want to initialize these once and use throughout the life of the
         * program, so declare here and read data in during rrtmgp_initialize().
         */
        CloudOptics cloud_optics_sw;
        CloudOptics cloud_optics_lw;

        /*
         * Gas concentrations. We want to initialize this once, since it will
         * determine the absorbing gases that we keep track of in the k-dist
         * objects. So we will initialize once, and then just update values
         * at runtime.
         */
        //GasConcs gas_concs;

        bool initialized = false;

        /*
         * The following routines provide a simple interface to RRTMGP. These
         * can be used as-is, but are intended to be wrapped by the SCREAM AD
         * interface to radiation.
         */
        void rrtmgp_initialize(GasConcs &gas_concs) {

            // If we've already initialized, just exit
            if (initialized) { 
                std::cout << "RRTMGP is already initialized; skipping\n";
                return; 
            }

            // Initialize YAKL
            if (!yakl::isInitialized()) {
                yakl::init();
            }

            // Load and initialize absorption coefficient data
            load_and_init(k_dist_sw, coefficients_file_sw, gas_concs);
            load_and_init(k_dist_lw, coefficients_file_lw, gas_concs);

            // Load and initialize cloud optical property look-up table information
            load_cld_lutcoeff(cloud_optics_sw, cloud_optics_file_sw);
            load_cld_lutcoeff(cloud_optics_lw, cloud_optics_file_lw);

            // We are now initialized!
            initialized = true;
        }

        void rrtmgp_finalize() {
            initialized = false;
            k_dist_sw.finalize();
            k_dist_lw.finalize();
            cloud_optics_sw.finalize(); //~CloudOptics();
            cloud_optics_lw.finalize(); //~CloudOptics();
        }

        void rrtmgp_main(
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev, 
                GasConcs &gas_concs, real2d &col_dry,
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0, 
                real2d &lwp, real2d &iwp, real2d &rel, real2d &rei,
                real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dn_dir,
                real2d &lw_flux_up, real2d &lw_flux_dn) {

            // Setup pointers to RRTMGP SW fluxes
            FluxesBroadband fluxes_sw;
            fluxes_sw.flux_up = sw_flux_up;
            fluxes_sw.flux_dn = sw_flux_dn;
            fluxes_sw.flux_dn_dir = sw_flux_dn_dir;

            // Setup pointers to RRTMGP LW fluxes
            FluxesBroadband fluxes_lw;
            fluxes_lw.flux_up = lw_flux_up;
            fluxes_lw.flux_dn = lw_flux_dn;

            // Convert cloud physical properties to optical properties for input to RRTMGP
            OpticalProps2str clouds_sw = get_cloud_optics_sw(cloud_optics_sw, k_dist_sw, p_lay, t_lay, lwp, iwp, rel, rei);
            OpticalProps1scl clouds_lw = get_cloud_optics_lw(cloud_optics_lw, k_dist_lw, p_lay, t_lay, lwp, iwp, rel, rei);

            // Do shortwave
            rrtmgp_sw(
                k_dist_sw, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry, 
                sfc_alb_dir, sfc_alb_dif, mu0, clouds_sw, fluxes_sw
            );

            // Do longwave
            rrtmgp_lw(
                k_dist_lw, p_lay, t_lay, p_lev, t_lev, gas_concs, col_dry,
                clouds_lw, fluxes_lw
            );
            
            // Calculate heating rates
        }

        OpticalProps2str get_cloud_optics_sw(
                CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist, 
                real2d &p_lay, real2d &t_lay, real2d &lwp, real2d &iwp, real2d &rel, real2d &rei) {

            // Problem sizes
            int ncol = t_lay.dimension[0];
            int nlay = t_lay.dimension[1];
 
            // Initialize optics
            OpticalProps2str clouds;
            clouds.init(kdist.get_band_lims_wavenumber());
            clouds.alloc_2str(ncol, nlay);

            // Needed for consistency with all-sky example problem?
            cloud_optics.set_ice_roughness(2);

            // Calculate cloud optics
            cloud_optics.cloud_optics(lwp, iwp, rel, rei, clouds);

            // Return optics
            return clouds;
        }


        OpticalProps1scl get_cloud_optics_lw(
                CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist, 
                real2d &p_lay, real2d &t_lay, real2d &lwp, real2d &iwp, real2d &rel, real2d &rei) {

            // Problem sizes
            int ncol = t_lay.dimension[0];
            int nlay = t_lay.dimension[1];

            // Initialize optics
            OpticalProps1scl clouds;
            clouds.init(kdist.get_band_lims_wavenumber());
            clouds.alloc_1scl(ncol, nlay);  // this is dumb, why do we need to init and alloc separately?!

            // Needed for consistency with all-sky example problem?
            cloud_optics.set_ice_roughness(2);

            // Calculate cloud optics
            cloud_optics.cloud_optics(lwp, iwp, rel, rei, clouds);

            // Return optics
            return clouds;
        }


        void rrtmgp_sw(
                GasOpticsRRTMGP &k_dist, 
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev, 
                GasConcs &gas_concs, real2d &col_dry,
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0, OpticalProps2str &clouds,
                FluxesBroadband &fluxes) {

            // Get problem sizes
            int nbnd = k_dist.get_nband();
            int ngpt = k_dist.get_ngpt();
            int ncol = p_lay.dimension[0];
            int nlay = p_lay.dimension[1];

            // Allocate space for optical properties
            OpticalProps2str optics;
            optics.alloc_2str(ncol, nlay, k_dist);

            // Do gas optics
            real2d toa_flux("toa_flux", ncol, ngpt);
            auto p_lay_host = p_lay.createHostCopy();
            bool top_at_1 = p_lay_host(1, 1) < p_lay_host(1, nlay);

            k_dist.gas_optics(top_at_1, p_lay, p_lev, t_lay, gas_concs, optics, toa_flux);

            // Combine gas and cloud optics
            clouds.delta_scale();
            clouds.increment(optics);

            // Compute fluxes
            rte_sw(optics, top_at_1, mu0, toa_flux, sfc_alb_dir, sfc_alb_dif, fluxes);
        }

        void rrtmgp_lw(
                GasOpticsRRTMGP &k_dist,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs, real2d &col_dry,
                OpticalProps1scl &clouds,
                FluxesBroadband &fluxes) {

            // Problem size
            int nbnd = k_dist.get_nband();
            int ngpt = k_dist.get_ngpt();
            int ncol = p_lay.dimension[0];
            int nlay = p_lay.dimension[1];

            // Allocate space for optical properties
            OpticalProps1scl optics;
            optics.alloc_1scl(ncol, nlay, k_dist);

            // Boundary conditions
            SourceFuncLW lw_sources;
            lw_sources.alloc(ncol, nlay, k_dist);
            real1d t_sfc   ("t_sfc"        ,ncol);
            real2d emis_sfc("emis_sfc",nbnd,ncol);

            // Surface temperature
            auto p_lay_host = p_lay.createHostCopy();
            bool top_at_1 = p_lay_host(1, 1) < p_lay_host(1, nlay);
            auto t_lev_host = t_lev.createHostCopy();
            memset( t_sfc    , t_lev_host(1, merge(nlay+1, 1, top_at_1)) );
            memset( emis_sfc , 0.98_wp                                   );

            // Do gas optics
            k_dist.gas_optics(top_at_1, p_lay, p_lev, t_lay, t_sfc, gas_concs, optics, lw_sources, real2d(), t_lev);

            // Combine gas and cloud optics
            clouds.increment(optics);

            // Get Gaussian quadrature weights
            // TODO: move this crap out of userland!
            // Weights and angle secants for first order (k=1) Gaussian quadrature.
            //   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
            //   after Abramowitz & Stegun 1972, page 921
            int constexpr max_gauss_pts = 4;
            realHost2d gauss_Ds_host ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
            gauss_Ds_host(1,1) = 1.66_wp      ; gauss_Ds_host(2,1) =         0._wp; gauss_Ds_host(3,1) =         0._wp; gauss_Ds_host(4,1) =         0._wp;
            gauss_Ds_host(1,2) = 1.18350343_wp; gauss_Ds_host(2,2) = 2.81649655_wp; gauss_Ds_host(3,2) =         0._wp; gauss_Ds_host(4,2) =         0._wp;
            gauss_Ds_host(1,3) = 1.09719858_wp; gauss_Ds_host(2,3) = 1.69338507_wp; gauss_Ds_host(3,3) = 4.70941630_wp; gauss_Ds_host(4,3) =         0._wp;
            gauss_Ds_host(1,4) = 1.06056257_wp; gauss_Ds_host(2,4) = 1.38282560_wp; gauss_Ds_host(3,4) = 2.40148179_wp; gauss_Ds_host(4,4) = 7.15513024_wp;

            realHost2d gauss_wts_host("gauss_wts",max_gauss_pts,max_gauss_pts);
            gauss_wts_host(1,1) = 0.5_wp         ; gauss_wts_host(2,1) = 0._wp          ; gauss_wts_host(3,1) = 0._wp          ; gauss_wts_host(4,1) = 0._wp          ;
            gauss_wts_host(1,2) = 0.3180413817_wp; gauss_wts_host(2,2) = 0.1819586183_wp; gauss_wts_host(3,2) = 0._wp          ; gauss_wts_host(4,2) = 0._wp          ;
            gauss_wts_host(1,3) = 0.2009319137_wp; gauss_wts_host(2,3) = 0.2292411064_wp; gauss_wts_host(3,3) = 0.0698269799_wp; gauss_wts_host(4,3) = 0._wp          ;
            gauss_wts_host(1,4) = 0.1355069134_wp; gauss_wts_host(2,4) = 0.2034645680_wp; gauss_wts_host(3,4) = 0.1298475476_wp; gauss_wts_host(4,4) = 0.0311809710_wp;

            real2d gauss_Ds ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
            real2d gauss_wts("gauss_wts",max_gauss_pts,max_gauss_pts);
            gauss_Ds_host .deep_copy_to(gauss_Ds );
            gauss_wts_host.deep_copy_to(gauss_wts);

            // Compute fluxes
            rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics, top_at_1, lw_sources, emis_sfc, fluxes);

        }

    }  // namespace rrtmgp
}  // namespace scream
