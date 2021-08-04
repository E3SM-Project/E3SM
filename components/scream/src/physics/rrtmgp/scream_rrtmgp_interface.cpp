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

        OpticalProps2str get_cloud_optics_sw(const int ncol, const int nlay, CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist, real2d &p_lay, real2d &t_lay, real2d &lwp, real2d &iwp, real2d &rel, real2d &rei);
        OpticalProps1scl get_cloud_optics_lw(const int ncol, const int nlay, CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist, real2d &p_lay, real2d &t_lay, real2d &lwp, real2d &iwp, real2d &rel, real2d &rei);

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

        void compute_band_by_band_surface_albedos(
                const int ncol, const int nswbands,
                real1d &sfc_alb_dir_vis, real1d &sfc_alb_dir_nir,
                real1d &sfc_alb_dif_vis, real1d &sfc_alb_dif_nir,
                real2d &sfc_alb_dir,     real2d &sfc_alb_dif) {

            EKAT_ASSERT_MSG(initialized, "Error! rrtmgp_initialize must be called before GasOpticsRRTMGP object can be used.");
            auto wavenumber_limits = k_dist_sw.get_band_lims_wavenumber();

            EKAT_ASSERT_MSG(size(wavenumber_limits, 1) == 2,
                            "Error! 1st dimension for wavenumber_limits should be 2.");
            EKAT_ASSERT_MSG(size(wavenumber_limits, 2) == nswbands,
                            "Error! 2nd dimension for wavenumber_limits should be " + std::to_string(nswbands) + " (nswbands).");

            // Loop over bands, and determine for each band whether it is broadly in the
            // visible or infrared part of the spectrum (visible or "not visible")
            parallel_for(Bounds<2>(nswbands, ncol), YAKL_LAMBDA(const int ibnd, const int icol) {

             // Threshold between visible and infrared is 0.7 micron, or 14286 cm^-1.
             const real visible_wavenumber_threshold = 14286;

             // Wavenumber is in the visible if it is above the visible wavenumber
             // threshold, and in the infrared if it is below the threshold
             const bool is_visible_wave1 = (wavenumber_limits(1, ibnd) > visible_wavenumber_threshold ? true : false);
             const bool is_visible_wave2 = (wavenumber_limits(2, ibnd) > visible_wavenumber_threshold ? true : false);

              if (is_visible_wave1 && is_visible_wave2) {

                // Entire band is in the visible
                sfc_alb_dir(icol,ibnd) = sfc_alb_dir_vis(icol);
                sfc_alb_dif(icol,ibnd) = sfc_alb_dif_vis(icol);

              } else if (!is_visible_wave1 && !is_visible_wave2) {

                // Entire band is in the longwave (near-infrared)
                sfc_alb_dir(icol,ibnd) = sfc_alb_dir_nir(icol);
                sfc_alb_dif(icol,ibnd) = sfc_alb_dif_nir(icol);

              } else {

                // Band straddles the visible to near-infrared transition, so we take
                // the albedo to be the average of the visible and near-infrared
                // broadband albedos
                sfc_alb_dir(icol,ibnd) = 0.5*(sfc_alb_dir_vis(icol) + sfc_alb_dir_nir(icol));
                sfc_alb_dif(icol,ibnd) = 0.5*(sfc_alb_dif_vis(icol) + sfc_alb_dif_nir(icol));

              }
            });
        }

        void rrtmgp_main(
                const int ncol, const int nlay,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev, 
                GasConcs &gas_concs,
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
            OpticalProps2str clouds_sw = get_cloud_optics_sw(ncol, nlay, cloud_optics_sw, k_dist_sw, p_lay, t_lay, lwp, iwp, rel, rei);
            OpticalProps1scl clouds_lw = get_cloud_optics_lw(ncol, nlay, cloud_optics_lw, k_dist_lw, p_lay, t_lay, lwp, iwp, rel, rei);        

            // Do shortwave
            rrtmgp_sw(
                ncol, nlay,
                k_dist_sw, p_lay, t_lay, p_lev, t_lev, gas_concs, 
                sfc_alb_dir, sfc_alb_dif, mu0, clouds_sw, fluxes_sw
            );

            // Do longwave
            rrtmgp_lw(
                ncol, nlay,
                k_dist_lw, p_lay, t_lay, p_lev, t_lev, gas_concs,
                clouds_lw, fluxes_lw
            );
            
            // Calculate heating rates
        }

        OpticalProps2str get_cloud_optics_sw(
                const int ncol, const int nlay,
                CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist,
                real2d &p_lay, real2d &t_lay, real2d &lwp, real2d &iwp, real2d &rel, real2d &rei) {
 
            // Initialize optics
            OpticalProps2str clouds;
            clouds.init(kdist.get_band_lims_wavenumber());
            clouds.alloc_2str(ncol, nlay);

            // Needed for consistency with all-sky example problem?
            cloud_optics.set_ice_roughness(2);
 
            // Limit effective radii to be within bounds of lookup table
            auto rel_limited = real2d("rel_limited", ncol, nlay);
            auto rei_limited = real2d("rei_limited", ncol, nlay);
            limit_to_bounds(rel, cloud_optics.radliq_lwr, cloud_optics.radliq_upr, rel_limited);
            limit_to_bounds(rei, cloud_optics.radice_lwr, cloud_optics.radice_upr, rei_limited);

            // Calculate cloud optics
            cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel_limited, rei_limited, clouds);

            // Return optics
            return clouds;
        }


        OpticalProps1scl get_cloud_optics_lw(
                const int ncol, const int nlay,
                CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist, 
                real2d &p_lay, real2d &t_lay, real2d &lwp, real2d &iwp, real2d &rel, real2d &rei) {

            // Initialize optics
            OpticalProps1scl clouds;
            clouds.init(kdist.get_band_lims_wavenumber());
            clouds.alloc_1scl(ncol, nlay);  // this is dumb, why do we need to init and alloc separately?!

            // Needed for consistency with all-sky example problem?
            cloud_optics.set_ice_roughness(2);

            // Limit effective radii to be within bounds of lookup table
            auto rel_limited = real2d("rel_limited", ncol, nlay);
            auto rei_limited = real2d("rei_limited", ncol, nlay);
            limit_to_bounds(rel, cloud_optics.radliq_lwr, cloud_optics.radliq_upr, rel_limited);
            limit_to_bounds(rei, cloud_optics.radice_lwr, cloud_optics.radice_upr, rei_limited);

            // Calculate cloud optics
            cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel_limited, rei_limited, clouds);

            // Return optics
            return clouds;
        }


        void rrtmgp_sw(
                const int ncol, const int nlay,
                GasOpticsRRTMGP &k_dist,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs,
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0, OpticalProps2str &clouds,
                FluxesBroadband &fluxes) {

            // Get problem sizes
            int nbnd = k_dist.get_nband();
            int ngpt = k_dist.get_ngpt();
            int ngas = gas_concs.get_num_gases();

            // Get daytime indices
            auto dayIndices = int1d("dayIndices", ncol);
            memset(dayIndices, -1);
            // Loop below has to be done on host, so create host copies
            // TODO: there is probably a way to do this on the device
            auto dayIndices_h = dayIndices.createHostCopy();
            auto mu0_h = mu0.createHostCopy();
            int nday = 0;
            for (int icol = 1; icol <= ncol; icol++) {
                if (mu0_h(icol) > 0) {
                    nday++;
                    dayIndices_h(nday) = icol;
                }
            }
            // Copy data back to the device
            dayIndices_h.deep_copy_to(dayIndices);
            if (nday == 0) { 
                std::cout << "WARNING: no daytime columns found for this chunk!\n";
                return;
            }

            // Subset mu0
            auto mu0_day = real1d("mu0_day", nday);
            parallel_for(Bounds<1>(nday), YAKL_LAMBDA(int iday) {
                mu0_day(iday) = mu0(dayIndices(iday));
            });

            // subset state variables
            auto p_lay_day = real2d("p_lay_day", nday, nlay);
            auto t_lay_day = real2d("t_lay_day", nday, nlay);
            parallel_for(Bounds<2>(nlay,nday), YAKL_LAMBDA(int ilay, int iday) {
                p_lay_day(iday,ilay) = p_lay(dayIndices(iday),ilay);
                t_lay_day(iday,ilay) = t_lay(dayIndices(iday),ilay);
            });
            auto p_lev_day = real2d("p_lev_day", nday, nlay+1);
            auto t_lev_day = real2d("t_lev_day", nday, nlay+1);
            parallel_for(Bounds<2>(nlay+1,nday), YAKL_LAMBDA(int ilev, int iday) {
                p_lev_day(iday,ilev) = p_lev(dayIndices(iday),ilev);
                t_lev_day(iday,ilev) = t_lev(dayIndices(iday),ilev);
            });

            // Subset gases
            auto gas_names = gas_concs.get_gas_names();
            GasConcs gas_concs_day;
            gas_concs_day.init(gas_names, nday, nlay);
            for (int igas = 1; igas <= ngas; igas++) {
                auto vmr_day = real2d("vmr_day", nday, nlay);
                auto vmr     = real2d("vmr"    , ncol, nlay);
                gas_concs.get_vmr(gas_names(igas), vmr);
                parallel_for(Bounds<2>(nlay,nday), YAKL_LAMBDA(int ilay, int iday) {
                    vmr_day(iday,ilay) = vmr(dayIndices(iday),ilay);
                });
                gas_concs_day.set_vmr(gas_names(igas), vmr_day);
            }

            // Subset cloud optics
            OpticalProps2str clouds_day;
            clouds_day.init(k_dist.get_band_lims_wavenumber());
            clouds_day.alloc_2str(nday, nlay);
            parallel_for(Bounds<3>(nbnd,nlay,nday), YAKL_LAMBDA(int ibnd, int ilay, int iday) {
                clouds_day.tau(iday,ilay,ibnd) = clouds.tau(dayIndices(iday),ilay,ibnd);
                clouds_day.ssa(iday,ilay,ibnd) = clouds.ssa(dayIndices(iday),ilay,ibnd);
                clouds_day.g  (iday,ilay,ibnd) = clouds.g  (dayIndices(iday),ilay,ibnd);
            });

            // RRTMGP assumes surface albedos have a screwy dimension ordering
            // for some strange reason, so we need to transpose these; also do
            // daytime subsetting in the same kernel
            real2d sfc_alb_dir_T("sfc_alb_dir", nbnd, nday);
            real2d sfc_alb_dif_T("sfc_alb_dif", nbnd, nday);
            parallel_for(Bounds<2>(nbnd,nday), YAKL_LAMBDA(int ibnd, int icol) {
                sfc_alb_dir_T(ibnd,icol) = sfc_alb_dir(dayIndices(icol),ibnd);
                sfc_alb_dif_T(ibnd,icol) = sfc_alb_dif(dayIndices(icol),ibnd);
            });

            // Allocate space for optical properties
            OpticalProps2str optics;
            optics.alloc_2str(nday, nlay, k_dist);

            // Do gas optics
            real2d toa_flux("toa_flux", nday, ngpt);
            auto p_lay_host = p_lay.createHostCopy();
            bool top_at_1 = p_lay_host(1, 1) < p_lay_host(1, nlay);

            k_dist.gas_optics(nday, nlay, top_at_1, p_lay_day, p_lev_day, t_lay_day, gas_concs_day, optics, toa_flux);

            // Combine gas and cloud optics
            clouds_day.delta_scale();
            clouds_day.increment(optics);

            // Compute fluxes on daytime columns
            auto flux_up_day = real2d("flux_up_day", nday, nlay+1);
            auto flux_dn_day = real2d("flux_dn_day", nday, nlay+1);
            auto flux_dn_dir_day = real2d("flux_dn_dir_day", nday, nlay+1);
            FluxesBroadband fluxes_day;
            fluxes_day.flux_up     = flux_up_day; //real2d("flux_up"    , nday,nlay+1);
            fluxes_day.flux_dn     = flux_dn_day; //real2d("flux_dn"    , nday,nlay+1);
            fluxes_day.flux_dn_dir = flux_dn_dir_day; //real2d("flux_dn_dir", nday,nlay+1);
            rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);

            // Zero out all fluxes before expanding daytime fluxes
            auto &flux_up = fluxes.flux_up;
            auto &flux_dn = fluxes.flux_dn;
            auto &flux_dn_dir = fluxes.flux_dn_dir;
            parallel_for(Bounds<2>(nlay+1,ncol), YAKL_LAMBDA(int ilev, int icol) {
                flux_up    (icol,ilev) = 0;
                flux_dn    (icol,ilev) = 0;
                flux_dn_dir(icol,ilev) = 0;
            });
            
            // Expand daytime fluxes to all columns
            parallel_for(Bounds<2>(nlay+1,nday), YAKL_LAMBDA(int ilev, int iday) {
                int icol = dayIndices(iday);
                flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
                flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
                flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
            });
        }

        void rrtmgp_lw(
                const int ncol, const int nlay,
                GasOpticsRRTMGP &k_dist,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs,
                OpticalProps1scl &clouds,
                FluxesBroadband &fluxes) {

            // Problem size
            int nbnd = k_dist.get_nband();

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
            k_dist.gas_optics(ncol, nlay, top_at_1, p_lay, p_lev, t_lay, t_sfc, gas_concs, optics, lw_sources, real2d(), t_lev);

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
