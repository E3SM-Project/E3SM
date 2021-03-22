#ifndef SCREAM_RRTMGP_INTERFACE_HPP
#define SCREAM_RRTMGP_INTERFACE_HPP

#include "cpp/rrtmgp/mo_gas_optics_rrtmgp.h"
#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"
#include "cpp/rte/mo_fluxes.h"

namespace scream {
    namespace rrtmgp {
        /* 
         * Objects containing k-distribution information need to be initialized
         * once and then persist throughout the life of the program, so we
         * declare them here within the rrtmgp namespace.
         */
        extern GasOpticsRRTMGP k_dist_sw;
        extern GasOpticsRRTMGP k_dist_lw;
        /*
         * Objects containing cloud optical property look-up table information.
         * We want to initialize these once and use throughout the life of the
         * program, so declare here and read data in during rrtmgp_initialize().
         */
        extern CloudOptics cloud_optics_sw;
        extern CloudOptics cloud_optics_lw;
        /*
         * Flag to indicate whether or not we have initialized RRTMGP
         */
        extern bool initialized;
        /*
         * Initialize data for RRTMGP driver
         */
        extern void rrtmgp_initialize(GasConcs &gas_concs);
        /*
         * Main driver code to run RRTMGP
         */
        extern void rrtmgp_main(
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev, 
                GasConcs &gas_concs,
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0, 
                real2d &lwp, real2d &iwp, real2d &real, real2d &rei,
                real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dn_dir,
                real2d &lw_flux_up, real2d &lw_flux_dn);
        /*
         * Perform any clean-up tasks
         */
        extern void rrtmgp_finalize();
        /*
         * Shortwave driver (called by rrtmgp_main)
         */
        extern void rrtmgp_sw(
                GasOpticsRRTMGP &k_dist, 
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev, 
                GasConcs &gas_concs, 
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0, OpticalProps2str &clouds,
                FluxesBroadband &fluxes);
        /*
         * Longwave driver (called by rrtmgp_main)
         */
        extern void rrtmgp_lw(
                GasOpticsRRTMGP &k_dist,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs,
                OpticalProps1scl &clouds,
                FluxesBroadband &fluxes);
    } // namespace rrtmgp
}  // namespace scream

#endif  // SCREAM_RRTMGP_INTERFACE_HPP
