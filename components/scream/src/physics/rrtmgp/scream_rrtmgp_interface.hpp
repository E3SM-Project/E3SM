#ifndef SCREAM_RRTMGP_INTERFACE_HPP
#define SCREAM_RRTMGP_INTERFACE_HPP

#include "cpp/rrtmgp/mo_gas_optics_rrtmgp.h"
#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"
#include "cpp/rte/mo_fluxes.h"
#include "cpp/const.h"
#include "physics/share/physics_constants.hpp"

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
                const int ncol, const int nlay,
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
                const int ncol, const int nlay,
                GasOpticsRRTMGP &k_dist, 
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev, 
                GasConcs &gas_concs, 
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0, OpticalProps2str &clouds,
                FluxesBroadband &fluxes);
        /*
         * Longwave driver (called by rrtmgp_main)
         */
        extern void rrtmgp_lw(
                const int ncol, const int nlay,
                GasOpticsRRTMGP &k_dist,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs,
                OpticalProps1scl &clouds,
                FluxesBroadband &fluxes);
        /* 
         * Provide a function to convert cloud (water and ice) mixing ratios to layer mass per unit area
         * (what E3SM refers to as "in-cloud water paths", a terminology we shun here to avoid confusion
         * with the standard practice of using "water path" to refer to the total column-integrated
         * quantities).
         */
        template<class T, int myMem, int myStyle> void mixing_ratio_to_cloud_mass(
                yakl::Array<T,2,myMem,myStyle> const &mixing_ratio, 
                yakl::Array<T,2,myMem,myStyle> const &cloud_fraction, 
                yakl::Array<T,2,myMem,myStyle> const &dp, 
                yakl::Array<T,2,myMem,myStyle>       &cloud_mass) {
            int ncol = mixing_ratio.dimension[0];
            int nlay = mixing_ratio.dimension[1];
            using physconst = scream::physics::Constants<Real>;
            parallel_for(Bounds<2>(nlay, ncol), YAKL_LAMBDA(int ilay, int icol) {
                // Compute in-cloud mixing ratio (mixing ratio of the cloudy part of the layer)
                // NOTE: these thresholds (from E3SM) seem arbitrary, but included here for consistency
                // This limits in-cloud mixing ratio to 0.005 kg/kg. According to note in cloud_diagnostics
                // in EAM, this is consistent with limits in MG2. Is this true for P3?
                // Note also that this limits cloud fraction in the computation to 0.0001, to avoid
                // dividing by zero. But if cloud fraction is zero, should we not just make incloud mixing
                // ration zero as well? In that case, we would probably expect mixing_ratio to also be zero.
                auto incloud_mixing_ratio = std::min(mixing_ratio(icol,ilay) / std::max(0.0001, cloud_fraction(icol,ilay)), 0.005);
                // Compute layer-integrated cloud mass (per unit area)
                cloud_mass(icol,ilay) = incloud_mixing_ratio * dp(icol,ilay) / physconst::gravit;
            });
        }

    } // namespace rrtmgp
}  // namespace scream

#endif  // SCREAM_RRTMGP_INTERFACE_HPP
