#ifndef SCREAM_RRTMGP_INTERFACE_HPP
#define SCREAM_RRTMGP_INTERFACE_HPP

#include "cpp/rrtmgp/mo_gas_optics_rrtmgp.h"
#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"
#include "cpp/extensions/fluxes_byband/mo_fluxes_byband.h"
#include "cpp/rrtmgp_const.h"

#include "physics/share/physics_constants.hpp"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/logging/ekat_logger.hpp"

namespace scream {
  void yakl_init ();
  void yakl_finalize();
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
        extern void rrtmgp_initialize(GasConcs &gas_concs,
                                      const std::string coefficients_file_sw, const std::string coefficients_file_lw,
                                      const std::string cloud_optics_file_sw, const std::string cloud_optics_file_lw,
                                      const std::shared_ptr<spdlog::logger>& logger);
        /*
         * Compute band-by-band surface albedos from broadband albedos.
         */
        extern void compute_band_by_band_surface_albedos(
                const int ncol, const int nswbands,
                real1d &sfc_alb_dir_vis, real1d &sfc_alb_dir_nir,
                real1d &sfc_alb_dif_vis, real1d &sfc_alb_dif_nir,
                real2d &sfc_alb_dir,     real2d &sfc_alb_dif);
        /*
         * Compute broadband visible/UV and near-infrared surface fluxes.
         */
        extern void compute_broadband_surface_fluxes(
                const int ncol, const int ktop, const int nswbands,
                real3d &sw_bnd_flux_dir , real3d &sw_bnd_flux_dif ,
                real1d &sfc_flux_dir_vis, real1d &sfc_flux_dir_nir,
                real1d &sfc_flux_dif_vis, real1d &sfc_flux_dif_nir);
        /*
         * Main driver code to run RRTMGP.
         * The input logger is in charge of outputing info to
         * screen and/or to file (or neither), depending on how it was set up.
         */
        extern void rrtmgp_main(
                const int ncol, const int nlay,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs,
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0,
                real2d &lwp, real2d &iwp, real2d &rel, real2d &rei, real2d &cldfrac,
                real3d &aer_tau_sw, real3d &aer_ssa_sw, real3d &aer_asm_sw, real3d &aer_tau_lw,
                real3d &cld_tau_sw_gpt, real3d &cld_tau_lw_gpt,
                real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dn_dir,
                real2d &lw_flux_up, real2d &lw_flux_dn,
                real2d &sw_clrsky_flux_up, real2d &sw_clrsky_flux_dn, real2d &sw_clrsky_flux_dn_dir,
                real2d &lw_clrsky_flux_up, real2d &lw_clrsky_flux_dn,
                real3d &sw_bnd_flux_up, real3d &sw_bnd_flux_dn, real3d &sw_bnd_flux_dn_dir,
                real3d &lw_bnd_flux_up, real3d &lw_bnd_flux_dn,
                const Real tsi_scaling,
                const std::shared_ptr<spdlog::logger>& logger);
        /*
         * Perform any clean-up tasks
         */
        extern void rrtmgp_finalize();
        /*
         * Shortwave driver (called by rrtmgp_main)
         */
        extern void rrtmgp_sw(const int ncol, const int nlay,
                GasOpticsRRTMGP &k_dist,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs,
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0,
                OpticalProps2str &aerosol, OpticalProps2str &clouds,
                FluxesByband &fluxes, FluxesBroadband &clrsky_fluxes, const Real tsi_scaling,
                const std::shared_ptr<spdlog::logger>& logger);
        /*
         * Longwave driver (called by rrtmgp_main)
         */
        extern void rrtmgp_lw(
                const int ncol, const int nlay,
                GasOpticsRRTMGP &k_dist,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs,
                OpticalProps1scl &aerosol, OpticalProps1scl &clouds,
                FluxesByband &fluxes, FluxesBroadband &clrsky_fluxes);
        /*
         * Return a subcolumn mask consistent with a specified overlap assumption
         */
        int3d get_subcolumn_mask(const int ncol, const int nlay, const int ngpt, real2d &cldf, const int overlap_option, int1d &seeds);
        /*
         * Compute cloud area from 3d subcol cloud property
         */
        void compute_cloud_area(
                int ncol, int nlay, int ngpt, Real pmin, Real pmax,
                const real2d& pmid, const real3d& cld_tau_gpt, real1d& cld_area);

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
            yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay, ncol), YAKL_LAMBDA(int ilay, int icol) {
                // Compute in-cloud mixing ratio (mixing ratio of the cloudy part of the layer)
                // NOTE: these thresholds (from E3SM) seem arbitrary, but included here for consistency
                // This limits in-cloud mixing ratio to 0.005 kg/kg. According to note in cloud_diagnostics
                // in EAM, this is consistent with limits in MG2. Is this true for P3?
                if (cloud_fraction(icol,ilay) > 0) {
                    // Compute layer-integrated cloud mass (per unit area)
                    auto incloud_mixing_ratio = std::min(mixing_ratio(icol,ilay) / std::max(0.0001, cloud_fraction(icol,ilay)), 0.005);
                    cloud_mass(icol,ilay) = incloud_mixing_ratio * dp(icol,ilay) / physconst::gravit;
                } else {
                    cloud_mass(icol,ilay) = 0;
                }
            });
        }

        /*
         * Routine to limit a quantity to set bounds. Used to make sure
         * effective radii are within the bounds of the cloud optical
         * property look-up tables, but could be used to limit other
         * fields as well.
         */
        template<class S, class T> void limit_to_bounds(S const &arr_in, T const lower, T const upper, S &arr_out) {
            yakl::c::parallel_for(arr_in.totElems(), YAKL_LAMBDA(int i) {
                arr_out.data()[i] = std::min(std::max(arr_in.data()[i], lower), upper);
            });
        }


    } // namespace rrtmgp
}  // namespace scream

#endif  // SCREAM_RRTMGP_INTERFACE_HPP
