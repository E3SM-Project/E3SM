#include "scream_rrtmgp_interface.hpp"
#include "rrtmgp_utils.hpp"
#include "mo_load_coefficients.h"
#include "mo_load_cloud_coefficients.h"
#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "cpp/rrtmgp/mo_gas_optics_rrtmgp.h"
#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"
#include "cpp/rte/mo_rte_sw.h"
#include "cpp/rte/mo_rte_lw.h"

namespace scream {
    void yakl_init ()
    {
      // Initialize yakl
      if(!yakl::isInitialized()) { yakl::init(); }
    }

    void yakl_finalize()
    {
      // Finalize YAKL
      yakl::finalize();
    }
    namespace rrtmgp {

        using yakl::fortran::parallel_for;
        using yakl::fortran::SimpleBounds;
        using yakl::intrinsics::merge;

        OpticalProps2str get_cloud_optics_sw(const int ncol, const int nlay, CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist, real2d &lwp, real2d &iwp, real2d &rel, real2d &rei);
        OpticalProps1scl get_cloud_optics_lw(const int ncol, const int nlay, CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist, real2d &lwp, real2d &iwp, real2d &rel, real2d &rei);
        OpticalProps2str get_subsampled_clouds(
                const int ncol, const int nlay, const int nbnd, const int ngpt,
                OpticalProps2str &cloud_optics, GasOpticsRRTMGP &kdist, real2d &cld, real2d &p_lay);
        OpticalProps1scl get_subsampled_clouds(
                const int ncol, const int nlay, const int nbnd, const int ngpt,
                OpticalProps1scl &cloud_optics, GasOpticsRRTMGP &kdist, real2d &cld, real2d &p_lay);

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
        void rrtmgp_initialize(GasConcs &gas_concs,
                               const std::string coefficients_file_sw, const std::string coefficients_file_lw,
                               const std::string cloud_optics_file_sw, const std::string cloud_optics_file_lw,
                               const std::shared_ptr<spdlog::logger>& logger) {

            // If we've already initialized, just exit
            if (initialized) {
                if (logger)
                  logger->info("RRTMGP is already initialized; skipping\n");
                return; 
            }

            // Initialize YAKL
            if (!yakl::isInitialized()) {  yakl::init(); }

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

            EKAT_ASSERT_MSG(yakl::intrinsics::size(wavenumber_limits, 1) == 2,
                            "Error! 1st dimension for wavenumber_limits should be 2.");
            EKAT_ASSERT_MSG(yakl::intrinsics::size(wavenumber_limits, 2) == nswbands,
                            "Error! 2nd dimension for wavenumber_limits should be " + std::to_string(nswbands) + " (nswbands).");

            // Loop over bands, and determine for each band whether it is broadly in the
            // visible or infrared part of the spectrum (visible or "not visible")
            parallel_for(SimpleBounds<2>(nswbands, ncol), YAKL_LAMBDA(const int ibnd, const int icol) {

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


        void compute_broadband_surface_fluxes(
                const int ncol, const int ktop, const int nswbands,
                real3d &sw_bnd_flux_dir , real3d &sw_bnd_flux_dif ,
                real1d &sfc_flux_dir_vis, real1d &sfc_flux_dir_nir,
                real1d &sfc_flux_dif_vis, real1d &sfc_flux_dif_nir) {
            // Band 10 straddles the near-IR and visible, so divide contributions from band 10 between both broadband sums
            // TODO: Hard-coding these band indices is really bad practice. If the bands ever were to change (like when
            // the RRTMG bands were re-ordered for RRTMGP), we would be using the wrong bands for the IR and UV/VIS. This
            // should be refactored to grab the correct bands by specifying appropriate wavenumber rather than index.
            //sfc_flux_dir_nir(i) = sum(sw_bnd_flux_dir(i+1,kbot,1:9))   + 0.5 * sw_bnd_flux_dir(i+1,kbot,10);
            //sfc_flux_dir_vis(i) = sum(sw_bnd_flux_dir(i+1,kbot,11:14)) + 0.5 * sw_bnd_flux_dir(i+1,kbot,10);
            //sfc_flux_dif_nir(i) = sum(sw_bnd_flux_dif(i+1,kbot,1:9))   + 0.5 * sw_bnd_flux_dif(i+1,kbot,10);
            //sfc_flux_dif_vis(i) = sum(sw_bnd_flux_dif(i+1,kbot,11:14)) + 0.5 * sw_bnd_flux_dif(i+1,kbot,10);

            // Initialize sums over bands
            memset(sfc_flux_dir_nir, 0);
            memset(sfc_flux_dir_vis, 0);
            memset(sfc_flux_dif_nir, 0);
            memset(sfc_flux_dif_vis, 0);
 
            // Threshold between visible and infrared is 0.7 micron, or 14286 cm^-1.
            const real visible_wavenumber_threshold = 14286;
            auto wavenumber_limits = k_dist_sw.get_band_lims_wavenumber();
            parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA(const int icol) {
              for (int ibnd = 1; ibnd <= nswbands; ++ibnd) {
                // Wavenumber is in the visible if it is above the visible wavenumber
                // threshold, and in the infrared if it is below the threshold
                const bool is_visible_wave1 = (wavenumber_limits(1, ibnd) > visible_wavenumber_threshold ? true : false);
                const bool is_visible_wave2 = (wavenumber_limits(2, ibnd) > visible_wavenumber_threshold ? true : false);

                if (is_visible_wave1 && is_visible_wave2) {

                    // Entire band is in the visible
                    sfc_flux_dir_vis(icol) += sw_bnd_flux_dir(icol,ktop,ibnd);
                    sfc_flux_dif_vis(icol) += sw_bnd_flux_dif(icol,ktop,ibnd);

                } else if (!is_visible_wave1 && !is_visible_wave2) {

                    // Entire band is in the longwave (near-infrared)
                    sfc_flux_dir_nir(icol) += sw_bnd_flux_dir(icol,ktop,ibnd);
                    sfc_flux_dif_nir(icol) += sw_bnd_flux_dif(icol,ktop,ibnd);

                } else {

                    // Band straddles the visible to near-infrared transition, so put half 
                    // the flux in visible and half in near-infrared fluxes
                    sfc_flux_dir_vis(icol) += 0.5 * sw_bnd_flux_dir(icol,ktop,ibnd);
                    sfc_flux_dif_vis(icol) += 0.5 * sw_bnd_flux_dif(icol,ktop,ibnd);
                    sfc_flux_dir_nir(icol) += 0.5 * sw_bnd_flux_dir(icol,ktop,ibnd);
                    sfc_flux_dif_nir(icol) += 0.5 * sw_bnd_flux_dif(icol,ktop,ibnd);
                }
              }
            });
        }

 
        void rrtmgp_main(
                const int ncol, const int nlay,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs,
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0,
                real2d &lwp, real2d &iwp, real2d &rel, real2d &rei, real2d &cldfrac,
                real3d &aer_tau_sw, real3d &aer_ssa_sw, real3d &aer_asm_sw, real3d &aer_tau_lw,
                real3d &cld_tau_sw_gpt,
                real3d &cld_tau_lw_gpt,
                real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dn_dir,
                real2d &lw_flux_up, real2d &lw_flux_dn,
                real2d &sw_clrsky_flux_up, real2d &sw_clrsky_flux_dn, real2d &sw_clrsky_flux_dn_dir,
                real2d &lw_clrsky_flux_up, real2d &lw_clrsky_flux_dn,
                real3d &sw_bnd_flux_up, real3d &sw_bnd_flux_dn, real3d &sw_bnd_flux_dn_dir,
                real3d &lw_bnd_flux_up, real3d &lw_bnd_flux_dn,
                const Real tsi_scaling,
                const std::shared_ptr<spdlog::logger>& logger) {

#ifdef SCREAM_RRTMGP_DEBUG
            // Sanity check inputs, and possibly repair
            check_range(t_lay      ,  k_dist_sw.get_temp_min(),         k_dist_sw.get_temp_max(), "rrtmgp_main::t_lay");
            check_range(t_lev      ,  k_dist_sw.get_temp_min(),         k_dist_sw.get_temp_max(), "rrtmgp_main::t_lev");
            check_range(p_lay      , k_dist_sw.get_press_min(),        k_dist_sw.get_press_max(), "rrtmgp_main::p_lay");
            check_range(p_lev      , k_dist_sw.get_press_min(),        k_dist_sw.get_press_max(), "rrtmgp_main::p_lev");
            check_range(sfc_alb_dir,                         0,                                1, "rrtmgp_main::sfc_alb_dir");
            check_range(sfc_alb_dif,                         0,                                1, "rrtmgp_main::sfc_alb_dif");
            check_range(mu0        ,                         0,                                1, "rrtmgp_main::mu0");
            check_range(lwp        ,                         0, std::numeric_limits<Real>::max(), "rrtmgp_main::lwp");
            check_range(iwp        ,                         0, std::numeric_limits<Real>::max(), "rrtmgp_main::iwp");
            check_range(rel        ,                         0, std::numeric_limits<Real>::max(), "rrtmgp_main::rel");
            check_range(rei        ,                         0, std::numeric_limits<Real>::max(), "rrtmgp_main::rei");
#endif

            // Setup pointers to RRTMGP SW fluxes
            FluxesByband fluxes_sw;
            fluxes_sw.flux_up = sw_flux_up;
            fluxes_sw.flux_dn = sw_flux_dn;
            fluxes_sw.flux_dn_dir = sw_flux_dn_dir;
            fluxes_sw.bnd_flux_up = sw_bnd_flux_up;
            fluxes_sw.bnd_flux_dn = sw_bnd_flux_dn;
            fluxes_sw.bnd_flux_dn_dir = sw_bnd_flux_dn_dir;
            // Clear-sky
            FluxesBroadband clrsky_fluxes_sw;
            clrsky_fluxes_sw.flux_up = sw_clrsky_flux_up;
            clrsky_fluxes_sw.flux_dn = sw_clrsky_flux_dn;
            clrsky_fluxes_sw.flux_dn_dir = sw_clrsky_flux_dn_dir;

            // Setup pointers to RRTMGP LW fluxes
            FluxesByband fluxes_lw;
            fluxes_lw.flux_up = lw_flux_up;
            fluxes_lw.flux_dn = lw_flux_dn;
            fluxes_lw.bnd_flux_up = lw_bnd_flux_up;
            fluxes_lw.bnd_flux_dn = lw_bnd_flux_dn;
            // Clear-sky
            FluxesBroadband clrsky_fluxes_lw;
            clrsky_fluxes_lw.flux_up = lw_clrsky_flux_up;
            clrsky_fluxes_lw.flux_dn = lw_clrsky_flux_dn;

            auto nswbands = k_dist_sw.get_nband();
            auto nlwbands = k_dist_lw.get_nband();

            // Setup aerosol optical properties
            OpticalProps2str aerosol_sw;
            OpticalProps1scl aerosol_lw;
            aerosol_sw.init(k_dist_sw.get_band_lims_wavenumber());
            aerosol_sw.alloc_2str(ncol, nlay);
            parallel_for(SimpleBounds<3>(nswbands,nlay,ncol) , YAKL_LAMBDA (int ibnd, int ilay, int icol) {
                aerosol_sw.tau(icol,ilay,ibnd) = aer_tau_sw(icol,ilay,ibnd);
                aerosol_sw.ssa(icol,ilay,ibnd) = aer_ssa_sw(icol,ilay,ibnd);
                aerosol_sw.g  (icol,ilay,ibnd) = aer_asm_sw(icol,ilay,ibnd);
            });
            aerosol_lw.init(k_dist_lw.get_band_lims_wavenumber());
            aerosol_lw.alloc_1scl(ncol, nlay);
            parallel_for(SimpleBounds<3>(nlwbands,nlay,ncol) , YAKL_LAMBDA (int ibnd, int ilay, int icol) {
                aerosol_lw.tau(icol,ilay,ibnd) = aer_tau_lw(icol,ilay,ibnd);
            });

#ifdef SCREAM_RRTMGP_DEBUG
            // Check aerosol optical properties
            // NOTE: these should already have been checked by precondition checks, but someday we might have
            // non-trivial aerosol optics, so this is still good to do here.
            check_range(aerosol_sw.tau,  0, 1e3, "rrtmgp_main:aerosol_sw.tau");
            check_range(aerosol_sw.ssa,  0,   1, "rrtmgp_main:aerosol_sw.ssa"); //, "aerosol_optics_sw.ssa");
            check_range(aerosol_sw.g  , -1,   1, "rrtmgp_main:aerosol_sw.g  "); //, "aerosol_optics_sw.g"  );
            check_range(aerosol_lw.tau,  0, 1e3, "rrtmgp_main:aerosol_lw.tau");
#endif


            // Convert cloud physical properties to optical properties for input to RRTMGP
            OpticalProps2str clouds_sw = get_cloud_optics_sw(ncol, nlay, cloud_optics_sw, k_dist_sw, lwp, iwp, rel, rei);
            OpticalProps1scl clouds_lw = get_cloud_optics_lw(ncol, nlay, cloud_optics_lw, k_dist_lw, lwp, iwp, rel, rei);        
            // Do subcolumn sampling to map bands -> gpoints based on cloud fraction and overlap assumption;
            // This implements the Monte Carlo Independing Column Approximation by mapping only a single 
            // subcolumn (cloud state) to each gpoint.
            auto nswgpts = k_dist_sw.get_ngpt();
            auto clouds_sw_gpt = get_subsampled_clouds(ncol, nlay, nswbands, nswgpts, clouds_sw, k_dist_sw, cldfrac, p_lay);
            // Longwave
            auto nlwgpts = k_dist_lw.get_ngpt();
            auto clouds_lw_gpt = get_subsampled_clouds(ncol, nlay, nlwbands, nlwgpts, clouds_lw, k_dist_lw, cldfrac, p_lay);

            // Copy cloud properties to outputs (is this needed, or can we just use pointers?)
            // Alternatively, just compute and output a subcolumn cloud mask
            parallel_for(SimpleBounds<3>(nswgpts, nlay, ncol), YAKL_LAMBDA (int igpt, int ilay, int icol) {
                cld_tau_sw_gpt(icol,ilay,igpt) = clouds_sw_gpt.tau(icol,ilay,igpt);
            });
            parallel_for(SimpleBounds<3>(nlwgpts, nlay, ncol), YAKL_LAMBDA (int igpt, int ilay, int icol) {
                cld_tau_lw_gpt(icol,ilay,igpt) = clouds_lw_gpt.tau(icol,ilay,igpt);
            });

#ifdef SCREAM_RRTMGP_DEBUG
            // Perform checks on optics; these would be caught by RRTMGP_EXPENSIVE_CHECKS in the RRTMGP code,
            // but we might want to provide additional debug info here. NOTE: we may actually want to move this
            // up higher in the code, I think optical props should go up higher since optical props are kind of
            // a parameterization of their own, and we might want to swap different choices. These checks go here
            // only because we need to run them on computed optical props, so if the optical props themselves get
            // computed up higher, then perform these checks higher as well
            check_range(clouds_sw.tau,  0, std::numeric_limits<Real>::max(), "rrtmgp_main:clouds_sw.tau");
            check_range(clouds_sw.ssa,  0,                                1, "rrtmgp_main:clouds_sw.ssa");
            check_range(clouds_sw.g  , -1,                                1, "rrtmgp_main:clouds_sw.g  ");
            check_range(clouds_sw.tau,  0, std::numeric_limits<Real>::max(), "rrtmgp_main:clouds_sw.tau");
#endif

            // Do shortwave
            rrtmgp_sw(
                ncol, nlay,
                k_dist_sw, p_lay, t_lay, p_lev, t_lev, gas_concs, 
                sfc_alb_dir, sfc_alb_dif, mu0, aerosol_sw, clouds_sw_gpt,
                fluxes_sw, clrsky_fluxes_sw,
                tsi_scaling, logger
            );

            // Do longwave
            rrtmgp_lw(
                ncol, nlay,
                k_dist_lw, p_lay, t_lay, p_lev, t_lev, gas_concs,
                aerosol_lw, clouds_lw_gpt,
                fluxes_lw, clrsky_fluxes_lw
            );
            
        }

        int3d get_subcolumn_mask(const int ncol, const int nlay, const int ngpt, real2d &cldf, const int overlap_option, int1d &seeds) {
            
            // Routine will return subcolumn mask with values of 0 indicating no cloud, 1 indicating cloud
            auto subcolumn_mask = int3d("subcolumn_mask", ncol, nlay, ngpt);

            // Subcolumn generators are a means for producing a variable x(i,j,k), where
            //
            //     c(i,j,k) = 1 for x(i,j,k) >  1 - cldf(i,j)
            //     c(i,j,k) = 0 for x(i,j,k) <= 1 - cldf(i,j)
            //
            // I am going to call this "cldx" to be just slightly less ambiguous
            auto cldx = real3d("cldx", ncol, nlay, ngpt);

            // Apply overlap assumption to set cldx
            if (overlap_option == 0) {  // Dummy mask, always cloudy
                memset(cldx, 1);
            } else {  // Default case, maximum-random overlap
                // Maximum-random overlap:
                // Uses essentially the algorithm described in eq (14) in Raisanen et al. 2004,
                // https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1256/qj.03.99. Also the same
                // algorithm used in RRTMG implementation of maximum-random overlap (see
                // https://github.com/AER-RC/RRTMG_SW/blob/master/src/mcica_subcol_gen_sw.f90)
                //
                // First, fill cldx with random numbers. Need to use a unique seed for each column!
                parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA(int icol) {
                    yakl::Random rand(seeds(icol));
                    for (int igpt = 1; igpt <= ngpt; igpt++) {
                        for (int ilay = 1; ilay <= nlay; ilay++) {
                            cldx(icol,ilay,igpt) = rand.genFP<Real>();
                        }
                    }
                });
                // Step down columns and apply algorithm from eq (14)
                parallel_for(SimpleBounds<2>(ngpt,ncol), YAKL_LAMBDA(int igpt, int icol) {
                    for (int ilay = 2; ilay <= nlay; ilay++) {
                        // Check cldx in level above and see if it satisfies conditions to create a cloudy subcolumn
                        if (cldx(icol,ilay-1,igpt) > 1.0 - cldf(icol,ilay-1)) {
                            // Cloudy subcolumn above, use same random number here so that clouds in these two adjacent
                            // layers are maximimally overlapped
                            cldx(icol,ilay,igpt) = cldx(icol,ilay-1,igpt);
                        } else {
                            // Cloud-less above, use new random number so that clouds are distributed
                            // randomly in this layer. Need to scale new random number to range
                            // [0, 1.0 - cldf(ilay-1)] because we have artifically changed the distribution
                            // of random numbers in this layer with the above branch of the conditional,
                            // which would otherwise inflate cloud fraction in this layer.
                            cldx(icol,ilay,igpt) = cldx(icol,ilay  ,igpt) * (1.0 - cldf(icol,ilay-1));
                        }
                    }
                });
            }

            // Use cldx array to create subcolumn mask
            parallel_for(SimpleBounds<3>(ngpt,nlay,ncol), YAKL_LAMBDA(int igpt, int ilay, int icol) {
                if (cldx(icol,ilay,igpt) > 1.0 - cldf(icol,ilay)) {
                    subcolumn_mask(icol,ilay,igpt) = 1;
                } else {
                    subcolumn_mask(icol,ilay,igpt) = 0;
                }
            });
            return subcolumn_mask;
        }


        OpticalProps2str get_subsampled_clouds(
                const int ncol, const int nlay, const int nbnd, const int ngpt,
                OpticalProps2str &cloud_optics, GasOpticsRRTMGP &kdist, real2d &cld, real2d &p_lay) {
            // Initialized subsampled optics
            OpticalProps2str subsampled_optics;
            subsampled_optics.init(kdist.get_band_lims_wavenumber(), kdist.get_band_lims_gpoint(), "subsampled_optics");
            subsampled_optics.alloc_2str(ncol, nlay);
            // Check that we do not have clouds with no optical properties; this would get corrected
            // when we assign optical props, but we want to use a "radiative cloud fraction"
            // for the subcolumn sampling too because otherwise we can get vertically-contiguous cloud
            // mask profiles with no actual cloud properties in between, which would just further overestimate
            // the vertical correlation of cloudy layers. I.e., cloudy layers might look maximally overlapped
            // even when separated by layers with no cloud properties, when in fact those layers should be
            // randomly overlapped.
            auto cldfrac_rad = real2d("cldfrac_rad", ncol, nlay);
            memset(cldfrac_rad, 0.0);  // Start with all zeros
            parallel_for(SimpleBounds<3>(nbnd,nlay,ncol), YAKL_LAMBDA (int ibnd, int ilay, int icol) {
                if (cloud_optics.tau(icol,ilay,ibnd) > 0) {
                    cldfrac_rad(icol,ilay) = cld(icol,ilay);
                } 
            });
            // Get subcolumn cloud mask; note that get_subcolumn_mask exposes overlap assumption as an option,
            // but the only currently supported options are 0 (trivial all-or-nothing cloud) or 1 (max-rand),
            // so overlap has not been exposed as an option beyond this subcolumn. In the future, we should
            // support generalized overlap as well, with parameters derived from DPSCREAM simulations with very
            // high resolution.
            int overlap = 1;
            // Get unique seeds for each column that are reproducible across different MPI rank layouts;
            // use decimal part of pressure for this, consistent with the implementation in EAM
            auto seeds = int1d("seeds", ncol);
            parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA(int icol) {
                seeds(icol) = 1e9 * (p_lay(icol,nlay) - int(p_lay(icol,nlay)));
            });
            auto cldmask = get_subcolumn_mask(ncol, nlay, ngpt, cldfrac_rad, overlap, seeds);
            // Assign optical properties to subcolumns (note this implements MCICA)
            auto gpoint_bands = kdist.get_gpoint_bands();
            parallel_for(SimpleBounds<3>(ngpt,nlay,ncol), YAKL_LAMBDA(int igpt, int ilay, int icol) {
                auto ibnd = gpoint_bands(igpt);
                if (cldmask(icol,ilay,igpt) == 1) {
                    subsampled_optics.tau(icol,ilay,igpt) = cloud_optics.tau(icol,ilay,ibnd);
                    subsampled_optics.ssa(icol,ilay,igpt) = cloud_optics.ssa(icol,ilay,ibnd);
                    subsampled_optics.g  (icol,ilay,igpt) = cloud_optics.g  (icol,ilay,ibnd);
                } else {
                    subsampled_optics.tau(icol,ilay,igpt) = 0;
                    subsampled_optics.ssa(icol,ilay,igpt) = 0;
                    subsampled_optics.g  (icol,ilay,igpt) = 0;
                } 
            });
            return subsampled_optics;
        }



        OpticalProps1scl get_subsampled_clouds(
                const int ncol, const int nlay, const int nbnd, const int ngpt,
                OpticalProps1scl &cloud_optics, GasOpticsRRTMGP &kdist, real2d &cld, real2d &p_lay) {
            // Initialized subsampled optics
            OpticalProps1scl subsampled_optics;
            subsampled_optics.init(kdist.get_band_lims_wavenumber(), kdist.get_band_lims_gpoint(), "subsampled_optics");
            subsampled_optics.alloc_1scl(ncol, nlay);
            // Check that we do not have clouds with no optical properties; this would get corrected
            // when we assign optical props, but we want to use a "radiative cloud fraction"
            // for the subcolumn sampling too because otherwise we can get vertically-contiguous cloud
            // mask profiles with no actual cloud properties in between, which would just further overestimate
            // the vertical correlation of cloudy layers. I.e., cloudy layers might look maximally overlapped
            // even when separated by layers with no cloud properties, when in fact those layers should be
            // randomly overlapped.
            auto cldfrac_rad = real2d("cldfrac_rad", ncol, nlay);
            memset(cldfrac_rad, 0.0);  // Start with all zeros
            parallel_for(SimpleBounds<3>(nbnd,nlay,ncol), YAKL_LAMBDA (int ibnd, int ilay, int icol) {
                if (cloud_optics.tau(icol,ilay,ibnd) > 0) {
                    cldfrac_rad(icol,ilay) = cld(icol,ilay);
                } 
            });
            // Get subcolumn cloud mask
            int overlap = 1;
            // Get unique seeds for each column that are reproducible across different MPI rank layouts;
            // use decimal part of pressure for this, consistent with the implementation in EAM; use different
            // seed values for longwave and shortwave
            auto seeds = int1d("seeds", ncol);
            parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA(int icol) {
                seeds(icol) = 1e9 * (p_lay(icol,nlay-1) - int(p_lay(icol,nlay-1)));
            });
            auto cldmask = get_subcolumn_mask(ncol, nlay, ngpt, cldfrac_rad, overlap, seeds);
            // Assign optical properties to subcolumns (note this implements MCICA)
            auto gpoint_bands = kdist.get_gpoint_bands();
            parallel_for(SimpleBounds<3>(ngpt,nlay,ncol), YAKL_LAMBDA(int igpt, int ilay, int icol) {
                auto ibnd = gpoint_bands(igpt);
                if (cldmask(icol,ilay,igpt) == 1) {
                    subsampled_optics.tau(icol,ilay,igpt) = cloud_optics.tau(icol,ilay,ibnd);
                } else {
                    subsampled_optics.tau(icol,ilay,igpt) = 0;
                } 
            });
            return subsampled_optics;
        }



        OpticalProps2str get_cloud_optics_sw(
                const int ncol, const int nlay,
                CloudOptics &cloud_optics, GasOpticsRRTMGP &kdist,
                real2d &lwp, real2d &iwp, real2d &rel, real2d &rei) {
 
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
                real2d &lwp, real2d &iwp, real2d &rel, real2d &rei) {

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
                real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0, 
                OpticalProps2str &aerosol, OpticalProps2str &clouds,
                FluxesByband &fluxes, FluxesBroadband &clrsky_fluxes,
                const Real tsi_scaling,
                const std::shared_ptr<spdlog::logger>& logger) {

            // Get problem sizes
            int nbnd = k_dist.get_nband();
            int ngpt = k_dist.get_ngpt();
            int ngas = gas_concs.get_num_gases();

            // Associate local pointers for fluxes
            auto &flux_up = fluxes.flux_up;
            auto &flux_dn = fluxes.flux_dn;
            auto &flux_dn_dir = fluxes.flux_dn_dir;
            auto &bnd_flux_up = fluxes.bnd_flux_up;
            auto &bnd_flux_dn = fluxes.bnd_flux_dn;
            auto &bnd_flux_dn_dir = fluxes.bnd_flux_dn_dir;
            auto &clrsky_flux_up = clrsky_fluxes.flux_up;
            auto &clrsky_flux_dn = clrsky_fluxes.flux_dn;
            auto &clrsky_flux_dn_dir = clrsky_fluxes.flux_dn_dir;

            // Reset fluxes to zero
            parallel_for(SimpleBounds<2>(nlay+1,ncol), YAKL_LAMBDA(int ilev, int icol) {
                flux_up    (icol,ilev) = 0;
                flux_dn    (icol,ilev) = 0;
                flux_dn_dir(icol,ilev) = 0;
                clrsky_flux_up    (icol,ilev) = 0;
                clrsky_flux_dn    (icol,ilev) = 0;
                clrsky_flux_dn_dir(icol,ilev) = 0;
            });
            parallel_for(SimpleBounds<3>(nbnd,nlay+1,ncol), YAKL_LAMBDA(int ibnd, int ilev, int icol) {
                bnd_flux_up    (icol,ilev,ibnd) = 0;
                bnd_flux_dn    (icol,ilev,ibnd) = 0;
                bnd_flux_dn_dir(icol,ilev,ibnd) = 0;
            });
 
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
                // No daytime columns in this chunk, skip the rest of this routine
                return;
            }

            // Subset mu0
            auto mu0_day = real1d("mu0_day", nday);
            parallel_for(SimpleBounds<1>(nday), YAKL_LAMBDA(int iday) {
                mu0_day(iday) = mu0(dayIndices(iday));
            });

            // subset state variables
            auto p_lay_day = real2d("p_lay_day", nday, nlay);
            auto t_lay_day = real2d("t_lay_day", nday, nlay);
            parallel_for(SimpleBounds<2>(nlay,nday), YAKL_LAMBDA(int ilay, int iday) {
                p_lay_day(iday,ilay) = p_lay(dayIndices(iday),ilay);
                t_lay_day(iday,ilay) = t_lay(dayIndices(iday),ilay);
            });
            auto p_lev_day = real2d("p_lev_day", nday, nlay+1);
            auto t_lev_day = real2d("t_lev_day", nday, nlay+1);
            parallel_for(SimpleBounds<2>(nlay+1,nday), YAKL_LAMBDA(int ilev, int iday) {
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
                parallel_for(SimpleBounds<2>(nlay,nday), YAKL_LAMBDA(int ilay, int iday) {
                    vmr_day(iday,ilay) = vmr(dayIndices(iday),ilay);
                });
                gas_concs_day.set_vmr(gas_names(igas), vmr_day);
            }

            // Subset aerosol optics
            OpticalProps2str aerosol_day;
            aerosol_day.init(k_dist.get_band_lims_wavenumber());
            aerosol_day.alloc_2str(nday, nlay);
            parallel_for(SimpleBounds<3>(nbnd,nlay,nday), YAKL_LAMBDA(int ibnd, int ilay, int iday) {
                aerosol_day.tau(iday,ilay,ibnd) = aerosol.tau(dayIndices(iday),ilay,ibnd);
                aerosol_day.ssa(iday,ilay,ibnd) = aerosol.ssa(dayIndices(iday),ilay,ibnd);
                aerosol_day.g  (iday,ilay,ibnd) = aerosol.g  (dayIndices(iday),ilay,ibnd);
            });

            // Subset cloud optics
            // TODO: nbnd -> ngpt once we pass sub-sampled cloud state
            OpticalProps2str clouds_day;
            clouds_day.init(k_dist.get_band_lims_wavenumber(), k_dist.get_band_lims_gpoint());
            clouds_day.alloc_2str(nday, nlay);
            parallel_for(SimpleBounds<3>(ngpt,nlay,nday), YAKL_LAMBDA(int igpt, int ilay, int iday) {
                clouds_day.tau(iday,ilay,igpt) = clouds.tau(dayIndices(iday),ilay,igpt);
                clouds_day.ssa(iday,ilay,igpt) = clouds.ssa(dayIndices(iday),ilay,igpt);
                clouds_day.g  (iday,ilay,igpt) = clouds.g  (dayIndices(iday),ilay,igpt);
            });

            // RRTMGP assumes surface albedos have a screwy dimension ordering
            // for some strange reason, so we need to transpose these; also do
            // daytime subsetting in the same kernel
            real2d sfc_alb_dir_T("sfc_alb_dir", nbnd, nday);
            real2d sfc_alb_dif_T("sfc_alb_dif", nbnd, nday);
            parallel_for(SimpleBounds<2>(nbnd,nday), YAKL_LAMBDA(int ibnd, int icol) {
                sfc_alb_dir_T(ibnd,icol) = sfc_alb_dir(dayIndices(icol),ibnd);
                sfc_alb_dif_T(ibnd,icol) = sfc_alb_dif(dayIndices(icol),ibnd);
            });

            // Temporaries we need for daytime-only fluxes
            auto flux_up_day = real2d("flux_up_day", nday, nlay+1);
            auto flux_dn_day = real2d("flux_dn_day", nday, nlay+1);
            auto flux_dn_dir_day = real2d("flux_dn_dir_day", nday, nlay+1);
            auto bnd_flux_up_day = real3d("bnd_flux_up_day", nday, nlay+1, nbnd);
            auto bnd_flux_dn_day = real3d("bnd_flux_dn_day", nday, nlay+1, nbnd);
            auto bnd_flux_dn_dir_day = real3d("bnd_flux_dn_dir_day", nday, nlay+1, nbnd);
            FluxesByband fluxes_day;
            fluxes_day.flux_up         = flux_up_day;
            fluxes_day.flux_dn         = flux_dn_day;
            fluxes_day.flux_dn_dir     = flux_dn_dir_day;
            fluxes_day.bnd_flux_up     = bnd_flux_up_day;
            fluxes_day.bnd_flux_dn     = bnd_flux_dn_day;
            fluxes_day.bnd_flux_dn_dir = bnd_flux_dn_dir_day;

            // Allocate space for optical properties
            OpticalProps2str optics;
            optics.alloc_2str(nday, nlay, k_dist);

            // Limit temperatures for gas optics look-up tables
            auto t_lay_limited = real2d("t_lay_limited", nday, nlay);
            limit_to_bounds(t_lay_day, k_dist_sw.get_temp_min(), k_dist_sw.get_temp_max(), t_lay_limited);

            // Do gas optics
            real2d toa_flux("toa_flux", nday, ngpt);
            auto p_lay_host = p_lay.createHostCopy();
            bool top_at_1 = p_lay_host(1, 1) < p_lay_host(1, nlay);

            k_dist.gas_optics(nday, nlay, top_at_1, p_lay_day, p_lev_day, t_lay_limited, gas_concs_day, optics, toa_flux);

#ifdef SCREAM_RRTMGP_DEBUG
            // Check gas optics
            check_range(optics.tau,  0, std::numeric_limits<Real>::max(), "rrtmgp_sw:optics.tau");
            check_range(optics.ssa,  0,                                1, "rrtmgp_sw:optics.ssa"); //, "optics.ssa");
            check_range(optics.g  , -1,                                1, "rrtmgp_sw:optics.g  "); //, "optics.g"  );
#endif

            // Apply tsi_scaling
            parallel_for(SimpleBounds<2>(ngpt,nday), YAKL_LAMBDA(int igpt, int iday) {
                toa_flux(iday,igpt) = tsi_scaling * toa_flux(iday,igpt);
            });

            // Combine gas and aerosol optics
            aerosol_day.delta_scale();
            aerosol_day.increment(optics);

            // Compute clearsky (gas + aerosol) fluxes on daytime columns
            rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);
            // Expand daytime fluxes to all columns
            parallel_for(SimpleBounds<2>(nlay+1,nday), YAKL_LAMBDA(int ilev, int iday) {
                int icol = dayIndices(iday);
                clrsky_flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
                clrsky_flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
                clrsky_flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
            });

            // Now merge in cloud optics and do allsky calculations

            // Combine gas and cloud optics
            clouds_day.delta_scale();
            clouds_day.increment(optics);
            // Compute fluxes on daytime columns
            rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);
            // Expand daytime fluxes to all columns
            parallel_for(SimpleBounds<2>(nlay+1,nday), YAKL_LAMBDA(int ilev, int iday) {
                int icol = dayIndices(iday);
                flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
                flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
                flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
            });
            parallel_for(SimpleBounds<3>(nbnd,nlay+1,nday), YAKL_LAMBDA(int ibnd, int ilev, int iday) {
                int icol = dayIndices(iday);
                bnd_flux_up    (icol,ilev,ibnd) = bnd_flux_up_day    (iday,ilev,ibnd);
                bnd_flux_dn    (icol,ilev,ibnd) = bnd_flux_dn_day    (iday,ilev,ibnd);
                bnd_flux_dn_dir(icol,ilev,ibnd) = bnd_flux_dn_dir_day(iday,ilev,ibnd);
            });
        }

        void rrtmgp_lw(
                const int ncol, const int nlay,
                GasOpticsRRTMGP &k_dist,
                real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
                GasConcs &gas_concs,
                OpticalProps1scl &aerosol,
                OpticalProps1scl &clouds,
                FluxesByband &fluxes, FluxesBroadband &clrsky_fluxes) {

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
            parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA(int icol) {
                t_sfc(icol) = t_lev(icol, merge(nlay+1, 1, top_at_1));
            });
            memset(emis_sfc , 0.98_wp);

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

            // Limit temperatures for gas optics look-up tables
            auto t_lay_limited = real2d("t_lay_limited", ncol, nlay);
            auto t_lev_limited = real2d("t_lev_limited", ncol, nlay+1);
            limit_to_bounds(t_lay, k_dist_lw.get_temp_min(), k_dist_lw.get_temp_max(), t_lay_limited);
            limit_to_bounds(t_lev, k_dist_lw.get_temp_min(), k_dist_lw.get_temp_max(), t_lev_limited);

            // Do gas optics
            k_dist.gas_optics(ncol, nlay, top_at_1, p_lay, p_lev, t_lay_limited, t_sfc, gas_concs, optics, lw_sources, real2d(), t_lev_limited);

#ifdef SCREAM_RRTMGP_DEBUG
            // Check gas optics
            check_range(optics.tau,  0, std::numeric_limits<Real>::max(), "rrtmgp_lw:optics.tau");
#endif

            // Combine gas and aerosol optics
            aerosol.increment(optics);

            // Compute clear-sky fluxes before we add in clouds
            rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics, top_at_1, lw_sources, emis_sfc, clrsky_fluxes);

            // Combine gas and cloud optics
            clouds.increment(optics);

            // Compute allsky fluxes
            rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics, top_at_1, lw_sources, emis_sfc, fluxes);

        }

        void compute_cloud_area(
                int ncol, int nlay, int ngpt, const Real pmin, const Real pmax, 
                const real2d& pmid, const real3d& cld_tau_gpt, real1d& cld_area) {
            // Subcolumn binary cld mask; if any layers with pressure between pmin and pmax are cloudy
            // then 2d subcol mask is 1, otherwise it is 0
            auto subcol_mask = real2d("subcol_mask", ncol, ngpt);
            memset(subcol_mask, 0);
            yakl::fortran::parallel_for(SimpleBounds<3>(ngpt, nlay, ncol), YAKL_LAMBDA(int igpt, int ilay, int icol) {
                // NOTE: using plev would need to assume level ordering (top to bottom or bottom to top), but
                // using play/pmid does not
                if (cld_tau_gpt(icol,ilay,igpt) > 0 && pmid(icol,ilay) >= pmin && pmid(icol,ilay) < pmax) {
                    subcol_mask(icol,igpt) = 1;
                }
            });
            // Compute average over subcols to get cloud area
            auto ngpt_inv = 1.0 / ngpt;
            memset(cld_area, 0);
            yakl::fortran::parallel_for(SimpleBounds<1>(ncol), YAKL_LAMBDA(int icol) {
                // This loop needs to be serial because of the atomic reduction
                for (int igpt = 1; igpt <= ngpt; ++igpt) {
                    cld_area(icol) += subcol_mask(icol,igpt) * ngpt_inv;
                }
            });
        }

    }  // namespace rrtmgp
}  // namespace scream
