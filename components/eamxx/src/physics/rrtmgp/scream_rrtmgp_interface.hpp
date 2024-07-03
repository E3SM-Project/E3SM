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

void init_kls ();
void finalize_kls();

namespace rrtmgp {

/*
 * Objects containing k-distribution information need to be initialized
 * once and then persist throughout the life of the program, so we
 * declare them here within the rrtmgp namespace.
 */
#ifdef RRTMGP_ENABLE_YAKL
extern GasOpticsRRTMGP k_dist_sw;
extern GasOpticsRRTMGP k_dist_lw;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
extern GasOpticsRRTMGPK k_dist_sw_k;
extern GasOpticsRRTMGPK k_dist_lw_k;
#endif

/*
 * Objects containing cloud optical property look-up table information.
 * We want to initialize these once and use throughout the life of the
 * program, so declare here and read data in during rrtmgp_initialize().
 */
#ifdef RRTMGP_ENABLE_YAKL
extern CloudOptics cloud_optics_sw;
extern CloudOptics cloud_optics_lw;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
extern CloudOpticsK cloud_optics_sw_k;
extern CloudOpticsK cloud_optics_lw_k;
#endif

/*
 * Flag to indicate whether or not we have initialized RRTMGP
 */
#ifdef RRTMGP_ENABLE_YAKL
extern bool initialized;
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
extern bool initialized_k;
#endif

/*
 * Initialize data for RRTMGP driver
 */
#ifdef RRTMGP_ENABLE_YAKL
extern void rrtmgp_initialize(
  GasConcs &gas_concs,
  const std::string& coefficients_file_sw, const std::string& coefficients_file_lw,
  const std::string& cloud_optics_file_sw, const std::string& cloud_optics_file_lw,
  const std::shared_ptr<spdlog::logger>& logger);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
extern void rrtmgp_initialize(
  GasConcsK &gas_concs,
  const std::string& coefficients_file_sw, const std::string& coefficients_file_lw,
  const std::string& cloud_optics_file_sw, const std::string& cloud_optics_file_lw,
  const std::shared_ptr<spdlog::logger>& logger);
#endif

/*
 * Compute band-by-band surface albedos from broadband albedos.
 */
#ifdef RRTMGP_ENABLE_YAKL
extern void compute_band_by_band_surface_albedos(
  const int ncol, const int nswbands,
  real1d &sfc_alb_dir_vis, real1d &sfc_alb_dir_nir,
  real1d &sfc_alb_dif_vis, real1d &sfc_alb_dif_nir,
  real2d &sfc_alb_dir,     real2d &sfc_alb_dif);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
extern void compute_band_by_band_surface_albedos(
  const int ncol, const int nswbands,
  real1dk &sfc_alb_dir_vis, real1dk &sfc_alb_dir_nir,
  real1dk &sfc_alb_dif_vis, real1dk &sfc_alb_dif_nir,
  real2dk &sfc_alb_dir,     real2dk &sfc_alb_dif);
#endif

/*
 * Compute broadband visible/UV and near-infrared surface fluxes.
 */
#ifdef RRTMGP_ENABLE_YAKL
extern void compute_broadband_surface_fluxes(
  const int ncol, const int ktop, const int nswbands,
  real3d &sw_bnd_flux_dir , real3d &sw_bnd_flux_dif ,
  real1d &sfc_flux_dir_vis, real1d &sfc_flux_dir_nir,
  real1d &sfc_flux_dif_vis, real1d &sfc_flux_dif_nir);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
extern void compute_broadband_surface_fluxes(
  const int ncol, const int ktop, const int nswbands,
  real3dk &sw_bnd_flux_dir , real3dk &sw_bnd_flux_dif ,
  real1dk &sfc_flux_dir_vis, real1dk &sfc_flux_dir_nir,
  real1dk &sfc_flux_dif_vis, real1dk &sfc_flux_dif_nir);
#endif

/*
 * Main driver code to run RRTMGP.
 * The input logger is in charge of outputing info to
 * screen and/or to file (or neither), depending on how it was set up.
 */
#ifdef RRTMGP_ENABLE_YAKL
extern void rrtmgp_main(
  const int ncol, const int nlay,
  real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
  GasConcs &gas_concs,
  real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0,
  real2d &lwp, real2d &iwp, real2d &rel, real2d &rei, real2d &cldfrac,
  real3d &aer_tau_sw, real3d &aer_ssa_sw, real3d &aer_asm_sw, real3d &aer_tau_lw,
  real3d &cld_tau_sw_bnd, real3d &cld_tau_lw_bnd,
  real3d &cld_tau_sw_gpt, real3d &cld_tau_lw_gpt,
  real2d &sw_flux_up, real2d &sw_flux_dn, real2d &sw_flux_dn_dir,
  real2d &lw_flux_up, real2d &lw_flux_dn,
  real2d &sw_clnclrsky_flux_up, real2d &sw_clnclrsky_flux_dn, real2d &sw_clnclrsky_flux_dn_dir,
  real2d &sw_clrsky_flux_up, real2d &sw_clrsky_flux_dn, real2d &sw_clrsky_flux_dn_dir,
  real2d &sw_clnsky_flux_up, real2d &sw_clnsky_flux_dn, real2d &sw_clnsky_flux_dn_dir,
  real2d &lw_clnclrsky_flux_up, real2d &lw_clnclrsky_flux_dn,
  real2d &lw_clrsky_flux_up, real2d &lw_clrsky_flux_dn,
  real2d &lw_clnsky_flux_up, real2d &lw_clnsky_flux_dn,
  real3d &sw_bnd_flux_up, real3d &sw_bnd_flux_dn, real3d &sw_bnd_flux_dn_dir,
  real3d &lw_bnd_flux_up, real3d &lw_bnd_flux_dn,
  const Real tsi_scaling,
  const std::shared_ptr<spdlog::logger>& logger,
  const bool extra_clnclrsky_diag = false, const bool extra_clnsky_diag = false);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
extern void rrtmgp_main(
  const int ncol, const int nlay,
  real2dk &p_lay, real2dk &t_lay, real2dk &p_lev, real2dk &t_lev,
  GasConcsK &gas_concs,
  real2dk &sfc_alb_dir, real2dk &sfc_alb_dif, real1dk &mu0,
  real2dk &lwp, real2dk &iwp, real2dk &rel, real2dk &rei, real2dk &cldfrac,
  real3dk &aer_tau_sw, real3dk &aer_ssa_sw, real3dk &aer_asm_sw, real3dk &aer_tau_lw,
  real3dk &cld_tau_sw_bnd, real3dk &cld_tau_lw_bnd,
  real3dk &cld_tau_sw_gpt, real3dk &cld_tau_lw_gpt,
  real2dk &sw_flux_up, real2dk &sw_flux_dn, real2dk &sw_flux_dn_dir,
  real2dk &lw_flux_up, real2dk &lw_flux_dn,
  real2dk &sw_clnclrsky_flux_up, real2dk &sw_clnclrsky_flux_dn, real2dk &sw_clnclrsky_flux_dn_dir,
  real2dk &sw_clrsky_flux_up, real2dk &sw_clrsky_flux_dn, real2dk &sw_clrsky_flux_dn_dir,
  real2dk &sw_clnsky_flux_up, real2dk &sw_clnsky_flux_dn, real2dk &sw_clnsky_flux_dn_dir,
  real2dk &lw_clnclrsky_flux_up, real2dk &lw_clnclrsky_flux_dn,
  real2dk &lw_clrsky_flux_up, real2dk &lw_clrsky_flux_dn,
  real2dk &lw_clnsky_flux_up, real2dk &lw_clnsky_flux_dn,
  real3dk &sw_bnd_flux_up, real3dk &sw_bnd_flux_dn, real3dk &sw_bnd_flux_dn_dir,
  real3dk &lw_bnd_flux_up, real3dk &lw_bnd_flux_dn,
  const Real tsi_scaling,
  const std::shared_ptr<spdlog::logger>& logger,
  const bool extra_clnclrsky_diag = false, const bool extra_clnsky_diag = false);
#endif

/*
 * Perform any clean-up tasks
 */
extern void rrtmgp_finalize();

/*
 * Shortwave driver (called by rrtmgp_main)
 */
#ifdef RRTMGP_ENABLE_YAKL
extern void rrtmgp_sw(
  const int ncol, const int nlay,
  GasOpticsRRTMGP &k_dist,
  real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
  GasConcs &gas_concs,
  real2d &sfc_alb_dir, real2d &sfc_alb_dif, real1d &mu0,
  OpticalProps2str &aerosol, OpticalProps2str &clouds,
  FluxesByband &fluxes, FluxesBroadband &clnclrsky_fluxes, FluxesBroadband &clrsky_fluxes, FluxesBroadband &clnsky_fluxes,
  const Real tsi_scaling,
  const std::shared_ptr<spdlog::logger>& logger,
  const bool extra_clnclrsky_diag, const bool extra_clnsky_diag);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
extern void rrtmgp_sw(
  const int ncol, const int nlay,
  GasOpticsRRTMGPK &k_dist,
  real2dk &p_lay, real2dk &t_lay, real2dk &p_lev, real2dk &t_lev,
  GasConcsK &gas_concs,
  real2dk &sfc_alb_dir, real2dk &sfc_alb_dif, real1dk &mu0,
  OpticalProps2strK &aerosol, OpticalProps2strK &clouds,
  FluxesBybandK &fluxes, FluxesBroadbandK &clnclrsky_fluxes, FluxesBroadbandK &clrsky_fluxes, FluxesBroadbandK &clnsky_fluxes,
  const Real tsi_scaling,
  const std::shared_ptr<spdlog::logger>& logger,
  const bool extra_clnclrsky_diag, const bool extra_clnsky_diag);
#endif

/*
 * Longwave driver (called by rrtmgp_main)
 */
#ifdef RRTMGP_ENABLE_YAKL
extern void rrtmgp_lw(
  const int ncol, const int nlay,
  GasOpticsRRTMGP &k_dist,
  real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
  GasConcs &gas_concs,
  OpticalProps1scl &aerosol, OpticalProps1scl &clouds,
  FluxesByband &fluxes, FluxesBroadband &clnclrsky_fluxes, FluxesBroadband &clrsky_fluxes, FluxesBroadband &clnsky_fluxes,
  const bool extra_clnclrsky_diag, const bool extra_clnsky_diag);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
extern void rrtmgp_lw(
  const int ncol, const int nlay,
  GasOpticsRRTMGPK &k_dist,
  real2dk &p_lay, real2dk &t_lay, real2dk &p_lev, real2dk &t_lev,
  GasConcsK &gas_concs,
  OpticalProps1sclK &aerosol, OpticalProps1sclK &clouds,
  FluxesBybandK &fluxes, FluxesBroadbandK &clnclrsky_fluxes, FluxesBroadbandK &clrsky_fluxes, FluxesBroadbandK &clnsky_fluxes,
  const bool extra_clnclrsky_diag, const bool extra_clnsky_diag);
#endif

/*
 * Return a subcolumn mask consistent with a specified overlap assumption
 */
#ifdef RRTMGP_ENABLE_YAKL
int3d get_subcolumn_mask(const int ncol, const int nlay, const int ngpt, real2d &cldf, const int overlap_option, int1d &seeds);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
int3dk get_subcolumn_mask(const int ncol, const int nlay, const int ngpt, real2dk &cldf, const int overlap_option, int1dk &seeds);
#endif

/*
 * Compute cloud area from 3d subcol cloud property
 */
#ifdef RRTMGP_ENABLE_YAKL
void compute_cloud_area(
  int ncol, int nlay, int ngpt, Real pmin, Real pmax,
  const real2d& pmid, const real3d& cld_tau_gpt, real1d& cld_area);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
void compute_cloud_area(
  int ncol, int nlay, int ngpt, Real pmin, Real pmax,
  const real2dk& pmid, const real3dk& cld_tau_gpt, real1dk& cld_area);
#endif

/*
 * Return select cloud-top diagnostics following AeroCom recommendation
 */
#ifdef RRTMGP_ENABLE_YAKL
void compute_aerocom_cloudtop(
  int ncol, int nlay, const real2d &tmid, const real2d &pmid,
  const real2d &p_del, const real2d &z_del, const real2d &qc,
  const real2d &qi, const real2d &rel, const real2d &rei,
  const real2d &cldfrac_tot, const real2d &nc,
  real1d &T_mid_at_cldtop, real1d &p_mid_at_cldtop,
  real1d &cldfrac_ice_at_cldtop, real1d &cldfrac_liq_at_cldtop,
  real1d &cldfrac_tot_at_cldtop, real1d &cdnc_at_cldtop,
  real1d &eff_radius_qc_at_cldtop, real1d &eff_radius_qi_at_cldtop);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
void compute_aerocom_cloudtop(
  int ncol, int nlay, const real2dk &tmid, const real2dk &pmid,
  const real2dk &p_del, const real2dk &z_del, const real2dk &qc,
  const real2dk &qi, const real2dk &rel, const real2dk &rei,
  const real2dk &cldfrac_tot, const real2dk &nc,
  real1dk &T_mid_at_cldtop, real1dk &p_mid_at_cldtop,
  real1dk &cldfrac_ice_at_cldtop, real1dk &cldfrac_liq_at_cldtop,
  real1dk &cldfrac_tot_at_cldtop, real1dk &cdnc_at_cldtop,
  real1dk &eff_radius_qc_at_cldtop, real1dk &eff_radius_qi_at_cldtop);
#endif

/*
 * Provide a function to convert cloud (water and ice) mixing ratios to layer mass per unit area
 * (what E3SM refers to as "in-cloud water paths", a terminology we shun here to avoid confusion
 * with the standard practice of using "water path" to refer to the total column-integrated
 * quantities).
 */
#ifdef RRTMGP_ENABLE_YAKL
template<class T, int myMem, int myStyle>
void mixing_ratio_to_cloud_mass(
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
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
template<class View1, class View2, class View3, class View4>
void mixing_ratio_to_cloud_mass(
  View1 const& mixing_ratio,
  View2 const& cloud_fraction,
  View3 const& dp,
  View4 const& cloud_mass)
{
  int ncol = mixing_ratio.extent(0);
  int nlay = mixing_ratio.extent(1);
  using physconst = scream::physics::Constants<Real>;
  Kokkos::parallel_for(conv::get_mdrp<2>({nlay, ncol}), KOKKOS_LAMBDA(int ilay, int icol) {
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
#endif

/*
 * Routine to limit a quantity to set bounds. Used to make sure
 * effective radii are within the bounds of the cloud optical
 * property look-up tables, but could be used to limit other
 * fields as well.
 */
#ifdef RRTMGP_ENABLE_YAKL
template<class S, class T>
void limit_to_bounds(S const &arr_in, T const lower, T const upper, S &arr_out) {
  yakl::c::parallel_for(arr_in.totElems(), YAKL_LAMBDA(int i) {
    arr_out.data()[i] = std::min(std::max(arr_in.data()[i], lower), upper);
  });
}
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
template<class S, class T>
void limit_to_bounds_k(S const &arr_in, T const lower, T const upper, S &arr_out) {
  Kokkos::parallel_for(arr_in.size(), KOKKOS_LAMBDA(int i) {
    arr_out.data()[i] = std::min(std::max(arr_in.data()[i], lower), upper);
  });
}
#endif

#ifdef RRTMGP_ENABLE_YAKL
int get_wavelength_index(OpticalProps &kdist, double wavelength);
int get_wavelength_index_sw(double wavelength);
int get_wavelength_index_lw(double wavelength);
#endif
#ifdef RRTMGP_ENABLE_KOKKOS
int get_wavelength_index(OpticalPropsK &kdist, double wavelength);
int get_wavelength_index_sw_k(double wavelength);
int get_wavelength_index_lw_k(double wavelength);
#endif

} // namespace rrtmgp
} // namespace scream

#endif  // SCREAM_RRTMGP_INTERFACE_HPP
