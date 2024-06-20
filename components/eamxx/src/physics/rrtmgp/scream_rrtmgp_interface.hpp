#ifndef SCREAM_RRTMGP_INTERFACE_HPP
#define SCREAM_RRTMGP_INTERFACE_HPP

#include "cpp/rrtmgp/mo_gas_optics_rrtmgp.h"
#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"
#include "cpp/extensions/fluxes_byband/mo_fluxes_byband.h"
#include "cpp/examples/mo_load_coefficients.h"
#include "cpp/rrtmgp_const.h"
#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "cpp/rrtmgp/mo_gas_optics_rrtmgp.h"
#include "cpp/extensions/cloud_optics/mo_cloud_optics.h"
#include "cpp/rte/mo_rte_sw.h"
#include "cpp/rte/mo_rte_lw.h"
#include "examples/all-sky/mo_load_cloud_coefficients.h"

#include "rrtmgp_utils.hpp"

#include "physics/share/physics_constants.hpp"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/logging/ekat_logger.hpp"
#include "ekat/util/ekat_math_utils.hpp"

#ifdef RRTMGP_ENABLE_KOKKOS
#include "Kokkos_Random.hpp"
#endif

namespace scream {

void init_kls ();
void finalize_kls();

namespace rrtmgp {

#ifdef RRTMGP_ENABLE_YAKL
extern GasOpticsRRTMGP k_dist_sw;
extern GasOpticsRRTMGP k_dist_lw;

extern CloudOptics cloud_optics_sw;
extern CloudOptics cloud_optics_lw;

extern bool initialized;

void rrtmgp_initialize(
  GasConcs &gas_concs,
  const std::string& coefficients_file_sw, const std::string& coefficients_file_lw,
  const std::string& cloud_optics_file_sw, const std::string& cloud_optics_file_lw,
  const std::shared_ptr<spdlog::logger>& logger);

void compute_band_by_band_surface_albedos(
  const int ncol, const int nswbands,
  real1d &sfc_alb_dir_vis, real1d &sfc_alb_dir_nir,
  real1d &sfc_alb_dif_vis, real1d &sfc_alb_dif_nir,
  real2d &sfc_alb_dir,     real2d &sfc_alb_dif);

void compute_broadband_surface_fluxes(
  const int ncol, const int ktop, const int nswbands,
  real3d &sw_bnd_flux_dir , real3d &sw_bnd_flux_dif ,
  real1d &sfc_flux_dir_vis, real1d &sfc_flux_dir_nir,
  real1d &sfc_flux_dif_vis, real1d &sfc_flux_dif_nir);

void rrtmgp_main(
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

void rrtmgp_finalize();

void rrtmgp_sw(
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

void rrtmgp_lw(
  const int ncol, const int nlay,
  GasOpticsRRTMGP &k_dist,
  real2d &p_lay, real2d &t_lay, real2d &p_lev, real2d &t_lev,
  GasConcs &gas_concs,
  OpticalProps1scl &aerosol, OpticalProps1scl &clouds,
  FluxesByband &fluxes, FluxesBroadband &clnclrsky_fluxes, FluxesBroadband &clrsky_fluxes, FluxesBroadband &clnsky_fluxes,
  const bool extra_clnclrsky_diag, const bool extra_clnsky_diag);

int3d get_subcolumn_mask(const int ncol, const int nlay, const int ngpt, real2d &cldf, const int overlap_option, int1d &seeds);

void compute_cloud_area(
  int ncol, int nlay, int ngpt, Real pmin, Real pmax,
  const real2d& pmid, const real3d& cld_tau_gpt, real1d& cld_area);

void compute_aerocom_cloudtop(
  int ncol, int nlay, const real2d &tmid, const real2d &pmid,
  const real2d &p_del, const real2d &z_del, const real2d &qc,
  const real2d &qi, const real2d &rel, const real2d &rei,
  const real2d &cldfrac_tot, const real2d &nc,
  real1d &T_mid_at_cldtop, real1d &p_mid_at_cldtop,
  real1d &cldfrac_ice_at_cldtop, real1d &cldfrac_liq_at_cldtop,
  real1d &cldfrac_tot_at_cldtop, real1d &cdnc_at_cldtop,
  real1d &eff_radius_qc_at_cldtop, real1d &eff_radius_qi_at_cldtop);

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

template<class S, class T>
void limit_to_bounds(S const &arr_in, T const lower, T const upper, S &arr_out) {
  yakl::c::parallel_for(arr_in.totElems(), YAKL_LAMBDA(int i) {
    arr_out.data()[i] = std::min(std::max(arr_in.data()[i], lower), upper);
  });
}

int get_wavelength_index(OpticalProps &kdist, double wavelength);
int get_wavelength_index_sw(double wavelength);
int get_wavelength_index_lw(double wavelength);
#endif // RRTMGP_ENABLE_YAKL

// New interface for Kokkos and flexible types
#ifdef RRTMGP_ENABLE_KOKKOS
template <typename RealT=Real, typename LayoutT=Kokkos::LayoutRight, typename DeviceT=DefaultDevice>
struct rrtmgp_interface {

using MDRP = typename conv::MDRP<LayoutT>;

template <typename T>
using view_t = Kokkos::View<T, LayoutT, DeviceT>;

template <typename T>
using oview_t = Kokkos::Experimental::OffsetView<T, LayoutT, DeviceT>;

template <typename T>
using hview_t = Kokkos::View<T, LayoutT, HostDevice>;

using pool_t = conv::MemPoolSingleton<RealT, DeviceT>;

/*
 * Objects containing k-distribution information need to be initialized
 * once and then persist throughout the life of the program, so we
 * declare them here within the rrtmgp namespace.
 */
static inline GasOpticsRRTMGPK<RealT, LayoutT, DeviceT> k_dist_sw_k;
static inline GasOpticsRRTMGPK<RealT, LayoutT, DeviceT> k_dist_lw_k;

/*
 * Objects containing cloud optical property look-up table information.
 * We want to initialize these once and use throughout the life of the
 * program, so declare here and read data in during rrtmgp_initialize().
 */
static inline CloudOpticsK<RealT, LayoutT, DeviceT> cloud_optics_sw_k;
static inline CloudOpticsK<RealT, LayoutT, DeviceT> cloud_optics_lw_k;

/*
 * Flag to indicate whether or not we have initialized RRTMGP
 */
static inline bool initialized_k = false;

/*
 * Initialize data for RRTMGP driver
 */
static void rrtmgp_initialize(
  GasConcsK<RealT, LayoutT, DeviceT> &gas_concs,
  const std::string& coefficients_file_sw, const std::string& coefficients_file_lw,
  const std::string& cloud_optics_file_sw, const std::string& cloud_optics_file_lw,
  const std::shared_ptr<spdlog::logger>& logger)
{
  // If we've already initialized, just exit
  if (initialized_k) {
    if (logger)
      logger->info("RRTMGP is already initialized; skipping\n");
    return;
  }

  // Initialize Kokkos
  if (!Kokkos::is_initialized()) {  Kokkos::initialize(); }

  // Load and initialize absorption coefficient data
  load_and_init(k_dist_sw_k, coefficients_file_sw, gas_concs);
  load_and_init(k_dist_lw_k, coefficients_file_lw, gas_concs);

  // Load and initialize cloud optical property look-up table information
  load_cld_lutcoeff(cloud_optics_sw_k, cloud_optics_file_sw);
  load_cld_lutcoeff(cloud_optics_lw_k, cloud_optics_file_lw);

  // initialize kokkos rrtmgp pool allocator
  const size_t base_ref = 18000;
  const size_t ncol = gas_concs.ncol;
  const size_t nlay = gas_concs.nlay;
  const size_t nlev = SCREAM_NUM_VERTICAL_LEV;
  const size_t my_size_ref = ncol * nlay * nlev;
  pool_t::init(2e6 * (float(my_size_ref) / base_ref));

  // We are now initialized!
  initialized_k = true;
}

/*
 * Compute band-by-band surface albedos from broadband albedos.
 */
template <typename SfcAlbDirVisT, typename SfcAlbDirNirT, typename SfcAlbDifVisT,
          typename SfcAlbDifNirT, typename SfcAlbDirT, typename SfcAlbDifT>
static void compute_band_by_band_surface_albedos(
  const int ncol, const int nswbands,
  SfcAlbDirVisT &sfc_alb_dir_vis, SfcAlbDirNirT &sfc_alb_dir_nir,
  SfcAlbDifVisT &sfc_alb_dif_vis, SfcAlbDifNirT &sfc_alb_dif_nir,
  SfcAlbDirT    &sfc_alb_dir,     SfcAlbDifT    &sfc_alb_dif)
{
  EKAT_ASSERT_MSG(initialized_k, "Error! rrtmgp_initialize must be called before GasOpticsRRTMGP object can be used.");
  auto wavenumber_limits = k_dist_sw_k.get_band_lims_wavenumber();

  EKAT_ASSERT_MSG(wavenumber_limits.extent(0) == 2,
                  "Error! 1st dimension for wavenumber_limits should be 2. It's " << wavenumber_limits.extent(0));
  EKAT_ASSERT_MSG(wavenumber_limits.extent(1) == static_cast<size_t>(nswbands),
                  "Error! 2nd dimension for wavenumber_limits should be " + std::to_string(nswbands) + " (nswbands).");

  // Loop over bands, and determine for each band whether it is broadly in the
  // visible or infrared part of the spectrum (visible or "not visible")
  Kokkos::parallel_for(MDRP::template get<2>({nswbands, ncol}), KOKKOS_LAMBDA(const int ibnd, const int icol) {

    // Threshold between visible and infrared is 0.7 micron, or 14286 cm^-1.
    const RealT visible_wavenumber_threshold = 14286;

    // Wavenumber is in the visible if it is above the visible wavenumber
    // threshold, and in the infrared if it is below the threshold
    const bool is_visible_wave1 = (wavenumber_limits(0, ibnd) > visible_wavenumber_threshold ? true : false);
    const bool is_visible_wave2 = (wavenumber_limits(1, ibnd) > visible_wavenumber_threshold ? true : false);

    if (is_visible_wave1 && is_visible_wave2) {
      // Entire band is in the visible
      sfc_alb_dir(icol,ibnd) = sfc_alb_dir_vis(icol);
      sfc_alb_dif(icol,ibnd) = sfc_alb_dif_vis(icol);
    }
    else if (!is_visible_wave1 && !is_visible_wave2) {
      // Entire band is in the longwave (near-infrared)
      sfc_alb_dir(icol,ibnd) = sfc_alb_dir_nir(icol);
      sfc_alb_dif(icol,ibnd) = sfc_alb_dif_nir(icol);
    }
    else {
      // Band straddles the visible to near-infrared transition, so we take
      // the albedo to be the average of the visible and near-infrared
      // broadband albedos
      sfc_alb_dir(icol,ibnd) = 0.5*(sfc_alb_dir_vis(icol) + sfc_alb_dir_nir(icol));
      sfc_alb_dif(icol,ibnd) = 0.5*(sfc_alb_dif_vis(icol) + sfc_alb_dif_nir(icol));
    }
  });
}

/*
 * Compute broadband visible/UV and near-infrared surface fluxes.
 */
template <typename SwBndFluxDirT, typename SwBndFluxDifT, typename SfcFluxDirVisT,
          typename SfcFluxDirNirT, typename SfcFluxDifVisT, typename SfcFluxDifNirT>
static void compute_broadband_surface_fluxes(
  const int ncol, const int ktop, const int nswbands,
  SwBndFluxDirT &sw_bnd_flux_dir , SwBndFluxDifT &sw_bnd_flux_dif ,
  SfcFluxDirVisT &sfc_flux_dir_vis, SfcFluxDirNirT &sfc_flux_dir_nir,
  SfcFluxDifVisT &sfc_flux_dif_vis, SfcFluxDifNirT &sfc_flux_dif_nir)
{
  // Band 10 straddles the near-IR and visible, so divide contributions from band 10 between both broadband sums
  // TODO: Hard-coding these band indices is really bad practice. If the bands ever were to change (like when
  // the RRTMG bands were re-ordered for RRTMGP), we would be using the wrong bands for the IR and UV/VIS. This
  // should be refactored to grab the correct bands by specifying appropriate wavenumber rather than index.
  //sfc_flux_dir_nir(i) = sum(sw_bnd_flux_dir(i+1,kbot,1:9))   + 0.5 * sw_bnd_flux_dir(i+1,kbot,10);
  //sfc_flux_dir_vis(i) = sum(sw_bnd_flux_dir(i+1,kbot,11:14)) + 0.5 * sw_bnd_flux_dir(i+1,kbot,10);
  //sfc_flux_dif_nir(i) = sum(sw_bnd_flux_dif(i+1,kbot,1:9))   + 0.5 * sw_bnd_flux_dif(i+1,kbot,10);
  //sfc_flux_dif_vis(i) = sum(sw_bnd_flux_dif(i+1,kbot,11:14)) + 0.5 * sw_bnd_flux_dif(i+1,kbot,10);

  // Initialize sums over bands
  Kokkos::deep_copy(sfc_flux_dir_nir, 0);
  Kokkos::deep_copy(sfc_flux_dir_vis, 0);
  Kokkos::deep_copy(sfc_flux_dif_nir, 0);
  Kokkos::deep_copy(sfc_flux_dif_vis, 0);

  // Threshold between visible and infrared is 0.7 micron, or 14286 cm^-1.
  const RealT visible_wavenumber_threshold = 14286;
  auto wavenumber_limits = k_dist_sw_k.get_band_lims_wavenumber();
  Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(const int icol) {
    for (int ibnd = 0; ibnd < nswbands; ++ibnd) {
      // Wavenumber is in the visible if it is above the visible wavenumber
      // threshold, and in the infrared if it is below the threshold
      const bool is_visible_wave1 = (wavenumber_limits(0, ibnd) > visible_wavenumber_threshold ? true : false);
      const bool is_visible_wave2 = (wavenumber_limits(1, ibnd) > visible_wavenumber_threshold ? true : false);

      if (is_visible_wave1 && is_visible_wave2) {
        // Entire band is in the visible
        sfc_flux_dir_vis(icol) += sw_bnd_flux_dir(icol,ktop,ibnd);
        sfc_flux_dif_vis(icol) += sw_bnd_flux_dif(icol,ktop,ibnd);
      }
      else if (!is_visible_wave1 && !is_visible_wave2) {
        // Entire band is in the longwave (near-infrared)
        sfc_flux_dir_nir(icol) += sw_bnd_flux_dir(icol,ktop,ibnd);
        sfc_flux_dif_nir(icol) += sw_bnd_flux_dif(icol,ktop,ibnd);
      }
      else {
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

/*
 * Main driver code to run RRTMGP.
 * The input logger is in charge of outputing info to
 * screen and/or to file (or neither), depending on how it was set up.
 */
template <typename PlayT, typename TlayT, typename PlevT, typename TlevT,
          typename SfcAlbDirT, typename SfcAlbDifT, typename Mu0T,
          typename LwpT, typename IwpT, typename RelT, typename ReiT, typename CldfracT,
          typename AerTauSwT, typename AerSsaSwT, typename AerAsmSwT, typename AerTauLwT,
          typename CldTauSwBndT, typename CldTauLwBndT,
          typename CldTauSwGptT, typename CldTauLwGptT,
          typename SwFluxUpT, typename SwFluxDnT, typename SwFluxDnDirT,
          typename LwFluxUpT, typename LwFluxDnT,
          typename SwClnclrskyFluxUpT, typename SwClnclrskyFluxDnT, typename SwClnclrskyFluxDnDirT,
          typename SwClrskyFluxUpT, typename SwClrskyFluxDnT, typename SwClrskyFluxDnDirT,
          typename SwClnskyFluxUpT, typename SwClnskyFluxDnT, typename SwClnskyFluxDnDirT,
          typename LwClnclrskyFluxUpT, typename LwClnclrskyFluxDnT,
          typename LwClrskyFluxUpT, typename LwClrskyFluxDnT,
          typename LwClnskyFluxUpT, typename LwClnskyFluxDnT,
          typename SwBndFluxUpT, typename SwBndFluxDnT, typename SwBndFluxDnDirT,
          typename LwBndFluxUpT, typename LwBndFluxDnT>
static void rrtmgp_main(
  const int ncol, const int nlay,
  PlayT &p_lay, TlayT &t_lay, PlevT &p_lev, TlevT &t_lev,
  GasConcsK<RealT, LayoutT, DeviceT> &gas_concs,
  SfcAlbDirT &sfc_alb_dir, SfcAlbDifT &sfc_alb_dif, Mu0T &mu0,
  LwpT &lwp, IwpT &iwp, RelT &rel, ReiT &rei, CldfracT &cldfrac,
  AerTauSwT &aer_tau_sw, AerSsaSwT &aer_ssa_sw, AerAsmSwT &aer_asm_sw, AerTauLwT &aer_tau_lw,
  CldTauSwBndT &cld_tau_sw_bnd, CldTauLwBndT &cld_tau_lw_bnd,
  CldTauSwGptT &cld_tau_sw_gpt, CldTauLwGptT &cld_tau_lw_gpt,
  SwFluxUpT &sw_flux_up, SwFluxDnT &sw_flux_dn, SwFluxDnDirT &sw_flux_dn_dir,
  LwFluxUpT &lw_flux_up, LwFluxDnT &lw_flux_dn,
  SwClnclrskyFluxUpT &sw_clnclrsky_flux_up, SwClnclrskyFluxDnT &sw_clnclrsky_flux_dn, SwClnclrskyFluxDnDirT &sw_clnclrsky_flux_dn_dir,
  SwClrskyFluxUpT &sw_clrsky_flux_up, SwClrskyFluxDnT &sw_clrsky_flux_dn, SwClrskyFluxDnDirT &sw_clrsky_flux_dn_dir,
  SwClnskyFluxUpT &sw_clnsky_flux_up, SwClnskyFluxDnT &sw_clnsky_flux_dn, SwClnskyFluxDnDirT &sw_clnsky_flux_dn_dir,
  LwClnclrskyFluxUpT &lw_clnclrsky_flux_up, LwClnclrskyFluxDnT &lw_clnclrsky_flux_dn,
  LwClrskyFluxUpT &lw_clrsky_flux_up, LwClrskyFluxDnT &lw_clrsky_flux_dn,
  LwClnskyFluxUpT &lw_clnsky_flux_up, LwClnskyFluxDnT &lw_clnsky_flux_dn,
  SwBndFluxUpT &sw_bnd_flux_up, SwBndFluxDnT &sw_bnd_flux_dn, SwBndFluxDnDirT &sw_bnd_flux_dn_dir,
  LwBndFluxUpT &lw_bnd_flux_up, LwBndFluxDnT &lw_bnd_flux_dn,
  const RealT tsi_scaling,
  const std::shared_ptr<spdlog::logger>& logger,
  const bool extra_clnclrsky_diag = false, const bool extra_clnsky_diag = false)
{
#ifdef SCREAM_RRTMGP_DEBUG
  // Sanity check inputs, and possibly repair
  check_range_k(t_lay      ,  k_dist_sw_k.get_temp_min(),         k_dist_sw_k.get_temp_max(), "rrtmgp_main::t_lay");
  check_range_k(t_lev      ,  k_dist_sw_k.get_temp_min(),         k_dist_sw_k.get_temp_max(), "rrtmgp_main::t_lev");
  check_range_k(p_lay      , k_dist_sw_k.get_press_min(),        k_dist_sw_k.get_press_max(), "rrtmgp_main::p_lay");
  check_range_k(p_lev      , k_dist_sw_k.get_press_min(),        k_dist_sw_k.get_press_max(), "rrtmgp_main::p_lev");
  check_range_k(sfc_alb_dir,                         0,                                1, "rrtmgp_main::sfc_alb_dir");
  check_range_k(sfc_alb_dif,                         0,                                1, "rrtmgp_main::sfc_alb_dif");
  check_range_k(mu0        ,                         0,                                1, "rrtmgp_main::mu0");
  check_range_k(lwp        ,                         0, std::numeric_limits<RealT>::max(), "rrtmgp_main::lwp");
  check_range_k(iwp        ,                         0, std::numeric_limits<RealT>::max(), "rrtmgp_main::iwp");
  check_range_k(rel        ,                         0, std::numeric_limits<RealT>::max(), "rrtmgp_main::rel");
  check_range_k(rei        ,                         0, std::numeric_limits<RealT>::max(), "rrtmgp_main::rei");
#endif

  // Setup pointers to RRTMGP SW fluxes
  FluxesBybandK<RealT, LayoutT, DeviceT> fluxes_sw;
  fluxes_sw.flux_up = sw_flux_up;
  fluxes_sw.flux_dn = sw_flux_dn;
  fluxes_sw.flux_dn_dir = sw_flux_dn_dir;
  fluxes_sw.bnd_flux_up = sw_bnd_flux_up;
  fluxes_sw.bnd_flux_dn = sw_bnd_flux_dn;
  fluxes_sw.bnd_flux_dn_dir = sw_bnd_flux_dn_dir;
  // Clean-clear-sky
  FluxesBroadbandK<RealT, LayoutT, DeviceT> clnclrsky_fluxes_sw;
  clnclrsky_fluxes_sw.flux_up = sw_clnclrsky_flux_up;
  clnclrsky_fluxes_sw.flux_dn = sw_clnclrsky_flux_dn;
  clnclrsky_fluxes_sw.flux_dn_dir = sw_clnclrsky_flux_dn_dir;
  // Clear-sky
  FluxesBroadbandK<RealT, LayoutT, DeviceT> clrsky_fluxes_sw;
  clrsky_fluxes_sw.flux_up = sw_clrsky_flux_up;
  clrsky_fluxes_sw.flux_dn = sw_clrsky_flux_dn;
  clrsky_fluxes_sw.flux_dn_dir = sw_clrsky_flux_dn_dir;
  // Clean-sky
  FluxesBroadbandK<RealT, LayoutT, DeviceT> clnsky_fluxes_sw;
  clnsky_fluxes_sw.flux_up = sw_clnsky_flux_up;
  clnsky_fluxes_sw.flux_dn = sw_clnsky_flux_dn;
  clnsky_fluxes_sw.flux_dn_dir = sw_clnsky_flux_dn_dir;

  // Setup pointers to RRTMGP LW fluxes
  FluxesBybandK<RealT, LayoutT, DeviceT> fluxes_lw;
  fluxes_lw.flux_up = lw_flux_up;
  fluxes_lw.flux_dn = lw_flux_dn;
  fluxes_lw.bnd_flux_up = lw_bnd_flux_up;
  fluxes_lw.bnd_flux_dn = lw_bnd_flux_dn;
  // Clean-clear-sky
  FluxesBroadbandK<RealT, LayoutT, DeviceT> clnclrsky_fluxes_lw;
  clnclrsky_fluxes_lw.flux_up = lw_clnclrsky_flux_up;
  clnclrsky_fluxes_lw.flux_dn = lw_clnclrsky_flux_dn;
  // Clear-sky
  FluxesBroadbandK<RealT, LayoutT, DeviceT> clrsky_fluxes_lw;
  clrsky_fluxes_lw.flux_up = lw_clrsky_flux_up;
  clrsky_fluxes_lw.flux_dn = lw_clrsky_flux_dn;
  // Clean-sky
  FluxesBroadbandK<RealT, LayoutT, DeviceT> clnsky_fluxes_lw;
  clnsky_fluxes_lw.flux_up = lw_clnsky_flux_up;
  clnsky_fluxes_lw.flux_dn = lw_clnsky_flux_dn;

  auto nswbands = k_dist_sw_k.get_nband();
  auto nlwbands = k_dist_lw_k.get_nband();

  // Setup aerosol optical properties
  OpticalProps2strK<RealT, LayoutT, DeviceT> aerosol_sw;
  OpticalProps1sclK<RealT, LayoutT, DeviceT> aerosol_lw;
  aerosol_sw.init(k_dist_sw_k.get_band_lims_wavenumber());
  aerosol_sw.alloc_2str(ncol, nlay);
  Kokkos::parallel_for(MDRP::template get<3>({nswbands,nlay,ncol}) , KOKKOS_LAMBDA (int ibnd, int ilay, int icol) {
    aerosol_sw.tau(icol,ilay,ibnd) = aer_tau_sw(icol,ilay,ibnd);
    aerosol_sw.ssa(icol,ilay,ibnd) = aer_ssa_sw(icol,ilay,ibnd);
    aerosol_sw.g  (icol,ilay,ibnd) = aer_asm_sw(icol,ilay,ibnd);
  });
  aerosol_lw.init(k_dist_lw_k.get_band_lims_wavenumber());
  aerosol_lw.alloc_1scl(ncol, nlay);
  Kokkos::parallel_for(MDRP::template get<3>({nlwbands,nlay,ncol}) , KOKKOS_LAMBDA (int ibnd, int ilay, int icol) {
    aerosol_lw.tau(icol,ilay,ibnd) = aer_tau_lw(icol,ilay,ibnd);
  });

#ifdef SCREAM_RRTMGP_DEBUG
  // Check aerosol optical properties
  // NOTE: these should already have been checked by precondition checks, but someday we might have
  // non-trivial aerosol optics, so this is still good to do here.
  check_range_k(aerosol_sw.tau,  0, 1e3, "rrtmgp_main:aerosol_sw.tau");
  check_range_k(aerosol_sw.ssa,  0,   1, "rrtmgp_main:aerosol_sw.ssa"); //, "aerosol_optics_sw.ssa");
  check_range_k(aerosol_sw.g  , -1,   1, "rrtmgp_main:aerosol_sw.g  "); //, "aerosol_optics_sw.g"  );
  check_range_k(aerosol_lw.tau,  0, 1e3, "rrtmgp_main:aerosol_lw.tau");
#endif

  // Convert cloud physical properties to optical properties for input to RRTMGP
  OpticalProps2strK<RealT, LayoutT, DeviceT> clouds_sw = get_cloud_optics_sw(ncol, nlay, cloud_optics_sw_k, k_dist_sw_k, lwp, iwp, rel, rei);
  OpticalProps1sclK<RealT, LayoutT, DeviceT> clouds_lw = get_cloud_optics_lw(ncol, nlay, cloud_optics_lw_k, k_dist_lw_k, lwp, iwp, rel, rei);
  Kokkos::deep_copy(cld_tau_sw_bnd, clouds_sw.tau);
  Kokkos::deep_copy(cld_tau_lw_bnd, clouds_lw.tau);

  // Do subcolumn sampling to map bands -> gpoints based on cloud fraction and overlap assumption;
  // This implements the Monte Carlo Independing Column Approximation by mapping only a single
  // subcolumn (cloud state) to each gpoint.
  auto nswgpts = k_dist_sw_k.get_ngpt();
  auto clouds_sw_gpt = get_subsampled_clouds(ncol, nlay, nswbands, nswgpts, clouds_sw, k_dist_sw_k, cldfrac, p_lay);
  // Longwave
  auto nlwgpts = k_dist_lw_k.get_ngpt();
  auto clouds_lw_gpt = get_subsampled_clouds(ncol, nlay, nlwbands, nlwgpts, clouds_lw, k_dist_lw_k, cldfrac, p_lay);

  // Copy cloud properties to outputs (is this needed, or can we just use pointers?)
  // Alternatively, just compute and output a subcolumn cloud mask
  Kokkos::parallel_for(MDRP::template get<3>({nswgpts, nlay, ncol}), KOKKOS_LAMBDA (int igpt, int ilay, int icol) {
    cld_tau_sw_gpt(icol,ilay,igpt) = clouds_sw_gpt.tau(icol,ilay,igpt);
  });
  Kokkos::parallel_for(MDRP::template get<3>({nlwgpts, nlay, ncol}), KOKKOS_LAMBDA (int igpt, int ilay, int icol) {
    cld_tau_lw_gpt(icol,ilay,igpt) = clouds_lw_gpt.tau(icol,ilay,igpt);
  });

#ifdef SCREAM_RRTMGP_DEBUG
  // Perform checks on optics; these would be caught by RRTMGP_EXPENSIVE_CHECKS in the RRTMGP code,
  // but we might want to provide additional debug info here. NOTE: we may actually want to move this
  // up higher in the code, I think optical props should go up higher since optical props are kind of
  // a parameterization of their own, and we might want to swap different choices. These checks go here
  // only because we need to run them on computed optical props, so if the optical props themselves get
  // computed up higher, then perform these checks higher as well
  check_range_k(clouds_sw.tau,  0, std::numeric_limits<RealT>::max(), "rrtmgp_main:clouds_sw.tau");
  check_range_k(clouds_sw.ssa,  0,                                1, "rrtmgp_main:clouds_sw.ssa");
  check_range_k(clouds_sw.g  , -1,                                1, "rrtmgp_main:clouds_sw.g  ");
  check_range_k(clouds_sw.tau,  0, std::numeric_limits<RealT>::max(), "rrtmgp_main:clouds_sw.tau");
#endif

  // Do shortwave
  rrtmgp_sw(
    ncol, nlay,
    k_dist_sw_k, p_lay, t_lay, p_lev, t_lev, gas_concs,
    sfc_alb_dir, sfc_alb_dif, mu0, aerosol_sw, clouds_sw_gpt,
    fluxes_sw, clnclrsky_fluxes_sw, clrsky_fluxes_sw, clnsky_fluxes_sw,
    tsi_scaling, logger,
    extra_clnclrsky_diag, extra_clnsky_diag
            );

  // Do longwave
  rrtmgp_lw(
    ncol, nlay,
    k_dist_lw_k, p_lay, t_lay, p_lev, t_lev, gas_concs,
    aerosol_lw, clouds_lw_gpt,
    fluxes_lw, clnclrsky_fluxes_lw, clrsky_fluxes_lw, clnsky_fluxes_lw,
    extra_clnclrsky_diag, extra_clnsky_diag
            );

}

/*
 * Perform any clean-up tasks
 */
static void rrtmgp_finalize()
{
  initialized_k = false;
  k_dist_sw_k.finalize();
  k_dist_lw_k.finalize();
  cloud_optics_sw_k.finalize(); //~CloudOptics();
  cloud_optics_lw_k.finalize(); //~CloudOptics();
  pool_t::finalize();
}

/*
 * Shortwave driver (called by rrtmgp_main)
 */
template <typename PlayT, typename TlayT, typename PlevT, typename TlevT,
          typename SfcAlbDirT, typename SfcAlbDifT, typename Mu0T>
static void rrtmgp_sw(
  const int ncol, const int nlay,
  GasOpticsRRTMGPK<RealT, LayoutT, DeviceT> &k_dist,
  PlayT &p_lay, TlayT &t_lay, PlevT &p_lev, TlevT &t_lev,
  GasConcsK<RealT, LayoutT, DeviceT> &gas_concs,
  SfcAlbDirT &sfc_alb_dir, SfcAlbDifT &sfc_alb_dif, Mu0T &mu0,
  OpticalProps2strK<RealT, LayoutT, DeviceT> &aerosol, OpticalProps2strK<RealT, LayoutT, DeviceT> &clouds,
  FluxesBybandK<RealT, LayoutT, DeviceT> &fluxes, FluxesBroadbandK<RealT, LayoutT, DeviceT> &clnclrsky_fluxes, FluxesBroadbandK<RealT, LayoutT, DeviceT> &clrsky_fluxes, FluxesBroadbandK<RealT, LayoutT, DeviceT> &clnsky_fluxes,
  const RealT tsi_scaling,
  const std::shared_ptr<spdlog::logger>& logger,
  const bool extra_clnclrsky_diag, const bool extra_clnsky_diag)
{
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
  auto &clnclrsky_flux_up = clnclrsky_fluxes.flux_up;
  auto &clnclrsky_flux_dn = clnclrsky_fluxes.flux_dn;
  auto &clnclrsky_flux_dn_dir = clnclrsky_fluxes.flux_dn_dir;
  auto &clrsky_flux_up = clrsky_fluxes.flux_up;
  auto &clrsky_flux_dn = clrsky_fluxes.flux_dn;
  auto &clrsky_flux_dn_dir = clrsky_fluxes.flux_dn_dir;
  auto &clnsky_flux_up = clnsky_fluxes.flux_up;
  auto &clnsky_flux_dn = clnsky_fluxes.flux_dn;
  auto &clnsky_flux_dn_dir = clnsky_fluxes.flux_dn_dir;

  // Reset fluxes to zero
  Kokkos::parallel_for(MDRP::template get<2>({nlay+1,ncol}), KOKKOS_LAMBDA(int ilev, int icol) {
    flux_up    (icol,ilev) = 0;
    flux_dn    (icol,ilev) = 0;
    flux_dn_dir(icol,ilev) = 0;
    clnclrsky_flux_up    (icol,ilev) = 0;
    clnclrsky_flux_dn    (icol,ilev) = 0;
    clnclrsky_flux_dn_dir(icol,ilev) = 0;
    clrsky_flux_up    (icol,ilev) = 0;
    clrsky_flux_dn    (icol,ilev) = 0;
    clrsky_flux_dn_dir(icol,ilev) = 0;
    clnsky_flux_up    (icol,ilev) = 0;
    clnsky_flux_dn    (icol,ilev) = 0;
    clnsky_flux_dn_dir(icol,ilev) = 0;
  });
  Kokkos::parallel_for(MDRP::template get<3>({nbnd,nlay+1,ncol}), KOKKOS_LAMBDA(int ibnd, int ilev, int icol) {
    bnd_flux_up    (icol,ilev,ibnd) = 0;
    bnd_flux_dn    (icol,ilev,ibnd) = 0;
    bnd_flux_dn_dir(icol,ilev,ibnd) = 0;
  });

  // Get daytime indices
  auto dayIndices = view_t<int*>("dayIndices", ncol);
  Kokkos::deep_copy(dayIndices, -1);

  int nday = 0;
  // Serialized for now.
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(int, int& nday_inner) {
    for (int icol = 0; icol < ncol; ++icol) {
      if (mu0(icol) > 0) {
        dayIndices(nday_inner++) = icol;
      }
    }
  }, Kokkos::Sum<int>(nday));

  if (nday == 0) {
    // No daytime columns in this chunk, skip the rest of this routine
    return;
  }

  // Subset mu0
  auto mu0_day = view_t<RealT*>("mu0_day", nday);
  Kokkos::parallel_for(nday, KOKKOS_LAMBDA(int iday) {
    mu0_day(iday) = mu0(dayIndices(iday));
  });

  // subset state variables
  auto p_lay_day = view_t<RealT**>("p_lay_day", nday, nlay);
  auto t_lay_day = view_t<RealT**>("t_lay_day", nday, nlay);
  Kokkos::parallel_for(MDRP::template get<2>({nlay,nday}), KOKKOS_LAMBDA(int ilay, int iday) {
    p_lay_day(iday,ilay) = p_lay(dayIndices(iday),ilay);
    t_lay_day(iday,ilay) = t_lay(dayIndices(iday),ilay);
  });
  auto p_lev_day = view_t<RealT**>("p_lev_day", nday, nlay+1);
  auto t_lev_day = view_t<RealT**>("t_lev_day", nday, nlay+1);
  Kokkos::parallel_for(MDRP::template get<2>({nlay+1,nday}), KOKKOS_LAMBDA(int ilev, int iday) {
    p_lev_day(iday,ilev) = p_lev(dayIndices(iday),ilev);
    t_lev_day(iday,ilev) = t_lev(dayIndices(iday),ilev);
  });

  // Subset gases
  auto gas_names = gas_concs.get_gas_names();
  GasConcsK<RealT, LayoutT, DeviceT> gas_concs_day;
  gas_concs_day.init(gas_names, nday, nlay);
  for (int igas = 0; igas < ngas; igas++) {
    auto vmr_day = view_t<RealT**>("vmr_day", nday, nlay);
    auto vmr     = view_t<RealT**>("vmr"    , ncol, nlay);
    gas_concs.get_vmr(gas_names[igas], vmr);
    Kokkos::parallel_for(MDRP::template get<2>({nlay,nday}), KOKKOS_LAMBDA(int ilay, int iday) {
      vmr_day(iday,ilay) = vmr(dayIndices(iday),ilay);
    });
    gas_concs_day.set_vmr(gas_names[igas], vmr_day);
  }

  // Subset aerosol optics
  OpticalProps2strK<RealT, LayoutT, DeviceT> aerosol_day;
  aerosol_day.init(k_dist.get_band_lims_wavenumber());
  aerosol_day.alloc_2str(nday, nlay);
  Kokkos::parallel_for(MDRP::template get<3>({nbnd,nlay,nday}), KOKKOS_LAMBDA(int ibnd, int ilay, int iday) {
    aerosol_day.tau(iday,ilay,ibnd) = aerosol.tau(dayIndices(iday),ilay,ibnd);
    aerosol_day.ssa(iday,ilay,ibnd) = aerosol.ssa(dayIndices(iday),ilay,ibnd);
    aerosol_day.g  (iday,ilay,ibnd) = aerosol.g  (dayIndices(iday),ilay,ibnd);
  });

  // Subset cloud optics
  // TODO: nbnd -> ngpt once we pass sub-sampled cloud state
  OpticalProps2strK<RealT, LayoutT, DeviceT> clouds_day;
  clouds_day.init(k_dist.get_band_lims_wavenumber(), k_dist.get_band_lims_gpoint());
  clouds_day.alloc_2str(nday, nlay);
  Kokkos::parallel_for(MDRP::template get<3>({ngpt,nlay,nday}), KOKKOS_LAMBDA(int igpt, int ilay, int iday) {
    clouds_day.tau(iday,ilay,igpt) = clouds.tau(dayIndices(iday),ilay,igpt);
    clouds_day.ssa(iday,ilay,igpt) = clouds.ssa(dayIndices(iday),ilay,igpt);
    clouds_day.g  (iday,ilay,igpt) = clouds.g  (dayIndices(iday),ilay,igpt);
  });

  // RRTMGP assumes surface albedos have a screwy dimension ordering
  // for some strange reason, so we need to transpose these; also do
  // daytime subsetting in the same kernel
  view_t<RealT**> sfc_alb_dir_T("sfc_alb_dir", nbnd, nday);
  view_t<RealT**> sfc_alb_dif_T("sfc_alb_dif", nbnd, nday);
  Kokkos::parallel_for(MDRP::template get<2>({nbnd,nday}), KOKKOS_LAMBDA(int ibnd, int icol) {
    sfc_alb_dir_T(ibnd,icol) = sfc_alb_dir(dayIndices(icol),ibnd);
    sfc_alb_dif_T(ibnd,icol) = sfc_alb_dif(dayIndices(icol),ibnd);
  });

  // Temporaries we need for daytime-only fluxes
  auto flux_up_day = view_t<RealT**>("flux_up_day", nday, nlay+1);
  auto flux_dn_day = view_t<RealT**>("flux_dn_day", nday, nlay+1);
  auto flux_dn_dir_day = view_t<RealT**>("flux_dn_dir_day", nday, nlay+1);
  auto bnd_flux_up_day = view_t<RealT***>("bnd_flux_up_day", nday, nlay+1, nbnd);
  auto bnd_flux_dn_day = view_t<RealT***>("bnd_flux_dn_day", nday, nlay+1, nbnd);
  auto bnd_flux_dn_dir_day = view_t<RealT***>("bnd_flux_dn_dir_day", nday, nlay+1, nbnd);
  FluxesBybandK<RealT, LayoutT, DeviceT> fluxes_day;
  fluxes_day.flux_up         = flux_up_day;
  fluxes_day.flux_dn         = flux_dn_day;
  fluxes_day.flux_dn_dir     = flux_dn_dir_day;
  fluxes_day.bnd_flux_up     = bnd_flux_up_day;
  fluxes_day.bnd_flux_dn     = bnd_flux_dn_day;
  fluxes_day.bnd_flux_dn_dir = bnd_flux_dn_dir_day;

  // Allocate space for optical properties
  OpticalProps2strK<RealT, LayoutT, DeviceT> optics;
  optics.alloc_2str(nday, nlay, k_dist);

  OpticalProps2strK<RealT, LayoutT, DeviceT> optics_no_aerosols;
  if (extra_clnsky_diag) {
    // Allocate space for optical properties (no aerosols)
    optics_no_aerosols.alloc_2str(nday, nlay, k_dist);
  }

  // Limit temperatures for gas optics look-up tables
  auto t_lay_limited = view_t<RealT**>("t_lay_limited", nday, nlay);
  limit_to_bounds_k(t_lay_day, k_dist_sw_k.get_temp_min(), k_dist_sw_k.get_temp_max(), t_lay_limited);

  // Do gas optics
  view_t<RealT**> toa_flux("toa_flux", nday, ngpt);
  bool top_at_1 = false;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(int, bool& val) {
    val |= p_lay(0, 0) < p_lay(0, nlay-1);
  }, Kokkos::LOr<bool>(top_at_1));

  oview_t<RealT***> col_gas("col_gas", std::make_pair(0, ncol-1), std::make_pair(0, nlay-1), std::make_pair(-1, k_dist.get_ngas()-1));

  k_dist.gas_optics(nday, nlay, top_at_1, p_lay_day, p_lev_day, t_lay_limited, gas_concs_day, col_gas, optics, toa_flux);
  if (extra_clnsky_diag) {
    k_dist.gas_optics(nday, nlay, top_at_1, p_lay_day, p_lev_day, t_lay_limited, gas_concs_day, col_gas, optics_no_aerosols, toa_flux);
  }

#ifdef SCREAM_RRTMGP_DEBUG
  // Check gas optics
  check_range_k(optics.tau,  0, std::numeric_limits<RealT>::max(), "rrtmgp_sw:optics.tau");
  check_range_k(optics.ssa,  0,                                1, "rrtmgp_sw:optics.ssa"); //, "optics.ssa");
  check_range_k(optics.g  , -1,                                1, "rrtmgp_sw:optics.g  "); //, "optics.g"  );
#endif

  // Apply tsi_scaling
  Kokkos::parallel_for(MDRP::template get<2>({ngpt,nday}), KOKKOS_LAMBDA(int igpt, int iday) {
    toa_flux(iday,igpt) = tsi_scaling * toa_flux(iday,igpt);
  });

  if (extra_clnclrsky_diag) {
    // Compute clear-clean-sky (just gas) fluxes on daytime columns
    rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);
    // Expand daytime fluxes to all columns
    Kokkos::parallel_for(MDRP::template get<2>({nlay+1,nday}), KOKKOS_LAMBDA(int ilev, int iday) {
      const int icol = dayIndices(iday);
      clnclrsky_flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
      clnclrsky_flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
      clnclrsky_flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
    });
  }

  // Combine gas and aerosol optics
  aerosol_day.delta_scale();
  aerosol_day.increment(optics);

  // Compute clearsky (gas + aerosol) fluxes on daytime columns
  rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);

  // Expand daytime fluxes to all columns
  Kokkos::parallel_for(MDRP::template get<2>({nlay+1,nday}), KOKKOS_LAMBDA(int ilev, int iday) {
    const int icol = dayIndices(iday);
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
  Kokkos::parallel_for(MDRP::template get<2>({nlay+1,nday}), KOKKOS_LAMBDA(int ilev, int iday) {
    const int icol = dayIndices(iday);
    flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
    flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
    flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
  });
  Kokkos::parallel_for(MDRP::template get<3>({nbnd,nlay+1,nday}), KOKKOS_LAMBDA(int ibnd, int ilev, int iday) {
    const int icol = dayIndices(iday);
    bnd_flux_up    (icol,ilev,ibnd) = bnd_flux_up_day    (iday,ilev,ibnd);
    bnd_flux_dn    (icol,ilev,ibnd) = bnd_flux_dn_day    (iday,ilev,ibnd);
    bnd_flux_dn_dir(icol,ilev,ibnd) = bnd_flux_dn_dir_day(iday,ilev,ibnd);
  });

  if (extra_clnsky_diag) {
    // First increment clouds in optics_no_aerosols
    clouds_day.increment(optics_no_aerosols);
    // Compute cleansky (gas + clouds) fluxes on daytime columns
    rte_sw(optics_no_aerosols, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);
    // Expand daytime fluxes to all columns
    Kokkos::parallel_for(MDRP::template get<2>({nlay+1,nday}), KOKKOS_LAMBDA(int ilev, int iday) {
      const int icol = dayIndices(iday);
      clnsky_flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
      clnsky_flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
      clnsky_flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
    });
  }
}

/*
 * Longwave driver (called by rrtmgp_main)
 */
template <typename PlayT, typename TlayT, typename PlevT, typename TlevT>
static void rrtmgp_lw(
  const int ncol, const int nlay,
  GasOpticsRRTMGPK<RealT, LayoutT, DeviceT> &k_dist,
  PlayT &p_lay, TlayT &t_lay, PlevT &p_lev, TlevT &t_lev,
  GasConcsK<RealT, LayoutT, DeviceT> &gas_concs,
  OpticalProps1sclK<RealT, LayoutT, DeviceT> &aerosol, OpticalProps1sclK<RealT, LayoutT, DeviceT> &clouds,
  FluxesBybandK<RealT, LayoutT, DeviceT> &fluxes, FluxesBroadbandK<RealT, LayoutT, DeviceT> &clnclrsky_fluxes, FluxesBroadbandK<RealT, LayoutT, DeviceT> &clrsky_fluxes, FluxesBroadbandK<RealT, LayoutT, DeviceT> &clnsky_fluxes,
  const bool extra_clnclrsky_diag, const bool extra_clnsky_diag)
{
  // Problem size
  int nbnd = k_dist.get_nband();

  // Associate local pointers for fluxes
  auto &flux_up           = fluxes.flux_up;
  auto &flux_dn           = fluxes.flux_dn;
  auto &bnd_flux_up       = fluxes.bnd_flux_up;
  auto &bnd_flux_dn       = fluxes.bnd_flux_dn;
  auto &clnclrsky_flux_up = clnclrsky_fluxes.flux_up;
  auto &clnclrsky_flux_dn = clnclrsky_fluxes.flux_dn;
  auto &clrsky_flux_up    = clrsky_fluxes.flux_up;
  auto &clrsky_flux_dn    = clrsky_fluxes.flux_dn;
  auto &clnsky_flux_up    = clnsky_fluxes.flux_up;
  auto &clnsky_flux_dn    = clnsky_fluxes.flux_dn;

  // Reset fluxes to zero
  Kokkos::parallel_for(
    MDRP::template get<2>({nlay + 1, ncol}), KOKKOS_LAMBDA(int ilev, int icol) {
      flux_up(icol, ilev)           = 0;
      flux_dn(icol, ilev)           = 0;
      clnclrsky_flux_up(icol, ilev) = 0;
      clnclrsky_flux_dn(icol, ilev) = 0;
      clrsky_flux_up(icol, ilev)    = 0;
      clrsky_flux_dn(icol, ilev)    = 0;
      clnsky_flux_up(icol, ilev)    = 0;
      clnsky_flux_dn(icol, ilev)    = 0;
    });
  Kokkos::parallel_for(
    MDRP::template get<3>({nbnd, nlay + 1, ncol}),
    KOKKOS_LAMBDA(int ibnd, int ilev, int icol) {
      bnd_flux_up(icol, ilev, ibnd) = 0;
      bnd_flux_dn(icol, ilev, ibnd) = 0;
    });

  // Allocate space for optical properties
  OpticalProps1sclK<RealT, LayoutT, DeviceT> optics;
  optics.alloc_1scl(ncol, nlay, k_dist);
  OpticalProps1sclK<RealT, LayoutT, DeviceT> optics_no_aerosols;
  if (extra_clnsky_diag) {
    // Allocate space for optical properties (no aerosols)
    optics_no_aerosols.alloc_1scl(ncol, nlay, k_dist);
  }

  // Boundary conditions
  SourceFuncLWK<RealT, LayoutT, DeviceT> lw_sources;
  lw_sources.alloc(ncol, nlay, k_dist);
  view_t<RealT*> t_sfc   ("t_sfc"        ,ncol);
  view_t<RealT**> emis_sfc("emis_sfc",nbnd,ncol);

  bool top_at_1 = false;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(int, bool& val) {
    val |= p_lay(0, 0) < p_lay(0, nlay-1);
  }, Kokkos::LOr<bool>(top_at_1));

  // Surface temperature
  Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    t_sfc(icol) = t_lev(icol, conv::merge(nlay, 0, top_at_1));
  });
  Kokkos::deep_copy(emis_sfc , 0.98);

  // Get Gaussian quadrature weights
  // TODO: move this crap out of userland!
  // Weights and angle secants for first order (k=1) Gaussian quadrature.
  //   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
  //   after Abramowitz & Stegun 1972, page 921
  int constexpr max_gauss_pts = 4;
  hview_t<RealT**> gauss_Ds_host ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
  gauss_Ds_host(0,0) = 1.66      ; gauss_Ds_host(1,0) =         0.; gauss_Ds_host(2,0) =         0.; gauss_Ds_host(3,0) =         0.;
  gauss_Ds_host(0,1) = 1.18350343; gauss_Ds_host(1,1) = 2.81649655; gauss_Ds_host(2,1) =         0.; gauss_Ds_host(3,1) =         0.;
  gauss_Ds_host(0,2) = 1.09719858; gauss_Ds_host(1,2) = 1.69338507; gauss_Ds_host(2,2) = 4.70941630; gauss_Ds_host(3,2) =         0.;
  gauss_Ds_host(0,3) = 1.06056257; gauss_Ds_host(1,3) = 1.38282560; gauss_Ds_host(2,3) = 2.40148179; gauss_Ds_host(3,3) = 7.15513024;

  hview_t<RealT**> gauss_wts_host("gauss_wts",max_gauss_pts,max_gauss_pts);
  gauss_wts_host(0,0) = 0.5         ; gauss_wts_host(1,0) = 0.          ; gauss_wts_host(2,0) = 0.          ; gauss_wts_host(3,0) = 0.          ;
  gauss_wts_host(0,1) = 0.3180413817; gauss_wts_host(1,1) = 0.1819586183; gauss_wts_host(2,1) = 0.          ; gauss_wts_host(3,1) = 0.          ;
  gauss_wts_host(0,2) = 0.2009319137; gauss_wts_host(1,2) = 0.2292411064; gauss_wts_host(2,2) = 0.0698269799; gauss_wts_host(3,2) = 0.          ;
  gauss_wts_host(0,3) = 0.1355069134; gauss_wts_host(1,3) = 0.2034645680; gauss_wts_host(2,3) = 0.1298475476; gauss_wts_host(3,3) = 0.0311809710;

  view_t<RealT**> gauss_Ds ("gauss_Ds" ,max_gauss_pts,max_gauss_pts);
  view_t<RealT**> gauss_wts("gauss_wts",max_gauss_pts,max_gauss_pts);
  Kokkos::deep_copy(gauss_Ds,  gauss_Ds_host);
  Kokkos::deep_copy(gauss_wts, gauss_wts_host);

  // Limit temperatures for gas optics look-up tables
  auto t_lay_limited = view_t<RealT**>("t_lay_limited", ncol, nlay);
  auto t_lev_limited = view_t<RealT**>("t_lev_limited", ncol, nlay+1);
  limit_to_bounds_k(t_lay, k_dist_lw_k.get_temp_min(), k_dist_lw_k.get_temp_max(), t_lay_limited);
  limit_to_bounds_k(t_lev, k_dist_lw_k.get_temp_min(), k_dist_lw_k.get_temp_max(), t_lev_limited);

  // Do gas optics
  oview_t<RealT***> col_gas("col_gas", std::make_pair(0, ncol-1), std::make_pair(0, nlay-1), std::make_pair(-1, k_dist.get_ngas()-1));
  k_dist.gas_optics(ncol, nlay, top_at_1, p_lay, p_lev, t_lay_limited, t_sfc, gas_concs, col_gas, optics, lw_sources, view_t<RealT**>(), t_lev_limited);
  if (extra_clnsky_diag) {
    k_dist.gas_optics(ncol, nlay, top_at_1, p_lay, p_lev, t_lay_limited, t_sfc, gas_concs, col_gas, optics_no_aerosols, lw_sources, view_t<RealT**>(), t_lev_limited);
  }

#ifdef SCREAM_RRTMGP_DEBUG
  // Check gas optics
  check_range_k(optics.tau,  0, std::numeric_limits<RealT>::max(), "rrtmgp_lw:optics.tau");
#endif

  if (extra_clnclrsky_diag) {
    // Compute clean-clear-sky fluxes before we add in aerosols and clouds
    rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics, top_at_1, lw_sources, emis_sfc, clnclrsky_fluxes);
  }

  // Combine gas and aerosol optics
  aerosol.increment(optics);

  // Compute clear-sky fluxes before we add in clouds
  rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics, top_at_1, lw_sources, emis_sfc, clrsky_fluxes);

  // Combine gas and cloud optics
  clouds.increment(optics);

  // Compute allsky fluxes
  rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics, top_at_1, lw_sources, emis_sfc, fluxes);

  if (extra_clnsky_diag) {
    // First increment clouds in optics_no_aerosols
    clouds.increment(optics_no_aerosols);
    // Compute clean-sky fluxes
    rte_lw(max_gauss_pts, gauss_Ds, gauss_wts, optics_no_aerosols, top_at_1, lw_sources, emis_sfc, clnsky_fluxes);
  }
}

/*
 * Return a subcolumn mask consistent with a specified overlap assumption
 */
template <typename CldfT, typename SeedsT>
static view_t<int***> get_subcolumn_mask(const int ncol, const int nlay, const int ngpt,
                                         CldfT &cldf, const int overlap_option, SeedsT &seeds)
{
  // Routine will return subcolumn mask with values of 0 indicating no cloud, 1 indicating cloud
  auto subcolumn_mask = view_t<int***>("subcolumn_mask", ncol, nlay, ngpt);

  // Subcolumn generators are a means for producing a variable x(i,j,k), where
  //
  //     c(i,j,k) = 1 for x(i,j,k) >  1 - cldf(i,j)
  //     c(i,j,k) = 0 for x(i,j,k) <= 1 - cldf(i,j)
  //
  // I am going to call this "cldx" to be just slightly less ambiguous
  auto cldx = view_t<RealT***>("cldx", ncol, nlay, ngpt);

  // Apply overlap assumption to set cldx
  if (overlap_option == 0) {  // Dummy mask, always cloudy
    Kokkos::deep_copy(cldx, 1);
  } else {  // Default case, maximum-random overlap
    // Maximum-random overlap:
    // Uses essentially the algorithm described in eq (14) in Raisanen et al. 2004,
    // https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1256/qj.03.99. Also the same
    // algorithm used in RRTMG implementation of maximum-random overlap (see
    // https://github.com/AER-RC/RRTMG_SW/blob/master/src/mcica_subcol_gen_sw.f90)
    //
    // First, fill cldx with random numbers. Need to use a unique seed for each column!
    // auto seeds_host = Kokkos::create_mirror_view(seeds);
    // Kokkos::deep_copy(seeds_host, seeds);
    // for (int icol = 0; icol < ncol; ++icol) {
    //   Kokkos::Random_XorShift64_Pool<> random_pool(seeds_host(icol));
    //   Kokkos::parallel_for(MDRP::template get<2>({ngpt, nlay}), KOKKOS_LAMBDA(int igpt, int ilay) {
    //     auto generator = random_pool.get_state();
    //     cldx(icol,ilay,igpt) = generator.drand(0., 1.);
    //     random_pool.free_state(generator);
    //   });
    // }
    Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
      conv::Random rand(seeds(icol));
      for (int igpt = 0; igpt < ngpt; igpt++) {
        for (int ilay = 0; ilay < nlay; ilay++) {
          cldx(icol,ilay,igpt) = rand.genFP<RealT>();
        }
      }
    });

    // Step down columns and apply algorithm from eq (14)
    Kokkos::parallel_for(MDRP::template get<2>({ngpt,ncol}), KOKKOS_LAMBDA(int igpt, int icol) {
      for (int ilay = 1; ilay < nlay; ilay++) {
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
  Kokkos::parallel_for(MDRP::template get<3>({ngpt,nlay,ncol}), KOKKOS_LAMBDA(int igpt, int ilay, int icol) {
    if (cldx(icol,ilay,igpt) > 1.0 - cldf(icol,ilay)) {
      subcolumn_mask(icol,ilay,igpt) = 1;
    } else {
      subcolumn_mask(icol,ilay,igpt) = 0;
    }
  });
  return subcolumn_mask;
}

/*
 * Compute cloud area from 3d subcol cloud property
 */
template <typename PmidT, typename CldTauGptT, typename CldAreaT>
static void compute_cloud_area(
  int ncol, int nlay, int ngpt, RealT pmin, RealT pmax,
  const PmidT& pmid, const CldTauGptT& cld_tau_gpt, CldAreaT& cld_area)
{
  // Subcolumn binary cld mask; if any layers with pressure between pmin and pmax are cloudy
  // then 2d subcol mask is 1, otherwise it is 0
  auto subcol_mask = view_t<RealT**>("subcol_mask", ncol, ngpt);
  Kokkos::parallel_for(MDRP::template get<3>({ngpt, nlay, ncol}), KOKKOS_LAMBDA(int igpt, int ilay, int icol) {
    // NOTE: using plev would need to assume level ordering (top to bottom or bottom to top), but
    // using play/pmid does not
    if (cld_tau_gpt(icol,ilay,igpt) > 0 && pmid(icol,ilay) >= pmin && pmid(icol,ilay) < pmax) {
      subcol_mask(icol,igpt) = 1;
    }
  });
  // Compute average over subcols to get cloud area
  auto ngpt_inv = 1.0 / ngpt;
  Kokkos::deep_copy(cld_area, 0);
  Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    // This loop needs to be serial because of the atomic reduction
    for (int igpt = 0; igpt < ngpt; ++igpt) {
      cld_area(icol) += subcol_mask(icol,igpt) * ngpt_inv;
    }
  });
}

/*
 * Return select cloud-top diagnostics following AeroCom recommendation
 */
template <typename TmidT, typename PmidT, typename PdelT, typename ZdelT, typename QcT,
          typename QiT, typename RelT, typename ReiT, typename CldfracTotT, typename NcT,
          typename TmidAtCldtopT, typename PmidAtCldtopT, typename CldfracIceAtCldtopT,
          typename CldfracLiqAtCldtopT, typename CldfracTotAtCldtopT, typename CdncAtCldtopT,
          typename EffRadiusQcAtCldtopT, typename EffRadiusQiAtCldTopT>
static void compute_aerocom_cloudtop(
  int ncol, int nlay, const TmidT &tmid, const PmidT &pmid,
  const PdelT &p_del, const ZdelT &z_del, const QcT &qc,
  const QiT &qi, const RelT &rel, const ReiT &rei,
  const CldfracTotT &cldfrac_tot, const NcT &nc,
  TmidAtCldtopT &T_mid_at_cldtop, PmidAtCldtopT &p_mid_at_cldtop,
  CldfracIceAtCldtopT &cldfrac_ice_at_cldtop, CldfracLiqAtCldtopT &cldfrac_liq_at_cldtop,
  CldfracTotAtCldtopT &cldfrac_tot_at_cldtop, CdncAtCldtopT &cdnc_at_cldtop,
  EffRadiusQcAtCldtopT &eff_radius_qc_at_cldtop, EffRadiusQiAtCldTopT &eff_radius_qi_at_cldtop)
{
  /* The goal of this routine is to calculate properties at cloud top
   * based on the AeroCom recommendation. See reference for routine
   * get_subcolumn_mask above, where equation 14 is used for the
   * maximum-random overlap assumption for subcolumn generation. We use
   * equation 13, the column counterpart.
   */
  // Set outputs to zero
  Kokkos::deep_copy(T_mid_at_cldtop, 0.0);
  Kokkos::deep_copy(p_mid_at_cldtop, 0.0);
  Kokkos::deep_copy(cldfrac_ice_at_cldtop, 0.0);
  Kokkos::deep_copy(cldfrac_liq_at_cldtop, 0.0);
  Kokkos::deep_copy(cldfrac_tot_at_cldtop, 0.0);
  Kokkos::deep_copy(cdnc_at_cldtop, 0.0);
  Kokkos::deep_copy(eff_radius_qc_at_cldtop, 0.0);
  Kokkos::deep_copy(eff_radius_qi_at_cldtop, 0.0);

  // Initialize the 1D "clear fraction" as 1 (totally clear)
  auto aerocom_clr = view_t<RealT*>("aerocom_clr", ncol);
  Kokkos::deep_copy(aerocom_clr, 1.0);

  // Get gravity acceleration constant from constants
  using physconst = scream::physics::Constants<RealT>;

  // TODO: move tunable constant to namelist
  constexpr RealT q_threshold = 0.0;  // BAD_CONSTANT!

  // TODO: move tunable constant to namelist
  constexpr RealT cldfrac_tot_threshold = 0.001;  // BAD_CONSTANT!

  // Loop over all columns in parallel
  Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    // Loop over all layers in serial (due to accumulative
    // product), starting at 2 (second highest) layer because the
    // highest is assumed to hav no clouds
    for(int ilay = 1; ilay < nlay; ++ilay) {
      // Only do the calculation if certain conditions are met
      if((qc(icol, ilay) + qi(icol, ilay)) > q_threshold &&
         (cldfrac_tot(icol, ilay) > cldfrac_tot_threshold)) {
        /* PART I: Probabilistically determining cloud top */
        // Populate aerocom_tmp as the clear-sky fraction
        // probability of this level, where aerocom_clr is that of
        // the previous level
        auto aerocom_tmp =
          aerocom_clr(icol) *
          (1.0 - ekat::impl::max(cldfrac_tot(icol, ilay - 1),
                                 cldfrac_tot(icol, ilay))) /
          (1.0 - ekat::impl::min(cldfrac_tot(icol, ilay - 1),
                                 1.0 - cldfrac_tot_threshold));
        // Temporary variable for probability "weights"
        auto aerocom_wts = aerocom_clr(icol) - aerocom_tmp;
        // Temporary variable for liquid "phase"
        auto aerocom_phi =
          qc(icol, ilay) / (qc(icol, ilay) + qi(icol, ilay));
        /* PART II: The inferred properties */
        /* In general, converting a 3D property X to a 2D cloud-top
         * counterpart x follows: x(i) += X(i,k) * weights * Phase
         * but X and Phase are not always needed */
        // T_mid_at_cldtop
        T_mid_at_cldtop(icol) += tmid(icol, ilay) * aerocom_wts;
        // p_mid_at_cldtop
        p_mid_at_cldtop(icol) += pmid(icol, ilay) * aerocom_wts;
        // cldfrac_ice_at_cldtop
        cldfrac_ice_at_cldtop(icol) +=
          (1.0 - aerocom_phi) * aerocom_wts;
        // cldfrac_liq_at_cldtop
        cldfrac_liq_at_cldtop(icol) += aerocom_phi * aerocom_wts;
        // cdnc_at_cldtop
        /* We need to convert nc from 1/mass to 1/volume first, and
         * from grid-mean to in-cloud, but after that, the
         * calculation follows the general logic */
        auto cdnc = nc(icol, ilay) * p_del(icol, ilay) /
          z_del(icol, ilay) / physconst::gravit /
          cldfrac_tot(icol, ilay);
        cdnc_at_cldtop(icol) += cdnc * aerocom_phi * aerocom_wts;
        // eff_radius_qc_at_cldtop
        eff_radius_qc_at_cldtop(icol) +=
          rel(icol, ilay) * aerocom_phi * aerocom_wts;
        // eff_radius_qi_at_cldtop
        eff_radius_qi_at_cldtop(icol) +=
          rei(icol, ilay) * (1.0 - aerocom_phi) * aerocom_wts;
        // Reset aerocom_clr to aerocom_tmp to accumulate
        aerocom_clr(icol) = aerocom_tmp;
      }
    }
    // After the serial loop over levels, the cloudy fraction is
    // defined as (1 - aerocom_clr). This is true because
    // aerocom_clr is the result of accumulative probabilities
    // (their products)
    cldfrac_tot_at_cldtop(icol) = 1.0 - aerocom_clr(icol);
  });
}

/*
 * Provide a function to convert cloud (water and ice) mixing ratios to layer mass per unit area
 * (what E3SM refers to as "in-cloud water paths", a terminology we shun here to avoid confusion
 * with the standard practice of using "water path" to refer to the total column-integrated
 * quantities).
 */
template<class View1, class View2, class View3, class View4>
static void mixing_ratio_to_cloud_mass(
  View1 const& mixing_ratio,
  View2 const& cloud_fraction,
  View3 const& dp,
  View4 const& cloud_mass)
{
  int ncol = mixing_ratio.extent(0);
  int nlay = mixing_ratio.extent(1);
  using physconst = scream::physics::Constants<Real>;
  Kokkos::parallel_for(MDRP::template get<2>({nlay, ncol}), KOKKOS_LAMBDA(int ilay, int icol) {
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
template<typename InT, typename T, typename OutT, typename std::enable_if<OutT::rank == 1>::type* = nullptr>
static void limit_to_bounds_k(InT const &arr_in, T const lower, T const upper, OutT &arr_out) {
  Kokkos::parallel_for(arr_out.size(), KOKKOS_LAMBDA(int i) {
    arr_out(i) = std::min(std::max(arr_in(i), lower), upper);
  });
}

template<typename InT, typename T, typename OutT, typename std::enable_if<OutT::rank == 2>::type* = nullptr>
static void limit_to_bounds_k(InT const &arr_in, T const lower, T const upper, OutT &arr_out) {
  Kokkos::parallel_for(MDRP::template get<2>({arr_out.extent(0), arr_out.extent(1)}), KOKKOS_LAMBDA(int i, int j) {
    arr_out(i, j) = std::min(std::max(arr_in(i, j), lower), upper);
  });
}


static int get_wavelength_index(OpticalPropsK<RealT, LayoutT, DeviceT> &kdist, RealT wavelength)
{
  // Get wavelength bounds for all wavelength bands
  auto wavelength_bounds = kdist.get_band_lims_wavelength();

  // Find the band index for the specified wavelength
  // Note that bands are stored in wavenumber space, units of cm-1, so if we are passed wavelength
  // in units of meters, we need a conversion factor of 10^2
  const int nbnds = kdist.get_nband();
  int band_index = -1;
  Kokkos::parallel_reduce(nbnds, KOKKOS_LAMBDA(int ibnd, int& band_index_inner) {
    if (wavelength_bounds(0,ibnd) < wavelength_bounds(1,ibnd)) {
      if (wavelength_bounds(0,ibnd) <= wavelength * 1e2 && wavelength * 1e2 <= wavelength_bounds(1,ibnd)) {
        band_index_inner = ibnd;
      }
    } else {
      if (wavelength_bounds(0,ibnd) >= wavelength * 1e2 && wavelength * 1e2 >= wavelength_bounds(1,ibnd)) {
        band_index_inner = ibnd;
      }
    }
  }, Kokkos::Max<int>(band_index));
  return band_index;
}

static inline int get_wavelength_index_sw_k(RealT wavelength)
{
  return get_wavelength_index(k_dist_sw_k, wavelength);
}

static inline int get_wavelength_index_lw_k(RealT wavelength)
{
  return get_wavelength_index(k_dist_lw_k, wavelength);
}

template <typename LwpT, typename IwpT, typename RelT, typename ReiT>
static OpticalProps2strK<RealT, LayoutT, DeviceT> get_cloud_optics_sw(
  const int ncol, const int nlay,
  CloudOpticsK<RealT, LayoutT, DeviceT> &cloud_optics, GasOpticsRRTMGPK<RealT, LayoutT, DeviceT> &kdist,
  LwpT &lwp, IwpT &iwp, RelT &rel, ReiT &rei) {

  // Initialize optics
  OpticalProps2strK<RealT, LayoutT, DeviceT> clouds;
  clouds.init(kdist.get_band_lims_wavenumber());
  clouds.alloc_2str(ncol, nlay);

  // Needed for consistency with all-sky example problem?
  cloud_optics.set_ice_roughness(2);

  // Limit effective radii to be within bounds of lookup table
  auto rel_limited = view_t<RealT**>("rel_limited", ncol, nlay);
  auto rei_limited = view_t<RealT**>("rei_limited", ncol, nlay);
  limit_to_bounds_k(rel, cloud_optics.radliq_lwr, cloud_optics.radliq_upr, rel_limited);
  limit_to_bounds_k(rei, cloud_optics.radice_lwr, cloud_optics.radice_upr, rei_limited);

  // Calculate cloud optics
  cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel_limited, rei_limited, clouds);

  // Return optics
  return clouds;
}

template <typename LwpT, typename IwpT, typename RelT, typename ReiT>
static OpticalProps1sclK<RealT, LayoutT, DeviceT> get_cloud_optics_lw(
  const int ncol, const int nlay,
  CloudOpticsK<RealT, LayoutT, DeviceT> &cloud_optics, GasOpticsRRTMGPK<RealT, LayoutT, DeviceT> &kdist,
  LwpT &lwp, IwpT &iwp, RelT &rel, ReiT &rei) {

  // Initialize optics
  OpticalProps1sclK<RealT, LayoutT, DeviceT> clouds;
  clouds.init(kdist.get_band_lims_wavenumber());
  clouds.alloc_1scl(ncol, nlay);  // this is dumb, why do we need to init and alloc separately?!

  // Needed for consistency with all-sky example problem?
  cloud_optics.set_ice_roughness(2);

  // Limit effective radii to be within bounds of lookup table
  auto rel_limited = view_t<RealT**>("rel_limited", ncol, nlay);
  auto rei_limited = view_t<RealT**>("rei_limited", ncol, nlay);
  limit_to_bounds_k(rel, cloud_optics.radliq_lwr, cloud_optics.radliq_upr, rel_limited);
  limit_to_bounds_k(rei, cloud_optics.radice_lwr, cloud_optics.radice_upr, rei_limited);

  // Calculate cloud optics
  cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel_limited, rei_limited, clouds);

  // Return optics
  return clouds;
}

template <typename CldT, typename PlayT>
static OpticalProps2strK<RealT, LayoutT, DeviceT> get_subsampled_clouds(
  const int ncol, const int nlay, const int nbnd, const int ngpt,
  OpticalProps2strK<RealT, LayoutT, DeviceT> &cloud_optics, GasOpticsRRTMGPK<RealT, LayoutT, DeviceT> &kdist, CldT &cld, PlayT &p_lay) {
  // Initialized subsampled optics
  OpticalProps2strK<RealT, LayoutT, DeviceT> subsampled_optics;
  subsampled_optics.init(kdist.get_band_lims_wavenumber(), kdist.get_band_lims_gpoint(), "subsampled_optics");
  subsampled_optics.alloc_2str(ncol, nlay);
  // Check that we do not have clouds with no optical properties; this would get corrected
  // when we assign optical props, but we want to use a "radiative cloud fraction"
  // for the subcolumn sampling too because otherwise we can get vertically-contiguous cloud
  // mask profiles with no actual cloud properties in between, which would just further overestimate
  // the vertical correlation of cloudy layers. I.e., cloudy layers might look maximally overlapped
  // even when separated by layers with no cloud properties, when in fact those layers should be
  // randomly overlapped.
  auto cldfrac_rad = view_t<RealT**>("cldfrac_rad", ncol, nlay);
  Kokkos::parallel_for(MDRP::template get<3>({nbnd,nlay,ncol}), KOKKOS_LAMBDA (int ibnd, int ilay, int icol) {
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
  auto seeds = view_t<int*>("seeds", ncol);
  Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    seeds(icol) = 1e9 * (p_lay(icol,nlay-1) - int(p_lay(icol,nlay-1)));
  });
  auto cldmask = get_subcolumn_mask(ncol, nlay, ngpt, cldfrac_rad, overlap, seeds);
  // Assign optical properties to subcolumns (note this implements MCICA)
  auto gpoint_bands = kdist.get_gpoint_bands();
  Kokkos::parallel_for(MDRP::template get<3>({ngpt,nlay,ncol}), KOKKOS_LAMBDA(int igpt, int ilay, int icol) {
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

template <typename CldT, typename PlayT>
static OpticalProps1sclK<RealT, LayoutT, DeviceT> get_subsampled_clouds(
  const int ncol, const int nlay, const int nbnd, const int ngpt,
  OpticalProps1sclK<RealT, LayoutT, DeviceT> &cloud_optics, GasOpticsRRTMGPK<RealT, LayoutT, DeviceT> &kdist,
  CldT &cld, PlayT &p_lay) {
  // Initialized subsampled optics
  OpticalProps1sclK<RealT, LayoutT, DeviceT> subsampled_optics;
  subsampled_optics.init(kdist.get_band_lims_wavenumber(), kdist.get_band_lims_gpoint(), "subsampled_optics");
  subsampled_optics.alloc_1scl(ncol, nlay);
  // Check that we do not have clouds with no optical properties; this would get corrected
  // when we assign optical props, but we want to use a "radiative cloud fraction"
  // for the subcolumn sampling too because otherwise we can get vertically-contiguous cloud
  // mask profiles with no actual cloud properties in between, which would just further overestimate
  // the vertical correlation of cloudy layers. I.e., cloudy layers might look maximally overlapped
  // even when separated by layers with no cloud properties, when in fact those layers should be
  // randomly overlapped.
  auto cldfrac_rad = view_t<RealT**>("cldfrac_rad", ncol, nlay);
  Kokkos::parallel_for(MDRP::template get<3>({nbnd,nlay,ncol}), KOKKOS_LAMBDA (int ibnd, int ilay, int icol) {
    if (cloud_optics.tau(icol,ilay,ibnd) > 0) {
      cldfrac_rad(icol,ilay) = cld(icol,ilay);
    }
  });
  // Get subcolumn cloud mask
  int overlap = 1;
  // Get unique seeds for each column that are reproducible across different MPI rank layouts;
  // use decimal part of pressure for this, consistent with the implementation in EAM; use different
  // seed values for longwave and shortwave
  auto seeds = view_t<int*>("seeds", ncol);
  Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    seeds(icol) = 1e9 * (p_lay(icol,nlay-2) - int(p_lay(icol,nlay-2)));
  });
  auto cldmask = get_subcolumn_mask(ncol, nlay, ngpt, cldfrac_rad, overlap, seeds);
  // Assign optical properties to subcolumns (note this implements MCICA)
  auto gpoint_bands = kdist.get_gpoint_bands();
  Kokkos::parallel_for(MDRP::template get<3>({ngpt,nlay,ncol}), KOKKOS_LAMBDA(int igpt, int ilay, int icol) {
      auto ibnd = gpoint_bands(igpt);
      if (cldmask(icol,ilay,igpt) == 1) {
        subsampled_optics.tau(icol,ilay,igpt) = cloud_optics.tau(icol,ilay,ibnd);
      } else {
        subsampled_optics.tau(icol,ilay,igpt) = 0;
      }
    });
  return subsampled_optics;
}

}; // struct rrtmgp_interface
#endif // RRTMGP_ENABLE_KOKKOS

} // namespace rrtmgp
} // namespace scream

#endif  // SCREAM_RRTMGP_INTERFACE_HPP
