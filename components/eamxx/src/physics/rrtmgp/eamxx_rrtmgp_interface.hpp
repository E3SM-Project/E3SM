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
  TIMED_KERNEL(yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay, ncol), YAKL_LAMBDA(int ilay, int icol) {
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
      }));
}

template<class S, class T>
void limit_to_bounds(S const &arr_in, T const lower, T const upper, S &arr_out) {
  TIMED_KERNEL(yakl::c::parallel_for(arr_in.totElems(), YAKL_LAMBDA(int i) {
    arr_out.data()[i] = std::min(std::max(arr_in.data()[i], lower), upper);
      }));
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
using hview_t = Kokkos::View<T, LayoutT, HostDevice>;

using pool_t = conv::MemPoolSingleton<RealT, LayoutT, DeviceT>;

using real1dk = view_t<RealT*>;
using real2dk = view_t<RealT**>;
using real3dk = view_t<RealT***>;
using creal1dk = view_t<const RealT*>;
using creal2dk = view_t<const RealT**>;
using creal3dk = view_t<const RealT***>;
using int1dk  = view_t<int*>;
using int2dk  = view_t<int**>;
using int3dk  = view_t<int***>;

using gas_optics_t   = GasOpticsRRTMGPK<RealT, LayoutT, DeviceT>;
using cloud_optics_t = CloudOpticsK<RealT, LayoutT, DeviceT>;
using gas_concs_t    = GasConcsK<RealT, LayoutT, DeviceT>;
using fluxes_t       = FluxesBybandK<RealT, LayoutT, DeviceT>;
using fluxes_broadband_t = FluxesBroadbandK<RealT, LayoutT, DeviceT>;
using optical_props_t = OpticalPropsK<RealT, LayoutT, DeviceT>;
using optical_props1_t = OpticalProps1sclK<RealT, LayoutT, DeviceT>;
using optical_props2_t = OpticalProps2strK<RealT, LayoutT, DeviceT>;
using source_func_t = SourceFuncLWK<RealT, LayoutT, DeviceT>;

/*
 * Objects containing k-distribution information need to be initialized
 * once and then persist throughout the life of the program, so we
 * declare them here within the rrtmgp namespace.
 */
static inline std::unique_ptr<gas_optics_t> k_dist_sw_k;
static inline std::unique_ptr<gas_optics_t> k_dist_lw_k;

/*
 * Objects containing cloud optical property look-up table information.
 * We want to initialize these once and use throughout the life of the
 * program, so declare here and read data in during rrtmgp_initialize().
 */
static inline std::unique_ptr<cloud_optics_t> cloud_optics_sw_k;
static inline std::unique_ptr<cloud_optics_t> cloud_optics_lw_k;

/*
 * Flag to indicate whether or not we have initialized RRTMGP
 */
static inline bool initialized_k = false;

/*
 * Initialize data for RRTMGP driver. Increase multiplier to allocate more pool space.
 */
static void rrtmgp_initialize(
  const gas_concs_t &gas_concs,
  const std::string& coefficients_file_sw, const std::string& coefficients_file_lw,
  const std::string& cloud_optics_file_sw, const std::string& cloud_optics_file_lw,
  const std::shared_ptr<spdlog::logger>& logger,
  const double multiplier = 1.0)
{
  // If we've already initialized, just exit
  if (initialized_k) {
    if (logger)
      logger->info("RRTMGP is already initialized; skipping\n");
    return;
  }

  // Initialize Kokkos
  if (!Kokkos::is_initialized()) {  Kokkos::initialize(); }

  // Create objects for static ptrs
  k_dist_sw_k = std::make_unique<gas_optics_t>();
  k_dist_lw_k = std::make_unique<gas_optics_t>();
  cloud_optics_sw_k = std::make_unique<cloud_optics_t>();
  cloud_optics_lw_k = std::make_unique<cloud_optics_t>();

  // Load and initialize absorption coefficient data
  load_and_init(*k_dist_sw_k, coefficients_file_sw, gas_concs);
  load_and_init(*k_dist_lw_k, coefficients_file_lw, gas_concs);

  // Load and initialize cloud optical property look-up table information
  load_cld_lutcoeff(*cloud_optics_sw_k, cloud_optics_file_sw);
  load_cld_lutcoeff(*cloud_optics_lw_k, cloud_optics_file_lw);

  // initialize kokkos rrtmgp pool allocator
  const size_t base_ref = 40000;
  const size_t ncol = gas_concs.ncol;
  const size_t nlay = gas_concs.nlay;
  const size_t nlev = SCREAM_NUM_VERTICAL_LEV;
  const size_t my_size_ref = ncol * nlay * nlev;
  pool_t::init(2e6 * (float(my_size_ref) / base_ref) * multiplier);

  // We are now initialized!
  initialized_k = true;
}

/*
 * Compute band-by-band surface albedos from broadband albedos.
 */
static void compute_band_by_band_surface_albedos(
  const int ncol, const int nswbands,
  const creal1dk &sfc_alb_dir_vis, const creal1dk &sfc_alb_dir_nir,
  const creal1dk &sfc_alb_dif_vis, const creal1dk &sfc_alb_dif_nir,
  const real2dk &sfc_alb_dir,     const real2dk &sfc_alb_dif)
{
  EKAT_ASSERT_MSG(initialized_k, "Error! rrtmgp_initialize must be called before GasOpticsRRTMGP object can be used.");
  auto wavenumber_limits = k_dist_sw_k->get_band_lims_wavenumber();

  EKAT_ASSERT_MSG(wavenumber_limits.extent(0) == 2,
                  "Error! 1st dimension for wavenumber_limits should be 2. It's " << wavenumber_limits.extent(0));
  EKAT_ASSERT_MSG(wavenumber_limits.extent(1) == static_cast<size_t>(nswbands),
                  "Error! 2nd dimension for wavenumber_limits should be " + std::to_string(nswbands) + " (nswbands).");

  // Loop over bands, and determine for each band whether it is broadly in the
  // visible or infrared part of the spectrum (visible or "not visible")
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nswbands, icol, ibnd,
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
  ));
}

/*
 * Compute broadband visible/UV and near-infrared surface fluxes.
 */
static void compute_broadband_surface_fluxes(
  const int ncol, const int ktop, const int nswbands,
  const real3dk &sw_bnd_flux_dir , const real3dk &sw_bnd_flux_dif ,
  const real1dk &sfc_flux_dir_vis, const real1dk &sfc_flux_dir_nir,
  const real1dk &sfc_flux_dif_vis, const real1dk &sfc_flux_dif_nir)
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
  auto wavenumber_limits = k_dist_sw_k->get_band_lims_wavenumber();
  TIMED_KERNEL(Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(const int icol) {
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
  }));
}

/*
 * Main driver code to run RRTMGP.
 * The input logger is in charge of outputing info to
 * screen and/or to file (or neither), depending on how it was set up.
 */
static void rrtmgp_main(
  const int ncol, const int nlay,
  const creal2dk &p_lay, const creal2dk &t_lay, const creal2dk &p_lev, const creal2dk &t_lev,
  gas_concs_t &gas_concs,
  const creal2dk &sfc_alb_dir, const creal2dk &sfc_alb_dif, const real1dk &mu0,
  const real2dk &lwp, const real2dk &iwp, const creal2dk &rel, const creal2dk &rei, const real2dk &cldfrac,
  const real3dk &aer_tau_sw, const real3dk &aer_ssa_sw, const real3dk &aer_asm_sw, const real3dk &aer_tau_lw,
  const real3dk &cld_tau_sw_bnd, const real3dk &cld_tau_lw_bnd,
  const real3dk &cld_tau_sw_gpt, const real3dk &cld_tau_lw_gpt,
  const real2dk &sw_flux_up, const real2dk &sw_flux_dn, const real2dk &sw_flux_dn_dir,
  const real2dk &lw_flux_up, const real2dk &lw_flux_dn,
  const real2dk &sw_clnclrsky_flux_up, const real2dk &sw_clnclrsky_flux_dn, const real2dk &sw_clnclrsky_flux_dn_dir,
  const real2dk &sw_clrsky_flux_up, const real2dk &sw_clrsky_flux_dn, const real2dk &sw_clrsky_flux_dn_dir,
  const real2dk &sw_clnsky_flux_up, const real2dk &sw_clnsky_flux_dn, const real2dk &sw_clnsky_flux_dn_dir,
  const real2dk &lw_clnclrsky_flux_up, const real2dk &lw_clnclrsky_flux_dn,
  const real2dk &lw_clrsky_flux_up, const real2dk &lw_clrsky_flux_dn,
  const real2dk &lw_clnsky_flux_up, const real2dk &lw_clnsky_flux_dn,
  const real3dk &sw_bnd_flux_up, const real3dk &sw_bnd_flux_dn, const real3dk &sw_bnd_flux_dn_dir,
  const real3dk &lw_bnd_flux_up, const real3dk &lw_bnd_flux_dn,
  const Real tsi_scaling,
  const std::shared_ptr<spdlog::logger>& logger,
  const bool extra_clnclrsky_diag = false, const bool extra_clnsky_diag = false)
{
  const int sw_nband = k_dist_sw_k->get_nband();
  const int lw_nband = k_dist_lw_k->get_nband();
  const int sw_ngpt = k_dist_sw_k->get_ngpt();
  const int lw_ngpt = k_dist_lw_k->get_ngpt();

#ifdef SCREAM_RRTMGP_DEBUG
  // Sanity check inputs, and possibly repair
  check_range_k(t_lay      ,  k_dist_sw_k->get_temp_min(),         k_dist_sw_k->get_temp_max(), "rrtmgp_main::t_lay");
  check_range_k(t_lev      ,  k_dist_sw_k->get_temp_min(),         k_dist_sw_k->get_temp_max(), "rrtmgp_main::t_lev");
  check_range_k(p_lay      , k_dist_sw_k->get_press_min(),        k_dist_sw_k->get_press_max(), "rrtmgp_main::p_lay");
  check_range_k(p_lev      , k_dist_sw_k->get_press_min(),        k_dist_sw_k->get_press_max(), "rrtmgp_main::p_lev");
  check_range_k(sfc_alb_dir,                         0,                                1, "rrtmgp_main::sfc_alb_dir");
  check_range_k(sfc_alb_dif,                         0,                                1, "rrtmgp_main::sfc_alb_dif");
  check_range_k(mu0        ,                         0,                                1, "rrtmgp_main::mu0");
  check_range_k(lwp        ,                         0, std::numeric_limits<RealT>::max(), "rrtmgp_main::lwp");
  check_range_k(iwp        ,                         0, std::numeric_limits<RealT>::max(), "rrtmgp_main::iwp");
  check_range_k(rel        ,                         0, std::numeric_limits<RealT>::max(), "rrtmgp_main::rel");
  check_range_k(rei        ,                         0, std::numeric_limits<RealT>::max(), "rrtmgp_main::rei");
#endif

  auto sw_band2gpt_mem = pool_t::template alloc<int>(2, sw_nband);
  auto sw_gpt2band_mem = pool_t::template alloc<int>(   sw_nband);
  auto lw_band2gpt_mem = pool_t::template alloc<int>(2, lw_nband);
  auto lw_gpt2band_mem = pool_t::template alloc<int>(   lw_nband);

  auto sw_cloud_band2gpt_mem = pool_t::template alloc<int>(2, sw_nband);
  auto sw_cloud_gpt2band_mem = pool_t::template alloc<int>(   sw_nband);
  auto lw_cloud_band2gpt_mem = pool_t::template alloc<int>(2, lw_nband);
  auto lw_cloud_gpt2band_mem = pool_t::template alloc<int>(   lw_nband);

  auto sw_subcloud_band2gpt_mem = pool_t::template alloc<int>(2, sw_nband);
  auto sw_subcloud_gpt2band_mem = pool_t::template alloc<int>(   sw_ngpt);
  auto lw_subcloud_band2gpt_mem = pool_t::template alloc<int>(2, lw_nband);
  auto lw_subcloud_gpt2band_mem = pool_t::template alloc<int>(   lw_ngpt);

  auto sw_tau_mem = pool_t::template alloc<RealT>(ncol, nlay, sw_nband);
  auto sw_ssa_mem = pool_t::template alloc<RealT>(ncol, nlay, sw_nband);
  auto sw_g_mem   = pool_t::template alloc<RealT>(ncol, nlay, sw_nband);
  auto lw_tau_mem = pool_t::template alloc<RealT>(ncol, nlay, lw_nband);

  auto sw_cloud_tau_mem = pool_t::template alloc<RealT>(ncol, nlay, sw_nband);
  auto sw_cloud_ssa_mem = pool_t::template alloc<RealT>(ncol, nlay, sw_nband);
  auto sw_cloud_g_mem   = pool_t::template alloc<RealT>(ncol, nlay, sw_nband);
  auto lw_cloud_tau_mem = pool_t::template alloc<RealT>(ncol, nlay, lw_nband);

  auto sw_subcloud_tau_mem = pool_t::template alloc<RealT>(ncol, nlay, sw_ngpt);
  auto sw_subcloud_ssa_mem = pool_t::template alloc<RealT>(ncol, nlay, sw_ngpt);
  auto sw_subcloud_g_mem   = pool_t::template alloc<RealT>(ncol, nlay, sw_ngpt);
  auto lw_subcloud_tau_mem = pool_t::template alloc<RealT>(ncol, nlay, lw_ngpt);

  // Setup pointers to RRTMGP SW fluxes
  fluxes_t fluxes_sw;
  fluxes_sw.flux_up = sw_flux_up;
  fluxes_sw.flux_dn = sw_flux_dn;
  fluxes_sw.flux_dn_dir = sw_flux_dn_dir;
  fluxes_sw.bnd_flux_up = sw_bnd_flux_up;
  fluxes_sw.bnd_flux_dn = sw_bnd_flux_dn;
  fluxes_sw.bnd_flux_dn_dir = sw_bnd_flux_dn_dir;
  // Clean-clear-sky
  fluxes_broadband_t clnclrsky_fluxes_sw;
  clnclrsky_fluxes_sw.flux_up = sw_clnclrsky_flux_up;
  clnclrsky_fluxes_sw.flux_dn = sw_clnclrsky_flux_dn;
  clnclrsky_fluxes_sw.flux_dn_dir = sw_clnclrsky_flux_dn_dir;
  // Clear-sky
  fluxes_broadband_t clrsky_fluxes_sw;
  clrsky_fluxes_sw.flux_up = sw_clrsky_flux_up;
  clrsky_fluxes_sw.flux_dn = sw_clrsky_flux_dn;
  clrsky_fluxes_sw.flux_dn_dir = sw_clrsky_flux_dn_dir;
  // Clean-sky
  fluxes_broadband_t clnsky_fluxes_sw;
  clnsky_fluxes_sw.flux_up = sw_clnsky_flux_up;
  clnsky_fluxes_sw.flux_dn = sw_clnsky_flux_dn;
  clnsky_fluxes_sw.flux_dn_dir = sw_clnsky_flux_dn_dir;

  // Setup pointers to RRTMGP LW fluxes
  fluxes_t fluxes_lw;
  fluxes_lw.flux_up = lw_flux_up;
  fluxes_lw.flux_dn = lw_flux_dn;
  fluxes_lw.bnd_flux_up = lw_bnd_flux_up;
  fluxes_lw.bnd_flux_dn = lw_bnd_flux_dn;
  // Clean-clear-sky
  fluxes_broadband_t clnclrsky_fluxes_lw;
  clnclrsky_fluxes_lw.flux_up = lw_clnclrsky_flux_up;
  clnclrsky_fluxes_lw.flux_dn = lw_clnclrsky_flux_dn;
  // Clear-sky
  fluxes_broadband_t clrsky_fluxes_lw;
  clrsky_fluxes_lw.flux_up = lw_clrsky_flux_up;
  clrsky_fluxes_lw.flux_dn = lw_clrsky_flux_dn;
  // Clean-sky
  fluxes_broadband_t clnsky_fluxes_lw;
  clnsky_fluxes_lw.flux_up = lw_clnsky_flux_up;
  clnsky_fluxes_lw.flux_dn = lw_clnsky_flux_dn;

  auto nswbands = k_dist_sw_k->get_nband();
  auto nlwbands = k_dist_lw_k->get_nband();

  // Setup aerosol optical properties
  optical_props2_t aerosol_sw;
  optical_props1_t aerosol_lw;
  aerosol_sw.init_no_alloc(k_dist_sw_k->get_band_lims_wavenumber(), sw_band2gpt_mem, sw_gpt2band_mem);
  aerosol_sw.alloc_2str_no_alloc(ncol, nlay, sw_tau_mem, sw_ssa_mem, sw_g_mem);
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol,nlay,nswbands, icol, ilay, ibnd,
    aerosol_sw.tau(icol,ilay,ibnd) = aer_tau_sw(icol,ilay,ibnd);
    aerosol_sw.ssa(icol,ilay,ibnd) = aer_ssa_sw(icol,ilay,ibnd);
    aerosol_sw.g  (icol,ilay,ibnd) = aer_asm_sw(icol,ilay,ibnd);
  ));
  aerosol_lw.init_no_alloc(k_dist_lw_k->get_band_lims_wavenumber(), lw_band2gpt_mem, lw_gpt2band_mem);
  aerosol_lw.alloc_1scl_no_alloc(ncol, nlay, lw_tau_mem);
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol,nlay,nlwbands, icol, ilay, ibnd,
    aerosol_lw.tau(icol,ilay,ibnd) = aer_tau_lw(icol,ilay,ibnd);
  ));

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
  optical_props2_t clouds_sw = get_cloud_optics_sw(ncol, nlay, *cloud_optics_sw_k, *k_dist_sw_k, lwp, iwp, rel, rei, sw_cloud_band2gpt_mem, sw_cloud_gpt2band_mem, sw_cloud_tau_mem, sw_cloud_ssa_mem, sw_cloud_g_mem);
  optical_props1_t clouds_lw = get_cloud_optics_lw(ncol, nlay, *cloud_optics_lw_k, *k_dist_lw_k, lwp, iwp, rel, rei, lw_cloud_band2gpt_mem, lw_cloud_gpt2band_mem, lw_cloud_tau_mem);
  Kokkos::deep_copy(cld_tau_sw_bnd, clouds_sw.tau);
  Kokkos::deep_copy(cld_tau_lw_bnd, clouds_lw.tau);

  // Do subcolumn sampling to map bands -> gpoints based on cloud fraction and overlap assumption;
  // This implements the Monte Carlo Independing Column Approximation by mapping only a single
  // subcolumn (cloud state) to each gpoint.
  auto nswgpts = k_dist_sw_k->get_ngpt();
  auto clouds_sw_gpt = get_subsampled_clouds(ncol, nlay, nswbands, nswgpts, clouds_sw, *k_dist_sw_k, cldfrac, p_lay, sw_subcloud_band2gpt_mem, sw_subcloud_gpt2band_mem, sw_subcloud_tau_mem, sw_subcloud_ssa_mem, sw_subcloud_g_mem);

  // Longwave
  auto nlwgpts = k_dist_lw_k->get_ngpt();
  auto clouds_lw_gpt = get_subsampled_clouds(ncol, nlay, nlwbands, nlwgpts, clouds_lw, *k_dist_lw_k, cldfrac, p_lay, lw_subcloud_band2gpt_mem, lw_subcloud_gpt2band_mem, lw_subcloud_tau_mem);

  // Copy cloud properties to outputs (is this needed, or can we just use pointers?)
  // Alternatively, just compute and output a subcolumn cloud mask
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol, nlay, nswgpts, icol, ilay, igpt,
    cld_tau_sw_gpt(icol,ilay,igpt) = clouds_sw_gpt.tau(icol,ilay,igpt);
  ));
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol, nlay, nlwgpts, icol, ilay, igpt,
    cld_tau_lw_gpt(icol,ilay,igpt) = clouds_lw_gpt.tau(icol,ilay,igpt);
  ));

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
    *k_dist_sw_k, p_lay, t_lay, p_lev, t_lev, gas_concs,
    sfc_alb_dir, sfc_alb_dif, mu0, aerosol_sw, clouds_sw_gpt,
    fluxes_sw, clnclrsky_fluxes_sw, clrsky_fluxes_sw, clnsky_fluxes_sw,
    tsi_scaling, logger,
    extra_clnclrsky_diag, extra_clnsky_diag
            );

  // Do longwave
  rrtmgp_lw(
    ncol, nlay,
    *k_dist_lw_k, p_lay, t_lay, p_lev, t_lev, gas_concs,
    aerosol_lw, clouds_lw_gpt,
    fluxes_lw, clnclrsky_fluxes_lw, clrsky_fluxes_lw, clnsky_fluxes_lw,
    extra_clnclrsky_diag, extra_clnsky_diag
            );

  pool_t::dealloc(sw_band2gpt_mem);
  pool_t::dealloc(sw_gpt2band_mem);
  pool_t::dealloc(lw_band2gpt_mem);
  pool_t::dealloc(lw_gpt2band_mem);

  pool_t::dealloc(sw_cloud_band2gpt_mem);
  pool_t::dealloc(sw_cloud_gpt2band_mem);
  pool_t::dealloc(lw_cloud_band2gpt_mem);
  pool_t::dealloc(lw_cloud_gpt2band_mem);

  pool_t::dealloc(sw_subcloud_band2gpt_mem);
  pool_t::dealloc(sw_subcloud_gpt2band_mem);
  pool_t::dealloc(lw_subcloud_band2gpt_mem);
  pool_t::dealloc(lw_subcloud_gpt2band_mem);

  pool_t::dealloc(sw_tau_mem);
  pool_t::dealloc(sw_ssa_mem);
  pool_t::dealloc(sw_g_mem);
  pool_t::dealloc(lw_tau_mem);

  pool_t::dealloc(sw_cloud_tau_mem);
  pool_t::dealloc(sw_cloud_ssa_mem);
  pool_t::dealloc(sw_cloud_g_mem);
  pool_t::dealloc(lw_cloud_tau_mem);

  pool_t::dealloc(sw_subcloud_tau_mem);
  pool_t::dealloc(sw_subcloud_ssa_mem);
  pool_t::dealloc(sw_subcloud_g_mem);
  pool_t::dealloc(lw_subcloud_tau_mem);
}

/*
 * Perform any clean-up tasks
 */
static void rrtmgp_finalize()
{
  initialized_k = false;
  k_dist_sw_k->finalize();
  k_dist_lw_k->finalize();
  cloud_optics_sw_k->finalize(); //~CloudOptics();
  cloud_optics_lw_k->finalize(); //~CloudOptics();
  k_dist_sw_k = nullptr;
  k_dist_lw_k = nullptr;
  cloud_optics_sw_k = nullptr;
  cloud_optics_lw_k = nullptr;
  pool_t::finalize();
}

/*
 * Shortwave driver (called by rrtmgp_main)
 */
static void rrtmgp_sw(
  const int ncol, const int nlay,
  gas_optics_t &k_dist,
  const creal2dk &p_lay, const creal2dk &t_lay, const creal2dk &p_lev, const creal2dk &t_lev,
  gas_concs_t &gas_concs,
  const creal2dk &sfc_alb_dir, const creal2dk &sfc_alb_dif, const real1dk &mu0,
  optical_props2_t &aerosol, optical_props2_t &clouds,
  fluxes_t &fluxes, fluxes_broadband_t &clnclrsky_fluxes, fluxes_broadband_t &clrsky_fluxes, fluxes_broadband_t &clnsky_fluxes,
  const Real tsi_scaling,
  const std::shared_ptr<spdlog::logger>& logger,
  const bool extra_clnclrsky_diag, const bool extra_clnsky_diag)
{
  // Get problem sizes
  const int nbnd = k_dist.get_nband();
  const int ngpt = k_dist.get_ngpt();
  const int ngas = gas_concs.get_num_gases();

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
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nlay+1, icol, ilev,
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
  ));
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol, nlay+1, nbnd, icol, ilev, ibnd,
    bnd_flux_up    (icol,ilev,ibnd) = 0;
    bnd_flux_dn    (icol,ilev,ibnd) = 0;
    bnd_flux_dn_dir(icol,ilev,ibnd) = 0;
  ));

  // Get daytime indices
  auto dayIndices = pool_t::template alloc<int>(ncol);
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
    pool_t::dealloc(dayIndices);
    return;
  }

  // Allocate temporaries from pool

  auto mu0_day = pool_t::template alloc<RealT>(nday);

  auto p_lay_day = pool_t::template alloc<RealT>(nday, nlay);
  auto t_lay_day = pool_t::template alloc<RealT>(nday, nlay);
  auto vmr_day = pool_t::template alloc<RealT>(nday, nlay);
  auto t_lay_limited = pool_t::template alloc<RealT>(nday, nlay);

  auto p_lev_day = pool_t::template alloc<RealT>(nday, nlay+1);
  auto t_lev_day = pool_t::template alloc<RealT>(nday, nlay+1);
  auto flux_up_day = pool_t::template alloc<RealT>(nday, nlay+1);
  auto flux_dn_day = pool_t::template alloc<RealT>(nday, nlay+1);
  auto flux_dn_dir_day = pool_t::template alloc<RealT>(nday, nlay+1);

  auto vmr = pool_t::template alloc<RealT>(ncol, nlay);

  auto sfc_alb_dir_T = pool_t::template alloc<RealT>(nbnd, nday);
  auto sfc_alb_dif_T = pool_t::template alloc<RealT>(nbnd, nday);

  auto toa_flux = pool_t::template alloc<RealT>(nday, ngpt);

  auto bnd_flux_up_day = pool_t::template alloc<RealT>(nday, nlay+1, nbnd);
  auto bnd_flux_dn_day = pool_t::template alloc<RealT>(nday, nlay+1, nbnd);
  auto bnd_flux_dn_dir_day = pool_t::template alloc<RealT>(nday, nlay+1, nbnd);

  auto col_gas = pool_t::template alloc<RealT>(ncol, nlay, k_dist.get_ngas()+1);

  auto concs_mem = pool_t::template alloc<RealT>(nday, nlay, ngas);

  auto sw_aero_tau_mem = pool_t::template alloc<RealT>(nday, nlay, nbnd);
  auto sw_aero_ssa_mem = pool_t::template alloc<RealT>(nday, nlay, nbnd);
  auto sw_aero_g_mem = pool_t::template alloc<RealT>(nday, nlay, nbnd);

  auto sw_cloud_tau_mem = pool_t::template alloc<RealT>(nday, nlay, ngpt);
  auto sw_cloud_ssa_mem = pool_t::template alloc<RealT>(nday, nlay, ngpt);
  auto sw_cloud_g_mem = pool_t::template alloc<RealT>(nday, nlay, ngpt);
  auto sw_optics_tau_mem = pool_t::template alloc<RealT>(nday, nlay, ngpt);
  auto sw_optics_ssa_mem = pool_t::template alloc<RealT>(nday, nlay, ngpt);
  auto sw_optics_g_mem = pool_t::template alloc<RealT>(nday, nlay, ngpt);
  auto sw_noaero_tau_mem = pool_t::template alloc<RealT>(nday, nlay, ngpt);
  auto sw_noaero_ssa_mem = pool_t::template alloc<RealT>(nday, nlay, ngpt);
  auto sw_noaero_g_mem = pool_t::template alloc<RealT>(nday, nlay, ngpt);

  auto sw_aero_band2gpt_mem = pool_t::template alloc<int>(2, nbnd);
  auto sw_aero_gpt2band_mem = pool_t::template alloc<int>(   nbnd);
  auto sw_cloud_band2gpt_mem = pool_t::template alloc<int>(2, nbnd);
  auto sw_cloud_gpt2band_mem = pool_t::template alloc<int>(   ngpt);
  auto sw_optics_band2gpt_mem = pool_t::template alloc<int>(2, nbnd);
  auto sw_optics_gpt2band_mem = pool_t::template alloc<int>(   ngpt);
  auto sw_noaero_band2gpt_mem = pool_t::template alloc<int>(2, nbnd);
  auto sw_noaero_gpt2band_mem = pool_t::template alloc<int>(   ngpt);

  // Subset mu0
  TIMED_KERNEL(Kokkos::parallel_for(nday, KOKKOS_LAMBDA(int iday) {
    mu0_day(iday) = mu0(dayIndices(iday));
  }));

  // subset state variables
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(nday, nlay, iday, ilay,
    p_lay_day(iday,ilay) = p_lay(dayIndices(iday),ilay);
    t_lay_day(iday,ilay) = t_lay(dayIndices(iday),ilay);
  ));
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(nday, nlay+1, iday, ilev,
    p_lev_day(iday,ilev) = p_lev(dayIndices(iday),ilev);
    t_lev_day(iday,ilev) = t_lev(dayIndices(iday),ilev);
  ));

  // Subset gases
  auto gas_names = gas_concs.get_gas_names();
  gas_concs_t gas_concs_day;
  gas_concs_day.init_no_alloc(gas_names, nday, nlay, concs_mem);
  for (int igas = 0; igas < ngas; igas++) {
    gas_concs.get_vmr(gas_names[igas], vmr);
    TIMED_KERNEL(FLATTEN_MD_KERNEL2(nday, nlay, iday, ilay,
      vmr_day(iday,ilay) = vmr(dayIndices(iday),ilay);
    ));
    gas_concs_day.set_vmr(gas_names[igas], vmr_day);
  }

  // Subset aerosol optics
  optical_props2_t aerosol_day;
  aerosol_day.init_no_alloc(k_dist.get_band_lims_wavenumber(), sw_aero_band2gpt_mem, sw_aero_gpt2band_mem);
  aerosol_day.alloc_2str_no_alloc(nday, nlay, sw_aero_tau_mem, sw_aero_ssa_mem, sw_aero_g_mem);
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(nday,nlay,nbnd, iday, ilay, ibnd,
    aerosol_day.tau(iday,ilay,ibnd) = aerosol.tau(dayIndices(iday),ilay,ibnd);
    aerosol_day.ssa(iday,ilay,ibnd) = aerosol.ssa(dayIndices(iday),ilay,ibnd);
    aerosol_day.g  (iday,ilay,ibnd) = aerosol.g  (dayIndices(iday),ilay,ibnd);
  ));

  // Subset cloud optics
  // TODO: nbnd -> ngpt once we pass sub-sampled cloud state
  optical_props2_t clouds_day;
  clouds_day.init_no_alloc(k_dist.get_band_lims_wavenumber(), k_dist.get_band_lims_gpoint(), sw_cloud_band2gpt_mem, sw_cloud_gpt2band_mem);
  clouds_day.alloc_2str_no_alloc(nday, nlay, sw_cloud_tau_mem, sw_cloud_ssa_mem, sw_cloud_g_mem);
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(nday,nlay,ngpt, iday, ilay, igpt,
    clouds_day.tau(iday,ilay,igpt) = clouds.tau(dayIndices(iday),ilay,igpt);
    clouds_day.ssa(iday,ilay,igpt) = clouds.ssa(dayIndices(iday),ilay,igpt);
    clouds_day.g  (iday,ilay,igpt) = clouds.g  (dayIndices(iday),ilay,igpt);
  ));

  // RRTMGP assumes surface albedos have a screwy dimension ordering
  // for some strange reason, so we need to transpose these; also do
  // daytime subsetting in the same kernel
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(nbnd,nday, ibnd, icol,
    sfc_alb_dir_T(ibnd,icol) = sfc_alb_dir(dayIndices(icol),ibnd);
    sfc_alb_dif_T(ibnd,icol) = sfc_alb_dif(dayIndices(icol),ibnd);
  ));

  // Temporaries we need for daytime-only fluxes
  fluxes_t fluxes_day;
  fluxes_day.flux_up         = flux_up_day;
  fluxes_day.flux_dn         = flux_dn_day;
  fluxes_day.flux_dn_dir     = flux_dn_dir_day;
  fluxes_day.bnd_flux_up     = bnd_flux_up_day;
  fluxes_day.bnd_flux_dn     = bnd_flux_dn_day;
  fluxes_day.bnd_flux_dn_dir = bnd_flux_dn_dir_day;

  // Allocate space for optical properties
  optical_props2_t optics;
  optics.alloc_2str_no_alloc(nday, nlay, k_dist, sw_optics_band2gpt_mem, sw_optics_gpt2band_mem, sw_optics_tau_mem, sw_optics_ssa_mem, sw_optics_g_mem);

  optical_props2_t optics_no_aerosols;
  if (extra_clnsky_diag) {
    // Allocate space for optical properties (no aerosols)
    optics_no_aerosols.alloc_2str_no_alloc(nday, nlay, k_dist, sw_noaero_band2gpt_mem, sw_noaero_gpt2band_mem, sw_noaero_tau_mem, sw_noaero_ssa_mem, sw_noaero_g_mem);
  }

  // Limit temperatures for gas optics look-up tables
  limit_to_bounds_k(t_lay_day, k_dist_sw_k->get_temp_min(), k_dist_sw_k->get_temp_max(), t_lay_limited);

  // Do gas optics
  bool top_at_1 = false;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(int, bool& val) {
    val |= p_lay(0, 0) < p_lay(0, nlay-1);
  }, Kokkos::LOr<bool>(top_at_1));

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
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(nday, ngpt, iday, igpt,
    toa_flux(iday,igpt) = tsi_scaling * toa_flux(iday,igpt);
  ));

  if (extra_clnclrsky_diag) {
    // Compute clear-clean-sky (just gas) fluxes on daytime columns
    rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);
    // Expand daytime fluxes to all columns
    TIMED_KERNEL(FLATTEN_MD_KERNEL2(nday, nlay+1, iday, ilev,
      const int icol = dayIndices(iday);
      clnclrsky_flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
      clnclrsky_flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
      clnclrsky_flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
    ));
  }

  // Combine gas and aerosol optics
  aerosol_day.delta_scale();
  aerosol_day.increment(optics);

  // Compute clearsky (gas + aerosol) fluxes on daytime columns
  rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);

  // Expand daytime fluxes to all columns
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(nday, nlay+1, iday, ilev,
    const int icol = dayIndices(iday);
    clrsky_flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
    clrsky_flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
    clrsky_flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
  ));

  // Now merge in cloud optics and do allsky calculations

  // Combine gas and cloud optics
  clouds_day.delta_scale();
  clouds_day.increment(optics);
  // Compute fluxes on daytime columns
  rte_sw(optics, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);
  // Expand daytime fluxes to all columns
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(nday, nlay+1, iday, ilev,
    const int icol = dayIndices(iday);
    flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
    flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
    flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
  ));
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(nday, nlay+1, nbnd, iday, ilev, ibnd,
    const int icol = dayIndices(iday);
    bnd_flux_up    (icol,ilev,ibnd) = bnd_flux_up_day    (iday,ilev,ibnd);
    bnd_flux_dn    (icol,ilev,ibnd) = bnd_flux_dn_day    (iday,ilev,ibnd);
    bnd_flux_dn_dir(icol,ilev,ibnd) = bnd_flux_dn_dir_day(iday,ilev,ibnd);
  ));

  if (extra_clnsky_diag) {
    // First increment clouds in optics_no_aerosols
    clouds_day.increment(optics_no_aerosols);
    // Compute cleansky (gas + clouds) fluxes on daytime columns
    rte_sw(optics_no_aerosols, top_at_1, mu0_day, toa_flux, sfc_alb_dir_T, sfc_alb_dif_T, fluxes_day);
    // Expand daytime fluxes to all columns
    TIMED_KERNEL(FLATTEN_MD_KERNEL2(nday, nlay+1, iday, ilev,
      const int icol = dayIndices(iday);
      clnsky_flux_up    (icol,ilev) = flux_up_day    (iday,ilev);
      clnsky_flux_dn    (icol,ilev) = flux_dn_day    (iday,ilev);
      clnsky_flux_dn_dir(icol,ilev) = flux_dn_dir_day(iday,ilev);
    ));
  }

  pool_t::dealloc(dayIndices);

  pool_t::dealloc(mu0_day);

  pool_t::dealloc(p_lay_day);
  pool_t::dealloc(t_lay_day);
  pool_t::dealloc(vmr_day);
  pool_t::dealloc(t_lay_limited);

  pool_t::dealloc(p_lev_day);
  pool_t::dealloc(t_lev_day);
  pool_t::dealloc(flux_up_day);
  pool_t::dealloc(flux_dn_day);
  pool_t::dealloc(flux_dn_dir_day);

  pool_t::dealloc(vmr);

  pool_t::dealloc(sfc_alb_dir_T);
  pool_t::dealloc(sfc_alb_dif_T);

  pool_t::dealloc(toa_flux);

  pool_t::dealloc(bnd_flux_up_day);
  pool_t::dealloc(bnd_flux_dn_day);
  pool_t::dealloc(bnd_flux_dn_dir_day);

  pool_t::dealloc(col_gas);

  pool_t::dealloc(concs_mem);

  pool_t::dealloc(sw_aero_tau_mem);
  pool_t::dealloc(sw_aero_ssa_mem);
  pool_t::dealloc(sw_aero_g_mem);

  pool_t::dealloc(sw_cloud_tau_mem);
  pool_t::dealloc(sw_cloud_ssa_mem);
  pool_t::dealloc(sw_cloud_g_mem);
  pool_t::dealloc(sw_optics_tau_mem);
  pool_t::dealloc(sw_optics_ssa_mem);
  pool_t::dealloc(sw_optics_g_mem);
  pool_t::dealloc(sw_noaero_tau_mem);
  pool_t::dealloc(sw_noaero_ssa_mem);
  pool_t::dealloc(sw_noaero_g_mem);

  pool_t::dealloc(sw_aero_band2gpt_mem);
  pool_t::dealloc(sw_aero_gpt2band_mem);
  pool_t::dealloc(sw_cloud_band2gpt_mem);
  pool_t::dealloc(sw_cloud_gpt2band_mem);
  pool_t::dealloc(sw_optics_band2gpt_mem);
  pool_t::dealloc(sw_optics_gpt2band_mem);
  pool_t::dealloc(sw_noaero_band2gpt_mem);
  pool_t::dealloc(sw_noaero_gpt2band_mem);
}

/*
 * Longwave driver (called by rrtmgp_main)
 */
static void rrtmgp_lw(
  const int ncol, const int nlay,
  gas_optics_t &k_dist,
  const creal2dk &p_lay, const creal2dk &t_lay, const creal2dk &p_lev, const creal2dk &t_lev,
  gas_concs_t &gas_concs,
  optical_props1_t &aerosol, optical_props1_t &clouds,
  fluxes_t &fluxes, fluxes_broadband_t &clnclrsky_fluxes, fluxes_broadband_t &clrsky_fluxes, fluxes_broadband_t &clnsky_fluxes,
  const bool extra_clnclrsky_diag, const bool extra_clnsky_diag)
{
  // Problem size
  int nbnd = k_dist.get_nband();
  int constexpr max_gauss_pts = 4;
  const int ngpt = k_dist.get_ngpt();

  auto t_sfc               = pool_t::template alloc<RealT>(ncol);
  auto emis_sfc            = pool_t::template alloc<RealT>(nbnd, ncol);
  auto gauss_Ds            = pool_t::template alloc<RealT>(max_gauss_pts, max_gauss_pts);
  auto gauss_wts           = pool_t::template alloc<RealT>(max_gauss_pts, max_gauss_pts);
  auto t_lay_limited       = pool_t::template alloc<RealT>(ncol, nlay);
  auto t_lev_limited       = pool_t::template alloc<RealT>(ncol, nlay+1);
  auto col_gas             = pool_t::template alloc<RealT>(ncol, nlay, k_dist.get_ngas()+1);
  auto lw_optics_tau_mem   = pool_t::template alloc<RealT>(ncol, nlay, ngpt);
  auto lw_noaero_tau_mem   = pool_t::template alloc<RealT>(ncol, nlay, ngpt);
  auto lay_source_mem      = pool_t::template alloc<RealT>(ncol, nlay, ngpt);
  auto lev_source_inc_mem  = pool_t::template alloc<RealT>(ncol, nlay, ngpt);
  auto lev_source_dec_mem  = pool_t::template alloc<RealT>(ncol, nlay, ngpt);
  auto sfc_source_mem      = pool_t::template alloc<RealT>(ncol, ngpt);

  auto lw_optics_band2gpt_mem = pool_t::template alloc<int>(2, nbnd);
  auto lw_optics_gpt2band_mem = pool_t::template alloc<int>(   ngpt);
  auto lw_noaero_band2gpt_mem = pool_t::template alloc<int>(2, nbnd);
  auto lw_noaero_gpt2band_mem = pool_t::template alloc<int>(   ngpt);
  auto lw_source_band2gpt_mem = pool_t::template alloc<int>(2, nbnd);
  auto lw_source_gpt2band_mem = pool_t::template alloc<int>(   ngpt);

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
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nlay + 1, icol, ilev,
      flux_up(icol, ilev)           = 0;
      flux_dn(icol, ilev)           = 0;
      clnclrsky_flux_up(icol, ilev) = 0;
      clnclrsky_flux_dn(icol, ilev) = 0;
      clrsky_flux_up(icol, ilev)    = 0;
      clrsky_flux_dn(icol, ilev)    = 0;
      clnsky_flux_up(icol, ilev)    = 0;
      clnsky_flux_dn(icol, ilev)    = 0;
  ));
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol, nlay + 1, nbnd, icol, ilev, ibnd,
      bnd_flux_up(icol, ilev, ibnd) = 0;
      bnd_flux_dn(icol, ilev, ibnd) = 0;
  ));

  // Allocate space for optical properties
  optical_props1_t optics;
  optics.alloc_1scl_no_alloc(ncol, nlay, k_dist, lw_optics_band2gpt_mem, lw_optics_gpt2band_mem, lw_optics_tau_mem);

  optical_props1_t optics_no_aerosols;
  if (extra_clnsky_diag) {
    // Allocate space for optical properties (no aerosols)
    optics_no_aerosols.alloc_1scl_no_alloc(ncol, nlay, k_dist, lw_noaero_band2gpt_mem, lw_noaero_gpt2band_mem, lw_noaero_tau_mem);
  }

  // Boundary conditions
  source_func_t lw_sources;
  lw_sources.alloc_no_alloc(ncol, nlay, k_dist, lw_source_band2gpt_mem, lw_source_gpt2band_mem, sfc_source_mem, lay_source_mem, lev_source_inc_mem, lev_source_dec_mem);

  bool top_at_1 = false;
  Kokkos::parallel_reduce(1, KOKKOS_LAMBDA(int, bool& val) {
    val |= p_lay(0, 0) < p_lay(0, nlay-1);
  }, Kokkos::LOr<bool>(top_at_1));

  // Surface temperature
  TIMED_KERNEL(Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    t_sfc(icol) = t_lev(icol, conv::merge(nlay, 0, top_at_1));
  }));
  Kokkos::deep_copy(emis_sfc , 0.98);

  // Get Gaussian quadrature weights
  // TODO: move this crap out of userland!
  // Weights and angle secants for first order (k=1) Gaussian quadrature.
  //   Values from Table 2, Clough et al, 1992, doi:10.1029/92JD01419
  //   after Abramowitz & Stegun 1972, page 921
  RealT gauss_Ds_host_raw[max_gauss_pts][max_gauss_pts] = {
    {1.66, 1.18350343, 1.09719858, 1.06056257},
    {0., 2.81649655, 1.69338507, 1.38282560},
    {0., 0., 4.70941630, 2.40148179},
    {0., 0., 0., 7.15513024}
  };
  hview_t<RealT**> gauss_Ds_host (&gauss_Ds_host_raw[0][0], max_gauss_pts, max_gauss_pts);

  RealT gauss_wts_host_raw[max_gauss_pts][max_gauss_pts] = {
    {0.5, 0.3180413817, 0.2009319137, 0.1355069134},
    {0., 0.1819586183, 0.2292411064, 0.2034645680},
    {0., 0., 0.0698269799, 0.1298475476},
    {0., 0., 0., 0.0311809710}
  };

  hview_t<RealT**> gauss_wts_host(&gauss_wts_host_raw[0][0],max_gauss_pts,max_gauss_pts);

  Kokkos::deep_copy(gauss_Ds,  gauss_Ds_host);
  Kokkos::deep_copy(gauss_wts, gauss_wts_host);

  // Limit temperatures for gas optics look-up tables
  limit_to_bounds_k(t_lay, k_dist_lw_k->get_temp_min(), k_dist_lw_k->get_temp_max(), t_lay_limited);
  limit_to_bounds_k(t_lev, k_dist_lw_k->get_temp_min(), k_dist_lw_k->get_temp_max(), t_lev_limited);

  // Do gas optics
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

  pool_t::dealloc(t_sfc);
  pool_t::dealloc(emis_sfc);
  pool_t::dealloc(gauss_Ds);
  pool_t::dealloc(gauss_wts);
  pool_t::dealloc(t_lay_limited);
  pool_t::dealloc(t_lev_limited);
  pool_t::dealloc(col_gas);
  pool_t::dealloc(lw_optics_tau_mem);
  pool_t::dealloc(lw_noaero_tau_mem);
  pool_t::dealloc(lay_source_mem);
  pool_t::dealloc(lev_source_inc_mem);
  pool_t::dealloc(lev_source_dec_mem);
  pool_t::dealloc(sfc_source_mem);

  pool_t::dealloc(lw_optics_band2gpt_mem);
  pool_t::dealloc(lw_optics_gpt2band_mem);
  pool_t::dealloc(lw_noaero_band2gpt_mem);
  pool_t::dealloc(lw_noaero_gpt2band_mem);
  pool_t::dealloc(lw_source_band2gpt_mem);
  pool_t::dealloc(lw_source_gpt2band_mem);
}

/*
 * Return a subcolumn mask consistent with a specified overlap assumption
 */
template <typename CldfT, typename SeedsT, typename SubcT>
static void get_subcolumn_mask(const int ncol, const int nlay, const int ngpt, const CldfT &cldf, const int overlap_option, const SeedsT &seeds, const SubcT& subcolumn_mask)
{
  // Subcolumn generators are a means for producing a variable x(i,j,k), where
  //
  //     c(i,j,k) = 1 for x(i,j,k) >  1 - cldf(i,j)
  //     c(i,j,k) = 0 for x(i,j,k) <= 1 - cldf(i,j)
  //
  // I am going to call this "cldx" to be just slightly less ambiguous
  auto cldx = pool_t::template alloc<RealT>(ncol, nlay, ngpt);

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
    //   TIMED_KERNEL(FLATTEN_MD_KERNEL2(nlay, ngpt, ilay, igpt,
    //     auto generator = random_pool.get_state();
    //     cldx(icol,ilay,igpt) = generator.drand(0., 1.);
    //     random_pool.free_state(generator);
    //   }));
    // }
    // TIMED_KERNEL(Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    //   conv::Random rand(seeds(icol));
    //   for (int igpt = 0; igpt < ngpt; igpt++) {
    //     for (int ilay = 0; ilay < nlay; ilay++) {
    //       cldx(icol,ilay,igpt) = rand.genFP<RealT>();
    //     }
    //   }
    // }));
    // Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    //     conv::Random rand(seeds(icol));
    //     for (int igpt = 0; igpt < ngpt; igpt++) {
    //       for (int ilay = 0; ilay < nlay; ilay++) {
    //         cldx(icol,ilay,igpt) = rand.genFP<RealT>();
    //       }
    //     }
    //   });
    TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol, nlay, ngpt, icol, ilay, igpt,
      conv::Random rand(seeds(icol) + ilay*ngpt + igpt);
      cldx(icol,ilay,igpt) = rand.genFP<RealT>();
    ));

    // Step down columns and apply algorithm from eq (14)
    TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol,ngpt, icol, igpt,
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
    ));
  }

  // Use cldx array to create subcolumn mask
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol,nlay,ngpt, icol, ilay, igpt,
    if (cldx(icol,ilay,igpt) > 1.0 - cldf(icol,ilay)) {
      subcolumn_mask(icol,ilay,igpt) = 1;
    } else {
      subcolumn_mask(icol,ilay,igpt) = 0;
    }
  ));

  pool_t::dealloc(cldx);
}

/*
 * Compute cloud area from 3d subcol cloud property
 */
static void compute_cloud_area(
  int ncol, int nlay, int ngpt, Real pmin, Real pmax,
  const creal2dk& pmid, const real3dk& cld_tau_gpt, const real1dk& cld_area)
{
  // Subcolumn binary cld mask; if any layers with pressure between pmin and pmax are cloudy
  // then 2d subcol mask is 1, otherwise it is 0
  auto subcol_mask = pool_t::template alloc_and_init<RealT>(ncol, ngpt);
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol, nlay, ngpt, icol, ilay, igpt,
    // NOTE: using plev would need to assume level ordering (top to bottom or bottom to top), but
    // using play/pmid does not
    if (cld_tau_gpt(icol,ilay,igpt) > 0 && pmid(icol,ilay) >= pmin && pmid(icol,ilay) < pmax) {
      subcol_mask(icol,igpt) = 1;
    }
  ));
  // Compute average over subcols to get cloud area
  auto ngpt_inv = 1.0 / ngpt;
  Kokkos::deep_copy(cld_area, 0);
  TIMED_KERNEL(Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    // This loop needs to be serial because of the atomic reduction
    for (int igpt = 0; igpt < ngpt; ++igpt) {
      cld_area(icol) += subcol_mask(icol,igpt) * ngpt_inv;
    }
  }));

  pool_t::dealloc(subcol_mask);
}

/*
 * Return select cloud-top diagnostics following AeroCom recommendation
 */
static void compute_aerocom_cloudtop(
  int ncol, int nlay, const creal2dk &tmid, const creal2dk &pmid,
  const creal2dk &p_del, const real2dk &z_del, const creal2dk &qc,
  const creal2dk &qi, const creal2dk &rel, const creal2dk &rei,
  const real2dk &cldfrac_tot, const creal2dk &nc,
  const real1dk &T_mid_at_cldtop, const real1dk &p_mid_at_cldtop,
  const real1dk &cldfrac_ice_at_cldtop, const real1dk &cldfrac_liq_at_cldtop,
  const real1dk &cldfrac_tot_at_cldtop, const real1dk &cdnc_at_cldtop,
  const real1dk &eff_radius_qc_at_cldtop, const real1dk &eff_radius_qi_at_cldtop)
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
  auto aerocom_clr = pool_t::template alloc<RealT>(ncol);
  Kokkos::deep_copy(aerocom_clr, 1.0);

  // Get gravity acceleration constant from constants
  using physconst = scream::physics::Constants<RealT>;

  // TODO: move tunable constant to namelist
  constexpr RealT q_threshold = 0.0;  // BAD_CONSTANT!

  // TODO: move tunable constant to namelist
  constexpr RealT cldfrac_tot_threshold = 0.001;  // BAD_CONSTANT!

  // Loop over all columns in parallel
  TIMED_KERNEL(Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
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
  }));

  pool_t::dealloc(aerocom_clr);
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
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ncol, nlay, icol, ilay,
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
  ));
}

/*
 * Routine to limit a quantity to set bounds. Used to make sure
 * effective radii are within the bounds of the cloud optical
 * property look-up tables, but could be used to limit other
 * fields as well.
 */
template<typename InT, typename T, typename OutT, typename std::enable_if<OutT::rank == 1>::type* dummy = nullptr>
static void limit_to_bounds_k(InT const &arr_in, T const lower, T const upper, OutT &arr_out) {
  TIMED_KERNEL(Kokkos::parallel_for(arr_out.size(), KOKKOS_LAMBDA(int i) {
    arr_out(i) = std::min(std::max(arr_in(i), lower), upper);
  }));
}

template<typename InT, typename T, typename OutT, typename std::enable_if<OutT::rank == 2>::type* dummy = nullptr>
static void limit_to_bounds_k(InT const &arr_in, T const lower, T const upper, OutT &arr_out) {
  const int ex0 = (int) arr_out.extent(0);
  const int ex1 = (int) arr_out.extent(1);
  TIMED_KERNEL(FLATTEN_MD_KERNEL2(ex0, ex1, i, j,
    arr_out(i, j) = std::min(std::max(arr_in(i, j), lower), upper);
  ));
}


static int get_wavelength_index(optical_props_t &kdist, RealT wavelength)
{
  auto band_lims_wvn = kdist.get_band_lims_wavenumber();
  auto wavelength_bounds = pool_t::template alloc<RealT>(band_lims_wvn.extent(0), band_lims_wvn.extent(1));
  // Get wavelength bounds for all wavelength bands
  kdist.get_band_lims_wavelength(wavelength_bounds);

  // Find the band index for the specified wavelength
  // Note that bands are stored in wavenumber space, units of cm-1, so if we are passed wavelength
  // in units of meters, we need a conversion factor of 10^2
  const int nbnds = kdist.get_nband();
  int band_index = -1;
  TIMED_KERNEL(Kokkos::parallel_reduce(nbnds, KOKKOS_LAMBDA(int ibnd, int& band_index_inner) {
    if (wavelength_bounds(0,ibnd) < wavelength_bounds(1,ibnd)) {
      if (wavelength_bounds(0,ibnd) <= wavelength * 1e2 && wavelength * 1e2 <= wavelength_bounds(1,ibnd)) {
        band_index_inner = ibnd;
      }
    } else {
      if (wavelength_bounds(0,ibnd) >= wavelength * 1e2 && wavelength * 1e2 >= wavelength_bounds(1,ibnd)) {
        band_index_inner = ibnd;
      }
    }
      }, Kokkos::Max<int>(band_index)));
  pool_t::dealloc(wavelength_bounds);
  return band_index;
}

static inline int get_wavelength_index_sw_k(RealT wavelength)
{
  return get_wavelength_index(*k_dist_sw_k, wavelength);
}

static inline int get_wavelength_index_lw_k(RealT wavelength)
{
  return get_wavelength_index(*k_dist_lw_k, wavelength);
}

static optical_props2_t get_cloud_optics_sw(
  const int ncol, const int nlay,
  cloud_optics_t &cloud_optics, gas_optics_t &kdist,
  const real2dk &lwp, const real2dk &iwp, const creal2dk &rel, const creal2dk &rei,
  const int2dk& sw_band2gpt_mem, const int1dk& sw_gpt2band_mem, const real3dk& sw_tau_mem, const real3dk& sw_ssa_mem, const real3dk& sw_g_mem) {

  // Initialize optics
  optical_props2_t clouds;
  clouds.init_no_alloc(kdist.get_band_lims_wavenumber(), sw_band2gpt_mem, sw_gpt2band_mem);
  clouds.alloc_2str_no_alloc(ncol, nlay, sw_tau_mem, sw_ssa_mem, sw_g_mem);

  // Needed for consistency with all-sky example problem?
  cloud_optics.set_ice_roughness(2);

  // Limit effective radii to be within bounds of lookup table
  auto rel_limited = pool_t::template alloc<RealT>(ncol, nlay);
  auto rei_limited = pool_t::template alloc<RealT>(ncol, nlay);
  limit_to_bounds_k(rel, cloud_optics.radliq_lwr, cloud_optics.radliq_upr, rel_limited);
  limit_to_bounds_k(rei, cloud_optics.radice_lwr, cloud_optics.radice_upr, rei_limited);

  // Calculate cloud optics
  cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel_limited, rei_limited, clouds);

  pool_t::dealloc(rel_limited);
  pool_t::dealloc(rei_limited);

  // Return optics
  return clouds;
}

static optical_props1_t get_cloud_optics_lw(
  const int ncol, const int nlay,
  cloud_optics_t &cloud_optics, gas_optics_t &kdist,
  const real2dk &lwp, const real2dk &iwp, const creal2dk &rel, const creal2dk &rei,
  const int2dk& lw_band2gpt_mem, const int1dk& lw_gpt2band_mem, const real3dk& lw_tau_mem) {

  // Initialize optics
  optical_props1_t clouds;
  clouds.init_no_alloc(kdist.get_band_lims_wavenumber(), lw_band2gpt_mem, lw_gpt2band_mem);
  clouds.alloc_1scl_no_alloc(ncol, nlay, lw_tau_mem);  // this is dumb, why do we need to init and alloc separately?!

  // Needed for consistency with all-sky example problem?
  cloud_optics.set_ice_roughness(2);

  // Limit effective radii to be within bounds of lookup table
  auto rel_limited = pool_t::template alloc<RealT>(ncol, nlay);
  auto rei_limited = pool_t::template alloc<RealT>(ncol, nlay);
  limit_to_bounds_k(rel, cloud_optics.radliq_lwr, cloud_optics.radliq_upr, rel_limited);
  limit_to_bounds_k(rei, cloud_optics.radice_lwr, cloud_optics.radice_upr, rei_limited);

  // Calculate cloud optics
  cloud_optics.cloud_optics(ncol, nlay, lwp, iwp, rel_limited, rei_limited, clouds);

  pool_t::dealloc(rel_limited);
  pool_t::dealloc(rei_limited);

  // Return optics
  return clouds;
}

static optical_props2_t get_subsampled_clouds(
  const int ncol, const int nlay, const int nbnd, const int ngpt,
  optical_props2_t &cloud_optics, gas_optics_t &kdist, const real2dk &cld, const creal2dk &p_lay,
  const int2dk& sw_band2gpt_mem, const int1dk& sw_gpt2band_mem, const real3dk& sw_tau_mem, const real3dk& sw_ssa_mem, const real3dk& sw_g_mem) {
  // Initialized subsampled optics
  optical_props2_t subsampled_optics;
  subsampled_optics.init_no_alloc(kdist.get_band_lims_wavenumber(), kdist.get_band_lims_gpoint(), sw_band2gpt_mem, sw_gpt2band_mem, "subsampled_optics");
  subsampled_optics.alloc_2str_no_alloc(ncol, nlay, sw_tau_mem, sw_ssa_mem, sw_g_mem);

  // Subcolumn mask with values of 0 indicating no cloud, 1 indicating cloud
  auto cldmask = pool_t::template alloc<int>(ncol, nlay, ngpt);

  // Check that we do not have clouds with no optical properties; this would get corrected
  // when we assign optical props, but we want to use a "radiative cloud fraction"
  // for the subcolumn sampling too because otherwise we can get vertically-contiguous cloud
  // mask profiles with no actual cloud properties in between, which would just further overestimate
  // the vertical correlation of cloudy layers. I.e., cloudy layers might look maximally overlapped
  // even when separated by layers with no cloud properties, when in fact those layers should be
  // randomly overlapped.
  auto cldfrac_rad = pool_t::template alloc_and_init<RealT>(ncol, nlay);
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol, nlay, nbnd, icol, ilay, ibnd,
    if (cloud_optics.tau(icol,ilay,ibnd) > 0) {
      cldfrac_rad(icol,ilay) = cld(icol,ilay);
    }
  ));

  // Get subcolumn cloud mask; note that get_subcolumn_mask exposes overlap assumption as an option,
  // but the only currently supported options are 0 (trivial all-or-nothing cloud) or 1 (max-rand),
  // so overlap has not been exposed as an option beyond this subcolumn. In the future, we should
  // support generalized overlap as well, with parameters derived from DPSCREAM simulations with very
  // high resolution.
  int overlap = 1;
  // Get unique seeds for each column that are reproducible across different MPI rank layouts;
  // use decimal part of pressure for this, consistent with the implementation in EAM
  auto seeds = pool_t::template alloc<int>(ncol);
  TIMED_KERNEL(Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    seeds(icol) = 1e9 * (p_lay(icol,nlay-1) - int(p_lay(icol,nlay-1)));
  }));
  get_subcolumn_mask(ncol, nlay, ngpt, cldfrac_rad, overlap, seeds, cldmask);

  // Assign optical properties to subcolumns (note this implements MCICA)
  auto gpoint_bands = kdist.get_gpoint_bands();
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol,nlay,ngpt, icol, ilay, igpt,
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
  ));

  pool_t::dealloc(cldmask);
  pool_t::dealloc(cldfrac_rad);
  pool_t::dealloc(seeds);

  return subsampled_optics;
}

static optical_props1_t get_subsampled_clouds(
  const int ncol, const int nlay, const int nbnd, const int ngpt,
  optical_props1_t &cloud_optics, gas_optics_t &kdist, const real2dk &cld, const creal2dk &p_lay,
  const int2dk& lw_band2gpt_mem, const int1dk& lw_gpt2band_mem, const real3dk& lw_tau_mem) {

  // Initialized subsampled optics
  optical_props1_t subsampled_optics;
  subsampled_optics.init_no_alloc(kdist.get_band_lims_wavenumber(), kdist.get_band_lims_gpoint(), lw_band2gpt_mem, lw_gpt2band_mem, "subsampled_optics");
  subsampled_optics.alloc_1scl_no_alloc(ncol, nlay, lw_tau_mem);

  // Subcolumn mask with values of 0 indicating no cloud, 1 indicating cloud
  auto cldmask = pool_t::template alloc<int>(ncol, nlay, ngpt);

  // Check that we do not have clouds with no optical properties; this would get corrected
  // when we assign optical props, but we want to use a "radiative cloud fraction"
  // for the subcolumn sampling too because otherwise we can get vertically-contiguous cloud
  // mask profiles with no actual cloud properties in between, which would just further overestimate
  // the vertical correlation of cloudy layers. I.e., cloudy layers might look maximally overlapped
  // even when separated by layers with no cloud properties, when in fact those layers should be
  // randomly overlapped.
  auto cldfrac_rad = pool_t::template alloc_and_init<RealT>(ncol, nlay);
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol,nlay,nbnd, icol, ilay, ibnd,
    if (cloud_optics.tau(icol,ilay,ibnd) > 0) {
      cldfrac_rad(icol,ilay) = cld(icol,ilay);
    }
  ));
  // Get subcolumn cloud mask
  int overlap = 1;
  // Get unique seeds for each column that are reproducible across different MPI rank layouts;
  // use decimal part of pressure for this, consistent with the implementation in EAM; use different
  // seed values for longwave and shortwave
  auto seeds = pool_t::template alloc<int>(ncol);
  TIMED_KERNEL(Kokkos::parallel_for(ncol, KOKKOS_LAMBDA(int icol) {
    seeds(icol) = 1e9 * (p_lay(icol,nlay-2) - int(p_lay(icol,nlay-2)));
  }));
  get_subcolumn_mask(ncol, nlay, ngpt, cldfrac_rad, overlap, seeds, cldmask);
  // Assign optical properties to subcolumns (note this implements MCICA)
  auto gpoint_bands = kdist.get_gpoint_bands();
  TIMED_KERNEL(FLATTEN_MD_KERNEL3(ncol,nlay,ngpt, icol, ilay, igpt,
    auto ibnd = gpoint_bands(igpt);
    if (cldmask(icol,ilay,igpt) == 1) {
      subsampled_optics.tau(icol,ilay,igpt) = cloud_optics.tau(icol,ilay,ibnd);
    } else {
      subsampled_optics.tau(icol,ilay,igpt) = 0;
    }
  ));

  pool_t::dealloc(cldmask);
  pool_t::dealloc(cldfrac_rad);
  pool_t::dealloc(seeds);

  return subsampled_optics;
}

}; // struct rrtmgp_interface
#endif // RRTMGP_ENABLE_KOKKOS

} // namespace rrtmgp
} // namespace scream

#endif  // SCREAM_RRTMGP_INTERFACE_HPP
