#ifndef GW_FUNCTIONS_HPP
#define GW_FUNCTIONS_HPP

#include "share/physics/physics_constants.hpp"

#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>
#include <ekat_parameter_list.hpp>

namespace scream {
namespace gw {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for gravity wave drag. We use the ETI pattern for
 * these functions.
 */

template <typename ScalarT, typename DeviceT>
struct Functions
{
  //
  // ------- Types --------
  //

  using Scalar = ScalarT;
  using Device = DeviceT;

  template <typename S>
  using BigPack = ekat::Pack<S,SCREAM_PACK_SIZE>;
  template <typename S>
  using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using IntSmallPack = SmallPack<Int>;
  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using Mask = ekat::Mask<BigPack<Scalar>::n>;
  using Smask = ekat::Mask<SmallPack<Scalar>::n>;

  using KT = KokkosTypes<Device>;

  using C = scream::physics::Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;

  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S>
  using uview_2d = typename ekat::template Unmanaged<view_2d<S> >;
  template <typename S>
  using uview_3d = typename ekat::template Unmanaged<view_3d<S> >;

  using MemberType = typename KT::MemberType;

  using WorkspaceManager = typename ekat::WorkspaceManager<Scalar, Device>;
  using Workspace        = typename WorkspaceManager::Workspace;

  //
  // ---------- GW constants ---------
  //
  struct GWC {
    // Index the cardinal directions.
    static inline constexpr int west = 0;
    static inline constexpr int east = 1;
    static inline constexpr int south = 2;
    static inline constexpr int north = 3;

    // rair/gravit
    static inline constexpr Real rog = C::Rair / C::gravit;

    // Background diffusivity.
    static inline constexpr Real dback = 0.05;

    // Minimum non-zero stress.
    static inline constexpr Real taumin = 1.e-10;

    // Maximum allowed change in u-c (before efficiency applied).
    static inline constexpr Real umcfac = 0.5;

    // Minimum value of (u-c)**2.
    static inline constexpr Real ubmc2mn = 0.01;

    // Useful to avoid promoting to 64-bit double in single prec
    static constexpr Real half=0.5;

    // Minimum value of Brunt-Vaisalla frequency squared.
    static constexpr Real n2min = 1.e-8;

    // Mysterious hardcoded value used in calculation of tndmax
    static constexpr Real temp1 = 86400;

    // spectrum averaging length [m]
    static constexpr Real tau_avg_length = 100e3;

    // Inverse Prandtl number.
    static constexpr Real prndl=0.25;

    // Integration interval to get bin average.
    static constexpr Real dca  = 0.1;

    // Width of gaussian in phase speed.
    static constexpr Real c0   = 30;

    // max altitude [m] to check for max heating
    static constexpr Real heating_altitude_max = 20e3;

    // min surface displacement height for orographic waves
    static constexpr Real orohmin = 10;

    // min wind speed for orographic waves
    static constexpr Real orovmin = 2;
  };

  //
  // --------- Structs for organizing data ---------
  //

  struct GwCommonInit {
    GwCommonInit() : initialized(false) {}

    // Tell us if initialize has been called
    bool initialized;

    // This flag preserves answers for vanilla CAM by making a few changes (e.g.
    // order of operations) when only orographic waves are on.
    bool orographic_only; // = .false.

    // Number of levels in the atmosphere.
    int pver;

    // Maximum number of waves allowed (i.e. wavenumbers are -pgwv:pgwv).
    int pgwv;

    // Bin width for spectrum.
    Real dc; // = huge(1._r8)

    // Reference speeds for the spectrum.
    view_1d<Real> cref;

    // Whether or not molecular diffusion is being done, and bottom level where
    // it is done.
    bool do_molec_diff; // = .false.
    int nbot_molec; // = huge(1)

    // Whether or not to enforce an upper boundary condition of tau = 0.
    bool tau_0_ubc; // = .false.

    // Critical Froude number.
    Real fcrit2; // = huge(1._r8)

    // Effective horizontal wave number.
    Real kwv; // = huge(1._r8)

    // 1/2 * horizontal wavenumber (for oro)
    Real oroko2;

    // Interface levels for gravity wave sources.
    int ktop; // = huge(1)
    int kbotbg; // = huge(1)

    // Effective wavenumber.
    Real effkwv; // = huge(1._r8)

    // Newtonian cooling coefficients.
    view_1d<Real> alpha;

    // Maximum wind tendency from stress divergence (before efficiency applied).
    Real tndmax; // = huge(1._r8)
  };

  struct GwConvectInit {
    GwConvectInit() : initialized(false) {}

    // Tell us if initialize has been called
    bool initialized;

    // Dimension for heating depth.
    int maxh;

    // Dimension for mean wind in heating.
    int maxuh;

    // Index for level for storm/steering flow (usually 700 mb)
    int k_src_wind;

    // Table of source spectra.
    view_3d<Real> mfcc;
  };

  struct GwFrontInit {
    GwFrontInit() : initialized(false) {}

    // Tell us if initialize has been called
    bool initialized;

    // Frontogenesis function critical threshold.
    Real frontgfc;

    // Level at which to check the frontogenesis function to determine when
    // waves will be launched.
    Int kfront;

    // Average value of gaussian over gravity wave spectrum bins, multiplied by
    // background source strength (taubgnd).
    view_1d<Real> fav;
  };

  //
  // --------- Util Functions ---------
  //

  // Pure function that interpolates the values of the input array along
  // dimension 2. This is obviously not a very generic routine, unlike, say,
  // CAM's lininterp. But it's used often enough that it seems worth providing
  // here.
  KOKKOS_INLINE_FUNCTION
  static void midpoint_interp(
    const MemberType& team,
    const uview_1d<const Real>& in,
    const uview_1d<Real>& interp)
  {
    EKAT_KERNEL_REQUIRE(in.size() == interp.size() + 1);
    Kokkos::parallel_for(
      Kokkos::TeamVectorRange(team, 0, in.extent(0)-1), [&] (const int k) {
        interp(k) = (in(k)+in(k+1)) / 2;
    });
  }

  // Take two components of a vector, and find the unit vector components and
  // total magnitude.
  KOKKOS_INLINE_FUNCTION
  static void get_unit_vector(const Real& u, const Real& v, Real& u_n, Real& v_n, Real& mag)
  {
    mag = sqrt(u*u + v*v);
    if (mag > 0) {
      u_n = u / mag;
      v_n = v / mag;
    }
    else {
      u_n = 0;
      v_n = 0;
    }
  }

  // Elemental version of a 2D dot product (since the intrinsic dot_product
  // is more suitable for arrays of contiguous vectors).
  KOKKOS_INLINE_FUNCTION
  static Real dot_2d(const Real& u1, const Real& v1, const Real& u2, const Real& v2)
  {
    return u1*u2 + v1*v2;
  }

  //
  // --------- Init/Finalize Functions ---------
  //
  static void gw_common_init(
    // Inputs
    const Int& pver_in,
    const Int& pgwv_in,
    const Real& dc_in,
    const uview_1d<const Real>& cref_in,
    const bool& orographic_only_in,
    const bool& do_molec_diff_in,
    const bool& tau_0_ubc_in,
    const Int& nbot_molec_in,
    const Int& ktop_in,
    const Int& kbotbg_in,
    const Real& fcrit2_in,
    const Real& kwv_in,
    const uview_1d<const Real>& alpha_in);

  static void gw_convect_init(
    // Inputs
    const Real& plev_src_wind,
    const uview_3d<const Real>& mfcc_in);

  static void gw_front_init(
    // Inputs
    const Real& taubgnd,
    const Real& frontgfc_in,
    const Int& kfront_in);

  static void gw_finalize()
  {
    s_common_init.cref  = decltype(s_common_init.cref)();
    s_common_init.alpha = decltype(s_common_init.alpha)();
    s_convect_init.mfcc = decltype(s_convect_init.mfcc)();
    s_front_init.fav    = decltype(s_front_init.fav)();
  }

  //
  // --------- Functions ---------
  //

  KOKKOS_FUNCTION
  static void gwd_compute_tendencies_from_stress_divergence(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const GwCommonInit& init,
    const Int& pver,
    const Int& pgwv,
    const bool& do_taper,
    const Real& dt,
    const Real& effgw,
    const Int& tend_level,
    const Int& max_level,
    const Real& lat,
    const uview_1d<const Real>& dpm,
    const uview_1d<const Real>& rdpm,
    const uview_1d<const Real>& c,
    const uview_1d<const Real>& ubm,
    const uview_1d<const Real>& t,
    const uview_1d<const Real>& nm,
    const Real& xv,
    const Real& yv,
    // Inputs/Outputs
    const uview_2d<Real>& tau,
    // Outputs
    const uview_2d<Real>& gwut,
    const uview_1d<Real>& utgw,
    const uview_1d<Real>& vtgw);

  /*
   * Compute profiles of background state quantities for the multiple
   * gravity wave drag parameterization.
   *
   * The parameterization is assumed to operate only where water vapor
   * concentrations are negligible in determining the density.
   */
  KOKKOS_FUNCTION
  static void gw_prof(
    // Inputs
    const MemberType& team,
    const Int& pver,
    const Real& cpair,
    const uview_1d<const Real>& t,
    const uview_1d<const Real>& pmid,
    const uview_1d<const Real>& pint,
    // Outputs
    const uview_1d<Real>& rhoi,
    const uview_1d<Real>& ti,
    const uview_1d<Real>& nm,
    const uview_1d<Real>& ni);

  KOKKOS_FUNCTION
  static void momentum_energy_conservation(
    // Inputs
    const MemberType& team,
    const Int& pver,
    const Int& tend_level,
    const Real& dt,
    const uview_2d<const Real>& taucd,
    const uview_1d<const Real>& pint,
    const uview_1d<const Real>& pdel,
    const uview_1d<const Real>& u,
    const uview_1d<const Real>& v,
    // Inputs/Outputs
    const uview_1d<Real>& dudt,
    const uview_1d<Real>& dvdt,
    const uview_1d<Real>& dsdt,
    const uview_1d<Real>& utgw,
    const uview_1d<Real>& vtgw,
    const uview_1d<Real>& ttgw);

  KOKKOS_FUNCTION
  static void gwd_compute_stress_profiles_and_diffusivities(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const GwCommonInit& init,
    const Int& pver,
    const Int& pgwv,
    const Int& src_level,
    const uview_1d<const Real>& ubi,
    const uview_1d<const Real>& c,
    const uview_1d<const Real>& rhoi,
    const uview_1d<const Real>& ni,
    const uview_1d<const Real>& kvtt,
    const uview_1d<const Real>& t,
    const uview_1d<const Real>& ti,
    const uview_1d<const Real>& piln,
    // Inputs/Outputs
    const uview_2d<Real>& tau);

  KOKKOS_FUNCTION
  static void gwd_project_tau(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const GwCommonInit& init,
    const Int& pver,
    const Int& pgwv,
    const Int& tend_level,
    const uview_2d<const Real>& tau,
    const uview_1d<const Real>& ubi,
    const uview_1d<const Real>& c,
    const Real& xv,
    const Real& yv,
    // Outputs
    const uview_2d<Real>& taucd);

  KOKKOS_FUNCTION
  static void gwd_precalc_rhoi(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const GwCommonInit& init,
    const Int& pver,
    const Int& pgwv,
    const Real& dt,
    const Int& tend_level,
    const uview_1d<const Real>& pmid,
    const uview_1d<const Real>& pint,
    const uview_1d<const Real>& t,
    const uview_2d<const Real>& gwut,
    const uview_1d<const Real>& ubm,
    const uview_1d<const Real>& nm,
    const uview_1d<const Real>& rdpm,
    const uview_1d<const Real>& c,
    const uview_2d<const Real>& q,
    const uview_1d<const Real>& dse,
    // Outputs
    const uview_1d<Real>& egwdffi,
    const uview_2d<Real>& qtgw,
    const uview_1d<Real>& dttdf,
    const uview_1d<Real>& dttke,
    const uview_1d<Real>& ttgw);

  KOKKOS_FUNCTION
  static void gw_drag_prof(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const GwCommonInit& init,
    const Int& pver,
    const Int& pgwv,
    const Int& src_level,
    const Int& max_level,
    const Int& tend_level,
    const bool& do_taper,
    const Real& dt,
    const Real& lat,
    const uview_1d<const Real>& t,
    const uview_1d<const Real>& ti,
    const uview_1d<const Real>& pmid,
    const uview_1d<const Real>& pint,
    const uview_1d<const Real>& dpm,
    const uview_1d<const Real>& rdpm,
    const uview_1d<const Real>& piln,
    const uview_1d<const Real>& rhoi,
    const uview_1d<const Real>& nm,
    const uview_1d<const Real>& ni,
    const uview_1d<const Real>& ubm,
    const uview_1d<const Real>& ubi,
    const Real& xv,
    const Real& yv,
    const Real& effgw,
    const uview_1d<const Real>& c,
    const uview_1d<const Real>& kvtt,
    const uview_2d<const Real>& q,
    const uview_1d<const Real>& dse,
    // Inputs/Outputs
    const uview_2d<Real>& tau,
    // Outputs
    const uview_1d<Real>& utgw,
    const uview_1d<Real>& vtgw,
    const uview_1d<Real>& ttgw,
    const uview_2d<Real>& qtgw,
    const uview_2d<Real>& taucd,
    const uview_1d<Real>& egwdffi,
    const uview_2d<Real>& gwut,
    const uview_1d<Real>& dttdf,
    const uview_1d<Real>& dttke);

  KOKKOS_FUNCTION
  static void gw_front_project_winds(
    // Inputs
    const MemberType& team,
    const Int& pver,
    const Int& kbot,
    const uview_1d<const Real>& u,
    const uview_1d<const Real>& v,
    // Outputs
    Real& xv,
    Real& yv,
    const uview_1d<Real>& ubm,
    const uview_1d<Real>& ubi);

  KOKKOS_FUNCTION
  static void gw_front_gw_sources(
    // Inputs
    const MemberType& team,
    const GwFrontInit& finit,
    const Int& pgwv,
    const Int& pver,
    const Int& kbot,
    const uview_1d<const Real>& frontgf,
    // Outputs
    const uview_2d<Real>& tau);

  KOKKOS_FUNCTION
  static void gw_cm_src(
    // Inputs
    const MemberType& team,
    const GwCommonInit& init,
    const GwFrontInit& finit,
    const Int& pver,
    const Int& pgwv,
    const Int& kbot,
    const uview_1d<const Real>& u,
    const uview_1d<const Real>& v,
    const uview_1d<const Real>& frontgf,
    // Outputs
    Int& src_level,
    Int& tend_level,
    const uview_2d<Real>& tau,
    const uview_1d<Real>& ubm,
    const uview_1d<Real>& ubi,
    Real& xv,
    Real& yv,
    const uview_1d<Real>& c);

  KOKKOS_FUNCTION
  static void gw_convect_project_winds(
    // Inputs
    const MemberType& team,
    const GwConvectInit& init,
    const Int& pver,
    const uview_1d<const Real>& u,
    const uview_1d<const Real>& v,
    // Outputs
    Real& xv,
    Real& yv,
    const uview_1d<Real>& ubm,
    const uview_1d<Real>& ubi);

  KOKKOS_FUNCTION
  static void gw_heating_depth(
    // Inputs
    const MemberType& team,
    const GwConvectInit& init,
    const Int& pver,
    const Real& maxq0_conversion_factor,
    const Real& hdepth_scaling_factor,
    const bool& use_gw_convect_old,
    const uview_1d<const Real>& zm,
    const uview_1d<const Real>& netdt,
    // Outputs
    Int& mini,
    Int& maxi,
    Real& hdepth,
    Real& maxq0_out,
    Real& maxq0);

  KOKKOS_FUNCTION
  static void gw_storm_speed(
    // Inputs
    const MemberType& team,
    const GwCommonInit& init,
    const GwConvectInit& cinit,
    const Int& pver,
    const Real& storm_speed_min,
    const uview_1d<const Real>& ubm,
    const Int& mini,
    const Int& maxi,
    // Outputs
    Int& storm_speed,
    Real& uh,
    Real& umin,
    Real& umax);

  KOKKOS_FUNCTION
  static void gw_convect_gw_sources(
    // Inputs
    const MemberType& team,
    const GwCommonInit& init,
    const GwConvectInit& cinit,
    const Int& pgwv,
    const Int& pver,
    const Real& lat,
    const Real& hdepth_min,
    const Real& hdepth,
    const Int& mini,
    const Int& maxi,
    const uview_1d<const Real>& netdt,
    const Real& uh,
    const Int& storm_speed,
    const Real& maxq0,
    const Real& umin,
    const Real& umax,
    // Outputs
    const uview_2d<Real>& tau);

  KOKKOS_FUNCTION
  static void gw_beres_src(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const GwCommonInit& init,
    const GwConvectInit& cinit,
    const Int& pver,
    const Int& pgwv,
    const Real& lat,
    const uview_1d<const Real>& u,
    const uview_1d<const Real>& v,
    const uview_1d<const Real>& netdt,
    const uview_1d<const Real>& zm,
    const Real& maxq0_conversion_factor,
    const Real& hdepth_scaling_factor,
    const Real& hdepth_min,
    const Real& storm_speed_min,
    const bool& use_gw_convect_old,
    // Outputs
    Int& src_level,
    Int& tend_level,
    const uview_2d<Real>& tau,
    const uview_1d<Real>& ubm,
    const uview_1d<Real>& ubi,
    Real& xv,
    Real& yv,
    const uview_1d<Real>& c,
    Real& hdepth,
    Real& maxq0_out);

  KOKKOS_FUNCTION
  static void gw_ediff(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const Int& pver,
    const Int& pgwv,
    const Int& kbot,
    const Int& ktop,
    const Int& tend_level,
    const Real& dt,
    const uview_2d<const Real>& gwut,
    const uview_1d<const Real>& ubm,
    const uview_1d<const Real>& nm,
    const uview_1d<const Real>& rho,
    const uview_1d<const Real>& pmid,
    const uview_1d<const Real>& rdpm,
    const uview_1d<const Real>& c,
    // Outputs
    const uview_1d<Real>& egwdffi,
    const uview_1d<Real>& decomp_ca,
    const uview_1d<Real>& decomp_cc,
    const uview_1d<Real>& decomp_dnom,
    const uview_1d<Real>& decomp_ze);

  KOKKOS_FUNCTION
  static void gw_diff_tend(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const Int& pver,
    const Int& kbot,
    const Int& ktop,
    const uview_1d<const Real>& q,
    const Real& dt,
    const uview_1d<const Real>& decomp_ca,
    const uview_1d<const Real>& decomp_cc,
    const uview_1d<const Real>& decomp_dnom,
    const uview_1d<const Real>& decomp_ze,
      // Outputs
    const uview_1d<Real>& dq);

  KOKKOS_FUNCTION
  static void gw_oro_src(
    // Inputs
    const MemberType& team,
    const GwCommonInit& init,
    const Int& pver,
    const Int& pgwv,
    const uview_1d<const Real>& u,
    const uview_1d<const Real>& v,
    const uview_1d<const Real>& t,
    const Real& sgh,
    const uview_1d<const Real>& pmid,
    const uview_1d<const Real>& pint,
    const uview_1d<const Real>& dpm,
    const uview_1d<const Real>& zm,
    const uview_1d<const Real>& nm,
    // Outputs
    Int& src_level,
    Int& tend_level,
    const uview_2d<Real>& tau,
    const uview_1d<Real>& ubm,
    const uview_1d<Real>& ubi,
    Real& xv,
    Real& yv,
    const uview_1d<Real>& c);

  KOKKOS_FUNCTION
  static void vd_lu_decomp(
    // Inputs
    const MemberType& team,
    const Int& pver,
    const Real& ksrf,
    const uview_1d<const Real>& kv,
    const uview_1d<const Real>& tmpi,
    const uview_1d<const Real>& rpdel,
    const Real& ztodt,
    const Real& cc_top,
    const Int& ntop,
    const Int& nbot,
    // Outputs
    const uview_1d<Real>& decomp_ca,
    const uview_1d<Real>& decomp_cc,
    const uview_1d<Real>& decomp_dnom,
    const uview_1d<Real>& decomp_ze);

  KOKKOS_FUNCTION
  static void vd_lu_solve(
    // Inputs
    const MemberType& team,
    const Workspace& workspace,
    const Int& pver,
    const uview_1d<const Real>& decomp_ca,
    const uview_1d<const Real>& decomp_cc,
    const uview_1d<const Real>& decomp_dnom,
    const uview_1d<const Real>& decomp_ze,
    const Int& ntop,
    const Int& nbot,
    const Real& cd_top,
    // Inputs/Outputs
    const uview_1d<Real>& q);

  //
  // --------- Members ---------
  //
  inline static GwCommonInit s_common_init;
  inline static GwConvectInit s_convect_init;
  inline static GwFrontInit s_front_init;

}; // struct Functions

} // namespace gw
} // namespace scream

// If a GPU build, without relocatable device code enabled, make all code available
// to the translation unit; otherwise, ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) \
                                && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
# include "impl/gw_compute_tendencies_from_stress_divergence_impl.hpp"
# include "impl/gw_prof_impl.hpp"
# include "impl/gw_momentum_energy_conservation_impl.hpp"
# include "impl/gw_compute_stress_profiles_and_diffusivities_impl.hpp"
# include "impl/gw_project_tau_impl.hpp"
# include "impl/gw_precalc_rhoi_impl.hpp"
# include "impl/gw_drag_prof_impl.hpp"
# include "impl/gw_front_project_winds_impl.hpp"
# include "impl/gw_front_gw_sources_impl.hpp"
# include "impl/gw_cm_src_impl.hpp"
# include "impl/gw_convect_project_winds_impl.hpp"
# include "impl/gw_heating_depth_impl.hpp"
# include "impl/gw_storm_speed_impl.hpp"
# include "impl/gw_convect_gw_sources_impl.hpp"
# include "impl/gw_beres_src_impl.hpp"
# include "impl/gw_ediff_impl.hpp"
# include "impl/gw_diff_tend_impl.hpp"
# include "impl/gw_oro_src_impl.hpp"
# include "impl/gw_common_init_impl.hpp"
# include "impl/gw_vd_lu_decomp_impl.hpp"
# include "impl/gw_vd_lu_solve_impl.hpp"
# include "impl/gw_convect_init_impl.hpp"
# include "impl/gw_front_init_impl.hpp"
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE
#endif // P3_FUNCTIONS_HPP
