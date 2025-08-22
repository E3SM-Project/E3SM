#ifndef GW_FUNCTIONS_HPP
#define GW_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"

#include "share/eamxx_types.hpp"

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

  static void gw_common_finalize()
  {
    s_common_init.cref  = decltype(s_common_init.cref)();
    s_common_init.alpha = decltype(s_common_init.alpha)();
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

  KOKKOS_FUNCTION
  static void gw_prof(
    // Inputs
    const Int& pver,
    const Int& ncol,
    const Spack& cpair,
    const uview_1d<const Spack>& t,
    const uview_1d<const Spack>& pmid,
    const uview_1d<const Spack>& pint,
    // Outputs
    const uview_1d<Spack>& rhoi,
    const uview_1d<Spack>& ti,
    const uview_1d<Spack>& nm,
    const uview_1d<Spack>& ni);

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
    const Int& pver,
    const Int& pgwv,
    const Int& ncol,
    const uview_1d<const Int>& tend_level,
    const uview_1d<const Spack>& tau,
    const uview_1d<const Spack>& ubi,
    const uview_1d<const Spack>& c,
    const uview_1d<const Spack>& xv,
    const uview_1d<const Spack>& yv,
    // Outputs
    const uview_1d<Spack>& taucd);

  KOKKOS_FUNCTION
  static void gwd_precalc_rhoi(
    // Inputs
    const Int& pver,
    const Int& pgwv,
    const Int& ncol,
    const Spack& dt,
    const uview_1d<const Int>& tend_level,
    const uview_1d<const Spack>& pmid,
    const uview_1d<const Spack>& pint,
    const uview_1d<const Spack>& t,
    const uview_1d<const Spack>& gwut,
    const uview_1d<const Spack>& ubm,
    const uview_1d<const Spack>& nm,
    const uview_1d<const Spack>& rdpm,
    const uview_1d<const Spack>& c,
    const uview_1d<const Spack>& q,
    const uview_1d<const Spack>& dse,
    // Outputs
    const uview_1d<Spack>& egwdffi,
    const uview_1d<Spack>& qtgw,
    const uview_1d<Spack>& dttdf,
    const uview_1d<Spack>& dttke,
    const uview_1d<Spack>& ttgw);

  KOKKOS_FUNCTION
  static void gw_drag_prof(
    // Inputs
    const Int& pver,
    const Int& pgwv,
    const Int& ncol,
    const uview_1d<const Int>& src_level,
    const uview_1d<const Int>& tend_level,
    const bool& do_taper,
    const Spack& dt,
    const uview_1d<const Spack>& lat,
    const uview_1d<const Spack>& t,
    const uview_1d<const Spack>& ti,
    const uview_1d<const Spack>& pmid,
    const uview_1d<const Spack>& pint,
    const uview_1d<const Spack>& dpm,
    const uview_1d<const Spack>& rdpm,
    const uview_1d<const Spack>& piln,
    const uview_1d<const Spack>& rhoi,
    const uview_1d<const Spack>& nm,
    const uview_1d<const Spack>& ni,
    const uview_1d<const Spack>& ubm,
    const uview_1d<const Spack>& ubi,
    const uview_1d<const Spack>& xv,
    const uview_1d<const Spack>& yv,
    const Spack& effgw,
    const uview_1d<const Spack>& c,
    const uview_1d<const Spack>& kvtt,
    const uview_1d<const Spack>& q,
    const uview_1d<const Spack>& dse,
    // Inputs/Outputs
    const uview_1d<Spack>& tau,
    // Outputs
    const uview_1d<Spack>& utgw,
    const uview_1d<Spack>& vtgw,
    const uview_1d<Spack>& ttgw,
    const uview_1d<Spack>& qtgw,
    const uview_1d<Spack>& taucd,
    const uview_1d<Spack>& egwdffi,
    const uview_1d<Spack>& gwut,
    const uview_1d<Spack>& dttdf,
    const uview_1d<Spack>& dttke);

  KOKKOS_FUNCTION
  static void gw_front_project_winds(
    // Inputs
    const Int& pver,
    const Int& ncol,
    const Int& kbot,
    const uview_1d<const Spack>& u,
    const uview_1d<const Spack>& v,
    // Outputs
    const uview_1d<Spack>& xv,
    const uview_1d<Spack>& yv,
    const uview_1d<Spack>& ubm,
    const uview_1d<Spack>& ubi);

  KOKKOS_FUNCTION
  static void gw_front_gw_sources(
    // Inputs
    const Int& pver,
    const Int& pgwv,
    const Int& ncol,
    const Int& kbot,
    const uview_1d<const Spack>& frontgf,
    // Outputs
    const uview_1d<Spack>& tau);

  KOKKOS_FUNCTION
  static void gw_cm_src(
    // Inputs
    const Int& pver,
    const Int& pgwv,
    const Int& ncol,
    const Int& kbot,
    const uview_1d<const Spack>& u,
    const uview_1d<const Spack>& v,
    const uview_1d<const Spack>& frontgf,
    // Outputs
    const uview_1d<Int>& src_level,
    const uview_1d<Int>& tend_level,
    const uview_1d<Spack>& tau,
    const uview_1d<Spack>& ubm,
    const uview_1d<Spack>& ubi,
    const uview_1d<Spack>& xv,
    const uview_1d<Spack>& yv,
    const uview_1d<Spack>& c);

  KOKKOS_FUNCTION
  static void gw_convect_project_winds(
    // Inputs
    const Int& pver,
    const Int& ncol,
    const uview_1d<const Spack>& u,
    const uview_1d<const Spack>& v,
    // Outputs
    const uview_1d<Spack>& xv,
    const uview_1d<Spack>& yv,
    const uview_1d<Spack>& ubm,
    const uview_1d<Spack>& ubi);

  KOKKOS_FUNCTION
  static void gw_heating_depth(
    // Inputs
    const Int& pver,
    const Int& ncol,
    const Spack& maxq0_conversion_factor,
    const Spack& hdepth_scaling_factor,
    const bool& use_gw_convect_old,
    const uview_1d<const Spack>& zm,
    const uview_1d<const Spack>& netdt,
    // Outputs
    const uview_1d<Int>& mini,
    const uview_1d<Int>& maxi,
    const uview_1d<Spack>& hdepth,
    const uview_1d<Spack>& maxq0_out,
    const uview_1d<Spack>& maxq0);

  KOKKOS_FUNCTION
  static void gw_storm_speed(
    // Inputs
    const Int& pver,
    const Int& ncol,
    const Spack& storm_speed_min,
    const uview_1d<const Spack>& ubm,
    const uview_1d<const Int>& mini,
    const uview_1d<const Int>& maxi,
    // Outputs
    const uview_1d<Int>& storm_speed,
    const uview_1d<Spack>& uh,
    const uview_1d<Spack>& umin,
    const uview_1d<Spack>& umax);

  KOKKOS_FUNCTION
  static void gw_convect_gw_sources(
    // Inputs
    const Int& pver,
    const Int& pgwv,
    const Int& ncol,
    const uview_1d<const Spack>& lat,
    const Spack& hdepth_min,
    const uview_1d<const Spack>& hdepth,
    const uview_1d<const Int>& mini,
    const uview_1d<const Int>& maxi,
    const uview_1d<const Spack>& netdt,
    const uview_1d<const Spack>& uh,
    const uview_1d<const Int>& storm_speed,
    const uview_1d<const Spack>& maxq0,
    const uview_1d<const Spack>& umin,
    const uview_1d<const Spack>& umax,
    // Outputs
    const uview_1d<Spack>& tau);

  KOKKOS_FUNCTION
  static void gw_beres_src(
    // Inputs
    const Int& pver,
    const Int& pgwv,
    const Int& ncol,
    const uview_1d<const Spack>& lat,
    const uview_1d<const Spack>& u,
    const uview_1d<const Spack>& v,
    const uview_1d<const Spack>& netdt,
    const uview_1d<const Spack>& zm,
    // Outputs
    const uview_1d<Int>& src_level,
    const uview_1d<Int>& tend_level,
    const uview_1d<Spack>& tau,
    const uview_1d<Spack>& ubm,
    const uview_1d<Spack>& ubi,
    const uview_1d<Spack>& xv,
    const uview_1d<Spack>& yv,
    const uview_1d<Spack>& c,
    const uview_1d<Spack>& hdepth,
    const uview_1d<Spack>& maxq0_out,
    // Inputs
    const Spack& maxq0_conversion_factor,
    const Spack& hdepth_scaling_factor,
    const Spack& hdepth_min,
    const Spack& storm_speed_min,
    const bool& use_gw_convect_old);

  KOKKOS_FUNCTION
  static void gw_ediff(
    // Inputs
    const Int& ncol,
    const Int& pver,
    const Int& kbot,
    const Int& ktop,
    const uview_1d<const Int>& tend_level,
    const uview_1d<const Spack>& gwut,
    const uview_1d<const Spack>& ubm,
    const uview_1d<const Spack>& nm,
    const uview_1d<const Spack>& rho,
    const Spack& dt,
    const Spack& gravit,
    const uview_1d<const Spack>& pmid,
    const uview_1d<const Spack>& rdpm,
    const uview_1d<const Spack>& c,
    // Outputs
    const uview_1d<Spack>& egwdffi);

  KOKKOS_FUNCTION
  static void gw_diff_tend(
    // Inputs
    const Int& ncol,
    const Int& pver,
    const Int& kbot,
    const Int& ktop,
    const uview_1d<const Spack>& q,
    const Spack& dt,
      // Outputs
    const uview_1d<Spack>& dq);

  KOKKOS_FUNCTION
  static void gw_oro_src(
    // Inputs
    const Int& pver,
    const Int& pgwv,
    const Int& ncol,
    const uview_1d<const Spack>& u,
    const uview_1d<const Spack>& v,
    const uview_1d<const Spack>& t,
    const uview_1d<const Spack>& sgh,
    const uview_1d<const Spack>& pmid,
    const uview_1d<const Spack>& pint,
    const uview_1d<const Spack>& dpm,
    const uview_1d<const Spack>& zm,
    const uview_1d<const Spack>& nm,
    // Outputs
    const uview_1d<Int>& src_level,
    const uview_1d<Int>& tend_level,
    const uview_1d<Spack>& tau,
    const uview_1d<Spack>& ubm,
    const uview_1d<Spack>& ubi,
    const uview_1d<Spack>& xv,
    const uview_1d<Spack>& yv,
    const uview_1d<Spack>& c);

  //
  // --------- Members ---------
  //
  inline static GwCommonInit s_common_init;

}; // struct Functions

} // namespace gw
} // namespace scream

// If a GPU build, without relocatable device code enabled, make all code available
// to the translation unit; otherwise, ETI is used.
#if defined(EAMXX_ENABLE_GPU) && !defined(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE) \
                                && !defined(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
# include "impl/gw_gwd_compute_tendencies_from_stress_divergence_impl.hpp"
# include "impl/gw_gw_prof_impl.hpp"
# include "impl/gw_momentum_energy_conservation_impl.hpp"
# include "impl/gw_gwd_compute_stress_profiles_and_diffusivities_impl.hpp"
# include "impl/gw_gwd_project_tau_impl.hpp"
# include "impl/gw_gwd_precalc_rhoi_impl.hpp"
# include "impl/gw_gw_drag_prof_impl.hpp"
# include "impl/gw_gw_front_project_winds_impl.hpp"
# include "impl/gw_gw_front_gw_sources_impl.hpp"
# include "impl/gw_gw_cm_src_impl.hpp"
# include "impl/gw_gw_convect_project_winds_impl.hpp"
# include "impl/gw_gw_heating_depth_impl.hpp"
# include "impl/gw_gw_storm_speed_impl.hpp"
# include "impl/gw_gw_convect_gw_sources_impl.hpp"
# include "impl/gw_gw_beres_src_impl.hpp"
# include "impl/gw_gw_ediff_impl.hpp"
# include "impl/gw_gw_diff_tend_impl.hpp"
# include "impl/gw_gw_oro_src_impl.hpp"
# include "impl/gw_gw_common_init_impl.hpp"
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE
#endif // P3_FUNCTIONS_HPP
