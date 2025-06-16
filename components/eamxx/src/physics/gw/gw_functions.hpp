#ifndef GW_FUNCTIONS_HPP
#define GW_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"

#include "share/eamxx_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"
#include "ekat/ekat_parameter_list.hpp"

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
  // ---------- GW constants ---------
  //
  struct GWC {
  };

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
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;
  template <typename S>
  using uview_2d = typename ekat::template Unmanaged<view_2d<S> >;

  using MemberType = typename KT::MemberType;

  using WorkspaceManager = typename ekat::WorkspaceManager<Spack, Device>;
  using Workspace        = typename WorkspaceManager::Workspace;

  //
  // --------- Functions ---------
  //

  KOKKOS_FUNCTION
  static void gwd_compute_tendencies_from_stress_divergence(
    // Inputs
    const Int& ncol,
    const Int& pver,
    const Int& pgwv,
    const Int& ngwv,
    const bool& do_taper,
    const Spack& dt,
    const Spack& effgw,
    const uview_1d<const Int>& tend_level,
    const uview_1d<const Spack>& lat,
    const uview_1d<const Spack>& dpm,
    const uview_1d<const Spack>& rdpm,
    const uview_1d<const Spack>& c,
    const uview_1d<const Spack>& ubm,
    const uview_1d<const Spack>& t,
    const uview_1d<const Spack>& nm,
    const uview_1d<const Spack>& xv,
    const uview_1d<const Spack>& yv,
    // Inputs/Outputs
    const uview_1d<Spack>& tau,
    // Outputs
    const uview_1d<Spack>& gwut,
    const uview_1d<Spack>& utgw,
    const uview_1d<Spack>& vtgw);

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
    const Int& pver,
    const Int& ncol,
    const uview_1d<const Int>& tend_level,
    const Spack& dt,
    const uview_1d<const Spack>& taucd,
    const uview_1d<const Spack>& pint,
    const uview_1d<const Spack>& pdel,
    const uview_1d<const Spack>& u,
    const uview_1d<const Spack>& v,
    // Inputs/Outputs
    const uview_1d<Spack>& dudt,
    const uview_1d<Spack>& dvdt,
    const uview_1d<Spack>& dsdt,
    const uview_1d<Spack>& utgw,
    const uview_1d<Spack>& vtgw,
    const uview_1d<Spack>& ttgw);

  KOKKOS_FUNCTION
  static void gwd_compute_stress_profiles_and_diffusivities(
    // Inputs
    const Int& pver,
    const Int& pgwv,
    const Int& ncol,
    const Int& ngwv,
    const uview_1d<const Int>& src_level,
    const uview_1d<const Spack>& ubi,
    const uview_1d<const Spack>& c,
    const uview_1d<const Spack>& rhoi,
    const uview_1d<const Spack>& ni,
    const uview_1d<const Spack>& kvtt,
    const uview_1d<const Spack>& t,
    const uview_1d<const Spack>& ti,
    const uview_1d<const Spack>& piln,
    // Inputs/Outputs
    const uview_1d<Spack>& tau);

  KOKKOS_FUNCTION
  static void gwd_project_tau(
    // Inputs
    const Int& pver,
    const Int& pgwv,
    const Int& ncol,
    const Int& ngwv,
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
    const Int& ngwv,
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
    const Int& ngwv,
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
    const Int& ngwv,
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
    const Int& ngwv,
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
#endif // GPU && !KOKKOS_ENABLE_*_RELOCATABLE_DEVICE_CODE
#endif // P3_FUNCTIONS_HPP
