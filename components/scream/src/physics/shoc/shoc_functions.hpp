#ifndef SHOC_FUNCTIONS_HPP
#define SHOC_FUNCTIONS_HPP

#include "physics/share/physics_constants.hpp"

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace shoc {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for SHOC. We use the ETI pattern for
 * these functions.
 *
 * SHOC assumptions:
 *  - Kokkos team policies have a vector length of 1
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

  using Mask  = ekat::Mask<Pack::n>;
  using Smask = ekat::Mask<Spack::n>;

  using KT = ekat::KokkosTypes<Device>;

  using C = physics::Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  template <typename S>
  using uview_2d = typename ekat::template Unmanaged<view_2d<S> >;

  using MemberType = typename KT::MemberType;

  using Workspace = typename ekat::WorkspaceManager<Spack, Device>::Workspace;

  //
  // --------- Functions ---------
  //
  KOKKOS_FUNCTION
  static void calc_shoc_varorcovar(
    const MemberType&            team,
    const Int&                   nlev, 
    const Scalar&                tunefac,
    const uview_1d<const Spack>& isotropy_zi, 
    const uview_1d<const Spack>& tkh_zi,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& invar1, 
    const uview_1d<const Spack>& invar2, 
    const uview_1d<Spack>&       varorcovar);

  KOKKOS_FUNCTION
  static void calc_shoc_vertflux(
    const MemberType& team,
    const Int& nlev,
    const uview_1d<const Spack>& tkh_zi,
    const uview_1d<const Spack>& dz_zi,
    const uview_1d<const Spack>& invar,
    const uview_1d<Spack>& vertflux);

  KOKKOS_FUNCTION
  static void shoc_diag_second_moments_srf(
    const Scalar& wthl_sfc, const Scalar& uw_sfc, const Scalar& vw_sfc,
    Scalar& ustar2, Scalar& wstar);

  KOKKOS_FUNCTION
  static void shoc_diag_second_moments_ubycond(
    Scalar& thl_sec, Scalar& qw_sec, Scalar& wthl_sec, Scalar& wqw_sec,
    Scalar& qwthl_sec, Scalar& uw_sec, Scalar& vw_sec, Scalar& wtke_sec);

  KOKKOS_FUNCTION
  static void update_host_dse(
    const MemberType& team,
    const Int& nlev,
    const uview_1d<const Spack>& thlm,
    const uview_1d<const Spack>& shoc_ql,
    const uview_1d<const Spack>& exner,
    const uview_1d<const Spack>& zt_grid,
    const Scalar& phis,
    const uview_1d<Spack>& host_dse);

  KOKKOS_FUNCTION
  static void shoc_pblintd_init_pot(
    const MemberType& team, const Int& nlev,
    const view_1d<const Spack>& thl, const view_1d<const Spack>& ql, const view_1d<const Spack>& q,
    const view_1d<Spack>& thv);

}; // struct Functions

} // namespace shoc
} // namespace scream

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
#ifdef KOKKOS_ENABLE_CUDA
# include "shoc_calc_shoc_varorcovar_impl.hpp"
# include "shoc_calc_shoc_vertflux_impl.hpp"
# include "shoc_diag_second_moments_srf_impl.hpp"
# include "shoc_diag_second_moments_ubycond_impl.hpp"
# include "shoc_update_host_dse_impl.hpp"
# include "shoc_pblintd_init_pot_impl.hpp"
#endif

#endif
