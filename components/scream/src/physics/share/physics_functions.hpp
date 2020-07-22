#ifndef PHYSICS_FUNCTIONS_HPP
#define PHYSICS_FUNCTIONS_HPP

#include "ekat/scream_types.hpp"
#include "ekat/scream_pack_kokkos.hpp"
#include "ekat/scream_workspace.hpp"
#include "physics_constants.hpp"

namespace scream {
namespace physics {

/*
 * Functions is a stateless struct used to encapsulate a
 * number of functions for p3. We use the ETI pattern for
 * these functions.
 *
 * P3 assumptions:
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
  using BigPack = scream::pack::BigPack<S>;
  template <typename S>
  using SmallPack = scream::pack::SmallPack<S>;
  using IntSmallPack = scream::pack::IntSmallPack;

  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  template <typename S>
  using Mask = scream::pack::Mask<BigPack<S>::n>;

  template <typename S>
  using SmallMask = scream::pack::Mask<SmallPack<S>::n>;

  using Smask = SmallMask<Scalar>;

  using KT = KokkosTypes<Device>;

  using C = scream::physics::Constants<Scalar>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;

  template <typename S, int N>
  using view_1d_ptr_array = typename KT::template view_1d_ptr_carray<S, N>;

  template <typename S>
  using uview_1d = typename ko::template Unmanaged<view_1d<S> >;
  template <typename S>
  using uview_2d = typename ko::template Unmanaged<view_2d<S> >;

  using MemberType = typename KT::MemberType;

  using Workspace = typename WorkspaceManager<Spack, Device>::Workspace;

  //
  // --------- Functions ---------
  //

  //  compute saturation vapor pressure
  //  polysvp1 returned in units of pa.
  //  t is input in units of k.
  //  ice refers to saturation with respect to liquid (false) or ice (true)
  KOKKOS_FUNCTION
  static Spack polysvp1(const Spack& t, const bool ice);

  //  compute saturation vapor pressure using Murphy and Koop(2005) formulation
  //  MurphyKoop_svp returned in units of pa.
  //  t is input in units of k.
  //  ice refers to saturation with respect to liquid (false) or ice (true)
  KOKKOS_FUNCTION
  static Spack MurphyKoop_svp(const Spack& t, const bool ice);

  // Calls a function to obtain the saturation vapor pressure, and then computes
  // and returns the saturation mixing ratio, with respect to either liquid or ice,
  // depending on value of 'ice'
  KOKKOS_FUNCTION
  static Spack qv_sat(const Spack& t_atm, const Spack& p_atm, const bool ice, const int& func_idx = 1);
};

} // namespace physics
} // namespace scream

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
#ifdef KOKKOS_ENABLE_CUDA
# include "physics_saturation_impl.hpp"
#endif

#endif
