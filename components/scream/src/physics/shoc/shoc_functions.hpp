#ifndef SHOC_FUNCTIONS_HPP
#define SHOC_FUNCTIONS_HPP

#include "ekat/scream_types.hpp"
#include "ekat/scream_pack_kokkos.hpp"
#include "ekat/scream_workspace.hpp"
#include "physics_constants.hpp"

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

  using MemberType = typename KT::MemberType;

  using Workspace = typename WorkspaceManager<Spack, Device>::Workspace;

  //
  // --------- Functions ---------
  //
}; // struct Functions

} // namespace shoc
} // namespace scream

// If a GPU build, make all code available to the translation unit; otherwise,
// ETI is used.
#ifdef KOKKOS_ENABLE_CUDA
//# include "shoc_functions_math_impl.hpp"

#endif

#endif
