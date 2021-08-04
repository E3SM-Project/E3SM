#ifndef SPA_FUNCTIONS_HPP
#define SPA_FUNCTIONS_HPP

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace spa {

template <typename ScalarT, typename DeviceT>
struct SPAFunctions
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

  using Pack = BigPack<Scalar>;
  using Spack = SmallPack<Scalar>;

  using KT = KokkosTypes<Device>;
  using MemberType = typename KT::MemberType;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  KOKKOS_FUNCTION
  static void reconstruct_pressure_profile(
    const Int ncols,
    const Int nlevs,
    const view_1d<const Spack>& hya,
    const view_1d<const Spack>& hyb,
    const view_1d<const Spack>& PS,
    const view_2d<Spack>&       pres);

  KOKKOS_FUNCTION
  static void aero_vertical_remap(
    const Int ncols,
    const Int nlevs_src,
    const Int nlevs_tgt,
    const view_2d<Spack>& pres_src,
    const view_2d<Spack>& pres_tgt,
    const view_2d<Spack>& aero_src,
    const view_2d<Spack>& aero_tgt);  //TODO, fix const for these when API in EKAT is fixed.

}; // struct Functions

} // namespace spa 
} // namespace scream

#endif // SPA_FUNCTIONS_HPP
