#ifndef CLDFRAC_FUNCTIONS_HPP
#define CLDFRAC_FUNCTIONS_HPP

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace cldfrac {

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
  using MemberType = typename KT::MemberType;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S> >;

  static void cldfraction_main(
    const Int nj, 
    const Int nk,
    const view_2d<const Pack>& qi, 
    const view_2d<const Pack>& alst, 
    const view_2d<Pack>& aist, 
    const view_2d<Pack>& ast);

  KOKKOS_FUNCTION
  static void cldfraction_calc_icefrac( 
    const MemberType& team,
    const Int& nk,
    const uview_1d<const Spack>& qi,
    const uview_1d<Spack>&       aist);

  KOKKOS_FUNCTION
  static void cldfraction_calc_totalfrac( 
    const MemberType& team,
    const Int& nk,
    const uview_1d<const Spack>& alst,
    const uview_1d<const Spack>& aist,
    const uview_1d<Spack>&       ast);

}; // struct Functions

} // namespace cldfrac
} // namespace scream

#endif // CLDFRAC_FUNCTIONS_HPP
