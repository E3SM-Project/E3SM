#ifndef CLD_FRAC_FUNCTIONS_HPP
#define CLD_FRAC_FUNCTIONS_HPP

#include "share/scream_types.hpp"

#include "ekat/ekat_pack_kokkos.hpp"
#include "ekat/ekat_workspace.hpp"

namespace scream {
namespace cld_fraction {

template <typename ScalarT, typename DeviceT>
struct CldFractionFunctions
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

  static void main(
    const Int nj, 
    const Int nk,
    const view_2d<const Pack>& qi, 
    const view_2d<const Pack>& liq_cld_frac, 
    const view_2d<Pack>& ice_cld_frac, 
    const view_2d<Pack>& tot_cld_frac);

  KOKKOS_FUNCTION
  static void calc_icefrac( 
    const MemberType& team,
    const Int& nk,
    const uview_1d<const Spack>& qi,
    const uview_1d<Spack>&       ice_cld_frac);

  KOKKOS_FUNCTION
  static void calc_totalfrac( 
    const MemberType& team,
    const Int& nk,
    const uview_1d<const Spack>& liq_cld_frac,
    const uview_1d<const Spack>& ice_cld_frac,
    const uview_1d<Spack>&       tot_cld_frac);

}; // struct Functions

} // namespace cld_fraction
} // namespace scream

#endif // CLD_FRAC_FUNCTIONS_HPP
