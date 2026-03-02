#ifndef CLD_FRAC_FUNCTIONS_HPP
#define CLD_FRAC_FUNCTIONS_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_workspace.hpp>

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

  using Pack = ekat::Pack<Scalar,SCREAM_PACK_SIZE>;

  using Mask = ekat::Mask<Pack::n>;

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
    const Real ice_threshold,
    const Real ice_4out_threshold,
    const view_2d<const Pack>& qi, 
    const view_2d<const Pack>& liq_cld_frac, 
    const view_2d<Pack>& ice_cld_frac, 
    const view_2d<Pack>& tot_cld_frac,
    const view_2d<Pack>& ice_cld_frac_4out, 
    const view_2d<Pack>& tot_cld_frac_4out);

  KOKKOS_FUNCTION
  static void calc_icefrac( 
    const MemberType& team,
    const Int& nk,
    const Real& threshold,
    const uview_1d<const Pack>& qi,
    const uview_1d<Pack>&       ice_cld_frac);

  KOKKOS_FUNCTION
  static void calc_totalfrac( 
    const MemberType& team,
    const Int& nk,
    const uview_1d<const Pack>& liq_cld_frac,
    const uview_1d<const Pack>& ice_cld_frac,
    const uview_1d<Pack>&       tot_cld_frac);

}; // struct Functions

} // namespace cld_fraction
} // namespace scream

#endif // CLD_FRAC_FUNCTIONS_HPP
