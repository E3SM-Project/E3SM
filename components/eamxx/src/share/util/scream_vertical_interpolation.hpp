#ifndef SCREAM_VERTICAL_INTERPOLATION_HPP
#define SCREAM_VERTICAL_INTERPOLATION_HPP

#include "share/scream_types.hpp"

#include "ekat/util/ekat_lin_interp.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace vinterp {

/*This utility can perform vertical interpolation for a particular variable
 *(for instance a field) from src levels to target levels
 *The relevant function is perform_vertical_interpolation().
 *A user can decide to use the function without a masked value
 *in which case the default masked value (masked_val) will be used.
 *Also a user can decide to provide a view_2d mask which will be filled
 *based on whether each point is masked or not.
 *Masking occurs when a value is out-of-bounds (i.e. requires an 
 *extrapolation).
 *Most of the functions call assume you are providing a 2d view (or lambda)
 *for the src, and input, and a 1d view for the target.
 *There is a function, perform_vertical_interpolation_impl_1d which 
 *where what is provided is all 1d views (or lambdas). However, in this 
 *case the user must provide the team and ekat::LinInterp as input as well.
 */

// ------- Types --------

template<typename T, int P>
using Pack = ekat::Pack<T,P>;

template<int P>
using Mask = ekat::Mask<P>;

using KT = KokkosTypes<DefaultDevice>;
using ExeSpace = typename KT::ExeSpace;
using ESU = ekat::ExeSpaceUtils<ExeSpace>;

using MemberType = typename KT::MemberType;

template<typename T, int P>
using LIV = ekat::LinInterp<T,P>;
 
template <typename T>
using view_1d = typename KT::template view_1d<T>;

template <typename T, int P>
using view_Nd = typename KT::template view_ND<T,P>;
template <typename T, int P>
using view_Nd_host = typename KT::template view_ND<T,P>::HostMirror;
    
constexpr Real masked_val = -std::numeric_limits<Real>::max();

//This is the generic function call where the user provides source levels, 
//target levels, the input to interpolate, the output that has been interpolated,
//a 2d-mask from the user which is filled based on what values
//should have a mask (i.e. are out-of-bounds and require an extrapolation),
//the number of source and target levels,
//and a mask value provided by the user (msk_val)
//The 2d-mask should have the same dimensions as the output
template<typename T, int P> 
void perform_vertical_interpolation(
  const view_Nd<Pack<T,P>,2>& x_src,
  const view_1d<Pack<T,P>>& x_tgt,
  const view_Nd<Pack<T,P>,2>& input,
  const view_Nd<Pack<T,P>,2>& output,
  const view_Nd<Mask<P>,2>& mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val);

// Tjhis is the same as the generic function above, but in the case where the user
// does not want to output the masked array.
template<typename T, int P> 
void perform_vertical_interpolation(
  const view_Nd<Pack<T,P>,2>& x_src,
  const view_1d<Pack<T,P>>& x_tgt,
  const view_Nd<Pack<T,P>,2>& input,
  const view_Nd<Pack<T,P>,2>& output,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val);

//This function call does not have a mask value provided by user
//so uses the default masked value (masked_val) in the
//scream_vertical_interpolation.hpp file
template<typename T, int P> 
void perform_vertical_interpolation(
  const view_Nd<Pack<T,P>,2>& x_src,
  const view_1d<Pack<T,P>>& x_tgt,
  const view_Nd<Pack<T,P>,2>& input,
  const view_Nd<Pack<T,P>,2>& output,
  const view_Nd<Mask<P>,2>& mask,
  const int nlevs_src,
  const int nlevs_tgt);

//This function call does not have a mask value provided by user
//so uses the default masked value (masked_val) defined 
//in the scream_vertical_interpolation.hpp file.
//In addition, this function does not require a 2d mask from the user
//so one is setup since this is what perform_vertical_interpolation_impl_2d
//requires  
template<typename T, int P> 
void perform_vertical_interpolation(
  const view_Nd<Pack<T,P>,2>& x_src,
  const view_1d<Pack<T,P>>& x_tgt,
  const view_Nd<Pack<T,P>,2>& input,
  const view_Nd<Pack<T,P>,2>& output,
  const int nlevs_src,
  const int nlevs_tgt);

template<typename T, int P> 
KOKKOS_FUNCTION
void perform_vertical_interpolation_impl_1d(
  const view_1d<Pack<T,P>>& x_src,
  const view_1d<Pack<T,P>>& x_tgt,
  const view_1d<Pack<T,P>>& input,
  const view_1d<Pack<T,P>>& output,
  const view_1d<Mask<P>>& mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const int icol,
  const Real msk_val,
  const MemberType& team,
  const LIV<T,P>& vert_interp);

template<typename T, int P, int N, int M> 
void perform_vertical_interpolation_impl_Nd(
  const view_Nd<const Pack<T,P>,M>& x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<Pack<T,P>,N>&       output,
  const view_Nd<Mask<P>,N>&         mask,
  const int                         nlevs_src,
  const int                         nlevs_tgt,
  const Real                        msk_val);

template<typename T, int P> 
void apply_interpolation(
  const                        int  num_levs,
  const                          T  mask_val,
  const                   LIV<T,P>& vert_interp,
  const view_Nd<const Pack<T,P>,2>& x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,2>& input,
  const view_Nd<      Pack<T,P>,2>& output,
  const view_Nd<        Mask<P>,2>& mask);

template<typename T, int P> 
void apply_interpolation(
  const                        int  num_levs,
  const                          T  mask_val,
  const                   LIV<T,P>& vert_interp,
  const view_Nd<const Pack<T,P>,2>& x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,3>& input,
  const view_Nd<      Pack<T,P>,3>& output,
  const view_Nd<        Mask<P>,3>& mask);

template<typename T, int P> 
void apply_interpolation(
  const                        int  num_levs,
  const                          T  mask_val,
  const                   LIV<T,P>& vert_interp,
  const view_Nd<const Pack<T,P>,2>& x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,4>& input,
  const view_Nd<      Pack<T,P>,4>& output,
  const view_Nd<        Mask<P>,4>& mask);
  
template<typename T, int P>
KOKKOS_FUNCTION
void apply_masking(
  const MemberType& team,
  const T           mask_val,
  const T           min_val,
  const T           max_val,
  const view_1d<const Pack<T,P>>& x_tgt,
  const view_1d<      Pack<T,P>>& out,
  const view_1d<      Mask<P>>&   mask);

} // namespace vinterp
} // namespace scream

#include "scream_vertical_interpolation_impl.hpp"

#endif // SCREAM_VERTICAL_INTERPOLATION_HPP
