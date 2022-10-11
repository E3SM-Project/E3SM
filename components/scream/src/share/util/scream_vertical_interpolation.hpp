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

template<typename T, int N>
using Pack = ekat::Pack<T, N>;

template<int N>
using Mask = ekat::Mask<N>;

using KT = KokkosTypes<DefaultDevice>;
using ExeSpace = typename KT::ExeSpace;
using ESU = ekat::ExeSpaceUtils<ExeSpace>;

using MemberType = typename KT::MemberType;

template<typename T, int N>
using LIV = ekat::LinInterp<T,N>;
 
template <typename S>
using view_1d = typename KT::template view_1d<S>;
template <typename S>
using view_2d = typename KT::template view_2d<S>;
    
const Real masked_val = -std::numeric_limits<Real>::max();

//This function call does not have a mask value provided by user
//so uses the default masked value (masked_val) defined 
//in the scream_vertical_interpolation.hpp file.
//In addition, this function does not require a 2d mask from the user
//so one is setup since this is what perform_vertical_interpolation_impl_2d
//requires  
template<typename Src, typename Tgt, typename Input, typename T, int N> 
void perform_vertical_interpolation(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_2d<Pack<T,N>>& output,
  const int nlevs_src,
  const int nlevs_tgt);

//This function call does not have a mask value provided by user
//so uses the default masked value (masked_val) in the
//scream_vertical_interpolation.hpp file
template<typename Src, typename Tgt, typename Input, typename T, int N> 
void perform_vertical_interpolation(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_2d<Pack<T,N>>& output,
  const view_2d<Mask<N>>& mask,
  const int nlevs_src,
  const int nlevs_tgt);

//This function call has a mask value provided by the user (msk_val)
//which is then used as the mask value instead of the default
//Also this function requires a 2d-mask from the user which is returned
//based on what values were required to be masked
template<typename Src, typename Tgt, typename Input, typename T, int N> 
void perform_vertical_interpolation(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_2d<Pack<T,N>>& output,
  const view_2d<Mask<N>>& mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real& msk_val);

template<typename Src, typename Tgt, typename Input, typename T, int N> 
void perform_vertical_interpolation_impl_2d(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_2d<Pack<T,N>>& output,
  const view_2d<Mask<N>>& mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real& msk_val);

template<typename Src, typename Tgt, typename Input, typename T, int N> 
void perform_vertical_interpolation_impl_1d(
  const Src& x_src,
  const Tgt& x_tgt,
  const Input& input,
  const view_1d<Pack<T,N>>& output,
  const view_1d<Mask<N>>& mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const int icol,
  const Real& msk_val,
  const MemberType& team,
  const LIV<T,N>& vert_interp);

} // namespace vinterp
} // namespace scream

#include "scream_vertical_interpolation_impl.hpp"

#endif // SCREAM_VERTICAL_INTERPOLATION_HPP
