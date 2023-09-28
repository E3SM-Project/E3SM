#ifndef SCREAM_VERTICAL_INTERPOLATION_HPP
#define SCREAM_VERTICAL_INTERPOLATION_HPP

#include "share/scream_types.hpp"

#include "ekat/util/ekat_lin_interp.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace vinterp {

/* This utility can perform vertical interpolation for a particular variable
 * (for instance a field) from src levels to target levels
 * The relevant function is perform_vertical_interpolation().
 * A user can decide to use the function without a masked value
 * in which case the default masked value (masked_val) will be used.
 * Also a user can decide to provide a view_Nd mask which will be filled
 * based on whether each point is masked or not.
 * Masking occurs when a value is out-of-bounds (i.e. requires an 
 * extrapolation).
 * The main function assumes that the user is providing as input:
 *   a 2D view for the source vertical levels to interpolate from
 *   a 1D view for the target vertical levels to interpolate onto
 *   a set of Nd views for the input data, output data and an optional mask
 *     where Nd is a N-dimensional view of dimension 2 or 3.
 * There is a function, perform_vertical_interpolation_impl_1d which 
 * is provided for 1d views (or lambdas). However, in this 
 * case the user must provide the team and ekat::LinInterp as input as well.
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

template <typename T>
using view_2d = typename KT::template view_2d<T>;

template <typename T>
using view_3d = typename KT::template view_3d<T>;

template <typename T, int N>
using view_Nd = typename KT::template view_ND<T,N>;
template <typename T, int N>
using view_Nd_host = typename KT::template view_ND<T,N>::HostMirror;
    
constexpr Real masked_val = -std::numeric_limits<Real>::max();

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val = masked_val);

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const view_Nd<      Mask<P>,N>&   mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val = masked_val);

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_2d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val = masked_val);

template<typename T, int P, int N> 
void perform_vertical_interpolation(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_2d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const view_Nd<      Mask<P>,N>&   mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const Real msk_val = masked_val);

/* ---------------------------------------------------------------------- 
 * Main interpolation routine that applies vertical interpolation to a
 * single vertical slice of data. 
 * ---------------------------------------------------------------------- */
template<typename T, int P> 
KOKKOS_FUNCTION
void apply_interpolation_impl_1d(
  const view_1d<const Pack<T,P>>& x_src,
  const view_1d<const Pack<T,P>>& x_tgt,
  const view_1d<const Pack<T,P>>& input,
  const view_1d<      Pack<T,P>>& output,
  const view_1d<      Mask<P>>&   mask,
  const int nlevs_src,
  const int nlevs_tgt,
  const int icol,
  const T msk_val,
  const MemberType& team,
  const LIV<T,P>& vert_interp);

/* ---------------------------------------------------------------------- 
 * Versions where x_tgt is a 2-D view
 * ---------------------------------------------------------------------- */
template<typename T, int P, int N> 
void perform_checks(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_2d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const int                         nlevs_src,
  const int                         nlevs_tgt);

template<typename T, int P> 
void apply_interpolation(
  const                      int  num_levs_src,
  const                      int  num_levs_tgt,
  const                        T  mask_val,
  const                 LIV<T,P>& vert_interp,
  const view_2d<const Pack<T,P>>& x_src,
  const view_2d<const Pack<T,P>>& x_tgt,
  const view_2d<const Pack<T,P>>& input,
  const view_2d<      Pack<T,P>>& output,
  const view_2d<        Mask<P>>& mask);

template<typename T, int P> 
void apply_interpolation(
  const                      int  num_levs_src,
  const                      int  num_levs_tgt,
  const                        T  mask_val,
  const                 LIV<T,P>& vert_interp,
  const view_2d<const Pack<T,P>>& x_src,
  const view_2d<const Pack<T,P>>& x_tgt,
  const view_3d<const Pack<T,P>>& input,
  const view_3d<      Pack<T,P>>& output,
  const view_3d<        Mask<P>>& mask);

/* ---------------------------------------------------------------------- 
 * Versions where x_tgt is a single 1-D vertical profile
 * ---------------------------------------------------------------------- */
template<typename T, int P, int N> 
void perform_checks(
  const view_2d<const Pack<T,P>>&   x_src,
  const view_1d<const Pack<T,P>>&   x_tgt,
  const view_Nd<const Pack<T,P>,N>& input,
  const view_Nd<      Pack<T,P>,N>& output,
  const int                         nlevs_src,
  const int                         nlevs_tgt);

template<typename T, int P> 
void apply_interpolation(
  const                      int  num_levs_src,
  const                      int  num_levs_tgt,
  const                        T  mask_val,
  const                 LIV<T,P>& vert_interp,
  const view_2d<const Pack<T,P>>& x_src,
  const view_1d<const Pack<T,P>>& x_tgt,
  const view_2d<const Pack<T,P>>& input,
  const view_2d<      Pack<T,P>>& output,
  const view_2d<        Mask<P>>& mask);

template<typename T, int P> 
void apply_interpolation(
  const                      int  num_levs_src,
  const                      int  num_levs_tgt,
  const                        T  mask_val,
  const                 LIV<T,P>& vert_interp,
  const view_2d<const Pack<T,P>>& x_src,
  const view_1d<const Pack<T,P>>& x_tgt,
  const view_3d<const Pack<T,P>>& input,
  const view_3d<      Pack<T,P>>& output,
  const view_3d<        Mask<P>>& mask);

// Helper function to allocate memory for an Nd mask on the fly.
template<int P, int N>
view_Nd<Mask<P>,N> allocate_mask(const std::vector<int>& extents);

} // namespace vinterp
} // namespace scream

#include "scream_vertical_interpolation_impl.hpp"

#endif // SCREAM_VERTICAL_INTERPOLATION_HPP
