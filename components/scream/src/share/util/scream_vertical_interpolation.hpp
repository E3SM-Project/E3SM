#ifndef SCREAM_VERTICAL_INTERPOLATION_HPP
#define SCREAM_VERTICAL_INTERPOLATION_HPP

#include "share/scream_types.hpp"

#include "ekat/util/ekat_lin_interp.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace vinterp {
  
//
// ------- Types --------
//
template<typename T, int N>
using Pack = ekat::Pack<T, N>;

template <typename S>
using SmallPack = Pack<S,SCREAM_SMALL_PACK_SIZE>;
using Spack = SmallPack<Real>;

template<int N>
using Pmask = ekat::Mask<N>;

using Smask = ekat::Mask<Spack::n>;

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
  
template<typename T, int N> 
void perform_vertical_interpolation(
  const view_2d<Pack<T,N>>& x_src,
  const view_1d<Pack<T,N>>& x_tgt,
  const view_2d<Pack<T,N>>& input,
  const view_2d<Pack<T,N>>& output,
  const int& nlevs_src,
  const int& nlevs_tgt);

template<typename T, int N> 
void perform_vertical_interpolation(
  const view_2d<Pack<T,N>>& x_src,
  const view_1d<Pack<T,N>>& x_tgt,
  const view_2d<Pack<T,N>>& input,
  const view_2d<Pack<T,N>>& output,
  const view_2d<Pmask<N>>& mask,
  const int& nlevs_src,
  const int& nlevs_tgt);

template<typename T, int N> 
void perform_vertical_interpolation(
  const view_2d<Pack<T,N>>& x_src,
  const view_1d<Pack<T,N>>& x_tgt,
  const view_2d<Pack<T,N>>& input,
  const view_2d<Pack<T,N>>& output,
  const view_2d<Pmask<N>>& mask,
  const int& nlevs_src,
  const int& nlevs_tgt,
  const Real& masked_val);

template<typename T, int N> 
void perform_vertical_interpolation_impl_2d(
  const view_2d<Pack<T,N>>& x_src,
  const view_1d<Pack<T,N>>& x_tgt,
  const view_2d<Pack<T,N>>& input,
  const view_2d<Pack<T,N>>& output,
  const view_2d<Pmask<N>>& mask,
  const int& nlevs_src,
  const int& nlevs_tgt,
  const Real& masked_val);

template<typename T, int N> 
void perform_vertical_interpolation_impl_1d(
  const view_1d<Pack<T,N>>& x_src,
  const view_1d<Pack<T,N>>& x_tgt,
  const view_1d<Pack<T,N>>& input,
  const view_1d<Pack<T,N>>& output,
  const view_1d<Pmask<N>>& mask,
  const int& nlevs_src,
  const int& nlevs_tgt,
  const int& icol,
  const Real& masked_val,
  const MemberType& team,
  const LIV<T,N>& vert_interp);

} // namespace vinterp
} // namespace scream

#include "scream_vertical_interpolation_impl.hpp"

#endif // SCREAM_VERTICAL_INTERPOLATION_HPP
