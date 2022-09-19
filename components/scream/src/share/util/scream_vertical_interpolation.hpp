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
template <typename S>
using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;
using Spack = SmallPack<Real>;

using Smask = ekat::Mask<Spack::n>;

using KT = KokkosTypes<DefaultDevice>;
using ExeSpace = typename KT::ExeSpace;
using ESU = ekat::ExeSpaceUtils<ExeSpace>;
using LIV = ekat::LinInterp<Real,Spack::n>;
 
template <typename S>
using view_1d = typename KT::template view_1d<S>;
template <typename S>
using view_2d = typename KT::template view_2d<S>;
    
const Real masked_val = -std::numeric_limits<Real>::max();
  
void perform_vertical_interpolation (const view_2d<const Spack>& x_src,
				     const view_1d<const Spack>& x_tgt,
				     const view_2d<const Spack>& input,
				     const view_2d<      Spack>& output,
				     const int& nlevs_src,
  				     const int& nlevs_tgt);

void perform_vertical_interpolation (const view_2d<const Spack>& x_src,
				     const view_1d<const Spack>& x_tgt,
				     const view_2d<const Spack>& input,
				     const view_2d<      Spack>& output,
				     const view_2d<      Smask>& mask,
				     const int& nlevs_src,
  				     const int& nlevs_tgt);

void perform_vertical_interpolation (const view_2d<const Spack>& x_src,
				     const view_1d<const Spack>& x_tgt,
				     const view_2d<const Spack>& input,
				     const view_2d<      Spack>& output,
				     const view_2d<      Smask>& mask,
				     const int& nlevs_src,
  				     const int& nlevs_tgt,
				     const Real& masked_val);

void perform_vertical_interpolation (const view_1d<const Spack>& x_src,
				     const view_1d<const Spack>& x_tgt,
				     const view_1d<const Spack>& input,
				     const view_1d<      Spack>& output,
				     const view_1d<      Smask>& mask,
				     const int& nlevs_src,
  				     const int& nlevs_tgt,
    				     const int& icol,
                                     const Real& masked_val,
      				     const LIV::MemberType& team,
        			     const LIV& vert_interp);

} // namespace vinterp
} // namespace scream

#endif // SCREAM_VERTICAL_INTERPOLATION_HPP
