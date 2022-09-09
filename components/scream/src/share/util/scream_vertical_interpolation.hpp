#ifndef SCREAM_VERTICAL_INTERPOLATION_HPP
#define SCREAM_VERTICAL_INTERPOLATION_HPP

#include "share/scream_types.hpp"

#include "ekat/util/ekat_lin_interp.hpp"
#include "ekat/ekat_pack_utils.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {

//
// ------- Types --------
//

template <typename S>
using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;
using Spack = SmallPack<Real>;

template <typename S>
using SmallPack = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;
using Sbool = SmallPack<bool>;

using Smask = ekat::Mask<Spack::n>;

using KT = KokkosTypes<DefaultDevice>;
using ExeSpace = typename KT::ExeSpace;
using ESU = ekat::ExeSpaceUtils<ExeSpace>;
using LIV = ekat::LinInterp<Real,Spack::n>;

 
template <typename S>
using view_1d = typename KT::template view_1d<S>;
template <typename S>
using view_2d = typename KT::template view_2d<S>;

const Real masked_val = -9.99e17;
  
void perform_vertical_interpolation (const view_2d<const Spack>&,
				     const view_1d<const Spack>&,
				     const view_2d<const Spack>&,
				     const view_2d<      Spack>&,
				     const int&,
  				     const int&);

void perform_vertical_interpolation (const view_2d<const Spack>&,
				     const view_1d<const Spack>&,
				     const view_2d<const Spack>&,
				     const view_2d<      Spack>&,
				     const view_2d<      Smask>&,
				     const int&,
  				     const int&);

void perform_vertical_interpolation (const view_2d<const Spack>&,
				     const view_1d<const Spack>&,
				     const view_2d<const Spack>&,
				     const view_2d<      Spack>&,
				     const view_2d<      Smask>&,
				     const int&,
  				     const int&,
				     const Real&);

void perform_vertical_interpolation (const view_1d<const Spack>&,
				     const view_1d<const Spack>&,
				     const view_1d<const Spack>&,
				     const view_1d<      Spack>&,
				     const view_1d<      Smask>&,
				     const int&,
  				     const int&,
    				     const int&,
                                     const Real&,
      				     const LIV::MemberType&,
        			     const LIV&);
  
} // namespace scream

#endif // SCREAM_VERTICAL_INTERPOLATION_HPP
