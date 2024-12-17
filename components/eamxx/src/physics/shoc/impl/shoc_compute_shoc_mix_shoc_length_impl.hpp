#ifndef SHOC_COMPUTE_SHOC_MIX_SHOC_LENGTH_IMPL_HPP
#define SHOC_COMPUTE_SHOC_MIX_SHOC_LENGTH_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::compute_shoc_mix_shoc_length(
  const MemberType&            team,
  const Int&                   nlev,
  const Scalar&                length_fac,
  const uview_1d<const Spack>& tke,
  const uview_1d<const Spack>& brunt,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& tk,  
  const Scalar&                l_inf,
  const uview_1d<Spack>&       shoc_mix)
{
  const Int nlev_pack = ekat::npack<Spack>(nlev);
  const auto maxlen = scream::shoc::Constants<Scalar>::maxlen;
  const auto vk = C::Karman;
  const bool tke_1p5_closure = scream::shoc::Constants<bool>::tke_1p5_closure;
  
  // Eddy turnover timescale
  const Scalar tscale = 400;

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    const Spack tkes = ekat::sqrt(tke(k));
    const Spack brunt2 = ekat::max(0, brunt(k));

    if (tke_1p5_closure){
       shoc_mix(k) = dz_zt(k); // may need to add possible stability correction
//       const auto stable_mask = brunt > 0;
//       if (stable_mask.any()){
//         shoc_mix(k) = ekat::min(dz_zt(k),ekat::max(0.1*dz_zt(k),ekat::sqrt(0.76
//	                        *tk(k)/0.1/ekat::sqrt(brunt+1.e-10))));
//       }
    }else{
       shoc_mix(k) = ekat::min(maxlen,
                            sp(2.8284)*(ekat::sqrt(1/((1/(tscale*tkes*vk*zt_grid(k)))
                            + (1/(tscale*tkes*l_inf))
                            + sp(0.01)*(brunt2/tke(k)))))/length_fac);
    }
  });
}

} // namespace shoc
} // namespace scream

#endif
