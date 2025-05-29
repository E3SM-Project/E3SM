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
  const bool&                  shoc_1p5tke,
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

  // Eddy turnover timescale
  const Scalar tscale = 400;

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    const Spack tkes = ekat::sqrt(tke(k));
    const Spack brunt2 = ekat::max(0, brunt(k));

    if (shoc_1p5tke){
      // If 1.5 TKE closure then set length scale to vertical grid spacing for
	    //   cells with unstable brunt vaisalla frequency.  Otherwise, overwrite the length
	    //   scale in stable cells with the new definition.

	    // Search for stable cells
      const auto stable_mask = brunt(k) > 0;

      // To avoid FPE when calculating sqrt(brunt), set brunt_tmp=0 in the case brunt<1.
      Spack brunt_tmp(stable_mask, brunt(k));

      // Define length scale for stable cells
      const auto length_tmp = ekat::sqrt(0.76*tk(k)/0.1/ekat::sqrt(brunt_tmp + 1.e-10));
      // Limit the stability corrected length scale between 0.1*dz and dz
	    const auto limited_len = ekat::min(dz_zt(k),ekat::max(0.1*dz_zt(k),length_tmp));

      // Set length scale to vertical grid if unstable, otherwise the stability adjusted value.
      shoc_mix(k).set(stable_mask, limited_len, dz_zt(k));
    } else{
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
