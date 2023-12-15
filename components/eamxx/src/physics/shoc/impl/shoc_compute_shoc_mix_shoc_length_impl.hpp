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
  const uview_1d<const Spack>& tke,
  const uview_1d<const Spack>& brunt,
  const uview_1d<const Spack>& zt_grid,
  const Scalar&                l_inf,
  const uview_1d<Spack>&       shoc_mix)
{
  const Int nlev_pack = ekat::npack<Spack>(nlev);
  const auto maxlen = scream::shoc::Constants<Scalar>::maxlen;
  const auto length_fac = scream::shoc::Constants<Scalar>::length_fac;
  const auto vk = C::Karman;
  
  // Eddy turnover timescale
  const Scalar tscale = 400;

  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    const Spack tkes = ekat::sqrt(tke(k));
    const Spack brunt2 = ekat::max(0, brunt(k));

    shoc_mix(k) = ekat::min(maxlen,
                            sp(2.8284)*(ekat::sqrt(1/((1/(tscale*tkes*vk*zt_grid(k)))
                            + (1/(tscale*tkes*l_inf))
                            + sp(0.01)*(brunt2/tke(k)))))/length_fac);
  });
}

} // namespace shoc
} // namespace scream

#endif
