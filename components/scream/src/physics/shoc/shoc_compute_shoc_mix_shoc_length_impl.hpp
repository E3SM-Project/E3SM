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
  const Scalar&                tscale,
  const uview_1d<const Spack>& zt_grid,
  const Scalar&                l_inf,
  const uview_1d<Spack>&       shoc_mix)
{
  const Int nlev_pack = ekat::npack<Spack>(nlev);
  const auto maxlen = scream::shoc::Constants<Scalar>::maxlen;
  const auto length_fac = scream::shoc::Constants<Scalar>::length_fac;
  const auto vk = C::Karman;

  Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_pack), [&] (const Int& k) {
    Spack tkes(ekat::sqrt(tke(k)));

    Spack brunt2(0);
    brunt2.set(brunt(k) >= 0, brunt(k));

    Spack return_val(maxlen);
    Spack comp_val(sp(2.8284)*(ekat::sqrt(1/((1/(tscale*tkes*vk*zt_grid(k)))
                                            +(1/(tscale*tkes*l_inf))
                                            +sp(0.01)*(brunt2/tke(k)))))/length_fac);
    return_val.set(return_val > comp_val, comp_val);

    shoc_mix(k) = return_val;
  });
}

} // namespace shoc
} // namespace scream

#endif
