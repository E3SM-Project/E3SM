#ifndef SHOC_COMPUTE_L_INF_SHOC_LENGTH_HPP
#define SHOC_COMPUTE_L_INF_SHOC_LENGTH_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::compute_l_inf_shoc_length(
  const MemberType&            team,
  const Int&                   nlev,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& tke,
  Scalar&                      l_inf)
{
  using ExeSpaceUtils = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  // Compute numerator
  Scalar numer = ExeSpaceUtils::view_reduction(team,0,nlev,
                                [&] (const int k) -> Spack {
    return ekat::sqrt(tke(k))*zt_grid(k)*dz_zt(k);
  });
  team.team_barrier(); // see comment in shoc_energy_integrals_impl.hpp

  // Compute denominator
  Scalar denom = ExeSpaceUtils::view_reduction(team,0,nlev,
                                [&] (const int k) -> Spack {
    return ekat::sqrt(tke(k))*dz_zt(k);
  });
  team.team_barrier();

  // Set l_inf
  l_inf = sp(0.1)*(numer/denom);
}

} // namespace shoc
} // namespace scream

#endif
