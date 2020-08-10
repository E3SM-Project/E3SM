#ifndef SHOC_CALC_SHOC_VERTFLUX_IMPL_HPP
#define SHOC_CALC_SHOC_VERTFLUX_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
void Functions<S,D>
::calc_shoc_vertflux(
  const MemberType& team,
  const Int& shcol, const Int& nlev, const Int& nlevi,
  const uview_2d<const Spack>& tkh_zi,
  const uview_2d<const Spack>& dz_zi,
  const uview_2d<const Spack>& invar,
  const uview_2d<Spack>& vertflux)
{
  for (Int k = 1; k < nlev; ++k) {
    const Int kt = k-1; // define upper grid point indicee
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, shcol), [&] (const Int& i) {
      const Real up_grid = invar(i,kt)[0];
      const Spack grid_dz = 1 / dz_zi(i,k); // vertical grid diff squared
      // Compute the vertical flux via downgradient diffusion
      vertflux(i,k) = -(tkh_zi(i,k) * grid_dz * (up_grid - invar(i,k)));
    });
  }
}

} // namespace shoc
} // namespace scream

#endif
