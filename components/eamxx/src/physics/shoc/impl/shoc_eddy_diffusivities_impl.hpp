#ifndef SHOC_EDDY_DIFFUSIVITIES_IMPL_HPP
#define SHOC_EDDY_DIFFUSIVITIES_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc eddy_diffusivities. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::eddy_diffusivities(
  const MemberType&            team,
  const Int&                   nlev,
  const Scalar&                 Ckh,
  const Scalar&                 Ckm,
  const Scalar&                pblh,
  const uview_1d<const Pack>& zt_grid,
  const uview_1d<const Pack>& tabs,
  const uview_1d<const Pack>& shoc_mix,
  const uview_1d<const Pack>& sterm_zt,
  const uview_1d<const Pack>& isotropy,
  const uview_1d<const Pack>& tke,
  const uview_1d<Pack>&       tkh,
  const uview_1d<Pack>&       tk)
{
  // Parameters

  // Minimum absolute temperature [K] to which apply extra mixing
  const Int tabs_crit = 182;
  // Transition depth [m] above PBL top to allow
  // stability diffusivities
  const Int pbl_trans = 200;
  // Dddy coefficients for stable PBL diffusivities
  const Scalar Ckh_s = 0.1;
  const Scalar Ckm_s = 0.1;

  const auto s_tabs = ekat::scalarize(tabs);

  const Int nlev_pack = ekat::npack<Pack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    // If surface layer temperature is running away, apply extra mixing
    //   based on traditional stable PBL diffusivities that are not damped
    //   by stability functions.
    const Mask condition = (zt_grid(k) < pblh+pbl_trans) && (s_tabs(nlev-1) < tabs_crit);
    tkh(k).set(condition, Ckh_s*ekat::square(shoc_mix(k))*ekat::sqrt(sterm_zt(k)));
    tk(k).set(condition,  Ckm_s*ekat::square(shoc_mix(k))*ekat::sqrt(sterm_zt(k)));

    tkh(k).set(!condition, Ckh*isotropy(k)*tke(k));
    tk(k).set(!condition,  Ckm*isotropy(k)*tke(k));
  });
}

} // namespace shoc
} // namespace scream

#endif
