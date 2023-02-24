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
  const Scalar&                obklen,
  const Scalar&                pblh,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& shoc_mix,
  const uview_1d<const Spack>& sterm_zt,
  const uview_1d<const Spack>& isotropy,
  const uview_1d<const Spack>& tke,
  const uview_1d<Spack>&       tkh,
  const uview_1d<Spack>&       tk)
{
  // Parameters

  // Critical value of dimensionless Monin-Obukhov length,
  // for which diffusivities are no longer damped
  const Int zL_crit_val = 100;
  // Transition depth [m] above PBL top to allow
  // stability diffusivities
  const Int pbl_trans = 200;
  // Turbulent coefficients
  const Scalar Ckh = 0.1;
  const Scalar Ckm = 0.1;
  // Maximum eddy coefficients for stable PBL diffusivities
  const Scalar Ckh_s_max = 0.1;
  const Scalar Ckm_s_max = 0.1;
  // Minimum allowable value for stability diffusivities
  const Scalar Ckh_s_min = 0.1;
  const Scalar Ckm_s_min = 0.1;

  const auto s_zt_grid = ekat::scalarize(zt_grid);

  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    // Dimensionless Okukhov length considering only
    // the lowest model grid layer height to scale
    const auto z_over_L = s_zt_grid(nlev-1)/obklen;

    // Compute diffusivity coefficient as function of dimensionless Obukhov,
    // given a critical value
    const Scalar Ckh_s = ekat::impl::max(Ckh_s_min,
                                         ekat::impl::min(Ckh_s_max,
                                                         z_over_L/zL_crit_val));
    const Scalar Ckm_s = ekat::impl::max(Ckm_s_min,
                                         ekat::impl::min(Ckm_s_max,
                                                         z_over_L/zL_crit_val));

    // If surface layer is stable, based on near surface dimensionless Monin-Obukov
    // use modified coefficients of tkh and tk that are primarily based on shear
    // production and SHOC length scale, to promote mixing within the PBL and to a
    // height slighty above to ensure smooth transition.
    const Smask condition = (zt_grid(k) < pblh+pbl_trans) && (z_over_L > 0);
    tkh(k).set(condition, Ckh_s*ekat::square(shoc_mix(k))*ekat::sqrt(sterm_zt(k)));
    tk(k).set(condition,  Ckm_s*ekat::square(shoc_mix(k))*ekat::sqrt(sterm_zt(k)));

    // Default definition of eddy diffusivity for heat and momentum
    tkh(k).set(!condition, Ckh*isotropy(k)*tke(k));
    tk(k).set(!condition,  Ckm*isotropy(k)*tke(k));
  });
}

} // namespace shoc
} // namespace scream

#endif
