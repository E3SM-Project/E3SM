#ifndef SHOC_DIAG_SECOND_MOMENTS_IMPL_HPP
#define SHOC_DIAG_SECOND_MOMENTS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc diag_second_moments. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::diag_second_moments(
  const MemberType& team, const Int& nlev, const Int& nlevi,
  const uview_1d<const Spack>& thetal, const uview_1d<const Spack>& qw, const uview_1d<const Spack>& u_wind,
  const uview_1d<const Spack>& v_wind, const uview_1d<const Spack>& tke, const uview_1d<const Spack>& isotropy,
  const uview_1d<const Spack>& tkh, const uview_1d<const Spack>& tk, const uview_1d<const Spack>& dz_zi,
  const uview_1d<const Spack>& zt_grid, const uview_1d<const Spack>& zi_grid, const uview_1d<const Spack>& shoc_mix,
  const uview_1d<Spack>& isotropy_zi, const uview_1d<Spack>& tkh_zi, const uview_1d<Spack>& tk_zi,
  const uview_1d<Spack>& thl_sec, const uview_1d<Spack>& qw_sec, const uview_1d<Spack>& wthl_sec, const uview_1d<Spack>& wqw_sec,
  const uview_1d<Spack>& qwthl_sec, const uview_1d<Spack>& uw_sec, const uview_1d<Spack>& vw_sec, const uview_1d<Spack>& wtke_sec,
  const uview_1d<Spack>& w_sec)
{
  // Purpose of this subroutine is to diagnose the second
  //  order moments needed for the SHOC parameterization.
  //  Namely these are variances of thetal, qw, and vertical
  //  velocity.  In addition the vertical fluxes of thetal, qw,
  //  u, v, TKE, and tracers are computed here as well as the
  //  correlation of qw and thetal.

  const auto thl2tune = 1;

  // moisture variance
  const auto qw2tune = 1;

  // temp moisture covariance
  const auto qwthl2tune = 1;

  // vertical velocity variance
  const auto w2tune = 1;

  // Interpolate some variables from the midpoint grid to the interface grid
  linear_interp(team, zt_grid, zi_grid, isotropy, isotropy_zi, nlev, nlevi, 0);
  linear_interp(team, zt_grid, zi_grid, tkh,      tkh_zi,      nlev, nlevi, 0);
  linear_interp(team, zt_grid, zi_grid, tk,       tk_zi,       nlev, nlevi, 0);
  team.team_barrier();

  // Vertical velocity variance is assumed to be propotional to the TKE
  const Int nlev_pack = ekat::npack<Spack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    w_sec(k) = w2tune*(sp(2.)/sp(3.))*tke(k);
  });

  // Calculate the temperature variance
  calc_shoc_varorcovar(team, nlev, thl2tune, isotropy_zi, tkh_zi, dz_zi, thetal, thetal, thl_sec);

  // Calculate the moisture variance
  calc_shoc_varorcovar(team, nlev ,qw2tune, isotropy_zi, tkh_zi, dz_zi, qw, qw, qw_sec);

  // Calculate the temperature and moisture covariance
  calc_shoc_varorcovar(team, nlev, qwthl2tune, isotropy_zi, tkh_zi, dz_zi, thetal, qw, qwthl_sec);

  // Calculate vertical flux for heat
  calc_shoc_vertflux(team, nlev, tkh_zi, dz_zi, thetal, wthl_sec);

  // Calculate vertical flux for moisture
  calc_shoc_vertflux(team, nlev, tkh_zi, dz_zi, qw, wqw_sec);

  // Calculate vertical flux for TKE
  calc_shoc_vertflux(team, nlev, tkh_zi, dz_zi, tke, wtke_sec);

  // Calculate vertical flux for momentum (zonal wind)
  calc_shoc_vertflux(team, nlev, tk_zi, dz_zi, u_wind, uw_sec);

  // Calculate vertical flux for momentum (meridional wind)
  calc_shoc_vertflux(team, nlev, tk_zi, dz_zi, v_wind, vw_sec);
}

} // namespace shoc
} // namespace scream

#endif
