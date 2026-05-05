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
  const Real& thl2tune, const Real& qw2tune, const Real& qwthl2tune, const Real& w2tune, const bool& shoc_1p5tke,
  const uview_1d<const Pack>& thetal, const uview_1d<const Pack>& qw, const uview_1d<const Pack>& u_wind,
  const uview_1d<const Pack>& v_wind, const uview_1d<const Pack>& tke, const uview_1d<const Pack>& isotropy,
  const uview_1d<const Pack>& tkh, const uview_1d<const Pack>& tk, const uview_1d<const Pack>& dz_zi,
  const uview_1d<const Pack>& zt_grid, const uview_1d<const Pack>& zi_grid, const uview_1d<const Pack>& shoc_mix,
  const uview_1d<Pack>& isotropy_zi, const uview_1d<Pack>& tkh_zi, const uview_1d<Pack>& tk_zi,
  const uview_1d<Pack>& thl_sec, const uview_1d<Pack>& qw_sec, const uview_1d<Pack>& wthl_sec, const uview_1d<Pack>& wqw_sec,
  const uview_1d<Pack>& qwthl_sec, const uview_1d<Pack>& uw_sec, const uview_1d<Pack>& vw_sec, const uview_1d<Pack>& wtke_sec,
  const uview_1d<Pack>& w_sec)
{
  // Purpose of this subroutine is to diagnose the second
  //  order moments needed for the SHOC parameterization.
  //  Namely these are variances of thetal, qw, and vertical
  //  velocity.  In addition the vertical fluxes of thetal, qw,
  //  u, v, TKE, and tracers are computed here as well as the
  //  correlation of qw and thetal.

  // Interpolate some variables from the midpoint grid to the interface grid
  linear_interp(team, zt_grid, zi_grid, isotropy, isotropy_zi, nlev, nlevi, 0);
  linear_interp(team, zt_grid, zi_grid, tkh,      tkh_zi,      nlev, nlevi, 0);
  linear_interp(team, zt_grid, zi_grid, tk,       tk_zi,       nlev, nlevi, 0);
  team.team_barrier();

  // Vertical velocity variance is assumed to be propotional to the TKE.
  //  If 1.5 TKE closure is activated then set to zero.
  const Int nlev_pack = ekat::npack<Pack>(nlev);
  Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
    w_sec(k) = shoc_1p5tke ? 0 : w2tune*(sp(2.)/sp(3.))*tke(k);
  });

  // For the following variances and covariance, if no SGS variability is desired then
  //  set these to zero.  Doing so, in conjuction with setting w3 and w2 (above) to zero
  //  will ensure that SHOC condensation reduces to an all-or-nothing scheme.
  if (shoc_1p5tke){
    const Int nlevi_pack = ekat::npack<Pack>(nlevi);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlevi_pack), [&] (const Int& k) {
      thl_sec(k) = 0;
      qw_sec(k) = 0;
      qwthl_sec(k) = 0;
    });
  }
  else{
    // Calculate the temperature variance
    calc_shoc_varorcovar(team, nlev, thl2tune, isotropy_zi, tkh_zi, dz_zi, thetal, thetal, thl_sec);

    // Calculate the moisture variance
    calc_shoc_varorcovar(team, nlev ,qw2tune, isotropy_zi, tkh_zi, dz_zi, qw, qw, qw_sec);

    // Calculate the temperature and moisture covariance
    calc_shoc_varorcovar(team, nlev, qwthl2tune, isotropy_zi, tkh_zi, dz_zi, thetal, qw, qwthl_sec);
  }

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
