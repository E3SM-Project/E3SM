#ifndef SHOC_TKE_IMPL_HPP
#define SHOC_TKE_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc shoc_tke. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

/*
 * Thus function advances the SGS
 * TKE equation due to shear production, buoyant
 * production, and dissipation processes.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::shoc_tke(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Scalar&                dtime,
  const Scalar&                lambda_low,
  const Scalar&                lambda_high,
  const Scalar&                lambda_slope,
  const Scalar&                lambda_thresh,
  const Scalar&                Ckh,
  const Scalar&                Ckm,
  const bool&                  shoc_1p5tke,
  const bool&                  do_3d_turb,
  const uview_1d<const Pack>& wthv_sec,
  const uview_1d<const Pack>& shear_strain3d,
  const uview_1d<const Pack>& shoc_mix,
  const uview_1d<const Pack>& dz_zi,
  const uview_1d<const Pack>& dz_zt,
  const uview_1d<const Pack>& pres,
  const uview_1d<const Pack>& tabs,
  const uview_1d<const Pack>& u_wind,
  const uview_1d<const Pack>& v_wind,
  const uview_1d<const Pack>& brunt,
  const uview_1d<const Pack>& zt_grid,
  const uview_1d<const Pack>& zi_grid,
  const Scalar&                pblh,
  const Workspace&             workspace,
  const uview_1d<Pack>&       tke,
  const uview_1d<Pack>&       tk,
  const uview_1d<Pack>&       tkh,
  const uview_1d<Pack>&       isotropy)
{
  // Define temporary variables
  uview_1d<Pack> sterm_zt, a_diss, sterm;
  workspace.template take_many_contiguous_unsafe<3>(
    {"sterm_zt", "a_diss", "sterm"},
    {&sterm_zt, &a_diss, &sterm});

  // Compute integrated column stability in lower troposphere
  Scalar brunt_int(0);
  integ_column_stability(team,nlev,dz_zt,pres,brunt,brunt_int);

  // If not using 3d turbulence then use the default 1D calculation for shear production
  if (!do_3d_turb){
    // Compute shear production term, which is on interface levels
    // This follows the methods of Bretheron and Park (2010)
    compute_shr_prod(team,nlevi,nlev,dz_zi,u_wind,v_wind,sterm);

    // Interpolate shear term from interface to thermo grid
    team.team_barrier();
    linear_interp(team,zi_grid,zt_grid,sterm,sterm_zt,nlevi,nlev,0);
  } else {
    // eddy_diffusivities still needs a midpoint shear magnitude for the
    // cold-surface fallback path. In 3D mode, use the SHOC-grid strain term
    // instead of leaving the old 1D shear workspace uninitialized.
    static constexpr Scalar Ck_sh = 0.1;
    const Int nlev_pack = ekat::npack<Pack>(nlev);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nlev_pack), [&] (const Int& k) {
      sterm_zt(k) = Ck_sh*ekat::max(Pack(0),shear_strain3d(k));
    });
    team.team_barrier();
  }

  // Advance sgs TKE
  adv_sgs_tke(team,nlev,dtime,shoc_1p5tke,do_3d_turb,shoc_mix,wthv_sec,sterm_zt,tk,brunt,shear_strain3d,tke,a_diss);

  // Compute isotropic time scale [s]
  isotropic_ts(team,nlev,lambda_low,lambda_high,lambda_slope,lambda_thresh,brunt_int,tke,a_diss,brunt,isotropy);

  // Compute eddy diffusivity for heat and momentum
  eddy_diffusivities(team,nlev,shoc_1p5tke,Ckh,Ckm,pblh,zt_grid,tabs,shoc_mix,sterm_zt,isotropy,tke,tkh,tk);

  // Release temporary variables from the workspace
  workspace.template release_many_contiguous<3>(
    {&sterm_zt, &a_diss, &sterm});
}

} // namespace shoc
} // namespace scream

#endif
