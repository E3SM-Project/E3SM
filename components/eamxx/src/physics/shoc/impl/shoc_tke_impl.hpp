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
  const uview_1d<const Spack>& wthv_sec,
  const uview_1d<const Spack>& shoc_mix,
  const uview_1d<const Spack>& dz_zi,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& u_wind,
  const uview_1d<const Spack>& v_wind,
  const uview_1d<const Spack>& brunt,
  const Scalar&                obklen,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const Scalar&                pblh,
  const Workspace&             workspace,
  const uview_1d<Spack>&       tke,
  const uview_1d<Spack>&       tk,
  const uview_1d<Spack>&       tkh,
  const uview_1d<Spack>&       isotropy)
{
  // Define temporary variables
  uview_1d<Spack> sterm_zt, a_diss, sterm;
  workspace.template take_many_contiguous_unsafe<3>(
    {"sterm_zt", "a_diss", "sterm"},
    {&sterm_zt, &a_diss, &sterm});

  // Compute integrated column stability in lower troposphere
  Scalar brunt_int(0);
  integ_column_stability(team,nlev,dz_zt,pres,brunt,brunt_int);

  // Compute shear production term, which is on interface levels
  // This follows the methods of Bretheron and Park (2010)
  compute_shr_prod(team,nlevi,nlev,dz_zi,u_wind,v_wind,sterm);

  // Interpolate shear term from interface to thermo grid
  team.team_barrier();
  linear_interp(team,zi_grid,zt_grid,sterm,sterm_zt,nlevi,nlev,0);

  // Advance sgs TKE
  adv_sgs_tke(team,nlev,dtime,shoc_mix,wthv_sec,sterm_zt,tk,tke,a_diss);

  // Compute isotropic time scale [s]
  isotropic_ts(team,nlev,brunt_int,tke,a_diss,brunt,isotropy);

  // Compute eddy diffusivity for heat and momentum
  eddy_diffusivities(team,nlev,obklen,pblh,zt_grid,shoc_mix,sterm_zt,isotropy,tke,tkh,tk);

  // Release temporary variables from the workspace
  workspace.template release_many_contiguous<3>(
    {&sterm_zt, &a_diss, &sterm});
}

} // namespace shoc
} // namespace scream

#endif
