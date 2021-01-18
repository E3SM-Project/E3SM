#ifndef SHOC_LENGTH_IMPL_HPP
#define SHOC_LENGTH_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::shoc_length(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Scalar&                host_dx,
  const Scalar&                host_dy,
  const Scalar&                pblh,
  const uview_1d<const Spack>& tke,
  const uview_1d<const Spack>& zt_grid,
  const uview_1d<const Spack>& zi_grid,
  const uview_1d<const Spack>& dz_zt,
  const uview_1d<const Spack>& wthv_sec,
  const uview_1d<const Spack>& thv,
  const uview_1d<Spack>&       thv_zi,
  const uview_1d<Spack>&       brunt,
  const uview_1d<Spack>&       shoc_mix)
{
  linear_interp(team,zt_grid,zi_grid,thv,thv_zi,nlev,nlevi,0);
  team.team_barrier();

  compute_brunt_shoc_length(team,nlev,nlevi,dz_zt,thv,thv_zi,brunt);
  team.team_barrier();

  Scalar l_inf = 0;
  compute_l_inf_shoc_length(team,nlev,zt_grid,dz_zt,tke,l_inf);

  Scalar conv_vel = 0;
  compute_conv_vel_shoc_length(team,nlev,pblh,zt_grid,dz_zt,thv,wthv_sec,conv_vel);
  team.team_barrier();

  Scalar tscale = 0;
  compute_conv_time_shoc_length(pblh,conv_vel,tscale);

  compute_shoc_mix_shoc_length(team,nlev,tke,brunt,tscale,zt_grid,l_inf,shoc_mix);
  team.team_barrier();

  check_length_scale_shoc_length(team,nlev,host_dx,host_dy,shoc_mix);
}

} // namespace shoc
} // namespace scream

#endif
