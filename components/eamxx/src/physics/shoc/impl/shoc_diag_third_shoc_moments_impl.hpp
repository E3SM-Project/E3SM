#ifndef SHOC_DIAG_THIRD_SHOC_MOMENTS_IMPL_HPP
#define SHOC_DIAG_THIRD_SHOC_MOMENTS_IMPL_HPP

#include "shoc_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace shoc {

/*
 * Implementation of shoc diag_third_shoc_moments. Clients should NOT
 * #include this file, but include shoc_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::diag_third_shoc_moments(
  const MemberType&            team,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Scalar&                c_diag_3rd_mom,
  const bool&                  shoc_1p5tke,
  const uview_1d<const Pack>& w_sec,
  const uview_1d<const Pack>& thl_sec,
  const uview_1d<const Pack>& wthl_sec,
  const uview_1d<const Pack>& isotropy,
  const uview_1d<const Pack>& brunt,
  const uview_1d<const Pack>& thetal,
  const uview_1d<const Pack>& tke,
  const uview_1d<const Pack>& dz_zt,
  const uview_1d<const Pack>& dz_zi,
  const uview_1d<const Pack>& zt_grid,
  const uview_1d<const Pack>& zi_grid,
  const Workspace&             workspace,
  const uview_1d<Pack>&       w3)
{
  // Define temporary variables
  uview_1d<Pack> isotropy_zi, w_sec_zi, brunt_zi, thetal_zi;
  workspace.template take_many_contiguous_unsafe<4>(
    {"isotropy_zi",  "w_sec_zi", "brunt_zi", "thetal_zi"},
    {&isotropy_zi, &w_sec_zi, &brunt_zi, &thetal_zi});

  // Constants
  const auto largeneg = SC::largeneg;
  const auto mintke = SC::mintke;

  // Interpolate variables onto the interface levels
  linear_interp(team,zt_grid,zi_grid,isotropy,isotropy_zi,nlev,nlevi,0);
  linear_interp(team,zt_grid,zi_grid,brunt,brunt_zi,nlev,nlevi,largeneg);
  linear_interp(team,zt_grid,zi_grid,w_sec,w_sec_zi,nlev,nlevi,sp(2.0/3.0)*mintke);
  linear_interp(team,zt_grid,zi_grid,thetal,thetal_zi,nlev,nlevi,0);
  team.team_barrier();

  // Diagnose the third moment of the vertical-velocity
  compute_diag_third_shoc_moment(team,nlev,nlevi,c_diag_3rd_mom,shoc_1p5tke,w_sec,thl_sec,wthl_sec,
                                 tke, dz_zt, dz_zi,isotropy_zi, brunt_zi,
                                 w_sec_zi,thetal_zi,w3);
  team.team_barrier();

  // Perform clipping to prevent unrealistically large values from occuring
  clipping_diag_third_shoc_moments(team,nlevi,w_sec_zi,w3);

  // Release temporary variables from the workspace
  workspace.template release_many_contiguous<4>(
    {&isotropy_zi, &w_sec_zi, &brunt_zi, &thetal_zi});
}

} // namespace shoc
} // namespace scream

#endif
