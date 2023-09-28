#include "shoc_functions.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace shoc {

template<>
void Functions<Real,DefaultDevice>
::update_prognostics_implicit_disp(
  const Int&                   shcol,
  const Int&                   nlev,
  const Int&                   nlevi,
  const Int&                   num_tracer,
  const Scalar&                dtime,
  const view_2d<const Spack>&  dz_zt,
  const view_2d<const Spack>&  dz_zi,
  const view_2d<const Spack>&  rho_zt,
  const view_2d<const Spack>&  zt_grid,
  const view_2d<const Spack>&  zi_grid,
  const view_2d<const Spack>&  tk,
  const view_2d<const Spack>&  tkh,
  const view_1d<const Scalar>& uw_sfc,
  const view_1d<const Scalar>& vw_sfc,
  const view_1d<const Scalar>& wthl_sfc,
  const view_1d<const Scalar>& wqw_sfc,
  const view_2d<const Spack>&  wtracer_sfc,
  const WorkspaceMgr&          workspace_mgr,
  const view_2d<Spack>&        thetal,
  const view_2d<Spack>&        qw,
  const view_3d<Spack>&        tracer,
  const view_2d<Spack>&        tke,
  const view_2d<Spack>&        u_wind,
  const view_2d<Spack>&        v_wind)
{
  using ExeSpace = typename KT::ExeSpace;

  const auto nlev_packs = ekat::npack<Spack>(nlev);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(shcol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    auto workspace  = workspace_mgr.get_workspace(team);

    update_prognostics_implicit(team, nlev, nlevi, num_tracer, dtime,
                                ekat::subview(dz_zt, i),
                                ekat::subview(dz_zi, i),
                                ekat::subview(rho_zt, i),
                                ekat::subview(zt_grid, i),
                                ekat::subview(zi_grid, i),
                                ekat::subview(tk, i),
                                ekat::subview(tkh, i),
                                uw_sfc(i),
                                vw_sfc(i),
                                wthl_sfc(i),
                                wqw_sfc(i),
                                ekat::subview(wtracer_sfc, i),
                                workspace,
                                ekat::subview(thetal, i),
                                ekat::subview(qw, i),
                                ekat::subview(tracer, i),
                                ekat::subview(tke, i),
                                ekat::subview(u_wind, i),
                                ekat::subview(v_wind, i));
  });
}

} // namespace shoc
} // namespace scream
