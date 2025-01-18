
#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 cloud sedimentation function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <>
void Functions<Real,DefaultDevice>
::cloud_sedimentation_disp(
    const uview_2d<Spack>& qc_incld,
    const uview_2d<const Spack>& rho,
    const uview_2d<const Spack>& inv_rho,
    const uview_2d<const Spack>& cld_frac_l,
    const uview_2d<const Spack>& acn,
    const uview_2d<const Spack>& inv_dz,
    const view_dnu_table& dnu,
    const WorkspaceManager& workspace_mgr,
    const Int& nj, const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& inv_dt, const bool& do_predict_nc,
    const uview_2d<Spack>& qc,
    const uview_2d<Spack>& nc,
    const uview_2d<Spack>& nc_incld,
    const uview_2d<Spack>& mu_c,
    const uview_2d<Spack>& lamc,
    const uview_2d<Spack>& qc_tend,
    const uview_2d<Spack>& nc_tend,
    const uview_1d<Scalar>& precip_liq_surf,
    const uview_1d<bool>& nucleationPossible,
    const uview_1d<bool>& hydrometeorsPresent)
{
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);
  // p3_cloud_sedimentation loop
  Kokkos::parallel_for(
    "p3_cloud_sedimentation",
    policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();
    auto workspace = workspace_mgr.get_workspace(team);
    if (!(nucleationPossible(i) || hydrometeorsPresent(i))) {
      return;
    }

    cloud_sedimentation(
      ekat::subview(qc_incld, i), ekat::subview(rho, i), ekat::subview(inv_rho, i), ekat::subview(cld_frac_l, i),
      ekat::subview(acn, i), ekat::subview(inv_dz, i), dnu, team, workspace,
      nk, ktop, kbot, kdir, dt, inv_dt, do_predict_nc,
      ekat::subview(qc, i), ekat::subview(nc, i), ekat::subview(nc_incld, i), ekat::subview(mu_c, i), ekat::subview(lamc, i), ekat::subview(qc_tend, i),
      ekat::subview(nc_tend, i),
      precip_liq_surf(i));
  });

}

} // namespace p3
} // namespace scream
