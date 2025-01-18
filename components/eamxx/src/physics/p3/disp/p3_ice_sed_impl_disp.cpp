
#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace p3 {

template <>
void Functions<Real,DefaultDevice>
::ice_sedimentation_disp(
  const uview_2d<const Spack>& rho,
  const uview_2d<const Spack>& inv_rho,
  const uview_2d<const Spack>& rhofaci,
  const uview_2d<const Spack>& cld_frac_i,
  const uview_2d<const Spack>& inv_dz,
  const WorkspaceManager& workspace_mgr,
  const Int& nj, const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& inv_dt,
  const uview_2d<Spack>& qi,
  const uview_2d<Spack>& qi_incld,
  const uview_2d<Spack>& ni,
  const uview_2d<Spack>& ni_incld,
  const uview_2d<Spack>& qm,
  const uview_2d<Spack>& qm_incld,
  const uview_2d<Spack>& bm,
  const uview_2d<Spack>& bm_incld,
  const uview_2d<Spack>& qi_tend,
  const uview_2d<Spack>& ni_tend,
  const view_ice_table& ice_table_vals,
  const uview_1d<Scalar>& precip_ice_surf,
  const uview_1d<bool>& nucleationPossible,
  const uview_1d<bool>& hydrometeorsPresent,
  const P3Runtime& runtime_options)
{
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);
  // p3_ice_sedimentation loop
  Kokkos::parallel_for("p3_ice_sedimentation",
    policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();
    if (!(nucleationPossible(i) || hydrometeorsPresent(i))) {
      return;
    }
    auto workspace = workspace_mgr.get_workspace(team);

    // Ice sedimentation:  (adaptive substepping)
    ice_sedimentation(
      ekat::subview(rho, i), ekat::subview(inv_rho, i), ekat::subview(rhofaci, i), ekat::subview(cld_frac_i, i),
      ekat::subview(inv_dz, i), team, workspace, nk, ktop, kbot, kdir, dt, inv_dt,
      ekat::subview(qi, i), ekat::subview(qi_incld, i), ekat::subview(ni, i), ekat::subview(ni_incld, i),
      ekat::subview(qm, i), ekat::subview(qm_incld, i), ekat::subview(bm, i), ekat::subview(bm_incld, i), ekat::subview(qi_tend, i), ekat::subview(ni_tend, i),
      ice_table_vals, precip_ice_surf(i), runtime_options);

 });
}

template <>
void Functions<Real,DefaultDevice>
::homogeneous_freezing_disp(
  const uview_2d<const Spack>& T_atm,
  const uview_2d<const Spack>& inv_exner,
  const Int& nj, const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir,
  const uview_2d<Spack>& qc,
  const uview_2d<Spack>& nc,
  const uview_2d<Spack>& qr,
  const uview_2d<Spack>& nr,
  const uview_2d<Spack>& qi,
  const uview_2d<Spack>& ni,
  const uview_2d<Spack>& qm,
  const uview_2d<Spack>& bm,
  const uview_2d<Spack>& th_atm,
  const uview_1d<bool>& nucleationPossible,
  const uview_1d<bool>& hydrometeorsPresent)
{
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);
  // p3_cloud_sedimentation loop
  Kokkos::parallel_for(
    "p3_homogeneous",
    policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();
    if (!(nucleationPossible(i) || hydrometeorsPresent(i))) {
      return;
    }

    // homogeneous freezing of cloud and rain
    homogeneous_freezing(
      ekat::subview(T_atm, i), ekat::subview(inv_exner, i), team, nk, ktop, kbot, kdir,
      ekat::subview(qc, i), ekat::subview(nc, i), ekat::subview(qr, i), ekat::subview(nr, i), ekat::subview(qi, i),
      ekat::subview(ni, i), ekat::subview(qm, i), ekat::subview(bm, i), ekat::subview(th_atm, i));

  });
}

} // namespace p3
} // namespace scream

