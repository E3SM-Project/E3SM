
#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace p3 {

template <>
void Functions<Real,DefaultDevice>
::rain_sedimentation_disp(
  const uview_2d<const Spack>& rho,
  const uview_2d<const Spack>& inv_rho,
  const uview_2d<const Spack>& rhofacr,
  const uview_2d<const Spack>& cld_frac_r,
  const uview_2d<const Spack>& inv_dz,
  const uview_2d<Spack>& qr_incld,
  const WorkspaceManager& workspace_mgr,
  const view_2d_table& vn_table_vals, const view_2d_table& vm_table_vals,
  const Int& nj, const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& inv_dt,
  const uview_2d<Spack>& qr,
  const uview_2d<Spack>& nr,
  const uview_2d<Spack>& nr_incld,
  const uview_2d<Spack>& mu_r,
  const uview_2d<Spack>& lamr,
  const uview_2d<Spack>& precip_liq_flux,
  const uview_2d<Spack>& qr_tend,
  const uview_2d<Spack>& nr_tend,
  const uview_1d<Scalar>& precip_liq_surf,
  const uview_1d<bool>& nucleationPossible,
  const uview_1d<bool>& hydrometeorsPresent,
  const P3Runtime& runtime_options)
{
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);
  // p3_rain_sedimentation loop
  Kokkos::parallel_for("p3_rain_sed_disp",
    policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();
    auto workspace = workspace_mgr.get_workspace(team);
    if (!(nucleationPossible(i) || hydrometeorsPresent(i))) {
      return;
    }

    // Rain sedimentation:  (adaptive substepping)
    rain_sedimentation(
      ekat::subview(rho, i), ekat::subview(inv_rho, i), ekat::subview(rhofacr, i), ekat::subview(cld_frac_r, i),
      ekat::subview(inv_dz, i), ekat::subview(qr_incld, i),
      team, workspace, vn_table_vals, vm_table_vals, nk, ktop, kbot, kdir, dt, inv_dt,
      ekat::subview(qr, i), ekat::subview(nr, i), ekat::subview(nr_incld, i), ekat::subview(mu_r, i),
      ekat::subview(lamr, i), ekat::subview(precip_liq_flux, i),
      ekat::subview(qr_tend, i), ekat::subview(nr_tend, i), precip_liq_surf(i), runtime_options);
  });

}
} // namespace p3
} // namespace scream
