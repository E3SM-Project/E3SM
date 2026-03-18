#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "share/physics/physics_functions.hpp" // also for ETI not on GPUs
#include "share/physics/physics_saturation_impl.hpp"

#include <ekat_subview_utils.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {
namespace p3 {

/*
 * Implementation of p3 main function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <>
void Functions<Real,DefaultDevice>
::p3_main_part3_disp(
  const Int& nj,
  const Int& nk_pack,
  const Scalar& max_total_ni,
  const view_dnu_table& dnu_table_vals,
  const view_ice_table& ice_table_vals,
  const uview_2d<const Pack>& inv_exner,
  const uview_2d<const Pack>& cld_frac_l,
  const uview_2d<const Pack>& cld_frac_r,
  const uview_2d<const Pack>& cld_frac_i,
  const uview_2d<Pack>& rho,
  const uview_2d<Pack>& inv_rho,
  const uview_2d<Pack>& rhofaci,
  const uview_2d<Pack>& qv,
  const uview_2d<Pack>& th_atm,
  const uview_2d<Pack>& qc,
  const uview_2d<Pack>& nc,
  const uview_2d<Pack>& qr,
  const uview_2d<Pack>& nr,
  const uview_2d<Pack>& qi,
  const uview_2d<Pack>& ni,
  const uview_2d<Pack>& qm,
  const uview_2d<Pack>& bm,
  const uview_2d<Pack>& mu_c,
  const uview_2d<Pack>& nu,
  const uview_2d<Pack>& lamc,
  const uview_2d<Pack>& mu_r,
  const uview_2d<Pack>& lamr,
  const uview_2d<Pack>& vap_liq_exchange,
  const uview_2d<Pack>& ze_rain,
  const uview_2d<Pack>& ze_ice,
  const uview_2d<Pack>& diag_vm_qi,
  const uview_2d<Pack>& diag_eff_radius_qi,
  const uview_2d<Pack>& diag_diam_qi,
  const uview_2d<Pack>& rho_qi,
  const uview_2d<Pack>& diag_equiv_reflectivity,
  const uview_2d<Pack>& diag_eff_radius_qc,
  const uview_2d<Pack>& diag_eff_radius_qr,
  const uview_1d<bool>& nucleationPossible,
  const uview_1d<bool>& hydrometeorsPresent,
  const P3Runtime& runtime_options)
{
  using ExeSpace = typename KT::ExeSpace;
  using TPF      = ekat::TeamPolicyFactory<ExeSpace>;

  const auto policy = TPF::get_default_team_policy(nj, nk_pack);
  // p3_cloud_sedimentation loop
  Kokkos::parallel_for(
    "p3_main_part3_disp",
    policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();
    if (!(nucleationPossible(i) || hydrometeorsPresent(i))) {
      return;
    }

    //
    // final checks to ensure consistency of mass/number
    // and compute diagnostic fields for output
    //
    p3_main_part3(
      team, nk_pack, max_total_ni, dnu_table_vals, ice_table_vals, ekat::subview(inv_exner, i), ekat::subview(cld_frac_l, i), ekat::subview(cld_frac_r, i),
      ekat::subview(cld_frac_i, i), ekat::subview(rho, i), ekat::subview(inv_rho, i), ekat::subview(rhofaci, i), ekat::subview(qv, i),
      ekat::subview(th_atm, i), ekat::subview(qc, i), ekat::subview(nc, i), ekat::subview(qr, i), ekat::subview(nr, i), ekat::subview(qi, i),
      ekat::subview(ni, i), ekat::subview(qm, i), ekat::subview(bm, i),
      ekat::subview(mu_c, i), ekat::subview(nu, i), ekat::subview(lamc, i), ekat::subview(mu_r, i), ekat::subview(lamr, i),
      ekat::subview(vap_liq_exchange, i), ekat::subview(ze_rain, i), ekat::subview(ze_ice, i), ekat::subview(diag_vm_qi, i), ekat::subview(diag_eff_radius_qi, i),
      ekat::subview(diag_diam_qi, i), ekat::subview(rho_qi, i), ekat::subview(diag_equiv_reflectivity, i), ekat::subview(diag_eff_radius_qc, i),
      ekat::subview(diag_eff_radius_qr, i), runtime_options);

  });
}

} // namespace p3
} // namespace scream

