#include "p3_functions.hpp" // for ETI only but harmless for GPU

#include <ekat_subview_utils.hpp>
#include <ekat_team_policy_utils.hpp>

namespace scream {
namespace p3 {

/*
 * Implementation of p3 main function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

#ifdef CLANGOPT_WORKAROUND
#pragma clang optimize off
#endif
template <>
void Functions<Real,DefaultDevice>
::p3_main_part2_disp(
  const Int& nj,
  const Int& nk,
  const Scalar& max_total_ni,
  const bool& predictNc,
  const bool& do_prescribed_CCN,
  const Scalar& dt,
  const Scalar& inv_dt,
  const uview_2d<const Pack>& hetfrz_immersion_nucleation_tend,
  const uview_2d<const Pack>& hetfrz_contact_nucleation_tend,
  const uview_2d<const Pack>& hetfrz_deposition_nucleation_tend,
  const view_dnu_table& dnu_table_vals,
  const view_ice_table& ice_table_vals,
  const view_collect_table& collect_table_vals,
  const view_2d_table& revap_table_vals,
  const uview_2d<const Pack>& pres,
  const uview_2d<const Pack>& dpres,
  const uview_2d<const Pack>& dz,
  const uview_2d<const Pack>& nc_nuceat_tend,
  const uview_2d<const Pack>& inv_exner,
  const uview_2d<const Pack>& exner,
  const uview_2d<const Pack>& inv_cld_frac_l,
  const uview_2d<const Pack>& inv_cld_frac_i,
  const uview_2d<const Pack>& inv_cld_frac_r,
  const uview_2d<const Pack>& ni_activated,
  const uview_2d<const Pack>& inv_qc_relvar,
  const uview_2d<const Pack>& cld_frac_i,
  const uview_2d<const Pack>& cld_frac_l,
  const uview_2d<const Pack>& cld_frac_r,
  const uview_2d<const Pack>& qv_prev,
  const uview_2d<const Pack>& t_prev,
//[shanyp 20260402
  const uview_2d<const Pack>& omega_mp,
//shanyp 20260402]
  const uview_2d<Pack>& T_atm,
  const uview_2d<Pack>& rho,
  const uview_2d<Pack>& inv_rho,
  const uview_2d<Pack>& qv_sat_l,
  const uview_2d<Pack>& qv_sat_i,
  const uview_2d<Pack>& qv_supersat_i,
  const uview_2d<Pack>& rhofacr,
  const uview_2d<Pack>& rhofaci,
  const uview_2d<Pack>& acn,
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
  const uview_2d<Pack>& qc_incld,
  const uview_2d<Pack>& qr_incld,
  const uview_2d<Pack>& qi_incld,
  const uview_2d<Pack>& qm_incld,
  const uview_2d<Pack>& nc_incld,
  const uview_2d<Pack>& nr_incld,
  const uview_2d<Pack>& ni_incld,
  const uview_2d<Pack>& bm_incld,
  const uview_2d<Pack>& mu_c,
  const uview_2d<Pack>& nu,
  const uview_2d<Pack>& lamc,
  const uview_2d<Pack>& cdist,
  const uview_2d<Pack>& cdist1,
  const uview_2d<Pack>& cdistr,
  const uview_2d<Pack>& mu_r,
  const uview_2d<Pack>& lamr,
  const uview_2d<Pack>& logn0r,
  const uview_2d<Pack>& qv2qi_depos_tend,
  const uview_2d<Pack>& precip_total_tend,
  const uview_2d<Pack>& nevapr,
  const uview_2d<Pack>& qr_evap_tend,
  const uview_2d<Pack>& vap_liq_exchange,
  const uview_2d<Pack>& vap_ice_exchange,
  const uview_2d<Pack>& liq_ice_exchange,
  const uview_2d<Pack>& qr2qv_evap,
  const uview_2d<Pack>& qi2qv_sublim,
  const uview_2d<Pack>& qc2qr_accret,
  const uview_2d<Pack>& qc2qr_autoconv,
  const uview_2d<Pack>& qv2qi_vapdep,
  const uview_2d<Pack>& qc2qi_berg,
  const uview_2d<Pack>& qc2qr_ice_shed,
  const uview_2d<Pack>& qc2qi_collect,
  const uview_2d<Pack>& qr2qi_collect,
  const uview_2d<Pack>& qc2qi_hetero_freeze,
  const uview_2d<Pack>& qr2qi_immers_freeze,
  const uview_2d<Pack>& qi2qr_melt,
  const uview_2d<Pack>& pratot,
  const uview_2d<Pack>& prctot,
  const uview_1d<bool>& nucleationPossible,
  const uview_1d<bool>& hydrometeorsPresent,
  const P3Runtime& runtime_options)
{
  using ExeSpace = typename KT::ExeSpace;
  using TPF      = ekat::TeamPolicyFactory<ExeSpace>;

  const Int nk_pack = ekat::npack<Pack>(nk);
  const auto policy = TPF::get_default_team_policy(nj, nk_pack);


  // p3_cloud_sedimentation loop
  Kokkos::parallel_for(
    "p3_main_part2_disp",
    policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();
    if (!(nucleationPossible(i) || hydrometeorsPresent(i))) {
      return;
    }

    // ------------------------------------------------------------------------------------------
    // main k-loop (for processes):
    p3_main_part2(
      team, nk_pack, max_total_ni, predictNc, do_prescribed_CCN, dt, inv_dt,
      ekat::subview(hetfrz_immersion_nucleation_tend, i),
      ekat::subview(hetfrz_contact_nucleation_tend, i),ekat::subview(hetfrz_deposition_nucleation_tend, i),
      dnu_table_vals, ice_table_vals, collect_table_vals, revap_table_vals,
      ekat::subview(pres, i), ekat::subview(dpres, i), ekat::subview(dz, i), ekat::subview(nc_nuceat_tend, i), ekat::subview(inv_exner, i),
      ekat::subview(exner, i), ekat::subview(inv_cld_frac_l, i), ekat::subview(inv_cld_frac_i, i), ekat::subview(inv_cld_frac_r, i),
      ekat::subview(ni_activated, i), ekat::subview(inv_qc_relvar, i), ekat::subview(cld_frac_i, i), ekat::subview(cld_frac_l, i),
//[shanyp 20260402
//      ekat::subview(cld_frac_r, i), ekat::subview(qv_prev, i), ekat::subview(t_prev, i), ekat::subview(T_atm, i), ekat::subview(rho, i),
      ekat::subview(cld_frac_r, i), ekat::subview(qv_prev, i), ekat::subview(t_prev, i), ekat::subview(omega_mp, i), ekat::subview(T_atm, i), ekat::subview(rho, i),
//shanyp 20260402]
      ekat::subview(inv_rho, i), ekat::subview(qv_sat_l, i), ekat::subview(qv_sat_i, i), ekat::subview(qv_supersat_i, i), ekat::subview(rhofacr, i),
      ekat::subview(rhofaci, i), ekat::subview(acn, i), ekat::subview(qv, i), ekat::subview(th_atm, i), ekat::subview(qc, i), ekat::subview(nc, i),
      ekat::subview(qr, i), ekat::subview(nr, i), ekat::subview(qi, i), ekat::subview(ni, i), ekat::subview(qm, i), ekat::subview(bm, i),
      ekat::subview(qc_incld, i),
      ekat::subview(qr_incld, i), ekat::subview(qi_incld, i), ekat::subview(qm_incld, i), ekat::subview(nc_incld, i), ekat::subview(nr_incld, i),
      ekat::subview(ni_incld, i), ekat::subview(bm_incld, i), ekat::subview(mu_c, i), ekat::subview(nu, i), ekat::subview(lamc, i), ekat::subview(cdist, i),
      ekat::subview(cdist1, i), ekat::subview(cdistr, i), ekat::subview(mu_r, i), ekat::subview(lamr, i), ekat::subview(logn0r, i),
      ekat::subview(qv2qi_depos_tend, i), ekat::subview(precip_total_tend, i),
      ekat::subview(nevapr, i), ekat::subview(qr_evap_tend, i), ekat::subview(vap_liq_exchange, i), ekat::subview(vap_ice_exchange, i), ekat::subview(liq_ice_exchange, i),
      ekat::subview(qr2qv_evap, i), ekat::subview(qi2qv_sublim, i), ekat::subview(qc2qr_accret, i), ekat::subview(qc2qr_autoconv, i),
      ekat::subview(qv2qi_vapdep, i), ekat::subview(qc2qi_berg, i), ekat::subview(qc2qr_ice_shed, i), ekat::subview(qc2qi_collect, i),
      ekat::subview(qr2qi_collect, i), ekat::subview(qc2qi_hetero_freeze, i), ekat::subview(qr2qi_immers_freeze, i),
      ekat::subview(qi2qr_melt, i),
      ekat::subview(pratot, i), ekat::subview(prctot, i), hydrometeorsPresent(i), nk, runtime_options);

    if (!hydrometeorsPresent(i)) return;
  });
}
#ifdef CLANGOPT_WORKAROUND
#pragma clang optimize on
#endif
} // namespace p3
} // namespace scream

