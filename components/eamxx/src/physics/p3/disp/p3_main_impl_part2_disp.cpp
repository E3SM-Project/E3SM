
#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 main function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

#ifdef SCREAM_SYSTEM_WORKAROUND_P3_PART2
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
  const uview_2d<const Spack>& hetfrz_immersion_nucleation_tend,
  const uview_2d<const Spack>& hetfrz_contact_nucleation_tend,
  const uview_2d<const Spack>& hetfrz_deposition_nucleation_tend,
  const view_dnu_table& dnu_table_vals,
  const view_ice_table& ice_table_vals,
  const view_collect_table& collect_table_vals,
  const view_2d_table& revap_table_vals,
  const uview_2d<const Spack>& pres,
  const uview_2d<const Spack>& dpres,
  const uview_2d<const Spack>& dz,
  const uview_2d<const Spack>& nc_nuceat_tend,
  const uview_2d<const Spack>& inv_exner,
  const uview_2d<const Spack>& exner,
  const uview_2d<const Spack>& inv_cld_frac_l,
  const uview_2d<const Spack>& inv_cld_frac_i,
  const uview_2d<const Spack>& inv_cld_frac_r,
  const uview_2d<const Spack>& ni_activated,
  const uview_2d<const Spack>& inv_qc_relvar,
  const uview_2d<const Spack>& cld_frac_i,
  const uview_2d<const Spack>& cld_frac_l,
  const uview_2d<const Spack>& cld_frac_r,
  const uview_2d<const Spack>& qv_prev,
  const uview_2d<const Spack>& t_prev,
  const uview_2d<Spack>& T_atm,
  const uview_2d<Spack>& rho,
  const uview_2d<Spack>& inv_rho,
  const uview_2d<Spack>& qv_sat_l,
  const uview_2d<Spack>& qv_sat_i,
  const uview_2d<Spack>& qv_supersat_i,
  const uview_2d<Spack>& rhofacr,
  const uview_2d<Spack>& rhofaci,
  const uview_2d<Spack>& acn,
  const uview_2d<Spack>& qv,
  const uview_2d<Spack>& th_atm,
  const uview_2d<Spack>& qc,
  const uview_2d<Spack>& nc,
  const uview_2d<Spack>& qr,
  const uview_2d<Spack>& nr,
  const uview_2d<Spack>& qi,
  const uview_2d<Spack>& ni,
  const uview_2d<Spack>& qm,
  const uview_2d<Spack>& bm,
  const uview_2d<Spack>& qc_incld,
  const uview_2d<Spack>& qr_incld,
  const uview_2d<Spack>& qi_incld,
  const uview_2d<Spack>& qm_incld,
  const uview_2d<Spack>& nc_incld,
  const uview_2d<Spack>& nr_incld,
  const uview_2d<Spack>& ni_incld,
  const uview_2d<Spack>& bm_incld,
  const uview_2d<Spack>& mu_c,
  const uview_2d<Spack>& nu,
  const uview_2d<Spack>& lamc,
  const uview_2d<Spack>& cdist,
  const uview_2d<Spack>& cdist1,
  const uview_2d<Spack>& cdistr,
  const uview_2d<Spack>& mu_r,
  const uview_2d<Spack>& lamr,
  const uview_2d<Spack>& logn0r,
  const uview_2d<Spack>& qv2qi_depos_tend,
  const uview_2d<Spack>& precip_total_tend,
  const uview_2d<Spack>& nevapr,
  const uview_2d<Spack>& qr_evap_tend,
  const uview_2d<Spack>& vap_liq_exchange,
  const uview_2d<Spack>& vap_ice_exchange,
  const uview_2d<Spack>& liq_ice_exchange,
  const uview_2d<Spack>& pratot,
  const uview_2d<Spack>& prctot,
  const uview_1d<bool>& nucleationPossible,
  const uview_1d<bool>& hydrometeorsPresent,
  const P3Runtime& runtime_options)
{
  using ExeSpace = typename KT::ExeSpace;
  const Int nk_pack = ekat::npack<Spack>(nk);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);


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
      ekat::subview(cld_frac_r, i), ekat::subview(qv_prev, i), ekat::subview(t_prev, i), ekat::subview(T_atm, i), ekat::subview(rho, i),
      ekat::subview(inv_rho, i), ekat::subview(qv_sat_l, i), ekat::subview(qv_sat_i, i), ekat::subview(qv_supersat_i, i), ekat::subview(rhofacr, i),
      ekat::subview(rhofaci, i), ekat::subview(acn, i), ekat::subview(qv, i), ekat::subview(th_atm, i), ekat::subview(qc, i), ekat::subview(nc, i),
      ekat::subview(qr, i), ekat::subview(nr, i), ekat::subview(qi, i), ekat::subview(ni, i), ekat::subview(qm, i), ekat::subview(bm, i),
      ekat::subview(qc_incld, i),
      ekat::subview(qr_incld, i), ekat::subview(qi_incld, i), ekat::subview(qm_incld, i), ekat::subview(nc_incld, i), ekat::subview(nr_incld, i),
      ekat::subview(ni_incld, i), ekat::subview(bm_incld, i), ekat::subview(mu_c, i), ekat::subview(nu, i), ekat::subview(lamc, i), ekat::subview(cdist, i),
      ekat::subview(cdist1, i), ekat::subview(cdistr, i), ekat::subview(mu_r, i), ekat::subview(lamr, i), ekat::subview(logn0r, i),
      ekat::subview(qv2qi_depos_tend, i), ekat::subview(precip_total_tend, i),
      ekat::subview(nevapr, i), ekat::subview(qr_evap_tend, i), ekat::subview(vap_liq_exchange, i), ekat::subview(vap_ice_exchange, i), ekat::subview(liq_ice_exchange, i),
      ekat::subview(pratot, i), ekat::subview(prctot, i), hydrometeorsPresent(i), nk, runtime_options);

    if (!hydrometeorsPresent(i)) return;
  });
}
#ifdef SCREAM_SYSTEM_WORKAROUND_P3_PART2
#pragma clang optimize on
#endif
} // namespace p3
} // namespace scream

