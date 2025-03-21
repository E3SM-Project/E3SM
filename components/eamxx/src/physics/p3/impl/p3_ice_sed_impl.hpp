#ifndef P3_ICE_SED_IMPL_HPP
#define P3_ICE_SED_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 ice sedimentation function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
typename Functions<S,D>::Spack
Functions<S,D>
::calc_bulk_rho_rime(
  const Spack& qi_tot, Spack& qi_rim, Spack& bi_rim,
  const P3Runtime& runtime_options,
  const Smask& context)
{
  constexpr Scalar bsmall       = C::BSMALL;
  constexpr Scalar qsmall       = C::QSMALL;
  const Scalar min_rime_rho     = runtime_options.min_rime_rho;
  const Scalar max_rime_rho     = runtime_options.max_rime_rho;

  Spack rho_rime(0);

  Smask bi_rim_gt_small = (bi_rim >= bsmall) && context;
  Smask bi_rim_lt_small = (bi_rim <  bsmall) && context;
  if (bi_rim_gt_small.any()) {
    rho_rime.set(bi_rim_gt_small, qi_rim / bi_rim);
  }

  Smask rho_rime_lt_min = rho_rime < min_rime_rho;
  Smask rho_rime_gt_max = rho_rime > max_rime_rho;

  // impose limits on rho_rime;  adjust bi_rim if needed
  rho_rime.set(bi_rim_gt_small && rho_rime_lt_min, min_rime_rho);
  rho_rime.set(bi_rim_gt_small && rho_rime_gt_max, max_rime_rho);
  Smask adjust = bi_rim_gt_small && (rho_rime_gt_max || rho_rime_lt_min);
  if (adjust.any()) {
    bi_rim.set(adjust, qi_rim / rho_rime);
  }

  qi_rim.set  (bi_rim_lt_small, 0);
  bi_rim.set  (bi_rim_lt_small, 0);
  rho_rime.set(bi_rim_lt_small, 0);

  // set upper constraint qi_rim <= qi_tot
  Smask qi_rim_gt_qi = (qi_rim > qi_tot) && (rho_rime > 0) && context;
  if (qi_rim_gt_qi.any()) {
    qi_rim.set(qi_rim_gt_qi, qi_tot);
    bi_rim.set(qi_rim_gt_qi, qi_rim / rho_rime);
  }

  // impose consistency
  Smask qi_rim_lt_small = (qi_rim < qsmall) && context;
  qi_rim.set(qi_rim_lt_small, 0);
  bi_rim.set(qi_rim_lt_small, 0);

  return rho_rime;
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::ice_sedimentation(
  const uview_1d<const Spack>& rho,
  const uview_1d<const Spack>& inv_rho,
  const uview_1d<const Spack>& rhofaci,
  const uview_1d<const Spack>& cld_frac_i,
  const uview_1d<const Spack>& inv_dz,
  const MemberType& team,
  const Workspace& workspace,
  const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& inv_dt,
  const uview_1d<Spack>& qi,
  const uview_1d<Spack>& qi_incld,
  const uview_1d<Spack>& ni,
  const uview_1d<Spack>& ni_incld,
  const uview_1d<Spack>& qm,
  const uview_1d<Spack>& qmr,
  const uview_1d<Spack>& qm_incld,
  const uview_1d<Spack>& bm,
  const uview_1d<Spack>& bm_incld,
  const uview_1d<Spack>& qi_tend,
  const uview_1d<Spack>& ni_tend,
  const view_ice_table& ice_table_vals,
  Scalar& precip_ice_surf,
  const P3Runtime& runtime_options)
{
  // Get temporary workspaces needed for the ice-sed calculation
  uview_1d<Spack> V_qit, V_nit, flux_nit, flux_bir, flux_qir, flux_qit, flux_qirr;
  workspace.template take_many_contiguous_unsafe<7>(
    {"V_qit", "V_nit", "flux_nit", "flux_bir", "flux_qir", "flux_qit","flux_qirr"},
    {&V_qit, &V_nit, &flux_nit, &flux_bir, &flux_qir, &flux_qit, &flux_qirr});

  const view_1d_ptr_array<Spack, 5>
    fluxes_ptr = {&flux_qit, &flux_nit, &flux_qir, &flux_bir,&flux_qirr},
    vs_ptr     = {&V_qit, &V_nit, &V_qit, &V_qit, &V_qit},
    qnr_ptr    = {&qi, &ni, &qm, &bm, &qmr};

  // find top, determine qxpresent
  const auto sqi = scalarize(qi);
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar nsmall = C::NSMALL;

  const Scalar ice_sedimentation_factor = runtime_options.ice_sedimentation_factor;

  bool log_qxpresent;
  const Int k_qxtop = find_top(team, sqi, qsmall, kbot, ktop, kdir, log_qxpresent);

  if (log_qxpresent) {
    Scalar dt_left   = dt;  // time remaining for sedi over full model (mp) time step
    Scalar prt_accum = 0.0; // precip rate for individual category

    // find bottom
    Int k_qxbot = find_bottom(team, sqi, qsmall, kbot, k_qxtop, kdir, log_qxpresent);

    while (dt_left > C::dt_left_tol) {
      Scalar Co_max = 0.0;
      Int kmin, kmax;
      const Int kmin_scalar = ( kdir == 1 ? k_qxbot : k_qxtop);
      const Int kmax_scalar = ( kdir == 1 ? k_qxtop : k_qxbot);

      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, V_qit.extent(0)), [&] (Int k) {
          V_qit(k) = 0;
          V_nit(k) = 0;
      });
      team.team_barrier();

      // Convert top/bot to pack indices
      ekat::impl::set_min_max(k_qxbot, k_qxtop, kmin, kmax, Spack::n);

      // compute Vq, Vn (get values from lookup table)
      Kokkos::parallel_reduce(
        Kokkos::TeamVectorRange(team, kmax-kmin+1), [&] (int pk_, Scalar& lmax) {

        const int pk = kmin + pk_;
        const auto range_pack = ekat::range<IntSmallPack>(pk*Spack::n);
        const auto range_mask = range_pack >= kmin_scalar && range_pack <= kmax_scalar;
        const auto qi_gt_small = range_mask && qi_incld(pk) > qsmall;
        if (qi_gt_small.any()) {
          // impose lower limits to prevent log(<0)
          ni_incld(pk).set(qi_gt_small, max(ni_incld(pk), nsmall));

          const auto rhop = calc_bulk_rho_rime(qi_incld(pk), qm_incld(pk), bm_incld(pk), runtime_options, qi_gt_small);
          qm(pk).set(qi_gt_small, qm_incld(pk)*cld_frac_i(pk) );
          bm(pk).set(qi_gt_small, bm_incld(pk)*cld_frac_i(pk) );

          TableIce tab;
          lookup_ice(qi_incld(pk), ni_incld(pk), qm_incld(pk), rhop, tab, qi_gt_small);

          const auto table_val_ni_fallspd = apply_table_ice(0, ice_table_vals, tab, qi_gt_small);
          const auto table_val_qi_fallspd = apply_table_ice(1, ice_table_vals, tab, qi_gt_small);
          const auto table_val_ni_lammax = apply_table_ice(6, ice_table_vals, tab, qi_gt_small);
          const auto table_val_ni_lammin = apply_table_ice(7, ice_table_vals, tab, qi_gt_small);

          // impose mean ice size bounds (i.e. apply lambda limiters)
          // note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
          ni_incld(pk).set(qi_gt_small, min(ni_incld(pk), table_val_ni_lammax * ni_incld(pk)));
          ni_incld(pk).set(qi_gt_small, max(ni_incld(pk), table_val_ni_lammin * ni_incld(pk)));
          ni(pk).set(qi_gt_small, ni_incld(pk) * cld_frac_i(pk));

          V_qit(pk).set(qi_gt_small, ice_sedimentation_factor * table_val_qi_fallspd * rhofaci(pk)); // mass-weighted   fall speed (with density factor)
          V_nit(pk).set(qi_gt_small, ice_sedimentation_factor * table_val_ni_fallspd * rhofaci(pk)); // number-weighted fall speed (with density factor)
        }
        const auto Co_max_local = max(qi_gt_small, 0,
                                      V_qit(pk) * dt_left * inv_dz(pk));
        if (Co_max_local > lmax) lmax = Co_max_local;
      }, Kokkos::Max<Scalar>(Co_max));
      team.team_barrier();

      generalized_sedimentation<5>(rho, inv_rho, inv_dz, team, nk, k_qxtop, k_qxbot, kbot, kdir, Co_max, dt_left, prt_accum, fluxes_ptr, vs_ptr, qnr_ptr);

      //Update _incld values with end-of-step cell-ave values
      //No prob w/ div by cld_frac_i because set to min of 1e-4 in interface.
      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, qi.extent(0)), [&] (int pk) {
	  qi_incld(pk)=qi(pk)/cld_frac_i(pk);
	  ni_incld(pk)=ni(pk)/cld_frac_i(pk);
	  qm_incld(pk)=qm(pk)/cld_frac_i(pk);
	  bm_incld(pk)=bm(pk)/cld_frac_i(pk);
	});


    } //end CFL substep loop

    Kokkos::single(
      Kokkos::PerTeam(team), [&] () {
        precip_ice_surf += prt_accum * C::INV_RHO_H2O * inv_dt;
    });
  }

  const Int nk_pack = ekat::npack<Spack>(nk);
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, nk_pack), [&] (int pk) {
      qi_tend(pk) = (qi(pk) - qi_tend(pk)) * inv_dt; // Liq. sedimentation tendency, measure
      ni_tend(pk) = (ni(pk) - ni_tend(pk)) * inv_dt; // Liq. # sedimentation tendency, measure
  });

  workspace.template release_many_contiguous<7>(
    {&V_qit, &V_nit, &flux_nit, &flux_bir, &flux_qir, &flux_qit, &flux_qirr});
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::homogeneous_freezing(
  const uview_1d<const Spack>& T_atm,
  const uview_1d<const Spack>& inv_exner,
  const MemberType& team,
  const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir,
  const uview_1d<Spack>& qc,
  const uview_1d<Spack>& nc,
  const uview_1d<Spack>& qr,
  const uview_1d<Spack>& nr,
  const uview_1d<Spack>& qi,
  const uview_1d<Spack>& ni,
  const uview_1d<Spack>& qm,
  const uview_1d<Spack>& qmr,
  const uview_1d<Spack>& bm,
  const uview_1d<Spack>& th_atm)
{
  constexpr Scalar qsmall          = C::QSMALL;
  constexpr Scalar nsmall          = C::NSMALL;
  constexpr Scalar T_homogfrz       = C::T_homogfrz;
  constexpr Scalar inv_rho_rimeMax = C::INV_RHO_RIMEMAX;
  constexpr Scalar inv_cp          = C::INV_CP;
  constexpr Scalar latice          = C::LatIce;

  const Int kmin_scalar = ( kdir == 1 ? kbot : ktop);
  const Int kmax_scalar = ( kdir == 1 ? ktop : kbot);

  // Convert top/bot to pack indices
  Int kmin, kmax;
  ekat::impl::set_min_max(kbot, ktop, kmin, kmax, Spack::n);

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, kmax-kmin+1), [&] (int pk_) {

    const int pk = kmin + pk_;

    // Set up masks
    const auto range_pack    = ekat::range<IntSmallPack>(pk*Spack::n);
    const auto range_mask    = range_pack >= kmin_scalar && range_pack <= kmax_scalar;
    const auto t_lt_homogf   = T_atm(pk) < T_homogfrz;
    const auto qc_ge_small   = range_mask && t_lt_homogf && qc(pk) >= qsmall;
    const auto qr_ge_small   = range_mask && t_lt_homogf && qr(pk) >= qsmall;

    Spack Qc_nuc(qc(pk)), Qr_nuc(qr(pk)), Nc_nuc(max(nc(pk), nsmall)), Nr_nuc(max(nr(pk), nsmall));

    qm(pk).set(qc_ge_small, qm(pk) + Qc_nuc);
    qi(pk).set(qc_ge_small, qi(pk) + Qc_nuc);
    bm(pk).set(qc_ge_small, bm(pk) + Qc_nuc*inv_rho_rimeMax);
    ni(pk).set(qc_ge_small, ni(pk) + Nc_nuc);
    th_atm(pk).set   (qc_ge_small, th_atm(pk) + inv_exner(pk)*Qc_nuc*latice*inv_cp);

    qm(pk).set(qr_ge_small, qm(pk) + Qr_nuc);
    qmr(pk).set(qr_ge_small, qmr(pk) + Qr_nuc);
    qi(pk).set(qr_ge_small, qi(pk) + Qr_nuc);
    bm(pk).set(qr_ge_small, bm(pk) + Qr_nuc*inv_rho_rimeMax);
    ni(pk).set(qr_ge_small, ni(pk) + Nr_nuc);
    th_atm(pk).set   (qr_ge_small, th_atm(pk) + inv_exner(pk)*Qr_nuc*latice*inv_cp);

    qc(pk).set(qc_ge_small, 0);
    nc(pk).set(qc_ge_small, 0);
    qr(pk).set(qr_ge_small, 0);
    nr(pk).set(qr_ge_small, 0);
  });
}

} // namespace p3
} // namespace scream

#endif // P3_ICE_SED_IMPL_HPP
