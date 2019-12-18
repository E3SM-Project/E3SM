#ifndef P3_FUNCTIONS_ICE_SED_IMPL_HPP
#define P3_FUNCTIONS_ICE_SED_IMPL_HPP

#include "p3_function.hpp" // for ETI only but harmless for GPU

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
  const Smask& qi_gt_small, const Spack& qi_tot, Spack& qi_rim, Spack& bi_rim)
{
  constexpr Scalar bsmall       = C::BSMALL;
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar rho_rime_min = C::rho_rimeMin;
  constexpr Scalar rho_rime_max = C::rho_rimeMax;

  Spack rho_rime(0);

  Smask bi_rim_gt_small = (bi_rim >= bsmall) && qi_gt_small;
  if (bi_rim_gt_small.any()) {
    rho_rime.set(qi_gt_small && bi_rim_gt_small, qi_rim / bi_rim);
  }

  Smask rho_rime_lt_min = rho_rime < rho_rime_min;
  Smask rho_rime_gt_max = rho_rime > rho_rime_max;

  // impose limits on rho_rime;  adjust bi_rim if needed
  rho_rime.set(bi_rim_gt_small && rho_rime_lt_min, rho_rime_min);
  rho_rime.set(bi_rim_gt_small && rho_rime_gt_max, rho_rime_max);
  Smask adjust = bi_rim_gt_small && (rho_rime_gt_max || rho_rime_lt_min);
  if (adjust.any()) {
    bi_rim.set(adjust, qi_rim / rho_rime);
  }

  qi_rim.set  (qi_gt_small && !bi_rim_gt_small, 0);
  bi_rim.set  (qi_gt_small && !bi_rim_gt_small, 0);
  rho_rime.set(qi_gt_small && !bi_rim_gt_small, 0);

  // set upper constraint qi_rim <= qi_tot
  Smask qi_rim_gt_qi = (qi_rim > qi_tot) && (rho_rime > 0) && qi_gt_small;
  if (qi_rim_gt_qi.any()) {
    qi_rim.set(qi_rim_gt_qi, qi_tot);
    bi_rim.set(qi_rim_gt_qi, qi_rim / rho_rime);
  }

  // impose consistency
  Smask qi_rim_lt_small = (qi_rim < qsmall) && qi_gt_small;
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
  const uview_1d<const Spack>& icldm,
  const uview_1d<const Spack>& inv_dzq,
  const MemberType& team,
  const Workspace& workspace,
  const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& odt,
  const uview_1d<Spack>& qitot,
  const uview_1d<Spack>& qitot_incld,
  const uview_1d<Spack>& nitot,
  const uview_1d<Spack>& nitot_incld,
  const uview_1d<Spack>& qirim,
  const uview_1d<Spack>& qirim_incld,
  const uview_1d<Spack>& birim,
  const uview_1d<Spack>& birim_incld,
  const uview_1d<Spack>& qi_tend,
  const uview_1d<Spack>& ni_tend,
  const view_itab_table& itab,
  Scalar& prt_sol)
{
  // Get temporary workspaces needed for the ice-sed calculation
  uview_1d<Spack> V_qit, V_nit, flux_nit, flux_bir, flux_qir, flux_qit;
  workspace.template take_many_contiguous_unsafe<6>(
    {"V_qit", "V_nit", "flux_nit", "flux_bir", "flux_qir", "flux_qit"},
    {&V_qit, &V_nit, &flux_nit, &flux_bir, &flux_qir, &flux_qit});

  const view_1d_ptr_array<Spack, 4>
    fluxes_ptr = {&flux_qit, &flux_nit, &flux_qir, &flux_bir},
    vs_ptr     = {&V_qit, &V_nit, &V_qit, &V_qit},
    qnr_ptr    = {&qitot, &nitot, &qirim, &birim};

  // find top, determine qxpresent
  const auto sqitot = scalarize(qitot);
  constexpr Scalar qsmall = C::QSMALL;
  constexpr Scalar nsmall = C::NSMALL;
  bool log_qxpresent;
  const Int k_qxtop = find_top(team, sqitot, qsmall, kbot, ktop, kdir, log_qxpresent);

  if (log_qxpresent) {
    Scalar dt_left   = dt;  // time remaining for sedi over full model (mp) time step
    Scalar prt_accum = 0.0; // precip rate for individual category

    // find bottom
    Int k_qxbot = find_bottom(team, sqitot, qsmall, kbot, k_qxtop, kdir, log_qxpresent);

    while (dt_left > C::dt_left_tol) {
      Scalar Co_max = 0.0;
      Int kmin, kmax;
      const Int kmin_scalar = ( kdir == 1 ? k_qxbot : k_qxtop);
      const Int kmax_scalar = ( kdir == 1 ? k_qxtop : k_qxbot);

      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, V_qit.extent(0)), [&] (Int k) {
          V_qit(k) = 0;
          V_nit(k) = 0;
      });
      team.team_barrier();

      // Convert top/bot to pack indices
      util::set_min_max(k_qxbot, k_qxtop, kmin, kmax, Spack::n);

      // compute Vq, Vn (get values from lookup table)
      Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, kmax-kmin+1), [&] (int pk_, Scalar& lmax) {

        const int pk = kmin + pk_;
        const auto range_pack = scream::pack::range<IntSmallPack>(pk*Spack::n);
        const auto range_mask = range_pack >= kmin_scalar && range_pack <= kmax_scalar;
        const auto qi_gt_small = range_mask && qitot_incld(pk) > qsmall;
        if (qi_gt_small.any()) {
          // impose lower limits to prevent log(<0)
          nitot_incld(pk).set(qi_gt_small, pack::max(nitot_incld(pk), nsmall));

          const auto rhop = calc_bulk_rho_rime(qi_gt_small, qitot_incld(pk), qirim_incld(pk), birim_incld(pk));

          TableIce t;
          lookup_ice(qi_gt_small, qitot_incld(pk), nitot_incld(pk), qirim_incld(pk), rhop, t);

          const auto f1pr01 = apply_table_ice(qi_gt_small, 0, itab, t);
          const auto f1pr02 = apply_table_ice(qi_gt_small, 1, itab, t);
          const auto f1pr09 = apply_table_ice(qi_gt_small, 6, itab, t);
          const auto f1pr10 = apply_table_ice(qi_gt_small, 7, itab, t);

          // impose mean ice size bounds (i.e. apply lambda limiters)
          // note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
          nitot_incld(pk).set(qi_gt_small, pack::min(nitot_incld(pk), f1pr09 * nitot_incld(pk)));
          nitot_incld(pk).set(qi_gt_small, pack::max(nitot_incld(pk), f1pr10 * nitot_incld(pk)));
          nitot(pk).set(qi_gt_small, nitot_incld(pk) * icldm(pk));

          V_qit(pk).set(qi_gt_small, f1pr02 * rhofaci(pk)); // mass-weighted   fall speed (with density factor)
          V_nit(pk).set(qi_gt_small, f1pr01 * rhofaci(pk)); // number-weighted fall speed (with density factor)
        }
        const auto Co_max_local = max(qi_gt_small, -1,
                                      V_qit(pk) * dt_left * inv_dzq(pk));
        if (Co_max_local > lmax) lmax = Co_max_local;
      }, Kokkos::Max<Scalar>(Co_max));
      team.team_barrier();

      generalized_sedimentation<4>(rho, inv_rho, inv_dzq, team, nk, k_qxtop, k_qxbot, kbot, kdir, Co_max, dt_left, prt_accum, fluxes_ptr, vs_ptr, qnr_ptr);
    }
    Kokkos::single(
      Kokkos::PerTeam(team), [&] () {
        prt_sol += prt_accum * C::INV_RHOW * odt;
    });
  }

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, qi_tend.extent(0)), [&] (int pk) {
      qi_tend(pk) = (qitot(pk) - qi_tend(pk)) * odt; // Liq. sedimentation tendency, measure
      ni_tend(pk) = (nitot(pk) - ni_tend(pk)) * odt; // Liq. # sedimentation tendency, measure
  });

  workspace.template release_many_contiguous<6>(
    {&V_qit, &V_nit, &flux_nit, &flux_bir, &flux_qir, &flux_qit});
}

} // namespace p3
} // namespace scream

#endif
