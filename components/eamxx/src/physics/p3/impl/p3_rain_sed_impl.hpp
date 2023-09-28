#ifndef P3_RAIN_SED_IMPL_HPP
#define P3_RAIN_SED_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace p3 {

/*
 * Implementation of p3 rain sedimentation functions. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::compute_rain_fall_velocity(
  const view_2d_table& vn_table_vals, const view_2d_table& vm_table_vals,
  const Spack& qr_incld, const Spack& rhofacr,
  Spack& nr_incld, Spack& mu_r, Spack& lamr, Spack& V_qr, Spack& V_nr,
  const Smask& context)
{
  Table3 table;
  Spack tmp1, tmp2; //ignore
  get_rain_dsd2(qr_incld, nr_incld, mu_r, lamr, tmp1, tmp2, context);

  if (context.any()) {
    lookup(mu_r, lamr, table, context);
    // mass-weighted fall speed:
    V_qr.set(context, apply_table(vm_table_vals, table) * rhofacr);
    // number-weighted fall speed:
    V_nr.set(context, apply_table(vn_table_vals, table) * rhofacr);
  }
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::rain_sedimentation(
  const uview_1d<const Spack>& rho,
  const uview_1d<const Spack>& inv_rho,
  const uview_1d<const Spack>& rhofacr,
  const uview_1d<const Spack>& cld_frac_r,
  const uview_1d<const Spack>& inv_dz,
  const uview_1d<Spack>& qr_incld,
  const MemberType& team,
  const Workspace& workspace,
  const view_2d_table& vn_table_vals, const view_2d_table& vm_table_vals,
  const Int& nk, const Int& ktop, const Int& kbot, const Int& kdir, const Scalar& dt, const Scalar& inv_dt,
  const uview_1d<Spack>& qr,
  const uview_1d<Spack>& nr,
  const uview_1d<Spack>& nr_incld,
  const uview_1d<Spack>& mu_r,
  const uview_1d<Spack>& lamr,
  const uview_1d<Spack>& precip_liq_flux,
  const uview_1d<Spack>& qr_tend,
  const uview_1d<Spack>& nr_tend,
  Scalar& precip_liq_surf)
{
  // Get temporary workspaces needed for the ice-sed calculation
  uview_1d<Spack> V_qr, V_nr, flux_qx, flux_nx;
  workspace.template take_many_contiguous_unsafe<4>(
    {"V_qr", "V_nr", "flux_qx", "flux_nx"},
    {&V_qr, &V_nr, &flux_qx, &flux_nx});

  const view_1d_ptr_array<Spack, 2>
    fluxes_ptr = {&flux_qx, &flux_nx},
    vs_ptr     = {&V_qr, &V_nr},
    qnr_ptr    = {&qr, &nr};

  const auto sflux_qx = scalarize(flux_qx);

  // find top, determine qxpresent
  const auto sqr = scalarize(qr);
  constexpr Scalar qsmall = C::QSMALL;
  bool log_qxpresent;
  const Int k_qxtop = find_top(team, sqr, qsmall, kbot, ktop, kdir, log_qxpresent);

  if (log_qxpresent) {
    Scalar dt_left   = dt;  // time remaining for sedi over full model (mp) time step
    Scalar prt_accum = 0.0; // precip rate for individual category

    // find bottom
    Int k_qxbot = find_bottom(team, sqr, qsmall, kbot, k_qxtop, kdir, log_qxpresent);

    while (dt_left > C::dt_left_tol) {
      Scalar Co_max = 0.0;
      Int kmin, kmax;
      Int kmin_scalar = ( kdir == 1 ? k_qxbot : k_qxtop);
      Int kmax_scalar = ( kdir == 1 ? k_qxtop : k_qxbot);

      Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, V_qr.extent(0)), [&] (Int k) {
        V_qr(k) = 0;
        V_nr(k) = 0;
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
        const auto qr_gt_small = range_mask && qr_incld(pk) > qsmall;
        if (qr_gt_small.any()) {
          compute_rain_fall_velocity(vn_table_vals, vm_table_vals,
                                     qr_incld(pk), rhofacr(pk),
                                     nr_incld(pk), mu_r(pk), lamr(pk),
				     V_qr(pk), V_nr(pk), qr_gt_small);

	  //in compute_rain_fall_velocity, get_rain_dsd2 keeps the drop-size
	  //distribution within reasonable bounds by modifying nr_incld.
	  //The next line maintains consistency between nr_incld and nr
	  nr(pk).set(qr_gt_small, nr_incld(pk)*cld_frac_r(pk));

        }
        const auto Co_max_local = max(qr_gt_small, 0,
                                      V_qr(pk) * dt_left * inv_dz(pk));
        if (Co_max_local > lmax) lmax = Co_max_local;
      }, Kokkos::Max<Scalar>(Co_max));
      team.team_barrier();

      generalized_sedimentation<2>(rho, inv_rho, inv_dz, team, nk, k_qxtop, k_qxbot, kbot, kdir, Co_max, dt_left, prt_accum, fluxes_ptr, vs_ptr, qnr_ptr);

      //Update _incld values with end-of-step cell-ave values
      //No prob w/ div by cld_frac_r because set to min of 1e-4 in interface.
      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, qr.extent(0)), [&] (int pk) {
	  qr_incld(pk)=qr(pk)/cld_frac_r(pk);
	  nr_incld(pk)=nr(pk)/cld_frac_r(pk);
	});

      // AaronDonahue, precip_liq_flux output
      kmin_scalar = ( kdir == 1 ? k_qxbot+1 : k_qxtop+1);
      kmax_scalar = ( kdir == 1 ? k_qxtop+1 : k_qxbot+1);
      ekat::impl::set_min_max(kmin_scalar, kmax_scalar, kmin, kmax, Spack::n);
      Kokkos::parallel_for(
       Kokkos::TeamVectorRange(team, kmax-kmin+1), [&] (int pk_) {
        const int pk = kmin + pk_;
        const auto range_pack = ekat::range<IntSmallPack>(pk*Spack::n);
        const auto range_mask = range_pack >= kmin_scalar && range_pack <= kmax_scalar;
        auto index_pack = range_pack-1;
        const auto lt_zero = index_pack < 0;
        index_pack.set(lt_zero, 0);
        const auto flux_qx_pk = index(sflux_qx, index_pack);
        precip_liq_flux(pk).set(range_mask, precip_liq_flux(pk) + flux_qx_pk);
      });
    }
    Kokkos::single(
      Kokkos::PerTeam(team), [&] () {
        precip_liq_surf += prt_accum * C::INV_RHO_H2O * inv_dt;
      });
  }

  const Int nk_pack = ekat::npack<Spack>(nk);
  Kokkos::parallel_for(
   Kokkos::TeamVectorRange(team, nk_pack), [&] (int pk) {
    qr_tend(pk) = (qr(pk) - qr_tend(pk)) * inv_dt; // Rain sedimentation tendency, measure
    nr_tend(pk) = (nr(pk) - nr_tend(pk)) * inv_dt; // Rain # sedimentation tendency, measure
  });

  workspace.template release_many_contiguous<4>(
    {&V_qr, &V_nr, &flux_qx, &flux_nx});
}

} // namespace p3
} // namespace scream

#endif
