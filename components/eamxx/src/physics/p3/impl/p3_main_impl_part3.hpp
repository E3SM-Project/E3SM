#ifndef P3_MAIN_IMPL_PART_3_HPP
#define P3_MAIN_IMPL_PART_3_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs
#include "physics/share/physics_saturation_impl.hpp"

#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 main function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main_part3(
  const MemberType& team,
  const Int& nk_pack,
  const Scalar& max_total_ni,
  const view_dnu_table& dnu,
  const view_ice_table& ice_table_vals,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<const Spack>& cld_frac_l,
  const uview_1d<const Spack>& cld_frac_r,
  const uview_1d<const Spack>& cld_frac_i,
  const uview_1d<Spack>& rho,
  const uview_1d<Spack>& inv_rho,
  const uview_1d<Spack>& rhofaci,
  const uview_1d<Spack>& qv,
  const uview_1d<Spack>& th_atm,
  const uview_1d<Spack>& qc,
  const uview_1d<Spack>& nc,
  const uview_1d<Spack>& qr,
  const uview_1d<Spack>& nr,
  const uview_1d<Spack>& qi,
  const uview_1d<Spack>& ni,
  const uview_1d<Spack>& qm,
  const uview_1d<Spack>& bm,
  const uview_1d<Spack>& mu_c,
  const uview_1d<Spack>& nu,
  const uview_1d<Spack>& lamc,
  const uview_1d<Spack>& mu_r,
  const uview_1d<Spack>& lamr,
  const uview_1d<Spack>& vap_liq_exchange,
  const uview_1d<Spack>& ze_rain,
  const uview_1d<Spack>& ze_ice,
  const uview_1d<Spack>& diag_vm_qi,
  const uview_1d<Spack>& diag_eff_radius_qi,
  const uview_1d<Spack>& diag_diam_qi,
  const uview_1d<Spack>& rho_qi,
  const uview_1d<Spack>& diag_equiv_reflectivity,
  const uview_1d<Spack>& diag_eff_radius_qc,
  const uview_1d<Spack>& diag_eff_radius_qr,
  const P3Runtime& runtime_options)
{
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar inv_cp       = C::INV_CP;
  constexpr Scalar nsmall       = C::NSMALL;
  constexpr Scalar latvap       = C::LatVap;
  constexpr Scalar latice       = C::LatIce;

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, nk_pack), [&] (Int k) {

    Spack
      ignore1  (0),
      ignore2  (0),
      table_val_qi_fallspd   (0),
      table_val_ice_eff_radius   (0),
      table_val_ni_lammax   (0),
      table_val_ni_lammin   (0),
      table_val_ice_reflectivity   (0),
      table_val_ice_mean_diam   (0),
      table_val_ice_bulk_dens   (0);

    // Cloud
    {
      const auto qc_gt_small = qc(k) >= qsmall;
      const auto qc_small    = !qc_gt_small;
      const auto qc_incld = qc(k)/cld_frac_l(k);
      auto nc_incld = nc(k)/cld_frac_l(k);

      get_cloud_dsd2(qc_incld, nc_incld, mu_c(k), rho(k), nu(k), dnu, lamc(k), ignore1, ignore2, qc_gt_small);

      if (qc_gt_small.any()) {
        nc(k).set(qc_gt_small,nc_incld*cld_frac_l(k)); //cld_dsd2 might have changed incld nc... need consistency.
        diag_eff_radius_qc(k).set(qc_gt_small, sp(0.5) * (mu_c(k) + 3) / lamc(k));
      }
      if (qc_small.any()) {
        qv(k)                .set(qc_small, qv(k)+qc(k));
        th_atm(k)            .set(qc_small, th_atm(k)-inv_exner(k)*qc(k)*latvap*inv_cp);
        vap_liq_exchange(k)  .set(qc_small, vap_liq_exchange(k) - qc(k));
        qc(k)                .set(qc_small, 0);
        nc(k)                .set(qc_small, 0);
      }
    }

    // Rain
    {
      const auto qr_gt_small = qr(k) >= qsmall;
      const auto qr_small    = !qr_gt_small;
      const auto qr_incld = qr(k)/cld_frac_r(k);
      auto nr_incld = nr(k)/cld_frac_r(k); //nr_incld is updated in get_rain_dsd2 but isn't used again

      get_rain_dsd2(
        qr_incld, nr_incld, mu_r(k), lamr(k), runtime_options, qr_gt_small);

      //Note that integrating over the drop-size PDF as done here should only be done to in-cloud
      //quantities but radar reflectivity is likely meant to be a cell ave. Thus nr in the next line
      //really should be cld_frac_r * nr/cld_frac_r. Not doing that since cld_frac_r cancels out.
      if (qr_gt_small.any()) {
        nr(k).set(qr_gt_small,nr_incld*cld_frac_r(k)); //rain_dsd2 might have changed incld nr... need consistency.
        ze_rain(k).set(qr_gt_small, nr(k)*(mu_r(k)+6)*(mu_r(k)+5)*(mu_r(k)+4)*
                       (mu_r(k)+3)*(mu_r(k)+2)*(mu_r(k)+1)/pow(lamr(k), sp(6.0))); // once f90 is gone, 6 can be int
        ze_rain(k).set(qr_gt_small, max(ze_rain(k), sp(1.e-22)));
        diag_eff_radius_qr(k).set(qr_gt_small, sp(1.5) / lamr(k));
      }

      if (qr_small.any()) {
        qv(k)              .set(qr_small, qv(k) + qr(k));
        th_atm(k)          .set(qr_small, th_atm(k) - inv_exner(k)*qr(k)*latvap*inv_cp);
        vap_liq_exchange(k).set(qr_small, vap_liq_exchange(k) - qr(k));
        qr(k)              .set(qr_small, 0);
        nr(k)              .set(qr_small, 0);
      }
    }

    // Ice
    {
      const auto qi_gt_small = qi(k) >= qsmall;
      const auto qi_small    = !qi_gt_small;

      // impose lower limits to prevent taking log of # < 0
      ni(k) = max(ni(k), nsmall);

      auto qi_incld = qi(k)/cld_frac_i(k);
      auto ni_incld = ni(k)/cld_frac_i(k);
      auto qm_incld = qm(k)/cld_frac_i(k);
      auto bm_incld = bm(k)/cld_frac_i(k);

      const auto rhop = calc_bulk_rho_rime(qi_incld, qm_incld, bm_incld, runtime_options, qi_gt_small);
      qm(k).set(qi_gt_small, qm_incld*cld_frac_i(k) );
      bm(k).set(qi_gt_small, bm_incld*cld_frac_i(k) );

      impose_max_total_ni(ni_incld, max_total_ni, inv_rho(k));

      TableIce table_ice;
      lookup_ice(qi_incld, ni_incld, qm_incld, rhop, table_ice, qi_gt_small);

      table_val_qi_fallspd.set(qi_gt_small, apply_table_ice(1,  ice_table_vals, table_ice, qi_gt_small));
      table_val_ice_eff_radius.set(qi_gt_small, apply_table_ice(5,  ice_table_vals, table_ice, qi_gt_small));
      table_val_ni_lammax.set(qi_gt_small, apply_table_ice(6,  ice_table_vals, table_ice, qi_gt_small));
      table_val_ni_lammin.set(qi_gt_small, apply_table_ice(7,  ice_table_vals, table_ice, qi_gt_small));
      table_val_ice_reflectivity.set(qi_gt_small, apply_table_ice(8,  ice_table_vals, table_ice, qi_gt_small));
      table_val_ice_mean_diam.set(qi_gt_small, apply_table_ice(10, ice_table_vals, table_ice, qi_gt_small));
      table_val_ice_bulk_dens.set(qi_gt_small, apply_table_ice(11, ice_table_vals, table_ice, qi_gt_small));

      // impose mean ice size bounds (i.e. apply lambda limiters)
      // note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
      ni_incld.set(qi_gt_small, min(ni_incld, table_val_ni_lammax * ni_incld));
      ni_incld.set(qi_gt_small, max(ni_incld, table_val_ni_lammin * ni_incld));
      ni(k) = ni_incld*cld_frac_i(k);

      // --this should already be done in s/r 'calc_bulkRhoRime'
      const auto qm_small = qm(k) < qsmall && qi_gt_small;
      qm(k).set(qm_small, 0);
      bm(k).set(qm_small, 0);

      // note that reflectivity from lookup table is normalized, so we need to multiply by N
      diag_vm_qi(k) .set(qi_gt_small, table_val_qi_fallspd * rhofaci(k));
      diag_eff_radius_qi(k).set(qi_gt_small, table_val_ice_eff_radius); // units are in m
      diag_diam_qi(k)  .set(qi_gt_small, table_val_ice_mean_diam);
      rho_qi(k).set(qi_gt_small, table_val_ice_bulk_dens);

      // note factor of air density below is to convert from m^6/kg to m^6/m^3
      ze_ice(k).set(qi_gt_small, ze_ice(k) + sp(0.1892)*table_val_ice_reflectivity*ni_incld*rho(k));   // sum contribution from each ice category (note: 0.1892 = 0.176/0.93);
      ze_ice(k).set(qi_gt_small, max(ze_ice(k), sp(1.e-22)));

      //above formula for ze only makes sense for in-cloud vals, but users expect cell-ave output.
      ze_ice(k).set(qi_gt_small, ze_ice(k)*cld_frac_i(k));

      qv(k).set(qi_small, qv(k) + qi(k));
      th_atm(k).set(qi_small, th_atm(k) - inv_exner(k)*qi(k)*(latvap+latice)*inv_cp);
      qi(k).set(qi_small, 0);
      ni(k).set(qi_small, 0);
      qm(k).set(qi_small, 0);
      bm(k).set(qi_small, 0);
      diag_diam_qi(k).set(qi_small, 0);
    }

    // sum ze components and convert to dBZ
    diag_equiv_reflectivity(k) = 10 * log10((ze_rain(k) + ze_ice(k))*sp(1.e18));

    // if qr is very small then set Nr to 0 (needs to be done here after call
    // to ice lookup table because a minimum Nr of nsmall will be set otherwise even if qr=0)
    nr(k).set(qr(k) < qsmall, 0);
  });
  team.team_barrier();
}

} // namespace p3
} // namespace scream

#endif // P3_MAIN_IMPL_PART_3_HPP
