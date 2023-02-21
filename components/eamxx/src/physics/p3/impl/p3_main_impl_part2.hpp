#ifndef P3_MAIN_IMPL_PART_2_HPP
#define P3_MAIN_IMPL_PART_2_HPP

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
::p3_main_part2(
  const MemberType& team,
  const Int& nk_pack,
  const bool& predictNc,
  const bool& do_prescribed_CCN,
  const Scalar& dt,
  const Scalar& inv_dt,
  const view_dnu_table& dnu,
  const view_ice_table& ice_table_vals,
  const view_collect_table& collect_table_vals,
  const view_2d_table& revap_table_vals,
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& dpres,
  const uview_1d<const Spack>& dz,
  const uview_1d<const Spack>& nc_nuceat_tend,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<const Spack>& exner,
  const uview_1d<const Spack>& inv_cld_frac_l,
  const uview_1d<const Spack>& inv_cld_frac_i,
  const uview_1d<const Spack>& inv_cld_frac_r,
  const uview_1d<const Spack>& ni_activated,
  const uview_1d<const Spack>& inv_qc_relvar,
  const uview_1d<const Spack>& cld_frac_i,
  const uview_1d<const Spack>& cld_frac_l,
  const uview_1d<const Spack>& cld_frac_r,
  const uview_1d<const Spack>& qv_prev,
  const uview_1d<const Spack>& t_prev,
  const uview_1d<Spack>& T_atm,
  const uview_1d<Spack>& rho,
  const uview_1d<Spack>& inv_rho,
  const uview_1d<Spack>& qv_sat_l,
  const uview_1d<Spack>& qv_sat_i,
  const uview_1d<Spack>& qv_supersat_i,
  const uview_1d<Spack>& rhofacr,
  const uview_1d<Spack>& rhofaci,
  const uview_1d<Spack>& acn,
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
  const uview_1d<Spack>& latent_heat_vapor,
  const uview_1d<Spack>& latent_heat_sublim,
  const uview_1d<Spack>& latent_heat_fusion,
  const uview_1d<Spack>& qc_incld,
  const uview_1d<Spack>& qr_incld,
  const uview_1d<Spack>& qi_incld,
  const uview_1d<Spack>& qm_incld,
  const uview_1d<Spack>& nc_incld,
  const uview_1d<Spack>& nr_incld,
  const uview_1d<Spack>& ni_incld,
  const uview_1d<Spack>& bm_incld,
  const uview_1d<Spack>& mu_c,
  const uview_1d<Spack>& nu,
  const uview_1d<Spack>& lamc,
  const uview_1d<Spack>& cdist,
  const uview_1d<Spack>& cdist1,
  const uview_1d<Spack>& cdistr,
  const uview_1d<Spack>& mu_r,
  const uview_1d<Spack>& lamr,
  const uview_1d<Spack>& logn0r,
  const uview_1d<Spack>& qv2qi_depos_tend,
  const uview_1d<Spack>& precip_total_tend,
  const uview_1d<Spack>& nevapr,
  const uview_1d<Spack>& qr_evap_tend,
  const uview_1d<Spack>& vap_liq_exchange,
  const uview_1d<Spack>& vap_ice_exchange,
  const uview_1d<Spack>& liq_ice_exchange,
  const uview_1d<Spack>& pratot,
  const uview_1d<Spack>& prctot,
  bool& hydrometeorsPresent, const Int& nk)
{
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar nsmall       = C::NSMALL;
  constexpr Scalar T_zerodegc   = C::T_zerodegc;
  constexpr Scalar max_total_ni = C::max_total_ni;
  constexpr Scalar f1r          = C::f1r;
  constexpr Scalar f2r          = C::f2r;
  constexpr Scalar nmltratio    = C::nmltratio;
  constexpr Scalar inv_cp       = C::INV_CP;

  team.team_barrier();
  hydrometeorsPresent = false;
  team.team_barrier();

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, nk_pack), [&] (Int k) {

    //compute mask to identify padded values in packs, which shouldn't be used in calculations
    const auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);
    const auto range_mask = range_pack < nk;
      
    // if relatively dry and no hydrometeors at this level, skip to end of k-loop (i.e. skip this level)
    const auto skip_all = ( !range_mask ||
        (qc(k)<qsmall && qr(k)<qsmall && qi(k)<qsmall &&
         T_atm(k)<T_zerodegc && qv_supersat_i(k)< -0.05) );
    
    if (skip_all.all()) {
      return; // skip all process rates
    }
    const auto not_skip_all = !skip_all;

    // All microphysics tendencies will be computed as IN-CLOUD, they will be mapped back to cell-average later.

    Spack
      // initialize warm-phase process rates
      qc2qr_accret_tend   (0), // cloud droplet accretion by rain
      qr2qv_evap_tend   (0), // rain evaporation
      qc2qr_autoconv_tend   (0), // cloud droplet autoconversion to rain
      nc_accret_tend   (0), // change in cloud droplet number from accretion by rain
      nc_selfcollect_tend   (0), // change in cloud droplet number from self-collection  (Not in paper?)
      nc2nr_autoconv_tend  (0), // change in cloud droplet number from autoconversion
      nr_selfcollect_tend   (0), // change in rain number from self-collection  (Not in paper?)
      nr_evap_tend   (0), // change in rain number from evaporation
      ncautr  (0), // change in rain number from autoconversion of cloud water

      // initialize ice-phase  process rates
      qi2qv_sublim_tend   (0), // sublimation of ice
      nr_ice_shed_tend  (0), // source for rain number from collision of rain/ice above freezing and shedding
      qc2qi_hetero_freeze_tend  (0), // immersion freezing droplets
      qr2qi_collect_tend   (0), // collection rain mass by ice
      qc2qr_ice_shed_tend   (0), // source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding
      qi2qr_melt_tend   (0), // melting of ice
      qc2qi_collect_tend   (0), // collection of cloud water by ice
      qr2qi_immers_freeze_tend  (0), // immersion freezing rain
      qv2qi_nucleat_tend   (0), // deposition/condensation freezing nuc
      ni2nr_melt_tend   (0), // melting of ice
      nc_collect_tend   (0), // change in cloud droplet number from collection by ice
      ncshdc  (0), // source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)
      nc2ni_immers_freeze_tend  (0), // immersion freezing droplets
      nr_collect_tend   (0), // change in rain number from collection by ice
      ni_selfcollect_tend   (0), // change in ice number from collection within a category (Not in paper?)
      ni_nucleat_tend   (0), // change in ice number from deposition/cond-freezing nucleation
      qv2qi_vapdep_tend   (0), // vapor deposition
      qc2qi_berg_tend  (0), // Bergeron process
      nr2ni_immers_freeze_tend  (0), // immersion freezing rain
      ni_sublim_tend   (0), // change in ice number from sublimation
      qc_growth_rate  (0), // wet growth rate

      // initialize time/space varying physical variables
      mu      (0), // TODO(doc)
      dv      (0), // TODO(doc)
      sc      (0), // TODO(doc)
      dqsdt   (0), // TODO(doc)
      dqsidt  (0), // TODO(doc)
      ab      (0), // TODO(doc)
      abi     (0), // TODO(doc)
      kap     (0), // TODO(doc)
      eii     (0), // temperature dependent aggregation efficiency

      // quantities related to process rates/parameters, interpolated from lookup tables:
      // For a more in depth reference to where these came from consult the file
      // "create_p3_lookupTable_1.F90-v4.1".  All line numbers below reference this
      // file.
      table_val_qi_fallspd(0), // mass-weighted fallspeed              See lines  731 -  808  ums
      table_val_ni_self_collect(0), // ice collection within a category     See lines  809 -  928  nagg
      table_val_qc2qi_collect(0), // collection of cloud water by ice     See lines  929 - 1009  nrwat
      table_val_qi2qr_melting(0), // melting                              See lines 1212 - 1279  vdep
      table_val_ice_eff_radius(0), // effective radius                     See lines 1281 - 1356  eff
      table_val_nr_collect(0), // collection of rain number by ice     See lines 1010 - 1209  nrrain
      table_val_qr2qi_collect(0), // collection of rain mass by ice       See lines 1010 - 1209  qrrain
      table_val_ni_lammax(0), // minimum ice number (lambda limiter)  See lines  704 -  705  nlarge
      table_val_ni_lammin(0), // maximum ice number (lambda limiter)  See lines  704 -  705  nsmall
      table_val_ice_reflectivity(0), // reflectivity                         See lines  731 -  808  refl
      table_val_qi2qr_vent_melt(0), // melting (ventilation term)           See lines 1212 - 1279  vdep1
      table_val_ice_mean_diam(0), // mass-weighted mean diameter          See lines 1212 - 1279  dmm
      table_val_ice_bulk_dens(0), // mass-weighted mean particle density  See lines 1212 - 1279  rhomm

      // TODO(doc)
      vtrmi1   (0),   // TODO(doc)
      rho_qm_cloud(400), // TODO(doc)
      epsi(0),        // TODO(doc)
      epsr(0),        // TODO(doc)
      epsc(0),        // TODO(doc)
      epsi_tot (0);   // inverse supersaturation relaxation timescale for combined ice categories

    Smask wetgrowth(false);

    // skip micro process calculations except nucleation/acvtivation if there no hydrometeors are present
    const auto skip_micro = skip_all || !(qc_incld(k) >= qsmall || qr_incld(k) >= qsmall || qi_incld(k) >= qsmall);
    const auto not_skip_micro = !skip_micro;

    if (not_skip_micro.any()) {
      // time/space varying physical variables
      get_time_space_phys_variables(
        T_atm(k), pres(k), rho(k), latent_heat_vapor(k), latent_heat_sublim(k), qv_sat_l(k), qv_sat_i(k),
        mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii, not_skip_micro);

      get_cloud_dsd2(qc_incld(k), nc_incld(k), mu_c(k), rho(k), nu(k), dnu,
                     lamc(k), cdist(k), cdist1(k), not_skip_micro);
      nc(k).set(not_skip_micro, nc_incld(k) * cld_frac_l(k));

      get_rain_dsd2(qr_incld(k), nr_incld(k), mu_r(k), lamr(k), cdistr(k), logn0r(k), not_skip_micro);
      nr(k).set(not_skip_micro, nr_incld(k) * cld_frac_r(k));

      impose_max_total_ni(ni_incld(k), max_total_ni, inv_rho(k), not_skip_micro);

      const auto qi_gt_small = qi_incld(k) >= qsmall && not_skip_micro;

      if (qi_gt_small.any()) {
        // impose lower limits to prevent taking log of # < 0
        ni_incld(k).set(qi_gt_small, max(ni_incld(k), nsmall));
        nr_incld(k).set(qi_gt_small, max(nr_incld(k), nsmall));

        const auto rhop = calc_bulk_rho_rime(qi_incld(k), qm_incld(k), bm_incld(k), qi_gt_small);
        qm(k).set(qi_gt_small, qm_incld(k)*cld_frac_i(k) );
        bm(k).set(qi_gt_small, bm_incld(k)*cld_frac_i(k) );

        TableIce table_ice;
        lookup_ice(qi_incld(k), ni_incld(k), qm_incld(k), rhop, table_ice, qi_gt_small);

        TableRain table_rain;
        lookup_rain(qr_incld(k), nr_incld(k), table_rain, qi_gt_small);

        // call to lookup table interpolation subroutines to get process rates
        table_val_qi_fallspd.set(qi_gt_small, apply_table_ice(1, ice_table_vals, table_ice, qi_gt_small));
        table_val_ni_self_collect.set(qi_gt_small, apply_table_ice(2, ice_table_vals, table_ice, qi_gt_small));
        table_val_qc2qi_collect.set(qi_gt_small, apply_table_ice(3, ice_table_vals, table_ice, qi_gt_small));
        table_val_qi2qr_melting.set(qi_gt_small, apply_table_ice(4, ice_table_vals, table_ice, qi_gt_small));
        table_val_ni_lammax.set(qi_gt_small, apply_table_ice(6, ice_table_vals, table_ice, qi_gt_small));
        table_val_ni_lammin.set(qi_gt_small, apply_table_ice(7, ice_table_vals, table_ice, qi_gt_small));
        table_val_qi2qr_vent_melt.set(qi_gt_small, apply_table_ice(9, ice_table_vals, table_ice, qi_gt_small));

        // ice-rain collection processes
        const auto qr_gt_small = qr_incld(k) >= qsmall && qi_gt_small;
        table_val_nr_collect.set(qr_gt_small, apply_table_coll(0, collect_table_vals, table_ice, table_rain, qi_gt_small));
        table_val_qr2qi_collect.set(qr_gt_small, apply_table_coll(1, collect_table_vals, table_ice, table_rain, qi_gt_small));

        // adjust Ni if needed to make sure mean size is in bounds (i.e. apply lambda limiters)
        // note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
        ni_incld(k).set(qi_gt_small, min(ni_incld(k), table_val_ni_lammax*ni_incld(k)));
        ni_incld(k).set(qi_gt_small, max(ni_incld(k), table_val_ni_lammin*ni_incld(k)));
      }

      // ----------------------------------------------------------------------
      // Begin calculations of microphysical processes
      // ----------------------------------------------------------------------

      // ......................................................................
      // ice processes
      // ......................................................................

      // collection of droplets
      ice_cldliq_collection(
        rho(k), T_atm(k), rhofaci(k), table_val_qc2qi_collect, qi_incld(k), qc_incld(k), ni_incld(k), nc_incld(k),
        qc2qi_collect_tend, nc_collect_tend, qc2qr_ice_shed_tend, ncshdc, not_skip_micro);

      // collection of rain
      ice_rain_collection(
        rho(k), T_atm(k), rhofaci(k), logn0r(k), table_val_nr_collect, table_val_qr2qi_collect, qi_incld(k), ni_incld(k), qr_incld(k),
        qr2qi_collect_tend, nr_collect_tend, not_skip_micro);

      // collection between ice categories

      // PMC nCat deleted lots of stuff here.

      // self-collection of ice
      ice_self_collection(
        rho(k), rhofaci(k), table_val_ni_self_collect, eii, qm_incld(k), qi_incld(k), ni_incld(k),
        ni_selfcollect_tend, not_skip_micro);

      // melting
      ice_melting(
        rho(k), T_atm(k), pres(k), rhofaci(k), table_val_qi2qr_melting, table_val_qi2qr_vent_melt, latent_heat_vapor(k), latent_heat_fusion(k), dv, sc, mu, kap, qv(k), qi_incld(k), ni_incld(k),
        qi2qr_melt_tend, ni2nr_melt_tend, not_skip_micro);

      // calculate wet growth
      ice_cldliq_wet_growth(
        rho(k), T_atm(k), pres(k), rhofaci(k), table_val_qi2qr_melting, table_val_qi2qr_vent_melt, latent_heat_vapor(k),
        latent_heat_fusion(k), dv, kap, mu, sc, qv(k), qc_incld(k), qi_incld(k), ni_incld(k), qr_incld(k),
        wetgrowth, qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend, not_skip_micro);

      // calculate total inverse ice relaxation timescale combined for all ice categories
      // note 'f1pr' values are normalized, so we need to multiply by N
      ice_relaxation_timescale(
        rho(k), T_atm(k), rhofaci(k), table_val_qi2qr_melting, table_val_qi2qr_vent_melt, dv, mu, sc, qi_incld(k), ni_incld(k),
        epsi, epsi_tot, not_skip_micro);

      // calculate rime density
      calc_rime_density(
        T_atm(k), rhofaci(k), table_val_qi_fallspd, acn(k), lamc(k), mu_c(k), qc_incld(k), qc2qi_collect_tend,
        vtrmi1, rho_qm_cloud, not_skip_micro);

      // contact and immersion freezing droplets
      cldliq_immersion_freezing(
        T_atm(k), lamc(k), mu_c(k), cdist1(k), qc_incld(k), inv_qc_relvar(k),
        qc2qi_hetero_freeze_tend, nc2ni_immers_freeze_tend, not_skip_micro);

      // for future: get rid of log statements below for rain freezing
      rain_immersion_freezing(
        T_atm(k), lamr(k), mu_r(k), cdistr(k), qr_incld(k),
        qr2qi_immers_freeze_tend, nr2ni_immers_freeze_tend, not_skip_micro);

      //  rime splintering (Hallet-Mossop 1974)
      // PMC comment: Morrison and Milbrandt 2015 part 1 and 2016 part 3 both say
      // that Hallet-Mossop should be neglected if 1 category to compensate for
      // artificial smearing out of ice DSD

      // ................................................
      //  condensation/evaporation/deposition/sublimation
      //    (use semi-analytic formulation)

      //  calculate rain evaporation including ventilation
      calc_liq_relaxation_timescale(
        revap_table_vals, rho(k), f1r, f2r, dv, mu, sc, mu_r(k), lamr(k), cdistr(k), cdist(k), qr_incld(k), qc_incld(k),
        epsr, epsc, not_skip_micro);

      evaporate_rain(qr_incld(k),qc_incld(k),nr_incld(k),qi_incld(k),
		     cld_frac_l(k),cld_frac_r(k),qv(k),qv_prev(k),qv_sat_l(k),qv_sat_i(k),
		     ab,abi,epsr,epsi_tot,T_atm(k),t_prev(k),latent_heat_sublim(k),dqsdt,dt,
		     qr2qv_evap_tend,nr_evap_tend, not_skip_micro);

      ice_deposition_sublimation(
	qi_incld(k), ni_incld(k), T_atm(k), qv_sat_l(k), qv_sat_i(k), epsi, abi, qv(k), inv_dt,
        qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend, not_skip_micro);
    }

    // deposition/condensation-freezing nucleation
    ice_nucleation(
      T_atm(k), inv_rho(k), ni(k), ni_activated(k), qv_supersat_i(k), inv_dt, predictNc, do_prescribed_CCN,
      qv2qi_nucleat_tend, ni_nucleat_tend, not_skip_all);

    // cloud water autoconversion
    // NOTE: cloud_water_autoconversion must be called before droplet_self_collection
    cloud_water_autoconversion(
      rho(k), qc_incld(k), nc_incld(k), inv_qc_relvar(k),
      qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, not_skip_all);

    // self-collection of droplets
    droplet_self_collection(
      rho(k), inv_rho(k), qc_incld(k),
      mu_c(k), nu(k), nc2nr_autoconv_tend, nc_selfcollect_tend, not_skip_all);

    // accretion of cloud by rain
    cloud_rain_accretion(
      rho(k), inv_rho(k), qc_incld(k), nc_incld(k), qr_incld(k), inv_qc_relvar(k),
      qc2qr_accret_tend, nc_accret_tend, not_skip_all);

    // self-collection and breakup of rain
    // (breakup following modified Verlinde and Cotton scheme)
    rain_self_collection(
      rho(k), qr_incld(k), nr_incld(k),
      nr_selfcollect_tend, not_skip_all);

    // Here we map the microphysics tendency rates back to CELL-AVERAGE quantities for updating
    // cell-average quantities.
    back_to_cell_average(
      cld_frac_l(k), cld_frac_r(k), cld_frac_i(k), qc2qr_accret_tend, qr2qv_evap_tend, qc2qr_autoconv_tend,
      nc_accret_tend, nc_selfcollect_tend, nc2nr_autoconv_tend, nr_selfcollect_tend, nr_evap_tend, ncautr, qi2qv_sublim_tend, nr_ice_shed_tend, qc2qi_hetero_freeze_tend,
      qr2qi_collect_tend, qc2qr_ice_shed_tend, qi2qr_melt_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, ni2nr_melt_tend, nc_collect_tend,
      ncshdc, nc2ni_immers_freeze_tend, nr_collect_tend, ni_selfcollect_tend,
      qv2qi_vapdep_tend, nr2ni_immers_freeze_tend, ni_sublim_tend, qv2qi_nucleat_tend, ni_nucleat_tend, qc2qi_berg_tend, not_skip_all);

    //
    // conservation of water
    //

    // don't bother checking vapor since all sinks already have limits imposed and the sum, therefore,
    // cannot possibly overdeplete qv [PMC: the above argument makes no sense to me. I think we don't
    // check qv because it is typically much greater than zero so seldom goes negative (and if it does
    // catastrophic failure is appropriate)]

    // cloud
    cloud_water_conservation(
      qc(k), dt,
      qc2qr_autoconv_tend, qc2qr_accret_tend, qc2qi_collect_tend, qc2qi_hetero_freeze_tend, qc2qr_ice_shed_tend, qc2qi_berg_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend, not_skip_all);

    // rain
    rain_water_conservation(
      qr(k), qc2qr_autoconv_tend, qc2qr_accret_tend, qi2qr_melt_tend, qc2qr_ice_shed_tend, dt,
      qr2qv_evap_tend, qr2qi_collect_tend, qr2qi_immers_freeze_tend, not_skip_all);

    // ice
    ice_water_conservation(
      qi(k), qv2qi_vapdep_tend, qv2qi_nucleat_tend, qc2qi_berg_tend, qr2qi_collect_tend, 
      qc2qi_collect_tend, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, dt,
      qi2qv_sublim_tend, qi2qr_melt_tend, not_skip_all);

    nc_conservation(nc(k), nc_selfcollect_tend, dt, nc_collect_tend, nc2ni_immers_freeze_tend,
                    nc_accret_tend, nc2nr_autoconv_tend, not_skip_all);
    nr_conservation(nr(k),ni2nr_melt_tend,nr_ice_shed_tend,ncshdc,nc2nr_autoconv_tend,dt,nmltratio,nr_collect_tend,
                    nr2ni_immers_freeze_tend,nr_selfcollect_tend,nr_evap_tend, not_skip_all);
    ni_conservation(ni(k),ni_nucleat_tend,nr2ni_immers_freeze_tend,nc2ni_immers_freeze_tend,dt,ni2nr_melt_tend,
                    ni_sublim_tend,ni_selfcollect_tend, not_skip_all);

    // make sure procs don't inappropriately push qv beyond ice saturation
    ice_supersat_conservation(qv2qi_vapdep_tend,qv2qi_nucleat_tend,cld_frac_i(k),qv(k),qv_sat_i(k),
			      latent_heat_sublim(k),th_atm(k)/inv_exner(k),dt,qi2qv_sublim_tend,qr2qv_evap_tend, not_skip_all);
    // make sure procs don't inappropriately push qv beyond liquid saturation
    prevent_liq_supersaturation(pres(k), T_atm(k), qv(k), latent_heat_vapor(k), latent_heat_sublim(k),dt,
				qv2qi_vapdep_tend, qv2qi_nucleat_tend, qi2qv_sublim_tend,qr2qv_evap_tend,
				not_skip_all);

    //---------------------------------------------------------------------------------
    // update prognostic microphysics and thermodynamics variables
    //---------------------------------------------------------------------------------

    //-- ice-phase dependent processes:
    update_prognostic_ice(
      qc2qi_hetero_freeze_tend, qc2qi_collect_tend, qc2qr_ice_shed_tend, nc_collect_tend, nc2ni_immers_freeze_tend, ncshdc, qr2qi_collect_tend, nr_collect_tend,  qr2qi_immers_freeze_tend,
      nr2ni_immers_freeze_tend, nr_ice_shed_tend, qi2qr_melt_tend, ni2nr_melt_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend, qv2qi_nucleat_tend, ni_nucleat_tend, ni_selfcollect_tend, ni_sublim_tend,
      qc2qi_berg_tend, inv_exner(k), latent_heat_sublim(k), latent_heat_fusion(k), predictNc, wetgrowth, dt, nmltratio,
      rho_qm_cloud, th_atm(k), qv(k), qi(k), ni(k), qm(k), bm(k), qc(k),
      nc(k), qr(k), nr(k), not_skip_all);

    //-- warm-phase only processes:
    update_prognostic_liquid(
      qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend, qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend,
      predictNc, do_prescribed_CCN, inv_rho(k), inv_exner(k), latent_heat_vapor(k), dt, th_atm(k), qv(k), qc(k), nc(k),
      qr(k), nr(k), not_skip_all);

    // AaronDonahue - Add extra variables needed from microphysics by E3SM:
    qv2qi_depos_tend(k)         .set(not_skip_all, qv2qi_vapdep_tend - qi2qv_sublim_tend + qv2qi_nucleat_tend);
    precip_total_tend(k)           .set(not_skip_all, qc2qr_accret_tend + qc2qr_autoconv_tend + qc2qr_ice_shed_tend + qc2qi_collect_tend);
    nevapr(k)          .set(not_skip_all, qi2qv_sublim_tend + qr2qv_evap_tend);
    qr_evap_tend(k)       .set(not_skip_all, qr2qv_evap_tend);
    vap_ice_exchange(k).set(not_skip_all, qv2qi_vapdep_tend - qi2qv_sublim_tend + qv2qi_nucleat_tend);
    vap_liq_exchange(k).set(not_skip_all, -qr2qv_evap_tend);
    liq_ice_exchange(k).set(not_skip_all, qc2qi_hetero_freeze_tend + qr2qi_immers_freeze_tend - qi2qr_melt_tend + qc2qi_berg_tend + qc2qi_collect_tend + qr2qi_collect_tend);

    // clipping for small hydrometeor values
    const auto qc_small    = qc(k) < qsmall    && not_skip_all;
    const auto qr_small    = qr(k) < qsmall    && not_skip_all;
    const auto qi_small = qi(k) < qsmall && not_skip_all;

    const auto qc_not_small    = qc(k) >= qsmall    && not_skip_all;
    const auto qr_not_small    = qr(k) >= qsmall    && not_skip_all;
    const auto qi_not_small = qi(k) >= qsmall && not_skip_all;

    qv(k).set(qc_small, qv(k) + qc(k));
    th_atm(k).set(qc_small, th_atm(k) - inv_exner(k) * qc(k) * latent_heat_vapor(k) * inv_cp);
    qc(k).set(qc_small, 0);
    nc(k).set(qc_small, 0);

    if (qc_not_small.any()) {
      hydrometeorsPresent = true;
    }

    qv(k).set(qr_small, qv(k) + qr(k));
    th_atm(k).set(qr_small, th_atm(k) - inv_exner(k) * qr(k) * latent_heat_vapor(k) * inv_cp);
    qr(k).set(qr_small, 0);
    nr(k).set(qr_small, 0);

    if (qr_not_small.any()) {
      hydrometeorsPresent = true;
    }

    qv(k).set(qi_small, qv(k) + qi(k));
    th_atm(k).set(qi_small, th_atm(k) - inv_exner(k) * qi(k) * latent_heat_sublim(k) * inv_cp);
    qi(k).set(qi_small, 0);
    ni(k).set(qi_small, 0);
    qm(k).set(qi_small, 0);
    bm(k).set(qi_small, 0);

    if (qi_not_small.any()) {
      hydrometeorsPresent = true;
    }

    // Outputs associated with aerocom comparison:
    pratot(k).set(not_skip_all, qc2qr_accret_tend); // cloud drop accretion by rain
    prctot(k).set(not_skip_all, qc2qr_autoconv_tend); // cloud drop autoconversion to rain

    //impose_max_total_ni is meant to operate on in-cloud vals. ni_incld is an output of
    //calculate_incloud_mixingratios below but we need to generate it earlier for impose_max_total_ni
    ni_incld(k).set(not_skip_all, ni(k)/cld_frac_i(k));
    impose_max_total_ni(ni_incld(k), max_total_ni, inv_rho(k), not_skip_all);
    ni(k).set(not_skip_all, ni_incld(k)*cld_frac_i(k));

    // Recalculate in-cloud values for sedimentation
    calculate_incloud_mixingratios(
      qc(k), qr(k), qi(k), qm(k), nc(k), nr(k), ni(k), bm(k), inv_cld_frac_l(k), inv_cld_frac_i(k), inv_cld_frac_r(k),
      qc_incld(k), qr_incld(k), qi_incld(k), qm_incld(k), nc_incld(k), nr_incld(k), ni_incld(k), bm_incld(k), not_skip_all);
  });
  team.team_barrier();
}

} // namespace p3
} // namespace scream

#endif // P3_MAIN_IMPL_PART_2_HPP
