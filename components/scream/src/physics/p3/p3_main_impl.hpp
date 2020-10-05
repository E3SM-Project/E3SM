#ifndef P3_MAIN_IMPL_HPP
#define P3_MAIN_IMPL_HPP

#include "physics/p3/p3_functions.hpp" // for ETI only but harmless for GPU
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs
#include "physics/share/physics_saturation_impl.hpp"

namespace scream {
namespace p3 {

/*
 * Implementation of p3 main function. Clients should NOT #include
 * this file, #include p3_functions.hpp instead.
 */

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main_init(
  const MemberType& team,
  const Int& nk_pack,
  const uview_1d<const Spack>& cld_frac_i,
  const uview_1d<const Spack>& cld_frac_l,
  const uview_1d<const Spack>& cld_frac_r,
  const uview_1d<const Spack>& exner,
  const uview_1d<const Spack>& th_atm,
  const uview_1d<const Spack>& dz,
  const uview_1d<Spack>& diag_equiv_reflectivity,
  const uview_1d<Spack>& ze_ice,
  const uview_1d<Spack>& ze_rain,
  const uview_1d<Spack>& diag_eff_radius_qc,
  const uview_1d<Spack>& diag_eff_radius_qi,
  const uview_1d<Spack>& inv_cld_frac_i,
  const uview_1d<Spack>& inv_cld_frac_l,
  const uview_1d<Spack>& inv_cld_frac_r,
  const uview_1d<Spack>& inv_exner,
  const uview_1d<Spack>& T_atm,
  const uview_1d<Spack>& qv,
  const uview_1d<Spack>& inv_dz,
  Scalar& precip_liq_surf,
  Scalar& precip_ice_surf,
  view_1d_ptr_array<Spack, 36>& zero_init)
{
  precip_liq_surf = 0;
  precip_ice_surf = 0;

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {

    diag_equiv_reflectivity(k)           = -99;
    ze_ice(k)            = 1.e-22;
    ze_rain(k)           = 1.e-22;
    diag_eff_radius_qc(k)         = 10.e-6;
    diag_eff_radius_qi(k)         = 25.e-6;
    inv_cld_frac_i(k)    = 1 / cld_frac_i(k);
    inv_cld_frac_l(k)    = 1 / cld_frac_l(k);
    inv_cld_frac_r(k)    = 1 / cld_frac_r(k);
    inv_exner(k)         = 1 / exner(k);
    T_atm(k)                 = th_atm(k) * inv_exner(k);
    qv(k)                = max(qv(k), 0);
    inv_dz(k)            = 1 / dz(k);

    for (size_t j = 0; j < zero_init.size(); ++j) {
      (*zero_init[j])(k) = 0;
    }
  });
  team.team_barrier();
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main_part1(
  const MemberType& team,
  const Int& nk,
  const bool& predictNc,
  const Scalar& dt,
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& dpres,
  const uview_1d<const Spack>& dz,
  const uview_1d<const Spack>& nc_nuceat_tend,
  const uview_1d<const Spack>& exner,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<const Spack>& inv_cld_frac_l,
  const uview_1d<const Spack>& inv_cld_frac_i,
  const uview_1d<const Spack>& inv_cld_frac_r,
  const uview_1d<const Spack>& latent_heat_vapor,
  const uview_1d<const Spack>& latent_heat_sublim,
  const uview_1d<const Spack>& latent_heat_fusion,
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
  const uview_1d<Spack>& qc_incld,
  const uview_1d<Spack>& qr_incld,
  const uview_1d<Spack>& qi_incld,
  const uview_1d<Spack>& qm_incld,
  const uview_1d<Spack>& nc_incld,
  const uview_1d<Spack>& nr_incld,
  const uview_1d<Spack>& ni_incld,
  const uview_1d<Spack>& bm_incld,
  bool& nucleationPossible,
  bool& hydrometeorsPresent)
{
  // Get access to saturation functions
  using physics = scream::physics::Functions<Scalar, Device>;

  // load constants into local vars
  constexpr Scalar g            = C::gravit;
  constexpr Scalar rho_1000mb   = C::RHO_1000MB;
  constexpr Scalar rho_600mb    = C::RHO_600MB;
  constexpr Scalar rho_h2o      = C::RHO_H2O;
  constexpr Scalar nccnst       = C::NCCNST;
  constexpr Scalar T_zerodegc     = C::T_zerodegc;
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar inv_cp       = C::INV_CP;

  nucleationPossible = false;
  hydrometeorsPresent = false;
  team.team_barrier();

  const Int nk_pack = ekat::npack<Spack>(nk);

  //
  // calculate some time-varying atmospheric variables
  // AaronDonahue - changed "rho" to be defined on nonhydrostatic
  // assumption, consistent with pressure based coordinate system
  //              - moved latent heat calculation to above.  Latent
  // heat is determined by calling a p3_util function so that it
  // can be made consistent with E3SM definition of latent heat
  //
  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {

    const auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);
    const auto range_mask = range_pack < nk;

    rho(k)          = dpres(k)/dz(k) / g;
    inv_rho(k)      = 1 / rho(k);
    qv_sat_l(k)     = physics::qv_sat(T_atm(k), pres(k), false, range_mask);
    qv_sat_i(k)     = physics::qv_sat(T_atm(k), pres(k), true,  range_mask);

    qv_supersat_i(k) = qv(k) / qv_sat_i(k) - 1;

    rhofacr(k) = pow(rho_1000mb * inv_rho(k), sp(.54));
    rhofaci(k) = pow(rho_600mb * inv_rho(k), sp(.54));
    Spack dum  = sp(1.496e-6) * pow(T_atm(k), sp(1.5)) / (T_atm(k) + 120); // this is mu
    acn(k)     = g * rho_h2o / (18 * dum); // 'a' parameter for droplet fallspeed (Stokes' law)

    if ( (T_atm(k) < T_zerodegc && qv_supersat_i(k) >= -0.05).any() ) {
      nucleationPossible = true;
    }

    // apply mass clipping if dry and mass is sufficiently small
    // (implying all mass is expected to evaporate/sublimate in one time step)
    auto drymass = qc(k) < qsmall;
    auto not_drymass = !drymass && range_mask;
    qv(k).set(drymass, qv(k) + qc(k));
    th_atm(k).set(drymass, th_atm(k) - exner(k) * qc(k) * latent_heat_vapor(k) * inv_cp);
    qc(k).set(drymass, 0);
    nc(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      hydrometeorsPresent = true; // updated further down
      // Apply droplet activation here (before other microphysical processes) for consistency with qc increase by saturation
      // adjustment already applied in macrophysics. If prescribed drop number is used, this is also a good place to
      // prescribe that value
      if (!predictNc) {
         nc(k).set(not_drymass, nccnst*inv_rho(k));
      } else {
         nc(k).set(not_drymass, max(nc(k) + nc_nuceat_tend(k) * dt, 0.0));
      }
    }

    drymass = qr(k) < qsmall;
    not_drymass = !drymass && range_mask;
    qv(k).set(drymass, qv(k) + qr(k));
    th_atm(k).set(drymass, th_atm(k) - exner(k) * qr(k) * latent_heat_vapor(k) * inv_cp);
    qr(k).set(drymass, 0);
    nr(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      hydrometeorsPresent = true; // updated further down
    }

    drymass = (qi(k) < qsmall || (qi(k) < 1.e-8 && qv_supersat_i(k) < -0.1));
    not_drymass = !drymass && range_mask;
    qv(k).set(drymass, qv(k) + qi(k));
    th_atm(k).set(drymass, th_atm(k) - exner(k) * qi(k) * latent_heat_sublim(k) * inv_cp);
    qi(k).set(drymass, 0);
    ni(k).set(drymass, 0);
    qm(k).set(drymass, 0);
    bm(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      hydrometeorsPresent = true; // final update
    }

    drymass = (qi(k) >= qsmall && qi(k) < 1.e-8 && T_atm(k) >= T_zerodegc);
    qr(k).set(drymass, qr(k) + qi(k));
    th_atm(k).set(drymass, th_atm(k) - exner(k) * qi(k) * latent_heat_fusion(k) * inv_cp);
    qi(k).set(drymass, 0);
    ni(k).set(drymass, 0);
    qm(k).set(drymass, 0);
    bm(k).set(drymass, 0);

    T_atm(k) = th_atm(k) * inv_exner(k);

    calculate_incloud_mixingratios(
      qc(k), qr(k), qi(k), qm(k), nc(k), nr(k), ni(k), bm(k),
      inv_cld_frac_l(k), inv_cld_frac_i(k), inv_cld_frac_r(k),
      qc_incld(k), qr_incld(k), qi_incld(k), qm_incld(k), nc_incld(k), nr_incld(k), ni_incld(k), bm_incld(k));
  });
  team.team_barrier();
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main_part2(
  const MemberType& team,
  const Int& nk_pack,
  const bool& predictNc,
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
  const uview_1d<const Spack>& exner,
  const uview_1d<const Spack>& inv_exner,
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
  constexpr Scalar T_zerodegc     = C::T_zerodegc;
  constexpr Scalar max_total_ni = C::max_total_ni;
  constexpr Scalar f1r          = C::f1r;
  constexpr Scalar f2r          = C::f2r;
  constexpr Scalar nmltratio    = C::nmltratio;
  constexpr Scalar inv_cp       = C::INV_CP;

  team.team_barrier();
  hydrometeorsPresent = false;
  team.team_barrier();

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {

    // if relatively dry and no hydrometeors at this level, skip to end of k-loop (i.e. skip this level)
    const auto skip_all = !(qc(k) >= qsmall || qr(k) >= qsmall || qi(k) >= qsmall) &&
      (T_atm(k) < T_zerodegc && qv_supersat_i(k) < -0.05);
    const auto not_skip_all = !skip_all;
    if (skip_all.all()) {
      return; // skip all process rates
    }

    //compute mask to identify padded values in packs, which are undefined
    const auto range_pack = ekat::range<IntSmallPack>(k*Spack::n);
    const auto range_mask = range_pack < nk;

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
        qi(k).set(qi_gt_small, qi_incld(k)*cld_frac_i(k) );
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
        qi2qr_melt_tend, ni2nr_melt_tend, range_mask, not_skip_micro);

      // calculate wet growth
      ice_cldliq_wet_growth(
        rho(k), T_atm(k), pres(k), rhofaci(k), table_val_qi2qr_melting, table_val_qi2qr_vent_melt, latent_heat_vapor(k), 
        latent_heat_fusion(k), dv, kap, mu, sc, qv(k), qc_incld(k), qi_incld(k), ni_incld(k), qr_incld(k),
        wetgrowth, qr2qi_collect_tend, qc2qi_collect_tend, qc_growth_rate, nr_ice_shed_tend, qc2qr_ice_shed_tend, range_mask, not_skip_micro);

      // calcualte total inverse ice relaxation timescale combined for all ice categories
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
        qi_incld(k), ni_incld(k), T_atm(k), qv_sat_l(k), qv_sat_i(k), epsi, abi, qv(k),
        qv2qi_vapdep_tend, qi2qv_sublim_tend, ni_sublim_tend, qc2qi_berg_tend, not_skip_micro);
    }

    // deposition/condensation-freezing nucleation
    ice_nucleation(
      T_atm(k), inv_rho(k), ni(k), ni_activated(k), qv_supersat_i(k), inv_dt, predictNc,
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

    // The microphysical process rates are computed above, based on the environmental conditions.
    // The rates are adjusted here (where necessary) such that the sum of the sinks of mass cannot
    // be greater than the sum of the sources, thereby resulting in overdepletion.
    //-- Limit ice process rates to prevent overdepletion of sources such that
    //   the subsequent adjustments are done with maximum possible rates for the
    //   time step.  (note: most ice rates are adjusted here since they must be done
    //   simultaneously (outside of iice-loops) to distribute reduction proportionally
    //   amongst categories.
    //PMC - might need to rethink above statement since only one category now.
    // AaronDonahue: Do we need the below checks for the new definition of
    // how qv2qi_vapdep_tend and qi2qv_sublim_tend are derived?
    // AaronDonahue: UPDATE, if we are following the implementation of MG
    // then the answer appears to be YES.  There is a similar check in MG
    // microphysics which limits qv2qi_vapdep_tend and qv2qi_nucleat_tend, but does not limit qi2qv_sublim_tend.
    // So similar but slightly different.  The overall answer would be that
    // qv2qi_vapdep_tend does need some limit.  The next questions are,
    //   1) Should we be taking qv2qi_nucleat_tend into consideration too?
    //   2) Is MG correct in NOT limiting qi2qv_sublim_tend?

    prevent_ice_overdepletion(pres(k), T_atm(k), qv(k), latent_heat_sublim(k), inv_dt, qv2qi_vapdep_tend, qi2qv_sublim_tend, range_mask, not_skip_all);

    // vapor -- not needed, since all sinks already have limits imposed and the sum, therefore,
    //          cannot possibly overdeplete qv

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
      qi(k), qv2qi_vapdep_tend, qv2qi_nucleat_tend, qc2qi_berg_tend, qr2qi_collect_tend, qc2qi_collect_tend, qr2qi_immers_freeze_tend, qc2qi_hetero_freeze_tend, dt,
      qi2qv_sublim_tend, qi2qr_melt_tend, not_skip_all);

    //---------------------------------------------------------------------------------
    // update prognostic microphysics and thermodynamics variables
    //---------------------------------------------------------------------------------

    //-- ice-phase dependent processes:
    update_prognostic_ice(
      qc2qi_hetero_freeze_tend, qc2qi_collect_tend, qc2qr_ice_shed_tend, nc_collect_tend, nc2ni_immers_freeze_tend, ncshdc, qr2qi_collect_tend, nr_collect_tend,  qr2qi_immers_freeze_tend,
      nr2ni_immers_freeze_tend, nr_ice_shed_tend, qi2qr_melt_tend, ni2nr_melt_tend, qi2qv_sublim_tend, qv2qi_vapdep_tend, qv2qi_nucleat_tend, ni_nucleat_tend, ni_selfcollect_tend, ni_sublim_tend,
      qc2qi_berg_tend, exner(k), latent_heat_sublim(k), latent_heat_fusion(k), predictNc, wetgrowth, dt, nmltratio,
      rho_qm_cloud, th_atm(k), qv(k), qi(k), ni(k), qm(k), bm(k), qc(k),
      nc(k), qr(k), nr(k), not_skip_all);

    //-- warm-phase only processes:
    update_prognostic_liquid(
      qc2qr_accret_tend, nc_accret_tend, qc2qr_autoconv_tend, nc2nr_autoconv_tend, ncautr, nc_selfcollect_tend, qr2qv_evap_tend, nr_evap_tend, nr_selfcollect_tend,
      predictNc, inv_rho(k), exner(k), latent_heat_vapor(k), dt, th_atm(k), qv(k), qc(k), nc(k),
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
    th_atm(k).set(qc_small, th_atm(k) - exner(k) * qc(k) * latent_heat_vapor(k) * inv_cp);
    qc(k).set(qc_small, 0);
    nc(k).set(qc_small, 0);

    if (qc_not_small.any()) {
      hydrometeorsPresent = true;
    }

    qv(k).set(qr_small, qv(k) + qr(k));
    th_atm(k).set(qr_small, th_atm(k) - exner(k) * qr(k) * latent_heat_vapor(k) * inv_cp);
    qr(k).set(qr_small, 0);
    nr(k).set(qr_small, 0);

    if (qr_not_small.any()) {
      hydrometeorsPresent = true;
    }

    qv(k).set(qi_small, qv(k) + qi(k));
    th_atm(k).set(qi_small, th_atm(k) - exner(k) * qi(k) * latent_heat_sublim(k) * inv_cp);
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

    impose_max_total_ni(ni(k), max_total_ni, inv_rho(k), not_skip_all);

    // Recalculate in-cloud values for sedimentation
    calculate_incloud_mixingratios(
      qc(k), qr(k), qi(k), qm(k), nc(k), nr(k), ni(k), bm(k), inv_cld_frac_l(k), inv_cld_frac_i(k), inv_cld_frac_r(k),
      qc_incld(k), qr_incld(k), qi_incld(k), qm_incld(k), nc_incld(k), nr_incld(k), ni_incld(k), bm_incld(k), not_skip_all);
  });
  team.team_barrier();
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main_part3(
  const MemberType& team,
  const Int& nk_pack,
  const view_dnu_table& dnu,
  const view_ice_table& ice_table_vals,
  const uview_1d<const Spack>& exner,
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
  const uview_1d<Spack>& latent_heat_vapor,
  const uview_1d<Spack>& latent_heat_sublim,
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
  const uview_1d<Spack>& diag_eff_radius_qc)
{
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar inv_cp       = C::INV_CP;
  constexpr Scalar max_total_ni = C::max_total_ni;
  constexpr Scalar nsmall       = C::NSMALL;

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {

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

      nc(k).set(qc_gt_small,nc_incld*cld_frac_l(k)); //cld_dsd2 might have changed incld nc... need consistency.
      diag_eff_radius_qc(k)       .set(qc_gt_small, sp(0.5) * (mu_c(k) + 3) / lamc(k));
      qv(k)              .set(qc_small, qv(k)+qc(k));
      th_atm(k)              .set(qc_small, th_atm(k)-exner(k)*qc(k)*latent_heat_vapor(k)*inv_cp);
      vap_liq_exchange(k).set(qc_small, vap_liq_exchange(k) - qc(k));
      qc(k)              .set(qc_small, 0);
      nc(k)              .set(qc_small, 0);
    }

    // Rain
    {
      const auto qr_gt_small = qr(k) >= qsmall;
      const auto qr_small    = !qr_gt_small;
      const auto qr_incld = qr(k)/cld_frac_r(k);
      auto nr_incld = nr(k)/cld_frac_r(k); //nr_incld is updated in get_rain_dsd2 but isn't used again

      get_rain_dsd2(
        qr_incld, nr_incld, mu_r(k), lamr(k), ignore1, ignore2, qr_gt_small);

      nr(k).set(qr_gt_small,nr_incld*cld_frac_r(k)); //rain_dsd2 might have changed incld nr... need consistency.

      //Note that integrating over the drop-size PDF as done here should only be done to in-cloud
      //quantities but radar reflectivity is likely meant to be a cell ave. Thus nr in the next line
      //really should be cld_frac_r * nr/cld_frac_r. Not doing that since cld_frac_r cancels out.
      ze_rain(k).set(qr_gt_small, nr(k)*(mu_r(k)+6)*(mu_r(k)+5)*(mu_r(k)+4)*
                     (mu_r(k)+3)*(mu_r(k)+2)*(mu_r(k)+1)/pow(lamr(k), sp(6.0))); // once f90 is gone, 6 can be int
      ze_rain(k).set(qr_gt_small, max(ze_rain(k), sp(1.e-22)));

      qv(k)              .set(qr_small, qv(k) + qr(k));
      th_atm(k)              .set(qr_small, th_atm(k) - exner(k)*qr(k)*latent_heat_vapor(k)*inv_cp);
      vap_liq_exchange(k).set(qr_small, vap_liq_exchange(k) - qr(k));
      qr(k)              .set(qr_small, 0);
      nr(k)              .set(qr_small, 0);
    }

    // Ice
    {
      //bugfix todo: this ice section should all use in-cloud values
      impose_max_total_ni(ni(k), max_total_ni, inv_rho(k));

      const auto qi_gt_small = qi(k) >= qsmall;
      const auto qi_small    = !qi_gt_small;

      // impose lower limits to prevent taking log of # < 0
      ni(k) = max(ni(k), nsmall);

      auto qi_incld = qi(k)/cld_frac_i(k);
      auto ni_incld = ni(k)/cld_frac_i(k);
      auto qm_incld = qm(k)/cld_frac_i(k);
      auto bm_incld = bm(k)/cld_frac_i(k);

      const auto rhop = calc_bulk_rho_rime(qi_incld, qm_incld, bm_incld, qi_gt_small);
      qi(k).set(qi_gt_small, qi_incld*cld_frac_i(k) );
      qm(k).set(qi_gt_small, qm_incld*cld_frac_i(k) );
      bm(k).set(qi_gt_small, bm_incld*cld_frac_i(k) );

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
      th_atm(k).set(qi_small, th_atm(k) - exner(k)*qi(k)*latent_heat_sublim(k)*inv_cp);
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

template <typename S, typename D>
void Functions<S,D>
::p3_main(
  const P3PrognosticState& prognostic_state,
  const P3DiagnosticInputs& diagnostic_inputs,
  const P3DiagnosticOutputs& diagnostic_outputs,
  const P3Infrastructure& infrastructure,
  const P3HistoryOnly& history_only,
  Int nj,
  Int nk)
{
  using ExeSpace = typename KT::ExeSpace;

  view_2d<Spack> latent_heat_sublim("latent_heat_sublim", nj, nk), latent_heat_vapor("latent_heat_vapor", nj, nk), latent_heat_fusion("latent_heat_fusion", nj, nk);

  get_latent_heat(nj, nk, latent_heat_vapor, latent_heat_sublim, latent_heat_fusion);

  const Int nk_pack = ekat::npack<Spack>(nk);
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);

  ekat::WorkspaceManager<Spack, Device> workspace_mgr(nk_pack, 47, policy);

  // load constants into local vars
  const     Scalar inv_dt          = 1 / infrastructure.dt;
  constexpr Int    kdir         = -1;
  const     Int    ktop         = kdir == -1 ? 0    : nk-1;
  const     Int    kbot         = kdir == -1 ? nk-1 : 0;
  constexpr bool   debug_ABORT  = false;

  // load tables
  view_1d_table mu_r_table_vals;
  view_2d_table vn_table_vals, vm_table_vals, revap_table_vals;
  view_ice_table ice_table_vals;
  view_collect_table collect_table_vals;
  view_dnu_table dnu;
  init_kokkos_ice_lookup_tables(ice_table_vals, collect_table_vals);
  init_kokkos_tables(vn_table_vals, vm_table_vals, revap_table_vals, mu_r_table_vals, dnu);

  // per-column bools
  view_2d<bool> bools("bools", nj, 2);

  // p3_main loop
  Kokkos::parallel_for(
    "p3 main loop",
    policy,
    KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();

    auto workspace = workspace_mgr.get_workspace(team);

    //
    // Get temporary workspaces needed for p3
    //
    uview_1d<Spack>
      mu_r,   // shape parameter of rain
      T_atm,      // temperature at the beginning of the microhpysics step [K]

      // 2D size distribution and fallspeed parameters
      lamr, logn0r, nu, cdist, cdist1, cdistr,

      // Variables needed for in-cloud calculations
      inv_cld_frac_i, inv_cld_frac_l, inv_cld_frac_r, // Inverse cloud fractions (1/cld)
      qc_incld, qr_incld, qi_incld, qm_incld, // In cloud mass-mixing ratios
      nc_incld, nr_incld, ni_incld, bm_incld, // In cloud number concentrations

      // Other
      inv_dz, inv_rho, ze_ice, ze_rain, prec, rho,
      rhofacr, rhofaci, acn, qv_sat_l, qv_sat_i, sup, qv_supersat_i,
      tmparr1, inv_exner, diag_equiv_reflectivity, diag_vm_qi, diag_diam_qi, pratot, prctot,

      // p3_tend_out, may not need these
      qtend_ignore, ntend_ignore;

    workspace.template take_many_and_reset<41>(
      {
        "mu_r", "T_atm", "lamr", "logn0r", "nu", "cdist", "cdist1", "cdistr",
        "inv_cld_frac_i", "inv_cld_frac_l", "inv_cld_frac_r", "qc_incld", "qr_incld", "qi_incld", "qm_incld",
        "nc_incld", "nr_incld", "ni_incld", "bm_incld",
        "inv_dz", "inv_rho", "ze_ice", "ze_rain", "prec", "rho",
        "rhofacr", "rhofaci", "acn", "qv_sat_l", "qv_sat_i", "sup", "qv_supersat_i",
        "tmparr1", "inv_exner", "diag_equiv_reflectivity", "diag_vm_qi", "diag_diam_qi",
        "pratot", "prctot", "qtend_ignore", "ntend_ignore"
      },
      {
        &mu_r, &T_atm, &lamr, &logn0r, &nu, &cdist, &cdist1, &cdistr,
        &inv_cld_frac_i, &inv_cld_frac_l, &inv_cld_frac_r, &qc_incld, &qr_incld, &qi_incld, &qm_incld,
        &nc_incld, &nr_incld, &ni_incld, &bm_incld,
        &inv_dz, &inv_rho, &ze_ice, &ze_rain, &prec, &rho,
        &rhofacr, &rhofaci, &acn, &qv_sat_l, &qv_sat_i, &sup, &qv_supersat_i,
        &tmparr1, &inv_exner, &diag_equiv_reflectivity, &diag_vm_qi, &diag_diam_qi,
        &pratot, &prctot, &qtend_ignore, &ntend_ignore
      });

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto opres               = ekat::subview(diagnostic_inputs.pres, i);
    const auto odz                 = ekat::subview(diagnostic_inputs.dz, i);
    const auto onc_nuceat_tend     = ekat::subview(diagnostic_inputs.nc_nuceat_tend, i);
    const auto oni_activated       = ekat::subview(diagnostic_inputs.ni_activated, i);
    const auto oinv_qc_relvar      = ekat::subview(diagnostic_inputs.inv_qc_relvar, i);
    const auto odpres              = ekat::subview(diagnostic_inputs.dpres, i);
    const auto oexner              = ekat::subview(diagnostic_inputs.exner, i);
    const auto ocld_frac_i         = ekat::subview(diagnostic_inputs.cld_frac_i, i);
    const auto ocld_frac_l         = ekat::subview(diagnostic_inputs.cld_frac_l, i);
    const auto ocld_frac_r         = ekat::subview(diagnostic_inputs.cld_frac_r, i);
    const auto ocol_location       = ekat::subview(infrastructure.col_location, i);
    const auto oqc                 = ekat::subview(prognostic_state.qc, i);
    const auto onc                 = ekat::subview(prognostic_state.nc, i);
    const auto oqr                 = ekat::subview(prognostic_state.qr, i);
    const auto onr                 = ekat::subview(prognostic_state.nr, i);
    const auto oqi                 = ekat::subview(prognostic_state.qi, i);
    const auto oqm                 = ekat::subview(prognostic_state.qm, i);
    const auto oni                 = ekat::subview(prognostic_state.ni, i);
    const auto obm                 = ekat::subview(prognostic_state.bm, i);
    const auto oqv                 = ekat::subview(prognostic_state.qv, i);
    const auto oth                 = ekat::subview(prognostic_state.th, i);
    const auto odiag_eff_radius_qc    = ekat::subview(diagnostic_outputs.diag_eff_radius_qc, i);
    const auto odiag_eff_radius_qi    = ekat::subview(diagnostic_outputs.diag_eff_radius_qi, i);
    const auto orho_qi             = ekat::subview(diagnostic_outputs.rho_qi, i);
    const auto omu_c               = ekat::subview(diagnostic_outputs.mu_c, i);
    const auto olamc               = ekat::subview(diagnostic_outputs.lamc, i);
    const auto oqv2qi_depos_tend   = ekat::subview(diagnostic_outputs.qv2qi_depos_tend, i);
    const auto oprecip_total_tend  = ekat::subview(diagnostic_outputs.precip_total_tend, i);
    const auto onevapr             = ekat::subview(diagnostic_outputs.nevapr, i);
    const auto oqr_evap_tend       = ekat::subview(diagnostic_outputs.qr_evap_tend, i);
    const auto oprecip_liq_flux    = ekat::subview(diagnostic_outputs.precip_liq_flux, i);
    const auto oprecip_ice_flux    = ekat::subview(diagnostic_outputs.precip_ice_flux, i);
    const auto oliq_ice_exchange   = ekat::subview(history_only.liq_ice_exchange, i);
    const auto ovap_liq_exchange   = ekat::subview(history_only.vap_liq_exchange, i);
    const auto ovap_ice_exchange   = ekat::subview(history_only.vap_ice_exchange, i);
    const auto olatent_heat_vapor  = ekat::subview(latent_heat_vapor, i);
    const auto olatent_heat_sublim = ekat::subview(latent_heat_sublim, i);
    const auto olatent_heat_fusion = ekat::subview(latent_heat_fusion, i);
    const auto oqv_prev            = ekat::subview(diagnostic_inputs.qv_prev, i);
    const auto ot_prev             = ekat::subview(diagnostic_inputs.t_prev, i);

    // Need to watch out for race conditions with these shared variables
    bool &nucleationPossible  = bools(i, 0);
    bool &hydrometeorsPresent = bools(i, 1);

    view_1d_ptr_array<Spack, 36> zero_init = {
      &mu_r, &lamr, &logn0r, &nu, &cdist, &cdist1, &cdistr,
      &qc_incld, &qr_incld, &qi_incld, &qm_incld,
      &nc_incld, &nr_incld, &ni_incld, &bm_incld,
      &inv_rho, &prec, &rho, &rhofacr, &rhofaci, &acn, &qv_sat_l, &qv_sat_i, &sup, &qv_supersat_i,
      &tmparr1, &qtend_ignore, &ntend_ignore,
      &omu_c, &olamc, &orho_qi, &oqv2qi_depos_tend, &oprecip_total_tend, &onevapr, &oprecip_liq_flux, &oprecip_ice_flux
    };

    // initialize
    p3_main_init(
      team, nk_pack,
      ocld_frac_i, ocld_frac_l, ocld_frac_r, oexner, oth, odz, diag_equiv_reflectivity,
      ze_ice, ze_rain, odiag_eff_radius_qc, odiag_eff_radius_qi, inv_cld_frac_i, inv_cld_frac_l,
      inv_cld_frac_r, inv_exner, T_atm, oqv, inv_dz,
      diagnostic_outputs.precip_liq_surf(i), diagnostic_outputs.precip_ice_surf(i), zero_init);

    p3_main_part1(
      team, nk, infrastructure.predictNc, infrastructure.dt,
      opres, odpres, odz, onc_nuceat_tend, oexner, inv_exner, inv_cld_frac_l, inv_cld_frac_i,
      inv_cld_frac_r, olatent_heat_vapor, olatent_heat_sublim, olatent_heat_fusion, T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr,
      rhofaci, acn, oqv, oth, oqc, onc, oqr, onr, oqi, oni, oqm,
      obm, qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld,
      ni_incld, bm_incld, nucleationPossible, hydrometeorsPresent);

    // There might not be any work to do for this team
    if (!(nucleationPossible || hydrometeorsPresent)) {
      return; // this is how you do a "continue" in a kokkos lambda
    }
    
    // ------------------------------------------------------------------------------------------
    // main k-loop (for processes):
    
    p3_main_part2(
      team, nk_pack, infrastructure.predictNc, infrastructure.dt, inv_dt,
      dnu, ice_table_vals, collect_table_vals, revap_table_vals, opres, odpres, odz, onc_nuceat_tend, oexner,
      inv_exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, oni_activated, oinv_qc_relvar, ocld_frac_i,
      ocld_frac_l, ocld_frac_r, oqv_prev, ot_prev, T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn,
      oqv, oth, oqc, onc, oqr, onr, oqi, oni, oqm, obm, olatent_heat_vapor,
      olatent_heat_sublim, olatent_heat_fusion, qc_incld, qr_incld, qi_incld, qm_incld, nc_incld,
      nr_incld, ni_incld, bm_incld, omu_c, nu, olamc, cdist, cdist1, cdistr,
      mu_r, lamr, logn0r, oqv2qi_depos_tend, oprecip_total_tend, onevapr, oqr_evap_tend,
      ovap_liq_exchange, ovap_ice_exchange, oliq_ice_exchange,
      pratot, prctot, hydrometeorsPresent, nk);

    //NOTE: At this point, it is possible to have negative (but small) nc, nr, ni.  This is not
    //      a problem; those values get clipped to zero in the sedimentation section (if necessary).
    //      (This is not done above simply for efficiency purposes.)
    
    if (!hydrometeorsPresent) return;

    // -----------------------------------------------------------------------------------------
    // End of main microphysical processes section
    // =========================================================================================

    // ==========================================================================================!
    // Sedimentation:

    // Cloud sedimentation:  (adaptive substepping)

    cloud_sedimentation(
      qc_incld, rho, inv_rho, ocld_frac_l, acn, inv_dz, dnu, team, workspace,
      nk, ktop, kbot, kdir, infrastructure.dt, inv_dt, infrastructure.predictNc,
      oqc, onc, nc_incld, omu_c, olamc, qtend_ignore, ntend_ignore,
      diagnostic_outputs.precip_liq_surf(i));

    // Rain sedimentation:  (adaptive substepping)
    rain_sedimentation(
      rho, inv_rho, rhofacr, ocld_frac_r, inv_dz, qr_incld, team, workspace,
      vn_table_vals, vm_table_vals, nk, ktop, kbot, kdir, infrastructure.dt, inv_dt, oqr,
      onr, nr_incld, mu_r, lamr, oprecip_liq_flux, qtend_ignore, ntend_ignore,
      diagnostic_outputs.precip_liq_surf(i));

    // Ice sedimentation:  (adaptive substepping)
    ice_sedimentation(
      rho, inv_rho, rhofaci, ocld_frac_i, inv_dz, team, workspace, nk, ktop, kbot,
      kdir, infrastructure.dt, inv_dt, oqi, qi_incld, oni, ni_incld,
      oqm, qm_incld, obm, bm_incld, qtend_ignore, ntend_ignore,
      ice_table_vals, diagnostic_outputs.precip_ice_surf(i));

    // homogeneous freezing of cloud and rain
    homogeneous_freezing(
      T_atm, oexner, olatent_heat_fusion, team, nk, ktop, kbot, kdir, oqc, onc, oqr, onr, oqi,
      oni, oqm, obm, oth);

    //
    // final checks to ensure consistency of mass/number
    // and compute diagnostic fields for output
    //
    p3_main_part3(
      team, nk_pack, dnu, ice_table_vals, oexner, ocld_frac_l, ocld_frac_r, ocld_frac_i, 
      rho, inv_rho, rhofaci, oqv, oth, oqc, onc, oqr, onr, oqi, oni,
      oqm, obm, olatent_heat_vapor, olatent_heat_sublim, omu_c, nu, olamc, mu_r, lamr,
      ovap_liq_exchange, ze_rain, ze_ice, diag_vm_qi, odiag_eff_radius_qi, diag_diam_qi,
      orho_qi, diag_equiv_reflectivity, odiag_eff_radius_qc);

    //
    // merge ice categories with similar properties

    //   note:  this should be relocated to above, such that the diagnostic
    //          ice properties are computed after merging

    // PMC nCat deleted nCat>1 stuff

#ifndef NDEBUG
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {
        tmparr1(k) = oth(k) * inv_exner(k);
    });

    check_values(oqv, tmparr1, ktop, kbot, infrastructure.it, debug_ABORT, 900,
                 team, ocol_location);
#endif
  });
}

} // namespace p3
} // namespace scream

#endif // P3_MAIN_IMPL_HPP
