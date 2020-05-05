#ifndef P3_FUNCTIONS_MAIN_IMPL_HPP
#define P3_FUNCTIONS_MAIN_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU

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
  const uview_1d<const Spack>& oicldm,
  const uview_1d<const Spack>& olcldm,
  const uview_1d<const Spack>& orcldm,
  const uview_1d<const Spack>& oexner,
  const uview_1d<const Spack>& oth,
  const uview_1d<Spack>& opratot,
  const uview_1d<Spack>& oprctot,
  const uview_1d<Spack>& prec,
  const uview_1d<Spack>& mu_r,
  const uview_1d<Spack>& odiag_ze,
  const uview_1d<Spack>& ze_ice,
  const uview_1d<Spack>& ze_rain,
  const uview_1d<Spack>& odiag_effc,
  const uview_1d<Spack>& odiag_effi,
  const uview_1d<Spack>& odiag_vmi,
  const uview_1d<Spack>& odiag_di,
  const uview_1d<Spack>& odiag_rhoi,
  const uview_1d<Spack>& ocmeiout,
  const uview_1d<Spack>& oprain,
  const uview_1d<Spack>& onevapr,
  const uview_1d<Spack>& orflx,
  const uview_1d<Spack>& osflx,
  const uview_1d<Spack>& inv_icldm,
  const uview_1d<Spack>& inv_lcldm,
  const uview_1d<Spack>& inv_rcldm,
  const uview_1d<Spack>& omu_c,
  const uview_1d<Spack>& olamc,
  const uview_1d<Spack>& inv_exner,
  const uview_1d<Spack>& t,
  const uview_1d<Spack>& oqv,
  Scalar& prt_liq,
  Scalar& prt_sol)
{
  prt_liq = 0;
  prt_sol = 0;

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {

    opratot(k)    = 0;
    oprctot(k)    = 0;
    prec(k)       = 0;
    mu_r(k)       = 0;
    odiag_ze(k)   = -99;
    ze_ice(k)     = 1.e-22;
    ze_rain(k)    = 1.e-22;
    odiag_effc(k) = 10.e-6;
    odiag_effi(k) = 25.e-6;
    odiag_vmi(k)  = 0;
    odiag_di(k)   = 0;
    odiag_rhoi(k) = 0;
    ocmeiout(k)   = 0;
    oprain(k)     = 0;
    onevapr(k)    = 0;
    orflx(k)      = 0;
    osflx(k)      = 0;
    inv_icldm(k)  = 1 / oicldm(k);
    inv_lcldm(k)  = 1 / olcldm(k);
    inv_rcldm(k)  = 1 / orcldm(k);
    omu_c(k)      = 0;
    olamc(k)      = 0;
    inv_exner(k)  = 1 / oexner(k);
    t(k)          = oth(k) * inv_exner(k);
    oqv(k)        = pack::max(oqv(k), 0);
  });
  team.team_barrier();
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main_pre_main_loop(
  const MemberType& team,
  const Int& nk,
  const bool& log_predictNc,
  const Scalar& dt,
  const uview_1d<const Spack>& opres,
  const uview_1d<const Spack>& opdel,
  const uview_1d<const Spack>& odzq,
  const uview_1d<const Spack>& onpccn,
  const uview_1d<const Spack>& oexner,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<const Spack>& inv_lcldm,
  const uview_1d<const Spack>& inv_icldm,
  const uview_1d<const Spack>& inv_rcldm,
  const uview_1d<const Spack>& oxxlv,
  const uview_1d<const Spack>& oxxls,
  const uview_1d<const Spack>& oxlf,
  const uview_1d<Spack>& t,
  const uview_1d<Spack>& rho,
  const uview_1d<Spack>& inv_rho,
  const uview_1d<Spack>& qvs,
  const uview_1d<Spack>& qvi,
  const uview_1d<Spack>& sup,
  const uview_1d<Spack>& supi,
  const uview_1d<Spack>& rhofacr,
  const uview_1d<Spack>& rhofaci,
  const uview_1d<Spack>& acn,
  const uview_1d<Spack>& oqv,
  const uview_1d<Spack>& oth,
  const uview_1d<Spack>& oqc,
  const uview_1d<Spack>& onc,
  const uview_1d<Spack>& oqr,
  const uview_1d<Spack>& onr,
  const uview_1d<Spack>& oqitot,
  const uview_1d<Spack>& onitot,
  const uview_1d<Spack>& oqirim,
  const uview_1d<Spack>& obirim,
  const uview_1d<Spack>& qc_incld,
  const uview_1d<Spack>& qr_incld,
  const uview_1d<Spack>& qitot_incld,
  const uview_1d<Spack>& qirim_incld,
  const uview_1d<Spack>& nc_incld,
  const uview_1d<Spack>& nr_incld,
  const uview_1d<Spack>& nitot_incld,
  const uview_1d<Spack>& birim_incld,
  bool& log_nucleationPossible,
  bool& log_hydrometeorsPresent)
{
  // load constants into local vars
  constexpr Scalar g            = C::gravit;
  constexpr Scalar rhosur       = C::RHOSUR;
  constexpr Scalar rhosui       = C::rhosui;
  constexpr Scalar rhow         = C::RHOW;
  constexpr Scalar nccnst       = C::NCCNST;
  constexpr Scalar zerodegc     = C::ZeroDegC;
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar inv_cp       = C::INV_CP;

  log_nucleationPossible = false;
  log_hydrometeorsPresent = false;
  team.team_barrier();

  const Int nk_pack = scream::pack::npack<Spack>(nk);

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

    const auto range_pack = scream::pack::range<IntSmallPack>(k*Spack::n);
    const auto range_mask = range_pack < nk;

    rho(k)     = opdel(k)/odzq(k) / g;
    inv_rho(k) = 1 / rho(k);
    qvs(k)     = qv_sat(t(k), opres(k), 0);
    qvi(k)     = qv_sat(t(k), opres(k), 1);

    sup(k)  = oqv(k) / qvs(k) - 1;
    supi(k) = oqv(k) / qvi(k) - 1;

    rhofacr(k) = pack::pow(rhosur * inv_rho(k), sp(.54));
    rhofaci(k) = pack::pow(rhosui * inv_rho(k), sp(.54));
    Spack dum  = sp(1.496e-6) * pack::pow(t(k), sp(1.5)) / (t(k) + 120); // this is mu
    acn(k)     = g * rhow / (18 * dum); // 'a' parameter for droplet fallspeed (Stokes' law)

    // specify cloud droplet number (for 1-moment version)
    if (!log_predictNc) {
      onc(k) = nccnst * inv_rho(k);
    }

    if ( ( (t(k) < zerodegc && supi(k) >= 0.05) ||
           (t(k) >= zerodegc && sup(k) >= 0.05) ).any() ) {
      log_nucleationPossible = true;
    }

    // apply mass clipping if dry and mass is sufficiently small
    // (implying all mass is expected to evaporate/sublimate in one time step)
    auto drymass = (oqc(k) < qsmall || (oqc(k) < 1.e-8 && sup(k) < -0.1));
    auto not_drymass = !drymass && range_mask;
    oqv(k).set(drymass, oqv(k) + oqc(k));
    oth(k).set(drymass, oth(k) - oexner(k) * oqc(k) * oxxlv(k) * inv_cp);
    oqc(k).set(drymass, 0);
    onc(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      log_hydrometeorsPresent = true; // updated further down
    }

    drymass = (oqr(k) < qsmall || (oqr(k) < 1.e-8 && sup(k) < -0.1));
    not_drymass = !drymass && range_mask;
    oqv(k).set(drymass, oqv(k) + oqr(k));
    oth(k).set(drymass, oth(k) - oexner(k) * oqr(k) * oxxlv(k) * inv_cp);
    oqr(k).set(drymass, 0);
    onr(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      log_hydrometeorsPresent = true; // updated further down
    }

    drymass = (oqitot(k) < qsmall || (oqitot(k) < 1.e-8 && supi(k) < -0.1));
    not_drymass = !drymass && range_mask;
    oqv(k).set(drymass, oqv(k) + oqitot(k));
    oth(k).set(drymass, oth(k) - oexner(k) * oqitot(k) * oxxls(k) * inv_cp);
    oqitot(k).set(drymass, 0);
    onitot(k).set(drymass, 0);
    oqirim(k).set(drymass, 0);
    obirim(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      log_hydrometeorsPresent = true; // final update
    }

    drymass = (oqitot(k) >= qsmall && oqitot(k) < 1.e-8 && t(k) >= zerodegc);
    oqr(k).set(drymass, oqr(k) + oqitot(k));
    oth(k).set(drymass, oth(k) - oexner(k) * oqitot(k) * oxlf(k) * inv_cp);
    oqitot(k).set(drymass, 0);
    onitot(k).set(drymass, 0);
    oqirim(k).set(drymass, 0);
    obirim(k).set(drymass, 0);

    t(k) = oth(k) * inv_exner(k);

    // Activation of cloud droplets
    if (log_predictNc) {
      onc(k) += onpccn(k) * dt;
    }

    calculate_incloud_mixingratios(
      oqc(k), oqr(k), oqitot(k), oqirim(k), onc(k), onr(k), onitot(k), obirim(k),
      inv_lcldm(k), inv_icldm(k), inv_rcldm(k),
      qc_incld(k), qr_incld(k), qitot_incld(k), qirim_incld(k), nc_incld(k), nr_incld(k), nitot_incld(k), birim_incld(k));
  });
  team.team_barrier();
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main_main_loop(
  const MemberType& team,
  const Int& nk_pack,
  const bool& log_predictNc,
  const Scalar& dt,
  const Scalar& odt,
  const view_dnu_table& dnu,
  const view_itab_table& itab,
  const view_itabcol_table& itabcol,
  const view_2d_table& revap_table,
  const uview_1d<const Spack>& opres,
  const uview_1d<const Spack>& opdel,
  const uview_1d<const Spack>& odzq,
  const uview_1d<const Spack>& onpccn,
  const uview_1d<const Spack>& oexner,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<const Spack>& inv_lcldm,
  const uview_1d<const Spack>& inv_icldm,
  const uview_1d<const Spack>& inv_rcldm,
  const uview_1d<const Spack>& onaai,
  const uview_1d<const Spack>& oicldm,
  const uview_1d<const Spack>& olcldm,
  const uview_1d<const Spack>& orcldm,
  const uview_1d<Spack>& t,
  const uview_1d<Spack>& rho,
  const uview_1d<Spack>& inv_rho,
  const uview_1d<Spack>& qvs,
  const uview_1d<Spack>& qvi,
  const uview_1d<Spack>& sup,
  const uview_1d<Spack>& supi,
  const uview_1d<Spack>& rhofacr,
  const uview_1d<Spack>& rhofaci,
  const uview_1d<Spack>& acn,
  const uview_1d<Spack>& oqv,
  const uview_1d<Spack>& oth,
  const uview_1d<Spack>& oqc,
  const uview_1d<Spack>& onc,
  const uview_1d<Spack>& oqr,
  const uview_1d<Spack>& onr,
  const uview_1d<Spack>& oqitot,
  const uview_1d<Spack>& onitot,
  const uview_1d<Spack>& oqirim,
  const uview_1d<Spack>& obirim,
  const uview_1d<Spack>& oxxlv,
  const uview_1d<Spack>& oxxls,
  const uview_1d<Spack>& oxlf,
  const uview_1d<Spack>& qc_incld,
  const uview_1d<Spack>& qr_incld,
  const uview_1d<Spack>& qitot_incld,
  const uview_1d<Spack>& qirim_incld,
  const uview_1d<Spack>& nc_incld,
  const uview_1d<Spack>& nr_incld,
  const uview_1d<Spack>& nitot_incld,
  const uview_1d<Spack>& birim_incld,
  const uview_1d<Spack>& omu_c,
  const uview_1d<Spack>& nu,
  const uview_1d<Spack>& olamc,
  const uview_1d<Spack>& cdist,
  const uview_1d<Spack>& cdist1,
  const uview_1d<Spack>& cdistr,
  const uview_1d<Spack>& mu_r,
  const uview_1d<Spack>& lamr,
  const uview_1d<Spack>& logn0r,
  const uview_1d<Spack>& ocmeiout,
  const uview_1d<Spack>& oprain,
  const uview_1d<Spack>& onevapr,
  const uview_1d<Spack>& oprer_evap,
  const uview_1d<Spack>& ovap_cld_exchange,
  const uview_1d<Spack>& ovap_liq_exchange,
  const uview_1d<Spack>& ovap_ice_exchange,
  const uview_1d<Spack>& oliq_ice_exchange,
  const uview_1d<Spack>& opratot,
  const uview_1d<Spack>& oprctot,
  bool& log_hydrometeorsPresent)
{
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar nsmall       = C::NSMALL;
  constexpr Scalar zerodegc     = C::ZeroDegC;
  constexpr Scalar max_total_Ni = C::max_total_Ni;
  constexpr Scalar f1r          = C::f1r;
  constexpr Scalar f2r          = C::f2r;
  constexpr Scalar nmltratio    = C::nmltratio;
  constexpr Scalar inv_cp       = C::INV_CP;

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {

    // if relatively dry and no hydrometeors at this level, skip to end of k-loop (i.e. skip this level)
    const auto skip_all = !(oqc(k) >= qsmall || oqr(k) >= qsmall || oqitot(k) >= qsmall) &&
      ( (t(k) < zerodegc && supi(k) < -0.05) || (t(k) >= zerodegc && sup(k) < -0.05) );
    const auto not_skip_all = !skip_all;
    if (skip_all.all()) {
      return; // skip all process rates
    }

    // All microphysics tendencies will be computed as IN-CLOUD, they will be mapped back to cell-average later.

    Spack
      // initialize warm-phase process rates
      qcacc   (0), // cloud droplet accretion by rain
      qrevp   (0), // rain evaporation
      qcaut   (0), // cloud droplet autoconversion to rain
      ncacc   (0), // change in cloud droplet number from accretion by rain
      ncnuc   (0), // change in cloud droplet number from activation of CCN
      ncslf   (0), // change in cloud droplet number from self-collection  (Not in paper?)
      ncautc  (0), // change in cloud droplet number from autoconversion
      qcnuc   (0), // activation of cloud droplets from CCN
      nrslf   (0), // change in rain number from self-collection  (Not in paper?)
      nrevp   (0), // change in rain number from evaporation
      ncautr  (0), // change in rain number from autoconversion of cloud water

      // initialize ice-phase  process rates
      qisub   (0), // sublimation of ice
      nrshdr  (0), // source for rain number from collision of rain/ice above freezing and shedding
      qcheti  (0), // immersion freezing droplets
      qrcol   (0), // collection rain mass by ice
      qcshd   (0), // source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding
      qimlt   (0), // melting of ice
      qccol   (0), // collection of cloud water by ice
      qrheti  (0), // immersion freezing rain
      qinuc   (0), // deposition/condensation freezing nuc
      nimlt   (0), // melting of ice
      nccol   (0), // change in cloud droplet number from collection by ice
      ncshdc  (0), // source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)
      ncheti  (0), // immersion freezing droplets
      nrcol   (0), // change in rain number from collection by ice
      nislf   (0), // change in ice number from collection within a category (Not in paper?)
      ninuc   (0), // change in ice number from deposition/cond-freezing nucleation
      qidep   (0), // vapor deposition
      qiberg  (0), // Bergeron process
      nrheti  (0), // immersion freezing rain
      nisub   (0), // change in ice number from sublimation
      qwgrth  (0), // wet growth rate

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
      f1pr02(0), // mass-weighted fallspeed              See lines  731 -  808  ums
      f1pr03(0), // ice collection within a category     See lines  809 -  928  nagg
      f1pr04(0), // collection of cloud water by ice     See lines  929 - 1009  nrwat
      f1pr05(0), // melting                              See lines 1212 - 1279  vdep
      f1pr06(0), // effective radius                     See lines 1281 - 1356  eff
      f1pr07(0), // collection of rain number by ice     See lines 1010 - 1209  nrrain
      f1pr08(0), // collection of rain mass by ice       See lines 1010 - 1209  qrrain
      f1pr09(0), // minimum ice number (lambda limiter)  See lines  704 -  705  nlarge
      f1pr10(0), // maximum ice number (lambda limiter)  See lines  704 -  705  nsmall
      f1pr13(0), // reflectivity                         See lines  731 -  808  refl
      f1pr14(0), // melting (ventilation term)           See lines 1212 - 1279  vdep1
      f1pr15(0), // mass-weighted mean diameter          See lines 1212 - 1279  dmm
      f1pr16(0), // mass-weighted mean particle density  See lines 1212 - 1279  rhomm

      // TODO(doc)
      vtrmi1   (0),   // TODO(doc)
      rhorime_c(400), // TODO(doc)
      epsi(0),        // TODO(doc)
      epsr(0),        // TODO(doc)
      epsc(0),        // TODO(doc)
      epsi_tot (0);   // inverse supersaturation relaxation timescale for combined ice categories

    Smask log_wetgrowth(false);

    // skip micro process calculations except nucleation/acvtivation if there no hydrometeors are present
    const auto skip_micro = skip_all || !(qc_incld(k) >= qsmall || qr_incld(k) >= qsmall || qitot_incld(k) >= qsmall);
    const auto not_skip_micro = !skip_micro;

    if (not_skip_micro.any()) {
      // time/space varying physical variables
      // TODO: needs smask protection (not_skip_micro)
      get_time_space_phys_variables(
        t(k), opres(k), rho(k), oxxlv(k), oxxls(k), qvs(k), qvi(k),
        mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii);

      get_cloud_dsd2(not_skip_micro, qc_incld(k), nc_incld(k), omu_c(k), rho(k), nu(k), dnu, olamc(k), cdist(k), cdist1(k), olcldm(k));

      get_rain_dsd2(not_skip_micro, qr_incld(k), nr_incld(k), mu_r(k), lamr(k), cdist(k), cdist1(k), orcldm(k));

      // TODO: needs smask protection
      impose_max_total_Ni(nitot_incld(k), max_total_Ni, inv_rho(k));

      const auto qitot_gt_small = qitot_incld(k) >= qsmall && not_skip_micro;

      if (qitot_gt_small.any()) {
        // impose lower limits to prevent taking log of # < 0
        nitot_incld(k).set(qitot_gt_small, pack::max(nitot_incld(k), nsmall));
        nr_incld(k).set(qitot_gt_small, pack::max(nr_incld(k), nsmall));

        const auto rhop = calc_bulk_rho_rime(qitot_gt_small, qitot_incld(k), qirim_incld(k), birim_incld(k));

        TableIce table_ice;
        lookup_ice(qitot_gt_small, qitot_incld(k), nitot_incld(k), qirim_incld(k), rhop, table_ice);

        TableRain table_rain;
        lookup_rain(qitot_gt_small, qr_incld(k), nr_incld(k), table_rain);

        // call to lookup table interpolation subroutines to get process rates
        f1pr02.set(qitot_gt_small, apply_table_ice(qitot_gt_small, 1,  itab, table_ice));
        f1pr03.set(qitot_gt_small, apply_table_ice(qitot_gt_small, 2,  itab, table_ice));
        f1pr04.set(qitot_gt_small, apply_table_ice(qitot_gt_small, 3,  itab, table_ice));
        f1pr05.set(qitot_gt_small, apply_table_ice(qitot_gt_small, 4,  itab, table_ice));
        f1pr09.set(qitot_gt_small, apply_table_ice(qitot_gt_small, 8,  itab, table_ice));
        f1pr10.set(qitot_gt_small, apply_table_ice(qitot_gt_small, 9,  itab, table_ice));
        f1pr14.set(qitot_gt_small, apply_table_ice(qitot_gt_small, 13, itab, table_ice));

        // ice-rain collection processes
        const auto qr_gt_small = qr_incld(k) >= qsmall && qitot_gt_small;
        f1pr07.set(qr_gt_small, apply_table_coll(qitot_gt_small, 0, itabcol, table_ice, table_rain));
        f1pr08.set(qr_gt_small, apply_table_coll(qitot_gt_small, 1, itabcol, table_ice, table_rain));

        // adjust Ni if needed to make sure mean size is in bounds (i.e. apply lambda limiters)
        // note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
        nitot_incld(k).set(qitot_gt_small, pack::min(nitot_incld(k), f1pr09*nitot_incld(k)));
        nitot_incld(k).set(qitot_gt_small, pack::max(nitot_incld(k), f1pr10*nitot_incld(k)));
      }

      // ----------------------------------------------------------------------
      // Begin calculations of microphysical processes
      // ----------------------------------------------------------------------

      // ......................................................................
      // ice processes
      // ......................................................................

      // TODO: needs smask protection (not_skip_micro)
      // collection of droplets
      ice_cldliq_collection(
        rho(k), t(k), rhofaci(k), f1pr04, qitot_incld(k), qc_incld(k), nitot_incld(k), nc_incld(k),
        qccol, nccol, qcshd, ncshdc);

      // TODO: needs smask protection
      // collection of rain
      ice_rain_collection(
        rho(k), t(k), rhofaci(k), logn0r(k), f1pr07, f1pr08, qitot_incld(k), nitot_incld(k), qr_incld(k),
        qrcol, nrcol);

      // collection between ice categories

      // PMC nCat deleted lots of stuff here.

      // TODO: needs smask protection
      // self-collection of ice
      ice_self_collection(
        rho(k), rhofaci(k), f1pr03, eii, qirim_incld(k), qitot_incld(k), nitot_incld(k),
        nislf);

      // TODO: needs smask protection
      // melting
      ice_melting(
        rho(k), t(k), opres(k), rhofaci(k), f1pr05, f1pr14, oxxlv(k), oxlf(k), dv, sc, mu, kap, oqv(k), qitot_incld(k), nitot_incld(k),
        qimlt, nimlt);

      // TODO: needs smask protection
      // calculate wet growth
      ice_cldliq_wet_growth(
        rho(k), t(k), opres(k), rhofaci(k), f1pr05, f1pr14, oxxlv(k), oxlf(k), dv, kap, mu, sc, oqv(k), qc_incld(k), qitot_incld(k), nitot_incld(k), qr_incld(k),
        log_wetgrowth, qrcol, qccol, qwgrth, nrshdr, qcshd);

      // TODO: needs smask protection
      // calcualte total inverse ice relaxation timescale combined for all ice categories
      // note 'f1pr' values are normalized, so we need to multiply by N
      ice_relaxation_timescale(
        rho(k), t(k), rhofaci(k), f1pr05, f1pr14, dv, mu, sc, qitot_incld(k), nitot_incld(k),
        epsi, epsi_tot);

      // TODO: needs smask protection
      // calculate rime density
      calc_rime_density(
        t(k), rhofaci(k), f1pr02, acn(k), olamc(k), omu_c(k), qc_incld(k), qccol,
        vtrmi1, rhorime_c);

      // TODO: needs smask protection
      // contact and immersion freezing droplets
      cldliq_immersion_freezing(
        t(k), olamc(k), omu_c(k), cdist1(k), qc_incld(k),
        qcheti, ncheti);

      // TODO: needs smask protection
      // for future: get rid of log statements below for rain freezing
      rain_immersion_freezing(
        t(k), lamr(k), mu_r(k), cdistr(k), qr_incld(k),
        qrheti, nrheti);

      //  rime splintering (Hallet-Mossop 1974)
      // PMC comment: Morrison and Milbrandt 2015 part 1 and 2016 part 3 both say
      // that Hallet-Mossop should be neglected if 1 category to compensate for
      // artificial smearing out of ice DSD

      // ................................................
      //  condensation/evaporation/deposition/sublimation
      //    (use semi-analytic formulation)

      // TODO: needs smask protection
      //  calculate rain evaporation including ventilation
      calc_liq_relaxation_timescale(
        revap_table, rho(k), f1r, f2r, dv, mu, sc, mu_r(k), lamr(k), cdistr(k), cdist(k), qr_incld(k), qc_incld(k),
        epsr, epsc);

      // TODO: needs smask protection
      evaporate_sublimate_precip(
        qr_incld(k), qc_incld(k), nr_incld(k), qitot_incld(k), olcldm(k), orcldm(k), qvs(k), ab, epsr, oqv(k),
        qrevp, nrevp);

      // TODO: needs smask protection
      ice_deposition_sublimation(
        qitot_incld(k), nitot_incld(k), t(k), qvs(k), qvi(k), epsi, abi, oqv(k),
        qidep, qisub, nisub, qiberg);
      }

    // TODO: needs smask protection (not_skip_all)
    // deposition/condensation-freezing nucleation
    ice_nucleation(
      t(k), inv_rho(k), onitot(k), onaai(k), supi(k), odt, log_predictNc,
      qinuc, ninuc);

    // TODO: needs smask protection
    // droplet activation
    droplet_activation(
      t(k), opres(k), oqv(k), oqc(k), inv_rho(k), sup(k), oxxlv(k), onpccn(k), log_predictNc, odt,
      qcnuc, ncnuc);

    // TODO: needs smask protection
    // cloud water autoconversion
    // NOTE: cloud_water_autoconversion must be called before droplet_self_collection
    cloud_water_autoconversion(
      rho(k), qc_incld(k), nc_incld(k),
      qcaut, ncautc, ncautr);

    // TODO: needs smask protection
    // self-collection of droplets
    droplet_self_collection(
      rho(k), inv_rho(k), qc_incld(k),
      omu_c(k), nu(k), ncautc, ncslf);

    // TODO: needs smask protection
    // accretion of cloud by rain
    cloud_rain_accretion(
      rho(k), inv_rho(k), qc_incld(k), nc_incld(k), qr_incld(k),
      qcacc, ncacc);

    // TODO: needs smask protection
    // self-collection and breakup of rain
    // (breakup following modified Verlinde and Cotton scheme)
    rain_self_collection(
      rho(k), qr_incld(k), nr_incld(k),
      nrslf);

    // TODO: needs smask protection
    // Here we map the microphysics tendency rates back to CELL-AVERAGE quantities for updating
    // cell-average quantities.
    back_to_cell_average(
      olcldm(k), orcldm(k), oicldm(k), qcacc, qrevp, qcaut,
      ncacc, ncslf, ncautc, nrslf, nrevp, ncautr, qcnuc, ncnuc, qisub, nrshdr, qcheti,
      qrcol, qcshd, qimlt, qccol, qrheti, nimlt, nccol, ncshdc, ncheti, nrcol, nislf,
      qidep, nrheti, nisub, qinuc, ninuc, qiberg);

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
    // how qidep and qisub are derived?
    // AaronDonahue: UPDATE, if we are following the implementation of MG
    // then the answer appears to be YES.  There is a similar check in MG
    // microphysics which limits qidep and qinuc, but does not limit qisub.
    // So similar but slightly different.  The overall answer would be that
    // qidep does need some limit.  The next questions are,
    //   1) Should we be taking qinuc into consideration too?
    //   2) Is MG correct in NOT limiting qisub?

    // TODO: needs smask protection
    prevent_ice_overdepletion(opres(k), t(k), oqv(k), oxxls(k), odt, qidep, qisub);

    // vapor -- not needed, since all sinks already have limits imposed and the sum, therefore,
    //          cannot possibly overdeplete qv

    // cloud
    // TODO: needs smask protection
    cloud_water_conservation(
      oqc(k), qcnuc, dt,
      qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep);

    // rain
    // TODO: needs smask protection
    rain_water_conservation(
      oqr(k), qcaut, qcacc, qimlt, qcshd, dt,
      qrevp, qrcol, qrheti);

    // ice
    // TODO: needs smask protection
    ice_water_conservation(
      oqitot(k), qidep, qinuc, qiberg, qrcol, qccol, qrheti, qcheti, dt,
      qisub, qimlt);

    //---------------------------------------------------------------------------------
    // update prognostic microphysics and thermodynamics variables
    //---------------------------------------------------------------------------------

    // TODO: needs smask protection
    //-- ice-phase dependent processes:
    update_prognostic_ice(
      qcheti, qccol, qcshd, nccol, ncheti, ncshdc, qrcol, nrcol,  qrheti, nrheti, nrshdr, qimlt, nimlt, qisub, qidep, qinuc, ninuc, nislf, nisub, qiberg, oexner(k), oxxls(k), oxlf(k), log_predictNc, log_wetgrowth, dt, nmltratio, rhorime_c,
      oth(k), oqv(k), oqitot(k), onitot(k), oqirim(k), obirim(k), oqc(k), onc(k), oqr(k), onr(k));

    // TODO: needs smask protection
    //-- warm-phase only processes:
    update_prognostic_liquid(
      qcacc, ncacc, qcaut, ncautc, qcnuc, ncautr, ncslf, qrevp, nrevp, nrslf, log_predictNc, inv_rho(k), oexner(k), oxxlv(k), dt,
      oth(k), oqv(k), oqc(k), onc(k), oqr(k), onr(k));

    // AaronDonahue - Add extra variables needed from microphysics by E3SM:
    ocmeiout(k)         .set(not_skip_all, qidep - qisub + qinuc);
    oprain(k)           .set(not_skip_all, qcacc + qcaut + qcshd + qccol);
    onevapr(k)          .set(not_skip_all, qisub + qrevp);
    oprer_evap(k)       .set(not_skip_all, qrevp);
    ovap_ice_exchange(k).set(not_skip_all, qidep - qisub + qinuc);
    ovap_liq_exchange(k).set(not_skip_all, -qrevp + qcnuc);
    oliq_ice_exchange(k).set(not_skip_all, qcheti + qrheti - qimlt + qiberg + qccol + qrcol);
    ovap_cld_exchange(k).set(not_skip_all, qcnuc);

    // clipping for small hydrometeor values
    const auto qc_small    = oqc(k) < qsmall    && not_skip_all;
    const auto qr_small    = oqr(k) < qsmall    && not_skip_all;
    const auto qitot_small = oqitot(k) < qsmall && not_skip_all;

    const auto qc_not_small    = oqc(k) >= qsmall    && not_skip_all;
    const auto qr_not_small    = oqr(k) >= qsmall    && not_skip_all;
    const auto qitot_not_small = oqitot(k) >= qsmall && not_skip_all;

    oqv(k).set(qc_small, oqv(k) + oqc(k));
    oth(k).set(qc_small, oth(k) - oexner(k) * oqc(k) * oxxlv(k) * inv_cp);
    oqc(k).set(qc_small, 0);
    onc(k).set(qc_small, 0);

    if (qc_not_small.any()) {
      log_hydrometeorsPresent = true;
    }

    oqv(k).set(qr_small, oqv(k) + oqr(k));
    oth(k).set(qr_small, oth(k) - oexner(k) * oqr(k) * oxxlv(k) * inv_cp);
    oqr(k).set(qr_small, 0);
    onr(k).set(qr_small, 0);

    if (qr_not_small.any()) {
      log_hydrometeorsPresent = true;
    }

    oqv(k).set(qitot_small, oqv(k) + oqitot(k));
    oth(k).set(qitot_small, oth(k) - oexner(k) * oqitot(k) * oxxls(k) * inv_cp);
    oqitot(k).set(qitot_small, 0);
    onitot(k).set(qitot_small, 0);
    oqirim(k).set(qitot_small, 0);
    obirim(k).set(qitot_small, 0);

    if (qitot_not_small.any()) {
      log_hydrometeorsPresent = true;
    }

    // TODO: needs smask protection
    impose_max_total_Ni(onitot(k), max_total_Ni, inv_rho(k));

    // Outputs associated with aerocom comparison:
    opratot(k) = qcacc; // cloud drop accretion by rain
    oprctot(k) = qcaut; // cloud drop autoconversion to rain

    // TODO: needs smask protection
    // Recalculate in-cloud values for sedimentation
    calculate_incloud_mixingratios(
      oqc(k), oqr(k), oqitot(k), oqirim(k), onc(k), onr(k), onitot(k), obirim(k), inv_lcldm(k), inv_icldm(k), inv_rcldm(k),
      qc_incld(k), qr_incld(k), qitot_incld(k), qirim_incld(k), nc_incld(k), nr_incld(k), nitot_incld(k), birim_incld(k));
    });
  team.team_barrier();
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main_post_main_loop(
  const MemberType& team,
  const Int& nk_pack,
  const bool& log_predictNc,
  const Scalar& dt,
  const Scalar& odt,
  const view_dnu_table& dnu,
  const view_itab_table& itab,
  const view_itabcol_table& itabcol,
  const view_2d_table& revap_table,
  const uview_1d<const Spack>& opres,
  const uview_1d<const Spack>& opdel,
  const uview_1d<const Spack>& odzq,
  const uview_1d<const Spack>& onpccn,
  const uview_1d<const Spack>& oexner,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<const Spack>& inv_lcldm,
  const uview_1d<const Spack>& inv_icldm,
  const uview_1d<const Spack>& inv_rcldm,
  const uview_1d<const Spack>& onaai,
  const uview_1d<const Spack>& oicldm,
  const uview_1d<const Spack>& olcldm,
  const uview_1d<const Spack>& orcldm,
  const uview_1d<Spack>& t,
  const uview_1d<Spack>& rho,
  const uview_1d<Spack>& inv_rho,
  const uview_1d<Spack>& qvs,
  const uview_1d<Spack>& qvi,
  const uview_1d<Spack>& sup,
  const uview_1d<Spack>& supi,
  const uview_1d<Spack>& rhofacr,
  const uview_1d<Spack>& rhofaci,
  const uview_1d<Spack>& acn,
  const uview_1d<Spack>& oqv,
  const uview_1d<Spack>& oth,
  const uview_1d<Spack>& oqc,
  const uview_1d<Spack>& onc,
  const uview_1d<Spack>& oqr,
  const uview_1d<Spack>& onr,
  const uview_1d<Spack>& oqitot,
  const uview_1d<Spack>& onitot,
  const uview_1d<Spack>& oqirim,
  const uview_1d<Spack>& obirim,
  const uview_1d<Spack>& oxxlv,
  const uview_1d<Spack>& oxxls,
  const uview_1d<Spack>& oxlf,
  const uview_1d<Spack>& qc_incld,
  const uview_1d<Spack>& qr_incld,
  const uview_1d<Spack>& qitot_incld,
  const uview_1d<Spack>& qirim_incld,
  const uview_1d<Spack>& nc_incld,
  const uview_1d<Spack>& nr_incld,
  const uview_1d<Spack>& nitot_incld,
  const uview_1d<Spack>& birim_incld,
  const uview_1d<Spack>& omu_c,
  const uview_1d<Spack>& nu,
  const uview_1d<Spack>& olamc,
  const uview_1d<Spack>& cdist,
  const uview_1d<Spack>& cdist1,
  const uview_1d<Spack>& cdistr,
  const uview_1d<Spack>& mu_r,
  const uview_1d<Spack>& lamr,
  const uview_1d<Spack>& logn0r,
  const uview_1d<Spack>& ocmeiout,
  const uview_1d<Spack>& oprain,
  const uview_1d<Spack>& onevapr,
  const uview_1d<Spack>& oprer_evap,
  const uview_1d<Spack>& ovap_cld_exchange,
  const uview_1d<Spack>& ovap_liq_exchange,
  const uview_1d<Spack>& ovap_ice_exchange,
  const uview_1d<Spack>& oliq_ice_exchange,
  const uview_1d<Spack>& opratot,
  const uview_1d<Spack>& oprctot,
  const uview_1d<Spack>& ze_rain,
  const uview_1d<Spack>& ze_ice,
  const uview_1d<Spack>& odiag_vmi,
  const uview_1d<Spack>& odiag_effi,
  const uview_1d<Spack>& odiag_di,
  const uview_1d<Spack>& odiag_rhoi,
  const uview_1d<Spack>& odiag_ze,
  const uview_1d<Spack>& tmparr1)
{
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar inv_cp       = C::INV_CP;
  constexpr Scalar max_total_Ni = C::max_total_Ni;
  constexpr Scalar nsmall       = C::NSMALL;

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {

    Spack
      ignore1  (0),
      ignore2  (0),
      f1pr02   (0),
      f1pr06   (0),
      f1pr09   (0),
      f1pr10   (0),
      f1pr13   (0),
      f1pr15   (0),
      f1pr16   (0);

    // Cloud
    {
      const auto qc_gt_small = oqc(k) >= qsmall;
      const auto qc_small    = !qc_gt_small;
      get_cloud_dsd2(
        qc_gt_small, oqc(k), onc(k), omu_c(k), rho(k), nu(k), dnu, olamc(k), ignore1, ignore2, olcldm(k));

      oqv(k)              .set(qc_small, oqv(k)+oqc(k));
      oth(k)              .set(qc_small, oth(k)-oexner(k)*oqc(k)*oxxlv(k)*inv_cp);
      ovap_liq_exchange(k).set(qc_small, ovap_liq_exchange(k) - oqc(k));
      ovap_cld_exchange(k).set(qc_small, ovap_cld_exchange(k) - oqc(k));
      oqc(k)              .set(qc_small, 0);
      onc(k)              .set(qc_small, 0);
    }

    // Rain
    {
      const auto qr_gt_small = oqr(k) >= qsmall;
      const auto qr_small    = !qr_gt_small;
      get_rain_dsd2(
        qr_gt_small, oqr(k), onr(k), mu_r(k), lamr(k), ignore1, ignore2, orcldm(k));

      ze_rain(k).set(qr_gt_small, onr(k)*(mu_r(k)+6)*(mu_r(k)+5)*(mu_r(k)+4)*
                     (mu_r(k)+3)*(mu_r(k)+2)*(mu_r(k)+1)/pow(lamr(k), 6));
      ze_rain(k).set(qr_gt_small, pack::max(ze_rain(k), 1.e-22));

      oqv(k)              .set(qr_small, oqv(k) + oqr(k));
      oth(k)              .set(qr_small, oth(k) - oexner(k)*oqr(k)*oxxlv(k)*inv_cp);
      ovap_liq_exchange(k).set(qr_small, ovap_liq_exchange(k) - oqr(k));
      oqr(k)              .set(qr_small, 0);
      onr(k)              .set(qr_small, 0);
    }

    // Ice
    {
      impose_max_total_Ni(onitot(k), max_total_Ni, inv_rho(k));

      const auto qi_gt_small = oqitot(k) >= qsmall;
      const auto qi_small    = !qi_gt_small;

      // impose lower limits to prevent taking log of # < 0
      onitot(k) = pack::max(onitot(k), nsmall);
      onr(k)    = pack::max(onr(k), nsmall);

      const auto rhop = calc_bulk_rho_rime(qi_gt_small, oqitot(k), oqirim(k), obirim(k));

      TableIce table_ice;
      lookup_ice(qi_gt_small, oqitot(k), onitot(k), oqirim(k), rhop, table_ice);

      f1pr02.set(qi_gt_small, apply_table_ice(qi_gt_small, 1,  itab, table_ice));
      f1pr06.set(qi_gt_small, apply_table_ice(qi_gt_small, 5,  itab, table_ice));
      f1pr09.set(qi_gt_small, apply_table_ice(qi_gt_small, 6,  itab, table_ice));
      f1pr10.set(qi_gt_small, apply_table_ice(qi_gt_small, 7,  itab, table_ice));
      f1pr13.set(qi_gt_small, apply_table_ice(qi_gt_small, 8,  itab, table_ice));
      f1pr15.set(qi_gt_small, apply_table_ice(qi_gt_small, 10, itab, table_ice));
      f1pr16.set(qi_gt_small, apply_table_ice(qi_gt_small, 11, itab, table_ice));

      // impose mean ice size bounds (i.e. apply lambda limiters)
      // note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
      onitot(k).set(qi_gt_small, pack::min(onitot(k), f1pr09 * onitot(k)));
      onitot(k).set(qi_gt_small, pack::max(onitot(k), f1pr10 * onitot(k)));

      // --this should already be done in s/r 'calc_bulkRhoRime'
      const auto qirim_small = oqirim(k) < qsmall && qi_gt_small;
      oqirim(k).set(qirim_small, 0);
      obirim(k).set(qirim_small, 0);

      // note that reflectivity from lookup table is normalized, so we need to multiply by N
      odiag_vmi(k) .set(qi_gt_small, f1pr02 * rhofaci(k));
      odiag_effi(k).set(qi_gt_small, f1pr06); // units are in m
      odiag_di(k)  .set(qi_gt_small, f1pr15);
      odiag_rhoi(k).set(qi_gt_small, f1pr16);

      // note factor of air density below is to convert from m^6/kg to m^6/m^3
      ze_ice(k).set(qi_gt_small, ze_ice(k) + sp(0.1892)*f1pr13*onitot(k)*rho(k)); // sum contribution from each ice category (note: 0.1892 = 0.176/0.93);
      ze_ice(k).set(qi_gt_small, pack::max(ze_ice(k), 1.e-22));

      oqv(k)     .set(qi_small, oqv(k) + oqitot(k));
      oth(k)     .set(qi_small, oth(k) - oexner(k)*oqitot(k)*oxxls(k)*inv_cp);
      oqitot(k)  .set(qi_small, 0);
      onitot(k)  .set(qi_small, 0);
      oqirim(k)  .set(qi_small, 0);
      obirim(k)  .set(qi_small, 0);
      odiag_di(k).set(qi_small, 0);
    }

    // sum ze components and convert to dBZ
    odiag_ze(k) = 10 * pack::log10((ze_rain(k) + ze_ice(k))*sp(1.e18));

    // if qr is very small then set Nr to 0 (needs to be done here after call
    // to ice lookup table because a minimum Nr of nsmall will be set otherwise even if qr=0)
    onr(k).set(oqr(k) < qsmall, 0);

#ifndef NDEBUG
    tmparr1(k) = oth(k) * inv_exner(k);
#endif
  });
  team.team_barrier();
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main(
  // inputs
  const view_2d<const Spack>& pres,          // pressure                             Pa
  const view_2d<const Spack>& dzq,           // vertical grid spacing                m
  const view_2d<const Spack>& npccn,         // IN ccn activated number tendency     kg-1 s-1
  const view_2d<const Spack>& naai,          // IN actived ice nuclei concentration  1/kg
  const Real&                 dt,            // model time step                      s
  const Int&                  ni,            // num columns
  const Int&                  nk,            // column size
  const Int&                  it,            // time step counter NOTE: starts at 1 for first time step
  const bool&                 log_predictNc, // .T. (.F.) for prediction (specification) of Nc
  const view_2d<const Spack>& pdel,          // pressure thickness                   Pa
  const view_2d<const Spack>& exner,         // Exner expression

  // inputs needed for PBUF variables used by other parameterizations
  const view_2d<const Spack>& icldm,         // ice cloud fraction
  const view_2d<const Spack>& lcldm,         // liquid cloud fraction
  const view_2d<const Spack>& rcldm,         // rain cloud fraction
  const view_2d<const Scalar>& col_location, // ni x 3

  // input/output  arguments
  const view_2d<Spack>& qc,    // cloud, mass mixing ratio         kg kg-1
  const view_2d<Spack>& nc,    // cloud, number mixing ratio       #  kg-1
  const view_2d<Spack>& qr,    // rain, mass mixing ratio          kg kg-1
  const view_2d<Spack>& nr,    // rain, number mixing ratio        #  kg-1
  const view_2d<Spack>& qitot, // ice, total mass mixing ratio     kg kg-1
  const view_2d<Spack>& qirim, // ice, rime mass mixing ratio      kg kg-1
  const view_2d<Spack>& nitot, // ice, total number mixing ratio   #  kg-1
  const view_2d<Spack>& birim, // ice, rime volume mixing ratio    m3 kg-1
  const view_2d<Spack>& qv,    // water vapor mixing ratio         kg kg-1
  const view_2d<Spack>& th,    // potential temperature            K

  // output arguments
  const view_1d<Scalar>& prt_liq,  // precipitation rate, liquid       m s-1
  const view_1d<Scalar>& prt_sol,  // precipitation rate, solid        m s-1
  const view_2d<Spack>& diag_ze,   // equivalent reflectivity          dBZ
  const view_2d<Spack>& diag_effc, // effective radius, cloud          m
  const view_2d<Spack>& diag_effi, // effective radius, ice            m
  const view_2d<Spack>& diag_vmi,  // mass-weighted fall speed of ice  m s-1
  const view_2d<Spack>& diag_di,   // mean diameter of ice             m
  const view_2d<Spack>& diag_rhoi, // bulk density of ice              kg m-3
  const view_2d<Spack>& mu_c,      // Size distribution shape parameter for radiation
  const view_2d<Spack>& lamc,      // Size distribution slope parameter for radiation

  // outputs for PBUF variables used by other parameterizations
  const view_2d<Spack>& cmeiout,          // qitend due to deposition/sublimation
  const view_2d<Spack>& prain,            // Total precipitation (rain + snow)
  const view_2d<Spack>& nevapr,           // evaporation of total precipitation (rain + snow)
  const view_2d<Spack>& prer_evap,        // evaporation of rain
  const view_2d<Spack>& rflx,             // grid-box average rain flux (kg m^-2 s^-1) pverp
  const view_2d<Spack>& sflx,             // grid-box average ice/snow flux (kg m^-2 s^-1) pverp
  const view_2d<Spack>& pratot,           // accretion of cloud by rain
  const view_2d<Spack>& prctot,           // autoconversion of cloud to rain
  const view_2d<Spack>& liq_ice_exchange, // sum of liq-ice phase change tendenices
  const view_2d<Spack>& vap_liq_exchange, // sum of vap-liq phase change tendenices
  const view_2d<Spack>& vap_ice_exchange, // sum of vap-ice phase change tendenices
  const view_2d<Spack>& vap_cld_exchange) // sum of vap-cld phase change tendenices
{
  using ExeSpace = typename KT::ExeSpace;

  view_2d<Spack> xxls("xxls", ni, nk), xxlv("xxlv", ni, nk), xlf("xlf", ni, nk);

  get_latent_heat(ni, nk, xxlv, xxls, xlf);

  const Int nk_pack = scream::pack::npack<Spack>(nk);
  const auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ni, nk_pack);

  WorkspaceManager<Spack, Device> workspace_mgr(nk_pack, 100, policy);

  // load constants into local vars
  constexpr Scalar qsmall       = C::QSMALL;
  constexpr Scalar nsmall       = C::NSMALL;
  constexpr Scalar inv_cp       = C::INV_CP;
  constexpr Scalar max_total_Ni = C::max_total_Ni;
  const     Scalar odt          = 1 / dt;
  constexpr Int    kdir         = -1;
  const     Int    ktop         = kdir == -1 ? 0    : nk-1;
  const     Int    kbot         = kdir == -1 ? nk-1 : 0;
  constexpr bool   debug_ABORT  = false;

  // load tables
  view_1d_table mu_r_table;
  view_2d_table vn_table, vm_table, revap_table;
  view_itab_table itab;
  view_itabcol_table itabcol;
  view_dnu_table dnu;
  init_kokkos_ice_lookup_tables(itab, itabcol);
  init_kokkos_tables(vn_table, vm_table, revap_table, mu_r_table, dnu);

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
      t,      // temperature at the beginning of the microhpysics step [K]

      // 2D size distribution and fallspeed parameters
      lamr, logn0r, nu, cdist, cdist1, cdistr,

      // Variables needed for in-cloud calculations
      inv_icldm, inv_lcldm, inv_rcldm, // Inverse cloud fractions (1/cld)
      qc_incld, qr_incld, qitot_incld, qirim_incld, // In cloud mass-mixing ratios
      nc_incld, nr_incld, nitot_incld, birim_incld, // In cloud number concentrations

      // Other
      inv_dzq, inv_rho, ze_ice, ze_rain, prec, rho,
      rhofacr, rhofaci, acn, qvs, qvi, sup, supi,
      tmparr1, inv_exner,

      // p3_tend_out, may not need these
      qc_tend, nc_tend, qr_tend, nr_tend, qi_tend, ni_tend;

    workspace.template take_many_and_reset<40>(
      {
        "mu_r", "t", "lamr", "logn0r", "nu", "cdist", "cdist1", "cdistr",
        "inv_icldm", "inv_lcldm", "inv_rcldm", "qc_incld", "qr_incld", "qitot_incld", "qirim_incld",
        "nc_incld", "nr_incld", "nitot_incld", "birim_incld",
        "inv_dzq", "inv_rho", "ze_ice", "ze_rain", "prec", "rho",
        "rhofacr", "rhofaci", "acn", "qvs", "qvi", "sup", "supi",
        "tmparr1", "inv_exner", "qc_tend", "nc_tend", "qr_tend", "nr_tend", "qi_tend", "ni_tend"
      },
      {
        &mu_r, &t, &lamr, &logn0r, &nu, &cdist, &cdist1, &cdistr,
        &inv_icldm, &inv_lcldm, &inv_rcldm, &qc_incld, &qr_incld, &qitot_incld, &qirim_incld,
        &nc_incld, &nr_incld, &nitot_incld, &birim_incld,
        &inv_dzq, &inv_rho, &ze_ice, &ze_rain, &prec, &rho,
        &rhofacr, &rhofaci, &acn, &qvs, &qvi, &sup, &supi,
        &tmparr1, &inv_exner, &qc_tend, &nc_tend, &qr_tend, &nr_tend, &qi_tend, &ni_tend
      });

    // Get single-column subviews of all inputs, shouldn't need any i-indexing
    // after this.
    const auto opres             = util::subview(pres, i);
    const auto odzq              = util::subview(dzq, i);
    const auto onpccn            = util::subview(npccn, i);
    const auto onaai             = util::subview(naai, i);
    const auto opdel             = util::subview(pdel, i);
    const auto oexner            = util::subview(exner, i);
    const auto oicldm            = util::subview(icldm, i);
    const auto olcldm            = util::subview(lcldm, i);
    const auto orcldm            = util::subview(rcldm, i);
    const auto ocol_location     = util::subview(col_location, i);
    const auto oqc               = util::subview(qc, i);
    const auto onc               = util::subview(nc, i);
    const auto oqr               = util::subview(qr, i);
    const auto onr               = util::subview(nr, i);
    const auto oqitot            = util::subview(qitot, i);
    const auto oqirim            = util::subview(qirim, i);
    const auto onitot            = util::subview(nitot, i);
    const auto obirim            = util::subview(birim, i);
    const auto oqv               = util::subview(qv, i);
    const auto oth               = util::subview(th, i);
    const auto odiag_ze          = util::subview(diag_ze, i);
    const auto odiag_effc        = util::subview(diag_effc, i);
    const auto odiag_effi        = util::subview(diag_effi, i);
    const auto odiag_vmi         = util::subview(diag_vmi, i);
    const auto odiag_di          = util::subview(diag_di, i);
    const auto odiag_rhoi        = util::subview(diag_rhoi, i);
    const auto omu_c             = util::subview(mu_c, i);
    const auto olamc             = util::subview(lamc, i);
    const auto ocmeiout          = util::subview(cmeiout, i);
    const auto oprain            = util::subview(prain, i);
    const auto onevapr           = util::subview(nevapr, i);
    const auto oprer_evap        = util::subview(prer_evap, i);
    const auto orflx             = util::subview(rflx, i);
    const auto osflx             = util::subview(sflx, i);
    const auto opratot           = util::subview(pratot, i);
    const auto oprctot           = util::subview(prctot, i);
    const auto oliq_ice_exchange = util::subview(liq_ice_exchange, i);
    const auto ovap_liq_exchange = util::subview(vap_liq_exchange, i);
    const auto ovap_ice_exchange = util::subview(vap_ice_exchange, i);
    const auto ovap_cld_exchange = util::subview(vap_cld_exchange, i);
    const auto oxxlv             = util::subview(xxlv, i);
    const auto oxxls             = util::subview(xxls, i);
    const auto oxlf              = util::subview(xlf, i);

    // Need to watch out for race conditions with these shared variables
    bool log_hydrometeorsPresent = false;
    bool log_nucleationPossible  = false;

    // initialize
    p3_main_init(
      team, nk_pack,
      oicldm, olcldm, orcldm, oexner, oth,
      opratot, oprctot, prec, mu_r, odiag_ze, ze_ice, ze_rain, odiag_effc, odiag_effi, odiag_vmi, odiag_di, odiag_rhoi, ocmeiout, oprain, onevapr, orflx, osflx, inv_icldm, inv_lcldm, inv_rcldm, omu_c, olamc, inv_exner, t, oqv,
      prt_liq(i), prt_sol(i));

    p3_main_pre_main_loop(
      team, nk, log_predictNc, dt,
      opres, opdel, odzq, onpccn, oexner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, oxxlv, oxxls, oxlf,
      t, rho, inv_rho, qvs, qvi, sup, supi, rhofacr, rhofaci, acn, oqv, oth, oqc, onc, oqr, onr, oqitot, onitot, oqirim, obirim, qc_incld, qr_incld, qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld,
      log_nucleationPossible, log_hydrometeorsPresent);

    // There might not be any work to do for this team
    if (!(log_nucleationPossible || log_hydrometeorsPresent)) {
      return; // this is how you do a "continue" in a kokkos lambda
    }

    team.team_barrier();
    log_hydrometeorsPresent = false; // reset value; used again below
    team.team_barrier();

    // ------------------------------------------------------------------------------------------
    // main k-loop (for processes):
    p3_main_main_loop(
      team, nk_pack, log_predictNc, dt, odt,
      dnu, itab, itabcol, revap_table,
      opres, opdel, odzq, onpccn, oexner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, onaai, oicldm, olcldm, orcldm,
      t, rho, inv_rho, qvs, qvi, sup, supi, rhofacr, rhofaci, acn, oqv, oth, oqc, onc, oqr, onr, oqitot, onitot, oqirim, obirim, oxxlv, oxxls, oxlf, qc_incld, qr_incld, qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld, omu_c, nu, olamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, ocmeiout, oprain, onevapr, oprer_evap, ovap_cld_exchange, ovap_liq_exchange, ovap_ice_exchange, oliq_ice_exchange, opratot, oprctot,
      log_hydrometeorsPresent);

    //NOTE: At this point, it is possible to have negative (but small) nc, nr, nitot.  This is not
    //      a problem; those values get clipped to zero in the sedimentation section (if necessary).
    //      (This is not done above simply for efficiency purposes.)

    if (!log_hydrometeorsPresent) return;

    // -----------------------------------------------------------------------------------------
    // End of main microphysical processes section
    // =========================================================================================

    // ==========================================================================================!
    // Sedimentation:

    // Cloud sedimentation:  (adaptive substepping)
    Kokkos::deep_copy(qc_tend, oqc);
    Kokkos::deep_copy(nc_tend, onc);
    cloud_sedimentation(
      qc_incld, rho, inv_rho, olcldm, acn, inv_dzq, dnu, team, workspace,
      nk, ktop, kbot, kdir, dt, odt, log_predictNc,
      oqc, onc, nc_incld, omu_c, olamc, qc_tend, nc_tend, prt_liq(i));

    // Rain sedimentation:  (adaptive substepping)
    Kokkos::deep_copy(qr_tend, oqr);
    Kokkos::deep_copy(nr_tend, onr);
    rain_sedimentation(
      rho, inv_rho, rhofacr, orcldm, inv_dzq, qr_incld, team, workspace, vn_table, vm_table,
      nk, ktop, kbot, kdir, dt, odt,
      oqr, onr, nr_incld, mu_r, lamr, orflx, qr_tend, nc_tend, prt_liq(i));

    // Ice sedimentation:  (adaptive substepping)
    Kokkos::deep_copy(qi_tend, oqitot);
    Kokkos::deep_copy(ni_tend, onitot);
    ice_sedimentation(
      rho, inv_rho, rhofaci, oicldm, inv_dzq, team, workspace,
      nk, ktop, kbot, kdir, dt, odt,
      oqitot, qitot_incld, onitot, nitot_incld, oqirim, qirim_incld, obirim, birim_incld, qi_tend, ni_tend, itab, prt_sol(i));

    // homogeneous freezing of cloud and rain
    homogeneous_freezing(
      t, oexner, oxlf, team, nk, ktop, kbot, kdir,
      oqc, onc, oqr, onr, oqitot, onitot, oqirim, obirim, oth);

    //
    // final checks to ensure consistency of mass/number
    // and compute diagnostic fields for output
    //
    p3_main_post_main_loop(
      team, nk_pack, log_predictNc, dt, odt,
      dnu, itab, itabcol, revap_table,
      opres, opdel, odzq, onpccn, oexner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, onaai, oicldm, olcldm, orcldm,
      t, rho, inv_rho, qvs, qvi, sup, supi, rhofacr, rhofaci, acn, oqv, oth, oqc, onc, oqr, onr, oqitot, onitot, oqirim, obirim, oxxlv, oxxls, oxlf, qc_incld, qr_incld, qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld, omu_c, nu, olamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, ocmeiout, oprain, onevapr, oprer_evap, ovap_cld_exchange, ovap_liq_exchange, ovap_ice_exchange, oliq_ice_exchange, opratot, oprctot, ze_rain, ze_ice, odiag_vmi, odiag_effi, odiag_di, odiag_rhoi, odiag_ze, tmparr1);

    //
    // merge ice categories with similar properties

    //   note:  this should be relocated to above, such that the diagnostic
    //          ice properties are computed after merging

    // PMC nCat deleted nCat>1 stuff

#ifndef NDEBUG
    check_values(oqv, tmparr1, ktop, kbot, it, debug_ABORT, 900, team, ocol_location);
#endif

  });
}

} // namespace p3
} // namespace scream

#endif
