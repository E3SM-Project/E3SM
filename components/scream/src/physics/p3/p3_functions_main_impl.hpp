#ifndef P3_FUNCTIONS_MAIN_IMPL_HPP
#define P3_FUNCTIONS_MAIN_IMPL_HPP

#include "p3_functions.hpp" // for ETI only but harmless for GPU
#include "physics_functions.hpp" // also for ETI not on GPUs
#include "physics_saturation_impl.hpp"

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
  const uview_1d<const Spack>& icldm,
  const uview_1d<const Spack>& lcldm,
  const uview_1d<const Spack>& rcldm,
  const uview_1d<const Spack>& exner,
  const uview_1d<const Spack>& th,
  const uview_1d<Spack>& pratot,
  const uview_1d<Spack>& prctot,
  const uview_1d<Spack>& prec,
  const uview_1d<Spack>& mu_r,
  const uview_1d<Spack>& diag_ze,
  const uview_1d<Spack>& ze_ice,
  const uview_1d<Spack>& ze_rain,
  const uview_1d<Spack>& diag_effc,
  const uview_1d<Spack>& diag_effi,
  const uview_1d<Spack>& diag_vmi,
  const uview_1d<Spack>& diag_di,
  const uview_1d<Spack>& diag_rhoi,
  const uview_1d<Spack>& cmeiout,
  const uview_1d<Spack>& prain,
  const uview_1d<Spack>& nevapr,
  const uview_1d<Spack>& rflx,
  const uview_1d<Spack>& sflx,
  const uview_1d<Spack>& inv_icldm,
  const uview_1d<Spack>& inv_lcldm,
  const uview_1d<Spack>& inv_rcldm,
  const uview_1d<Spack>& mu_c,
  const uview_1d<Spack>& lamc,
  const uview_1d<Spack>& inv_exner,
  const uview_1d<Spack>& t,
  const uview_1d<Spack>& qv,
  Scalar& prt_liq,
  Scalar& prt_sol)
{
  prt_liq = 0;
  prt_sol = 0;

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {

    pratot(k)    = 0;
    prctot(k)    = 0;
    prec(k)      = 0;
    mu_r(k)      = 0;
    diag_ze(k)   = -99;
    ze_ice(k)    = 1.e-22;
    ze_rain(k)   = 1.e-22;
    diag_effc(k) = 10.e-6;
    diag_effi(k) = 25.e-6;
    diag_vmi(k)  = 0;
    diag_di(k)   = 0;
    diag_rhoi(k) = 0;
    cmeiout(k)   = 0;
    prain(k)     = 0;
    nevapr(k)    = 0;
    rflx(k)      = 0;
    sflx(k)      = 0;
    inv_icldm(k) = 1 / icldm(k);
    inv_lcldm(k) = 1 / lcldm(k);
    inv_rcldm(k) = 1 / rcldm(k);
    mu_c(k)      = 0;
    lamc(k)      = 0;
    inv_exner(k) = 1 / exner(k);
    t(k)         = th(k) * inv_exner(k);
    qv(k)        = pack::max(qv(k), 0);
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
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& pdel,
  const uview_1d<const Spack>& dzq,
  const uview_1d<const Spack>& ncnuc,
  const uview_1d<const Spack>& exner,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<const Spack>& inv_lcldm,
  const uview_1d<const Spack>& inv_icldm,
  const uview_1d<const Spack>& inv_rcldm,
  const uview_1d<const Spack>& xxlv,
  const uview_1d<const Spack>& xxls,
  const uview_1d<const Spack>& xlf,
  const uview_1d<Spack>& t,
  const uview_1d<Spack>& rho,
  const uview_1d<Spack>& inv_rho,
  const uview_1d<Spack>& qvs,
  const uview_1d<Spack>& qvi,
  const uview_1d<Spack>& supi,
  const uview_1d<Spack>& rhofacr,
  const uview_1d<Spack>& rhofaci,
  const uview_1d<Spack>& acn,
  const uview_1d<Spack>& qv,
  const uview_1d<Spack>& th,
  const uview_1d<Spack>& qc,
  const uview_1d<Spack>& nc,
  const uview_1d<Spack>& qr,
  const uview_1d<Spack>& nr,
  const uview_1d<Spack>& qitot,
  const uview_1d<Spack>& nitot,
  const uview_1d<Spack>& qirim,
  const uview_1d<Spack>& birim,
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
  // Get access to saturation functions
  using physics = scream::physics::Functions<Scalar, Device>;

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

    rho(k)     = pdel(k)/dzq(k) / g;
    inv_rho(k) = 1 / rho(k);
    qvs(k)     = physics::qv_sat(t(k), pres(k), 0);
    qvi(k)     = physics::qv_sat(t(k), pres(k), 1);

    supi(k) = qv(k) / qvi(k) - 1;

    rhofacr(k) = pack::pow(rhosur * inv_rho(k), sp(.54));
    rhofaci(k) = pack::pow(rhosui * inv_rho(k), sp(.54));
    Spack dum  = sp(1.496e-6) * pack::pow(t(k), sp(1.5)) / (t(k) + 120); // this is mu
    acn(k)     = g * rhow / (18 * dum); // 'a' parameter for droplet fallspeed (Stokes' law)

    if ( (t(k) < zerodegc && supi(k) >= -0.05).any() ) {
      log_nucleationPossible = true;
    }

    // apply mass clipping if dry and mass is sufficiently small
    // (implying all mass is expected to evaporate/sublimate in one time step)
    auto drymass = qc(k) < qsmall;
    auto not_drymass = !drymass && range_mask;
    qv(k).set(drymass, qv(k) + qc(k));
    th(k).set(drymass, th(k) - exner(k) * qc(k) * xxlv(k) * inv_cp);
    qc(k).set(drymass, 0);
    nc(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      log_hydrometeorsPresent = true; // updated further down
      // Apply droplet activation here (before other microphysical processes) for consistency with qc increase by saturation
      // adjustment already applied in macrophysics. If prescribed drop number is used, this is also a good place to
      // prescribe that value
      if (!log_predictNc) {
         nc(k).set(not_drymass, nccnst*inv_rho(k));
      } else {
         nc(k).set(not_drymass, pack::max(nc(k) + ncnuc(k) * dt, 0.0));
      }
    }

    drymass = qr(k) < qsmall;
    not_drymass = !drymass && range_mask;
    qv(k).set(drymass, qv(k) + qr(k));
    th(k).set(drymass, th(k) - exner(k) * qr(k) * xxlv(k) * inv_cp);
    qr(k).set(drymass, 0);
    nr(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      log_hydrometeorsPresent = true; // updated further down
    }

    drymass = (qitot(k) < qsmall || (qitot(k) < 1.e-8 && supi(k) < -0.1));
    not_drymass = !drymass && range_mask;
    qv(k).set(drymass, qv(k) + qitot(k));
    th(k).set(drymass, th(k) - exner(k) * qitot(k) * xxls(k) * inv_cp);
    qitot(k).set(drymass, 0);
    nitot(k).set(drymass, 0);
    qirim(k).set(drymass, 0);
    birim(k).set(drymass, 0);
    if ( not_drymass.any() ) {
      log_hydrometeorsPresent = true; // final update
    }

    drymass = (qitot(k) >= qsmall && qitot(k) < 1.e-8 && t(k) >= zerodegc);
    qr(k).set(drymass, qr(k) + qitot(k));
    th(k).set(drymass, th(k) - exner(k) * qitot(k) * xlf(k) * inv_cp);
    qitot(k).set(drymass, 0);
    nitot(k).set(drymass, 0);
    qirim(k).set(drymass, 0);
    birim(k).set(drymass, 0);

    t(k) = th(k) * inv_exner(k);

    calculate_incloud_mixingratios(
      qc(k), qr(k), qitot(k), qirim(k), nc(k), nr(k), nitot(k), birim(k),
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
  const uview_1d<const Spack>& pres,
  const uview_1d<const Spack>& pdel,
  const uview_1d<const Spack>& dzq,
  const uview_1d<const Spack>& ncnuc,
  const uview_1d<const Spack>& exner,
  const uview_1d<const Spack>& inv_exner,
  const uview_1d<const Spack>& inv_lcldm,
  const uview_1d<const Spack>& inv_icldm,
  const uview_1d<const Spack>& inv_rcldm,
  const uview_1d<const Spack>& naai,
  const uview_1d<const Spack>& qc_relvar,
  const uview_1d<const Spack>& icldm,
  const uview_1d<const Spack>& lcldm,
  const uview_1d<const Spack>& rcldm,
  const uview_1d<Spack>& t,
  const uview_1d<Spack>& rho,
  const uview_1d<Spack>& inv_rho,
  const uview_1d<Spack>& qvs,
  const uview_1d<Spack>& qvi,
  const uview_1d<Spack>& supi,
  const uview_1d<Spack>& rhofacr,
  const uview_1d<Spack>& rhofaci,
  const uview_1d<Spack>& acn,
  const uview_1d<Spack>& qv,
  const uview_1d<Spack>& th,
  const uview_1d<Spack>& qc,
  const uview_1d<Spack>& nc,
  const uview_1d<Spack>& qr,
  const uview_1d<Spack>& nr,
  const uview_1d<Spack>& qitot,
  const uview_1d<Spack>& nitot,
  const uview_1d<Spack>& qirim,
  const uview_1d<Spack>& birim,
  const uview_1d<Spack>& xxlv,
  const uview_1d<Spack>& xxls,
  const uview_1d<Spack>& xlf,
  const uview_1d<Spack>& qc_incld,
  const uview_1d<Spack>& qr_incld,
  const uview_1d<Spack>& qitot_incld,
  const uview_1d<Spack>& qirim_incld,
  const uview_1d<Spack>& nc_incld,
  const uview_1d<Spack>& nr_incld,
  const uview_1d<Spack>& nitot_incld,
  const uview_1d<Spack>& birim_incld,
  const uview_1d<Spack>& mu_c,
  const uview_1d<Spack>& nu,
  const uview_1d<Spack>& lamc,
  const uview_1d<Spack>& cdist,
  const uview_1d<Spack>& cdist1,
  const uview_1d<Spack>& cdistr,
  const uview_1d<Spack>& mu_r,
  const uview_1d<Spack>& lamr,
  const uview_1d<Spack>& logn0r,
  const uview_1d<Spack>& cmeiout,
  const uview_1d<Spack>& prain,
  const uview_1d<Spack>& nevapr,
  const uview_1d<Spack>& prer_evap,
  const uview_1d<Spack>& vap_liq_exchange,
  const uview_1d<Spack>& vap_ice_exchange,
  const uview_1d<Spack>& liq_ice_exchange,
  const uview_1d<Spack>& pratot,
  const uview_1d<Spack>& prctot,
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

  team.team_barrier();
  log_hydrometeorsPresent = false;
  team.team_barrier();

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(team, nk_pack), [&] (Int k) {

    // if relatively dry and no hydrometeors at this level, skip to end of k-loop (i.e. skip this level)
    const auto skip_all = !(qc(k) >= qsmall || qr(k) >= qsmall || qitot(k) >= qsmall) &&
      (t(k) < zerodegc && supi(k) < -0.05);
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
      ncslf   (0), // change in cloud droplet number from self-collection  (Not in paper?)
      ncautc  (0), // change in cloud droplet number from autoconversion
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
      get_time_space_phys_variables(
        t(k), pres(k), rho(k), xxlv(k), xxls(k), qvs(k), qvi(k),
        mu, dv, sc, dqsdt, dqsidt, ab, abi, kap, eii, not_skip_micro);

      get_cloud_dsd2(qc_incld(k), nc_incld(k), mu_c(k), rho(k), nu(k), dnu, lamc(k), cdist(k), cdist1(k), lcldm(k), not_skip_micro);
      nc(k).set(not_skip_micro, nc_incld(k) * lcldm(k));

      get_rain_dsd2(qr_incld(k), nr_incld(k), mu_r(k), lamr(k), cdistr(k), logn0r(k), rcldm(k), not_skip_micro);
      nr(k).set(not_skip_micro, nr_incld(k) * rcldm(k));

      impose_max_total_Ni(nitot_incld(k), max_total_Ni, inv_rho(k), not_skip_micro);

      const auto qitot_gt_small = qitot_incld(k) >= qsmall && not_skip_micro;

      if (qitot_gt_small.any()) {
        // impose lower limits to prevent taking log of # < 0
        nitot_incld(k).set(qitot_gt_small, pack::max(nitot_incld(k), nsmall));
        nr_incld(k).set(qitot_gt_small, pack::max(nr_incld(k), nsmall));

        const auto rhop = calc_bulk_rho_rime(qitot_incld(k), qirim_incld(k), birim_incld(k), qitot_gt_small);

        TableIce table_ice;
        lookup_ice(qitot_incld(k), nitot_incld(k), qirim_incld(k), rhop, table_ice, qitot_gt_small);

        TableRain table_rain;
        lookup_rain(qr_incld(k), nr_incld(k), table_rain, qitot_gt_small);

        // call to lookup table interpolation subroutines to get process rates
        f1pr02.set(qitot_gt_small, apply_table_ice(1, itab, table_ice, qitot_gt_small));
        f1pr03.set(qitot_gt_small, apply_table_ice(2, itab, table_ice, qitot_gt_small));
        f1pr04.set(qitot_gt_small, apply_table_ice(3, itab, table_ice, qitot_gt_small));
        f1pr05.set(qitot_gt_small, apply_table_ice(4, itab, table_ice, qitot_gt_small));
        f1pr09.set(qitot_gt_small, apply_table_ice(6, itab, table_ice, qitot_gt_small));
        f1pr10.set(qitot_gt_small, apply_table_ice(7, itab, table_ice, qitot_gt_small));
        f1pr14.set(qitot_gt_small, apply_table_ice(9, itab, table_ice, qitot_gt_small));

        // ice-rain collection processes
        const auto qr_gt_small = qr_incld(k) >= qsmall && qitot_gt_small;
        f1pr07.set(qr_gt_small, apply_table_coll(0, itabcol, table_ice, table_rain, qitot_gt_small));
        f1pr08.set(qr_gt_small, apply_table_coll(1, itabcol, table_ice, table_rain, qitot_gt_small));

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

      // collection of droplets
      ice_cldliq_collection(
        rho(k), t(k), rhofaci(k), f1pr04, qitot_incld(k), qc_incld(k), nitot_incld(k), nc_incld(k),
        qccol, nccol, qcshd, ncshdc, not_skip_micro);

      // collection of rain
      ice_rain_collection(
        rho(k), t(k), rhofaci(k), logn0r(k), f1pr07, f1pr08, qitot_incld(k), nitot_incld(k), qr_incld(k),
        qrcol, nrcol, not_skip_micro);

      // collection between ice categories

      // PMC nCat deleted lots of stuff here.

      // self-collection of ice
      ice_self_collection(
        rho(k), rhofaci(k), f1pr03, eii, qirim_incld(k), qitot_incld(k), nitot_incld(k),
        nislf, not_skip_micro);

      // melting
      ice_melting(
        rho(k), t(k), pres(k), rhofaci(k), f1pr05, f1pr14, xxlv(k), xlf(k), dv, sc, mu, kap, qv(k), qitot_incld(k), nitot_incld(k),
        qimlt, nimlt, not_skip_micro);

      // calculate wet growth
      ice_cldliq_wet_growth(
        rho(k), t(k), pres(k), rhofaci(k), f1pr05, f1pr14, xxlv(k), xlf(k), dv, kap, mu, sc, qv(k), qc_incld(k), qitot_incld(k), nitot_incld(k), qr_incld(k),
        log_wetgrowth, qrcol, qccol, qwgrth, nrshdr, qcshd, not_skip_micro);

      // calcualte total inverse ice relaxation timescale combined for all ice categories
      // note 'f1pr' values are normalized, so we need to multiply by N
      ice_relaxation_timescale(
        rho(k), t(k), rhofaci(k), f1pr05, f1pr14, dv, mu, sc, qitot_incld(k), nitot_incld(k),
        epsi, epsi_tot, not_skip_micro);

      // calculate rime density
      calc_rime_density(
        t(k), rhofaci(k), f1pr02, acn(k), lamc(k), mu_c(k), qc_incld(k), qccol,
        vtrmi1, rhorime_c, not_skip_micro);

      // contact and immersion freezing droplets
      cldliq_immersion_freezing(
	t(k), lamc(k), mu_c(k), cdist1(k), qc_incld(k), qc_relvar(k),
        qcheti, ncheti, not_skip_micro);

      // for future: get rid of log statements below for rain freezing
      rain_immersion_freezing(
        t(k), lamr(k), mu_r(k), cdistr(k), qr_incld(k),
        qrheti, nrheti, not_skip_micro);

      //  rime splintering (Hallet-Mossop 1974)
      // PMC comment: Morrison and Milbrandt 2015 part 1 and 2016 part 3 both say
      // that Hallet-Mossop should be neglected if 1 category to compensate for
      // artificial smearing out of ice DSD

      // ................................................
      //  condensation/evaporation/deposition/sublimation
      //    (use semi-analytic formulation)

      //  calculate rain evaporation including ventilation
      calc_liq_relaxation_timescale(
        revap_table, rho(k), f1r, f2r, dv, mu, sc, mu_r(k), lamr(k), cdistr(k), cdist(k), qr_incld(k), qc_incld(k),
        epsr, epsc, not_skip_micro);

      evaporate_sublimate_precip(
        qr_incld(k), qc_incld(k), nr_incld(k), qitot_incld(k), lcldm(k), rcldm(k), qvs(k), ab, epsr, qv(k),
        qrevp, nrevp, not_skip_micro);

      ice_deposition_sublimation(
        qitot_incld(k), nitot_incld(k), t(k), qvs(k), qvi(k), epsi, abi, qv(k),
        qidep, qisub, nisub, qiberg, not_skip_micro);
    }

    // deposition/condensation-freezing nucleation
    ice_nucleation(
      t(k), inv_rho(k), nitot(k), naai(k), supi(k), odt, log_predictNc,
      qinuc, ninuc, not_skip_all);

    // cloud water autoconversion
    // NOTE: cloud_water_autoconversion must be called before droplet_self_collection
    cloud_water_autoconversion(
      rho(k), qc_incld(k), nc_incld(k), qc_relvar(k),
      qcaut, ncautc, ncautr, not_skip_all);

    // self-collection of droplets
    droplet_self_collection(
      rho(k), inv_rho(k), qc_incld(k),
      mu_c(k), nu(k), ncautc, ncslf, not_skip_all);

    // accretion of cloud by rain
    cloud_rain_accretion(
      rho(k), inv_rho(k), qc_incld(k), nc_incld(k), qr_incld(k), qc_relvar(k),
      qcacc, ncacc, not_skip_all);

    // self-collection and breakup of rain
    // (breakup following modified Verlinde and Cotton scheme)
    rain_self_collection(
      rho(k), qr_incld(k), nr_incld(k),
      nrslf, not_skip_all);

    // Here we map the microphysics tendency rates back to CELL-AVERAGE quantities for updating
    // cell-average quantities.
    back_to_cell_average(
      lcldm(k), rcldm(k), icldm(k), qcacc, qrevp, qcaut,
      ncacc, ncslf, ncautc, nrslf, nrevp, ncautr, qisub, nrshdr, qcheti,
      qrcol, qcshd, qimlt, qccol, qrheti, nimlt, nccol, ncshdc, ncheti, nrcol, nislf,
      qidep, nrheti, nisub, qinuc, ninuc, qiberg, not_skip_all);

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

    prevent_ice_overdepletion(pres(k), t(k), qv(k), xxls(k), odt, qidep, qisub, not_skip_all);

    // vapor -- not needed, since all sinks already have limits imposed and the sum, therefore,
    //          cannot possibly overdeplete qv

    // cloud
    cloud_water_conservation(
      qc(k), dt,
      qcaut, qcacc, qccol, qcheti, qcshd, qiberg, qisub, qidep, not_skip_all);

    // rain
    rain_water_conservation(
      qr(k), qcaut, qcacc, qimlt, qcshd, dt,
      qrevp, qrcol, qrheti, not_skip_all);

    // ice
    ice_water_conservation(
      qitot(k), qidep, qinuc, qiberg, qrcol, qccol, qrheti, qcheti, dt,
      qisub, qimlt, not_skip_all);

    //---------------------------------------------------------------------------------
    // update prognostic microphysics and thermodynamics variables
    //---------------------------------------------------------------------------------

    //-- ice-phase dependent processes:
    update_prognostic_ice(
      qcheti, qccol, qcshd, nccol, ncheti, ncshdc, qrcol, nrcol,  qrheti, nrheti, nrshdr, qimlt, nimlt, qisub, qidep, qinuc, ninuc, nislf, nisub, qiberg, exner(k), xxls(k), xlf(k), log_predictNc, log_wetgrowth, dt, nmltratio, rhorime_c,
      th(k), qv(k), qitot(k), nitot(k), qirim(k), birim(k), qc(k), nc(k), qr(k), nr(k), not_skip_all);

    //-- warm-phase only processes:
    update_prognostic_liquid(
      qcacc, ncacc, qcaut, ncautc, ncautr, ncslf, qrevp, nrevp, nrslf, log_predictNc, inv_rho(k), exner(k), xxlv(k), dt,
      th(k), qv(k), qc(k), nc(k), qr(k), nr(k), not_skip_all);

    // AaronDonahue - Add extra variables needed from microphysics by E3SM:
    cmeiout(k)         .set(not_skip_all, qidep - qisub + qinuc);
    prain(k)           .set(not_skip_all, qcacc + qcaut + qcshd + qccol);
    nevapr(k)          .set(not_skip_all, qisub + qrevp);
    prer_evap(k)       .set(not_skip_all, qrevp);
    vap_ice_exchange(k).set(not_skip_all, qidep - qisub + qinuc);
    vap_liq_exchange(k).set(not_skip_all, -qrevp);
    liq_ice_exchange(k).set(not_skip_all, qcheti + qrheti - qimlt + qiberg + qccol + qrcol);

    // clipping for small hydrometeor values
    const auto qc_small    = qc(k) < qsmall    && not_skip_all;
    const auto qr_small    = qr(k) < qsmall    && not_skip_all;
    const auto qitot_small = qitot(k) < qsmall && not_skip_all;

    const auto qc_not_small    = qc(k) >= qsmall    && not_skip_all;
    const auto qr_not_small    = qr(k) >= qsmall    && not_skip_all;
    const auto qitot_not_small = qitot(k) >= qsmall && not_skip_all;

    qv(k).set(qc_small, qv(k) + qc(k));
    th(k).set(qc_small, th(k) - exner(k) * qc(k) * xxlv(k) * inv_cp);
    qc(k).set(qc_small, 0);
    nc(k).set(qc_small, 0);

    if (qc_not_small.any()) {
      log_hydrometeorsPresent = true;
    }

    qv(k).set(qr_small, qv(k) + qr(k));
    th(k).set(qr_small, th(k) - exner(k) * qr(k) * xxlv(k) * inv_cp);
    qr(k).set(qr_small, 0);
    nr(k).set(qr_small, 0);

    if (qr_not_small.any()) {
      log_hydrometeorsPresent = true;
    }

    qv(k).set(qitot_small, qv(k) + qitot(k));
    th(k).set(qitot_small, th(k) - exner(k) * qitot(k) * xxls(k) * inv_cp);
    qitot(k).set(qitot_small, 0);
    nitot(k).set(qitot_small, 0);
    qirim(k).set(qitot_small, 0);
    birim(k).set(qitot_small, 0);

    if (qitot_not_small.any()) {
      log_hydrometeorsPresent = true;
    }

    impose_max_total_Ni(nitot(k), max_total_Ni, inv_rho(k), not_skip_all);

    // Outputs associated with aerocom comparison:
    pratot(k).set(not_skip_all, qcacc); // cloud drop accretion by rain
    prctot(k).set(not_skip_all, qcaut); // cloud drop autoconversion to rain

    // Recalculate in-cloud values for sedimentation
    calculate_incloud_mixingratios(
      qc(k), qr(k), qitot(k), qirim(k), nc(k), nr(k), nitot(k), birim(k), inv_lcldm(k), inv_icldm(k), inv_rcldm(k),
      qc_incld(k), qr_incld(k), qitot_incld(k), qirim_incld(k), nc_incld(k), nr_incld(k), nitot_incld(k), birim_incld(k), not_skip_all);

  });
  team.team_barrier();
}

template <typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>
::p3_main_post_main_loop(
  const MemberType& team,
  const Int& nk_pack,
  const view_dnu_table& dnu,
  const view_itab_table& itab,
  const uview_1d<const Spack>& exner,
  const uview_1d<const Spack>& lcldm,
  const uview_1d<const Spack>& rcldm,
  const uview_1d<Spack>& rho,
  const uview_1d<Spack>& inv_rho,
  const uview_1d<Spack>& rhofaci,
  const uview_1d<Spack>& qv,
  const uview_1d<Spack>& th,
  const uview_1d<Spack>& qc,
  const uview_1d<Spack>& nc,
  const uview_1d<Spack>& qr,
  const uview_1d<Spack>& nr,
  const uview_1d<Spack>& qitot,
  const uview_1d<Spack>& nitot,
  const uview_1d<Spack>& qirim,
  const uview_1d<Spack>& birim,
  const uview_1d<Spack>& xxlv,
  const uview_1d<Spack>& xxls,
  const uview_1d<Spack>& mu_c,
  const uview_1d<Spack>& nu,
  const uview_1d<Spack>& lamc,
  const uview_1d<Spack>& mu_r,
  const uview_1d<Spack>& lamr,
  const uview_1d<Spack>& vap_liq_exchange,
  const uview_1d<Spack>& ze_rain,
  const uview_1d<Spack>& ze_ice,
  const uview_1d<Spack>& diag_vmi,
  const uview_1d<Spack>& diag_effi,
  const uview_1d<Spack>& diag_di,
  const uview_1d<Spack>& diag_rhoi,
  const uview_1d<Spack>& diag_ze,
  const uview_1d<Spack>& diag_effc)
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
      const auto qc_gt_small = qc(k) >= qsmall;
      const auto qc_small    = !qc_gt_small;
      get_cloud_dsd2(qc(k), nc(k), mu_c(k), rho(k), nu(k), dnu, lamc(k), ignore1, ignore2, lcldm(k), qc_gt_small);

      diag_effc(k)       .set(qc_gt_small, sp(0.5) * (mu_c(k) + 3) / lamc(k));
      qv(k)              .set(qc_small, qv(k)+qc(k));
      th(k)              .set(qc_small, th(k)-exner(k)*qc(k)*xxlv(k)*inv_cp);
      vap_liq_exchange(k).set(qc_small, vap_liq_exchange(k) - qc(k));
      qc(k)              .set(qc_small, 0);
      nc(k)              .set(qc_small, 0);
    }

    // Rain
    {
      const auto qr_gt_small = qr(k) >= qsmall;
      const auto qr_small    = !qr_gt_small;
      get_rain_dsd2(
        qr(k), nr(k), mu_r(k), lamr(k), ignore1, ignore2, rcldm(k), qr_gt_small);

      ze_rain(k).set(qr_gt_small, nr(k)*(mu_r(k)+6)*(mu_r(k)+5)*(mu_r(k)+4)*
                     (mu_r(k)+3)*(mu_r(k)+2)*(mu_r(k)+1)/pow(lamr(k), sp(6.0))); // once f90 is gone, 6 can be int
      ze_rain(k).set(qr_gt_small, pack::max(ze_rain(k), sp(1.e-22)));

      qv(k)              .set(qr_small, qv(k) + qr(k));
      th(k)              .set(qr_small, th(k) - exner(k)*qr(k)*xxlv(k)*inv_cp);
      vap_liq_exchange(k).set(qr_small, vap_liq_exchange(k) - qr(k));
      qr(k)              .set(qr_small, 0);
      nr(k)              .set(qr_small, 0);
    }

    // Ice
    {
      impose_max_total_Ni(nitot(k), max_total_Ni, inv_rho(k));

      const auto qi_gt_small = qitot(k) >= qsmall;
      const auto qi_small    = !qi_gt_small;

      // impose lower limits to prevent taking log of # < 0
      nitot(k) = pack::max(nitot(k), nsmall);
      nr(k)    = pack::max(nr(k), nsmall);

      const auto rhop = calc_bulk_rho_rime(qitot(k), qirim(k), birim(k), qi_gt_small);

      TableIce table_ice;
      lookup_ice(qitot(k), nitot(k), qirim(k), rhop, table_ice, qi_gt_small);

      f1pr02.set(qi_gt_small, apply_table_ice(1,  itab, table_ice, qi_gt_small));
      f1pr06.set(qi_gt_small, apply_table_ice(5,  itab, table_ice, qi_gt_small));
      f1pr09.set(qi_gt_small, apply_table_ice(6,  itab, table_ice, qi_gt_small));
      f1pr10.set(qi_gt_small, apply_table_ice(7,  itab, table_ice, qi_gt_small));
      f1pr13.set(qi_gt_small, apply_table_ice(8,  itab, table_ice, qi_gt_small));
      f1pr15.set(qi_gt_small, apply_table_ice(10, itab, table_ice, qi_gt_small));
      f1pr16.set(qi_gt_small, apply_table_ice(11, itab, table_ice, qi_gt_small));

      // impose mean ice size bounds (i.e. apply lambda limiters)
      // note that the Nmax and Nmin are normalized and thus need to be multiplied by existing N
      nitot(k).set(qi_gt_small, pack::min(nitot(k), f1pr09 * nitot(k)));
      nitot(k).set(qi_gt_small, pack::max(nitot(k), f1pr10 * nitot(k)));

      // --this should already be done in s/r 'calc_bulkRhoRime'
      const auto qirim_small = qirim(k) < qsmall && qi_gt_small;
      qirim(k).set(qirim_small, 0);
      birim(k).set(qirim_small, 0);

      // note that reflectivity from lookup table is normalized, so we need to multiply by N
      diag_vmi(k) .set(qi_gt_small, f1pr02 * rhofaci(k));
      diag_effi(k).set(qi_gt_small, f1pr06); // units are in m
      diag_di(k)  .set(qi_gt_small, f1pr15);
      diag_rhoi(k).set(qi_gt_small, f1pr16);

      // note factor of air density below is to convert from m^6/kg to m^6/m^3
      ze_ice(k).set(qi_gt_small, ze_ice(k) + sp(0.1892)*f1pr13*nitot(k)*rho(k)); // sum contribution from each ice category (note: 0.1892 = 0.176/0.93);
      ze_ice(k).set(qi_gt_small, pack::max(ze_ice(k), sp(1.e-22)));

      qv(k)     .set(qi_small, qv(k) + qitot(k));
      th(k)     .set(qi_small, th(k) - exner(k)*qitot(k)*xxls(k)*inv_cp);
      qitot(k)  .set(qi_small, 0);
      nitot(k)  .set(qi_small, 0);
      qirim(k)  .set(qi_small, 0);
      birim(k)  .set(qi_small, 0);
      diag_di(k).set(qi_small, 0);
    }

    // sum ze components and convert to dBZ
    diag_ze(k) = 10 * pack::log10((ze_rain(k) + ze_ice(k))*sp(1.e18));

    // if qr is very small then set Nr to 0 (needs to be done here after call
    // to ice lookup table because a minimum Nr of nsmall will be set otherwise even if qr=0)
    nr(k).set(qr(k) < qsmall, 0);
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
  const view_2d<const Spack>& ncnuc,         // IN ccn activated number tendency     kg-1 s-1
  const view_2d<const Spack>& naai,          // IN actived ice nuclei concentration  1/kg
  const view_2d<const Spack>& qc_relvar,     // Assumed SGS 1/(var(qc)/mean(qc))     kg2/kg2
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
  const view_2d<Spack>& vap_ice_exchange) // sum of vap-ice phase change tendenices
{
  using ExeSpace = typename KT::ExeSpace;

  view_2d<Spack> xxls("xxls", ni, nk), xxlv("xxlv", ni, nk), xlf("xlf", ni, nk);

  get_latent_heat(ni, nk, xxlv, xxls, xlf);

  const Int nk_pack = scream::pack::npack<Spack>(nk);
  const auto policy = util::ExeSpaceUtils<ExeSpace>::get_default_team_policy(ni, nk_pack);

  WorkspaceManager<Spack, Device> workspace_mgr(nk_pack, 100, policy);

  // load constants into local vars
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

  // per-column bools
  view_2d<bool> bools("bools", ni, 2);

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
    const auto oncnuc            = util::subview(ncnuc, i);
    const auto onaai             = util::subview(naai, i);
    const auto oqc_relvar        = util::subview(qc_relvar, i);
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
    const auto oxxlv             = util::subview(xxlv, i);
    const auto oxxls             = util::subview(xxls, i);
    const auto oxlf              = util::subview(xlf, i);

    // Need to watch out for race conditions with these shared variables
    bool &log_nucleationPossible  = bools(i, 0);
    bool &log_hydrometeorsPresent = bools(i, 1);

    // initialize
    p3_main_init(
      team, nk_pack,
      oicldm, olcldm, orcldm, oexner, oth,
      opratot, oprctot, prec, mu_r, odiag_ze, ze_ice, ze_rain, odiag_effc, odiag_effi, odiag_vmi, odiag_di, odiag_rhoi, ocmeiout, oprain, onevapr, orflx, osflx, inv_icldm, inv_lcldm, inv_rcldm, omu_c, olamc, inv_exner, t, oqv,
      prt_liq(i), prt_sol(i));

    p3_main_pre_main_loop(
      team, nk, log_predictNc, dt,
      opres, opdel, odzq, oncnuc, oexner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, oxxlv, oxxls, oxlf,
      t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, acn, oqv, oth, oqc, onc, oqr, onr, oqitot, onitot, oqirim, obirim, qc_incld, qr_incld, qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld,
      log_nucleationPossible, log_hydrometeorsPresent);

    // There might not be any work to do for this team
    if (!(log_nucleationPossible || log_hydrometeorsPresent)) {
      return; // this is how you do a "continue" in a kokkos lambda
    }

    // ------------------------------------------------------------------------------------------
    // main k-loop (for processes):
    p3_main_main_loop(
      team, nk_pack, log_predictNc, dt, odt,
      dnu, itab, itabcol, revap_table,
      opres, opdel, odzq, oncnuc, oexner, inv_exner, inv_lcldm, inv_icldm, inv_rcldm, onaai, oqc_relvar, oicldm, olcldm, orcldm,
      t, rho, inv_rho, qvs, qvi, supi, rhofacr, rhofaci, acn, oqv, oth, oqc, onc, oqr, onr, oqitot, onitot, oqirim, obirim, oxxlv, oxxls, oxlf, qc_incld, qr_incld, qitot_incld, qirim_incld, nc_incld, nr_incld, nitot_incld, birim_incld, omu_c, nu, olamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r, ocmeiout, oprain, onevapr, oprer_evap, ovap_liq_exchange, ovap_ice_exchange, oliq_ice_exchange, opratot, oprctot,
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
      team, nk_pack, dnu, itab,
      oexner, olcldm, orcldm,
      rho, inv_rho, rhofaci, oqv, oth, oqc, onc, oqr, onr, oqitot, onitot, oqirim, obirim, oxxlv, oxxls, omu_c, nu, olamc, mu_r, lamr, ovap_liq_exchange, ze_rain, ze_ice, odiag_vmi, odiag_effi, odiag_di, odiag_rhoi, odiag_ze, odiag_effc);

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

    check_values(oqv, tmparr1, ktop, kbot, it, debug_ABORT, 900, team, ocol_location);
#endif

  });
}

} // namespace p3
} // namespace scream

#endif
