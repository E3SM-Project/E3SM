
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

template <>
void Functions<Real,DefaultDevice>
::p3_main_init_disp(
  const Int& nj, const Int& nk_pack, const uview_2d<const Spack>& cld_frac_i, const uview_2d<const Spack>& cld_frac_l,
  const uview_2d<const Spack>& cld_frac_r, const uview_2d<const Spack>& inv_exner, const uview_2d<const Spack>& th_atm,
  const uview_2d<const Spack>& dz, const uview_2d<Spack>& diag_equiv_reflectivity, const uview_2d<Spack>& ze_ice,
  const uview_2d<Spack>& ze_rain, const uview_2d<Spack>& diag_eff_radius_qc, const uview_2d<Spack>& diag_eff_radius_qi,
  const uview_2d<Spack>& diag_eff_radius_qr,
  const uview_2d<Spack>& inv_cld_frac_i, const uview_2d<Spack>& inv_cld_frac_l, const uview_2d<Spack>& inv_cld_frac_r,
  const uview_2d<Spack>& exner, const uview_2d<Spack>& T_atm, const uview_2d<Spack>& qv, const uview_2d<Spack>& inv_dz,
  const uview_1d<Scalar>& precip_liq_surf, const uview_1d<Scalar>& precip_ice_surf,
  const uview_2d<Spack>& mu_r, const uview_2d<Spack>& lamr, const uview_2d<Spack>& logn0r, const uview_2d<Spack>& nu,
  const uview_2d<Spack>& cdist, const uview_2d<Spack>& cdist1, const uview_2d<Spack>& cdistr,
  const uview_2d<Spack>& qc_incld, const uview_2d<Spack>& qr_incld, const uview_2d<Spack>& qi_incld,
  const uview_2d<Spack>& qm_incld, const uview_2d<Spack>& nc_incld, const uview_2d<Spack>& nr_incld, const uview_2d<Spack>& ni_incld,
  const uview_2d<Spack>& bm_incld, const uview_2d<Spack>& inv_rho, const uview_2d<Spack>& prec, const  uview_2d<Spack>& rho, const uview_2d<Spack>& rhofacr,
  const uview_2d<Spack>& rhofaci, const uview_2d<Spack>& acn, const uview_2d<Spack>& qv_sat_l, const uview_2d<Spack>& qv_sat_i, const uview_2d<Spack>& sup,
  const uview_2d<Spack>& qv_supersat_i, const uview_2d<Spack>& qtend_ignore, const uview_2d<Spack>& ntend_ignore, const uview_2d<Spack>& mu_c,
  const uview_2d<Spack>& lamc, const uview_2d<Spack>& rho_qi, const uview_2d<Spack>& qv2qi_depos_tend, const uview_2d<Spack>& precip_total_tend,
  const uview_2d<Spack>& nevapr, const uview_2d<Spack>& precip_liq_flux, const uview_2d<Spack>& precip_ice_flux)
{
  using ExeSpace = typename KT::ExeSpace;
  const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(nj, nk_pack);

  Kokkos::parallel_for("p3_main_init",
         policy, KOKKOS_LAMBDA(const MemberType& team) {

    const Int i = team.league_rank();
    precip_liq_surf(i) = 0;
    precip_ice_surf(i) = 0;

    Kokkos::parallel_for(
      Kokkos::TeamVectorRange(team, nk_pack), [&] (Int k) {
        diag_equiv_reflectivity(i,k) = -99;
        ze_ice(i,k)            = 1.e-22;
        ze_rain(i,k)           = 1.e-22;
        diag_eff_radius_qc(i,k)         = 10.e-6;
        diag_eff_radius_qi(i,k)         = 25.e-6;
        diag_eff_radius_qr(i,k)         = 500.e-6;
        inv_cld_frac_i(i,k)    = 1 / cld_frac_i(i,k);
        inv_cld_frac_l(i,k)    = 1 / cld_frac_l(i,k);
        inv_cld_frac_r(i,k)    = 1 / cld_frac_r(i,k);
        exner(i,k)         = 1 / inv_exner(i,k);
        T_atm(i,k)                 = th_atm(i,k) * exner(i,k);
        qv(i,k)                = max(qv(i,k), 0);
        inv_dz(i,k)            = 1 / dz(i,k);
        mu_r(i,k)               = 0.;
        lamr(i,k)               = 0.;
        logn0r(i,k)             = 0.;
        nu(i,k)                 = 0.;
        cdist(i,k)              = 0.;
        cdist1(i,k)             = 0.;
        cdistr(i,k)             = 0.;
        qc_incld(i,k)           = 0.;
        qr_incld(i,k)           = 0.;
        qi_incld(i,k)           = 0.;
        qm_incld(i,k)           = 0.;
        nc_incld(i,k)           = 0.;
        nr_incld(i,k)           = 0.;
        ni_incld(i,k)           = 0.;
        bm_incld(i,k)           = 0.;
        inv_rho(i,k)            = 0.;
        prec(i,k)               = 0.;
        rho(i,k)                = 0.;
        rhofacr(i,k)            = 0.;
        rhofaci(i,k)            = 0.;
        acn(i,k)                = 0.;
        qv_sat_l(i,k)           = 0.;
        qv_sat_i(i,k)           = 0.;
        sup(i,k)                = 0.;
        qv_supersat_i(i,k)      = 0.;
        qtend_ignore(i,k)       = 0.;
        ntend_ignore(i,k)       = 0.;
        mu_c(i,k)               = 0.;
        lamc(i,k)               = 0.;
        rho_qi(i,k)             = 0.;
        qv2qi_depos_tend(i,k)   = 0.;
        precip_total_tend(i,k)  = 0.;
        nevapr(i,k)             = 0.;
        precip_liq_flux(i,k)    = 0.;
        precip_ice_flux(i,k)    = 0.;
   });
 });
}

template <>
Int Functions<Real,DefaultDevice>
::p3_main_internal_disp(
  const P3Runtime& runtime_options,
  const P3PrognosticState& prognostic_state,
  const P3DiagnosticInputs& diagnostic_inputs,
  const P3DiagnosticOutputs& diagnostic_outputs,
  const P3Infrastructure& infrastructure,
  const P3HistoryOnly& history_only,
  const P3LookupTables& lookup_tables,
  const P3Temporaries& temporaries,
  const WorkspaceManager& workspace_mgr,
  Int nj,
  Int nk)
{
  using ExeSpace = typename KT::ExeSpace;

  const Int nk_pack = ekat::npack<Spack>(nk);

  // load constants into local vars
  const     Scalar inv_dt          = 1 / infrastructure.dt;
  constexpr Int    kdir         = -1;
  const     Int    ktop         = kdir == -1 ? 0    : nk-1;
  const     Int    kbot         = kdir == -1 ? nk-1 : 0;
  constexpr bool   debug_ABORT  = false;

  const bool do_ice_production = runtime_options.do_ice_production;

  // per-column bools
  view_1d<bool> nucleationPossible("nucleationPossible", nj);
  view_1d<bool> hydrometeorsPresent("hydrometeorsPresent", nj);

  // Get views of all inputs
  auto pres                    = diagnostic_inputs.pres;
  auto dz                      = diagnostic_inputs.dz;
  auto nc_nuceat_tend          = diagnostic_inputs.nc_nuceat_tend;
  auto nccn_prescribed         = diagnostic_inputs.nccn;
  auto ni_activated            = diagnostic_inputs.ni_activated;
  auto inv_qc_relvar           = diagnostic_inputs.inv_qc_relvar;
  auto dpres                   = diagnostic_inputs.dpres;
  auto inv_exner               = diagnostic_inputs.inv_exner;
  auto cld_frac_i              = diagnostic_inputs.cld_frac_i;
  auto cld_frac_l              = diagnostic_inputs.cld_frac_l;
  auto cld_frac_r              = diagnostic_inputs.cld_frac_r;
  auto col_location            = infrastructure.col_location;
  auto qc                      = prognostic_state.qc;
  auto nc                      = prognostic_state.nc;
  auto qr                      = prognostic_state.qr;
  auto nr                      = prognostic_state.nr;
  auto qi                      = prognostic_state.qi;
  auto qm                      = prognostic_state.qm;
  auto ni                      = prognostic_state.ni;
  auto bm                      = prognostic_state.bm;
  auto qv                      = prognostic_state.qv;
  auto th                      = prognostic_state.th;
  auto diag_eff_radius_qc      = diagnostic_outputs.diag_eff_radius_qc;
  auto diag_eff_radius_qi      = diagnostic_outputs.diag_eff_radius_qi;
  auto diag_eff_radius_qr      = diagnostic_outputs.diag_eff_radius_qr;
  auto qv2qi_depos_tend        = diagnostic_outputs.qv2qi_depos_tend;
  auto rho_qi                  = diagnostic_outputs.rho_qi;
  auto precip_liq_flux         = diagnostic_outputs.precip_liq_flux;
  auto precip_ice_flux         = diagnostic_outputs.precip_ice_flux;
  auto precip_total_tend       = diagnostic_outputs.precip_total_tend;
  auto nevapr                  = diagnostic_outputs.nevapr;
  auto qv_prev                 = diagnostic_inputs.qv_prev;
  auto t_prev                  = diagnostic_inputs.t_prev;
  // Inputs for the heteogeneous freezing
  auto hetfrz_immersion_nucleation_tend  = diagnostic_inputs.hetfrz_immersion_nucleation_tend;
  auto hetfrz_contact_nucleation_tend    = diagnostic_inputs.hetfrz_contact_nucleation_tend;
  auto hetfrz_deposition_nucleation_tend = diagnostic_inputs.hetfrz_deposition_nucleation_tend;

  auto liq_ice_exchange        = history_only.liq_ice_exchange;
  auto vap_liq_exchange        = history_only.vap_liq_exchange;
  auto vap_ice_exchange        = history_only.vap_ice_exchange;
  auto mu_r                    = temporaries.mu_r;
  auto T_atm                   = temporaries.T_atm;
  auto lamr                    = temporaries.lamr;
  auto logn0r                  = temporaries.logn0r;
  auto nu                      = temporaries.nu;
  auto cdist                   = temporaries.cdist;
  auto cdist1                  = temporaries.cdist1;
  auto cdistr                  = temporaries.cdistr;
  auto inv_cld_frac_i          = temporaries.inv_cld_frac_i;
  auto inv_cld_frac_l          = temporaries.inv_cld_frac_l;
  auto inv_cld_frac_r          = temporaries.inv_cld_frac_r;
  auto qc_incld                = temporaries.qc_incld;
  auto qr_incld                = temporaries.qr_incld;
  auto qi_incld                = temporaries.qi_incld;
  auto qm_incld                = temporaries.qm_incld;
  auto nc_incld                = temporaries.nc_incld;
  auto nr_incld                = temporaries.nr_incld;
  auto ni_incld                = temporaries.ni_incld;
  auto bm_incld                = temporaries.bm_incld;
  auto inv_dz                  = temporaries.inv_dz;
  auto inv_rho                 = temporaries.inv_rho;
  auto ze_ice                  = temporaries.ze_ice;
  auto ze_rain                 = temporaries.ze_rain;
  auto prec                    = temporaries.prec;
  auto rho                     = temporaries.rho;
  auto rhofacr                 = temporaries.rhofacr;
  auto rhofaci                 = temporaries.rhofaci;
  auto acn                     = temporaries.acn;
  auto qv_sat_l                = temporaries.qv_sat_l;
  auto qv_sat_i                = temporaries.qv_sat_i;
  auto sup                     = temporaries.sup;
  auto qv_supersat_i           = temporaries.qv_supersat_i;
  auto tmparr2                 = temporaries.tmparr2;
  auto exner                   = temporaries.exner;
  auto diag_equiv_reflectivity = temporaries.diag_equiv_reflectivity;
  auto diag_vm_qi              = temporaries.diag_vm_qi;
  auto diag_diam_qi            = temporaries.diag_diam_qi;
  auto pratot                  = temporaries.pratot;
  auto prctot                  = temporaries.prctot;
  auto qtend_ignore            = temporaries.qtend_ignore;
  auto ntend_ignore            = temporaries.ntend_ignore;
  auto mu_c                    = temporaries.mu_c;
  auto lamc                    = temporaries.lamc;
  auto qr_evap_tend            = temporaries.qr_evap_tend;
  auto v_qc                    = temporaries.v_qc;
  auto v_nc                    = temporaries.v_nc;
  auto flux_qx                 = temporaries.flux_qx;
  auto flux_nx                 = temporaries.flux_nx;
  auto v_qit                   = temporaries.v_qit;
  auto v_nit                   = temporaries.v_nit;
  auto flux_nit                = temporaries.flux_nit;
  auto flux_bir                = temporaries.flux_bir;
  auto flux_qir                = temporaries.flux_qir;
  auto flux_qit                = temporaries.flux_qit;
  auto v_qr                    = temporaries.v_qr;
  auto v_nr                    = temporaries.v_nr;

  // we do not want to measure init stuff
  auto start = std::chrono::steady_clock::now();

  // initialize
  p3_main_init_disp(
      nj, nk_pack, cld_frac_i, cld_frac_l, cld_frac_r, inv_exner, th, dz, diag_equiv_reflectivity,
      ze_ice, ze_rain, diag_eff_radius_qc, diag_eff_radius_qi, diag_eff_radius_qr,
      inv_cld_frac_i, inv_cld_frac_l, inv_cld_frac_r, exner, T_atm, qv, inv_dz,
      diagnostic_outputs.precip_liq_surf, diagnostic_outputs.precip_ice_surf,
      mu_r, lamr, logn0r, nu, cdist, cdist1, cdistr,
      qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld,
      inv_rho, prec, rho, rhofacr, rhofaci, acn, qv_sat_l, qv_sat_i, sup, qv_supersat_i,
      qtend_ignore, ntend_ignore, mu_c, lamc, rho_qi, qv2qi_depos_tend, precip_total_tend,
      nevapr, precip_liq_flux, precip_ice_flux);

  p3_main_part1_disp(
      nj, nk, infrastructure.predictNc, infrastructure.prescribedCCN, infrastructure.dt,
      pres, dpres, dz, nc_nuceat_tend, nccn_prescribed, inv_exner, exner, inv_cld_frac_l, inv_cld_frac_i,
      inv_cld_frac_r,
      T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr,
      rhofaci, acn, qv, th, qc, nc, qr, nr, qi, ni, qm,
      bm, qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld,
      ni_incld, bm_incld, nucleationPossible, hydrometeorsPresent, runtime_options);

  // ------------------------------------------------------------------------------------------
  // main k-loop (for processes):

  p3_main_part2_disp(
      nj, nk, runtime_options.max_total_ni, infrastructure.predictNc, infrastructure.prescribedCCN, infrastructure.dt, inv_dt,
      hetfrz_immersion_nucleation_tend, hetfrz_contact_nucleation_tend, hetfrz_deposition_nucleation_tend,
      lookup_tables.dnu_table_vals, lookup_tables.ice_table_vals, lookup_tables.collect_table_vals,
      lookup_tables.revap_table_vals, pres, dpres, dz, nc_nuceat_tend, inv_exner,
      exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r, ni_activated, inv_qc_relvar, cld_frac_i,
      cld_frac_l, cld_frac_r, qv_prev, t_prev, T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn,
      qv, th, qc, nc, qr, nr, qi, ni, qm, bm, qc_incld, qr_incld, qi_incld, qm_incld, nc_incld,
      nr_incld, ni_incld, bm_incld, mu_c, nu, lamc, cdist, cdist1, cdistr,
      mu_r, lamr, logn0r, qv2qi_depos_tend, precip_total_tend, nevapr, qr_evap_tend,
      vap_liq_exchange, vap_ice_exchange, liq_ice_exchange,
      pratot, prctot, nucleationPossible, hydrometeorsPresent, runtime_options);

  //NOTE: At this point, it is possible to have negative (but small) nc, nr, ni.  This is not
  //      a problem; those values get clipped to zero in the sedimentation section (if necessary).
  //      (This is not done above simply for efficiency purposes.)

  // -----------------------------------------------------------------------------------------
  // End of main microphysical processes section
  // =========================================================================================

  // ==========================================================================================!
  // Sedimentation:

  // Cloud sedimentation:  (adaptive substepping)
  cloud_sedimentation_disp(
      qc_incld, rho, inv_rho, cld_frac_l, acn, inv_dz, lookup_tables.dnu_table_vals, workspace_mgr,
      nj, nk, ktop, kbot, kdir, infrastructure.dt, inv_dt, infrastructure.predictNc,
      qc, nc, nc_incld, mu_c, lamc, qtend_ignore, ntend_ignore,
      diagnostic_outputs.precip_liq_surf, nucleationPossible, hydrometeorsPresent);


  // Rain sedimentation:  (adaptive substepping)
  rain_sedimentation_disp(
      rho, inv_rho, rhofacr, cld_frac_r, inv_dz, qr_incld, workspace_mgr,
      lookup_tables.vn_table_vals, lookup_tables.vm_table_vals, nj, nk, ktop, kbot, kdir, infrastructure.dt, inv_dt, qr,
      nr, nr_incld, mu_r, lamr, precip_liq_flux, qtend_ignore, ntend_ignore,
      diagnostic_outputs.precip_liq_surf, nucleationPossible, hydrometeorsPresent, runtime_options);

  // Ice sedimentation:  (adaptive substepping)
  ice_sedimentation_disp(
      rho, inv_rho, rhofaci, cld_frac_i, inv_dz, workspace_mgr, nj, nk, ktop, kbot,
      kdir, infrastructure.dt, inv_dt, qi, qi_incld, ni, ni_incld,
      qm, qm_incld, bm, bm_incld, qtend_ignore, ntend_ignore,
      lookup_tables.ice_table_vals, diagnostic_outputs.precip_ice_surf, nucleationPossible, hydrometeorsPresent, runtime_options);

  // homogeneous freezing f cloud and rain
  if(do_ice_production) {
    homogeneous_freezing_disp(T_atm, inv_exner, nj, nk, ktop, kbot, kdir, qc,
                              nc, qr, nr, qi, ni, qm, bm, th,
                              nucleationPossible, hydrometeorsPresent);
  }

  //
  // final checks to ensure consistency of mass/number
  // and compute diagnostic fields for output
  //
  p3_main_part3_disp(
      nj, nk_pack, runtime_options.max_total_ni, lookup_tables.dnu_table_vals, lookup_tables.ice_table_vals, inv_exner, cld_frac_l, cld_frac_r, cld_frac_i,
      rho, inv_rho, rhofaci, qv, th, qc, nc, qr, nr, qi, ni,
      qm, bm, mu_c, nu, lamc, mu_r, lamr,
      vap_liq_exchange, ze_rain, ze_ice, diag_vm_qi, diag_eff_radius_qi, diag_diam_qi,
      rho_qi, diag_equiv_reflectivity, diag_eff_radius_qc, diag_eff_radius_qr, nucleationPossible, hydrometeorsPresent,
      runtime_options);

  //
  // merge ice categories with similar properties

  //   note:  this should be relocated to above, such that the diagnostic
  //          ice properties are computed after merging

  // PMC nCat deleted nCat>1 stuff

#ifndef NDEBUG
  Kokkos::parallel_for(
      Kokkos::MDRangePolicy<ExeSpace, Kokkos::Rank<2>>({0, 0}, {nj, nk_pack}), KOKKOS_LAMBDA (int i, int k) {
      tmparr2(i,k) = th(i,k) * exner(i,k);
  });
  check_values_disp(qv, tmparr2, ktop, kbot, infrastructure.it, debug_ABORT, 900, col_location, nj, nk);
#endif
  Kokkos::fence();

  auto finish = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(finish - start);
  return duration.count();
}

} // namespace p3
} // namespace scream
