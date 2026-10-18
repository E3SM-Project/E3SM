#include "catch2/catch.hpp"
//[unittest
#include "p3_functions.hpp"
#include "p3_unit_tests_common.hpp"

#include "share/physics/physics_functions.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_team_policy_utils.hpp>

namespace scream {
namespace p3 {
namespace unit_test {

/*
 * Unit tests for the new "prognostic supersaturation" pathways added to
 * p3_main_part1 (aerosol activation, marked "[shanyp") and p3_main_part2
 * (condensation/evaporation, marked "[shanyp") in:
 *   impl/p3_main_impl_part1.hpp
 *   impl/p3_main_impl_part2.hpp
 *
 * These tests exercise p3_main_part1/p3_main_part2 directly (single column,
 * one Kokkos pack) rather than through the legacy F90-comparison ("_host")
 * bridges, since those bridges (tests/infra/p3_test_data.cpp) have not yet
 * been updated to pass the new arguments (qv2qc_aeroact, omega_mp, tke_mp,
 * qv2qc_condevap) added to p3_main_part1/p3_main_part2's signatures.
 *
 * Rather than re-deriving the (lookup-table-dependent) numeric formulas used
 * internally, these are "differential" property tests: for the same inputs,
 * we run the kernel once with the new feature switched off (p3_super_sat =
 * false, matching pre-existing behavior) and once with it switched on, and
 * check that:
 *   (a) the new diagnostic outputs (qv2qc_aeroact / qv2qc_condevap) are
 *       nonzero only when the feature is active and its trigger condition
 *       is met,
 *   (b) qc/qv (and, where applicable, th_atm) are updated consistently with
 *       those diagnostics (conservation of total water / consistent latent
 *       heating), and
 *   (c) when the trigger condition is not met, or the feature is off, the
 *       state is bit-for-bit identical to the "off" run.
 */
template <typename D>
struct UnitWrap::UnitTest<D>::TestPrognosticSupersat : public UnitWrap::UnitTest<D>::Base {

  using TPF     = ekat::TeamPolicyFactory<ExeSpace>;
  using physics = scream::physics::Functions<Scalar, Device>;

  // Bag of results pulled back from a single-column p3_main_part1 call.
  struct Part1Result {
    Real qv, qc, nc, th_atm, qv2qc_aeroact;
    bool is_nucleat_possible, is_hydromet_present;
  };

  // Bag of results pulled back from a single-column p3_main_part2 call.
  struct Part2Result {
    Real qv, qc, th_atm, qv2qc_condevap;
    bool is_hydromet_present;
  };

  // Build a device view of length n, every pack lane set to val.
  static view_1d<Pack> make_view(Real val, Int n) {
    view_1d<Pack> v("v", n);
    auto vh = Kokkos::create_mirror_view(v);
    for (Int k = 0; k < n; ++k) vh(k) = Pack(val);
    Kokkos::deep_copy(v, vh);
    return v;
  }

  //-----------------------------------------------------------------------
  // Run p3_main_part1 for a single column/level (broadcast across the pack)
  // and return the resulting state. aeroact_sentinel lets us pre-seed the
  // qv2qc_aeroact diagnostic with a non-zero value to probe whether it gets
  // properly reset when the activation branch is not taken.
  //-----------------------------------------------------------------------
  Part1Result run_part1(bool p3_super_sat, bool shoc_enable_condensation,
                         Real T_atm_in, Real pres_in, Real qv_in, Real qc_in, Real nc_in,
                         Real nccn_prescribed_in, Real aeroact_sentinel = 0)
  {
    constexpr Scalar g   = C::gravit.value;
    const Int nk         = Pack::n;
    const Int nk_pack     = 1;

    const Scalar dt = 1800.0;

    typename Functions::P3Runtime runtime_options; // value-initialized -> all new bools default false
    runtime_options.p3_super_sat            = p3_super_sat;
    runtime_options.shoc_enable_condensation = shoc_enable_condensation;

    auto pres            = make_view(pres_in, nk_pack);
    auto dpres           = make_view(1000.0, nk_pack);
    auto dz              = make_view(1000.0 / (1.0 * g), nk_pack); // -> rho ~ 1.0
    auto nc_nuceat_tend  = make_view(0.0, nk_pack);
    auto nccn_prescribed = make_view(nccn_prescribed_in, nk_pack);
    auto inv_exner       = make_view(1.0, nk_pack);
    auto exner           = make_view(1.0, nk_pack);
    auto inv_cld_frac_l  = make_view(1.0, nk_pack); // fully cloudy column
    auto inv_cld_frac_i  = make_view(1.0, nk_pack);
    auto inv_cld_frac_r  = make_view(1.0, nk_pack);

    auto T_atm     = make_view(T_atm_in, nk_pack);
    auto rho       = make_view(0.0, nk_pack);
    auto inv_rho   = make_view(0.0, nk_pack);
    auto qv_sat_l  = make_view(0.0, nk_pack);
    auto qv_sat_i  = make_view(0.0, nk_pack);
    auto qv_supersat_i = make_view(0.0, nk_pack);
    auto rhofacr   = make_view(0.0, nk_pack);
    auto rhofaci   = make_view(0.0, nk_pack);
    auto acn       = make_view(0.0, nk_pack);
    auto qv        = make_view(qv_in, nk_pack);
    auto th_atm    = make_view(300.0, nk_pack);
    auto qc        = make_view(qc_in, nk_pack);
    auto nc        = make_view(nc_in, nk_pack);
    auto qr        = make_view(0.0, nk_pack);
    auto nr        = make_view(0.0, nk_pack);
    auto qi        = make_view(0.0, nk_pack);
    auto ni        = make_view(0.0, nk_pack);
    auto qm        = make_view(0.0, nk_pack);
    auto bm        = make_view(0.0, nk_pack);
    auto qc_incld  = make_view(0.0, nk_pack);
    auto qr_incld  = make_view(0.0, nk_pack);
    auto qi_incld  = make_view(0.0, nk_pack);
    auto qm_incld  = make_view(0.0, nk_pack);
    auto nc_incld  = make_view(0.0, nk_pack);
    auto nr_incld  = make_view(0.0, nk_pack);
    auto ni_incld  = make_view(0.0, nk_pack);
    auto bm_incld  = make_view(0.0, nk_pack);
    auto qv2qc_aeroact = make_view(aeroact_sentinel, nk_pack);

    view_1d<bool> bools_d("bools", 2);

    // do_predict_nc = false, do_prescribed_CCN = false: droplet number in the
    // (unmodified) "else" leg of p3_main_part1 comes from a fixed nccnst, so
    // that leg's behavior is completely independent of p3_super_sat and lets
    // us isolate the effect of the new activation block cleanly.
    const bool do_predict_nc     = false;
    const bool do_prescribed_CCN = false;

    auto policy = TPF::get_default_team_policy(1, nk_pack);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      Functions::p3_main_part1(
        team, nk, do_predict_nc, do_prescribed_CCN, dt,
        pres, dpres, dz, nc_nuceat_tend, nccn_prescribed, inv_exner, exner,
        inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
        T_atm, rho, inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci,
        acn, qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm,
        qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld,
        qv2qc_aeroact, bools_d(0), bools_d(1), runtime_options);
    });

    Part1Result result;
    auto get0 = [](const view_1d<Pack>& v) {
      auto vh = Kokkos::create_mirror_view(v);
      Kokkos::deep_copy(vh, v);
      return vh(0)[0];
    };
    result.qv            = get0(qv);
    result.qc            = get0(qc);
    result.nc            = get0(nc);
    result.th_atm        = get0(th_atm);
    result.qv2qc_aeroact  = get0(qv2qc_aeroact);

    auto bools_h = Kokkos::create_mirror_view(bools_d);
    Kokkos::deep_copy(bools_h, bools_d);
    result.is_nucleat_possible = bools_h(0);
    result.is_hydromet_present = bools_h(1);

    return result;
  }

  //-----------------------------------------------------------------------
  // Run p3_main_part2 for a single column/level and return the resulting
  // state. condevap_sentinel lets us pre-seed qv2qc_condevap to probe reset
  // behavior when the condensation/evaporation branch is not taken.
  //-----------------------------------------------------------------------
  Part2Result run_part2(bool p3_super_sat, bool shoc_enable_condensation,
                        Real T_atm_in, Real pres_in, Real qv_in, Real qc_in, Real nc_in,
                        Real omega_mp_in, Real tke_mp_in, Real condevap_sentinel = 0)
  {
    constexpr Scalar g   = C::gravit.value;
    const Int nk_pack     = 1;
    const Int nk          = Pack::n;

    const Scalar dt     = 1800.0;
    const Scalar inv_dt = 1.0 / dt;
    constexpr Scalar max_total_ni = 740.0e3;

    typename Functions::P3Runtime runtime_options;
    runtime_options.p3_super_sat             = p3_super_sat;
    runtime_options.shoc_enable_condensation  = shoc_enable_condensation;
    runtime_options.do_ice_production         = false; // isolate from ice-phase code paths

    auto hetfrz_immersion  = make_view(0.0, nk_pack);
    auto hetfrz_contact    = make_view(0.0, nk_pack);
    auto hetfrz_deposition = make_view(0.0, nk_pack);

    const auto& tables = Base::lookup_tables;
    const auto dnu                = tables.dnu_table_vals;
    const auto ice_table_vals     = tables.ice_table_vals;
    const auto collect_table_vals = tables.collect_table_vals;
    const auto revap_table_vals   = tables.revap_table_vals;

    auto pres            = make_view(pres_in, nk_pack);
    auto dpres           = make_view(1000.0, nk_pack);
    auto dz              = make_view(1000.0 / (1.0 * g), nk_pack); // -> rho ~ 1.0
    auto nc_nuceat_tend  = make_view(0.0, nk_pack);
    auto inv_exner       = make_view(1.0, nk_pack);
    auto exner           = make_view(1.0, nk_pack);
    auto inv_cld_frac_l  = make_view(1.0, nk_pack); // fully cloudy column
    auto inv_cld_frac_i  = make_view(1.0, nk_pack);
    auto inv_cld_frac_r  = make_view(1.0, nk_pack);
    auto ni_activated    = make_view(0.0, nk_pack);
    auto inv_qc_relvar   = make_view(1.0, nk_pack);
    auto cld_frac_i      = make_view(1.0, nk_pack);
    auto cld_frac_l      = make_view(1.0, nk_pack);
    auto cld_frac_r      = make_view(1.0, nk_pack);
    auto qv_prev         = make_view(qv_in, nk_pack);
    auto t_prev          = make_view(T_atm_in, nk_pack);
    auto omega_mp        = make_view(omega_mp_in, nk_pack);
    auto tke_mp          = make_view(tke_mp_in, nk_pack);
    auto qv2qc_condevap  = make_view(condevap_sentinel, nk_pack);

    auto T_atm     = make_view(T_atm_in, nk_pack);
    // Mirror p3_main_part1's atmospheric-state pre-computation (rho,
    // inv_rho, rhofacr, rhofaci, acn, qv_sat_l) since part2 only READS
    // these fields; they are normally populated by part1 on the same
    // column earlier in the same timestep.
    constexpr Scalar rho_1000mb = C::RHO_1000MB.value;
    constexpr Scalar rho_600mb  = C::RHO_600MB;
    constexpr Scalar rho_h2o    = C::RHO_H2O.value;
    const Real rho_val     = 1.0; // matches dpres=1000, dz=1000/g used below
    const Real inv_rho_val = 1.0 / rho_val;
    const Real rhofacr_val = std::pow(rho_1000mb * inv_rho_val, 0.54);
    const Real rhofaci_val = std::pow(rho_600mb  * inv_rho_val, 0.54);
    const Real mu_dum_val  = 1.496e-6 * std::pow(T_atm_in, 1.5) / (T_atm_in + 120.0);
    const Real acn_val     = g * rho_h2o / (18.0 * mu_dum_val);
    auto rho       = make_view(rho_val, nk_pack);
    auto inv_rho   = make_view(inv_rho_val, nk_pack);
    // p3_main_part2 only ever READS qv_sat_l/qv_sat_i (they are computed and
    // populated earlier, in p3_main_part1, on the same column). Since we
    // call part2 standalone here, replicate that pre-computation ourselves
    // so ssatqv_l = qv - qv_sat_l reflects the actual saturation state.
    const Mask context(true);
    const Real qv_sat_l_val = physics::qv_sat_dry(Pack(T_atm_in), Pack(pres_in), false, context)[0];
    auto qv_sat_l  = make_view(qv_sat_l_val, nk_pack);
    auto qv_sat_i  = make_view(0.0, nk_pack);
    auto qv_supersat_i = make_view(0.0, nk_pack);
    auto rhofacr   = make_view(rhofacr_val, nk_pack);
    auto rhofaci   = make_view(rhofaci_val, nk_pack);
    auto acn       = make_view(acn_val, nk_pack);
    auto qv        = make_view(qv_in, nk_pack);
    auto th_atm    = make_view(300.0, nk_pack);
    auto qc        = make_view(qc_in, nk_pack);
    auto nc        = make_view(nc_in, nk_pack);
    auto qr        = make_view(0.0, nk_pack);
    auto nr        = make_view(0.0, nk_pack);
    auto qi        = make_view(0.0, nk_pack);
    auto ni        = make_view(0.0, nk_pack);
    auto qm        = make_view(0.0, nk_pack);
    auto bm        = make_view(0.0, nk_pack);
    // Fully cloudy column (cld_frac_l = 1) => in-cloud == grid-mean.
    auto qc_incld  = make_view(qc_in, nk_pack);
    auto qr_incld  = make_view(0.0, nk_pack);
    auto qi_incld  = make_view(0.0, nk_pack);
    auto qm_incld  = make_view(0.0, nk_pack);
    auto nc_incld  = make_view(nc_in, nk_pack);
    auto nr_incld  = make_view(0.0, nk_pack);
    auto ni_incld  = make_view(0.0, nk_pack);
    auto bm_incld  = make_view(0.0, nk_pack);

    auto mu_c    = make_view(1.0, nk_pack);
    auto nu      = make_view(1.0, nk_pack);
    auto lamc    = make_view(1.0e5, nk_pack);
    auto cdist   = make_view(0.0, nk_pack);
    auto cdist1  = make_view(0.0, nk_pack);
    auto cdistr  = make_view(0.0, nk_pack);
    auto mu_r    = make_view(1.0, nk_pack);
    auto lamr    = make_view(1.0e4, nk_pack);
    auto logn0r  = make_view(0.0, nk_pack);
    auto qv2qi_depos_tend    = make_view(0.0, nk_pack);
    auto precip_total_tend  = make_view(0.0, nk_pack);
    auto nevapr             = make_view(0.0, nk_pack);
    auto qr_evap_tend       = make_view(0.0, nk_pack);
    auto vap_liq_exchange   = make_view(0.0, nk_pack);
    auto vap_ice_exchange   = make_view(0.0, nk_pack);
    auto liq_ice_exchange   = make_view(0.0, nk_pack);
    auto qr2qv_evap         = make_view(0.0, nk_pack);
    auto qi2qv_sublim       = make_view(0.0, nk_pack);
    auto qc2qr_accret       = make_view(0.0, nk_pack);
    auto qc2qr_autoconv     = make_view(0.0, nk_pack);
    auto qv2qi_vapdep       = make_view(0.0, nk_pack);
    auto qc2qi_berg         = make_view(0.0, nk_pack);
    auto qc2qr_ice_shed     = make_view(0.0, nk_pack);
    auto qc2qi_collect      = make_view(0.0, nk_pack);
    auto qr2qi_collect      = make_view(0.0, nk_pack);
    auto qc2qi_hetero_freeze = make_view(0.0, nk_pack);
    auto qr2qi_immers_freeze = make_view(0.0, nk_pack);
    auto qi2qr_melt         = make_view(0.0, nk_pack);
    auto pratot             = make_view(0.0, nk_pack);
    auto prctot             = make_view(0.0, nk_pack);

    view_1d<bool> bools_d("bools", 1);

    const bool do_predict_nc     = false;
    const bool do_prescribed_CCN = false;

    auto policy = TPF::get_default_team_policy(1, nk_pack);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      Functions::p3_main_part2(
        team, nk_pack, max_total_ni, do_predict_nc, do_prescribed_CCN, dt, inv_dt,
        hetfrz_immersion, hetfrz_contact, hetfrz_deposition,
        dnu, ice_table_vals, collect_table_vals, revap_table_vals,
        pres, dpres, dz, nc_nuceat_tend, inv_exner, exner, inv_cld_frac_l, inv_cld_frac_i, inv_cld_frac_r,
        ni_activated, inv_qc_relvar, cld_frac_i, cld_frac_l, cld_frac_r, qv_prev,
        t_prev, omega_mp, tke_mp, qv2qc_condevap, T_atm, rho,
        inv_rho, qv_sat_l, qv_sat_i, qv_supersat_i, rhofacr, rhofaci, acn,
        qv, th_atm, qc, nc, qr, nr, qi, ni, qm, bm,
        qc_incld, qr_incld, qi_incld, qm_incld, nc_incld, nr_incld, ni_incld, bm_incld,
        mu_c, nu, lamc, cdist, cdist1, cdistr, mu_r, lamr, logn0r,
        qv2qi_depos_tend, precip_total_tend, nevapr, qr_evap_tend,
        vap_liq_exchange, vap_ice_exchange, liq_ice_exchange,
        qr2qv_evap, qi2qv_sublim, qc2qr_accret, qc2qr_autoconv,
        qv2qi_vapdep, qc2qi_berg, qc2qr_ice_shed, qc2qi_collect, qr2qi_collect,
        qc2qi_hetero_freeze, qr2qi_immers_freeze, qi2qr_melt,
        pratot, prctot, bools_d(0), nk, runtime_options);
    });

    Part2Result result;
    auto get0 = [](const view_1d<Pack>& v) {
      auto vh = Kokkos::create_mirror_view(v);
      Kokkos::deep_copy(vh, v);
      return vh(0)[0];
    };
    result.qv             = get0(qv);
    result.qc             = get0(qc);
    result.th_atm         = get0(th_atm);
    result.qv2qc_condevap = get0(qv2qc_condevap);

    auto bools_h = Kokkos::create_mirror_view(bools_d);
    Kokkos::deep_copy(bools_h, bools_d);
    result.is_hydromet_present = bools_h(0);

    return result;
  }

  //=======================================================================
  // p3_main_part1: aerosol activation ("[shanyp 20260816" qv2qc_aeroact)
  //=======================================================================
  void run_property_part1()
  {
    constexpr Scalar latvap = C::LatVap.value;
    constexpr Scalar inv_cp = C::INV_CP.value;
    const Scalar dt = 1800.0;

    const Real T_atm_in = 285.0;
    const Real pres_in  = 90000.0;
    const Real qc_in    = 1.0e-4; // well above qsmall, avoids the "drymass" clip
    const Real nc_in    = 1.0e7;
    const Real nccn_prescribed_in = 5.0e8; // forces nc to increase upon activation

    // Determine qv_sat_l for the chosen T_atm/pres so we can construct
    // clearly super- and sub-saturated qv values (qv_sat_l is recomputed by
    // p3_main_part1 itself from T_atm/pres, so this must match exactly).
    const Mask context(true);
    Pack qsl = physics::qv_sat_dry(Pack(T_atm_in), Pack(pres_in), false, context);
    const Real qv_sat_l_val = qsl[0];

    const Real qv_supersat = qv_sat_l_val * 1.10; // 10% supersaturated
    const Real qv_subsat    = qv_sat_l_val * 0.90; // 10% subsaturated

    // --- Feature completely off: new diagnostic must be exactly zero ---
    {
      auto off = run_part1(/*p3_super_sat=*/false, /*shoc_enable_condensation=*/false,
                            T_atm_in, pres_in, qv_supersat, qc_in, nc_in, nccn_prescribed_in);
      REQUIRE(off.qv2qc_aeroact == 0.0);
    }

    // --- Feature on, but shoc_enable_condensation=true disables it (same
    //     gating "!shoc_enable_condensation && p3_super_sat" as the
    //     condensation/evaporation block) ---
    {
      auto shoc_on = run_part1(/*p3_super_sat=*/true, /*shoc_enable_condensation=*/true,
                                T_atm_in, pres_in, qv_supersat, qc_in, nc_in, nccn_prescribed_in);
      REQUIRE(shoc_on.qv2qc_aeroact == 0.0);
    }

    // --- Feature on + sub-saturated: activation should not trigger, and
    //     the resulting state must match the "feature off" run exactly
    //     (do_predict_nc=do_prescribed_CCN=false makes the rest of
    //     p3_main_part1 independent of p3_super_sat). ---
    {
      auto off        = run_part1(false, false, T_atm_in, pres_in, qv_subsat, qc_in, nc_in, nccn_prescribed_in);
      auto on_subsat   = run_part1(true,  false, T_atm_in, pres_in, qv_subsat, qc_in, nc_in, nccn_prescribed_in);
      REQUIRE(on_subsat.qv2qc_aeroact == 0.0);
      REQUIRE(on_subsat.qc     == off.qc);
      REQUIRE(on_subsat.qv     == off.qv);
      REQUIRE(on_subsat.nc     == off.nc);
      REQUIRE(on_subsat.th_atm == off.th_atm);
    }

    // --- Feature on + super-saturated: activation should trigger. Compare
    //     against the "off" run (same inputs) to isolate exactly the delta
    //     introduced by the new block. ---
    {
      auto off = run_part1(false, false, T_atm_in, pres_in, qv_supersat, qc_in, nc_in, nccn_prescribed_in);
      auto on  = run_part1(true,  false, T_atm_in, pres_in, qv_supersat, qc_in, nc_in, nccn_prescribed_in);

      REQUIRE(off.qv2qc_aeroact == 0.0);
      REQUIRE(on.qv2qc_aeroact > 0.0); // supersaturated + nccn_prescribed >> nc => nc increases => qc gained

      const Real delta_qc = on.qc - off.qc;
      const Real delta_qv = on.qv - off.qv;
      const Real delta_th = on.th_atm - off.th_atm;

      // qv2qc_aeroact(k) = deltaqc/dt, and qc += deltaqc while qv -= deltaqc
      REQUIRE(delta_qc == Approx(on.qv2qc_aeroact * dt).epsilon(1e-10));
      // total water (qv+qc) is conserved by the activation step
      REQUIRE(delta_qv == Approx(-delta_qc).epsilon(1e-10));
      // th_atm change is consistent with condensational latent heating
      // (exner = 1 in run_part1's inputs, so th_atm delta reduces to
      // latvap*inv_cp*delta_qc)
      REQUIRE(delta_th == Approx(latvap * inv_cp * delta_qc).epsilon(1e-8));

      // activation can only add cloud water/number, never remove it
      REQUIRE(on.nc >= off.nc);
      REQUIRE(on.qc >= off.qc);
    }
  }

  //=======================================================================
  // p3_main_part2: condensation/evaporation ("[shanyp 20260814/20260816"
  // qv2qc_condevap)
  //=======================================================================
  void run_property_part2()
  {
    const Scalar dt = 1800.0;
    constexpr Scalar latvap = C::LatVap.value;
    constexpr Scalar inv_cp = C::INV_CP.value;

    const Real T_atm_in = 285.0;
    const Real pres_in  = 90000.0;
    const Real qc_in    = 1.0e-4; // well above qsmall so condmask is satisfied
    const Real nc_in    = 1.0e7;
    // omega_mp = tke_mp = 0 zeroes out the vertical-velocity term (w0d) in
    // the condensation formula, which otherwise depends on lookup-table
    // derived quantities (epsc, ab) that are impractical to hand-derive here.
    const Real omega_mp_in = 0.0;
    const Real tke_mp_in   = 0.0;

    const Mask context(true);
    Pack qsl = physics::qv_sat_dry(Pack(T_atm_in), Pack(pres_in), false, context);
    const Real qv_sat_l_val = qsl[0];

    const Real qv_supersat = qv_sat_l_val * 1.10;
    const Real qv_subsat    = qv_sat_l_val * 0.90;

    // --- Feature off: diagnostic must be exactly zero ---
    {
      auto off = run_part2(false, false, T_atm_in, pres_in, qv_supersat, qc_in, nc_in, omega_mp_in, tke_mp_in);
      REQUIRE(off.qv2qc_condevap == 0.0);
    }

    // --- shoc_enable_condensation=true disables the new block even with
    //     p3_super_sat=true ---
    {
      auto shoc_on = run_part2(true, true, T_atm_in, pres_in, qv_supersat, qc_in, nc_in, omega_mp_in, tke_mp_in);
      REQUIRE(shoc_on.qv2qc_condevap == 0.0);
    }

    // --- Supersaturated => condensation (qv2qc_condevap > 0). qv is
    //     decremented by exactly qv2qc_condevap*dt directly in the new
    //     block (a reliable check); qc, however, is *not* guaranteed to
    //     increase by that same amount, because later in the same
    //     p3_main_part2 call, warm-rain autoconversion/accretion
    //     (update_prognostic_liquid) also consumes qc (converting some of
    //     it to qr) using the *post-condensation* qc/nc state -- so the
    //     "on" run autoconverts more mass to rain than the "off" run does.
    //     We therefore only check the qv side exactly, and check qc/qr
    //     directionally. ---
    {
      auto on = run_part2(true, false, T_atm_in, pres_in, qv_supersat, qc_in, nc_in, omega_mp_in, tke_mp_in);
      REQUIRE(on.qv2qc_condevap > 0.0);

      auto off = run_part2(false, false, T_atm_in, pres_in, qv_supersat, qc_in, nc_in, omega_mp_in, tke_mp_in);
      const Real delta_qv = on.qv - off.qv;

      // The qv decrement happens directly in the new block, unaffected by
      // any downstream physics, so this holds exactly.
      REQUIRE(delta_qv == Approx(-(on.qv2qc_condevap * dt)).epsilon(1e-8));
      // Net effect: turning on condensation can only add (not remove)
      // cloud water/number relative to the "off" run, even after
      // downstream autoconversion/accretion partially consumes it.
      REQUIRE(on.qc >= off.qc);

      // NOTE: this exercises the fix for a previously-uninitialized local
      // "pcc" variable in p3_main_impl_part2.hpp's condensation block,
      // which used to make this th_atm update a permanent no-op. th_atm
      // is now updated using qv2qc_condevap (the same quantity used for
      // qc/qv), so the latent-heat release is exact given exner=1,
      // cld_frac_l=1 in this setup.
      const Real delta_th_atm = on.th_atm - off.th_atm;
      REQUIRE(delta_th_atm == Approx(latvap * inv_cp * on.qv2qc_condevap * dt).epsilon(1e-8));
    }

    // --- Subsaturated: qc_incld >= qsmall still makes condmask true (it
    //     only checks qc_incld against qsmall, not saturation state), so
    //     the branch still executes but should now evaporate cloud water
    //     (qv2qc_condevap <= 0). As above, only the qv side is checked
    //     exactly; qc is also affected by downstream autoconversion using
    //     the post-evaporation qc, so we only check its direction. ---
    {
      auto on = run_part2(true, false, T_atm_in, pres_in, qv_subsat, qc_in, nc_in, omega_mp_in, tke_mp_in);
      REQUIRE(on.qv2qc_condevap <= 0.0);

      auto off = run_part2(false, false, T_atm_in, pres_in, qv_subsat, qc_in, nc_in, omega_mp_in, tke_mp_in);
      const Real delta_qv = on.qv - off.qv;
      REQUIRE(delta_qv == Approx(-(on.qv2qc_condevap * dt)).epsilon(1e-8));
      REQUIRE(on.qc <= off.qc);
    }
  }

  //=======================================================================
  // Diagnostic-array reset behavior. qv2qc_aeroact/qv2qc_condevap are only
  // ever written via masked .set() calls guarded by the activation /
  // condensation masks; neither p3_main_part1 nor p3_main_part2 explicitly
  // zeroes them out elsewhere. This test documents that, with the feature
  // off (or its trigger condition unmet), a caller-provided non-zero
  // "sentinel" value currently survives the call unchanged -- i.e. these
  // diagnostics are NOT self-resetting and must be zeroed by the caller
  // every timestep before invoking p3_main_part1/p3_main_part2.
  //=======================================================================
  void run_property_diagnostic_reset()
  {
    const Real T_atm_in = 285.0;
    const Real pres_in  = 90000.0;
    const Real qc_in    = 1.0e-4;
    const Real nc_in    = 1.0e7;
    const Real sentinel = 12345.0;

    const Mask context(true);
    Pack qsl = physics::qv_sat_dry(Pack(T_atm_in), Pack(pres_in), false, context);
    const Real qv_subsat = qsl[0] * 0.90;

    // p3_super_sat = false => qv2qc_aeroact block is skipped entirely.
    auto p1 = run_part1(false, false, T_atm_in, pres_in, qv_subsat, qc_in, nc_in, 5.0e8, sentinel);
    REQUIRE(p1.qv2qc_aeroact == sentinel);

    auto p2 = run_part2(false, false, T_atm_in, pres_in, qv_subsat, qc_in, nc_in, 0.0, 0.0, sentinel);
    REQUIRE(p2.qv2qc_condevap == sentinel);
  }

  void run_property()
  {
    run_property_part1();
    run_property_part2();
    run_property_diagnostic_reset();
  }
};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("p3_prognostic_supersat_part1_activation", "[p3]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPrognosticSupersat;

  T t;
  t.run_property_part1();
}

TEST_CASE("p3_prognostic_supersat_part2_condevap", "[p3]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPrognosticSupersat;

  T t;
  t.run_property_part2();
}

TEST_CASE("p3_prognostic_supersat_diagnostic_reset", "[p3]")
{
  using T = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPrognosticSupersat;

  T t;
  t.run_property_diagnostic_reset();
}

} // empty namespace
//unittest]
