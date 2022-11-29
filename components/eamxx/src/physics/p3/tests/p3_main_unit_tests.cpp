#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestP3Main {

static void run_phys_p3_main_part1()
{
  // TODO
}

static void run_phys_p3_main_part2()
{
  // TODO
}

static void run_phys_p3_main_part3()
{
  // TODO
}

static void run_phys_p3_main()
{
  // TODO
}

static void run_phys()
{
  run_phys_p3_main_part1();
  run_phys_p3_main_part2();
  run_phys_p3_main_part3();
  run_phys_p3_main();
}

static void run_bfb_p3_main_part1()
{
  auto engine = setup_random_test();

  constexpr Scalar qsmall = C::QSMALL; //PMC wouldn't it make more sense to define qsmall at a higher level since used in part1, part2, and part3?
  constexpr Scalar T_zerodegc   = C::T_zerodegc;
  constexpr Scalar sup_upper = -0.05;
  constexpr Scalar sup_lower = -0.1;

  P3MainPart1Data isds_fortran[] = {
    //            kts, kte, ktop, kbot, kdir, do_predict_nc, do_prescribed_CCN,       dt
    P3MainPart1Data(1,  72,    1,   72,    1, false,          true,           1.800E+03),
    P3MainPart1Data(1,  72,    1,   72,    1, true,           true,           1.800E+03),
    P3MainPart1Data(1,  72,   72,    1,   -1, false,          false,          1.800E+03),
    P3MainPart1Data(1,  72,   72,    1,   -1, true,           false,          1.800E+03),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(P3MainPart1Data);

  for (auto& d : isds_fortran) {
    const auto qsmall_r = std::make_pair(0, qsmall*2); //PMC this range seems inappropriately small
    d.randomize(engine, {
        {d.T_atm, {T_zerodegc - 10, T_zerodegc + 10}},
        {d.qv_supersat_i, {sup_lower -.05, sup_upper + .05}},
        {d.qc, qsmall_r}, {d.qr, qsmall_r}, {d.qi, qsmall_r} });
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  P3MainPart1Data isds_cxx[num_runs] = {
    P3MainPart1Data(isds_fortran[0]),
    P3MainPart1Data(isds_fortran[1]),
    P3MainPart1Data(isds_fortran[2]),
    P3MainPart1Data(isds_fortran[3]),
  };

    // Get data from fortran
  for (auto& d : isds_fortran) {
    p3_main_part1(d);
  }

  // Get data from cxx
  for (auto& d : isds_cxx) {
    p3_main_part1_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir, d.do_predict_nc, d.do_prescribed_CCN, d.dt,
                    d.pres, d.dpres, d.dz, d.nc_nuceat_tend, d.nccn_prescribed, d.inv_exner, d.exner, d.inv_cld_frac_l, d.inv_cld_frac_i,
                    d.inv_cld_frac_r, d.latent_heat_vapor, d.latent_heat_sublim, d.latent_heat_fusion,
                    d.T_atm, d.rho, d.inv_rho, d.qv_sat_l, d.qv_sat_i, d.qv_supersat_i, d.rhofacr, d.rhofaci,
                    d.acn, d.qv, d.th_atm, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni, d.qm, d.bm, d.qc_incld, d.qr_incld, d.qi_incld,
                    d.qm_incld, d.nc_incld, d.nr_incld, d.ni_incld, d.bm_incld,
                    &d.is_nucleat_possible, &d.is_hydromet_present);
  }

  if (SCREAM_BFB_TESTING) {
    for (Int i = 0; i < num_runs; ++i) {
      Int start = std::min(isds_fortran[i].kbot, isds_fortran[i].ktop) - 1; // 0-based indx
      Int end   = std::max(isds_fortran[i].kbot, isds_fortran[i].ktop);     // 0-based indx
      for (Int k = start; k < end; ++k) {
        REQUIRE(isds_fortran[i].T_atm[k]         == isds_cxx[i].T_atm[k]);
        REQUIRE(isds_fortran[i].rho[k]           == isds_cxx[i].rho[k]);
        REQUIRE(isds_fortran[i].inv_rho[k]       == isds_cxx[i].inv_rho[k]);
        REQUIRE(isds_fortran[i].qv_sat_l[k]      == isds_cxx[i].qv_sat_l[k]);
        REQUIRE(isds_fortran[i].qv_sat_i[k]      == isds_cxx[i].qv_sat_i[k]);
        REQUIRE(isds_fortran[i].qv_supersat_i[k] == isds_cxx[i].qv_supersat_i[k]);
        REQUIRE(isds_fortran[i].rhofacr[k]       == isds_cxx[i].rhofacr[k]);
        REQUIRE(isds_fortran[i].rhofaci[k]       == isds_cxx[i].rhofaci[k]);
        REQUIRE(isds_fortran[i].acn[k]           == isds_cxx[i].acn[k]);
        REQUIRE(isds_fortran[i].qv[k]            == isds_cxx[i].qv[k]);
        REQUIRE(isds_fortran[i].th_atm[k]        == isds_cxx[i].th_atm[k]);
        REQUIRE(isds_fortran[i].qc[k]            == isds_cxx[i].qc[k]);
        REQUIRE(isds_fortran[i].nc[k]            == isds_cxx[i].nc[k]);
        REQUIRE(isds_fortran[i].qr[k]            == isds_cxx[i].qr[k]);
        REQUIRE(isds_fortran[i].nr[k]            == isds_cxx[i].nr[k]);
        REQUIRE(isds_fortran[i].qi[k]            == isds_cxx[i].qi[k]);
        REQUIRE(isds_fortran[i].ni[k]            == isds_cxx[i].ni[k]);
        REQUIRE(isds_fortran[i].qm[k]            == isds_cxx[i].qm[k]);
        REQUIRE(isds_fortran[i].bm[k]            == isds_cxx[i].bm[k]);
        REQUIRE(isds_fortran[i].qc_incld[k]      == isds_cxx[i].qc_incld[k]);
        REQUIRE(isds_fortran[i].qr_incld[k]      == isds_cxx[i].qr_incld[k]);
        REQUIRE(isds_fortran[i].qi_incld[k]      == isds_cxx[i].qi_incld[k]);
        REQUIRE(isds_fortran[i].qm_incld[k]      == isds_cxx[i].qm_incld[k]);
        REQUIRE(isds_fortran[i].nc_incld[k]      == isds_cxx[i].nc_incld[k]);
        REQUIRE(isds_fortran[i].nr_incld[k]      == isds_cxx[i].nr_incld[k]);
        REQUIRE(isds_fortran[i].ni_incld[k]      == isds_cxx[i].ni_incld[k]);
        REQUIRE(isds_fortran[i].bm_incld[k]      == isds_cxx[i].bm_incld[k]);
      }
      REQUIRE( isds_fortran[i].is_hydromet_present == isds_cxx[i].is_hydromet_present );
      REQUIRE( isds_fortran[i].is_nucleat_possible == isds_cxx[i].is_nucleat_possible );
    }
  }
}

static void run_bfb_p3_main_part2()
{
  auto engine = setup_random_test();

  constexpr Scalar qsmall     = C::QSMALL;
  constexpr Scalar T_zerodegc   = C::T_zerodegc;
  constexpr Scalar sup_upper = -0.05;
  constexpr Scalar sup_lower = -0.1;

  P3MainPart2Data isds_fortran[] = {
    //            kts, kte, ktop, kbot, kdir, do_predict_nc, do_prescribed_CCN,       dt
    P3MainPart2Data(1,  72,    1,   72,    1, false,         true,        1.800E+03),
    P3MainPart2Data(1,  72,    1,   72,    1, true,          true,        1.800E+03),
    P3MainPart2Data(1,  72,   72,    1,   -1, false,         false,       1.800E+03),
    P3MainPart2Data(1,  72,   72,    1,   -1, true,          false,       1.800E+03),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(P3MainPart2Data);

  for (auto& d : isds_fortran) {
    const auto qsmall_r = std::make_pair(0, qsmall*2);
    d.randomize(engine, {
        {d.T_atm,  {T_zerodegc - 10, T_zerodegc + 10}},
        {d.t_prev, {T_zerodegc - 10, T_zerodegc + 10}},
        {d.qv_supersat_i, {sup_lower -.05, sup_upper + .05}},
        {d.qc, qsmall_r}, {d.qr, qsmall_r}, {d.qi, qsmall_r} });
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  P3MainPart2Data isds_cxx[num_runs] = {
    P3MainPart2Data(isds_fortran[0]),
    P3MainPart2Data(isds_fortran[1]),
    P3MainPart2Data(isds_fortran[2]),
    P3MainPart2Data(isds_fortran[3]),
  };

  // Get data from fortran
  for (auto& d : isds_fortran) {
    p3_main_part2(d);
  }

  // Get data from cxx
  for (auto& d : isds_cxx) {
    p3_main_part2_f(
      d.kts, d.kte, d.kbot, d.ktop, d.kdir, d.do_predict_nc, d.do_prescribed_CCN, d.dt, d.inv_dt,
      d.pres, d.dpres, d.dz, d.nc_nuceat_tend, d.inv_exner, d.exner, d.inv_cld_frac_l, d.inv_cld_frac_i,
      d.inv_cld_frac_r, d.ni_activated, d.inv_qc_relvar, d.cld_frac_i, d.cld_frac_l, d.cld_frac_r, d.qv_prev, d.t_prev,
      d.T_atm, d.rho, d.inv_rho, d.qv_sat_l, d.qv_sat_i, d.qv_supersat_i, d.rhofacr, d.rhofaci, d.acn, d.qv, d.th_atm, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni,
      d.qm, d.bm, d.latent_heat_vapor, d.latent_heat_sublim, d.latent_heat_fusion, d.qc_incld, d.qr_incld, d.qi_incld, d.qm_incld, d.nc_incld, d.nr_incld,
      d.ni_incld, d.bm_incld, d.mu_c, d.nu, d.lamc, d.cdist, d.cdist1, d.cdistr, d.mu_r, d.lamr, d.logn0r, d.qv2qi_depos_tend, d.precip_total_tend,
      d.nevapr, d.qr_evap_tend, d.vap_liq_exchange, d.vap_ice_exchange, d.liq_ice_exchange, d.pratot,
      d.prctot, &d.is_hydromet_present);
  }

  if (SCREAM_BFB_TESTING) {
    for (Int i = 0; i < num_runs; ++i) {
      Int start = std::min(isds_fortran[i].kbot, isds_fortran[i].ktop) - 1; // 0-based indx
      Int end   = std::max(isds_fortran[i].kbot, isds_fortran[i].ktop);     // 0-based indx
      for (Int k = start; k < end; ++k) {
        REQUIRE(isds_fortran[i].T_atm[k]              == isds_cxx[i].T_atm[k]);
        REQUIRE(isds_fortran[i].rho[k]                == isds_cxx[i].rho[k]);
        REQUIRE(isds_fortran[i].inv_rho[k]            == isds_cxx[i].inv_rho[k]);
        REQUIRE(isds_fortran[i].qv_sat_l[k]           == isds_cxx[i].qv_sat_l[k]);
        REQUIRE(isds_fortran[i].qv_sat_i[k]           == isds_cxx[i].qv_sat_i[k]);
        REQUIRE(isds_fortran[i].qv_supersat_i[k]      == isds_cxx[i].qv_supersat_i[k]);
        REQUIRE(isds_fortran[i].rhofacr[k]            == isds_cxx[i].rhofacr[k]);
        REQUIRE(isds_fortran[i].rhofaci[k]            == isds_cxx[i].rhofaci[k]);
        REQUIRE(isds_fortran[i].acn[k]                == isds_cxx[i].acn[k]);
        REQUIRE(isds_fortran[i].qv[k]                 == isds_cxx[i].qv[k]);
        REQUIRE(isds_fortran[i].th_atm[k]             == isds_cxx[i].th_atm[k]);
        REQUIRE(isds_fortran[i].qc[k]                 == isds_cxx[i].qc[k]);
        REQUIRE(isds_fortran[i].nc[k]                 == isds_cxx[i].nc[k]);
        REQUIRE(isds_fortran[i].qr[k]                 == isds_cxx[i].qr[k]);
        REQUIRE(isds_fortran[i].nr[k]                 == isds_cxx[i].nr[k]);
        REQUIRE(isds_fortran[i].qi[k]                 == isds_cxx[i].qi[k]);
        REQUIRE(isds_fortran[i].ni[k]                 == isds_cxx[i].ni[k]);
        REQUIRE(isds_fortran[i].qm[k]                 == isds_cxx[i].qm[k]);
        REQUIRE(isds_fortran[i].bm[k]                 == isds_cxx[i].bm[k]);
        REQUIRE(isds_fortran[i].latent_heat_vapor[k]  == isds_cxx[i].latent_heat_vapor[k]);
        REQUIRE(isds_fortran[i].latent_heat_sublim[k] == isds_cxx[i].latent_heat_sublim[k]);
        REQUIRE(isds_fortran[i].latent_heat_fusion[k] == isds_cxx[i].latent_heat_fusion[k]);
        REQUIRE(isds_fortran[i].qc_incld[k]           == isds_cxx[i].qc_incld[k]);
        REQUIRE(isds_fortran[i].qr_incld[k]           == isds_cxx[i].qr_incld[k]);
        REQUIRE(isds_fortran[i].qi_incld[k]           == isds_cxx[i].qi_incld[k]);
        REQUIRE(isds_fortran[i].qm_incld[k]           == isds_cxx[i].qm_incld[k]);
        REQUIRE(isds_fortran[i].nc_incld[k]           == isds_cxx[i].nc_incld[k]);
        REQUIRE(isds_fortran[i].nr_incld[k]           == isds_cxx[i].nr_incld[k]);
        REQUIRE(isds_fortran[i].ni_incld[k]           == isds_cxx[i].ni_incld[k]);
        REQUIRE(isds_fortran[i].bm_incld[k]           == isds_cxx[i].bm_incld[k]);
        REQUIRE(isds_fortran[i].mu_c[k]               == isds_cxx[i].mu_c[k]);
        REQUIRE(isds_fortran[i].nu[k]                 == isds_cxx[i].nu[k]);
        REQUIRE(isds_fortran[i].lamc[k]               == isds_cxx[i].lamc[k]);
        REQUIRE(isds_fortran[i].cdist[k]              == isds_cxx[i].cdist[k]);
        REQUIRE(isds_fortran[i].cdist1[k]             == isds_cxx[i].cdist1[k]);
        REQUIRE(isds_fortran[i].cdistr[k]             == isds_cxx[i].cdistr[k]);
        REQUIRE(isds_fortran[i].mu_r[k]               == isds_cxx[i].mu_r[k]);
        REQUIRE(isds_fortran[i].lamr[k]               == isds_cxx[i].lamr[k]);
        REQUIRE(isds_fortran[i].logn0r[k]             == isds_cxx[i].logn0r[k]);
        REQUIRE(isds_fortran[i].qv2qi_depos_tend[k]            == isds_cxx[i].qv2qi_depos_tend[k]);
        REQUIRE(isds_fortran[i].precip_total_tend[k]  == isds_cxx[i].precip_total_tend[k]);
        REQUIRE(isds_fortran[i].nevapr[k]             == isds_cxx[i].nevapr[k]);
        REQUIRE(isds_fortran[i].qr_evap_tend[k]       == isds_cxx[i].qr_evap_tend[k]);
        REQUIRE(isds_fortran[i].vap_liq_exchange[k]   == isds_cxx[i].vap_liq_exchange[k]);
        REQUIRE(isds_fortran[i].vap_ice_exchange[k]   == isds_cxx[i].vap_ice_exchange[k]);
        REQUIRE(isds_fortran[i].liq_ice_exchange[k]   == isds_cxx[i].liq_ice_exchange[k]);
        REQUIRE(isds_fortran[i].pratot[k]             == isds_cxx[i].pratot[k]);
        REQUIRE(isds_fortran[i].prctot[k]             == isds_cxx[i].prctot[k]);
      }
      REQUIRE( isds_fortran[i].is_hydromet_present == isds_cxx[i].is_hydromet_present );
    }
  }
}

static void run_bfb_p3_main_part3()
{
  auto engine = setup_random_test();

  constexpr Scalar qsmall     = C::QSMALL;

  P3MainPart3Data isds_fortran[] = {
    //            kts, kte, ktop, kbot, kdir
    P3MainPart3Data(1,  72,    1,   72,    1),
    P3MainPart3Data(1,  72,    1,   72,    1),
    P3MainPart3Data(1,  72,   72,    1,   -1),
    P3MainPart3Data(1,  72,   72,    1,   -1),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(P3MainPart3Data);

  for (auto& d : isds_fortran) {
    const auto qsmall_r = std::make_pair(0, qsmall*2);
    d.randomize(engine, { {d.qc, qsmall_r}, {d.qr, qsmall_r}, {d.qi, qsmall_r} });
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  P3MainPart3Data isds_cxx[num_runs] = {
    P3MainPart3Data(isds_fortran[0]),
    P3MainPart3Data(isds_fortran[1]),
    P3MainPart3Data(isds_fortran[2]),
    P3MainPart3Data(isds_fortran[3]),
  };

  // Get data from fortran
  for (auto& d : isds_fortran) {
    p3_main_part3(d);
  }

  // Get data from cxx
  for (auto& d : isds_cxx) {
    p3_main_part3_f(
      d.kts, d.kte, d.kbot, d.ktop, d.kdir,
      d.inv_exner, d.cld_frac_l, d.cld_frac_r, d.cld_frac_i,
      d.rho, d.inv_rho, d.rhofaci, d.qv, d.th_atm, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni, d.qm, d.bm, d.latent_heat_vapor, d.latent_heat_sublim,
      d.mu_c, d.nu, d.lamc, d.mu_r, d.lamr, d.vap_liq_exchange,
      d. ze_rain, d.ze_ice, d.diag_vm_qi, d.diag_eff_radius_qi, d.diag_diam_qi, d.rho_qi, d.diag_equiv_reflectivity, d.diag_eff_radius_qc);
  }

  if (SCREAM_BFB_TESTING) {
    for (Int i = 0; i < num_runs; ++i) {
      Int start = std::min(isds_fortran[i].kbot, isds_fortran[i].ktop) - 1; // 0-based indx
      Int end   = std::max(isds_fortran[i].kbot, isds_fortran[i].ktop);     // 0-based indx
      for (Int k = start; k < end; ++k) {
        REQUIRE(isds_fortran[i].rho[k]                     == isds_cxx[i].rho[k]);
        REQUIRE(isds_fortran[i].inv_rho[k]                 == isds_cxx[i].inv_rho[k]);
        REQUIRE(isds_fortran[i].rhofaci[k]                 == isds_cxx[i].rhofaci[k]);
        REQUIRE(isds_fortran[i].qv[k]                      == isds_cxx[i].qv[k]);
        REQUIRE(isds_fortran[i].th_atm[k]                  == isds_cxx[i].th_atm[k]);
        REQUIRE(isds_fortran[i].qc[k]                      == isds_cxx[i].qc[k]);
        REQUIRE(isds_fortran[i].nc[k]                      == isds_cxx[i].nc[k]);
        REQUIRE(isds_fortran[i].qr[k]                      == isds_cxx[i].qr[k]);
        REQUIRE(isds_fortran[i].nr[k]                      == isds_cxx[i].nr[k]);
        REQUIRE(isds_fortran[i].qi[k]                      == isds_cxx[i].qi[k]);
        REQUIRE(isds_fortran[i].ni[k]                      == isds_cxx[i].ni[k]);
        REQUIRE(isds_fortran[i].qm[k]                      == isds_cxx[i].qm[k]);
        REQUIRE(isds_fortran[i].bm[k]                      == isds_cxx[i].bm[k]);
        REQUIRE(isds_fortran[i].latent_heat_vapor[k]       == isds_cxx[i].latent_heat_vapor[k]);
        REQUIRE(isds_fortran[i].latent_heat_sublim[k]      == isds_cxx[i].latent_heat_sublim[k]);
        REQUIRE(isds_fortran[i].mu_c[k]                    == isds_cxx[i].mu_c[k]);
        REQUIRE(isds_fortran[i].nu[k]                      == isds_cxx[i].nu[k]);
        REQUIRE(isds_fortran[i].lamc[k]                    == isds_cxx[i].lamc[k]);
        REQUIRE(isds_fortran[i].mu_r[k]                    == isds_cxx[i].mu_r[k]);
        REQUIRE(isds_fortran[i].lamr[k]                    == isds_cxx[i].lamr[k]);
        REQUIRE(isds_fortran[i].vap_liq_exchange[k]        == isds_cxx[i].vap_liq_exchange[k]);
        REQUIRE(isds_fortran[i].ze_rain[k]                 == isds_cxx[i].ze_rain[k]);
        REQUIRE(isds_fortran[i].ze_ice[k]                  == isds_cxx[i].ze_ice[k]);
        REQUIRE(isds_fortran[i].diag_vm_qi[k]              == isds_cxx[i].diag_vm_qi[k]);
        REQUIRE(isds_fortran[i].diag_eff_radius_qi[k]         == isds_cxx[i].diag_eff_radius_qi[k]);
        REQUIRE(isds_fortran[i].diag_diam_qi[k]            == isds_cxx[i].diag_diam_qi[k]);
        REQUIRE(isds_fortran[i].rho_qi[k]                  == isds_cxx[i].rho_qi[k]);
        REQUIRE(isds_fortran[i].diag_equiv_reflectivity[k] == isds_cxx[i].diag_equiv_reflectivity[k]);
        REQUIRE(isds_fortran[i].diag_eff_radius_qc[k]         == isds_cxx[i].diag_eff_radius_qc[k]);
      }
    }
  }
}

static void run_bfb_p3_main()
{
  auto engine = setup_random_test();

  P3MainData isds_fortran[] = {
    //      its, ite, kts, kte,   it,        dt, do_predict_nc, do_prescribed_CCN
    P3MainData(1, 10,   1,  72,    1, 1.800E+03, false, true),
    P3MainData(1, 10,   1,  72,    1, 1.800E+03, true,  false),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(P3MainData);

  for (auto& d : isds_fortran) {
    d.randomize(engine, {
        {d.pres           , {1.00000000E+02 , 9.87111111E+04}},
        {d.dz             , {1.22776609E+02 , 3.49039167E+04}},
        {d.nc_nuceat_tend , {0              , 0}},
        {d.nccn_prescribed, {0              , 0}},
        {d.ni_activated   , {0              , 0}},
        {d.dpres          , {1.37888889E+03, 1.39888889E+03}},
        {d.inv_exner      , {1.00371345E+00, 3.19721007E+00}},
        {d.cld_frac_i     , {1              , 1}},
        {d.cld_frac_l     , {1              , 1}},
        {d.cld_frac_r     , {1              , 1}},
        {d.inv_qc_relvar  , {1              , 1}},
        {d.qc             , {0              , 1.00000000E-04}},
        {d.nc             , {1.00000000E+06 , 1.00000000E+06}},
        {d.qr             , {0              , 1.00000000E-05}},
        {d.nr             , {1.00000000E+06 , 1.00000000E+06}},
        {d.qi             , {0              , 1.00000000E-04}},
        {d.qm             , {0              , 1.00000000E-04}},
        {d.ni             , {1.00000000E+06 , 1.00000000E+06}},
        {d.bm             , {0              , 1.00000000E-02}},
        {d.qv             , {0              , 5.00000000E-02}},
        {d.qv_prev        , {0              , 5.00000000E-02}},
        {d.th_atm         , {6.72653866E+02 , 1.07954335E+03}}, //PMC - this range seems insane
        {d.t_prev         , {1.50000000E+02 , 3.50000000E+02}},
    });
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  P3MainData isds_cxx[num_runs] = {
    P3MainData(isds_fortran[0]),
    P3MainData(isds_fortran[1]),
  };

  // Get data from fortran
  for (auto& d : isds_fortran) {
    p3_main(d);
  }

  // Get data from cxx
  for (auto& d : isds_cxx) {
    d.transpose<ekat::TransposeDirection::c2f>();
    p3_main_f(
      d.qc, d.nc, d.qr, d.nr, d.th_atm, d.qv, d.dt, d.qi, d.qm, d.ni,
      d.bm, d.pres, d.dz, d.nc_nuceat_tend, d.nccn_prescribed, d.ni_activated, d.inv_qc_relvar, d.it, d.precip_liq_surf,
      d.precip_ice_surf, d.its, d.ite, d.kts, d.kte, d.diag_eff_radius_qc, d.diag_eff_radius_qi,
      d.rho_qi, d.do_predict_nc, d.do_prescribed_CCN, d.dpres, d.inv_exner, d.qv2qi_depos_tend,
      d.precip_liq_flux, d.precip_ice_flux, d.cld_frac_r, d.cld_frac_l, d.cld_frac_i, 
      d.liq_ice_exchange, d.vap_liq_exchange, d.vap_ice_exchange, d.qv_prev, d.t_prev);
    d.transpose<ekat::TransposeDirection::f2c>();
  }

  if (SCREAM_BFB_TESTING) {
    for (Int i = 0; i < num_runs; ++i) {
      const auto& df90 = isds_fortran[i];
      const auto& dcxx = isds_fortran[i];
      const auto tot = isds_fortran[i].total(df90.qc);
      for (Int t = 0; t < tot; ++t) {
        REQUIRE(df90.qc[t]                == dcxx.qc[t]);
        REQUIRE(df90.nc[t]                == dcxx.nc[t]);
        REQUIRE(df90.qr[t]                == dcxx.qr[t]);
        REQUIRE(df90.nr[t]                == dcxx.nr[t]);
        REQUIRE(df90.qi[t]                == dcxx.qi[t]);
        REQUIRE(df90.qm[t]                == dcxx.qm[t]);
        REQUIRE(df90.ni[t]                == dcxx.ni[t]);
        REQUIRE(df90.bm[t]                == dcxx.bm[t]);
        REQUIRE(df90.qv[t]                == dcxx.qv[t]);
        REQUIRE(df90.th_atm[t]            == dcxx.th_atm[t]);
        REQUIRE(df90.diag_eff_radius_qc[t]         == dcxx.diag_eff_radius_qc[t]);
        REQUIRE(df90.diag_eff_radius_qi[t]         == dcxx.diag_eff_radius_qi[t]);
        REQUIRE(df90.rho_qi[t]            == dcxx.rho_qi[t]);
        REQUIRE(df90.mu_c[t]              == dcxx.mu_c[t]);
        REQUIRE(df90.lamc[t]              == dcxx.lamc[t]);
        REQUIRE(df90.qv2qi_depos_tend[t]           == dcxx.qv2qi_depos_tend[t]);
        REQUIRE(df90.precip_total_tend[t] == dcxx.precip_total_tend[t]);
        REQUIRE(df90.nevapr[t]            == dcxx.nevapr[t]);
        REQUIRE(df90.qr_evap_tend[t]      == dcxx.qr_evap_tend[t]);
        REQUIRE(df90.liq_ice_exchange[t]  == dcxx.liq_ice_exchange[t]);
        REQUIRE(df90.vap_liq_exchange[t]  == dcxx.vap_liq_exchange[t]);
        REQUIRE(df90.vap_ice_exchange[t]  == dcxx.vap_ice_exchange[t]);
        REQUIRE(df90.precip_liq_flux[t]   == dcxx.precip_liq_flux[t]);
        REQUIRE(df90.precip_ice_flux[t]   == dcxx.precip_ice_flux[t]);
        REQUIRE(df90.precip_liq_surf[t]   == dcxx.precip_liq_surf[t]);
        REQUIRE(df90.precip_ice_surf[t]   == dcxx.precip_ice_surf[t]);
      }
      REQUIRE(df90.precip_liq_flux[tot]   == dcxx.precip_liq_flux[tot]);
      REQUIRE(df90.precip_ice_flux[tot]   == dcxx.precip_ice_flux[tot]);
      REQUIRE(df90.precip_liq_surf[tot]   == dcxx.precip_liq_surf[tot]);
      REQUIRE(df90.precip_ice_surf[tot]   == dcxx.precip_ice_surf[tot]);
    }
  }
}

static void run_bfb()
{
  run_bfb_p3_main_part1();
  run_bfb_p3_main_part2();
  run_bfb_p3_main_part3();
  run_bfb_p3_main();
}

};

}
}
}

namespace {

TEST_CASE("p3_main", "[p3_functions]")
{
  using TP3 = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestP3Main;

  TP3::run_phys();
  TP3::run_bfb();

  scream::p3::P3GlobalForFortran::deinit();
}

} // namespace
