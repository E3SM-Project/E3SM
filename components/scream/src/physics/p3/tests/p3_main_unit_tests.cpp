#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"

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
  constexpr Scalar qsmall     = C::QSMALL;
  constexpr Scalar zerodegc   = C::ZeroDegC;
  constexpr Scalar sup_upper = -0.05;
  constexpr Scalar sup_lower = -0.1;

  const std::array< std::pair<Real, Real>, P3MainPart1Data::NUM_ARRAYS > ranges = {
    std::make_pair(0, 1), // pres
    std::make_pair(0, 1), // dpres
    std::make_pair(0, 1), // dz
    std::make_pair(0, 1), // nc_nuceat_tend
    std::make_pair(0, 1), // exner
    std::make_pair(0, 1), // inv_exner
    std::make_pair(0, 1), // inv_cld_frac_l
    std::make_pair(0, 1), // inv_cld_frac_i
    std::make_pair(0, 1), // inv_cld_frac_r
    std::make_pair(0, 1), // latent_heat_vapor
    std::make_pair(0, 1), // latent_heat_sublim
    std::make_pair(0, 1), // latent_heat_fusion
    std::make_pair(zerodegc - 10, zerodegc + 10), // t
    std::make_pair(0, 1), // rho
    std::make_pair(0, 1), // inv_rho
    std::make_pair(0, 1), // qv_sat_l
    std::make_pair(0, 1), // qv_sat_i
    std::make_pair(sup_lower -.05, sup_upper + .05), // qv_supersat_i
    std::make_pair(0, 1), // rhofacr
    std::make_pair(0, 1), // rhofaci
    std::make_pair(0, 1), // acn
    std::make_pair(0, 1), // qv
    std::make_pair(0, 1), // th
    std::make_pair(0, qsmall * 2), // qc
    std::make_pair(0, 1), // nc
    std::make_pair(0, qsmall * 2), // qr
    std::make_pair(0, 1), // nr
    std::make_pair(0, qsmall * 2), // qi
    std::make_pair(0, 1), // ni
    std::make_pair(0, 1), // qm
    std::make_pair(0, 1), // bm
    std::make_pair(0, 1), // qc_incld
    std::make_pair(0, 1), // qr_incld
    std::make_pair(0, 1), // qi_incld
    std::make_pair(0, 1), // qm_incld
    std::make_pair(0, 1), // nc_incld
    std::make_pair(0, 1), // nr_incld
    std::make_pair(0, 1), // ni_incld
    std::make_pair(0, 1), // bm_incld
  };

  P3MainPart1Data isds_fortran[] = {
    //              kts, kte, ktop, kbot, kdir, do_predict_nc,        dt, ranges
    P3MainPart1Data(1,  72,    1,   72,    1, false,         1.800E+03, ranges),
    P3MainPart1Data(1,  72,    1,   72,    1, true,          1.800E+03, ranges),
    P3MainPart1Data(1,  72,   72,    1,   -1, false,         1.800E+03, ranges),
    P3MainPart1Data(1,  72,   72,    1,   -1, true,          1.800E+03, ranges),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(P3MainPart1Data);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  P3MainPart1Data isds_cxx[num_runs] = {
    P3MainPart1Data(isds_fortran[0]),
    P3MainPart1Data(isds_fortran[1]),
    P3MainPart1Data(isds_fortran[2]),
    P3MainPart1Data(isds_fortran[3]),
  };

    // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    p3_main_part1(isds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    P3MainPart1Data& d = isds_cxx[i];
    p3_main_part1_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir, d.do_predict_nc, d.dt,
                            d.pres, d.dpres, d.dz, d.nc_nuceat_tend, d.exner, d.inv_exner, d.inv_cld_frac_l, d.inv_cld_frac_i, 
                            d.inv_cld_frac_r, d.latent_heat_vapor, d.latent_heat_sublim, d.latent_heat_fusion,
                            d.t, d.rho, d.inv_rho, d.qv_sat_l, d.qv_sat_i, d.qv_supersat_i, d.rhofacr, d.rhofaci,
                            d.acn, d.qv, d.th, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni, d.qm, d.bm, d.qc_incld, d.qr_incld, d.qi_incld,
                            d.qm_incld, d.nc_incld, d.nr_incld, d.ni_incld, d.bm_incld,
                            &d.is_nucleat_possible, &d.is_hydromet_present);
  }

  for (Int i = 0; i < num_runs; ++i) {
    Int start = std::min(isds_fortran[i].kbot, isds_fortran[i].ktop) - 1; // 0-based indx
    Int end   = std::max(isds_fortran[i].kbot, isds_fortran[i].ktop);     // 0-based indx
    for (Int k = start; k < end; ++k) {
      REQUIRE(isds_fortran[i].t[k]             == isds_cxx[i].t[k]);
      REQUIRE(isds_fortran[i].rho[k]           == isds_cxx[i].rho[k]);
      REQUIRE(isds_fortran[i].inv_rho[k]       == isds_cxx[i].inv_rho[k]);
      REQUIRE(isds_fortran[i].qv_sat_l[k]      == isds_cxx[i].qv_sat_l[k]);
      REQUIRE(isds_fortran[i].qv_sat_i[k]      == isds_cxx[i].qv_sat_i[k]);
      REQUIRE(isds_fortran[i].qv_supersat_i[k] == isds_cxx[i].qv_supersat_i[k]);
      REQUIRE(isds_fortran[i].rhofacr[k]       == isds_cxx[i].rhofacr[k]);
      REQUIRE(isds_fortran[i].rhofaci[k]       == isds_cxx[i].rhofaci[k]);
      REQUIRE(isds_fortran[i].acn[k]           == isds_cxx[i].acn[k]);
      REQUIRE(isds_fortran[i].qv[k]            == isds_cxx[i].qv[k]);
      REQUIRE(isds_fortran[i].th[k]            == isds_cxx[i].th[k]);
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

static void run_bfb_p3_main_part2()
{
  constexpr Scalar qsmall     = C::QSMALL;
  constexpr Scalar zerodegc   = C::ZeroDegC;
  constexpr Scalar sup_upper = -0.05;
  constexpr Scalar sup_lower = -0.1;

  const std::array< std::pair<Real, Real>, P3MainPart2Data::NUM_ARRAYS > ranges = {
    std::make_pair(0, 1), // pres
    std::make_pair(0, 1), // dpres
    std::make_pair(0, 1), // dz
    std::make_pair(0, 1), // npccn
    std::make_pair(0, 1), // exner
    std::make_pair(0, 1), // inv_exner
    std::make_pair(0, 1), // inv_cld_frac_l
    std::make_pair(0, 1), // inv_cld_frac_i
    std::make_pair(0, 1), // inv_cld_frac_r
    std::make_pair(0, 1), // ni_activated
    std::make_pair(0, 1), // inv_qc_relvar
    std::make_pair(0, 1), // cld_frac_i
    std::make_pair(0, 1), // cld_frac_l
    std::make_pair(0, 1), // cld_frac_r
    std::make_pair(zerodegc - 10, zerodegc + 10), // t
    std::make_pair(0, 1), // rho
    std::make_pair(0, 1), // inv_rho
    std::make_pair(0, 1), // qv_sat_l
    std::make_pair(0, 1), // qv_sat_i
    std::make_pair(sup_lower -.05, sup_upper + .05), // qv_supersat_i
    std::make_pair(0, 1), // rhofacr
    std::make_pair(0, 1), // rhofaci
    std::make_pair(0, 1), // acn
    std::make_pair(0, 1), // qv
    std::make_pair(0, 1), // th
    std::make_pair(0, qsmall * 2), // qc
    std::make_pair(0, 1), // nc
    std::make_pair(0, qsmall * 2), // qr
    std::make_pair(0, 1), // nr
    std::make_pair(0, qsmall * 2), // qi
    std::make_pair(0, 1), // ni
    std::make_pair(0, 1), // qm
    std::make_pair(0, 1), // bm
    std::make_pair(0, 1), // latent_heat_vapor
    std::make_pair(0, 1), // latent_heat_sublim
    std::make_pair(0, 1), // latent_heat_fusion
    std::make_pair(0, 1), // qc_incld
    std::make_pair(0, 1), // qr_incld
    std::make_pair(0, 1), // qi_incld
    std::make_pair(0, 1), // qm_incld
    std::make_pair(0, 1), // nc_incld
    std::make_pair(0, 1), // nr_incld
    std::make_pair(0, 1), // ni_incld
    std::make_pair(0, 1), // bm_incld
    std::make_pair(0, 1), // mu_c
    std::make_pair(0, 1), // nu
    std::make_pair(0, 1), // lamc
    std::make_pair(0, 1), // cdist
    std::make_pair(0, 1), // cdist1
    std::make_pair(0, 1), // cdistr
    std::make_pair(0, 1), // mu_r
    std::make_pair(0, 1), // lamr
    std::make_pair(0, 1), // logn0r
    std::make_pair(0, 1), // cmeiout
    std::make_pair(0, 1), // precip_total_tend
    std::make_pair(0, 1), // nevapr
    std::make_pair(0, 1), // qr_evap_tend
    std::make_pair(0, 1), // vap_liq_exchange
    std::make_pair(0, 1), // vap_ice_exchange
    std::make_pair(0, 1), // liq_ice_exchange
    std::make_pair(0, 1), // pratot
    std::make_pair(0, 1)  // prctot
  };

  P3MainPart2Data isds_fortran[] = {
    //              kts, kte, ktop, kbot, kdir, do_predict_nc,        dt, ranges
    P3MainPart2Data(1,  72,    1,   72,    1, false,         1.800E+03, ranges),
    P3MainPart2Data(1,  72,    1,   72,    1, true,          1.800E+03, ranges),
    P3MainPart2Data(1,  72,   72,    1,   -1, false,         1.800E+03, ranges),
    P3MainPart2Data(1,  72,   72,    1,   -1, true,          1.800E+03, ranges),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(P3MainPart2Data);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  P3MainPart2Data isds_cxx[num_runs] = {
    P3MainPart2Data(isds_fortran[0]),
    P3MainPart2Data(isds_fortran[1]),
    P3MainPart2Data(isds_fortran[2]),
    P3MainPart2Data(isds_fortran[3]),
  };

    // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    p3_main_part2(isds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    P3MainPart2Data& d = isds_cxx[i];
    p3_main_part2_f(
      d.kts, d.kte, d.kbot, d.ktop, d.kdir, d.do_predict_nc, d.dt, d.inv_dt,
      d.pres, d.dpres, d.dz, d.nc_nuceat_tend, d.exner, d.inv_exner, d.inv_cld_frac_l, d.inv_cld_frac_i, 
      d.inv_cld_frac_r, d.ni_activated, d.inv_qc_relvar, d.cld_frac_i, d.cld_frac_l, d.cld_frac_r,
      d.t, d.rho, d.inv_rho, d.qv_sat_l, d.qv_sat_i, d.qv_supersat_i, d.rhofacr, d.rhofaci, d.acn, d.qv, d.th, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni,
      d.qm, d.bm, d.latent_heat_vapor, d.latent_heat_sublim, d.latent_heat_fusion, d.qc_incld, d.qr_incld, d.qi_incld, d.qm_incld, d.nc_incld, d.nr_incld,
      d.ni_incld, d.bm_incld, d.mu_c, d.nu, d.lamc, d.cdist, d.cdist1, d.cdistr, d.mu_r, d.lamr, d.logn0r, d.cmeiout, d.precip_total_tend,
      d.nevapr, d.qr_evap_tend, d.vap_liq_exchange, d.vap_ice_exchange, d.liq_ice_exchange, d.pratot,
      d.prctot, &d.is_hydromet_present);
  }

  for (Int i = 0; i < num_runs; ++i) {
    Int start = std::min(isds_fortran[i].kbot, isds_fortran[i].ktop) - 1; // 0-based indx
    Int end   = std::max(isds_fortran[i].kbot, isds_fortran[i].ktop);     // 0-based indx
    for (Int k = start; k < end; ++k) {
      REQUIRE(isds_fortran[i].t[k]                  == isds_cxx[i].t[k]);
      REQUIRE(isds_fortran[i].rho[k]                == isds_cxx[i].rho[k]);
      REQUIRE(isds_fortran[i].inv_rho[k]            == isds_cxx[i].inv_rho[k]);
      REQUIRE(isds_fortran[i].qv_sat_l[k]           == isds_cxx[i].qv_sat_l[k]);
      REQUIRE(isds_fortran[i].qv_sat_i[k]           == isds_cxx[i].qv_sat_i[k]);
      REQUIRE(isds_fortran[i].qv_supersat_i[k]      == isds_cxx[i].qv_supersat_i[k]);
      REQUIRE(isds_fortran[i].rhofacr[k]            == isds_cxx[i].rhofacr[k]);
      REQUIRE(isds_fortran[i].rhofaci[k]            == isds_cxx[i].rhofaci[k]);
      REQUIRE(isds_fortran[i].acn[k]                == isds_cxx[i].acn[k]);
      REQUIRE(isds_fortran[i].qv[k]                 == isds_cxx[i].qv[k]);
      REQUIRE(isds_fortran[i].th[k]                 == isds_cxx[i].th[k]);
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
      REQUIRE(isds_fortran[i].cmeiout[k]            == isds_cxx[i].cmeiout[k]);
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

static void run_bfb_p3_main_part3()
{
  constexpr Scalar qsmall     = C::QSMALL;

  const std::array< std::pair<Real, Real>, P3MainPart3Data::NUM_ARRAYS > ranges = {
    std::make_pair(0, 1), // exner
    std::make_pair(0, 1), // cld_frac_l
    std::make_pair(0, 1), // cld_frac_r
    std::make_pair(0, 1), // rho
    std::make_pair(0, 1), // inv_rho
    std::make_pair(0, 1), // rhofaci
    std::make_pair(0, 1), // qv
    std::make_pair(0, 1), // th
    std::make_pair(0, qsmall * 2), // qc
    std::make_pair(0, 1), // nc
    std::make_pair(0, qsmall * 2), // qr
    std::make_pair(0, 1), // nr
    std::make_pair(0, qsmall * 2), // qi
    std::make_pair(0, 1), // ni
    std::make_pair(0, 1), // qm
    std::make_pair(0, 1), // bm
    std::make_pair(0, 1), // latent_heat_vapor
    std::make_pair(0, 1), // latent_heat_sublim
    std::make_pair(0, 1), // mu_c
    std::make_pair(0, 1), // nu
    std::make_pair(0, 1), // lamc
    std::make_pair(0, 1), // mu_r
    std::make_pair(0, 1), // lamr
    std::make_pair(0, 1), // vap_liq_exchange
    std::make_pair(0, 1), // ze_rain
    std::make_pair(0, 1), // ze_ice
    std::make_pair(0, 1), // diag_vmi
    std::make_pair(0, 1), // diag_effi
    std::make_pair(0, 1), // diag_di
    std::make_pair(0, 1), // rho_qi
    std::make_pair(0, 1), // diag_ze
    std::make_pair(0, 1), // diag_effc
  };

  P3MainPart3Data isds_fortran[] = {
    //               kts, kte, ktop, kbot, kdir, ranges
    P3MainPart3Data(1,  72,    1,   72,    1, ranges),
    P3MainPart3Data(1,  72,    1,   72,    1, ranges),
    P3MainPart3Data(1,  72,   72,    1,   -1, ranges),
    P3MainPart3Data(1,  72,   72,    1,   -1, ranges),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(P3MainPart3Data);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  P3MainPart3Data isds_cxx[num_runs] = {
    P3MainPart3Data(isds_fortran[0]),
    P3MainPart3Data(isds_fortran[1]),
    P3MainPart3Data(isds_fortran[2]),
    P3MainPart3Data(isds_fortran[3]),
  };

    // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    p3_main_part3(isds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    P3MainPart3Data& d = isds_cxx[i];
    p3_main_part3_f(
      d.kts, d.kte, d.kbot, d.ktop, d.kdir,
      d.exner, d.cld_frac_l, d.cld_frac_r,
      d.rho, d.inv_rho, d.rhofaci, d.qv, d.th, d.qc, d.nc, d.qr, d.nr, d.qi, d.ni, d.qm, d.bm, d.latent_heat_vapor, d.latent_heat_sublim,
      d.mu_c, d.nu, d.lamc, d.mu_r, d.lamr, d.vap_liq_exchange,
      d. ze_rain, d.ze_ice, d.diag_vmi, d.diag_effi, d.diag_di, d.rho_qi, d.diag_ze, d.diag_effc);
  }

  for (Int i = 0; i < num_runs; ++i) {
    Int start = std::min(isds_fortran[i].kbot, isds_fortran[i].ktop) - 1; // 0-based indx
    Int end   = std::max(isds_fortran[i].kbot, isds_fortran[i].ktop);     // 0-based indx
    for (Int k = start; k < end; ++k) {
      REQUIRE(isds_fortran[i].rho[k]                == isds_cxx[i].rho[k]);
      REQUIRE(isds_fortran[i].inv_rho[k]            == isds_cxx[i].inv_rho[k]);
      REQUIRE(isds_fortran[i].rhofaci[k]            == isds_cxx[i].rhofaci[k]);
      REQUIRE(isds_fortran[i].qv[k]                 == isds_cxx[i].qv[k]);
      REQUIRE(isds_fortran[i].th[k]                 == isds_cxx[i].th[k]);
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
      REQUIRE(isds_fortran[i].mu_c[k]               == isds_cxx[i].mu_c[k]);
      REQUIRE(isds_fortran[i].nu[k]                 == isds_cxx[i].nu[k]);
      REQUIRE(isds_fortran[i].lamc[k]               == isds_cxx[i].lamc[k]);
      REQUIRE(isds_fortran[i].mu_r[k]               == isds_cxx[i].mu_r[k]);
      REQUIRE(isds_fortran[i].lamr[k]               == isds_cxx[i].lamr[k]);
      REQUIRE(isds_fortran[i].vap_liq_exchange[k]   == isds_cxx[i].vap_liq_exchange[k]);
      REQUIRE(isds_fortran[i].ze_rain[k]            == isds_cxx[i].ze_rain[k]);
      REQUIRE(isds_fortran[i].ze_ice[k]             == isds_cxx[i].ze_ice[k]);
      REQUIRE(isds_fortran[i].diag_vmi[k]           == isds_cxx[i].diag_vmi[k]);
      REQUIRE(isds_fortran[i].diag_effi[k]          == isds_cxx[i].diag_effi[k]);
      REQUIRE(isds_fortran[i].diag_di[k]            == isds_cxx[i].diag_di[k]);
      REQUIRE(isds_fortran[i].rho_qi[k]             == isds_cxx[i].rho_qi[k]);
      REQUIRE(isds_fortran[i].diag_ze[k]            == isds_cxx[i].diag_ze[k]);
      REQUIRE(isds_fortran[i].diag_effc[k]          == isds_cxx[i].diag_effc[k]);
    }
  }
}

static void run_bfb_p3_main()
{
  const std::array< std::pair<Real, Real>, P3MainData::NUM_INPUT_ARRAYS > ranges = {
    std::make_pair(1.00000000E+02 , 9.87111111E+04), // pres
    std::make_pair(1.22776609E+02 , 3.49039167E+04), // dz
    std::make_pair(0              , 0), // nc_nuceat_tend
    std::make_pair(0              , 0), // ni_activated
    std::make_pair(1.37888889E+03 , 1.39888889E+03), // dpres
    std::make_pair(1.00371345E+00 , 3.19721007E+00), // exner
    std::make_pair(1              , 1), // cld_frac_i
    std::make_pair(1              , 1), // cld_frac_l
    std::make_pair(1              , 1), // cld_frac_r
    std::make_pair(1              , 1), // inv_qc_relvar
    std::make_pair(0              , 1.00000000E-04), // qc
    std::make_pair(1.00000000E+06 , 1.00000000E+06), // nc
    std::make_pair(0              , 1.00000000E-05), // qr
    std::make_pair(1.00000000E+06 , 1.00000000E+06), // nr
    std::make_pair(0              , 1.00000000E-04), // qi
    std::make_pair(0              , 1.00000000E-04), // qm
    std::make_pair(1.00000000E+06 , 1.00000000E+06), // ni
    std::make_pair(0              , 1.00000000E-02), // bm
    std::make_pair(0              , 5.00000000E-02), // qv
    std::make_pair(6.72653866E+02 , 1.07954335E+03), // th
  };

  P3MainData isds_fortran[] = {
    //      its,  ite, kts, kte,   it,        dt, do_predict_nc, ranges
    P3MainData(1, 10,   1,  72,    1, 1.800E+03, false, ranges),
    P3MainData(1, 10,   1,  72,    1, 1.800E+03, true,  ranges),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(P3MainData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  P3MainData isds_cxx[num_runs] = {
    P3MainData(isds_fortran[0]),
    P3MainData(isds_fortran[1]),
  };

    // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    p3_main(isds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    P3MainData& d = isds_cxx[i];
    d.transpose<ekat::util::TransposeDirection::c2f>();
    p3_main_f(
      d.qc, d.nc, d.qr, d.nr, d.th, d.qv, d.dt, d.qi, d.qm, d.ni,
      d.bm, d.pres, d.dz, d.nc_nuceat_tend, d.ni_activated, d.inv_qc_relvar, d.it, d.precip_liq_surf,
      d.precip_ice_surf, d.its, d.ite, d.kts, d.kte, d.diag_effc, d.diag_effi,
      d.rho_qi, d.do_predict_nc, d.dpres, d.exner, d.cmeiout, d.precip_total_tend,
      d.nevapr, d.qr_evap_tend, d.precip_liq_flux, d.precip_ice_flux, d.cld_frac_r, d.cld_frac_l, d.cld_frac_i, d.mu_c,
      d.lamc, d.liq_ice_exchange, d.vap_liq_exchange, d.vap_ice_exchange);
    d.transpose<ekat::util::TransposeDirection::f2c>();
  }

  for (Int i = 0; i < num_runs; ++i) {
    for (Int t = 0; t < isds_fortran[i].nt(); ++t) {
      REQUIRE(isds_fortran[i].qc[t]                == isds_cxx[i].qc[t]);
      REQUIRE(isds_fortran[i].nc[t]                == isds_cxx[i].nc[t]);
      REQUIRE(isds_fortran[i].qr[t]                == isds_cxx[i].qr[t]);
      REQUIRE(isds_fortran[i].nr[t]                == isds_cxx[i].nr[t]);
      REQUIRE(isds_fortran[i].qi[t]                == isds_cxx[i].qi[t]);
      REQUIRE(isds_fortran[i].qm[t]                == isds_cxx[i].qm[t]);
      REQUIRE(isds_fortran[i].ni[t]                == isds_cxx[i].ni[t]);
      REQUIRE(isds_fortran[i].bm[t]                == isds_cxx[i].bm[t]);
      REQUIRE(isds_fortran[i].qv[t]                == isds_cxx[i].qv[t]);
      REQUIRE(isds_fortran[i].th[t]                == isds_cxx[i].th[t]);
      REQUIRE(isds_fortran[i].diag_effc[t]         == isds_cxx[i].diag_effc[t]);
      REQUIRE(isds_fortran[i].diag_effi[t]         == isds_cxx[i].diag_effi[t]);
      REQUIRE(isds_fortran[i].rho_qi[t]            == isds_cxx[i].rho_qi[t]);
      REQUIRE(isds_fortran[i].mu_c[t]              == isds_cxx[i].mu_c[t]);
      REQUIRE(isds_fortran[i].lamc[t]              == isds_cxx[i].lamc[t]);
      REQUIRE(isds_fortran[i].cmeiout[t]           == isds_cxx[i].cmeiout[t]);
      REQUIRE(isds_fortran[i].precip_total_tend[t] == isds_cxx[i].precip_total_tend[t]);
      REQUIRE(isds_fortran[i].nevapr[t]            == isds_cxx[i].nevapr[t]);
      REQUIRE(isds_fortran[i].qr_evap_tend[t]      == isds_cxx[i].qr_evap_tend[t]);
      REQUIRE(isds_fortran[i].liq_ice_exchange[t]  == isds_cxx[i].liq_ice_exchange[t]);
      REQUIRE(isds_fortran[i].vap_liq_exchange[t]  == isds_cxx[i].vap_liq_exchange[t]);
      REQUIRE(isds_fortran[i].vap_ice_exchange[t]  == isds_cxx[i].vap_ice_exchange[t]);
      REQUIRE(isds_fortran[i].precip_liq_flux[t]   == isds_cxx[i].precip_liq_flux[t]);
      REQUIRE(isds_fortran[i].precip_ice_flux[t]   == isds_cxx[i].precip_ice_flux[t]);
      REQUIRE(isds_fortran[i].precip_liq_surf[t]   == isds_cxx[i].precip_liq_surf[t]);
      REQUIRE(isds_fortran[i].precip_ice_surf[t]   == isds_cxx[i].precip_ice_surf[t]);
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
