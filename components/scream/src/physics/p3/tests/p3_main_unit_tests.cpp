#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"

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

static void run_phys_p3_main_pre_loop()
{
  // TODO
}

static void run_phys()
{
  run_phys_p3_main_pre_loop();
}

static void run_bfb_p3_main_pre_loop()
{
  constexpr Scalar qsmall     = C::QSMALL;
  constexpr Scalar zerodegc   = C::ZeroDegC;
  constexpr Scalar sup_upper = -0.05;
  constexpr Scalar sup_lower = -0.1;

  const std::array< std::pair<Real, Real>, P3MainPreLoopData::NUM_ARRAYS > ranges = {
    std::make_pair(0, 1), // pres
    std::make_pair(0, 1), // pdel
    std::make_pair(0, 1), // dzq
    std::make_pair(0, 1), // npccn
    std::make_pair(0, 1), // exner
    std::make_pair(0, 1), // inv_exner
    std::make_pair(0, 1), // inv_lcldm
    std::make_pair(0, 1), // inv_icldm
    std::make_pair(0, 1), // inv_rcldm
    std::make_pair(0, 1), // xxlv
    std::make_pair(0, 1), // xxls
    std::make_pair(0, 1), // xlf
    std::make_pair(zerodegc - 10, zerodegc + 10), // t
    std::make_pair(0, 1), // rho
    std::make_pair(0, 1), // inv_rho
    std::make_pair(0, 1), // qvs
    std::make_pair(0, 1), // qvi
    std::make_pair(sup_lower -.05, sup_upper + .05), // sup
    std::make_pair(sup_lower -.05, sup_upper + .05), // supi
    std::make_pair(0, 1), // rhofacr
    std::make_pair(0, 1), // rhofaci
    std::make_pair(0, 1), // acn
    std::make_pair(0, 1), // qv
    std::make_pair(0, 1), // th
    std::make_pair(0, qsmall * 2), // qc
    std::make_pair(0, 1), // nc
    std::make_pair(0, qsmall * 2), // qr
    std::make_pair(0, 1), // nr
    std::make_pair(0, qsmall * 2), // qitot
    std::make_pair(0, 1), // nitot
    std::make_pair(0, 1), // qirim
    std::make_pair(0, 1), // birim
    std::make_pair(0, 1), // qc_incld
    std::make_pair(0, 1), // qr_incld
    std::make_pair(0, 1), // qitot_incld
    std::make_pair(0, 1), // qirim_incld
    std::make_pair(0, 1), // nc_incld
    std::make_pair(0, 1), // nr_incld
    std::make_pair(0, 1), // nitot_incld
    std::make_pair(0, 1), // birim_incld
  };

  P3MainPreLoopData isds_fortran[] = {
    //              kts, kte, ktop, kbot, kdir, log_predictNc,        dt, ranges
    P3MainPreLoopData(1,  1,    1,   72,   -1, false,         1.800E+03, ranges),
    P3MainPreLoopData(1,  1,    1,   72,   -1, true,          1.800E+03, ranges),
    P3MainPreLoopData(1,  1,   72,    1,    1, false,         1.800E+03, ranges),
    P3MainPreLoopData(1,  1,   72,    1,    1, true,          1.800E+03, ranges),
  };

  static constexpr Int num_runs = sizeof(isds_fortran) / sizeof(P3MainPreLoopData);

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  P3MainPreLoopData isds_cxx[num_runs] = {
    P3MainPreLoopData(isds_fortran[0]),
    P3MainPreLoopData(isds_fortran[1]),
    P3MainPreLoopData(isds_fortran[2]),
    P3MainPreLoopData(isds_fortran[3]),
  };

    // Get data from fortran
  for (Int i = 0; i < num_runs; ++i) {
    p3_main_pre_main_loop(isds_fortran[i]);
  }

  // Get data from cxx
  for (Int i = 0; i < num_runs; ++i) {
    P3MainPreLoopData& d = isds_cxx[i];
    p3_main_pre_main_loop_f(d.kts, d.kte, d.ktop, d.kbot, d.kdir, d.log_predictNc, d.dt,
                            d.pres, d.pdel, d.dzq, d.npccn, d.exner, d.inv_exner, d.inv_lcldm, d.inv_icldm, d.inv_rcldm, d.xxlv, d.xxls, d.xlf,
                            d.t, d.rho, d.inv_rho, d.qvs, d.qvi, d.sup, d.supi, d.rhofacr, d.rhofaci,
                            d.acn, d.qv, d.th, d.qc, d.nc, d.qr, d.nr, d.qitot, d.nitot, d.qirim, d.birim, d.qc_incld, d.qr_incld, d.qitot_incld,
                            d.qirim_incld, d.nc_incld, d.nr_incld, d.nitot_incld, d.birim_incld,
                            &d.log_nucleationPossible, &d.log_hydrometeorsPresent);
  }

  for (Int i = 0; i < num_runs; ++i) {
    Int start = std::min(isds_fortran[i].kbot, isds_fortran[i].ktop) - 1; // 0-based indx
    Int end   = std::max(isds_fortran[i].kbot, isds_fortran[i].ktop);     // 0-based indx
    for (Int k = start; k < end; ++k) {
      REQUIRE(isds_fortran[i].t[k]           == isds_cxx[i].t[k]);
      REQUIRE(isds_fortran[i].rho[k]         == isds_cxx[i].rho[k]);
      REQUIRE(isds_fortran[i].inv_rho[k]     == isds_cxx[i].inv_rho[k]);
      REQUIRE(isds_fortran[i].qvs[k]         == isds_cxx[i].qvs[k]);
      REQUIRE(isds_fortran[i].qvi[k]         == isds_cxx[i].qvi[k]);
      REQUIRE(isds_fortran[i].sup[k]         == isds_cxx[i].sup[k]);
      REQUIRE(isds_fortran[i].supi[k]        == isds_cxx[i].supi[k]);
      REQUIRE(isds_fortran[i].rhofacr[k]     == isds_cxx[i].rhofacr[k]);
      REQUIRE(isds_fortran[i].rhofaci[k]     == isds_cxx[i].rhofaci[k]);
      REQUIRE(isds_fortran[i].acn[k]         == isds_cxx[i].acn[k]);
      REQUIRE(isds_fortran[i].qv[k]          == isds_cxx[i].qv[k]);
      REQUIRE(isds_fortran[i].th[k]          == isds_cxx[i].th[k]);
      REQUIRE(isds_fortran[i].qc[k]          == isds_cxx[i].qc[k]);
      REQUIRE(isds_fortran[i].nc[k]          == isds_cxx[i].nc[k]);
      REQUIRE(isds_fortran[i].qr[k]          == isds_cxx[i].qr[k]);
      REQUIRE(isds_fortran[i].nr[k]          == isds_cxx[i].nr[k]);
      REQUIRE(isds_fortran[i].qitot[k]       == isds_cxx[i].qitot[k]);
      REQUIRE(isds_fortran[i].nitot[k]       == isds_cxx[i].nitot[k]);
      REQUIRE(isds_fortran[i].qirim[k]       == isds_cxx[i].qirim[k]);
      REQUIRE(isds_fortran[i].birim[k]       == isds_cxx[i].birim[k]);
      REQUIRE(isds_fortran[i].qc_incld[k]    == isds_cxx[i].qc_incld[k]);
      REQUIRE(isds_fortran[i].qr_incld[k]    == isds_cxx[i].qr_incld[k]);
      REQUIRE(isds_fortran[i].qitot_incld[k] == isds_cxx[i].qitot_incld[k]);
      REQUIRE(isds_fortran[i].qirim_incld[k] == isds_cxx[i].qirim_incld[k]);
      REQUIRE(isds_fortran[i].nc_incld[k]    == isds_cxx[i].nc_incld[k]);
      REQUIRE(isds_fortran[i].nr_incld[k]    == isds_cxx[i].nr_incld[k]);
      REQUIRE(isds_fortran[i].nitot_incld[k] == isds_cxx[i].nitot_incld[k]);
      REQUIRE(isds_fortran[i].birim_incld[k] == isds_cxx[i].birim_incld[k]);
    }
    REQUIRE(isds_fortran[i].log_hydrometeorsPresent == isds_cxx[i].log_hydrometeorsPresent);
    REQUIRE(isds_fortran[i].log_nucleationPossible  == isds_cxx[i].log_nucleationPossible);
  }
}

static void run_bfb()
{
  run_bfb_p3_main_pre_loop();
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
