#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestShocMain {

  static void run_bfb()
  {
    ShocMainData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(ShocMainData);

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ShocMainData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      shoc_main(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      shoc_main_f(d.shcol, d.nlev, d.nlevi, d.dtime, d.nadv, d.host_dx, d.host_dy, d.thv, d.zt_grid, d.zi_grid, d.pres, d.presi, d.pdel, d.wthl_sfc, d.wqw_sfc, d.uw_sfc, d.vw_sfc, d.wtracer_sfc, d.num_qtracers, d.w_field, d.exner, d.phis, d.host_dse, d.tke, d.thetal, d.qw, d.u_wind, d.v_wind, d.qtracers, d.wthv_sec, d.tkh, d.tk, d.shoc_ql, d.shoc_cldfrac, d.pblh, d.shoc_mix, d.isotropy, d.w_sec, d.thl_sec, d.qw_sec, d.qwthl_sec, d.wthl_sec, d.wqw_sec, d.wtke_sec, d.uw_sec, d.vw_sec, d.w3, d.wqls_sec, d.brunt, d.shoc_ql2);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      ShocMainData& d_f90 = f90_data[i];
      ShocMainData& d_cxx = cxx_data[i];
      for (Int k = 0; k < d_f90.total(d_f90.host_dse); ++k) {
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.host_dse));
        REQUIRE(d_f90.host_dse[k] == d_cxx.host_dse[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.tke));
        REQUIRE(d_f90.tke[k] == d_cxx.tke[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.thetal));
        REQUIRE(d_f90.thetal[k] == d_cxx.thetal[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.qw));
        REQUIRE(d_f90.qw[k] == d_cxx.qw[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.u_wind));
        REQUIRE(d_f90.u_wind[k] == d_cxx.u_wind[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.v_wind));
        REQUIRE(d_f90.v_wind[k] == d_cxx.v_wind[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.wthv_sec));
        REQUIRE(d_f90.wthv_sec[k] == d_cxx.wthv_sec[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.tkh));
        REQUIRE(d_f90.tkh[k] == d_cxx.tkh[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.tk));
        REQUIRE(d_f90.tk[k] == d_cxx.tk[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.shoc_ql));
        REQUIRE(d_f90.shoc_ql[k] == d_cxx.shoc_ql[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.shoc_cldfrac));
        REQUIRE(d_f90.shoc_cldfrac[k] == d_cxx.shoc_cldfrac[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.shoc_mix));
        REQUIRE(d_f90.shoc_mix[k] == d_cxx.shoc_mix[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.isotropy));
        REQUIRE(d_f90.isotropy[k] == d_cxx.isotropy[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.w_sec));
        REQUIRE(d_f90.w_sec[k] == d_cxx.w_sec[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.wqls_sec));
        REQUIRE(d_f90.wqls_sec[k] == d_cxx.wqls_sec[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.brunt));
        REQUIRE(d_f90.brunt[k] == d_cxx.brunt[k]);
        REQUIRE(d_f90.total(d_f90.host_dse) == d_cxx.total(d_cxx.shoc_ql2));
        REQUIRE(d_f90.shoc_ql2[k] == d_cxx.shoc_ql2[k]);
      }
      for (Int k = 0; k < d_f90.total(d_f90.qtracers); ++k) {
        REQUIRE(d_f90.total(d_f90.qtracers) == d_cxx.total(d_cxx.qtracers));
        REQUIRE(d_f90.qtracers[k] == d_cxx.qtracers[k]);
      }
      for (Int k = 0; k < d_f90.total(d_f90.pblh); ++k) {
        REQUIRE(d_f90.total(d_f90.pblh) == d_cxx.total(d_cxx.pblh));
        REQUIRE(d_f90.pblh[k] == d_cxx.pblh[k]);
      }
      for (Int k = 0; k < d_f90.total(d_f90.thl_sec); ++k) {
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.thl_sec));
        REQUIRE(d_f90.thl_sec[k] == d_cxx.thl_sec[k]);
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.qw_sec));
        REQUIRE(d_f90.qw_sec[k] == d_cxx.qw_sec[k]);
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.qwthl_sec));
        REQUIRE(d_f90.qwthl_sec[k] == d_cxx.qwthl_sec[k]);
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.wthl_sec));
        REQUIRE(d_f90.wthl_sec[k] == d_cxx.wthl_sec[k]);
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.wqw_sec));
        REQUIRE(d_f90.wqw_sec[k] == d_cxx.wqw_sec[k]);
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.wtke_sec));
        REQUIRE(d_f90.wtke_sec[k] == d_cxx.wtke_sec[k]);
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.uw_sec));
        REQUIRE(d_f90.uw_sec[k] == d_cxx.uw_sec[k]);
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.vw_sec));
        REQUIRE(d_f90.vw_sec[k] == d_cxx.vw_sec[k]);
        REQUIRE(d_f90.total(d_f90.thl_sec) == d_cxx.total(d_cxx.w3));
        REQUIRE(d_f90.w3[k] == d_cxx.w3[k]);
      }

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("shoc_main_bfb", "[shoc]")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestShocMain;

  TestStruct::run_bfb();
}

} // empty namespace
