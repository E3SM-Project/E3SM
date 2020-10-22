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
struct UnitWrap::UnitTest<D>::TestUpdatePrognosticsImplicit {

  static void run_bfb()
  {
    UpdatePrognosticsImplicitData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(UpdatePrognosticsImplicitData);

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize();
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    UpdatePrognosticsImplicitData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      update_prognostics_implicit(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      update_prognostics_implicit_f(d.shcol, d.nlev, d.nlevi, d.num_tracer, d.dtime, d.dz_zt, d.dz_zi, d.rho_zt, d.zt_grid, d.zi_grid, d.tk, d.tkh, d.uw_sfc, d.vw_sfc, d.wthl_sfc, d.wqw_sfc, d.wtracer_sfc, d.thetal, d.qw, d.tracer, d.tke, d.u_wind, d.v_wind);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      UpdatePrognosticsImplicitData& d_f90 = f90_data[i];
      UpdatePrognosticsImplicitData& d_cxx = cxx_data[i];
      for (Int k = 0; k < d_f90.total(d_f90.thetal); ++k) {
        REQUIRE(d_f90.total(d_f90.thetal) == d_cxx.total(d_cxx.thetal));
        REQUIRE(d_f90.thetal[k] == d_cxx.thetal[k]);
        REQUIRE(d_f90.total(d_f90.thetal) == d_cxx.total(d_cxx.qw));
        REQUIRE(d_f90.qw[k] == d_cxx.qw[k]);
        REQUIRE(d_f90.total(d_f90.thetal) == d_cxx.total(d_cxx.tke));
        REQUIRE(d_f90.tke[k] == d_cxx.tke[k]);
        REQUIRE(d_f90.total(d_f90.thetal) == d_cxx.total(d_cxx.u_wind));
        REQUIRE(d_f90.u_wind[k] == d_cxx.u_wind[k]);
        REQUIRE(d_f90.total(d_f90.thetal) == d_cxx.total(d_cxx.v_wind));
        REQUIRE(d_f90.v_wind[k] == d_cxx.v_wind[k]);
      }
      for (Int k = 0; k < d_f90.total(d_f90.tracer); ++k) {
        REQUIRE(d_f90.total(d_f90.tracer) == d_cxx.total(d_cxx.tracer));
        REQUIRE(d_f90.tracer[k] == d_cxx.tracer[k]);
      }

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("update_prognostics_implicit_bfb", "[shoc]")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestUpdatePrognosticsImplicit;

  TestStruct::run_bfb();
}

} // empty namespace