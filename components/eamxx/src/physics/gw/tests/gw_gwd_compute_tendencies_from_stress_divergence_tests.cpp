#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/gw/gw_functions.hpp"
#include "physics/gw/gw_functions_f90.hpp"

#include "gw_unit_tests_common.hpp"

namespace scream {
namespace gw {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestGwdComputeTendenciesFromStressDivergence {

  static void run_bfb()
  {
    auto engine = setup_random_test();

    GwdComputeTendenciesFromStressDivergenceData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(GwdComputeTendenciesFromStressDivergenceData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    GwdComputeTendenciesFromStressDivergenceData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      gwd_compute_tendencies_from_stress_divergence(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      gwd_compute_tendencies_from_stress_divergence_f(d.pver, d.-pgwv:pgwv, d.0:pver, d.-ngwv:ngwv, d.ncol, d.ngwv, d.do_taper, d.dt, d.effgw, d.tend_level, d.lat, d.dpm, d.rdpm, d.c, d.ubm, d.t, d.nm, d.xv, d.yv, d.tau, d.gwut, d.utgw, d.vtgw);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      for (Int i = 0; i < num_runs; ++i) {
        GwdComputeTendenciesFromStressDivergenceData& d_f90 = f90_data[i];
        GwdComputeTendenciesFromStressDivergenceData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.tau); ++k) {
          REQUIRE(d_f90.total(d_f90.tau) == d_cxx.total(d_cxx.tau));
          REQUIRE(d_f90.tau[k] == d_cxx.tau[k]);
        }
        for (Int k = 0; k < d_f90.total(d_f90.gwut); ++k) {
          REQUIRE(d_f90.total(d_f90.gwut) == d_cxx.total(d_cxx.gwut));
          REQUIRE(d_f90.gwut[k] == d_cxx.gwut[k]);
        }
        for (Int k = 0; k < d_f90.total(d_f90.utgw); ++k) {
          REQUIRE(d_f90.total(d_f90.utgw) == d_cxx.total(d_cxx.utgw));
          REQUIRE(d_f90.utgw[k] == d_cxx.utgw[k]);
          REQUIRE(d_f90.total(d_f90.utgw) == d_cxx.total(d_cxx.vtgw));
          REQUIRE(d_f90.vtgw[k] == d_cxx.vtgw[k]);
        }

      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace gw
} // namespace scream

namespace {

TEST_CASE("gwd_compute_tendencies_from_stress_divergence_bfb", "[gw]")
{
  using TestStruct = scream::gw::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestGwdComputeTendenciesFromStressDivergence;

  TestStruct::run_bfb();
}

} // empty namespace
