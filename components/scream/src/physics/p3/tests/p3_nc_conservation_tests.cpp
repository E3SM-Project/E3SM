#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"

#include "p3_unit_tests_common.hpp"

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestNcConservation {

  static void run_bfb()
  {
    NcConservationData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(NcConservationData);

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    NcConservationData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      nc_conservation(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      nc_conservation_f(d.nc, d.nc_selfcollect_tend, d.dt, &d.nc_collect_tend, &d.nc2ni_immers_freeze_tend, &d.nc_accret_tend, &d.nc2nr_autoconv_tend);
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      NcConservationData& d_f90 = f90_data[i];
      NcConservationData& d_cxx = cxx_data[i];
      REQUIRE(d_f90.nc_collect_tend == d_cxx.nc_collect_tend);
      REQUIRE(d_f90.nc2ni_immers_freeze_tend == d_cxx.nc2ni_immers_freeze_tend);
      REQUIRE(d_f90.nc_accret_tend == d_cxx.nc_accret_tend);
      REQUIRE(d_f90.nc2nr_autoconv_tend == d_cxx.nc2nr_autoconv_tend);

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("nc_conservation_bfb", "[p3]")
{
  using TestStruct = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestNcConservation;

  TestStruct::run_bfb();
}

} // empty namespace
