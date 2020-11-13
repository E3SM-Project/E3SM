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
struct UnitWrap::UnitTest<D>::TestNrConservation {

  static void run_bfb()
  {
    NrConservationData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(NrConservationData);

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    NrConservationData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      nr_conservation(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      nr_conservation_f(d.nr, d.ni2nr_melt_tend, d.nr_ice_shed_tend, d.ncshdc, d.nc2nr_autoconv_tend, d.dt, &d.nr_collect_tend, &d.nr2ni_immers_freeze_tend, &d.nr_selfcollect_tend, &d.nr_evap_tend);
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      NrConservationData& d_f90 = f90_data[i];
      NrConservationData& d_cxx = cxx_data[i];
      REQUIRE(d_f90.nr_collect_tend == d_cxx.nr_collect_tend);
      REQUIRE(d_f90.nr2ni_immers_freeze_tend == d_cxx.nr2ni_immers_freeze_tend);
      REQUIRE(d_f90.nr_selfcollect_tend == d_cxx.nr_selfcollect_tend);
      REQUIRE(d_f90.nr_evap_tend == d_cxx.nr_evap_tend);

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("nr_conservation_bfb", "[p3]")
{
  using TestStruct = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestNrConservation;

  TestStruct::run_bfb();
}

} // empty namespace
