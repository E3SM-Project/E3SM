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
struct UnitWrap::UnitTest<D>::TestNiConservation {

  static void run_bfb()
  {
    NiConservationData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(NiConservationData);

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    NiConservationData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      ni_conservation(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      ni_conservation_f(d.ni, d.ni_nucleat_tend, d.nr2ni_immers_freeze_tend, d.nc2ni_immers_freeze_tend, d.dt, &d.ni2nr_melt_tend, &d.ni_sublim_tend, &d.ni_selfcollect_tend);
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      NiConservationData& d_f90 = f90_data[i];
      NiConservationData& d_cxx = cxx_data[i];
      REQUIRE(d_f90.ni2nr_melt_tend == d_cxx.ni2nr_melt_tend);
      REQUIRE(d_f90.ni_sublim_tend == d_cxx.ni_sublim_tend);
      REQUIRE(d_f90.ni_selfcollect_tend == d_cxx.ni_selfcollect_tend);

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("ni_conservation_bfb", "[p3]")
{
  using TestStruct = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestNiConservation;

  TestStruct::run_bfb();
}

} // empty namespace
