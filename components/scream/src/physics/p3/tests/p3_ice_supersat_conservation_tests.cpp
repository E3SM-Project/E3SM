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
struct UnitWrap::UnitTest<D>::TestIceSupersatConservation {

  static void run_bfb()
  {
    IceSupersatConservationData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(IceSupersatConservationData);

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    IceSupersatConservationData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      ice_supersat_conservation(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      ice_supersat_conservation_f(&d.qidep, &d.qinuc, d.cld_frac_i, d.qv, d.qv_sat_i, d.latent_heat_sublim, d.t_atm, d.dt);
    }

    // Verify BFB results, all data should be in C layout
    for (Int i = 0; i < num_runs; ++i) {
      IceSupersatConservationData& d_f90 = f90_data[i];
      IceSupersatConservationData& d_cxx = cxx_data[i];
      REQUIRE(d_f90.qidep == d_cxx.qidep);
      REQUIRE(d_f90.qinuc == d_cxx.qinuc);

    }
  } // run_bfb

};

} // namespace unit_test
} // namespace p3
} // namespace scream

namespace {

TEST_CASE("ice_supersat_conservation_bfb", "[p3]")
{
  using TestStruct = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIceSupersatConservation;

  TestStruct::run_bfb();
}

} // empty namespace
