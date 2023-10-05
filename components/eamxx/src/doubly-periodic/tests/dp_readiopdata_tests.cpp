#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "doubly-periodic/dp_functions.hpp"
#include "doubly-periodic/dp_functions_f90.hpp"

#include "dp_unit_tests_common.hpp"

namespace scream {
namespace dp {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestReadiopdata {

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ReadiopdataData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(ReadiopdataData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ReadiopdataData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      readiopdata(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      readiopdata_f(d.plev, d.iop_update_phase1, d.hyam, d.hybm);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      for (Int i = 0; i < num_runs; ++i) {
        ReadiopdataData& d_f90 = f90_data[i];
        ReadiopdataData& d_cxx = cxx_data[i];

      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace dp
} // namespace scream

namespace {

TEST_CASE("readiopdata_bfb", "[dp]")
{
  using TestStruct = scream::dp::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestReadiopdata;

  TestStruct::run_bfb();
}

} // empty namespace
