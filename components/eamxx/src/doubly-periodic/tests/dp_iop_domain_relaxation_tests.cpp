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
struct UnitWrap::UnitTest<D>::TestIopDomainRelaxation {

  static void run_bfb()
  {
    auto engine = setup_random_test();

    IopDomainRelaxationData f90_data[] = {
      // TODO
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(IopDomainRelaxationData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    IopDomainRelaxationData cxx_data[] = {
      // TODO
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      iop_domain_relaxation(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      iop_domain_relaxation_f(d.nelemd, d.np, d.nlev, d.elem, d.hvcoord, d.hybrid, d.t1, d.dp, d.nelemd_todo, d.np_todo, d.dt);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      for (Int i = 0; i < num_runs; ++i) {
        IopDomainRelaxationData& d_f90 = f90_data[i];
        IopDomainRelaxationData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.dp); ++k) {
          REQUIRE(d_f90.total(d_f90.dp) == d_cxx.total(d_cxx.dp));
          REQUIRE(d_f90.dp[k] == d_cxx.dp[k]);
        }

      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace dp
} // namespace scream

namespace {

TEST_CASE("iop_domain_relaxation_bfb", "[dp]")
{
  using TestStruct = scream::dp::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestIopDomainRelaxation;

  TestStruct::run_bfb();
}

} // empty namespace
