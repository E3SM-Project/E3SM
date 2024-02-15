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
struct UnitWrap::UnitTest<D>::TestAdvanceIopNudging {

  static void run_bfb()
  {
    auto engine = setup_random_test();

    AdvanceIopNudgingData f90_data[] = {
      //                    plev, scm_dt, ps_in
      AdvanceIopNudgingData(72,   0.1,    1000),
      AdvanceIopNudgingData(27,   0.1,    1000),
      AdvanceIopNudgingData(32,   0.1,    1000),
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(AdvanceIopNudgingData);

    // Generate random input data
    // Alternatively, you can use the f90_data construtors/initializer lists to hardcode data
    for (auto& d : f90_data) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    AdvanceIopNudgingData cxx_data[] = {
      AdvanceIopNudgingData(f90_data[0]),
      AdvanceIopNudgingData(f90_data[1]),
      AdvanceIopNudgingData(f90_data[2]),
   };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      advance_iop_nudging(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      advance_iop_nudging_f(d.plev, d.scm_dt, d.ps_in, d.t_in, d.q_in, d.tobs, d.qobs,
                            d.hyai, d.hyam, d.hybi, d.hybm,
                            d.t_update, d.q_update, d.relaxt, d.relaxq);
    }

    // We can't call into fortran. Due to all the dependencies it has, it's not possible
    // to build it in standalone eamxx. Without fortran, we cannot do BFB tests.
#if 0
    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      for (Int i = 0; i < num_runs; ++i) {
        AdvanceIopNudgingData& d_f90 = f90_data[i];
        AdvanceIopNudgingData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.t_update); ++k) {
          REQUIRE(d_f90.total(d_f90.t_update) == d_cxx.total(d_cxx.t_update));
          REQUIRE(d_f90.t_update[k] == d_cxx.t_update[k]);
          REQUIRE(d_f90.total(d_f90.t_update) == d_cxx.total(d_cxx.q_update));
          REQUIRE(d_f90.q_update[k] == d_cxx.q_update[k]);
          REQUIRE(d_f90.total(d_f90.t_update) == d_cxx.total(d_cxx.relaxt));
          REQUIRE(d_f90.relaxt[k] == d_cxx.relaxt[k]);
          REQUIRE(d_f90.total(d_f90.t_update) == d_cxx.total(d_cxx.relaxq));
          REQUIRE(d_f90.relaxq[k] == d_cxx.relaxq[k]);
        }

      }
    }
#endif
  } // run_bfb
};

} // namespace unit_test
} // namespace dp
} // namespace scream

namespace {

TEST_CASE("advance_iop_nudging_bfb", "[dp]")
{
  using TestStruct = scream::dp::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestAdvanceIopNudging;

  TestStruct::run_bfb();
}

} // empty namespace
