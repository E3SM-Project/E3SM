#include "catch2/catch.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestVdShocDecompandSolve : public UnitWrap::UnitTest<D>::Base {

  void run_bfb()
  {
    auto engine = Base::get_engine();

    VdShocDecompandSolveData baseline_data[] = {
      // shcol, nlev, nlevi, dtime, n_rhs
      VdShocDecompandSolveData(10, 71, 72, 5, 19),
      VdShocDecompandSolveData(10, 12, 13, 2.5, 7),
      VdShocDecompandSolveData(7, 16, 17, 1, 2),
      VdShocDecompandSolveData(2, 7, 8, 1, 1)
    };

    static constexpr Int num_runs = sizeof(baseline_data) / sizeof(VdShocDecompandSolveData);

    // Generate random input data. Diagonals in solver data will be overwritten
    // after results of decomp routine.
    for (Int i = 0; i < num_runs; ++i) {
      VdShocDecompandSolveData& d_baseline = baseline_data[i];
      d_baseline.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    VdShocDecompandSolveData cxx_data[] = {
      VdShocDecompandSolveData(baseline_data[0]),
      VdShocDecompandSolveData(baseline_data[1]),
      VdShocDecompandSolveData(baseline_data[2]),
      VdShocDecompandSolveData(baseline_data[3])
    };

    // Assume all data is in C layout

    // Read baseline data.
    if (this->m_baseline_action == COMPARE) {
      for (auto& d: baseline_data) {
        d.read(Base::m_fid);
      }
    }

    // Get data from cxx
    for (Int i = 0; i < num_runs; ++i) {
      VdShocDecompandSolveData& d_cxx = cxx_data[i];
      vd_shoc_decomp_and_solve(d_cxx);
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        VdShocDecompandSolveData& d_baseline = baseline_data[i];
        VdShocDecompandSolveData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_baseline.total(d_baseline.var); ++k) {
          REQUIRE(d_baseline.total(d_baseline.var) == d_cxx.total(d_cxx.var));
          REQUIRE(d_baseline.var[k] == d_cxx.var[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        cxx_data[i].write(Base::m_fid);
      }
    }
  } // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("vd_shoc_solve_bfb", "[shoc]")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestVdShocDecompandSolve;

  TestStruct().run_bfb();
}

} // empty namespace
