#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "shoc_unit_tests_common.hpp"

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestVdShocDecompandSolve {

  static void run_bfb()
  {
    auto engine = setup_random_test();

    VdShocDecompandSolveData f90_data[] = {
      // shcol, nlev, nlevi, dtime, n_rhs
      VdShocDecompandSolveData(10, 71, 72, 5, 19),
      VdShocDecompandSolveData(10, 12, 13, 2.5, 7),
      VdShocDecompandSolveData(7, 16, 17, 1, 2),
      VdShocDecompandSolveData(2, 7, 8, 1, 1)
    };

    static constexpr Int num_runs = sizeof(f90_data) / sizeof(VdShocDecompandSolveData);

    // Generate random input data. Diagonals in solver data will be overwritten
    // after results of decomp routine.
    for (Int i = 0; i < num_runs; ++i) {
      VdShocDecompandSolveData& d_f90 = f90_data[i];
      d_f90.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    VdShocDecompandSolveData cxx_data[] = {
      VdShocDecompandSolveData(f90_data[0]),
      VdShocDecompandSolveData(f90_data[1]),
      VdShocDecompandSolveData(f90_data[2]),
      VdShocDecompandSolveData(f90_data[3])
    };

    // Assume all data is in C layout

    // Get data from fortran.
    for (Int i = 0; i < num_runs; ++i) {
      VdShocDecompandSolveData& d_f90 = f90_data[i];
      // expects data in C layout
      vd_shoc_decomp_and_solve(d_f90);
    }

    // Get data from cxx
    for (Int i = 0; i < num_runs; ++i) {
      VdShocDecompandSolveData& d_cxx = cxx_data[i];

      d_cxx.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      vd_shoc_decomp_and_solve_f(d_cxx.shcol, d_cxx.nlev, d_cxx.nlevi, d_cxx.n_rhs,
                                 d_cxx.kv_term, d_cxx.tmpi, d_cxx.rdp_zt,
                                 d_cxx.dtime, d_cxx.flux, d_cxx.var);
      d_cxx.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      for (Int i = 0; i < num_runs; ++i) {
        VdShocDecompandSolveData& d_f90 = f90_data[i];
        VdShocDecompandSolveData& d_cxx = cxx_data[i];
        for (Int k = 0; k < d_f90.total(d_f90.var); ++k) {
          REQUIRE(d_f90.total(d_f90.var) == d_cxx.total(d_cxx.var));
          REQUIRE(d_f90.var[k] == d_cxx.var[k]);
        }
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

  TestStruct::run_bfb();
}

} // empty namespace
