#include "catch2/catch.hpp"

#include "tms_unit_tests_common.hpp"

#include "share/eamxx_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "tms_functions.hpp"
#include "tms_functions_f90.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

namespace scream {
namespace tms {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestComputeTMS {

  static void run_property()
  {
    // Should property tests be created?
  } // run_property

  static void run_bfb()
  {
    auto engine = setup_random_test();

    ComputeTMSData f90_data[] = {
      //             ncols, nlevs
      ComputeTMSData(12,    72),
      ComputeTMSData(8,     12),
      ComputeTMSData(7,     16),
      ComputeTMSData(2,     7)
    };

    // Generate random input data
    for (auto& d : f90_data) {
      d.randomize(engine, { {d.sgh, {0.5, 1.5}} });
    }

    // Create copies of data for use by cxx. Needs to happen before fortran calls so that
    // inout data is in original state
    ComputeTMSData cxx_data[] = {
      ComputeTMSData(f90_data[0]),
      ComputeTMSData(f90_data[1]),
      ComputeTMSData(f90_data[2]),
      ComputeTMSData(f90_data[3])
    };

    // Assume all data is in C layout

    // Get data from fortran
    for (auto& d : f90_data) {
      // expects data in C layout
      compute_tms(d);
    }

    // Get data from cxx
    for (auto& d : cxx_data) {
      d.transpose<ekat::TransposeDirection::c2f>(); // _f expects data in fortran layout
      compute_tms_f(d.ncols, d.nlevs,
                    d.u_wind, d.v_wind, d.t_mid, d.p_mid, d.exner,
                    d.z_mid, d.sgh, d.landfrac, d.ksrf, d.taux, d.tauy);
      d.transpose<ekat::TransposeDirection::f2c>(); // go back to C layout
    }

    // Verify BFB results, all data should be in C layout
    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(f90_data) / sizeof(ComputeTMSData);

      for (int r = 0; r<num_runs; ++r) {
        ComputeTMSData& d_f90 = f90_data[r];
        ComputeTMSData& d_cxx = cxx_data[r];
        REQUIRE(d_f90.total(d_f90.ksrf) == d_cxx.total(d_cxx.ksrf));
        REQUIRE(d_f90.total(d_f90.taux) == d_cxx.total(d_cxx.taux));
        REQUIRE(d_f90.total(d_f90.tauy) == d_cxx.total(d_cxx.tauy));

        for (int i=0; i<d_f90.total(d_f90.ksrf); ++i) {
          REQUIRE(d_f90.ksrf[i] == d_cxx.ksrf[i]);
          REQUIRE(d_f90.taux[i] == d_cxx.taux[i]);
          REQUIRE(d_f90.tauy[i] == d_cxx.tauy[i]);
        }
      }
    }
  } // run_bfb
};

} // namespace unit_test
} // namespace tms
} // namespace scream

namespace {

//TEST_CASE("compute_tms_property", "tms")
//{
//  using TestStruct = scream::tms::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeTMS;
//  TestStruct::run_property();
//}

TEST_CASE("compute_tms_bfb", "tms")
{
  using TestStruct = scream::tms::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestComputeTMS;
  TestStruct::run_bfb();
}

} // empty namespace
