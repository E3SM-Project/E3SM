#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "shoc_functions.hpp"
#include "shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"
#include "share/util/scream_setup_random_test.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestSecondMomSrf {

  static void run_property()
  {
    // Property test for the SHOC subroutine:
    //  diag_second_moments_srf

    static constexpr Int shcol    = 5;

    // Given several columns worth of data with varying conditions,
    //  verify that the output is as expected.

    // Define the surface heat flux [K m/s]
    static constexpr Real wthl_sfc[shcol] = {-0.03, 0.02, 0, -0.2, 0.04};
    // Define the zonal wind flux [m2/s2]
    static constexpr Real uw_sfc[shcol] = {-0.02, 0.2, -0.9, 0.3, 0};
    // Define the meridional wind flux [m2/s2]
    static constexpr Real vw_sfc[shcol] = {-0.4, 0.1, 0.01, -0.2, 0};

    // Initialize data structure for bridging to F90
    DiagSecondMomentsSrfData SDS(shcol);

    // Load up the data
    for (Int s = 0; s < shcol; ++s){
      SDS.wthl_sfc[s] = wthl_sfc[s];
      SDS.uw_sfc[s] = uw_sfc[s];
      SDS.vw_sfc[s] = vw_sfc[s];
    }

    // Call the C++ implementation
    shoc_diag_second_moments_srf_f(SDS.shcol, SDS.wthl_sfc, SDS.uw_sfc, SDS.vw_sfc, SDS.ustar2, SDS.wstar);

    // Verify the output
    for (Int s = 0; s < shcol; ++s){
      // if surface heat flux is less than or equal to zero then
      //  wstar should be zero, otherwise, it should be greater
      //  than zero.
      if (SDS.wthl_sfc[s] <= 0){
        REQUIRE(SDS.wstar[s] == 0);
      }
      else {
        REQUIRE(SDS.wstar[s] > 0);
        // Check result is within reasonable bounds
        REQUIRE(SDS.wstar[s] < 1);
      }

      // if both input surface wind fluxes are zero then
      //   ustar2 should be zero, else it should always be positive
      if (SDS.uw_sfc[s] == 0 && SDS.vw_sfc[s] == 0){
        REQUIRE(SDS.ustar2[s] == 0);
      }
      else{
        REQUIRE(SDS.ustar2[s] > 0);
        // Check result is within reasonable bounds
        REQUIRE(SDS.ustar2[s] < 10);
      }
    }

  }

  static void run_bfb()
  {
  #if 0
    auto engine = setup_random_test();

    SHOCSecondMomentSrfData mom_srf_data_f90[] = {
      //                      shcol
      SHOCSecondMomentSrfData(36),
      SHOCSecondMomentSrfData(72),
      SHOCSecondMomentSrfData(128),
      SHOCSecondMomentSrfData(256),
    };

    for (auto& d : mom_srf_data_f90) {
      d.randomize(engine, { {d.wthl, {-1, 1}} });
    }

    SHOCSecondMomentSrfData mom_srf_data_cxx[] = {
      SHOCSecondMomentSrfData(mom_srf_data_f90[0]),
      SHOCSecondMomentSrfData(mom_srf_data_f90[1]),
      SHOCSecondMomentSrfData(mom_srf_data_f90[2]),
      SHOCSecondMomentSrfData(mom_srf_data_f90[3]),
    };

    for (auto& d : mom_srf_data_f90) {
      // expects data in C layout
      shoc_diag_second_moments_srf(d);
    }

    for (auto& d : mom_srf_data_cxx) {
      shoc_diag_second_moments_srf_f(d.shcol, d.wthl_sfc, d.uw_sfc, d.vw_sfc, d.ustar2, d.wstar);
    }

    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(mom_srf_data_f90) / sizeof(SHOCSecondMomentSrfData);
      for (Int i = 0; i < num_runs; ++i) {
        Int shcol = mom_srf_data_cxx[i].shcol;
        for (Int k = 0; k < shcol; ++k) {
          REQUIRE(mom_srf_data_f90[i].ustar2[k] == mom_srf_data_cxx[i].ustar2[k]);
          REQUIRE(mom_srf_data_f90[i].wstar[k]  == mom_srf_data_cxx[i].wstar[k]);
        }
      }
    }
#endif
  }

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("shoc_second_moments_srf_property", "shoc") {
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSecondMomSrf;

  TestStruct::run_property();
}

TEST_CASE("shoc_second_moments_srf_bfb", "shoc") {
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSecondMomSrf;

  TestStruct::run_bfb();
}

} // namespace
