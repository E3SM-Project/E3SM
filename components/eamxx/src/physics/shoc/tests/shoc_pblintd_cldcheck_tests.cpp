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
struct UnitWrap::UnitTest<D>::TestPblintdCldCheck {

  static void run_property()
  {
    static constexpr Int shcol = 5;
    static constexpr Int nlev = 3;
    static constexpr Int nlevi = nlev+1;

    // Tests for the subroutine pblintd_cldcheck

    // Define interface height [m]
    static constexpr Real zi[nlevi] = {1500, 1000, 500, 0};
    // Define cloud fraction [-] (surface value will be modified later)
    static constexpr Real cldn[nlev] = {0.5, 0.4, 0.0};
    // PBL height [m] (feed small value for all columns)
    static constexpr Real pblh = 10;

    // Initialize data structure for bridging to F90
    PblintdCldcheckData SDS(shcol, nlev, nlevi);

    // Test that the inputs are reasonable
    REQUIRE( (SDS.shcol == shcol && SDS.nlev == nlev && SDS.nlevi == nlevi) );
    REQUIRE(shcol > 0);
    REQUIRE(nlevi == nlev+1);

    // Fill in test data
    for(Int s = 0; s < shcol; ++s) {
      SDS.pblh[s] = pblh;
      for(Int n = 0; n < nlev; ++n) {
        const auto offset = n + s * nlev;
        SDS.cldn[offset] = cldn[n];
        // Feed lowest level 0 for first column and non
        //  zero for subsequent columns
        if (n == nlev-1){
          SDS.cldn[offset] = std::min(1.0,s*0.1);
        }
      }
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
        SDS.zi[offset] = zi[n];
      }
    }

    // check to make sure the input data make sense
    for(Int s = 0; s < shcol; ++s) {
      REQUIRE(SDS.pblh[s] > 0);
      for(Int n = 0; n < nlev; ++n){
        const auto offset = n + s * nlev;
        REQUIRE(SDS.cldn[offset] >= 0);
        REQUIRE(SDS.cldn[offset] <= 1);
      }
      for(Int n = 0; n < nlevi; ++n) {
        const auto offset = n + s * nlevi;
        REQUIRE(SDS.zi[offset] >= 0);
      }
    }

    // Call the fortran implementation
    pblintd_cldcheck(SDS);

    // Check the result
    for(Int s = 0; s < shcol; ++s) {
      const auto lowlev = (nlev-1) + s * nlev;
      // If surface cloud exists verify that pbl value was adjusted
      if (SDS.cldn[lowlev] > 0){
        REQUIRE(SDS.pblh[s] > pblh);
      }
    }

  } // run_property

  static void run_bfb()
  {
    auto engine = setup_random_test();

    PblintdCldcheckData cldcheck_data_f90[] = {
      //                      shcol, nlev, nlevi
      PblintdCldcheckData(36,  128, 129),
      PblintdCldcheckData(72,  128, 129),
      PblintdCldcheckData(128, 128, 129),
      PblintdCldcheckData(256, 128, 129),
    };

    for (auto& d : cldcheck_data_f90) {
      d.randomize(engine);
    }

    PblintdCldcheckData cldcheck_data_cxx[] = {
      PblintdCldcheckData(cldcheck_data_f90[0]),
      PblintdCldcheckData(cldcheck_data_f90[1]),
      PblintdCldcheckData(cldcheck_data_f90[2]),
      PblintdCldcheckData(cldcheck_data_f90[3]),
    };

    for (auto& d : cldcheck_data_f90) {
      // expects data in C layout
      pblintd_cldcheck(d);
    }

    for (auto& d : cldcheck_data_cxx) {
      d.transpose<ekat::TransposeDirection::c2f>();
      shoc_pblintd_cldcheck_f(d.shcol, d.nlev, d.nlevi, d.zi, d.cldn, d.pblh);
    }

    if (SCREAM_BFB_TESTING) {
      static constexpr Int num_runs = sizeof(cldcheck_data_f90) / sizeof(PblintdCldcheckData);
      for (Int i = 0; i < num_runs; ++i) {
        const Int shcol = cldcheck_data_cxx[i].shcol;
        for (Int k = 0; k < shcol; ++k) {
          REQUIRE(cldcheck_data_f90[i].pblh[k]  == cldcheck_data_cxx[i].pblh[k]);
        }
      }
    }
  }  // run_bfb

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("shoc_pblintd_cldcheck_property", "shoc") {
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdCldCheck;

  TestStruct::run_property();

}

TEST_CASE("shoc_pblintd_cldcheck_bfb", "shoc") {
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdCldCheck;

  TestStruct::run_bfb();

}

}  // empty namespace
