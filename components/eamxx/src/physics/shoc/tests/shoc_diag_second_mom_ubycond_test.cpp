#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "shoc_functions.hpp"
#include "shoc_test_data.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/eamxx_types.hpp"
#include "share/util/eamxx_setup_random_test.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

//#include "share/eamxx_types.hpp"
#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestSecondMomUbycond : public UnitWrap::UnitTest<D>::Base {

  void run_property()
  {
    // Property test for SHOC subroutine:
    //   diag_second_moments_ubycond

    static constexpr Int shcol    = 2;

    // Note that this subroutine does not have any inputs,
    //  only outputs.  All outputs should be zero.

    // Initialize data structure for bridging to F90
    DiagSecondMomentsUbycondData SDS(shcol);

    // Test that the inputs are reasonable
    REQUIRE(SDS.shcol == shcol);
    REQUIRE(shcol > 0);

    // Call the C++ implementation
    diag_second_moments_ubycond(SDS);

    // Verify the result
    //  all output should be zero.

    for (Int s = 0; s < shcol; ++s){
      REQUIRE(SDS.thl_sec[s] == 0);
      REQUIRE(SDS.qw_sec[s] == 0);
      REQUIRE(SDS.qwthl_sec[s] == 0);
      REQUIRE(SDS.wthl_sec[s] == 0);
      REQUIRE(SDS.wqw_sec[s] == 0);
      REQUIRE(SDS.uw_sec[s] == 0);
      REQUIRE(SDS.vw_sec[s] == 0);
      REQUIRE(SDS.wtke_sec[s] == 0);
    }
  }

  void run_bfb()
  {
    auto engine = Base::get_engine();

    DiagSecondMomentsUbycondData uby_fortran[] = {
      // shcol
      DiagSecondMomentsUbycondData(128),
      DiagSecondMomentsUbycondData(128),
      DiagSecondMomentsUbycondData(128),
      DiagSecondMomentsUbycondData(128),
    };

    static constexpr Int num_runs = sizeof(uby_fortran) / sizeof(DiagSecondMomentsUbycondData);

    for (auto& d : uby_fortran) {
      d.randomize(engine);
    }

    // Create copies of data for use by cxx. Needs to happen before reads so that
    // inout data is in original state
    DiagSecondMomentsUbycondData uby_cxx[num_runs] = {
      DiagSecondMomentsUbycondData(uby_fortran[0]),
      DiagSecondMomentsUbycondData(uby_fortran[1]),
      DiagSecondMomentsUbycondData(uby_fortran[2]),
      DiagSecondMomentsUbycondData(uby_fortran[3]),
    };

    // Read baseline data
    if (this->m_baseline_action == COMPARE) {
      for (auto& d : uby_fortran) {
        d.read(Base::m_fid);
      }
    }

    for (auto& d : uby_cxx) {
      diag_second_moments_ubycond(d);
    }

    if (SCREAM_BFB_TESTING && this->m_baseline_action == COMPARE) {
      for (Int i = 0; i < num_runs; ++i) {
        const Int shcol = uby_cxx[i].shcol;
        for (Int k = 0; k < shcol; ++k) {
          REQUIRE(uby_fortran[i].thl_sec[k]      == uby_cxx[i].thl_sec[k]);
          REQUIRE(uby_fortran[i].qw_sec[k]       == uby_cxx[i].qw_sec[k]);
          REQUIRE(uby_fortran[i].qwthl_sec[k]    == uby_cxx[i].qwthl_sec[k]);
          REQUIRE(uby_fortran[i].wthl_sec[k]     == uby_cxx[i].wthl_sec[k]);
          REQUIRE(uby_fortran[i].wqw_sec[k]      == uby_cxx[i].wqw_sec[k]);
          REQUIRE(uby_fortran[i].uw_sec[k]       == uby_cxx[i].uw_sec[k]);
          REQUIRE(uby_fortran[i].vw_sec[k]       == uby_cxx[i].vw_sec[k]);
          REQUIRE(uby_fortran[i].wtke_sec[k]     == uby_cxx[i].wtke_sec[k]);
        }
      }
    } // SCREAM_BFB_TESTING
    else if (this->m_baseline_action == GENERATE) {
      for (Int i = 0; i < num_runs; ++i) {
        uby_cxx[i].write(Base::m_fid);
      }
    }
  }

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("second_mom_uby_property", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSecondMomUbycond;

  TestStruct().run_property();
}

TEST_CASE("second_mom_uby_bfb", "shoc")
{
  using TestStruct = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSecondMomUbycond;

  TestStruct().run_bfb();
}

} // namespace
