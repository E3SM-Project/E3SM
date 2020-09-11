#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"

//#include "share/scream_types.hpp"
#include <algorithm>
#include <array>
#include <random>
#include <thread>

namespace scream {
namespace shoc {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestSecondMomUbycond {

static void run_second_mom_ubycond_bfb()
{
  SHOCSecondMomentUbycondData uby_fortran[] = {
    // shcol
    SHOCSecondMomentUbycondData(128),
    SHOCSecondMomentUbycondData(128),
    SHOCSecondMomentUbycondData(128),
    SHOCSecondMomentUbycondData(128),
  };

  static constexpr Int num_runs = sizeof(uby_fortran) / sizeof(SHOCSecondMomentUbycondData);

  for (auto& d : uby_fortran) {
    d.randomize();
  }

  // Create copies of data for use by cxx. Needs to happen before fortran calls so that
  // inout data is in original state
  SHOCSecondMomentUbycondData uby_cxx[num_runs] = {
    SHOCSecondMomentUbycondData(uby_fortran[0]),
    SHOCSecondMomentUbycondData(uby_fortran[1]),
    SHOCSecondMomentUbycondData(uby_fortran[2]),
    SHOCSecondMomentUbycondData(uby_fortran[3]),
  };

  // Get data from fortran
  for (auto& d : uby_fortran) {
    shoc_diag_second_moments_ubycond(d);
  }

  for (auto& d : uby_cxx) {
    shoc_diag_second_moments_ubycond_f(d.shcol(), d.thl, d.qw, d.qwthl, d.wthl, d.wqw, d.uw, d.vw, d.wtke);
  }


  for (Int i = 0; i < num_runs; ++i) {
    Int shcol = uby_cxx[i].shcol();
    for (Int k = 0; k < shcol; ++k) {
      REQUIRE(uby_fortran[i].thl[k]      == uby_cxx[i].thl[k]);
      REQUIRE(uby_fortran[i].qw[k]       == uby_cxx[i].qw[k]);
      REQUIRE(uby_fortran[i].qwthl[k]    == uby_cxx[i].qwthl[k]);
      REQUIRE(uby_fortran[i].wthl[k]     == uby_cxx[i].wthl[k]);
      REQUIRE(uby_fortran[i].wqw[k]      == uby_cxx[i].wqw[k]);
      REQUIRE(uby_fortran[i].uw[k]       == uby_cxx[i].uw[k]);
      REQUIRE(uby_fortran[i].vw[k]       == uby_cxx[i].vw[k]);
      REQUIRE(uby_fortran[i].wtke[k]     == uby_cxx[i].wtke[k]);
    }
  }
}

static void run_second_mom_ubycond_phys()
{
    // TODO
}
};

}
}
}

namespace {
TEST_CASE("second_mom_uby", "shoc") {
  using TRS = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSecondMomUbycond;
  TRS::run_second_mom_ubycond_phys();
  TRS::run_second_mom_ubycond_bfb();
}
}
