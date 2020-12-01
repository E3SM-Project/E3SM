#include "catch2/catch.hpp"

#include "shoc_unit_tests_common.hpp"

#include "physics/shoc/shoc_functions.hpp"
#include "physics/shoc/shoc_functions_f90.hpp"
#include "physics/share/physics_constants.hpp"
#include "share/scream_types.hpp"

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

static void run_pblintd_cldcheck_bfb()
{
  PblintdCldcheckData cldcheck_data_f90[] = {
    //                      shcol, nlev, nlevi
    PblintdCldcheckData(36,  128, 129),
    PblintdCldcheckData(72,  128, 129),
    PblintdCldcheckData(128, 128, 129),
    PblintdCldcheckData(256, 128, 129),
  };

  static constexpr Int num_runs = sizeof(cldcheck_data_f90) / sizeof(PblintdCldcheckData);

  for (auto& d : cldcheck_data_f90) {
    d.randomize();
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

  for (Int i = 0; i < num_runs; ++i) {
    const Int shcol = cldcheck_data_cxx[i].shcol;
    for (Int k = 0; k < shcol; ++k) {
      REQUIRE(cldcheck_data_f90[i].pblh[k]  == cldcheck_data_cxx[i].pblh[k]);
    }
  }
}

static void run_pblintd_cldcheck_phys()
{
  // TODO
}

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("shoc_pblintd_cldcheck", "shoc") {
  using TRS = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdCldCheck;

  TRS::run_pblintd_cldcheck_phys();
  TRS::run_pblintd_cldcheck_bfb();
}
}
