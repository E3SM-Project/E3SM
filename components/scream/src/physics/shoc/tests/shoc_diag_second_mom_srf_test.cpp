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
struct UnitWrap::UnitTest<D>::TestSecondMomSrf {

static void run_second_mom_srf_bfb()
{
  SHOCSecondMomentSrfData mom_srf_data_f90[] = {
    //                      shcol
    SHOCSecondMomentSrfData(36),
    SHOCSecondMomentSrfData(72),
    SHOCSecondMomentSrfData(128),
    SHOCSecondMomentSrfData(256),
  };

  static constexpr Int num_runs = sizeof(mom_srf_data_f90) / sizeof(SHOCSecondMomentSrfData);

  for (Int i = 0; i < num_runs; ++i) {
    mom_srf_data_f90[i].randomize({}, {}, { {-1, 1} });
  }

  SHOCSecondMomentSrfData mom_srf_data_cxx[] = {
    SHOCSecondMomentSrfData(mom_srf_data_f90[0]),
    SHOCSecondMomentSrfData(mom_srf_data_f90[1]),
    SHOCSecondMomentSrfData(mom_srf_data_f90[2]),
    SHOCSecondMomentSrfData(mom_srf_data_f90[3]),
  };

  for (Int i = 0; i < num_runs; ++i) {
    // expects data in C layout
    shoc_diag_second_moments_srf(mom_srf_data_f90[i]);
  }

  for (Int i = 0; i < num_runs; ++i) {
    SHOCSecondMomentSrfData& d = mom_srf_data_cxx[i];
    shoc_diag_second_moments_srf_f(d.shcol, d.wthl, d.uw, d.vw, d.ustar2, d.wstar);
  }

  for (Int i = 0; i < num_runs; ++i) {
    Int shcol = mom_srf_data_cxx[i].shcol;
    for (Int k = 0; k < shcol; ++k) {
      REQUIRE(mom_srf_data_f90[i].ustar2[k] == mom_srf_data_cxx[i].ustar2[k]);
      REQUIRE(mom_srf_data_f90[i].wstar[k]  == mom_srf_data_cxx[i].wstar[k]);
    }
  }
}

static void run_second_mom_srf_phys()
{
  // TODO
}

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("shoc_second_moments_srf_property", "shoc") {
  using TRS = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSecondMomSrf;

  TRS::run_second_mom_srf_phys();
}

TEST_CASE("shoc_second_moments_srf_bfb", "shoc") {
  using TRS = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestSecondMomSrf;

  TRS::run_second_mom_srf_bfb();
}

}
