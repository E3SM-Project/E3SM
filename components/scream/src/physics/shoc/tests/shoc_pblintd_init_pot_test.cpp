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
struct UnitWrap::UnitTest<D>::TestPblintdInitPot {

static void run_pblintd_init_pot_bfb()
{
  PblintdInitPotData pblintd_init_pot_data_f90[] = {
    //                     shcol, nlev
    PblintdInitPotData(36,  72),
    PblintdInitPotData(72,  72),
    PblintdInitPotData(128, 72),
    PblintdInitPotData(256, 72),
  };

  static constexpr Int num_runs = sizeof(pblintd_init_pot_data_f90) / sizeof(PblintdInitPotData);

  for (auto& d : pblintd_init_pot_data_f90) {
    d.randomize();
  }

  PblintdInitPotData pblintd_init_pot_data_cxx[] = {
    PblintdInitPotData(pblintd_init_pot_data_f90[0]),
    PblintdInitPotData(pblintd_init_pot_data_f90[1]),
    PblintdInitPotData(pblintd_init_pot_data_f90[2]),
    PblintdInitPotData(pblintd_init_pot_data_f90[3]),
  };

  for (auto& d : pblintd_init_pot_data_f90) {
    // expects data in C layout
    pblintd_init_pot(d);
  }

  for (auto& d : pblintd_init_pot_data_cxx) {
    shoc_pblintd_init_pot_f(d.shcol, d.nlev, d.thl, d.ql, d.q, d.thv);
  }

  for (Int i = 0; i < num_runs; ++i) {
    const Int shcol = pblintd_init_pot_data_cxx[i].shcol;
    const Int nlev  = pblintd_init_pot_data_cxx[i].nlev;
    for (Int j = 0; j < shcol; ++j ) {
      for (Int k = 0; k < nlev; ++k) {
        REQUIRE(pblintd_init_pot_data_f90[i].thv[j*k] == pblintd_init_pot_data_cxx[i].thv[j*k]);
      }
    }
  }
}

static void run_pblintd_init_pot_phys()
{
  // TODO
}

};

} // namespace unit_test
} // namespace shoc
} // namespace scream

namespace {

TEST_CASE("shoc_pblintd_init_pot_property", "shoc") {
  using TRS = scream::shoc::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestPblintdInitPot;

  TRS::run_pblintd_init_pot_phys();
  TRS::run_pblintd_init_pot_bfb();
}
}
