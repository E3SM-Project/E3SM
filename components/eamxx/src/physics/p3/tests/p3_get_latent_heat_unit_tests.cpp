#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "ekat/ekat_pack.hpp"
#include "ekat/kokkos/ekat_kokkos_utils.hpp"
#include "p3_functions.hpp"
#include "p3_functions_f90.hpp"
#include "p3_f90.hpp"

#include "p3_unit_tests_common.hpp"

#include <thread>
#include <array>
#include <algorithm>
#include <random>
#include <iomanip>      // std::setprecision

namespace scream {
namespace p3 {
namespace unit_test {

template <typename D>
struct UnitWrap::UnitTest<D>::TestLatentHeat {

  static void run_latent_heat_bfb()
  {
    LatentHeatData latent_fortran[] = {
      //           its, ite, kts, kte
      LatentHeatData(1,   7,   1,  10),
      LatentHeatData(1,   7,   1,  10),
      LatentHeatData(1,   7,   1,  10),
      LatentHeatData(1,   7,   1,  10),
    };

    static constexpr Int num_runs = sizeof(latent_fortran) / sizeof(LatentHeatData);

    LatentHeatData latent_cxx[num_runs] = {
      LatentHeatData(latent_fortran[0]),
      LatentHeatData(latent_fortran[1]),
      LatentHeatData(latent_fortran[2]),
      LatentHeatData(latent_fortran[3]),
    };

    for (Int i = 0; i < num_runs; ++i) {
      LatentHeatData& h = latent_fortran[i];
      get_latent_heat(h);

      LatentHeatData& d = latent_cxx[i];
      get_latent_heat_f(d.its, d.ite, d.kts, d.kte, d.v, d.s, d.f);

      if (SCREAM_BFB_TESTING) {
        REQUIRE(h.total(h.v) == d.total(d.v));

        for (Int j = 0; j < h.total(h.v); ++j) {
          REQUIRE(d.v[j] == h.v[j]);
          REQUIRE(d.s[j] == h.s[j]);
          REQUIRE(d.f[j] == h.f[j]);
        }
      }
    }
  }

  static void run_latent_heat_phys()
  {
    // TODO
  }
};

}
}
}

namespace {

TEST_CASE("p3_latent_heat", "[p3_functions]")
{
  using TD = scream::p3::unit_test::UnitWrap::UnitTest<scream::DefaultDevice>::TestLatentHeat;

  TD::run_latent_heat_phys();
  TD::run_latent_heat_bfb();
}

}
