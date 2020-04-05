#include "catch2/catch.hpp"

#include "share/scream_types.hpp"
#include "share/util/scream_utils.hpp"
#include "share/scream_kokkos.hpp"
#include "share/scream_pack.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/p3/p3_functions_f90.hpp"
#include "physics/p3/p3_f90.hpp"
#include "share/util/scream_kokkos_utils.hpp"

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
    const std::array< std::pair<Real, Real>, LatentHeatData::NUM_ARRAYS > ranges = {
      std::make_pair(4.056E-01, 1.153E+00), // v
      std::make_pair(1.000E+02, 5.000E+02), // s
      std::make_pair(1.000E+01, 2.000E+01), // f
    };

    LatentHeatData latent_fortran[] = {
     // its, ite, kts, kte, ranges
      LatentHeatData(1,  7,   1,  10,  ranges),
      LatentHeatData(1,  7,   1,  10,  ranges),
      LatentHeatData(1,  7,   1,  10,  ranges),
      LatentHeatData(1,  7,   1,  10,  ranges),
    };

    static constexpr Int num_runs = sizeof(latent_fortran) / sizeof(LatentHeatData);

    LatentHeatData latent_cxx[num_runs] = {
      LatentHeatData(latent_fortran[0]),
      LatentHeatData(latent_fortran[1]),
      LatentHeatData(latent_fortran[2]),
      LatentHeatData(latent_fortran[3]),
    };


    for (Int i = 0; i < num_runs; ++i) {
      get_latent_heat(latent_fortran[i]);
    }

    for (Int i = 0; i < num_runs; ++i) {
      LatentHeatData& d = latent_cxx[i];
      get_latent_heat_f(d.its, d.ite, d.kts, d.kte, d.v, d.s, d.f);

      LatentHeatData& h = latent_fortran[i];
      Int ni = d.ite - d.its + 1;
      Int nk = d.kte - d.kts + 1;
      Int n = 0;
      for (Int k = 0; k < nk; ++k) {
        for (Int j = 0; j < ni; ++j) {
printf("n= %d, d.s= %f, h.s= %f ",n,d.s[n],h.s[n]);
          REQUIRE(d.s[n] == h.s[n]);
          REQUIRE(d.v[n] == h.v[n]);
          REQUIRE(d.f[n] == h.f[n]);
          ++n;
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
