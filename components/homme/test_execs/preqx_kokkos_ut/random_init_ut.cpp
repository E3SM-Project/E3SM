#include <catch2/catch.hpp>

#include "ElementsState.hpp"
#include "ElementsGeometry.hpp"
#include "PhysicalConstants.hpp"
#include "Types.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"

#include <iostream>
#include <random>

using namespace Homme;

// ====================== RANDOM INITIALIZATION ====================== //

TEST_CASE("dp3d_intervals", "Testing Elements::random_init") {
  constexpr int num_elems = 5;
  constexpr Real max_pressure = 32.0;
  constexpr Real ps0 = 32.0/1000.0;
  constexpr Real hyai0 = 0.001;
  constexpr Real rel_threshold = 128.0 * std::numeric_limits<Real>::epsilon();
  ElementsState state;

  std::random_device rd;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  state.init(num_elems);
  state.randomize(seed, max_pressure, ps0, hyai0);
  HostViewManaged<Scalar * [NUM_TIME_LEVELS][NP][NP][NUM_LEV]> dp3d("host dp3d",
                                                                    num_elems);
  Kokkos::deep_copy(dp3d, state.m_dp3d);
  for (int ie = 0; ie < num_elems; ++ie) {
    for (int tl = 0; tl < NUM_TIME_LEVELS; ++tl) {
      for (int igp = 0; igp < NP; ++igp) {
        for (int jgp = 0; jgp < NP; ++jgp) {
          HostViewUnmanaged<Scalar[NUM_LEV]> dp3d_col =
              Homme::subview(dp3d, ie, tl, igp, jgp);
          Real sum = 0.0;
          for (int level = 0; level < NUM_PHYSICAL_LEV; ++level) {
            const int ilev = level / VECTOR_SIZE;
            const int vec_lev = level % VECTOR_SIZE;
            sum += dp3d_col(ilev)[vec_lev];
            REQUIRE(dp3d_col(ilev)[vec_lev] > 0.0);
          }
          Real rel_error = compare_answers(max_pressure, sum+ps0);
          REQUIRE(rel_threshold >= rel_error);
        }
      }
    }
  }
}

TEST_CASE("d_dinv_check", "Testing Elements::random_init") {
  constexpr int num_elems = 5;
  constexpr Real rel_threshold = 1e-10;
  std::random_device rd;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");

  ElementsGeometry geometry;
  geometry.init(num_elems,false,true,PhysicalConstants::rearth0);
  geometry.randomize(seed);

  HostViewManaged<Real * [2][2][NP][NP]> d("host d", num_elems);
  HostViewManaged<Real * [2][2][NP][NP]> dinv("host dinv", num_elems);
  Kokkos::deep_copy(d, geometry.m_d);
  Kokkos::deep_copy(dinv, geometry.m_dinv);
  for (int ie = 0; ie < num_elems; ++ie) {
    for (int igp = 0; igp < NP; ++igp) {
      for (int jgp = 0; jgp < NP; ++jgp) {
        const Real det_1 = d(ie, 0, 0, igp, jgp) * d(ie, 1, 1, igp, jgp) -
                           d(ie, 0, 1, igp, jgp) * d(ie, 1, 0, igp, jgp);
        REQUIRE(det_1 > 0.0);
        const Real det_2 = dinv(ie, 0, 0, igp, jgp) * dinv(ie, 1, 1, igp, jgp) -
                           dinv(ie, 0, 1, igp, jgp) * dinv(ie, 1, 0, igp, jgp);
        REQUIRE(det_2 > 0.0);
        for (int i = 0; i < 2; ++i) {
          for (int j = 0; j < 2; ++j) {
            Real pt_product = 0.0;
            for (int k = 0; k < 2; ++k) {
              pt_product += d(ie, i, k, igp, jgp) * dinv(ie, k, j, igp, jgp);
            }
            const Real expected = (i == j) ? 1.0 : 0.0;
            const Real rel_error = compare_answers(expected, pt_product);
            REQUIRE(rel_threshold >= rel_error);
          }
        }
      }
    }
  }
}
