#include <catch2/catch.hpp>

#include <iostream>

#include "Elements.hpp"
#include "Tracers.hpp"
#include "Types.hpp"
#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"

using namespace Homme;

// ====================== RANDOM INITIALIZATION ====================== //

TEST_CASE("dp3d_intervals", "Testing Elements::random_init") {
  constexpr int num_elems = 5;
  constexpr Real max_pressure = 32.0;
  constexpr Real rel_threshold = 128.0 * std::numeric_limits<Real>::epsilon();
  Elements elements;
  elements.random_init(num_elems, max_pressure);
  HostViewManaged<Scalar * [NUM_TIME_LEVELS][NP][NP][NUM_LEV]> dp3d("host dp3d",
                                                                    num_elems);
  Kokkos::deep_copy(dp3d, elements.m_dp3d);
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
          Real rel_error = compare_answers(max_pressure, sum);
          REQUIRE(rel_threshold >= rel_error);
        }
      }
    }
  }
}

TEST_CASE("d_dinv_check", "Testing Elements::random_init") {
  constexpr int num_elems = 5;
  constexpr Real rel_threshold = 1e-10;
  Elements elements;
  elements.random_init(num_elems);
  HostViewManaged<Real * [2][2][NP][NP]> d("host d", num_elems);
  HostViewManaged<Real * [2][2][NP][NP]> dinv("host dinv", num_elems);
  Kokkos::deep_copy(d, elements.m_d);
  Kokkos::deep_copy(dinv, elements.m_dinv);
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

#if 0
TEST_CASE("tracers_check", "Testing Tracers::Tracers(int, int)") {
  // Ensures three things - that genRandArray results in an array starting and
  // ending with values between the specified bounds, that it does not exceed
  // the bounds specified, and that the tracers aren't accidentally overwriting
  // each other
  constexpr int num_elems = 3;
  constexpr int num_tracers = 5;
  constexpr Real min_val = 5.3, max_val = 9.3;
  constexpr Real signature = min_val - 1.0;
  std::random_device rd;
  std::mt19937_64 engine(rd());
  std::uniform_real_distribution<Real> dist(min_val, max_val);
  Tracers tracers(num_elems, num_tracers);
  for (int ie = 0; ie < num_elems; ++ie) {
    for (int iq = 0; iq < num_tracers; ++iq) {
      Tracers::Tracer t = tracers.tracer(ie, iq);
      genRandArray(t.qtens, engine, dist);
      REQUIRE(t.qtens(0, 0, 0)[0] >= min_val);
      REQUIRE(t.qtens(0, 0, 0)[0] <= max_val);
      REQUIRE(t.qtens(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] >= min_val);
      REQUIRE(t.qtens(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] <= max_val);
      t.qtens(0, 0, 0)[0] = signature;
      t.qtens(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] = signature;

      genRandArray(t.vstar_qdp, engine, dist);
      REQUIRE(t.vstar_qdp(0, 0, 0, 0)[0] >= min_val);
      REQUIRE(t.vstar_qdp(0, 0, 0, 0)[0] <= max_val);
      REQUIRE(t.vstar_qdp(1, NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] >=
              min_val);
      REQUIRE(t.vstar_qdp(1, NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] <=
              max_val);
      t.vstar_qdp(0, 0, 0, 0)[0] = signature;
      t.vstar_qdp(1, NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] = signature;

      genRandArray(t.qlim, engine, dist);
      REQUIRE(t.qlim(0, 0)[0] >= min_val);
      REQUIRE(t.qlim(0, 0)[0] <= max_val);
      REQUIRE(t.qlim(1, NUM_LEV - 1)[VECTOR_SIZE - 1] >= min_val);
      REQUIRE(t.qlim(1, NUM_LEV - 1)[VECTOR_SIZE - 1] <= max_val);
      t.qlim(0, 0)[0] = signature;
      t.qlim(1, NUM_LEV - 1)[VECTOR_SIZE - 1] = signature;

      genRandArray(t.q, engine, dist);
      REQUIRE(t.q(0, 0, 0)[0] >= min_val);
      REQUIRE(t.q(0, 0, 0)[0] <= max_val);
      REQUIRE(t.q(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] >= min_val);
      REQUIRE(t.q(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] <= max_val);
      t.q(0, 0, 0)[0] = signature;
      t.q(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] = signature;

      genRandArray(t.qtens_biharmonic, engine, dist);
      REQUIRE(t.qtens_biharmonic(0, 0, 0)[0] >= min_val);
      REQUIRE(t.qtens_biharmonic(0, 0, 0)[0] <= max_val);
      REQUIRE(t.qtens_biharmonic(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] >= min_val);
      REQUIRE(t.qtens_biharmonic(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] <= max_val);
      t.qtens_biharmonic(0, 0, 0)[0] = signature;
      t.qtens_biharmonic(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] = signature;
    }
  }
  for (int ie = 0; ie < num_elems; ++ie) {
    for (int iq = 0; iq < num_tracers; ++iq) {
      Tracers::Tracer t = tracers.tracer(ie, iq);
      REQUIRE(t.qtens(0, 0, 0)[0] == signature);
      REQUIRE(t.qtens(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] ==
              signature);
      REQUIRE(t.vstar_qdp(0, 0, 0, 0)[0] == signature);
      REQUIRE(t.vstar_qdp(1, NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] ==
              signature);
      REQUIRE(t.qlim(0, 0)[0] == signature);
      REQUIRE(t.qlim(1, NUM_LEV - 1)[VECTOR_SIZE - 1] == signature);
      REQUIRE(t.q(0, 0, 0)[0] == signature);
      REQUIRE(t.q(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] == signature);
      REQUIRE(t.qtens_biharmonic(0, 0, 0)[0] == signature);
      REQUIRE(t.qtens_biharmonic(NP - 1, NP - 1, NUM_LEV - 1)[VECTOR_SIZE - 1] == signature);
    }
  }
}
#endif
