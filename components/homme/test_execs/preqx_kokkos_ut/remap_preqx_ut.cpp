#include <catch2/catch.hpp>

#include <limits>

#include "RemapFunctor.hpp"
#include "PpmRemap.hpp"
#include "PhysicalConstants.hpp"

#include "utilities/SubviewUtils.hpp"
#include "utilities/TestUtils.hpp"

#include <random>

TEST_CASE("remap_interface", "vertical remap") {

  using namespace Homme;
  using rngAlg = std::mt19937_64;
  using namespace Remap;
  using namespace Remap::Ppm;

  constexpr int num_elems = 4;
  std::random_device rd;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");

  Elements elements;
  elements.init(num_elems,seed, /*alloc_gradphis = */ false,
                PhysicalConstants::rearth0);
  elements.randomize(seed);

  Tracers tracers;
  tracers.init(num_elems,QSIZE_D);
  tracers.randomize(seed+1);

  HybridVCoord hvcoord;
  hvcoord.random_init(seed+2);

  Real dp3d_min;
  {
    auto h_dp3d = Kokkos::create_mirror_view(elements.m_state.m_dp3d);
    Kokkos::deep_copy(h_dp3d,elements.m_state.m_dp3d);
    Real* ptr = reinterpret_cast<Real*>(h_dp3d.data());
    dp3d_min = *std::min_element(ptr,ptr+h_dp3d.size()*VECTOR_SIZE);
  }

  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]> eta_dot_dpdn ("",num_elems);
  ExecViewManaged<Scalar*[Q_NUM_TIME_LEVELS][QSIZE_D][NP][NP][NUM_LEV]> qdp("",num_elems);

  // TODO: make dt random
  constexpr int np1 = 0;
  constexpr int n0_qdp = 0;
  // Note: the bounds on the distribution for dt are strictly linked to how ps_v is (randomly)
  //       init-ed in ElementsState, and how eta_dot_dpdn is (randomly) init-ed in ElementsDerivedState
  //       (here we are initing eta_dot_dpdn in the same way). In particular, this interval *should* ensure that
  //       dp3d[k] + dt*(eta_dot_dpdn[k+1]-eta_dot_dpdn[k]) > 0, which is needed to pass the test
  std::uniform_real_distribution<Real> eta_pdf(0.01*dp3d_min, 0.1*dp3d_min);
  std::uniform_real_distribution<Real> dt_pdf(0.01, 10);
  rngAlg engine(seed+3);
  const Real dt = dt_pdf(engine);
  genRandArray(eta_dot_dpdn,engine,eta_pdf);
  genRandArray(qdp,engine,dt_pdf);

  SECTION("states_only") {
    constexpr bool rsplit_non_zero = true;
    constexpr int qsize = 0;
    using RF = RemapFunctor<rsplit_non_zero, PpmVertRemap<PpmMirrored>>;
    RF remap(qsize, elements, tracers, hvcoord);
    REQUIRE_NOTHROW(remap.run_remap(np1, n0_qdp, dt));
  }
  SECTION("tracers_only") {
    constexpr bool rsplit_non_zero = false;
    constexpr int qsize = QSIZE_D;
    using RF = RemapFunctor<rsplit_non_zero, PpmVertRemap<PpmMirrored>>;
    RF remap(qsize, elements, tracers, hvcoord);
    REQUIRE_NOTHROW(remap.run_remap(np1, n0_qdp, dt));
  }
  SECTION("states_tracers") {
    constexpr bool rsplit_non_zero = true;
    constexpr int qsize = QSIZE_D;
    using RF = RemapFunctor<rsplit_non_zero, PpmVertRemap<PpmMirrored>>;
    RF remap(qsize, elements, tracers, hvcoord);
    REQUIRE_NOTHROW(remap.run_remap(np1, n0_qdp, dt));
  }
}
