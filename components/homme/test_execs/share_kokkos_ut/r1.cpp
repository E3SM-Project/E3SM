#include <catch2/catch.hpp>

#include <limits>

#include "RemapFunctor.hpp"
#include "PpmRemap.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/TestUtils.hpp"

#include <random>

using namespace Homme;
using namespace Remap;
using namespace Ppm;
  
using rngAlg = std::mt19937_64;

/* This object is meant for testing different configurations of the PPM vertical
 * remap method.
 * boundary_cond needs to be one of the PPM boundary condition objects,
 * which provide the indexes to loop over
 */
template <typename boundary_cond>
class ppm_remap_functor_test {
  static_assert(std::is_base_of<PpmBoundaryConditions, boundary_cond>::value,
                "PPM Remap test must have a supported boundary condition");

public:
  ppm_remap_functor_test(const int num_elems, const int num_remap)
      : ne(num_elems), num_remap(num_remap), remap(num_elems, num_remap),
        src_layer_thickness_kokkos("source layer thickness", num_elems),
        tgt_layer_thickness_kokkos("target layer thickness", num_elems),
        remap_vals("values to remap", num_elems, num_remap) {}

  struct TagGridTest {};
  struct TagPPMTest {};
  struct TagRemapTest {};

  static bool nan_boundaries(
      HostViewUnmanaged<Real * [NP][NP][_ppm_consts::DPO_PHYSICAL_LEV]> host) {
    for (int ie = 0; ie < host.extent_int(0); ++ie) {
      for (int igp = 0; igp < host.extent_int(1); ++igp) {
        for (int jgp = 0; jgp < host.extent_int(2); ++jgp) {
          for (int k = 0; k < _ppm_consts::INITIAL_PADDING - _ppm_consts::gs;
               ++k) {
            host(ie, igp, jgp, k) = std::numeric_limits<Real>::quiet_NaN();
          }
          for (int k = _ppm_consts::INITIAL_PADDING + NUM_PHYSICAL_LEV +
                       _ppm_consts::gs;
               k < _ppm_consts::DPO_PHYSICAL_LEV; ++k) {
            host(ie, igp, jgp, k) = std::numeric_limits<Real>::quiet_NaN();
          }
        }
      }
    }
    return true;
  }

  void test_grid() {
    std::random_device rd;
    const unsigned int catchRngSeed = Catch::rngSeed();
    const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
    std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
    rngAlg engine(seed);
    genRandArray(remap.m_dpo, engine,
                 std::uniform_real_distribution<Real>(0.125, 1000),
                 nan_boundaries);
    Kokkos::parallel_for(
        Homme::get_default_team_policy<ExecSpace, TagGridTest>(ne), *this);
    Kokkos::fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagGridTest &, const TeamMember& team) const {
    KernelVariables kv(team);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      remap.compute_grids(kv, Homme::subview(remap.m_dpo, kv.ie, igp, jgp),
                          Homme::subview(remap.m_ppmdx, kv.ie, igp, jgp));
    });
  }

  void test_ppm() {
    // Hack to get the right number of outputs
    remap.m_ao = ExecViewManaged<Real * [NP][NP][_ppm_consts::AO_PHYSICAL_LEV]>(
        "ao", num_remap * ne);
    remap.m_dma = ExecViewManaged<Real * [NP][NP][_ppm_consts::DMA_PHYSICAL_LEV]>(
        "dma", num_remap * ne);
    remap.m_ai = ExecViewManaged<Real * [NP][NP][_ppm_consts::AI_PHYSICAL_LEV]>(
        "ai", num_remap * ne);
    remap.m_parabola_coeffs = ExecViewManaged<Real * [NP][NP][3][NUM_PHYSICAL_LEV]>(
        "parabola coeffs", num_remap * ne);

    std::random_device rd;
    const unsigned int catchRngSeed = Catch::rngSeed();
    const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
    std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
    rngAlg engine(seed);
    genRandArray(remap.m_ppmdx, engine,
                 std::uniform_real_distribution<Real>(0.125, 1000));
    for (int i = 0; i < num_remap; ++i) {
      genRandArray(remap.m_ao, engine,
                   std::uniform_real_distribution<Real>(0.125, 1000),
                   nan_boundaries);
    }
    Kokkos::parallel_for(
        Homme::get_default_team_policy<ExecSpace, TagPPMTest>(ne), *this);
    Kokkos::fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagPPMTest &, const TeamMember& team) const {
    KernelVariables kv(team);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, num_remap * NP * NP),
                         [&](const int &loop_idx) {
      const int var = loop_idx / NP / NP;
      const int igp = (loop_idx / NP) % NP;
      const int jgp = loop_idx % NP;
      remap.compute_ppm(
          kv,
          Homme::subview(remap.m_ao, kv.ie * num_remap + var, igp, jgp),
          Homme::subview(remap.m_ppmdx, kv.ie, igp, jgp),
          Homme::subview(remap.m_dma, kv.ie * num_remap + var, igp, jgp),
          Homme::subview(remap.m_ai, kv.ie * num_remap + var, igp, jgp),
          Homme::subview(remap.m_parabola_coeffs, kv.ie * num_remap + var, igp, jgp));
    });
  }

  // This ensures that the represented grid is not degenerate
  // The grid is represented by the interval widths, so they must all be
  // positive. The top of the grid must also be a fixed value,
  // so the sum of the intervals must be top, and the bottom is assumed to be 0
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]>
  generate_grid_intervals(rngAlg &engine, const Real &top, std::string name) {
    HostViewManaged<Scalar * [NP][NP][NUM_LEV]> grid("grid", ne);
    genRandArray(grid, engine,
                 std::uniform_real_distribution<Real>(0.0625, top));
    constexpr int last_vector = (NUM_PHYSICAL_LEV - 1) % VECTOR_SIZE;
    for (int ie = 0; ie < ne; ++ie) {
      for (int igp = 0; igp < NP; ++igp) {
        for (int jgp = 0; jgp < NP; ++jgp) {
          auto grid_slice = Homme::subview(grid, ie, igp, jgp);
          auto start = reinterpret_cast<Real *>(grid_slice.data());
          auto end = start + NUM_PHYSICAL_LEV;
          std::sort(start, end);
          grid_slice(NUM_LEV - 1)[last_vector] = top;

          // Changing grid[i] from absolute value to incremental value
          // compared to grid[i-1]. Note: there should be NPL+1 abs
          // values and NPL incremental values. We only generated NPL
          // abs values. We implicitly assume that there is an extra
          // grid value of 0 below the 1st generated one, so that
          // grid_slice(0)[0] ends up being the increment between the
          // first (the 0 which we never generated) and second level
          for (int k = NUM_PHYSICAL_LEV - 1; k > 0; --k) {
            const int vector_level = k / VECTOR_SIZE;
            const int vector = k % VECTOR_SIZE;
            const int lower_vector_level = (k - 1) / VECTOR_SIZE;
            const int lower_vector = (k - 1) % VECTOR_SIZE;
            grid_slice(vector_level)[vector] -=
                grid_slice(lower_vector_level)[lower_vector];
          }
        }
      }
    }

    ExecViewManaged<Scalar * [NP][NP][NUM_LEV]> intervals(name, ne);
    Kokkos::deep_copy(intervals, grid);
    return intervals;
  }

  void initialize_layers(rngAlg &engine) {
    // Note that these must have the property that
    // sum(src_layer_thickness) = sum(tgt_layer_thickness)
    // To do this, we generate two grids and compute the interval lengths

    const Real top = std::uniform_real_distribution<Real>(1.0, 1024.0)(engine);

    src_layer_thickness_kokkos =
        generate_grid_intervals(engine, top, "kokkos source layer thickness");
    tgt_layer_thickness_kokkos =
        generate_grid_intervals(engine, top, "kokkos target layer thickness");
  }

  void test_remap() {
    std::random_device rd;
    const unsigned int catchRngSeed = Catch::rngSeed();
    const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
    std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
    rngAlg engine(seed);
    std::uniform_real_distribution<Real> dist(0.125, 1000.0);
    genRandArray(remap_vals, engine, dist);

    REQUIRE(remap_vals.extent_int(0) == ne);
    REQUIRE(remap_vals.extent_int(1) == num_remap);
    REQUIRE(remap_vals.extent_int(2) == NP);
    REQUIRE(remap_vals.extent_int(3) == NP);
    REQUIRE(remap_vals.extent_int(4) == NUM_LEV);
    ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]> tmp =
        Homme::subview(remap_vals, 0, 0);
    REQUIRE(tmp.extent_int(0) == NP);
    REQUIRE(tmp.extent_int(1) == NP);
    REQUIRE(tmp.extent_int(2) == NUM_LEV);

    initialize_layers(engine);

    Kokkos::parallel_for(
        Homme::get_default_team_policy<ExecSpace, TagRemapTest>(ne), *this);
    Kokkos::fence();
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagRemapTest &, const TeamMember& team) const {
    KernelVariables kv(team);
    remap.compute_grids_phase(
        kv, Homme::subview(src_layer_thickness_kokkos, kv.ie),
        Homme::subview(tgt_layer_thickness_kokkos, kv.ie));
    for (int var = 0; var < num_remap; ++var) {
      remap.compute_remap_phase(kv, Homme::subview(remap_vals, kv.ie, var));
    }
  }

  const int ne, num_remap;
  PpmVertRemap<boundary_cond> remap;
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]> src_layer_thickness_kokkos;
  ExecViewManaged<Scalar * [NP][NP][NUM_LEV]> tgt_layer_thickness_kokkos;
  ExecViewManaged<Scalar * * [NP][NP][NUM_LEV]> remap_vals;
};

TEST_CASE("ppm_mirrored", "vertical remap") {
  constexpr int num_elems = 2;
  constexpr int num_remap = 3;
  ppm_remap_functor_test<PpmMirrored> remap_test_mirrored(num_elems, num_remap);
  SECTION("grid") { remap_test_mirrored.test_grid(); }
  SECTION("ppm") { remap_test_mirrored.test_ppm(); }
  SECTION("remap") { remap_test_mirrored.test_remap(); }
}


TEST_CASE("binary_search","binary_search")
{
  constexpr int length = 100;
  constexpr int num_tests = 5;

  std::random_device rd;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  rngAlg engine(seed);
  std::uniform_real_distribution<Real> pdf(0.01, 10);

  // Testing the binary search utility
  HostViewManaged<Real[length+1]> pin("");
  HostViewManaged<Real[length+2]> pio("");

  for (int itest=0; itest<num_tests; ++itest) {
    // Generate random pin, but sorted
    genRandArray(pin,engine,pdf);
    std::sort(pin.data(),pin.data()+pin.size());
    pin(0) = 0;

    // Generate random pio, but sorted
    genRandArray(pio,engine,pdf);
    std::sort(pio.data(),pio.data()+pio.size());
    pio(0) = 0;

    // Make pio(length+1)>pio(length), and
    // make pio(length)=pin(length)
    if (pio(length)>=pin(length)) {
      pin(length) = pio(length);
    } else {
      pio(length) = pin(length);
    }
    pio(length+1) = pio(length)+1.0;

    // For each k, run binary search and check that the result maches the
    // post-condition.
    for (int k=0; k<length; ++k) {
      int kk = k;
      binary_search(pio,pin(k+1),kk);
      if (kk==(length+1)) {
        kk = length;
      }
      REQUIRE ( (pio(kk) <= pin(k+1) &&
                pio(kk+1) >= pin(k+1) &&
                !(pio(kk+1) == pin(k+1) && k!=length)) );
    }
  }
}
