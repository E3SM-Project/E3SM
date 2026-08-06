#include <catch2/catch.hpp>

#include <random>

#include "ElementOps.hpp"
#include "Types.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

using namespace Homme;

extern "C" {
void init_f90 (const Real* hyai_ptr, const Real& ps0);
void compute_r_star_f90(const int& num_elems,
                        const bool& moist,
                        const Real*& Q,
                        Real*& R_star);
} // extern "C"

// ============= ELEMENT OPS ================ //

TEST_CASE("elem_ops", "elem_ops") {

  constexpr int num_elems = 2;
  auto policy = Homme::get_default_team_policy<ExecSpace>(num_elems);

  std::random_device rd;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  using rngAlg = std::mt19937_64;
  rngAlg engine(seed);
  std::uniform_real_distribution<Real> pdf(0.01, 1.0);

  HybridVCoord hvcoord;
  hvcoord.random_init(seed);

  decltype(hvcoord.hybrid_ai)::HostMirror hyai = Kokkos::create_mirror_view(hvcoord.hybrid_ai);
  Kokkos::deep_copy(hyai,hvcoord.hybrid_ai);
  const Real* hyai_ptr = hyai.data();
  init_f90(hyai_ptr,hvcoord.ps0);

  // F90 and CXX views (and host mirrors of cxx views)
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  dp_f90("",num_elems);
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  theta_ref_f90("",num_elems);

  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   dp_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   theta_ref_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   p_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> p_i_cxx("",num_elems);

  auto h_theta_ref = Kokkos::create_mirror_view(theta_ref_cxx);

  ElementOps elem_ops;
  elem_ops.init(hvcoord);

  SECTION("r_star") {
    for (bool moist : {false, true}) {

      // Use some existing buffers
      auto Q_cxx = theta_ref_cxx;
      auto Q_f90 = theta_ref_f90;
      auto R_f90 = dp_f90;
      genRandArray(Q_cxx,engine,pdf);
      sync_to_host(Q_cxx,Q_f90);

      // Run cxx version
      Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const TeamMember& team) {
        KernelVariables kv(team);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                             [&](const int idx) {
          const int igp = idx / NP;
          const int jgp = idx % NP;

          auto Q = Homme::subview(Q_cxx,kv.ie,igp,jgp);
          auto R = Homme::subview(dp_cxx,kv.ie,igp,jgp);

          elem_ops.get_R_star(kv,moist,Q,R);
        });
      });

      // Run f90 version
      const Real* Q_ptr = Q_f90.data();
      Real* R_ptr = R_f90.data();
      compute_r_star_f90(num_elems,moist,Q_ptr,R_ptr);

      // Compare answers
      auto R_cxx = Kokkos::create_mirror_view(dp_cxx);
      Kokkos::deep_copy(R_cxx,dp_cxx);

      for (int ie=0; ie<num_elems; ++ie) {
        for (int igp=0; igp<NP; ++igp) {
          for (int jgp=0; jgp<NP; ++jgp) {
            for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
              const int ilev = k / VECTOR_SIZE;
              const int ivec = k % VECTOR_SIZE;
              REQUIRE (R_cxx(ie,igp,jgp,ilev)[ivec] == R_f90(ie,k,igp,jgp));
            }
          }
        }
      }
    }
  }

}
