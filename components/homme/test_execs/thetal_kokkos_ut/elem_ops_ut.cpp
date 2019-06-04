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
void set_theta_ref_f90(const int& num_elems,
                       const Real*& dp,
                       Real*& theta_ref);
} // extern "C"

// ============= ELEMENT OPS ================ //

TEST_CASE("elem_ops", "elem_ops") {

  constexpr int num_elems = 2;

  std::random_device rd;
  using rngAlg = std::mt19937_64;
  rngAlg engine(rd());
  std::uniform_real_distribution<Real> pdf(0.01, 1.0);

  HybridVCoord hvcoord;
  hvcoord.random_init(rd());

  decltype(hvcoord.hybrid_ai)::HostMirror hyai = Kokkos::create_mirror_view(hvcoord.hybrid_ai);
  Kokkos::deep_copy(hyai,hvcoord.hybrid_ai);
  const Real* hyai_ptr = hyai.data();

  // F90 and CXX views (and host mirrors of cxx views)
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  dp_f90("",num_elems);
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  theta_ref_f90("",num_elems);

  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   dp_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   theta_ref_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   p_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> p_i_cxx("",num_elems);

  auto h_theta_ref = Kokkos::create_mirror_view(theta_ref_cxx);

  ElementOps elem_ops;
  ColumnOps  col_ops;

  SECTION("theta_ref") {
    // Create random inputs
    genRandArray(theta_ref_cxx,engine,pdf);
    genRandArray(dp_cxx,engine,pdf);

    sync_to_host(theta_ref_cxx,theta_ref_f90);
    sync_to_host(dp_cxx,dp_f90);

    // Run fortran version
    init_f90(hyai_ptr,hvcoord.ps0);
    const Real* dp_f90_ptr = dp_f90.data();
    Real* theta_ref_f90_ptr = theta_ref_f90.data();
    set_theta_ref_f90(num_elems,dp_f90_ptr,theta_ref_f90_ptr);

    // Run cxx versino
    elem_ops.init(hvcoord);

    auto policy = Homme::get_default_team_policy<ExecSpace>(num_elems);
    Kokkos::parallel_for(policy,KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        auto dp  = Homme::subview(dp_cxx,kv.ie,igp,jgp);
        auto p   = Homme::subview(p_cxx,kv.ie,igp,jgp);
        auto p_i = Homme::subview(p_i_cxx,kv.ie,igp,jgp);
        p_i(0)[0] = hvcoord.hybrid_ai0*hvcoord.ps0;

        col_ops.column_scan_mid_to_int<true>(kv,dp,p_i);
        col_ops.compute_midpoint_values(kv,p_i,p);

        auto theta = Homme::subview(theta_ref_cxx,kv.ie,igp,jgp);
        elem_ops.compute_theta_ref(kv,p,theta);
      });
    });

    // Compare answers
    Kokkos::deep_copy(h_theta_ref,theta_ref_cxx);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<num_elems; ++igp) {
        for (int jgp=0; jgp<num_elems; ++jgp) {
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            const int ilev = k / VECTOR_SIZE;
            const int ivec = k % VECTOR_SIZE;
            REQUIRE (h_theta_ref(ie,igp,jgp,ilev)[ivec] == theta_ref_f90(ie,k,igp,jgp));
          }
        }
      }
    }
  }
}
