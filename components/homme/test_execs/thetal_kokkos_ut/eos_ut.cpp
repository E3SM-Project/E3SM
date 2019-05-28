#include <catch2/catch.hpp>

#include <iostream>
#include <random>

#include "EquationOfState.hpp"
#include "Types.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"

using namespace Homme;

// ============= EQUATION OF STATE ================ //

extern "C" {
void init_f90 (const Real* hyai_ptr, const Real& ps0, const bool& hydrostatic);
void pnh_and_exner_from_eos_f90(const int& num_elems,
                                const Real*& vtheta_dp,
                                const Real*& dp3d,
                                const Real*& phi_i,
                                      Real*& pnh,
                                      Real*& exner,
                                      Real*& dpnh_dp_i);
} // extern "C"

TEST_CASE("pnh_and_exner_from_eos", "pnh_and_exner_from_eos") {

  constexpr auto LAST_LEV_P = ColInfo<NUM_INTERFACE_LEV>::LastPack;
  constexpr auto LAST_INTERFACE_VEC_IDX = ColInfo<NUM_INTERFACE_LEV>::LastVecEnd;

  constexpr int num_elems = 10;

  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  vtheta_dp_f90("",num_elems);
  HostViewManaged<Real*[NUM_INTERFACE_LEV][NP][NP]> phi_i_f90("",num_elems);
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  dp3d_f90("",num_elems);
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  pnh_f90("",num_elems);
  HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>  exner_f90("",num_elems);
  HostViewManaged<Real*[NUM_INTERFACE_LEV][NP][NP]> dpnh_dp_i_f90("",num_elems);

  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   vtheta_dp_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> phi_i_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   dp3d_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> dp3d_i_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   pi_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   pnh_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   exner_cxx("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> dpnh_dp_i_cxx("",num_elems);

  std::random_device rd;
  using rngAlg = std::mt19937_64;
  rngAlg engine(rd());
  std::uniform_real_distribution<Real> pdf(0.01, 1.0);

  genRandArray(vtheta_dp_cxx,engine,pdf);
  genRandArray(dp3d_cxx,engine,pdf);

  // Note: to avoid errors in pnh_and_exner_from_eos, we need phi to be increasing.
  //       Rather than using a constraint (which may call the function many times,
  //       we simply ask that there are no duplicates, then we sort it later.
  auto sort_and_chek = [](const ExecViewManaged<Scalar[NUM_LEV_P]>::HostMirror v)->bool {
    Real* start = reinterpret_cast<Real*>(v.data());
    Real* end   = reinterpret_cast<Real*>(v.data()) + NUM_LEV_P*VECTOR_SIZE;
    std::sort(start,end);
    std::reverse(start,end);
    auto it = std::unique(start,end);
    return it==end;
  };
  for (int ie=0; ie<num_elems; ++ie) {
    for (int igp=0; igp<NP; ++igp) {
      for (int jgp=0; jgp<NP; ++ jgp) {
        genRandArray(Homme::subview(phi_i_cxx,ie,igp,jgp),engine,pdf,sort_and_chek);
      }
    }
  }

  // Copy to f90 views

  auto h_vtheta_dp = Kokkos::create_mirror_view(vtheta_dp_cxx);
  auto h_phi_i = Kokkos::create_mirror_view(phi_i_cxx);
  auto h_dp3d = Kokkos::create_mirror_view(dp3d_cxx);
  for (int ie=0; ie<num_elems; ++ie) {
    for (int igp=0; igp<NP; ++igp) {
      for (int jgp=0; jgp<NP; ++ jgp) {
        for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
          const int ilev = k / VECTOR_SIZE;
          const int ivec = k % VECTOR_SIZE;

          vtheta_dp_f90(ie,k,igp,jgp) = h_vtheta_dp(ie,igp,jgp,ilev)[ivec];
          phi_i_f90(ie,k,igp,jgp) = h_phi_i(ie,igp,jgp,ilev)[ivec];
          dp3d_f90(ie,k,igp,jgp) = h_dp3d(ie,igp,jgp,ilev)[ivec];
        }
        phi_i_f90(ie,NUM_INTERFACE_LEV-1,igp,jgp) = h_phi_i(ie,igp,jgp,LAST_LEV_P)[LAST_INTERFACE_VEC_IDX];
      }
    }
  }

  bool hydrostatic = pdf(engine) > 0.5;

  HybridVCoord hvcoord;
  hvcoord.random_init(rd());

  // Init f90
  decltype(hvcoord.hybrid_ai)::HostMirror hyai = Kokkos::create_mirror_view(hvcoord.hybrid_ai);
  Kokkos::deep_copy(hyai,hvcoord.hybrid_ai);
  const Real* hyai_ptr = hyai.data();
  init_f90(hyai_ptr,hvcoord.ps0,hydrostatic);

  // Run f90
  const Real* vtheta_dp_ptr = vtheta_dp_f90.data();
  const Real* dp3d_ptr = dp3d_f90.data();
  const Real* phi_i_ptr = phi_i_f90.data();
  Real* pnh_ptr = pnh_f90.data();
  Real* exner_ptr = exner_f90.data();
  Real* dpnh_dp_i_ptr = dpnh_dp_i_f90.data();
  pnh_and_exner_from_eos_f90(num_elems,vtheta_dp_ptr,dp3d_ptr,phi_i_ptr,pnh_ptr,exner_ptr,dpnh_dp_i_ptr);

  ElementOps elem_ops;
  EquationOfState eos;
  eos.init(hydrostatic,hvcoord);

  // Compute pi from dp3d and compute dp3d_i
  Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                       KOKKOS_LAMBDA(const TeamMember& team) {
    KernelVariables kv(team);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      elem_ops.compute_interface_values(kv,Homme::subview(dp3d_cxx,kv.ie,igp,jgp),
                                           Homme::subview(dp3d_i_cxx,kv.ie,igp,jgp));

      Kokkos::single(Kokkos::PerThread(kv.team),[&](){
        auto pi = viewAsReal(Homme::subview(pi_cxx,kv.ie,igp,jgp));
        auto dp3d = viewAsReal(Homme::subview(dp3d_cxx,kv.ie,igp,jgp));

        pi(0) = hvcoord.hybrid_ai(0)*hvcoord.ps0;
        for (int k=0; k<NUM_PHYSICAL_LEV-1;++k) {
          pi(k+1) = pi(k) + dp3d(k);
          pi(k) += dp3d(k)/2.0;
        }
        pi(NUM_PHYSICAL_LEV-1) += dp3d(NUM_PHYSICAL_LEV-1)/2.0;
      });
    });
  });
  Kokkos::fence();

  // Compute pnh and exner
  Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                       KOKKOS_LAMBDA(const TeamMember& team) {
    KernelVariables kv(team);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto vtheta_dp = Homme::subview(vtheta_dp_cxx,kv.ie,igp,jgp);
      auto phi_i = Homme::subview(phi_i_cxx,kv.ie,igp,jgp);
      auto pi = Homme::subview(pi_cxx,kv.ie,igp,jgp);
      auto pnh = Homme::subview(pnh_cxx,kv.ie,igp,jgp);
      auto exner = Homme::subview(exner_cxx,kv.ie,igp,jgp);

      eos.compute_pnh_and_exner(kv,vtheta_dp,phi_i,pi,pnh,exner);
    });
    kv.team_barrier();
  });
  Kokkos::fence();

  // Compute dpnh_dp_i
  Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                       KOKKOS_LAMBDA(const TeamMember& team) {
    KernelVariables kv(team);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      auto dp_i = Homme::subview(dp3d_i_cxx,kv.ie,igp,jgp);
      auto pnh = Homme::subview(pnh_cxx,kv.ie,igp,jgp);
      auto dpnh_dp_i = Homme::subview(dpnh_dp_i_cxx,kv.ie,igp,jgp);

      eos.compute_dpnh_dp_i(kv,pnh,dp_i,dpnh_dp_i);
    });
    kv.team_barrier();
  });
  Kokkos::fence();

  // Now, compare results
  auto h_exner = Kokkos::create_mirror_view(exner_cxx);
  auto h_pnh = Kokkos::create_mirror_view(pnh_cxx);
  auto h_dpnh_dp_i = Kokkos::create_mirror_view(dpnh_dp_i_cxx);
  Kokkos::deep_copy(h_exner,exner_cxx);
  Kokkos::deep_copy(h_pnh,pnh_cxx);
  Kokkos::deep_copy(h_dpnh_dp_i,dpnh_dp_i_cxx);

  for (int ie=0; ie<num_elems; ++ie) {
    for (int igp=0; igp<NP; ++igp) {
      for (int jgp=0; jgp<NP; ++jgp) {
        for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
          const int ilev = k / VECTOR_SIZE;
          const int ivec = k % VECTOR_SIZE;

          // REQUIRE (h_exner(ie,igp,jgp,ilev)[ivec] == exner_f90(ie,k,igp,jgp));
          REQUIRE (h_pnh(ie,igp,jgp,ilev)[ivec] == pnh_f90(ie,k,igp,jgp));
          REQUIRE (h_dpnh_dp_i(ie,igp,jgp,ilev)[ivec] == dpnh_dp_i_f90(ie,k,igp,jgp));
        }
        REQUIRE (h_dpnh_dp_i(ie,igp,jgp,LAST_LEV_P)[LAST_INTERFACE_VEC_IDX] == dpnh_dp_i_f90(ie,NUM_INTERFACE_LEV-1,igp,jgp));
      }
    }
  }
}
