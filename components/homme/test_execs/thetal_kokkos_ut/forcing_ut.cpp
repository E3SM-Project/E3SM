#include <catch2/catch.hpp>

#include <iostream>
#include <random>

#include "Context.hpp"
#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "ElementsForcing.hpp"
#include "HybridVCoord.hpp"
#include "ForcingFunctor.hpp"
#include "PhysicalConstants.hpp"
#include "SimulationParams.hpp"
#include "Tracers.hpp"
#include "Types.hpp"

#include "mpi/Comm.hpp"

#include "utilities/MathUtils.hpp"
#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"

using namespace Homme;

template<typename T>
using HVM = HostViewManaged<T>;


// ============= THETA MODEL FORCING ================ //

extern "C" {
void init_f90 (const int& num_elems,
               const Real* hyai_ptr, const Real* hybi_ptr,
               const Real* hyam_ptr, const Real* hybm_ptr,
               const Real* gradphis,
               const Real& ps0, const int& qsize);
void set_forcing_pointers_f90 (Real*& q_ptr, Real*& fq_ptr, Real*& qdp_ptr,
                               Real*& v_ptr, Real*& w_ptr, Real*& vtheta_ptr, 
                               Real*& dp_ptr, Real*& phinh_ptr, Real*& ps_ptr,
                               Real*& fm_ptr, Real*& ft_ptr,
                               Real*& fvtheta_ptr, Real*& fphi_ptr);
void tracers_forcing_f90 (const Real& dt, const int& np1, const int& np1_qdp, const bool& hydrostatic, const bool& moist, const bool& adjustment);
void dynamics_forcing_f90 (const Real& dt, const int& np1);
void cleanup_f90();
} // extern "C"

TEST_CASE("forcing", "forcing") {
  using rngAlg = std::mt19937_64;
  using ipdf = std::uniform_int_distribution<int>;
  using dpdf = std::uniform_real_distribution<double>;

  std::random_device rd;
  constexpr int num_elems = 10;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  rngAlg engine(seed);

  // Init everything through singleton, which is what happens in normal runs
  auto& c = Context::singleton();
  auto& p = c.create<SimulationParams>();

  auto& hv = c.create<HybridVCoord>();
  hv.random_init(seed);

  auto& geo     = c.create<ElementsGeometry>();
  geo.init(num_elems,true, /* alloc_gradphis = */ true,
           PhysicalConstants::rearth0);
  geo.randomize(seed);

  auto& state   = c.create<ElementsState>();
  state.init(num_elems);

  auto& forcing = c.create<ElementsForcing>();
  forcing.init(num_elems);

  auto& tracers = c.create<Tracers>();
  p.qsize = ipdf(1,QSIZE_D)(engine);
  tracers.init(num_elems,p.qsize);
  // Kokkos::deep_copy(tracers.fq,0.0);

  const Real dt = dpdf(0.1,10.0)(engine);
  const int np1 = ipdf(0,2)(engine);
  const int np1_qdp = ipdf(0,1)(engine);

  // Init the f90 side
  auto h_hyai = Kokkos::create_mirror_view(hv.hybrid_ai);
  auto h_hybi = Kokkos::create_mirror_view(hv.hybrid_bi);
  auto h_hyam = Kokkos::create_mirror_view(hv.hybrid_am);
  auto h_hybm = Kokkos::create_mirror_view(hv.hybrid_bm);
  auto h_gradphis = Kokkos::create_mirror_view(geo.m_gradphis);
  Kokkos::deep_copy(h_hyai,hv.hybrid_ai);
  Kokkos::deep_copy(h_hybi,hv.hybrid_bi);
  Kokkos::deep_copy(h_hyam,hv.hybrid_am);
  Kokkos::deep_copy(h_hybm,hv.hybrid_bm);
  Kokkos::deep_copy(h_gradphis,geo.m_gradphis);
  init_f90(num_elems,
           h_hyai.data(), h_hybi.data(),
           reinterpret_cast<Real*>(h_hyam.data()),
           reinterpret_cast<Real*>(h_hybm.data()),
           h_gradphis.data(), hv.ps0, p.qsize);

  // Create f90-layout views
  HVM<Real*[QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]>                    q_f90 ("",num_elems);
  HVM<Real*[QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]>                    fq_f90 ("",num_elems);
  HVM<Real*[Q_NUM_TIME_LEVELS][QSIZE_D][NUM_PHYSICAL_LEV][NP][NP]> qdp_f90("",num_elems);

  HVM<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]> v_f90("",num_elems);
  HVM<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>   w_f90("",num_elems);
  HVM<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>    dp_f90("",num_elems);
  HVM<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>    vtheta_f90("",num_elems);
  HVM<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>   phinh_f90("",num_elems);
  HVM<Real*[NUM_TIME_LEVELS][NP][NP]>                      ps_f90("",num_elems);

  HVM<Real*[NUM_PHYSICAL_LEV][3][NP][NP]> fm_f90("",num_elems);
  HVM<Real*[NUM_PHYSICAL_LEV][NP][NP]>    ft_f90("",num_elems);
  HVM<Real*[NUM_PHYSICAL_LEV][NP][NP]>    fvtheta_f90("",num_elems);
  HVM<Real*[NUM_INTERFACE_LEV][NP][NP]>   fphi_f90("",num_elems);

  // Get pointers, and init f90
  Real* q_ptr      = q_f90.data();
  Real* fq_ptr     = fq_f90.data();
  Real* qdp_ptr    = qdp_f90.data();

  Real* v_ptr      = v_f90.data();
  Real* w_ptr      = w_f90.data();
  Real* vtheta_ptr = vtheta_f90.data();
  Real* dp_ptr     = dp_f90.data();
  Real* phinh_ptr  = phinh_f90.data();
  Real* ps_ptr     = ps_f90.data();

  Real* fm_ptr      = fm_f90.data();
  Real* ft_ptr      = ft_f90.data();
  Real* fvtheta_ptr = fvtheta_f90.data();
  Real* fphi_ptr    = fphi_f90.data();

  set_forcing_pointers_f90 (q_ptr, fq_ptr, qdp_ptr,
                            v_ptr, w_ptr, vtheta_ptr, dp_ptr, phinh_ptr, ps_ptr,
                            fm_ptr, ft_ptr, fvtheta_ptr, fphi_ptr);

  // Host mirrors of cxx views (for results checking)
  auto h_q       = Kokkos::create_mirror_view(tracers.Q);
  auto h_fq      = Kokkos::create_mirror_view(tracers.fq);
  auto h_qdp     = Kokkos::create_mirror_view(tracers.qdp);

  auto h_dp      = Kokkos::create_mirror_view(state.m_dp3d);
  auto h_v       = Kokkos::create_mirror_view(state.m_v);
  auto h_w       = Kokkos::create_mirror_view(state.m_w_i);
  auto h_vtheta  = Kokkos::create_mirror_view(state.m_vtheta_dp);
  auto h_phi     = Kokkos::create_mirror_view(state.m_phinh_i);

  auto h_fm      = Kokkos::create_mirror_view(forcing.m_fm);
  auto h_ft      = Kokkos::create_mirror_view(forcing.m_ft);
  auto h_fvtheta = Kokkos::create_mirror_view(forcing.m_fvtheta);
  auto h_fphi    = Kokkos::create_mirror_view(forcing.m_fphi);

  SECTION("tracers") {
    std::cout << "Testing tracers forcing.\n";
    for (const bool hydrostatic : {true,false}) {
      std::cout << " -> hydrostatic mode: " << (hydrostatic ? "true" : "false") << "\n";
      for (const MoistDry moisture : {MoistDry::DRY,MoistDry::MOIST}) {
        std::cout << "   -> moisture: " << (moisture==MoistDry::MOIST ? "moist" : "dry") << "\n";
        for (const bool adjustment : {true,false}) {
          std::cout << "     -> adjustment: " << (adjustment ? "true" : "false") << "\n";

          // Reset state, tracers, and forcing to the original random values
          state.randomize(seed, 10*hv.ps0, hv.ps0, hv.hybrid_ai0);
          tracers.randomize(seed,-1e-3,1e-3);
          forcing.randomize(seed);

          // Sync views
          sync_to_host(tracers.Q,q_f90);
          sync_to_host(tracers.fq,fq_f90);
          sync_to_host(tracers.qdp,qdp_f90);

          sync_to_host(state.m_v,v_f90);
          sync_to_host(state.m_w_i,w_f90);
          sync_to_host(state.m_vtheta_dp,vtheta_f90);
          sync_to_host(state.m_dp3d,dp_f90);
          sync_to_host(state.m_phinh_i,phinh_f90);
          // ps has same layout in cxx and f90
          Kokkos::deep_copy(ps_f90,state.m_ps_v);

          sync_to_host<3>(forcing.m_fm,fm_f90);
          sync_to_host(forcing.m_ft,ft_f90);
          sync_to_host(forcing.m_fvtheta,fvtheta_f90);
          sync_to_host(forcing.m_fphi,fphi_f90);

          p.theta_hydrostatic_mode = hydrostatic;

          // Create the functor
          FunctorsBuffersManager fbm;
          ForcingFunctor ff;

          // Create and set buffers in the forcing functor
          fbm.request_size(ff.requested_buffer_size());
          fbm.allocate();
          ff.init_buffers(fbm);

          // Run tracers forcing (cxx and f90)
          ff.tracers_forcing(dt,np1,np1_qdp,adjustment,moisture);
          tracers_forcing_f90(dt,np1+1,np1_qdp+1,hydrostatic,moisture==MoistDry::MOIST,adjustment);

          // Compare answers
          Kokkos::deep_copy(h_dp,state.m_dp3d);
          Kokkos::deep_copy(h_fphi,forcing.m_fphi);
          Kokkos::deep_copy(h_fvtheta,forcing.m_fvtheta);
          Kokkos::deep_copy(h_fq,tracers.fq);
          Kokkos::deep_copy(h_qdp,tracers.qdp);
          Kokkos::deep_copy(h_q,tracers.Q);
          for (int ie=0; ie<num_elems; ++ie) {
            for (int igp=0; igp<NP; ++igp) {
              for (int jgp=0; jgp<NP; ++jgp) {
                for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
                  const int ilev = k / VECTOR_SIZE;
                  const int ivec = k % VECTOR_SIZE;

                  if(h_dp(ie,np1,igp,jgp,ilev)[ivec]!=dp_f90(ie,np1,k,igp,jgp)) {
                    printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf ("dp_cxx: %3.16f\n",h_dp(ie,np1,igp,jgp,ilev)[ivec]);
                    printf ("dp_f90: %3.16f\n",dp_f90(ie,np1,k,igp,jgp));
                  }
                  REQUIRE(h_dp(ie,np1,igp,jgp,ilev)[ivec]==dp_f90(ie,np1,k,igp,jgp));

                  if(h_fphi(ie,igp,jgp,ilev)[ivec]!=fphi_f90(ie,k,igp,jgp)) {
                    printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf ("fphi_cxx: %3.16f\n",h_fphi(ie,igp,jgp,ilev)[ivec]);
                    printf ("fphi_f90: %3.16f\n",fphi_f90(ie,k,igp,jgp));
                  }
                  REQUIRE(h_fphi(ie,igp,jgp,ilev)[ivec]==fphi_f90(ie,k,igp,jgp));

                  if(h_fvtheta(ie,igp,jgp,ilev)[ivec]!=fvtheta_f90(ie,k,igp,jgp)) {
                    printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf ("fvtheta_cxx: %3.16f\n",h_fvtheta(ie,igp,jgp,ilev)[ivec]);
                    printf ("fvtheta_f90: %3.16f\n",fvtheta_f90(ie,k,igp,jgp));
                  }
                  REQUIRE(h_fvtheta(ie,igp,jgp,ilev)[ivec]==fvtheta_f90(ie,k,igp,jgp));

                  for (int iq=0; iq<p.qsize; ++iq) {
                    if(h_fq(ie,iq,igp,jgp,ilev)[ivec]!=fq_f90(ie,iq,k,igp,jgp)) {
                      printf ("ie,iq,k,igp,jgp: %d, %d, %d, %d, %d\n",ie,iq,k,igp,jgp);
                      printf ("fq_cxx: %3.16f\n",h_fq(ie,iq,igp,jgp,ilev)[ivec]);
                      printf ("fq_f90: %3.16f\n",fq_f90(ie,iq,k,igp,jgp));
                    }
                    REQUIRE(h_fq(ie,iq,igp,jgp,ilev)[ivec]==fq_f90(ie,iq,k,igp,jgp));

                    if(h_qdp(ie,np1_qdp,iq,igp,jgp,ilev)[ivec]!=qdp_f90(ie,np1_qdp,iq,k,igp,jgp)) {
                      printf ("ie,iq,k,igp,jgp: %d, %d, %d, %d, %d\n",ie,iq,k,igp,jgp);
                      printf ("qdp_cxx: %3.16f\n",h_qdp(ie,np1_qdp,iq,igp,jgp,ilev)[ivec]);
                      printf ("qdp_f90: %3.16f\n",qdp_f90(ie,np1_qdp,iq,k,igp,jgp));
                    }
                    REQUIRE(h_qdp(ie,np1_qdp,iq,igp,jgp,ilev)[ivec]==qdp_f90(ie,np1_qdp,iq,k,igp,jgp));

                    if(h_q(ie,iq,igp,jgp,ilev)[ivec]!=q_f90(ie,iq,k,igp,jgp)) {
                      printf ("ie,iq,k,igp,jgp: %d, %d, %d, %d, %d\n",ie,iq,k,igp,jgp);
                      printf ("q_cxx: %3.16f\n",h_q(ie,iq,igp,jgp,ilev)[ivec]);
                      printf ("q_f90: %3.16f\n",q_f90(ie,iq,k,igp,jgp));
                    }
                    REQUIRE(h_q(ie,iq,igp,jgp,ilev)[ivec]==q_f90(ie,iq,k,igp,jgp));
                  }
                }
                int k = NUM_INTERFACE_LEV-1;
                const int ilev = k / VECTOR_SIZE;
                const int ivec = k % VECTOR_SIZE;

                if(h_fphi(ie,igp,jgp,ilev)[ivec]!=fphi_f90(ie,k,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("fphi_cxx: %3.16f\n",h_fphi(ie,igp,jgp,ilev)[ivec]);
                  printf ("fphi_f90: %3.16f\n",fphi_f90(ie,k,igp,jgp));
                }
                REQUIRE(h_fphi(ie,igp,jgp,ilev)[ivec]==fphi_f90(ie,k,igp,jgp));
              }
            }
          }
        }
      }
    }
  }

  SECTION ("states") {

    // Set theta_hydrostatic_mode in simulation parameters, since it is
    // used later in the ForcingFunctor setup. If you don't set it, the FF
    // will have conditional jumps depending on uninited memory
    auto& p = Context::singleton().get<SimulationParams>();
    p.theta_hydrostatic_mode = false;

    // Reset state and forcing to the original random values
    std::cout << "Testing dynamics forcing.\n";

    state.randomize(seed, 10*hv.ps0, hv.ps0, hv.hybrid_ai0);
    forcing.randomize(seed);

    // Sync views
    sync_to_host(state.m_v,v_f90);
    sync_to_host(state.m_w_i,w_f90);
    sync_to_host(state.m_vtheta_dp,vtheta_f90);
    sync_to_host(state.m_dp3d,dp_f90);
    sync_to_host(state.m_phinh_i,phinh_f90);

    sync_to_host<3>(forcing.m_fm,fm_f90);
    sync_to_host(forcing.m_ft,ft_f90);
    sync_to_host(forcing.m_fvtheta,fvtheta_f90);
    sync_to_host(forcing.m_fphi,fphi_f90);

    // Create the functor
    FunctorsBuffersManager fbm;
    ForcingFunctor ff;

    // Create and set buffers in the forcing functor
    fbm.request_size(ff.requested_buffer_size());
    fbm.allocate();
    ff.init_buffers(fbm);

    // Run states forcing (cxx and f90)
    ff.states_forcing(dt,np1);
    dynamics_forcing_f90(dt,np1+1);

    // Sync results back to host
    Kokkos::deep_copy(h_v,state.m_v);
    Kokkos::deep_copy(h_w,state.m_w_i);
    Kokkos::deep_copy(h_phi,state.m_phinh_i);
    Kokkos::deep_copy(h_vtheta,state.m_vtheta_dp);

    // Compare answers
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            const int ilev = k / VECTOR_SIZE;
            const int ivec = k % VECTOR_SIZE;

            if(h_v(ie,np1,0,igp,jgp,ilev)[ivec]!=v_f90(ie,np1,k,0,igp,jgp)) {
              printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
              printf ("u_cxx: %3.16f\n",h_v(ie,np1,0,igp,jgp,ilev)[ivec]);
              printf ("u_f90: %3.16f\n",v_f90(ie,np1,0,k,igp,jgp));
            }
            REQUIRE(h_v(ie,np1,0,igp,jgp,ilev)[ivec]==v_f90(ie,np1,k,0,igp,jgp));

            if(h_v(ie,np1,1,igp,jgp,ilev)[ivec]!=v_f90(ie,np1,k,1,igp,jgp)) {
              printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
              printf ("v_cxx: %3.16f\n",h_v(ie,np1,1,igp,jgp,ilev)[ivec]);
              printf ("v_f90: %3.16f\n",v_f90(ie,np1,1,k,igp,jgp));
            }
            REQUIRE(h_v(ie,np1,1,igp,jgp,ilev)[ivec]==v_f90(ie,np1,k,1,igp,jgp));

            if(h_w(ie,np1,igp,jgp,ilev)[ivec]!=w_f90(ie,np1,k,igp,jgp)) {
              printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
              printf ("w_cxx: %3.16f\n",h_w(ie,np1,igp,jgp,ilev)[ivec]);
              printf ("w_f90: %3.16f\n",w_f90(ie,np1,k,igp,jgp));
            }
            REQUIRE(h_w(ie,np1,igp,jgp,ilev)[ivec]==w_f90(ie,np1,k,igp,jgp));

            if(h_phi(ie,np1,igp,jgp,ilev)[ivec]!=phinh_f90(ie,np1,k,igp,jgp)) {
              printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
              printf ("phi_cxx: %3.16f\n",h_phi(ie,np1,igp,jgp,ilev)[ivec]);
              printf ("phi_f90: %3.16f\n",phinh_f90(ie,np1,k,igp,jgp));
            }
            REQUIRE(h_phi(ie,np1,igp,jgp,ilev)[ivec]==phinh_f90(ie,np1,k,igp,jgp));

            if(h_vtheta(ie,np1,igp,jgp,ilev)[ivec]!=vtheta_f90(ie,np1,k,igp,jgp)) {
              printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
              printf ("u_cxx: %3.16f\n",h_vtheta(ie,np1,igp,jgp,ilev)[ivec]);
              printf ("u_f90: %3.16f\n",vtheta_f90(ie,np1,k,igp,jgp));
            }
            REQUIRE(h_vtheta(ie,np1,igp,jgp,ilev)[ivec]==vtheta_f90(ie,np1,k,igp,jgp));
          }

          const int k = NUM_PHYSICAL_LEV;
          const int ilev = k / VECTOR_SIZE;
          const int ivec = k % VECTOR_SIZE;
          REQUIRE(h_w(ie,np1,igp,jgp,ilev)[ivec]==w_f90(ie,np1,k,igp,jgp));
        }
      }
    }
  }

  // The tester.cpp file (where the 'main' is), inits the comm in
  // the context. When there are multiple test_cases/sections, we
  // need to make sure the context is returned in the same status
  // that it was found in. The comm is the only structure that
  // is present when we enter a test_case, so just copy it,
  // finalize the singleton, then re-set the (same) comm in the context.
  auto old_comm = c.get_ptr<Comm>();
  c.finalize_singleton();
  auto& new_comm = c.create<Comm>();
  new_comm = *old_comm;
  cleanup_f90();
}
