#include <catch2/catch.hpp>

#include <random>

#include "Types.hpp"
#include "Context.hpp"
#include "VerticalRemapManager.hpp"
#include "SimulationParams.hpp"
#include "Elements.hpp"
#include "HybridVCoord.hpp"
#include "Tracers.hpp"
#include "mpi/Connectivity.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"
#include "utilities/MathUtils.hpp"

using namespace Homme;

extern "C" {
void init_remap_f90 (const int& ne,
                     const Real* hyai_ptr, const Real* hybi_ptr,
                     const Real* hyam_ptr, const Real* hybm_ptr,
                     Real* dvv, Real* mp,
                     const Real& ps0);
void init_phis_f90 (const Real*& phis_ptr, const Real*& gradphis_ptr);

void run_remap_f90 (const int& np1, const int& np1_qdp, const Real& dt,
                    const int& rsplit, const int& qsize, const int& remap_alg,
                    Real*& dp_ptr, Real*& vtheta_dp_ptr, Real*& w_i_ptr,
                    Real*& phi_i_ptr, Real*& v_ptr, Real*& ps_ptr);
void cleanup_f90();
} // extern "C"


TEST_CASE("remap", "remap_testing") {

  // Catch runs these blocks of code multiple times, namely once per each
  // session within each test case. This is problematic for Context, which
  // is a static singleton.
  // We cannot call 'create' unless we are sure the object is not already stored
  // in the context. One solution is to call 'create_if_not_there', but that's not what
  // happens in mpi_cxx_f90_interface, which is called by the geometry_interface
  // fortran module.
  // Two solutions:
  //  - cleaning up the context at the end of TEST_CASE: this would also delete
  //    the comm object in the context, so you have to re-create it.
  //    Notice, however, that the comm would *already be there* when this block
  //    of code is executed for the first time (is created in tester.cpp),
  //    so you need to check if it's there first.
  //  - change mpi_cxx_f90_interface, to create the Connectivity only if not
  //    already present.
  //
  // Among the two, the former seems cleaner, since it does not affect the
  // src folder of Homme, only the test one. So I'm going with that.
  // More precisely, I'm getting a copy of the existing Comm from the context,
  // and reset it back in it after the cleanup

  constexpr int ne = 2;

  // The random numbers generator
  std::random_device rd;
  using rngAlg = std::mt19937_64;
  const int seed = 1; // Change to the following line after debugging
  // const int seed = rd();
  rngAlg engine(seed);
  using RPDF = std::uniform_real_distribution<Real>;
  using IPDF = std::uniform_int_distribution<int>;

  // Use stuff from Context, to increase similarity with actual runs
  auto& c = Context::singleton();

  // Init parameters
  auto& params = c.create<SimulationParams>();
  params.params_set = true;

  // Create and init hvcoord and ref_elem, needed to init the fortran interface
  auto& hvcoord = c.create<HybridVCoord>();
  hvcoord.random_init(seed);

  auto hyai = Kokkos::create_mirror_view(hvcoord.hybrid_ai);
  auto hybi = Kokkos::create_mirror_view(hvcoord.hybrid_bi);
  auto hyam = Kokkos::create_mirror_view(hvcoord.hybrid_am);
  auto hybm = Kokkos::create_mirror_view(hvcoord.hybrid_bm);
  Kokkos::deep_copy(hyai,hvcoord.hybrid_ai);
  Kokkos::deep_copy(hybi,hvcoord.hybrid_bi);
  Kokkos::deep_copy(hyam,hvcoord.hybrid_am);
  Kokkos::deep_copy(hybm,hvcoord.hybrid_bm);
  const Real* hyai_ptr  = hyai.data();
  const Real* hybi_ptr  = hybi.data();
  const Real* hyam_ptr  = reinterpret_cast<Real*>(hyam.data());
  const Real* hybm_ptr  = reinterpret_cast<Real*>(hybm.data());

  std::vector<Real> dvv(NP*NP);
  std::vector<Real> mp(NP*NP);

  // This will also init the c connectivity.
  init_remap_f90(ne,hyai_ptr,hybi_ptr,hyam_ptr,hybm_ptr,dvv.data(),mp.data(),hvcoord.ps0);
  const int num_elems = c.get<Connectivity>().get_num_local_elements();
  params.qsize = std::max(int(QSIZE_D-1),0);

  // Create and init elements/tracers
  auto& elems = c.create<Elements>();
  elems.init(num_elems,false,true);
  const auto max_pressure = 1000.0 + hvcoord.ps0; // This ensures max_p > ps0
  auto& geo = elems.m_geometry;
  elems.m_geometry.randomize(seed); // Only needed for phis and gradphis

  // Do not use all tracers, for better testing 
  auto& tracers = c.create<Tracers>();
  tracers.init(elems.num_elems(),params.qsize);
  tracers.randomize(seed);

  // Init f90
  auto d        = Kokkos::create_mirror_view(geo.m_d);
  auto dinv     = Kokkos::create_mirror_view(geo.m_dinv);
  auto phis     = Kokkos::create_mirror_view(geo.m_phis);
  auto gradphis = Kokkos::create_mirror_view(geo.m_gradphis);
  auto fcor     = Kokkos::create_mirror_view(geo.m_fcor);
  auto spmp     = Kokkos::create_mirror_view(geo.m_spheremp);
  auto rspmp    = Kokkos::create_mirror_view(geo.m_rspheremp);
  auto tVisc    = Kokkos::create_mirror_view(geo.m_tensorvisc);
  auto sph2c    = Kokkos::create_mirror_view(geo.m_vec_sph2cart);
  auto mdet     = Kokkos::create_mirror_view(geo.m_metdet);
  auto minv     = Kokkos::create_mirror_view(geo.m_metinv);

  const Real* phis_ptr     = phis.data();
  const Real* gradphis_ptr = gradphis.data();

  // For remap, we only need phis and gradphis. Copy from cxx to f90.
  Kokkos::deep_copy(phis,geo.m_phis);
  Kokkos::deep_copy(gradphis,geo.m_gradphis);
  init_phis_f90(phis_ptr,gradphis_ptr);

  using Scalar2dF90  = HostViewManaged<Real*[NUM_TIME_LEVELS][NP][NP]>;
  using ScalarF90    = HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>;
  using ScalarIntF90 = HostViewManaged<Real*[NUM_INTERFACE_LEV][NP][NP]>;
  using VectorF90    = HostViewManaged<Real*[NUM_PHYSICAL_LEV][2][NP][NP]>;
  using ScalarStateF90    = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>;
  using ScalarStateIntF90 = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>;
  using VectorStateF90    = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]>;

  Scalar2dF90       ps_f90("",elems.num_elems());
  ScalarStateF90    dp3d_f90("",elems.num_elems());
  ScalarStateF90    vtheta_dp_f90("",elems.num_elems());
  ScalarStateIntF90 w_i_f90("",elems.num_elems());
  ScalarStateIntF90 phinh_i_f90("",elems.num_elems());
  VectorStateF90    v_f90("",elems.num_elems());

  ScalarF90    omega_p_f90("",elems.num_elems());
  ScalarIntF90 eta_dot_dpdn_f90("",elems.num_elems());

  // Lambda to compute min of dp3d, to give meaningful initial value to derived.m_eta_dot_dpdn
  auto dp3d_min = [] (decltype(elems.m_state.m_dp3d) dp3d) -> Real {
    Real the_min = std::numeric_limits<Real>::max();
    for (int ie = 0; ie < dp3d.extent_int(0); ++ie) {
      // Because this constraint is difficult to satisfy for all of the tensors,
      // incrementally generate the view
      for (int tl = 0; tl < NUM_TIME_LEVELS; ++tl) {
        for (int igp = 0; igp < NP; ++igp) {
          for (int jgp = 0; jgp < NP; ++jgp) {
            ExecViewUnmanaged<Scalar[NUM_LEV]> pt_dp3d =
                Homme::subview(dp3d, ie, tl, igp, jgp);
            auto h_dp3d = Kokkos::create_mirror_view(pt_dp3d);
            Kokkos::deep_copy(h_dp3d,pt_dp3d);
            for (int ilev=0; ilev<NUM_LEV; ++ilev) {
              for (int iv=0; iv<VECTOR_SIZE; ++iv) {
                the_min = std::min(the_min,h_dp3d(ilev)[iv]);
              }
            }
          }
        }
      }
    }
    return the_min;
  };

  SECTION ("run_remap") {
    auto remap_algs = {RemapAlg::PPM_FIXED_PARABOLA, RemapAlg::PPM_MIRRORED, RemapAlg::PPM_FIXED_MEANS};
    auto remap_alg_f90 = [](const RemapAlg alg)->int {
      if (alg==RemapAlg::PPM_MIRRORED) {
        return 1;
      } else if (alg==RemapAlg::PPM_FIXED_PARABOLA) {
        return 2;
      } else if (alg==RemapAlg::PPM_FIXED_MEANS) {
        return 3;
      }
      return -1;
    };
    for (const bool hydrostatic : {true, false}) {
      std::cout << " -> " << (hydrostatic ? "hydrostatic" : "non-hydrostatic") << "\n";
      for (const int rsplit : {3,0}) {
        std::cout << "   -> rsplit = " << rsplit << "\n";
        for (auto alg : remap_algs) {
          std::cout << "     -> remap alg = " << remapAlg2str(alg) << "\n";
          // Set the parameters
          params.rsplit = rsplit;
          params.remap_alg = alg;
          params.theta_hydrostatic_mode = hydrostatic;

          // Generate timestep stage data
          const Real dt      = RPDF(1.0,100.0)(engine);
          const int  np1     = IPDF(0,NUM_TIME_LEVELS-1)(engine);
          const int  np1_qdp = IPDF(0,Q_NUM_TIME_LEVELS-1)(engine);

          // Randomize state/derived
          elems.m_state.randomize(seed,max_pressure,hvcoord.ps0,hvcoord.hybrid_ai0,geo.m_phis);
          elems.m_derived.randomize(seed,dp3d_min(elems.m_state.m_dp3d));

          // Note: to avoid errors in the equation of state, we need phi to be increasing.
          //       Moreover, f90 caar 'assumes' phi already stores phis, while cxx caar uses phis.

          // Copy initial values to f90
          sync_to_host(elems.m_state.m_dp3d, dp3d_f90);
          sync_to_host(elems.m_state.m_vtheta_dp, vtheta_dp_f90);
          sync_to_host(elems.m_state.m_w_i, w_i_f90);
          sync_to_host(elems.m_state.m_phinh_i, phinh_i_f90);
          sync_to_host(elems.m_state.m_v, v_f90);
          // The view for ps_v has the same layout in memory, so use Kokkos::deep_copy
          Kokkos::deep_copy(ps_f90,elems.m_state.m_ps_v);

          sync_to_host_p2i(elems.m_derived.m_eta_dot_dpdn, eta_dot_dpdn_f90);
          sync_to_host(elems.m_derived.m_omega_p, omega_p_f90);

          // Create the remap functor
          // Note: ALL the options must be set in params *before* creating the vrm.
          VerticalRemapManager vrm;
          vrm.run_remap(np1,np1_qdp,dt);

          auto dp3d_ptr = dp3d_f90.data();
          auto vtheta_dp_ptr = vtheta_dp_f90.data();
          auto w_i_ptr = w_i_f90.data();
          auto phinh_i_ptr = phinh_i_f90.data();
          auto v_ptr = v_f90.data();
          auto ps_ptr = ps_f90.data();
          run_remap_f90 (np1+1, np1_qdp+1, dt,
                         rsplit, params.qsize, remap_alg_f90(alg),
                         dp3d_ptr, vtheta_dp_ptr, w_i_ptr,
                         phinh_i_ptr, v_ptr, ps_ptr);

          // Compare answers
          auto h_dp3d      = Kokkos::create_mirror_view(elems.m_state.m_dp3d);
          auto h_vtheta_dp = Kokkos::create_mirror_view(elems.m_state.m_vtheta_dp);
          auto h_w_i       = Kokkos::create_mirror_view(elems.m_state.m_w_i);
          auto h_phinh_i   = Kokkos::create_mirror_view(elems.m_state.m_phinh_i);
          auto h_v         = Kokkos::create_mirror_view(elems.m_state.m_v);

          auto h_eta_dot_dpdn = Kokkos::create_mirror_view(elems.m_derived.m_eta_dot_dpdn);
          auto h_omega_p      = Kokkos::create_mirror_view(elems.m_derived.m_omega_p);

          Kokkos::deep_copy(h_dp3d     , elems.m_state.m_dp3d);
          Kokkos::deep_copy(h_vtheta_dp, elems.m_state.m_vtheta_dp);
          Kokkos::deep_copy(h_w_i      , elems.m_state.m_w_i);
          Kokkos::deep_copy(h_phinh_i  , elems.m_state.m_phinh_i);
          Kokkos::deep_copy(h_v        , elems.m_state.m_v);

          Kokkos::deep_copy(h_eta_dot_dpdn, elems.m_derived.m_eta_dot_dpdn);
          Kokkos::deep_copy(h_omega_p     , elems.m_derived.m_omega_p);

          for (int ie=0; ie<num_elems; ++ie) {
            auto dp3d_cxx      = viewAsReal(Homme::subview(h_dp3d,ie,np1));
            auto vtheta_dp_cxx = viewAsReal(Homme::subview(h_vtheta_dp,ie,np1));
            auto w_i_cxx       = viewAsReal(Homme::subview(h_w_i,ie,np1));
            auto phinh_i_cxx   = viewAsReal(Homme::subview(h_phinh_i,ie,np1));
            auto v_cxx         = viewAsReal(Homme::subview(h_v,ie,np1));

            auto eta_dot_dpdn_cxx = viewAsReal(Homme::subview(h_eta_dot_dpdn,ie));
            auto omega_p_cxx      = viewAsReal(Homme::subview(h_omega_p,ie));

            for (int igp=0; igp<NP; ++igp) {
              for (int jgp=0; jgp<NP; ++jgp) {
                for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
                  // dp3d
                  if(dp3d_cxx(igp,jgp,k)!=dp3d_f90(ie,np1,k,igp,jgp)) {
                    printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf("dp3d cxx: %3.40f\n",dp3d_cxx(igp,jgp,k));
                    printf("dp3d f90: %3.40f\n",dp3d_f90(ie,np1,k,igp,jgp));
                  }
                  REQUIRE(dp3d_cxx(igp,jgp,k)==dp3d_f90(ie,np1,k,igp,jgp));

                  // vtheta_dp
                  if(vtheta_dp_cxx(igp,jgp,k)!=vtheta_dp_f90(ie,np1,k,igp,jgp)) {
                    printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf("vtheta_dp cxx: %3.40f\n",vtheta_dp_cxx(igp,jgp,k));
                    printf("vtheta_dp f90: %3.40f\n",vtheta_dp_f90(ie,np1,k,igp,jgp));
                  }
                  REQUIRE(vtheta_dp_cxx(igp,jgp,k)==vtheta_dp_f90(ie,np1,k,igp,jgp));

                  // w_i
                  if(w_i_cxx(igp,jgp,k)!=w_i_f90(ie,np1,k,igp,jgp)) {
                    printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf("w_i cxx: %3.40f\n",w_i_cxx(igp,jgp,k));
                    printf("w_i f90: %3.40f\n",w_i_f90(ie,np1,k,igp,jgp));
                  }
                  REQUIRE(w_i_cxx(igp,jgp,k)==w_i_f90(ie,np1,k,igp,jgp));

                  // phinh_i
                  if(phinh_i_cxx(igp,jgp,k)!=phinh_i_f90(ie,np1,k,igp,jgp)) {
                    printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf("phinh_i cxx: %3.40f\n",phinh_i_cxx(igp,jgp,k));
                    printf("phinh_i f90: %3.40f\n",phinh_i_f90(ie,np1,k,igp,jgp));
                  }
                  REQUIRE(phinh_i_cxx(igp,jgp,k)==phinh_i_f90(ie,np1,k,igp,jgp));

                  // u
                  if(v_cxx(0,igp,jgp,k)!=v_f90(ie,np1,k,0,igp,jgp)) {
                    printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf("u cxx: %3.40f\n",v_cxx(0,igp,jgp,k));
                    printf("u f90: %3.40f\n",v_f90(ie,np1,k,0,igp,jgp));
                  }
                  REQUIRE(v_cxx(0,igp,jgp,k)==v_f90(ie,np1,k,0,igp,jgp));

                  // v
                  if(v_cxx(1,igp,jgp,k)!=v_f90(ie,np1,k,1,igp,jgp)) {
                    printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf("v cxx: %3.40f\n",v_cxx(1,igp,jgp,k));
                    printf("v f90: %3.40f\n",v_f90(ie,np1,k,1,igp,jgp));
                  }
                  REQUIRE(v_cxx(1,igp,jgp,k)==v_f90(ie,np1,k,1,igp,jgp));
                }

                // Check last interface for w_i and phinh_i
                int k = NUM_PHYSICAL_LEV;
                if(w_i_cxx(igp,jgp,k)!=w_i_f90(ie,np1,k,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("w_i cxx: %3.40f\n",w_i_cxx(igp,jgp,k));
                  printf("w_i f90: %3.40f\n",w_i_f90(ie,np1,k,igp,jgp));
                }
                REQUIRE(w_i_cxx(igp,jgp,k)==w_i_f90(ie,np1,k,igp,jgp));
                if(phinh_i_cxx(igp,jgp,k)!=phinh_i_f90(ie,np1,k,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("phinh_i cxx: %3.40f\n",phinh_i_cxx(igp,jgp,k));
                  printf("phinh_i f90: %3.40f\n",phinh_i_f90(ie,np1,k,igp,jgp));
                }
                REQUIRE(phinh_i_cxx(igp,jgp,k)==phinh_i_f90(ie,np1,k,igp,jgp));
            }}
          }
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
