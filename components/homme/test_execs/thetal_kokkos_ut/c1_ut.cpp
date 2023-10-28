#include <catch2/catch.hpp>

#include <random>

#include "Types.hpp"
#include "Context.hpp"
#include "CaarFunctorImpl.hpp"
#include "SimulationParams.hpp"
#include "Tracers.hpp"
#include "PhysicalConstants.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

#define NNE 10
#define HOWMANY 20


using namespace Homme;

#if 0
extern "C" {
void init_caar_f90 (const int& ne,
               const Real* hyai_ptr, const Real* hybi_ptr,
               const Real* hyam_ptr, const Real* hybm_ptr,
               Real* dvv, Real* mp,
               const Real& ps0);
void init_geo_views_f90 (Real*& d_ptr,Real*& dinv_ptr,
               const Real*& phis_ptr, const Real*& gradphis_ptr,
               Real*& fcor_ptr,
               Real*& sphmp_ptr, Real*& rspmp_ptr,
               Real*& tVisc_ptr, Real*& sph2c_ptr,
               Real*& metdet_ptr, Real*& metinv_ptr);

void run_limiter_f90 (const int& np1, Real*& dp_ptr, Real*& vtheta_dp_ptr);
void cleanup_f90();
} // extern "C"
#endif


TEST_CASE("caar", "caar_testing") {

  constexpr int ne = NNE;

  // The random numbers generator
  std::random_device rd;
  using rngAlg = std::mt19937_64;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  rngAlg engine(seed);
  using RPDF = std::uniform_real_distribution<Real>;
  using IPDF = std::uniform_int_distribution<int>;

  // Use stuff from Context, to increase similarity with actual runs
  auto& c = Context::singleton();

  // Init parameters
  auto& params = c.create<SimulationParams>();
  params.params_set = true;
  //since init_params... is not called here and thus F setup is not transferred, manually 
  //modify thresholds values as in control mod
  params.dp3d_thresh = 0.125;
  params.vtheta_thresh = 100.0;

  // Create and init hvcoord and ref_elem, needed to init the fortran interface
  auto& hvcoord = c.create<HybridVCoord>();
  auto& ref_FE  = c.create<ReferenceElement>();
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
  //init_caar_f90(ne,hyai_ptr,hybi_ptr,hyam_ptr,hybm_ptr,dvv.data(),mp.data(),hvcoord.ps0);
  ref_FE.init_mass(mp.data());
  ref_FE.init_deriv(dvv.data());

  // Create and init elements
  const int num_elems = c.get<Connectivity>().get_num_local_elements();

  auto& elems = c.create<Elements>();
  elems.init(num_elems,false,true,PhysicalConstants::rearth0);
  const auto max_pressure = 1000.0 + hvcoord.ps0; // This ensures max_p > ps0
  auto& geo = elems.m_geometry;
  elems.m_geometry.randomize(seed); // Only needed for phis and gradphis

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
  Kokkos::deep_copy(phis,geo.m_phis);
  Kokkos::deep_copy(gradphis,geo.m_gradphis);

  Real* d_ptr        = d.data();
  Real* dinv_ptr     = dinv.data();
  Real* spmp_ptr     = spmp.data();
  Real* rspmp_ptr    = rspmp.data();
  Real* tVisc_ptr    = tVisc.data();
  Real* sph2c_ptr    = sph2c.data();
  Real* mdet_ptr     = mdet.data();
  Real* minv_ptr     = minv.data();
  const Real* phis_ptr     = phis.data();
  const Real* gradphis_ptr = gradphis.data();
  Real* fcor_ptr     = fcor.data();

  // Get the f90 values for geometric views.
//  init_geo_views_f90(d_ptr,dinv_ptr,phis_ptr,gradphis_ptr,fcor_ptr,
//                     spmp_ptr,rspmp_ptr,tVisc_ptr,
//                     sph2c_ptr,mdet_ptr,minv_ptr);

  Kokkos::deep_copy(geo.m_d,d);
  Kokkos::deep_copy(geo.m_dinv,dinv);
  Kokkos::deep_copy(geo.m_spheremp,spmp);
  Kokkos::deep_copy(geo.m_rspheremp,rspmp);
  Kokkos::deep_copy(geo.m_tensorvisc,tVisc);
  Kokkos::deep_copy(geo.m_vec_sph2cart,sph2c);
  Kokkos::deep_copy(geo.m_metdet,mdet);
  Kokkos::deep_copy(geo.m_metinv,minv);
  Kokkos::deep_copy(geo.m_fcor,fcor);

  // Get or create and init other structures needed by HVF
  auto& bm = c.create<MpiBuffersManager>();
  auto& sphop = c.create<SphereOperators>();
  auto& tracers = c.create<Tracers>();
  auto& limiter = c.create<LimiterFunctor>(elems,hvcoord,params);

  sphop.setup(geo,ref_FE);
  if (!bm.is_connectivity_set ()) {
    bm.set_connectivity(c.get_ptr<Connectivity>());
  }

//  using ScalarF90    = HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>;
//  using ScalarIntF90 = HostViewManaged<Real*[NUM_INTERFACE_LEV][NP][NP]>;
//  using VectorF90    = HostViewManaged<Real*[NUM_PHYSICAL_LEV][2][NP][NP]>;
//  using ScalarStateF90    = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>;
//  using ScalarStateIntF90 = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>;
//  using VectorStateF90    = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]>;

//  ScalarStateF90    dp3d_f90("",elems.num_elems());
//  ScalarStateF90    vtheta_dp_f90("",elems.num_elems());
//  ScalarStateIntF90 w_i_f90("",elems.num_elems());
//  ScalarStateIntF90 phinh_i_f90("",elems.num_elems());
//  VectorStateF90    v_f90("",elems.num_elems());

//  ScalarF90    omega_p_f90("",elems.num_elems());
//  ScalarIntF90 eta_dot_dpdn_f90("",elems.num_elems());
//  VectorF90    vn0_f90("",elems.num_elems());

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

  auto& comm = c.get<Comm>();
  const int rank = comm.rank();

  SECTION ("caar_run") {
    //for (const bool hydrostatic : {true,false}) {
    for (const bool hydrostatic : {true}) {
      if (comm.root()) {
        std::cout << " -> " << (hydrostatic ? "Hydrostatic\n" : "Non-Hydrostatic\n");
      }
      //auto adv_forms = {AdvectionForm::Conservative, AdvectionForm::NonConservative};
      auto adv_forms = {AdvectionForm::NonConservative};
      for (const AdvectionForm adv_form : adv_forms) {
        if (comm.root()) {
          std::cout << "  -> " << (adv_form==AdvectionForm::Conservative ? "Conservative" : "Non-Conservative") << " theta advection\n";
        }
        //for (int rsplit : {3,0}) {
        for (int rsplit : {3}) {
          if (comm.root()) {
            std::cout << "   -> rsplit = " << rsplit << "\n";
          }
          // Set the parameters
          params.theta_hydrostatic_mode = hydrostatic;
          params.theta_adv_form = adv_form;
          params.rsplit = rsplit;

          // Generate RK stage data
          Real dt = RPDF(1.0,10.0)(engine);
          Real eta_ave_w = RPDF(0.1,1.0)(engine);
          Real scale1 = RPDF(1.0,2.0)(engine);
          Real scale2 = RPDF(1.0,2.0)(engine);
          Real scale3 = RPDF(1.0,2.0)(engine);

          int  np1 = IPDF(0,2)(engine);

#if 0	  
          // Sync scalars across ranks (only np1 is *really* necessary, but might as well...)
          auto mpi_comm = Context::singleton().get<Comm>().mpi_comm();
          MPI_Bcast(&dt,1,MPI_DOUBLE,0,mpi_comm);
          MPI_Bcast(&scale1,1,MPI_DOUBLE,0,mpi_comm);
          MPI_Bcast(&scale2,1,MPI_DOUBLE,0,mpi_comm);
          MPI_Bcast(&scale3,1,MPI_DOUBLE,0,mpi_comm);
          MPI_Bcast(&eta_ave_w,1,MPI_DOUBLE,0,mpi_comm);
          MPI_Bcast(&np1,1,MPI_INT,0,mpi_comm);
#endif

          const int  n0  = (np1+1)%3;
          const int  nm1 = (np1+2)%3;

          RKStageData data (nm1, n0, np1, 0, dt, eta_ave_w, scale1, scale2, scale3);

          // Randomize state/derived
          elems.m_state.randomize(seed,max_pressure,hvcoord.ps0,hvcoord.hybrid_ai0,geo.m_phis);
          elems.m_derived.randomize(seed,dp3d_min(elems.m_state.m_dp3d));

#if 0	  
          // Copy initial values to f90
          sync_to_host(elems.m_state.m_dp3d, dp3d_f90);
          sync_to_host(elems.m_state.m_vtheta_dp, vtheta_dp_f90);
          sync_to_host(elems.m_state.m_w_i, w_i_f90);
          sync_to_host(elems.m_state.m_phinh_i, phinh_i_f90);
          sync_to_host(elems.m_state.m_v, v_f90);

          sync_to_host<2>(elems.m_derived.m_vn0, vn0_f90);
          sync_to_host(elems.m_derived.m_eta_dot_dpdn, eta_dot_dpdn_f90);
          sync_to_host(elems.m_derived.m_omega_p, omega_p_f90);
#endif

          // Create the Caar functor
          CaarFunctorImpl caar(elems,tracers,ref_FE,hvcoord,sphop,params);
          FunctorsBuffersManager fbm;
          fbm.request_size( caar.requested_buffer_size() );
          fbm.request_size(limiter.requested_buffer_size());
          fbm.allocate();
          caar.init_buffers(fbm);
          limiter.init_buffers(fbm);
          caar.init_boundary_exchanges(c.get_ptr<MpiBuffersManager>());

          // Run cxx
	  for(int ii = 0; ii <  HOWMANY ; ii++)  caar.run(data);
	}
      }
    }
  }

  // Cleanup (see comment at the top for explanation of the treatment of Comm)
//  auto old_comm = c.get_ptr<Comm>();
//  c.finalize_singleton();
//  auto& new_comm = c.create<Comm>();
//  new_comm = *old_comm;

//  cleanup_f90();
}
