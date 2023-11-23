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

#define NNE 2
#define HOWMANY 3


using namespace Homme;


TEST_CASE("caar", "caar_testing") {

  constexpr int ne = NNE;

  // The random numbers generator
  std::random_device rd;
  using rngAlg = std::mt19937_64;
  const unsigned int catchRngSeed = Catch::rngSeed();
  //const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  const unsigned int seed = 1;
  //std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  std::cout << "seed: " << seed << "\n";
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
  //ref_FE.init_mass(mp.data());
  //ref_FE.init_deriv(dvv.data());

  ref_FE.random_init(seed);

  // Create and init elements
  const int num_elems = NNE ; //c.get<Connectivity>().get_num_local_elements();

  auto& elems = c.create<Elements>();
  elems.init(num_elems,false,true,PhysicalConstants::rearth0);
  const auto max_pressure = 1000.0 + hvcoord.ps0; // This ensures max_p > ps0
  auto& geo = elems.m_geometry;
  elems.m_geometry.randomize(seed); // Only needed for phis and gradphis

  // Get or create and init other structures needed by HVF
  auto& sphop = c.create<SphereOperators>();
  auto& tracers = c.create<Tracers>();
  auto& limiter = c.create<LimiterFunctor>(elems,hvcoord,params);

  sphop.setup(geo,ref_FE);

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


  SECTION ("caar_run") {

    //for (const bool hydrostatic : {true,false}) {
    for (const bool hydrostatic : {true}) {
        std::cout << " -> " << (hydrostatic ? "Hydrostatic\n" : "Non-Hydrostatic\n");
      //auto adv_forms = {AdvectionForm::Conservative, AdvectionForm::NonConservative};
      auto adv_forms = {AdvectionForm::NonConservative};
      for (const AdvectionForm adv_form : adv_forms) {
          std::cout << "  -> " << (adv_form==AdvectionForm::Conservative ? "Conservative" : "Non-Conservative") << " theta advection\n";
        //for (int rsplit : {3,0}) {
        for (int rsplit : {3}) {
            std::cout << "   -> rsplit = " << rsplit << "\n";
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

          const int  n0  = (np1+1)%3;
          const int  nm1 = (np1+2)%3;

          RKStageData data (nm1, n0, np1, 0, dt, eta_ave_w, scale1, scale2, scale3);

          // Randomize state/derived
          elems.m_state.randomize(seed,max_pressure,hvcoord.ps0,hvcoord.hybrid_ai0,geo.m_phis);
          elems.m_derived.randomize(seed,dp3d_min(elems.m_state.m_dp3d));

          // Create the Caar functor

          CaarFunctorImpl caar(elems,tracers,ref_FE,hvcoord,sphop,params);
          FunctorsBuffersManager fbm;
          fbm.request_size( caar.requested_buffer_size() );
          fbm.request_size(limiter.requested_buffer_size());
          fbm.allocate();
          caar.init_buffers(fbm);
          limiter.init_buffers(fbm);

          // Run cxx
	  for(int ii = 0; ii <  HOWMANY ; ii++)  caar.run(data);
	}
      }
    }
  }

}
