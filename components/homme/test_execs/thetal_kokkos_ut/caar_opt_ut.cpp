#include <catch2/catch.hpp>

#include "Hommexx_config.h"

#ifdef HOMMEXX_ENABLE_CAAR_OPT

#include <random>
#include <cmath>

#include "Types.hpp"
#include "Context.hpp"
#include "CaarFunctorImpl-caar-opt.hpp"
#include "SimulationParams.hpp"
#include "Tracers.hpp"
#include "PhysicalConstants.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

using namespace Homme;

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
void run_caar_f90 (const int& nm1, const int& n0, const int& np1,
                   const Real& dt, const Real& eta_ave_w,
                   const Real& scale1, const Real& scale2, const Real& scale3,
                   const bool& hydrostatic, const bool& adv_conservative,
                   const int& rsplit, const int& pgrad,
                   Real*& dp_ptr, Real*& vtheta_dp_ptr, Real*& w_i_ptr,
                   Real*& phi_i_ptr, Real*& v_ptr,
                   Real*& vn0_ptr, Real*& etadot_dpdn_ptr, Real*& omega_p_ptr);

void run_limiter_f90 (const int& np1, Real*& dp_ptr, Real*& vtheta_dp_ptr);
void cleanup_f90();
} // extern "C"

TEST_CASE("caar_opt", "caar_opt_testing") {

  constexpr int ne = 2;
  constexpr Real tol = 1e-12;

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
  // Manually modify threshold values as in control mod
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
  init_caar_f90(ne,hyai_ptr,hybi_ptr,hyam_ptr,hybm_ptr,dvv.data(),mp.data(),hvcoord.ps0);
  ref_FE.init_mass(mp.data());
  ref_FE.init_deriv(dvv.data());

  // Create and init elements
  const int num_elems = c.get<Connectivity>().get_num_local_elements();

  auto& elems = c.create<Elements>();
  elems.init(num_elems,true,PhysicalConstants::rearth0);
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
  init_geo_views_f90(d_ptr,dinv_ptr,phis_ptr,gradphis_ptr,fcor_ptr,
                     spmp_ptr,rspmp_ptr,tVisc_ptr,
                     sph2c_ptr,mdet_ptr,minv_ptr);

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

  using ScalarF90         = HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>;
  using ScalarIntF90      = HostViewManaged<Real*[NUM_INTERFACE_LEV][NP][NP]>;
  using VectorF90         = HostViewManaged<Real*[NUM_PHYSICAL_LEV][2][NP][NP]>;
  using ScalarStateF90    = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>;
  using ScalarStateIntF90 = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>;
  using VectorStateF90    = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]>;

  ScalarStateF90    dp3d_f90("",elems.num_elems());
  ScalarStateF90    vtheta_dp_f90("",elems.num_elems());
  ScalarStateIntF90 w_i_f90("",elems.num_elems());
  ScalarStateIntF90 phinh_i_f90("",elems.num_elems());
  VectorStateF90    v_f90("",elems.num_elems());

  ScalarF90    omega_p_f90("",elems.num_elems());
  ScalarIntF90 eta_dot_dpdn_f90("",elems.num_elems());
  VectorF90    vn0_f90("",elems.num_elems());

  // Lambda to compute min of dp3d, to give meaningful initial value to derived.m_eta_dot_dpdn
  auto dp3d_min = [] (decltype(elems.m_state.m_dp3d) dp3d) -> Real {
    Real the_min = std::numeric_limits<Real>::max();
    for (int ie = 0; ie < dp3d.extent_int(0); ++ie) {
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

  auto check_val = [&](const char* name, Real cxx_val, Real f90_val, int ie, int k, int igp, int jgp) {
    REQUIRE(!std::isnan(cxx_val));
    const Real rel_diff = compare_answers(f90_val, cxx_val);
    if (rel_diff >= tol) {
      printf("DIFF [%s] rank,ie,k,igp,jgp: %d, %d, %d, %d, %d\n", name, rank, ie, k, igp, jgp);
      printf("  cxx: %3.40f\n", cxx_val);
      printf("  f90: %3.40f\n", f90_val);
      printf("  rel_diff: %e (tol: %e)\n", rel_diff, tol);
    }
    REQUIRE(rel_diff < tol);
  };

  SECTION ("caar_opt_run") {
    for (const bool hydrostatic : {true,false}) {
      if (comm.root()) {
        std::cout << " -> " << (hydrostatic ? "Hydrostatic\n" : "Non-Hydrostatic\n");
      }
      auto adv_forms = {AdvectionForm::Conservative, AdvectionForm::NonConservative};
      for (const AdvectionForm adv_form : adv_forms) {
        if (comm.root()) {
          std::cout << "  -> " << (adv_form==AdvectionForm::Conservative ? "Conservative" : "Non-Conservative") << " theta advection\n";
        }
        for (int rsplit : {3,0}) {
          if (comm.root()) {
            std::cout << "   -> rsplit = " << rsplit << "\n";
          }
          for (const int pgrad : {1,0}) {
            if (comm.root()) {
              std::cout << "    -> pgrad_correction = " << pgrad << "\n";
            }
            // Set the parameters
            params.theta_hydrostatic_mode = hydrostatic;
            params.theta_adv_form = adv_form;
            params.rsplit = rsplit;
            params.pgrad_correction = (pgrad != 0);

            // Generate RK stage data
            Real dt = RPDF(1.0,10.0)(engine);
            Real eta_ave_w = RPDF(0.1,1.0)(engine);
            Real scale1 = RPDF(1.0,2.0)(engine);
            Real scale2 = RPDF(1.0,2.0)(engine);
            Real scale3 = RPDF(1.0,2.0)(engine);

            int  np1 = IPDF(0,2)(engine);

            // Sync scalars across ranks
            auto mpi_comm = Context::singleton().get<Comm>().mpi_comm();
            MPI_Bcast(&dt,1,MPI_DOUBLE,0,mpi_comm);
            MPI_Bcast(&scale1,1,MPI_DOUBLE,0,mpi_comm);
            MPI_Bcast(&scale2,1,MPI_DOUBLE,0,mpi_comm);
            MPI_Bcast(&scale3,1,MPI_DOUBLE,0,mpi_comm);
            MPI_Bcast(&eta_ave_w,1,MPI_DOUBLE,0,mpi_comm);
            MPI_Bcast(&np1,1,MPI_INT,0,mpi_comm);

            const int  n0  = (np1+1)%3;
            const int  nm1 = (np1+2)%3;

            RKStageData data (nm1, n0, np1, 0, dt, eta_ave_w, scale1, scale2, scale3);

            // Randomize state/derived
            elems.m_state.randomize(seed,max_pressure,hvcoord.ps0,hvcoord.hybrid_ai0,geo.m_phis);
            elems.m_derived.randomize(seed,dp3d_min(elems.m_state.m_dp3d));

            // Copy initial values to f90
            sync_to_host(elems.m_state.m_dp3d, dp3d_f90);
            sync_to_host(elems.m_state.m_vtheta_dp, vtheta_dp_f90);
            sync_to_host(elems.m_state.m_w_i, w_i_f90);
            sync_to_host(elems.m_state.m_phinh_i, phinh_i_f90);
            sync_to_host(elems.m_state.m_v, v_f90);

            sync_to_host<2>(elems.m_derived.m_vn0, vn0_f90);
            sync_to_host(elems.m_derived.m_eta_dot_dpdn, eta_dot_dpdn_f90);
            sync_to_host(elems.m_derived.m_omega_p, omega_p_f90);

            // Create the Caar functor with CAAR-opt
            CaarFunctorImpl caar(elems,tracers,ref_FE,hvcoord,sphop,params);
            FunctorsBuffersManager fbm;
            fbm.request_size( caar.requested_buffer_size() );
            fbm.request_size( limiter.requested_buffer_size() );
            fbm.allocate();
            caar.init_buffers(fbm);
            limiter.init_buffers(fbm);
            caar.init_boundary_exchanges(c.get_ptr<MpiBuffersManager>());

            // Run cxx
            caar.run(data);

            auto dp3d_ptr = dp3d_f90.data();
            auto vtheta_dp_ptr = vtheta_dp_f90.data();
            auto w_i_ptr = w_i_f90.data();
            auto phinh_i_ptr = phinh_i_f90.data();
            auto v_ptr = v_f90.data();
            auto vn0_ptr = vn0_f90.data();
            auto eta_dot_dpdn_ptr = eta_dot_dpdn_f90.data();
            auto omega_p_ptr = omega_p_f90.data();
            run_caar_f90 (nm1+1, n0+1, np1+1, dt, eta_ave_w, scale1, scale2, scale3,
                          hydrostatic, adv_form==AdvectionForm::Conservative, rsplit,
                          pgrad, dp3d_ptr, vtheta_dp_ptr, w_i_ptr,
                          phinh_i_ptr, v_ptr,
                          vn0_ptr, eta_dot_dpdn_ptr, omega_p_ptr);

            // Compare answers
            auto h_dp3d      = Kokkos::create_mirror_view(elems.m_state.m_dp3d);
            auto h_vtheta_dp = Kokkos::create_mirror_view(elems.m_state.m_vtheta_dp);
            auto h_w_i       = Kokkos::create_mirror_view(elems.m_state.m_w_i);
            auto h_phinh_i   = Kokkos::create_mirror_view(elems.m_state.m_phinh_i);
            auto h_v         = Kokkos::create_mirror_view(elems.m_state.m_v);

            auto h_vn0          = Kokkos::create_mirror_view(elems.m_derived.m_vn0);
            auto h_eta_dot_dpdn = Kokkos::create_mirror_view(elems.m_derived.m_eta_dot_dpdn);
            auto h_omega_p      = Kokkos::create_mirror_view(elems.m_derived.m_omega_p);

            Kokkos::deep_copy(h_dp3d     , elems.m_state.m_dp3d);
            Kokkos::deep_copy(h_vtheta_dp, elems.m_state.m_vtheta_dp);
            Kokkos::deep_copy(h_w_i      , elems.m_state.m_w_i);
            Kokkos::deep_copy(h_phinh_i  , elems.m_state.m_phinh_i);
            Kokkos::deep_copy(h_v        , elems.m_state.m_v);

            Kokkos::deep_copy(h_vn0         , elems.m_derived.m_vn0);
            Kokkos::deep_copy(h_eta_dot_dpdn, elems.m_derived.m_eta_dot_dpdn);
            Kokkos::deep_copy(h_omega_p     , elems.m_derived.m_omega_p);

            for (int ie=0; ie<num_elems; ++ie) {
              auto dp3d_cxx      = viewAsReal(Homme::subview(h_dp3d,ie,np1));
              auto vtheta_dp_cxx = viewAsReal(Homme::subview(h_vtheta_dp,ie,np1));
              auto w_i_cxx       = viewAsReal(Homme::subview(h_w_i,ie,np1));
              auto phinh_i_cxx   = viewAsReal(Homme::subview(h_phinh_i,ie,np1));
              auto v_cxx         = viewAsReal(Homme::subview(h_v,ie,np1));

              auto vn0_cxx          = viewAsReal(Homme::subview(h_vn0,ie));
              auto eta_dot_dpdn_cxx = viewAsReal(Homme::subview(h_eta_dot_dpdn,ie));
              auto omega_p_cxx      = viewAsReal(Homme::subview(h_omega_p,ie));

              for (int igp=0; igp<NP; ++igp) {
                for (int jgp=0; jgp<NP; ++jgp) {
                  for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
                    check_val("dp3d", dp3d_cxx(igp,jgp,k), dp3d_f90(ie,np1,k,igp,jgp), ie, k, igp, jgp);
                    check_val("vtheta_dp", vtheta_dp_cxx(igp,jgp,k), vtheta_dp_f90(ie,np1,k,igp,jgp), ie, k, igp, jgp);
                    check_val("w_i", w_i_cxx(igp,jgp,k), w_i_f90(ie,np1,k,igp,jgp), ie, k, igp, jgp);
                    check_val("phinh_i", phinh_i_cxx(igp,jgp,k), phinh_i_f90(ie,np1,k,igp,jgp), ie, k, igp, jgp);
                    check_val("u", v_cxx(0,igp,jgp,k), v_f90(ie,np1,k,0,igp,jgp), ie, k, igp, jgp);
                    check_val("v", v_cxx(1,igp,jgp,k), v_f90(ie,np1,k,1,igp,jgp), ie, k, igp, jgp);
                    check_val("un0", vn0_cxx(0,igp,jgp,k), vn0_f90(ie,k,0,igp,jgp), ie, k, igp, jgp);
                    check_val("vn0", vn0_cxx(1,igp,jgp,k), vn0_f90(ie,k,1,igp,jgp), ie, k, igp, jgp);
                    check_val("eta_dot_dpdn", eta_dot_dpdn_cxx(igp,jgp,k), eta_dot_dpdn_f90(ie,k,igp,jgp), ie, k, igp, jgp);
                    check_val("omega_p", omega_p_cxx(igp,jgp,k), omega_p_f90(ie,k,igp,jgp), ie, k, igp, jgp);
                  }

                  // Check last interface for w_i and phinh_i
                  const int k = NUM_PHYSICAL_LEV;
                  check_val("w_i_surf", w_i_cxx(igp,jgp,k), w_i_f90(ie,np1,k,igp,jgp), ie, k, igp, jgp);
                  check_val("phinh_i_surf", phinh_i_cxx(igp,jgp,k), phinh_i_f90(ie,np1,k,igp,jgp), ie, k, igp, jgp);
                }
              }
            }
          }
        }
      }
    }
  }

  SECTION ("caar_opt_split_form") {
    // Tests split-form advection (theta_advect_form = 2) in CAAR-opt path
    for (const bool hydrostatic : {true,false}) {
      for (int rsplit : {3,0}) {
        for (const int pgrad : {1,0}) {
          params.theta_hydrostatic_mode = hydrostatic;
          params.theta_adv_form = AdvectionForm::Split;
          params.rsplit = rsplit;
          params.pgrad_correction = (pgrad != 0);

          Real dt = RPDF(1.0,10.0)(engine);
          Real eta_ave_w = RPDF(0.1,1.0)(engine);
          Real scale1 = RPDF(1.0,2.0)(engine);
          Real scale2 = RPDF(1.0,2.0)(engine);
          Real scale3 = RPDF(1.0,2.0)(engine);
          int  np1 = IPDF(0,2)(engine);

          auto mpi_comm = Context::singleton().get<Comm>().mpi_comm();
          MPI_Bcast(&dt,1,MPI_DOUBLE,0,mpi_comm);
          MPI_Bcast(&scale1,1,MPI_DOUBLE,0,mpi_comm);
          MPI_Bcast(&scale2,1,MPI_DOUBLE,0,mpi_comm);
          MPI_Bcast(&scale3,1,MPI_DOUBLE,0,mpi_comm);
          MPI_Bcast(&eta_ave_w,1,MPI_DOUBLE,0,mpi_comm);
          MPI_Bcast(&np1,1,MPI_INT,0,mpi_comm);

          const int  n0  = (np1+1)%3;
          const int  nm1 = (np1+2)%3;

          RKStageData data (nm1, n0, np1, 0, dt, eta_ave_w, scale1, scale2, scale3);

          elems.m_state.randomize(seed,max_pressure,hvcoord.ps0,hvcoord.hybrid_ai0,geo.m_phis);
          elems.m_derived.randomize(seed,dp3d_min(elems.m_state.m_dp3d));

          CaarFunctorImpl caar(elems,tracers,ref_FE,hvcoord,sphop,params);
          FunctorsBuffersManager fbm;
          fbm.request_size( caar.requested_buffer_size() );
          fbm.request_size( limiter.requested_buffer_size() );
          fbm.allocate();
          caar.init_buffers(fbm);
          limiter.init_buffers(fbm);
          caar.init_boundary_exchanges(c.get_ptr<MpiBuffersManager>());

          caar.run(data);

          auto h_dp3d      = Kokkos::create_mirror_view(elems.m_state.m_dp3d);
          auto h_vtheta_dp = Kokkos::create_mirror_view(elems.m_state.m_vtheta_dp);
          auto h_w_i       = Kokkos::create_mirror_view(elems.m_state.m_w_i);
          auto h_phinh_i   = Kokkos::create_mirror_view(elems.m_state.m_phinh_i);
          auto h_v         = Kokkos::create_mirror_view(elems.m_state.m_v);

          Kokkos::deep_copy(h_dp3d     , elems.m_state.m_dp3d);
          Kokkos::deep_copy(h_vtheta_dp, elems.m_state.m_vtheta_dp);
          Kokkos::deep_copy(h_w_i      , elems.m_state.m_w_i);
          Kokkos::deep_copy(h_phinh_i  , elems.m_state.m_phinh_i);
          Kokkos::deep_copy(h_v        , elems.m_state.m_v);

          for (int ie=0; ie<num_elems; ++ie) {
            auto dp3d_cxx      = viewAsReal(Homme::subview(h_dp3d,ie,np1));
            auto vtheta_dp_cxx = viewAsReal(Homme::subview(h_vtheta_dp,ie,np1));
            auto w_i_cxx       = viewAsReal(Homme::subview(h_w_i,ie,np1));
            auto phinh_i_cxx   = viewAsReal(Homme::subview(h_phinh_i,ie,np1));
            auto v_cxx         = viewAsReal(Homme::subview(h_v,ie,np1));

            for (int igp=0; igp<NP; ++igp) {
              for (int jgp=0; jgp<NP; ++jgp) {
                for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
                  REQUIRE(!std::isnan(dp3d_cxx(igp,jgp,k)));
                  REQUIRE(!std::isnan(vtheta_dp_cxx(igp,jgp,k)));
                  REQUIRE(!std::isnan(v_cxx(0,igp,jgp,k)));
                  REQUIRE(!std::isnan(v_cxx(1,igp,jgp,k)));
                  REQUIRE(!std::isnan(w_i_cxx(igp,jgp,k)));
                  REQUIRE(!std::isnan(phinh_i_cxx(igp,jgp,k)));
                  REQUIRE(std::isfinite(dp3d_cxx(igp,jgp,k)));
                  REQUIRE(std::isfinite(vtheta_dp_cxx(igp,jgp,k)));
                  REQUIRE(std::isfinite(v_cxx(0,igp,jgp,k)));
                  REQUIRE(std::isfinite(v_cxx(1,igp,jgp,k)));
                }
              }
            }
          }
        }
      }
    }
  }

  SECTION ("limiter_dp3d") {

    params.rsplit = 1;
    params.theta_hydrostatic_mode = false;

    elems.m_state.randomize(seed,max_pressure,hvcoord.ps0,hvcoord.hybrid_ai0,geo.m_phis);

    sync_to_host(elems.m_state.m_dp3d, dp3d_f90);
    sync_to_host(elems.m_state.m_vtheta_dp, vtheta_dp_f90);

    auto dp3d_ptr = dp3d_f90.data();
    auto vtheta_dp_ptr = vtheta_dp_f90.data();

    CaarFunctorImpl caar(elems,tracers,ref_FE,hvcoord,sphop,params);
    FunctorsBuffersManager fbm;
    fbm.request_size( caar.requested_buffer_size() );
    fbm.request_size( limiter.requested_buffer_size() );
    fbm.allocate();
    caar.init_buffers(fbm);
    limiter.init_buffers(fbm);

    int  np1 = IPDF(0,2)(engine);
    RKStageData data;
    data.np1 = np1;
    caar.set_rk_stage_data(data);

    limiter.run(np1);

    run_limiter_f90(np1+1, dp3d_ptr, vtheta_dp_ptr);

    auto h_dp3d      = Kokkos::create_mirror_view(elems.m_state.m_dp3d);
    auto h_vtheta_dp = Kokkos::create_mirror_view(elems.m_state.m_vtheta_dp);

    Kokkos::deep_copy(h_dp3d     , elems.m_state.m_dp3d);
    Kokkos::deep_copy(h_vtheta_dp, elems.m_state.m_vtheta_dp);
    for (int ie=0; ie<num_elems; ++ie) {
      auto dp3d_cxx      = viewAsReal(Homme::subview(h_dp3d,ie,np1));
      auto vtheta_dp_cxx = viewAsReal(Homme::subview(h_vtheta_dp,ie,np1));

      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            check_val("limiter_dp3d", dp3d_cxx(igp,jgp,k), dp3d_f90(ie,np1,k,igp,jgp), ie, k, igp, jgp);
            check_val("limiter_vtheta_dp", vtheta_dp_cxx(igp,jgp,k), vtheta_dp_f90(ie,np1,k,igp,jgp), ie, k, igp, jgp);
          }
        }
      }
    }
  }

  // Cleanup
  auto old_comm = c.get_ptr<Comm>();
  c.finalize_singleton();
  auto& new_comm = c.create<Comm>();
  new_comm = *old_comm;

  cleanup_f90();
}

#endif // HOMMEXX_ENABLE_CAAR_OPT
