#include <catch2/catch.hpp>

#include <random>

#include "Types.hpp"
#include "Context.hpp"
#include "ElementsGeometry.hpp"
#include "ElementsState.hpp"
#include "ElementsDerivedState.hpp"
#include "FunctorsBuffersManager.hpp"
#include "HybridVCoord.hpp"
#include "HyperviscosityFunctorImpl.hpp"
#include "SimulationParams.hpp"
#include "SphereOperators.hpp"
#include "mpi/MpiBuffersManager.hpp"
#include "mpi/Connectivity.hpp"
#include "PhysicalConstants.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

using namespace Homme;

extern "C" {
void init_hv_f90 (const int& ne,
               const Real* hyai_ptr, const Real* hybi_ptr,
               const Real* hyam_ptr, const Real* hybm_ptr,
               Real* dvv, Real* mp,
               const Real& ps0, const int& hypervis_subcycle,
               const Real& nu, const Real& nu_div, const Real& nu_top,
               const Real& nu_p, const Real& nu_s);
void init_geo_views_f90 (Real*& d_ptr,Real*& dinv_ptr,
               const Real*& phis_ptr, const Real*& gradphis_ptr,
               Real*& fcor,
               Real*& sphmp_ptr, Real*& rspmp_ptr,
               Real*& tVisc_ptr, Real*& sph2c_ptr,
               Real*& metdet_ptr, Real*& metinv_ptr);
void initialize_reference_states_f90 (const Real*& phis,
                                      const Real*& dp_ref,
                                      const Real*& theta_ref,
                                      const Real*& phi_ref);
void biharmonic_wk_theta_f90 (const int& np1, const Real& hv_scaling, const bool& hydrostatic,
                              const Real*& dp, const Real*& vtheta_dp,
                              const Real*& w,  const Real*& phi, const Real*& v,
                              Real*& dptens, Real*& ttens, Real*& wtens,
                              Real*& phitens, Real*& vtens);
void advance_hypervis_f90 (const int& np1, const Real& dt, const Real& eta_ave_w,
                           const Real& hv_scaling, const bool& hydrostatic,
                           const Real*& dp_ref_ptr, const Real*& theta_ref_ptr, const Real*& phi_ref_ptr,
                           Real*& v_state, Real*& w_state, Real*& vtheta_state,
                           Real*& dp_state, Real*& phinh_state);
void cleanup_f90();
} // extern "C"

// This class is basically a HyperviscosityFunctorImpl, but
// using this instead of HVF we can access protected members
// of HVF (e.g., for initialization) without exposing them in HVF.

class HVFTester : public HyperviscosityFunctorImpl {
public:
  HVFTester (const SimulationParams&     params,
             const ElementsGeometry&     geometry,
             const ElementsState&        state,
             const ElementsDerivedState& derived)
   : HyperviscosityFunctorImpl(params,geometry,state,derived)
  {
    // Nothing to do here
  }

  ~HVFTester () = default;

  void set_timestep_data (const int np1, const Real dt, const Real eta_ave_w)
  {
    m_data.np1 = np1;
    m_data.dt = dt;
    m_data.dt_hvs = (m_data.hypervis_subcycle > 0 ) ? dt/m_data.hypervis_subcycle : -1.0;
    m_data.dt_hvs_tom = -1.0;// set to dt/m_data.hypervis_subcycle_tom;
    m_data.eta_ave_w = eta_ave_w;
  }

  void set_hv_data (const Real hv_scaling, const Real nu_ratio1, const Real nu_ratio2)
  {
    m_data.nu_ratio1 = nu_ratio1;
    m_data.nu_ratio2 = nu_ratio2;
    m_data.consthv = (hv_scaling==0.0);
  }

  using ScalarTens = ExecViewUnmanaged<Scalar*   [NP][NP][NUM_LEV]>;
#ifdef XX_NONBFB_COMING
  using ScalarTensInt = ExecViewUnmanaged<Scalar*   [NP][NP][NUM_LEV_P]>;
#else
  using ScalarTensInt = ScalarTens;
#endif
  using VectorTens = ExecViewUnmanaged<Scalar*[2][NP][NP][NUM_LEV]>;

  ScalarTens get_dptens () const { return m_buffers.dptens; }
  ScalarTens get_ttens ()  const { return m_buffers.ttens; }
  ScalarTensInt get_wtens ()  const { return m_buffers.wtens; }
  ScalarTens get_phitens ()  const { return m_buffers.phitens; }
  VectorTens get_vtens ()  const { return m_buffers.vtens; }

  bool process_nh_vars () const { return m_process_nh_vars; }
};

TEST_CASE("hvf", "biharmonic") {

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
  //keep nu_top=0 till nu_scale_top is set here
  params.nu_top            = 0.0; //RPDF(1e-6,1e-3)(engine);
  params.nu                = RPDF(1e-1,1e3)(engine);
  params.nu_p              = RPDF(1e-6,1e-3)(engine);
  params.nu_s              = RPDF(1e-6,1e-3)(engine);
  //do not set to 0 in the test, won't work
  params.nu_div            = RPDF(1e-6,1e-3)(engine);
  params.hypervis_scaling  = RPDF(0.1,1.0)(engine);
  params.hypervis_subcycle = IPDF(1,3)(engine);
  params.hypervis_subcycle_tom = 0;
  params.params_set = true;

  // Sync params across ranks
  MPI_Bcast(&params.nu_top,1,MPI_DOUBLE,0,c.get<Comm>().mpi_comm());
  MPI_Bcast(&params.nu,1,MPI_DOUBLE,0,c.get<Comm>().mpi_comm());
  MPI_Bcast(&params.nu_p,1,MPI_DOUBLE,0,c.get<Comm>().mpi_comm());
  MPI_Bcast(&params.nu_s,1,MPI_DOUBLE,0,c.get<Comm>().mpi_comm());
  MPI_Bcast(&params.nu_div,1,MPI_DOUBLE,0,c.get<Comm>().mpi_comm());
  //reset below, not bcasted
  MPI_Bcast(&params.hypervis_scaling,1,MPI_DOUBLE,0,c.get<Comm>().mpi_comm());
  MPI_Bcast(&params.hypervis_subcycle,1,MPI_INT,0,c.get<Comm>().mpi_comm());

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
// nu is set differently for tensor than for const hv, move this call down
// or reset nu only
  init_hv_f90(ne,hyai_ptr,hybi_ptr,hyam_ptr,hybm_ptr,dvv.data(),mp.data(),
              hvcoord.ps0,params.hypervis_subcycle,
              params.nu, params.nu_div, params.nu_top,
              params.nu_p, params.nu_s);
  ref_FE.init_mass(mp.data());
  ref_FE.init_deriv(dvv.data());

  // Create and init elements
  const int num_elems = c.get<Connectivity>().get_num_local_elements();

  auto& geo = c.create<ElementsGeometry>();
  geo.init(num_elems,false,true,PhysicalConstants::rearth0);
  geo.randomize(seed);

  auto& state = c.create<ElementsState>();
  state.init(num_elems);
  const auto max_pressure = 1000.0 + hvcoord.ps0; // This ensures max_p > ps0

  auto& derived = c.create<ElementsDerivedState>();
  derived.init(num_elems);
  derived.randomize(seed,RPDF(0.1,1.0)(engine));

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
  Real* fcor_ptr     = fcor.data();
  Real* spmp_ptr     = spmp.data();
  Real* rspmp_ptr    = rspmp.data();
  Real* tVisc_ptr    = tVisc.data();
  Real* sph2c_ptr    = sph2c.data();
  Real* mdet_ptr     = mdet.data();
  Real* minv_ptr     = minv.data();
  const Real* phis_ptr     = phis.data();
  const Real* gradphis_ptr = gradphis.data();

  // This will also init the c connectivity.
  init_geo_views_f90(d_ptr,dinv_ptr,phis_ptr,gradphis_ptr,
                     fcor_ptr, spmp_ptr,rspmp_ptr,tVisc_ptr,
                     sph2c_ptr,mdet_ptr,minv_ptr);

  Kokkos::deep_copy(geo.m_d,d);
  Kokkos::deep_copy(geo.m_dinv,dinv);
  Kokkos::deep_copy(geo.m_spheremp,spmp);
  Kokkos::deep_copy(geo.m_rspheremp,rspmp);
  Kokkos::deep_copy(geo.m_tensorvisc,tVisc);
  Kokkos::deep_copy(geo.m_vec_sph2cart,sph2c);
  Kokkos::deep_copy(geo.m_metdet,mdet);
  Kokkos::deep_copy(geo.m_metinv,minv);

  // Get or create and init other structures needed by HVF
  auto& bmm = c.create<MpiBuffersManagerMap>();
  auto& sphop = c.create<SphereOperators>();

  sphop.setup(geo,ref_FE);
  if (!bmm.is_connectivity_set ()) {
    bmm.set_connectivity(c.get_ptr<Connectivity>());
  }

  SECTION ("biharmonic_wk_theta") {
    std::cout << "Biharmonic wk theta test:\n";
    for (const bool hydrostatic : {true, false}) {
      std::cout << " -> " << (hydrostatic ? "hydrostatic" : "non-hydrostatic") << "\n";

      for (Real hv_scaling : {0.0, RPDF(0.5,5.0)(engine)}) {
        std::cout << "   -> hypervis scaling = " << hv_scaling << "\n";
        params.theta_hydrostatic_mode = hydrostatic;

        // Generate timestep settings
        const Real dt = RPDF(1.0,10.0)(engine);
        const Real eta_ave_w = RPDF(0.1,10.0)(engine);
        int np1 = IPDF(0,2)(engine);
        // Sync np1 across ranks. If they are not synced, we may get stuck in an mpi wait
        MPI_Bcast(&np1,1,MPI_INT,0,c.get<Comm>().mpi_comm());

        // Create the HVF tester
        HVFTester hvf(params,geo,state,derived);

        FunctorsBuffersManager fbm;
        fbm.request_size( hvf.requested_buffer_size() );
        fbm.allocate();
        hvf.init_buffers(fbm);

        hvf.set_timestep_data(np1,dt,eta_ave_w);

        // The be needs to be inited after the hydrostatic option has been set
        hvf.init_boundary_exchanges();

        // Update hv settings
        params.hypervis_scaling = hv_scaling;
        if (params.nu != params.nu_div) {
          Real ratio = params.nu_div / params.nu;
          if (params.hypervis_scaling != 0.0) {
            params.nu_ratio1 = ratio;
            params.nu_ratio2 = 1.0;
          }else{
            params.nu_ratio1 = ratio;
            params.nu_ratio2 = 1.0;
          }
        }else{
          params.nu_ratio1 = 1.0;
          params.nu_ratio2 = 1.0;
        }

        // Randomize inputs
        state.randomize(seed,max_pressure,hvcoord.ps0,hvcoord.hybrid_ai0,geo.m_phis);

        // Set the hv scaling
        hvf.set_hv_data(hv_scaling,params.nu_ratio1,params.nu_ratio2);

        // Run kokkos version
        hvf.biharmonic_wk_theta();

        // Run fortran version
        using ScalarStateF90    = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>;
        using ScalarStateIntF90 = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>;
        using VectorStateF90    = HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]>;

        using ScalarViewF90 = HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]>;
        using VectorViewF90 = HostViewManaged<Real*[NUM_PHYSICAL_LEV][2][NP][NP]>;

        ScalarStateF90 dp3d_f90("",num_elems);
        ScalarStateF90 vtheta_f90("",num_elems);
        ScalarStateIntF90 w_f90("",num_elems);
        ScalarStateIntF90 phi_f90("",num_elems);
        VectorStateF90 v_f90("",num_elems);

        sync_to_host(state.m_dp3d,dp3d_f90);
        sync_to_host(state.m_vtheta_dp,vtheta_f90);
        sync_to_host(state.m_w_i,w_f90);
        sync_to_host(state.m_phinh_i,phi_f90);
        sync_to_host(state.m_v,v_f90);

        const Real* dp_ptr     = reinterpret_cast<const Real*>(dp3d_f90.data());
        const Real* vtheta_ptr = reinterpret_cast<const Real*>(vtheta_f90.data());
        const Real* w_ptr      = reinterpret_cast<const Real*>(w_f90.data());
        const Real* phi_ptr    = reinterpret_cast<const Real*>(phi_f90.data());
        const Real* v_ptr      = reinterpret_cast<const Real*>(v_f90.data());

        ScalarViewF90 dptens_f90("",num_elems);
        ScalarViewF90 ttens_f90("",num_elems);
        ScalarViewF90 wtens_f90("",num_elems);
        ScalarViewF90 phitens_f90("",num_elems);
        VectorViewF90 vtens_f90("",num_elems);

        auto dptens_ptr  = dptens_f90.data();
        auto ttens_ptr   = ttens_f90.data();
        auto wtens_ptr   = wtens_f90.data();
        auto phitens_ptr = phitens_f90.data();
        auto vtens_ptr   = vtens_f90.data();
        biharmonic_wk_theta_f90(np1+1, params.hypervis_scaling, hydrostatic,
                                dp_ptr,vtheta_ptr,w_ptr,phi_ptr,v_ptr,
                                dptens_ptr,ttens_ptr,wtens_ptr, phitens_ptr,vtens_ptr);

        // Compare answers
        auto h_dptens = Kokkos::create_mirror_view(hvf.get_dptens());
        auto h_ttens = Kokkos::create_mirror_view(hvf.get_ttens());
        auto h_wtens = Kokkos::create_mirror_view(hvf.get_wtens());
        auto h_phitens = Kokkos::create_mirror_view(hvf.get_phitens());
        auto h_vtens = Kokkos::create_mirror_view(hvf.get_vtens());

        Kokkos::deep_copy(h_dptens,hvf.get_dptens());
        Kokkos::deep_copy(h_ttens,hvf.get_ttens());
        Kokkos::deep_copy(h_wtens,hvf.get_wtens());
        Kokkos::deep_copy(h_phitens,hvf.get_phitens());
        Kokkos::deep_copy(h_vtens,hvf.get_vtens());
        for (int ie=0; ie<num_elems; ++ie) {
          auto dptens_cxx = viewAsReal(Homme::subview(h_dptens,ie));
          auto ttens_cxx = viewAsReal(Homme::subview(h_ttens,ie));
          auto wtens_cxx = viewAsReal(Homme::subview(h_wtens,ie));
          auto phitens_cxx = viewAsReal(Homme::subview(h_phitens,ie));
          auto vtens_cxx = viewAsReal(Homme::subview(h_vtens,ie));
          for (int igp=0; igp<NP; ++igp) {
            for (int jgp=0; jgp<NP; ++jgp) {
              for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
                if(dptens_cxx(igp,jgp,k)!=dptens_f90(ie,k,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("dptens cxx: %3.40f\n",dptens_cxx(igp,jgp,k));
                  printf("dptens f90: %3.40f\n",dptens_f90(ie,k,igp,jgp));
                }
                REQUIRE(dptens_cxx(igp,jgp,k)==dptens_f90(ie,k,igp,jgp));
                if(ttens_cxx(igp,jgp,k)!=ttens_f90(ie,k,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("ttens cxx: %3.17f\n",ttens_cxx(igp,jgp,k));
                  printf("ttens f90: %3.17f\n",ttens_f90(ie,k,igp,jgp));
                }
                REQUIRE(ttens_cxx(igp,jgp,k)==ttens_f90(ie,k,igp,jgp));
                if(vtens_cxx(0,igp,jgp,k)!=vtens_f90(ie,k,0,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("vtens cxx: %3.40f\n",vtens_cxx(0,igp,jgp,k));
                  printf("vtens f90: %3.40f\n",vtens_f90(ie,k,0,igp,jgp));
                }
                REQUIRE(vtens_cxx(0,igp,jgp,k)==vtens_f90(ie,k,0,igp,jgp));
                if(vtens_cxx(1,igp,jgp,k)!=vtens_f90(ie,k,1,igp,jgp)) {
                  printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf("vtens cxx: %3.40f\n",vtens_cxx(1,igp,jgp,k));
                  printf("vtens f90: %3.40f\n",vtens_f90(ie,k,1,igp,jgp));
                }
                REQUIRE(vtens_cxx(1,igp,jgp,k)==vtens_f90(ie,k,1,igp,jgp));

                if (hvf.process_nh_vars()) {
                  if(wtens_cxx(igp,jgp,k)!=wtens_f90(ie,k,igp,jgp)) {
                    printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf("wtens cxx: %3.17f\n",wtens_cxx(igp,jgp,k));
                    printf("wtens f90: %3.17f\n",wtens_f90(ie,k,igp,jgp));
                  }
                  REQUIRE(wtens_cxx(igp,jgp,k)==wtens_f90(ie,k,igp,jgp));
                  if(phitens_cxx(igp,jgp,k)!=phitens_f90(ie,k,igp,jgp)) {
                    printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf("phitens cxx: %3.17f\n",phitens_cxx(igp,jgp,k));
                    printf("phitens f90: %3.17f\n",phitens_f90(ie,k,igp,jgp));
                  }
                  REQUIRE(phitens_cxx(igp,jgp,k)==phitens_f90(ie,k,igp,jgp));
                }
              }
            }
          }
        }
      }
    }
  }

  SECTION ("hypervis") {
    std::cout << "Hypervis test:\n";

    HostViewManaged<Real*[NP][NP]> phis_f90("",num_elems);
    HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]> dp_ref_f90("",num_elems);
    HostViewManaged<Real*[NUM_PHYSICAL_LEV][NP][NP]> theta_ref_f90("",num_elems);
    HostViewManaged<Real*[NUM_INTERFACE_LEV][NP][NP]> phi_ref_f90("",num_elems);

    sync_to_host(geo.m_phis, phis_f90);

    const Real* phis_ptr      = phis_f90.data();
    const Real* dp_ref_ptr    = dp_ref_f90.data();
    const Real* theta_ref_ptr = theta_ref_f90.data();
    const Real* phi_ref_ptr   = phi_ref_f90.data();

    // Have F90 compute reference states and initialize
    // in C++ RefStates. These are not dependent
    // on choice for hydrostatic mode or hv_scaling,
    // nor do the values change when state is randomized.
    initialize_reference_states_f90(phis_ptr,
                                    dp_ref_ptr,
                                    theta_ref_ptr,
                                    phi_ref_ptr);

    for (const bool hydrostatic : {true, false}) {
      std::cout << " -> " << (hydrostatic ? "hydrostatic" : "non-hydrostatic") << "\n";

      for (Real hv_scaling : {0.0, 1.2345}) {
        std::cout << "   -> hypervis scaling = " << hv_scaling << "\n";
        params.theta_hydrostatic_mode = hydrostatic;

        // Generate timestep settings
        const Real dt = RPDF(1e-5,1e-3)(engine);
        //randomize it? also, dpdiss is not tested in here
        const Real eta_ave_w = 1.0;
        int  np1 = IPDF(0,2)(engine);
        // Sync np1 across ranks. If they are not synced, we may get stuck in an mpi wait
        MPI_Bcast(&np1,1,MPI_INT,0,c.get<Comm>().mpi_comm());

        // Create the HVF tester
        HVFTester hvf(params,geo,state,derived);

        FunctorsBuffersManager fbm;
        fbm.request_size( hvf.requested_buffer_size() );
        fbm.allocate();
        hvf.init_buffers(fbm);

        hvf.set_timestep_data(np1,dt,eta_ave_w);

        // Generate random states
        state.randomize(seed);

        // The HV functor as a whole is more delicate than biharmonic_wk.
        // In particular, the EOS is used a couple of times. This means
        // that inputs *must* satisfy some minimum requirements, like
        // dp>0, vtheta>0, and d(phi)>0. This is very unlikely with random
        // inputs coming from state.randomize(seed), so we generate data
        // as "realistic" as possible, and perturb it.
        // This computation mimics that of
        // src/theta-l/share/element_ops.F90:initialize_reference_states().
        using PDF = std::uniform_real_distribution<Real>;
        ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> perturb("",num_elems);

        static constexpr Real T1 =
          PhysicalConstants::Tref_lapse_rate*PhysicalConstants::Tref*PhysicalConstants::cp/PhysicalConstants::g;
        static constexpr Real T0 = PhysicalConstants::Tref-T1;

        constexpr Real noise_lvl = 0.05;
        genRandArray(perturb,engine,PDF(-noise_lvl,noise_lvl));
        EquationOfState eos;
        eos.init(hydrostatic,hvcoord);

        ElementOps elem_ops;
        elem_ops.init(hvcoord);

        ExecViewManaged<Scalar[NUM_LEV]> buf_m("");
        ExecViewManaged<Scalar[NUM_LEV_P]> buf_i("");
        Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                             KOKKOS_LAMBDA(const TeamMember& team){
          KernelVariables kv(team);
          Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                               [&](const int idx){
            const int igp = idx / NP;
            const int jgp = idx % NP;

            auto noise = Homme::subview(perturb,kv.ie,igp,jgp);
            auto dp = Homme::subview(state.m_dp3d,kv.ie,np1,igp,jgp);
            auto theta = Homme::subview(state.m_vtheta_dp,kv.ie,np1,igp,jgp);
            auto phi = Homme::subview(state.m_phinh_i,kv.ie,np1,igp,jgp);

            // First, compute dp = dp_ref+noise
            hvcoord.compute_dp_ref(kv,state.m_ps_v(kv.ie,np1,igp,jgp),dp);
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                                 [&](const int ilev){
              dp(ilev) *= 1.0 + noise(ilev);
            });
            // Compute pressure
            elem_ops.compute_hydrostatic_p(kv,dp,buf_i,buf_m);

            // Compute vtheta_dp = theta_ref*dp, where
            // theta_ref = T0/exner + T1, exner = (p/p0)^k
            // theta_ref mimics computation in src/theta-l/share/element_ops.F90:set_theta_ref()
            Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,NUM_LEV),
                                 [&](const int ilev){
              theta(ilev) = pow(buf_m(ilev)/PhysicalConstants::p0,PhysicalConstants::kappa);
              theta(ilev) = T0/theta(ilev) + T1;
              theta(ilev) *= dp(ilev);
            });

            // Compute phi
            eos.compute_phi_i(kv,geo.m_phis(kv.ie,igp,jgp),
                              theta,buf_m,phi);
          });
        });

        // The be needs to be inited after the hydrostatic option has been set
        hvf.init_boundary_exchanges();

        // Copy states into f90 pointers
        HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][2][NP][NP]> v_f90("",num_elems);
        HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>   w_f90("",num_elems);
        HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>    dp_f90("",num_elems);
        HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_PHYSICAL_LEV][NP][NP]>    vtheta_f90("",num_elems);
        HostViewManaged<Real*[NUM_TIME_LEVELS][NUM_INTERFACE_LEV][NP][NP]>   phinh_f90("",num_elems);

        sync_to_host(state.m_v,v_f90);
        sync_to_host(state.m_w_i,w_f90);
        sync_to_host(state.m_dp3d,dp_f90);
        sync_to_host(state.m_vtheta_dp,vtheta_f90);
        sync_to_host(state.m_phinh_i,phinh_f90);

        Real* v_f90_ptr      = v_f90.data();
        Real* w_f90_ptr      = w_f90.data();
        Real* dp_f90_ptr     = dp_f90.data();
        Real* vtheta_f90_ptr = vtheta_f90.data();
        Real* phinh_f90_ptr  = phinh_f90.data();

        // Update hv settings
        params.hypervis_scaling = hv_scaling;
        if (params.nu != params.nu_div) {
          Real ratio = params.nu_div / params.nu;
          if (params.hypervis_scaling != 0.0) {
            params.nu_ratio1 = ratio;
            params.nu_ratio2 = 1.0;
          }else{
            params.nu_ratio1 = ratio;
            params.nu_ratio2 = 1.0;
          }
        }else{
          params.nu_ratio1 = 1.0;
          params.nu_ratio2 = 1.0;
        }

        // Set the viscosity params
        hvf.set_hv_data(hv_scaling,params.nu_ratio1,params.nu_ratio2);

        // Run the cxx functor
        hvf.run(np1,dt,eta_ave_w);

        // Run the f90 functor
        advance_hypervis_f90(np1+1,dt,eta_ave_w, hv_scaling, hydrostatic,
                             dp_ref_ptr, theta_ref_ptr, phi_ref_ptr,
                             v_f90_ptr, w_f90_ptr, vtheta_f90_ptr, dp_f90_ptr, phinh_f90_ptr);

        // Compare answers
        auto v_cxx      = Kokkos::create_mirror_view(state.m_v);
        auto w_cxx      = Kokkos::create_mirror_view(state.m_w_i);
        auto vtheta_cxx = Kokkos::create_mirror_view(state.m_vtheta_dp);
        auto dp_cxx     = Kokkos::create_mirror_view(state.m_dp3d);
        auto phinh_cxx  = Kokkos::create_mirror_view(state.m_phinh_i);

        Kokkos::deep_copy(v_cxx,      state.m_v);
        Kokkos::deep_copy(w_cxx,      state.m_w_i);
        Kokkos::deep_copy(vtheta_cxx, state.m_vtheta_dp);
        Kokkos::deep_copy(dp_cxx,     state.m_dp3d);
        Kokkos::deep_copy(phinh_cxx,  state.m_phinh_i);

        for (int ie=0; ie<num_elems; ++ie) {
          for (int igp=0; igp<NP; ++igp) {
            for (int jgp=0; jgp<NP; ++jgp) {
              for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
                const int ilev = k / VECTOR_SIZE;
                const int ivec = k % VECTOR_SIZE;

                if (v_cxx(ie,np1,0,igp,jgp,ilev)[ivec]!=v_f90(ie,np1,k,0,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("v_cxx: %3.40f\n",v_cxx(ie,np1,0,igp,jgp,ilev)[ivec]);
                  printf ("v_f90: %3.40f\n",v_f90(ie,np1,k,0,igp,jgp));
                }
                REQUIRE (v_cxx(ie,np1,0,igp,jgp,ilev)[ivec]==v_f90(ie,np1,k,0,igp,jgp));

                if (v_cxx(ie,np1,1,igp,jgp,ilev)[ivec]!=v_f90(ie,np1,k,1,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("v_cxx: %3.40f\n",v_cxx(ie,np1,1,igp,jgp,ilev)[ivec]);
                  printf ("v_f90: %3.40f\n",v_f90(ie,np1,k,1,igp,jgp));
                }
                REQUIRE (v_cxx(ie,np1,1,igp,jgp,ilev)[ivec]==v_f90(ie,np1,k,1,igp,jgp));

                if (dp_cxx(ie,np1,igp,jgp,ilev)[ivec]!=dp_f90(ie,np1,k,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("dp_cxx: %3.16f\n",dp_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                  printf ("dp_f90: %3.16f\n",dp_f90(ie,np1,k,igp,jgp));
                }
                REQUIRE (dp_cxx(ie,np1,igp,jgp,ilev)[ivec]==dp_f90(ie,np1,k,igp,jgp));

                if (vtheta_cxx(ie,np1,igp,jgp,ilev)[ivec]!=vtheta_f90(ie,np1,k,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("vtheta_cxx: %3.16f\n",vtheta_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                  printf ("vtheta_f90: %3.16f\n",vtheta_f90(ie,np1,k,igp,jgp));
                }
                REQUIRE (vtheta_cxx(ie,np1,igp,jgp,ilev)[ivec]==vtheta_f90(ie,np1,k,igp,jgp));

                if (hvf.process_nh_vars()) {
                  if (w_cxx(ie,np1,igp,jgp,ilev)[ivec]!=w_f90(ie,np1,k,igp,jgp)) {
                    printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf ("w_cxx: %3.16f\n",w_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                    printf ("w_f90: %3.16f\n",w_f90(ie,np1,k,igp,jgp));
                  }
                  REQUIRE (w_cxx(ie,np1,igp,jgp,ilev)[ivec]==w_f90(ie,np1,k,igp,jgp));

                  if (phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]!=phinh_f90(ie,np1,k,igp,jgp)) {
                    printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                    printf ("phinh_cxx: %3.16f\n",phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                    printf ("phinh_f90: %3.16f\n",phinh_f90(ie,np1,k,igp,jgp));
                  }
                  REQUIRE (phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]==phinh_f90(ie,np1,k,igp,jgp));
                }
              }

              if (hvf.process_nh_vars()) {
                // Last interface
                const int k = NUM_INTERFACE_LEV-1;
                const int ilev = ColInfo<NUM_INTERFACE_LEV>::LastPack;
                const int ivec = ColInfo<NUM_INTERFACE_LEV>::LastPackEnd;

                if (phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]!=phinh_f90(ie,np1,k,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("phinh_cxx: %3.16f\n",phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                  printf ("phinh_f90: %3.16f\n",phinh_f90(ie,np1,k,igp,jgp));
                }
                REQUIRE (phinh_cxx(ie,np1,igp,jgp,ilev)[ivec]==phinh_f90(ie,np1,k,igp,jgp));

                if (w_cxx(ie,np1,igp,jgp,ilev)[ivec]!=w_f90(ie,np1,k,igp,jgp)) {
                  printf ("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                  printf ("w_cxx: %3.16f\n",w_cxx(ie,np1,igp,jgp,ilev)[ivec]);
                  printf ("w_f90: %3.16f\n",w_f90(ie,np1,k,igp,jgp));
                }
                REQUIRE (w_cxx(ie,np1,igp,jgp,ilev)[ivec]==w_f90(ie,np1,k,igp,jgp));
              }
            }
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
