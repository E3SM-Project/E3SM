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

#include "utilities/TestUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

using namespace Homme;

extern "C" {
void init_f90 (const int& ne, const Real* hyai_ptr, const Real* dvv, const Real* mp,
               const Real& ps0, const Real& hypervis_subcycle,
               const Real& nu, const Real& nu_div, const bool& hydrostatic);
void init_geo_views_f90 (const Real*& d_ptr,const Real*& dinv_ptr,
               const Real*& sphmp_ptr, const Real*& rspmp_ptr,
               const Real*& tVisc_ptr, const Real*& sph2c_ptr,
               const Real*& metdet_ptr, const Real*& metinv_ptr);
void biharmonic_wk_theta_f90 (const int& np1, const Real& hv_scaling,
                              const Real*& dp, const Real*& vtheta_dp,
                              const Real*& w,  const Real*& phi, const Real*& v,
                              Real*& dptens, Real*& ttens, Real*& wtens,
                              Real*& phitens, Real*& vtens);
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
  {}

  ~HVFTester () = default;

  void set_timestep_data (const int np1, const Real dt, const Real eta_ave_w)
  {
    m_data.np1 = np1;
    m_data.dt = dt/m_data.hypervis_subcycle;
    m_data.eta_ave_w = eta_ave_w;
  }

  void set_hv_data (const Real hv_scaling, const Real nu_ratio1, const Real nu_ratio2)
  {
    m_data.hypervis_scaling = hv_scaling;
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
  const int seed = 1; // Change to the following line after debugging
  // const int seed = rd();
  rngAlg engine(seed);
  std::uniform_real_distribution<Real> rpdf(0.01, 1.0);
  std::uniform_int_distribution<int> nuexp(10, 17);
  std::uniform_int_distribution<int> ipdf(1, 3);

  // Use stuff from Context, to increase similarity with actual runs
  auto& c = Context::singleton();
  auto old_comm = c.get_ptr<Comm>();

  // Init parameters
  auto& params = c.create<SimulationParams>();
  params.nu_top            = rpdf(engine)*std::pow(10.0,nuexp(engine));
  params.nu                = rpdf(engine)*std::pow(10.0,nuexp(engine));
  params.nu_p              = rpdf(engine)*std::pow(10.0,nuexp(engine));
  params.nu_s              = rpdf(engine)*std::pow(10.0,nuexp(engine));
  params.nu_div            = rpdf(engine)*std::pow(10.0,nuexp(engine));
  params.hypervis_scaling  = ipdf(engine);
  params.hypervis_subcycle = ipdf(engine);
  params.theta_hydrostatic_mode = false; //(rpdf(engine) > 0.5);
  params.params_set = true;

  // Create and init hvcoord and ref_elem, needed to init the fortran interface
  auto& hvcoord = c.create<HybridVCoord>();
  auto& ref_FE  = c.create<ReferenceElement>();
  hvcoord.random_init(seed);
  ref_FE.random_init(seed);

  auto hyai = Kokkos::create_mirror_view(hvcoord.hybrid_ai);
  auto dvv = Kokkos::create_mirror_view(ref_FE.get_deriv());
  auto mp = Kokkos::create_mirror_view(ref_FE.get_mass());
  Kokkos::deep_copy(hyai,hvcoord.hybrid_ai);
  Kokkos::deep_copy(dvv,ref_FE.get_deriv());
  Kokkos::deep_copy(mp,ref_FE.get_mass());
  const Real* hyai_ptr  = hyai.data();
  const Real* dvv_ptr   = dvv.data();
  const Real* mp_ptr    = mp.data();

  // This will also init the c connectivity.
  init_f90(ne,hyai_ptr,dvv_ptr,mp_ptr,
           hvcoord.ps0,params.hypervis_subcycle,
           params.nu, params.nu_div, params.theta_hydrostatic_mode);

  // Create and init elements
  const int num_elems = c.get<Connectivity>().get_num_local_elements();

  auto& geo = c.create<ElementsGeometry>();
  geo.init(num_elems,false, /* alloc_gradphis = */ true);
  geo.randomize(seed);

  auto& state = c.create<ElementsState>();
  state.init(num_elems);
  const auto max_pressure = rpdf(engine) + hvcoord.ps0; // This ensures max_p > ps0
  state.randomize(seed,max_pressure,hvcoord.ps0);

  auto& derived = c.create<ElementsDerivedState>();
  derived.init(num_elems);
  derived.randomize(seed,rpdf(engine));

  // Init f90
  auto d     = Kokkos::create_mirror_view(geo.m_d);
  auto dinv  = Kokkos::create_mirror_view(geo.m_dinv);
  auto spmp = Kokkos::create_mirror_view(geo.m_spheremp);
  auto rspmp = Kokkos::create_mirror_view(geo.m_rspheremp);
  auto tVisc = Kokkos::create_mirror_view(geo.m_tensorvisc);
  auto sph2c = Kokkos::create_mirror_view(geo.m_vec_sph2cart);
  auto mdet = Kokkos::create_mirror_view(geo.m_metdet);
  auto minv = Kokkos::create_mirror_view(geo.m_metinv);

  Kokkos::deep_copy(d,geo.m_d);
  Kokkos::deep_copy(dinv,geo.m_dinv);
  Kokkos::deep_copy(spmp,geo.m_spheremp);
  Kokkos::deep_copy(rspmp,geo.m_rspheremp);
  Kokkos::deep_copy(tVisc,geo.m_tensorvisc);
  Kokkos::deep_copy(sph2c,geo.m_vec_sph2cart);
  Kokkos::deep_copy(mdet,geo.m_metdet);
  Kokkos::deep_copy(minv,geo.m_metinv);

  const Real* d_ptr     = d.data();
  const Real* dinv_ptr  = dinv.data();
  const Real* spmp_ptr  = spmp.data();
  const Real* rspmp_ptr = rspmp.data();
  const Real* tVisc_ptr = tVisc.data();
  const Real* sph2c_ptr = sph2c.data();
  const Real* mdet_ptr  = mdet.data();
  const Real* minv_ptr  = minv.data();

  // This will also init the c connectivity.
  init_geo_views_f90(d_ptr,dinv_ptr,
                     spmp_ptr,rspmp_ptr,tVisc_ptr,
                     sph2c_ptr,mdet_ptr,minv_ptr);

  // Get or create and init other structures needed by HVF
  auto& bmm = c.create<MpiBuffersManagerMap>();
  auto& sphop = c.create<SphereOperators>();
  FunctorsBuffersManager fbm;

  sphop.setup(geo,ref_FE);
  if (!bmm.is_connectivity_set ()) {
    bmm.set_connectivity(c.get_ptr<Connectivity>());
  }

  // Create the HVF tester
  HVFTester hvf(params,geo,state,derived);
  fbm.request_size( hvf.requested_buffer_size() );
  fbm.allocate();
  hvf.init_buffers(fbm);
  hvf.init_boundary_exchanges();

  SECTION ("biharmonic_wk_theta") {
    for (Real hv_scaling : {0.0, rpdf(engine)}) {

      // Generate timestep settings
      const Real dt = 1.0;//rpdf(engine);
      const Real eta_ave_w = 1.0;//rpdf(engine);
      const int  np1 = 1;//ipdf(engine);
      hvf.set_timestep_data(np1-1,dt,eta_ave_w);

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
      biharmonic_wk_theta_f90(np1, params.hypervis_scaling,
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
                printf("hv_scaling: %3.17f\n",hv_scaling);
                printf("dptens cxx: %3.40f\n",dptens_cxx(igp,jgp,k));
                printf("dptens f90: %3.40f\n",dptens_f90(ie,k,igp,jgp));
              }
              REQUIRE(dptens_cxx(igp,jgp,k)==dptens_f90(ie,k,igp,jgp));
              if(ttens_cxx(igp,jgp,k)!=ttens_f90(ie,k,igp,jgp)) {
                printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                printf("hv_scaling: %3.17f\n",hv_scaling);
                printf("ttens cxx: %3.17f\n",ttens_cxx(igp,jgp,k));
                printf("ttens f90: %3.17f\n",ttens_f90(ie,k,igp,jgp));
              }
              REQUIRE(ttens_cxx(igp,jgp,k)==ttens_f90(ie,k,igp,jgp));
              if(wtens_cxx(igp,jgp,k)!=wtens_f90(ie,k,igp,jgp)) {
                printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                printf("hv_scaling: %3.17f\n",hv_scaling);
                printf("wtens cxx: %3.17f\n",wtens_cxx(igp,jgp,k));
                printf("wtens f90: %3.17f\n",wtens_f90(ie,k,igp,jgp));
              }
              REQUIRE(wtens_cxx(igp,jgp,k)==wtens_f90(ie,k,igp,jgp));
              if(phitens_cxx(igp,jgp,k)!=phitens_f90(ie,k,igp,jgp)) {
                printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                printf("hv_scaling: %3.17f\n",hv_scaling);
                printf("phitens cxx: %3.17f\n",phitens_cxx(igp,jgp,k));
                printf("phitens f90: %3.17f\n",phitens_f90(ie,k,igp,jgp));
              }
              REQUIRE(phitens_cxx(igp,jgp,k)==phitens_f90(ie,k,igp,jgp));

              if(vtens_cxx(0,igp,jgp,k)!=vtens_f90(ie,k,0,igp,jgp)) {
                printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                printf("hv_scaling: %3.17f\n",hv_scaling);
                printf("vtens cxx: %3.17f\n",vtens_cxx(0,igp,jgp,k));
                printf("vtens f90: %3.17f\n",vtens_f90(ie,k,0,igp,jgp));
              }
              REQUIRE(vtens_cxx(0,igp,jgp,k)==vtens_f90(ie,k,0,igp,jgp));
              if(vtens_cxx(1,igp,jgp,k)!=vtens_f90(ie,k,1,igp,jgp)) {
                printf("ie,k,igp,jgp: %d, %d, %d, %d\n",ie,k,igp,jgp);
                printf("hv_scaling: %3.17f\n",hv_scaling);
                printf("vtens cxx: %3.17f\n",vtens_cxx(1,igp,jgp,k));
                printf("vtens f90: %3.17f\n",vtens_f90(ie,k,1,igp,jgp));
              }
              REQUIRE(vtens_cxx(1,igp,jgp,k)==vtens_f90(ie,k,1,igp,jgp));
            }
          }
        }
      }
    }
  }

  // TODO
  // SECTION ("run") {
  //   // Store the current state as 'initial' state, then re-init the state
  //   elements.m_state.save_state();
  // }

  c.finalize_singleton();
  auto& new_comm = c.create<Comm>();
  new_comm = *old_comm;
  cleanup_f90();
}
