/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "HyperviscosityFunctorImpl.hpp"

#include "Context.hpp"
#include "FunctorsBuffersManager.hpp"
#include "profiling.hpp"

#include "mpi/BoundaryExchange.hpp"
#include "mpi/MpiBuffersManager.hpp"
#include "mpi/Connectivity.hpp"

namespace Homme
{

HyperviscosityFunctorImpl::
HyperviscosityFunctorImpl (const SimulationParams&     params,
                           const ElementsGeometry&     geometry,
                           const ElementsState&        state,
                           const ElementsDerivedState& derived)
 : m_num_elems(state.num_elems())
 , m_data (params.hypervis_subcycle,params.hypervis_subcycle_tom,
		       params.nu_ratio1,params.nu_ratio2,params.nu_top,params.nu,
		       params.nu_p,params.nu_s,params.hypervis_scaling)
 , m_state   (state)
 , m_derived (derived)
 , m_geometry (geometry)
 , m_sphere_ops (Context::singleton().get<SphereOperators>())
 , m_hvcoord (Context::singleton().get<HybridVCoord>())
 , m_policy_update_states (Homme::get_default_team_policy<ExecSpace,TagUpdateStates>(m_num_elems))
 , m_policy_first_laplace (Homme::get_default_team_policy<ExecSpace,TagFirstLaplaceHV>(m_num_elems))
 , m_policy_pre_exchange (Homme::get_default_team_policy<ExecSpace, TagHyperPreExchange>(m_num_elems))
 , m_policy_nutop_laplace (Homme::get_default_team_policy<ExecSpace, TagNutopLaplace>(m_num_elems))
 , m_policy_nutop_update_states (Homme::get_default_team_policy<ExecSpace,TagNutopUpdateStates>(m_num_elems))
 , m_tu(m_policy_update_states)
{
  init_params(params);

  // Make sure the sphere operators have buffers large enough to accommodate this functor's needs
  m_sphere_ops.allocate_buffers(m_tu);
}

HyperviscosityFunctorImpl::
HyperviscosityFunctorImpl (const int num_elems, const SimulationParams &params)
  : m_num_elems(num_elems)
  , m_data (params.hypervis_subcycle,params.hypervis_subcycle_tom,
		        params.nu_ratio1,params.nu_ratio2,params.nu_top,params.nu,
		        params.nu_p,params.nu_s,params.hypervis_scaling)
  , m_hvcoord (Context::singleton().get<HybridVCoord>())
  , m_policy_update_states (Homme::get_default_team_policy<ExecSpace,TagUpdateStates>(m_num_elems))
  , m_policy_first_laplace (Homme::get_default_team_policy<ExecSpace,TagFirstLaplaceHV>(m_num_elems))
  , m_policy_pre_exchange (Homme::get_default_team_policy<ExecSpace, TagHyperPreExchange>(m_num_elems))
  , m_policy_nutop_laplace (Homme::get_default_team_policy<ExecSpace, TagNutopLaplace>(m_num_elems))
  , m_policy_nutop_update_states (Homme::get_default_team_policy<ExecSpace,TagNutopUpdateStates>(m_num_elems))
  , m_tu(m_policy_update_states)
{
  init_params(params);
}

void HyperviscosityFunctorImpl::init_params(const SimulationParams& params)
{
  // Sanity check
  assert(params.params_set);

  if (m_data.nu_top>0) {

    m_nu_scale_top = ExecViewManaged<Scalar[NUM_LEV]>("nu_scale_top");
    ExecViewManaged<Scalar[NUM_LEV]>::HostMirror h_nu_scale_top;
    h_nu_scale_top = Kokkos::create_mirror_view(m_nu_scale_top);

    const auto etai_h = Kokkos::create_mirror_view(m_hvcoord.etai);
    const auto etam_h = Kokkos::create_mirror_view(m_hvcoord.etam);
    Kokkos::deep_copy(etai_h, m_hvcoord.etai);
    Kokkos::deep_copy(etam_h, m_hvcoord.etam);

    for (int phys_lev=0; phys_lev < NUM_LEV*VECTOR_SIZE; ++phys_lev) {
      const int ilev = phys_lev / VECTOR_SIZE;
      const int ivec = phys_lev % VECTOR_SIZE;

      Real ptop_over_press;

      //prevent padding of nu_scale to get nans to avoid last interface levels
      //of w, phi to be nans, too
      if( phys_lev < NUM_PHYSICAL_LEV ){
        //etai is num_interface_lev, that is, 129 or 73
        //etam is num_lev, so packs
        if ( etai_h(0) == 0.0) {
          ptop_over_press = etam_h(0)[0] / etam_h(ilev)[ivec];
        }else{
          ptop_over_press = etai_h(0) / etam_h(ilev)[ivec];
        }
      }else{
          ptop_over_press = 0.0;
      }

      auto val = 16.0*ptop_over_press*ptop_over_press / (ptop_over_press*ptop_over_press + 1);
      if ( val < 0.15 ) val = 0.0;
      h_nu_scale_top(ilev)[ivec] = val;

      // This is the equivalent of nlev_tom in the F90 code.
      if (val != 0) m_nu_scale_top_ilev_pack_lim = phys_lev + 1;
    }
    Kokkos::deep_copy(m_nu_scale_top, h_nu_scale_top);

    // Convert to pack index.
    m_nu_scale_top_ilev_pack_lim = ((m_nu_scale_top_ilev_pack_lim + VECTOR_SIZE - 1) /
                                    VECTOR_SIZE);
  }

  // Init ElementOps
  m_elem_ops.init(m_hvcoord);

  // Init Equation of state
  m_eos.init(params.theta_hydrostatic_mode,m_hvcoord);

#ifdef HOMMEXX_BFB_TESTING
  m_process_nh_vars = true;
#else
  m_process_nh_vars = !params.theta_hydrostatic_mode;
#endif
}

void HyperviscosityFunctorImpl::setup(const ElementsGeometry&     geometry,
                                      const ElementsState&        state,
                                      const ElementsDerivedState& derived)
{
  m_state = state;
  assert(m_num_elems == m_state.num_elems()); //Sanity check
  m_derived = derived;
  m_geometry = geometry;
  m_sphere_ops = Context::singleton().get<SphereOperators>();

  // Make sure the sphere operators have buffers large enough to accommodate this functor's needs
  m_sphere_ops.allocate_buffers(m_tu);
}

int HyperviscosityFunctorImpl::requested_buffer_size () const {
  constexpr int size_mid_scalar =   NP*NP*NUM_LEV*VECTOR_SIZE;
  constexpr int size_mid_vector = 2*NP*NP*NUM_LEV*VECTOR_SIZE;
  constexpr int size_int_scalar =   NP*NP*NUM_LEV_P*VECTOR_SIZE;

  // Number of scalar/vector int/mid buffers needed, with size nelems
  const int mid_vectors_nelems = 1;
  const int int_scalars_nelems = 0;
  const int mid_scalars_nelems = 2 + (m_process_nh_vars ? 2 : 0);

  const int size = m_num_elems*(mid_scalars_nelems*size_mid_scalar +
                                mid_vectors_nelems*size_mid_vector +
                                int_scalars_nelems*size_int_scalar);

  return size;
}

void HyperviscosityFunctorImpl::init_buffers (const FunctorsBuffersManager& fbm) {
  Errors::runtime_check(fbm.allocated_size()>=requested_buffer_size(), "Error! Buffers size not sufficient.\n");

  constexpr int size_mid_scalar =   NP*NP*NUM_LEV;
  constexpr int size_mid_vector = 2*NP*NP*NUM_LEV;

  auto mem_in = fbm.get_memory();
  Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());
  const int nelems = m_geometry.num_elems();

  // Tens quantities (persistent views => nelems)
  m_buffers.dptens = decltype(m_buffers.dptens)(mem,nelems);
  mem += size_mid_scalar*nelems;

  m_buffers.ttens = decltype(m_buffers.ttens)(mem,nelems);
  mem += size_mid_scalar*nelems;

  if (m_process_nh_vars) {
    m_buffers.wtens = decltype(m_buffers.wtens)(mem,nelems);
    mem += size_mid_scalar*nelems;

    m_buffers.phitens = decltype(m_buffers.phitens)(mem,nelems);
    mem += size_mid_scalar*nelems;
  }

  m_buffers.vtens = decltype(m_buffers.vtens)(mem,nelems);
  mem += size_mid_vector*nelems;

  const int used_mem = reinterpret_cast<Real*>(mem)-mem_in;
  if (used_mem < requested_buffer_size()) {
    printf("[HyperviscosityFunctorImpl] Warning! We used less memory than we said we would: %d instead of %d\n",
           used_mem, requested_buffer_size());
  } else if (used_mem > requested_buffer_size()) {
    std::string msg = "[HyperviscosityFunctorImpl] Error! We used more memory than we said we would: "
                    + std::to_string(used_mem) + " instead of "
                    + std::to_string(requested_buffer_size()) + "\n";
    Errors::runtime_abort(msg);
  }
}

void HyperviscosityFunctorImpl::init_boundary_exchanges () {
  auto bm_exchange = Context::singleton().get<MpiBuffersManagerMap>()[MPI_EXCHANGE];
  m_be = std::make_shared<BoundaryExchange>();
  m_be_tom = std::make_shared<BoundaryExchange>();
  std::shared_ptr<BoundaryExchange> bes[] = {m_be, m_be_tom};
  const int nlevs[] = {NUM_LEV, m_nu_scale_top_ilev_pack_lim};
  for (int i = 0; i < 2; ++i) {
    if (i == 1 && m_data.nu_top <= 0) continue;
    auto be = bes[i];
    const auto nlev = nlevs[i];
    be->set_buffers_manager(bm_exchange);
    if (m_process_nh_vars) {
      be->set_num_fields(0, 0, 6);
    } else {
      be->set_num_fields(0, 0, 4);
    }
    be->register_field(m_buffers.dptens, nlev);
    be->register_field(m_buffers.ttens, nlev);
    if (m_process_nh_vars) {
      be->register_field(m_buffers.wtens, nlev);
      be->register_field(m_buffers.phitens, nlev);
    }
    be->register_field(m_buffers.vtens, 2, 0, nlev);
    be->registration_completed();
  }
}//initBE

void HyperviscosityFunctorImpl::run (const int np1, const Real dt, const Real eta_ave_w)
{
  m_data.np1 = np1;

  m_data.dt = dt;
  if (m_data.hypervis_subcycle > 0) { 
    m_data.dt_hvs = dt/m_data.hypervis_subcycle;
  }else{
    //won't be used
    m_data.dt_hvs = -1.0;
  }
  if (m_data.hypervis_subcycle_tom > 0) { 
    m_data.dt_hvs_tom = dt/m_data.hypervis_subcycle_tom;
  }else{
    //won't be used
    m_data.dt_hvs_tom = -1.0;
  }
  m_data.eta_ave_w = eta_ave_w;

  // Convert vtheta_dp -> theta
  auto state = m_state;
  Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(state.num_elems()),
                       KOKKOS_LAMBDA(const TeamMember& team) {
    const int ie = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      // theta->vtheta
      auto vtheta = Homme::subview(state.m_vtheta_dp,ie,np1,igp,jgp);
      auto dp = Homme::subview(state.m_dp3d,ie,np1,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,NUM_LEV),
                           [&](const int ilev) {
        vtheta(ilev) /= dp(ilev);
      });
    });
  });
  Kokkos::fence();

  for (int icycle = 0; icycle < m_data.hypervis_subcycle; ++icycle) {
    GPTLstart("hvf-bhwk");
    biharmonic_wk_theta ();
    GPTLstop("hvf-bhwk");

    Kokkos::parallel_for(m_policy_pre_exchange, *this);
    Kokkos::fence();

    // Exchange
    assert (m_be->is_registration_completed());
    GPTLstart("hvf-bexch");
    m_be->exchange();
    GPTLstop("hvf-bexch");

    // Update states
    Kokkos::parallel_for(m_policy_update_states, *this);
    Kokkos::fence();
  } //subcycle

  // Convert theta back to vtheta, and adjust w at surface
  auto geo = m_geometry;
  auto process_nh_vars = m_process_nh_vars;
  Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(state.num_elems()),
                       KOKKOS_LAMBDA(const TeamMember& team) {
    const int ie = team.league_rank();
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,NP*NP),
                         [&](const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      // theta->vtheta
      auto vtheta = Homme::subview(state.m_vtheta_dp,ie,np1,igp,jgp);
      auto dp = Homme::subview(state.m_dp3d,ie,np1,igp,jgp);
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,NUM_LEV),
                           [&](const int ilev) {
        vtheta(ilev) *= dp(ilev);
      });

      // Fix w at surface:
      // Adjust w_i at the surface, since velocity has changed
      if (process_nh_vars) {
        Kokkos::single(Kokkos::PerThread(team),[&](){
          using InfoI = ColInfo<NUM_INTERFACE_LEV>;
          using InfoM = ColInfo<NUM_PHYSICAL_LEV>;
          constexpr int LAST_MID_PACK     = InfoM::LastPack;
          constexpr int LAST_MID_PACK_END = InfoM::LastPackEnd;
          constexpr int LAST_INT_PACK     = InfoI::LastPack;
          constexpr int LAST_INT_PACK_END = InfoI::LastPackEnd;
          constexpr Real g = PhysicalConstants::g;

          const auto& grad_x = geo.m_gradphis(ie,0,igp,jgp);
          const auto& grad_y = geo.m_gradphis(ie,1,igp,jgp);
          const auto& u = state.m_v(ie,np1,0,igp,jgp,LAST_MID_PACK)[LAST_MID_PACK_END];
          const auto& v = state.m_v(ie,np1,1,igp,jgp,LAST_MID_PACK)[LAST_MID_PACK_END];

          auto& w = state.m_w_i(ie,np1,igp,jgp,LAST_INT_PACK)[LAST_INT_PACK_END];

          w = (u*grad_x+v*grad_y) / g;
        });
      }
    });
  });//conversion back to vtheta

  Kokkos::fence();

  // sponge layer 
  if (m_data.nu_top > 0) {
    for (int icycle = 0; icycle < m_data.hypervis_subcycle_tom; ++icycle) {
      // laplace(fields) --> ttens, etc.
      Kokkos::parallel_for(m_policy_nutop_laplace, *this);
      Kokkos::fence();

      // exchange is done on ttens, dptens, vtens, etc.
      assert (m_be->is_registration_completed());
      GPTLstart("hvf-bexch");
      m_be_tom->exchange();
      GPTLstop("hvf-bexch");

      Kokkos::parallel_for(m_policy_nutop_update_states, *this);
      Kokkos::fence();
    }
  } // for sponge layer
} // run()

void HyperviscosityFunctorImpl::biharmonic_wk_theta() const
{
  // For the first laplacian we use a differnt kernel, which uses directly the states
  // at timelevel np1 as inputs, and subtracts the reference states.
  // This way we avoid copying the states to *tens buffers.
  Kokkos::parallel_for(m_policy_first_laplace, *this);
  Kokkos::fence();

  // Exchange
  assert (m_be->is_registration_completed());
  GPTLstart("hvf-bexch");
  m_be->exchange(m_geometry.m_rspheremp);
  GPTLstop("hvf-bexch");

  // Compute second laplacian, tensor or const hv
  const int ne = m_geometry.num_elems();
  if ( m_data.consthv ) {
    auto policy = Homme::get_default_team_policy<ExecSpace,TagSecondLaplaceConstHV>(ne);
    Kokkos::parallel_for(policy, *this);
  }else{
    auto policy = Homme::get_default_team_policy<ExecSpace,TagSecondLaplaceTensorHV>(ne);
    Kokkos::parallel_for(policy, *this);
  }
  Kokkos::fence();
} //biharmonic

// Laplace for nu_top
KOKKOS_INLINE_FUNCTION
void HyperviscosityFunctorImpl::operator() (const TagNutopLaplace&, const TeamMember& team) const {
  KernelVariables kv(team, m_tu);

  using MidColumn = decltype(Homme::subview(m_buffers.wtens,0,0,0));

  // Laplacian of layer thickness
  m_sphere_ops.laplace_simple(kv,
                              Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1),
                              Homme::subview(m_buffers.dptens,kv.ie),
                              m_nu_scale_top_ilev_pack_lim);
  // Laplacian of theta
  m_sphere_ops.laplace_simple(kv,
                              Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1),
                              Homme::subview(m_buffers.ttens,kv.ie),
                              m_nu_scale_top_ilev_pack_lim);

  if (m_process_nh_vars) {
    // Laplacian of vertical velocity (do not compute last interface)
    m_sphere_ops.laplace_simple<NUM_LEV,NUM_LEV_P>(kv,
                                                   Homme::subview(m_state.m_w_i,kv.ie,m_data.np1),
                                                   Homme::subview(m_buffers.wtens,kv.ie),
                                                   m_nu_scale_top_ilev_pack_lim);
    // Laplacian of geopotential (do not compute last interface)
    m_sphere_ops.laplace_simple<NUM_LEV,NUM_LEV_P>(kv,
                                                   Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1),
                                                   Homme::subview(m_buffers.phitens,kv.ie),
                                                   m_nu_scale_top_ilev_pack_lim);
  }

  // Laplacian of velocity
  m_sphere_ops.vlaplace_sphere_wk_contra(kv, m_data.nu_ratio1,
                                         Homme::subview(m_state.m_v,kv.ie,m_data.np1),
                                         Homme::subview(m_buffers.vtens,kv.ie),
                                         m_nu_scale_top_ilev_pack_lim);

  kv.team_barrier();

  Kokkos::parallel_for(
    Kokkos::TeamThreadRange(kv.team,NP*NP),
    [&] (const int idx) {
      const int igp = idx / NP;
      const int jgp = idx % NP;

      const auto utens  = Homme::subview(m_buffers.vtens,kv.ie,0,igp,jgp);
      const auto vtens  = Homme::subview(m_buffers.vtens,kv.ie,1,igp,jgp);
      const auto ttens  = Homme::subview(m_buffers.ttens,kv.ie,igp,jgp);
      const auto dptens = Homme::subview(m_buffers.dptens,kv.ie,igp,jgp);
     
      MidColumn wtens, phitens;
      if (m_process_nh_vars) {
        wtens   = Homme::subview(m_buffers.wtens,kv.ie,igp,jgp);
        phitens = Homme::subview(m_buffers.phitens,kv.ie,igp,jgp);
      }

      Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(kv.team, m_nu_scale_top_ilev_pack_lim),
        [&] (const int ilev) {
          
          const auto xf = m_data.dt_hvs_tom  * m_nu_scale_top(ilev) * m_data.nu_top;
          utens(ilev)  *= xf;
          vtens(ilev)  *= xf;
          ttens(ilev)  *= xf;
          dptens(ilev) *= xf;

          if (m_process_nh_vars) {
            wtens(ilev)   *= xf;
            phitens(ilev) *= xf;
          }

        }); // threadvectorrange
    }); // teamthreadrange
} // TagNutopLaplace

KOKKOS_INLINE_FUNCTION
void HyperviscosityFunctorImpl::operator() (const TagNutopUpdateStates&, const TeamMember& team) const {
  KernelVariables kv(team, m_tu);

  using MidColumn = decltype(Homme::subview(m_buffers.wtens,0,0,0));
  using IntColumn = decltype(Homme::subview(m_state.m_w_i,0,0,0,0));

  Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                       [&](const int idx) {
    const int igp = idx / NP;
    const int jgp = idx % NP;

    // Add Xtens quantities back to the states, except for vtheta
    auto u = Homme::subview(m_state.m_v,kv.ie,m_data.np1,0,igp,jgp);
    auto v = Homme::subview(m_state.m_v,kv.ie,m_data.np1,1,igp,jgp);
    auto vtheta = Homme::subview(m_state.m_vtheta_dp,kv.ie,m_data.np1,igp,jgp);
    auto dp     = Homme::subview(m_state.m_dp3d,kv.ie,m_data.np1,igp,jgp);

    auto utens   = Homme::subview(m_buffers.vtens,kv.ie,0,igp,jgp);
    auto vtens   = Homme::subview(m_buffers.vtens,kv.ie,1,igp,jgp);
    auto ttens   = Homme::subview(m_buffers.ttens,kv.ie,igp,jgp);
    auto dptens  = Homme::subview(m_buffers.dptens,kv.ie,igp,jgp);
    const auto& rspheremp = m_geometry.m_rspheremp(kv.ie,igp,jgp);

    MidColumn wtens, phitens;
    IntColumn w, phi_i;

    if (m_process_nh_vars) {
      wtens   = Homme::subview(m_buffers.wtens,kv.ie,igp,jgp);
      phitens = Homme::subview(m_buffers.phitens,kv.ie,igp,jgp);
      w       = Homme::subview(m_state.m_w_i,kv.ie,m_data.np1,igp,jgp);
      phi_i   = Homme::subview(m_state.m_phinh_i,kv.ie,m_data.np1,igp,jgp);
    }

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, m_nu_scale_top_ilev_pack_lim),
                         [&](const int ilev) {
      utens(ilev)   *= rspheremp;
      vtens(ilev)   *= rspheremp;
      ttens(ilev)   *= rspheremp;
      dptens(ilev)  *= rspheremp;
      u(ilev)      += utens(ilev);
      v(ilev)      += vtens(ilev);
      vtheta(ilev) += ttens(ilev);
      dp(ilev)     += dptens(ilev);

      if (m_process_nh_vars) {
        wtens(ilev)   *= rspheremp;
        phitens(ilev) *= rspheremp;
        w(ilev)     += wtens(ilev);
        phi_i(ilev) += phitens(ilev);
      }
    }); // threadvectorrange
  }); // threadteamrange
} // tagUpdateStates2

} // namespace Homme
