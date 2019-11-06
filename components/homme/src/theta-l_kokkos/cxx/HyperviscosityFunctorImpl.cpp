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
 : m_data (params.hypervis_subcycle,params.nu_ratio1,params.nu_ratio2,params.nu_top,params.nu,params.nu_p,params.nu_s,params.hypervis_scaling)
 , m_state   (state)
 , m_derived (derived)
 , m_geometry (geometry)
 , m_sphere_ops (Context::singleton().get<SphereOperators>())
 , m_hvcoord (Context::singleton().get<HybridVCoord>())
 , m_policy_ref_states (Homme::get_default_team_policy<ExecSpace,TagRefStates>(state.num_elems()))
 , m_policy_update_states (Homme::get_default_team_policy<ExecSpace,TagUpdateStates>(state.num_elems()))
 , m_policy_first_laplace (Homme::get_default_team_policy<ExecSpace,TagFirstLaplaceHV>(state.num_elems()))
 , m_policy_pre_exchange (Homme::get_default_team_policy<ExecSpace, TagHyperPreExchange>(state.num_elems()))
{
  // Sanity check
  assert(params.params_set);

  if (m_data.nu_top>0) {

    m_nu_scale_top = ExecViewManaged<Scalar[NUM_LEV]>("nu_scale_top");
    ExecViewManaged<Scalar[NUM_LEV]>::HostMirror h_nu_scale_top;
    h_nu_scale_top = Kokkos::create_mirror_view(m_nu_scale_top);

    Kokkos::Array<Real,NUM_BIHARMONIC_PHYSICAL_LEVELS> lev_nu_scale_top = { 4.0, 2.0, 1.0 };
    for (int phys_lev=0; phys_lev<NUM_BIHARMONIC_PHYSICAL_LEVELS; ++phys_lev) {
      const int ilev = phys_lev / VECTOR_SIZE;
      const int ivec = phys_lev % VECTOR_SIZE;
      h_nu_scale_top(ilev)[ivec] = lev_nu_scale_top[phys_lev]*m_data.nu_top;
    }
    Kokkos::deep_copy(m_nu_scale_top, h_nu_scale_top);
  }

  // Init ElementOps
  m_elem_ops.init(m_hvcoord);

  // Init Equation of state
  m_theta_hydrostatic_mode = params.theta_hydrostatic_mode;
  m_eos.init(m_theta_hydrostatic_mode,m_hvcoord);

  // Make sure the sphere operators have buffers large enough to accommodate this functor's needs
  m_sphere_ops.allocate_buffers(Homme::get_default_team_policy<ExecSpace>(m_state.num_elems()));
}

int HyperviscosityFunctorImpl::requested_buffer_size () const {
  constexpr int size_mid_scalar =   NP*NP*NUM_LEV*VECTOR_SIZE;
  constexpr int size_mid_vector = 2*NP*NP*NUM_LEV*VECTOR_SIZE;
  constexpr int size_int_scalar =   NP*NP*NUM_LEV_P*VECTOR_SIZE;
  constexpr int size_bhm_scalar =   NP*NP*NUM_BIHARMONIC_LEV*VECTOR_SIZE;
  constexpr int size_bhm_vector = 2*NP*NP*NUM_BIHARMONIC_LEV*VECTOR_SIZE;

  const int nelems = m_geometry.num_elems();
  const int nteams = get_num_concurrent_teams(m_policy_pre_exchange); 

  // Number of scalar/vector int/mid/bhm buffers needed, with size nteams/nelems
  constexpr int mid_vectors_nelems = 1;
  int int_scalars_nelems = 1;
  int mid_scalars_nelems = 5;
#ifdef XX_NONBFB_COMING
  // wtens on interfaces
  ++int_scalars_nelems;
  --mid_scalars_nelems;
#endif

#ifdef HV_USE_THETA_REF
  ++mid_scalars_nelems;
#endif

  constexpr int bhm_scalars_nteams = 4;
  constexpr int bhm_vectors_nteams = 1;
  constexpr int int_scalars_nteams = 1;
  constexpr int mid_scalars_nteams = 1;

  const int size = nelems*(mid_scalars_nelems*size_mid_scalar +
                           mid_vectors_nelems*size_mid_vector +
                           int_scalars_nelems*size_int_scalar) +
                   nteams*(bhm_scalars_nteams*size_bhm_scalar +
                           bhm_vectors_nteams*size_bhm_vector +
                           mid_scalars_nteams*size_mid_scalar +
                           int_scalars_nteams*size_int_scalar);

  return size;
}

void HyperviscosityFunctorImpl::init_buffers (const FunctorsBuffersManager& fbm) {
  Errors::runtime_check(fbm.allocated_size()>=requested_buffer_size(), "Error! Buffers size not sufficient.\n");

  constexpr int size_mid_scalar =   NP*NP*NUM_LEV;
  constexpr int size_mid_vector = 2*NP*NP*NUM_LEV;
  constexpr int size_int_scalar =   NP*NP*NUM_LEV_P;
  constexpr int size_bhm_scalar =   NP*NP*NUM_BIHARMONIC_LEV;
  constexpr int size_bhm_vector = 2*NP*NP*NUM_BIHARMONIC_LEV;

  auto mem_in = fbm.get_memory();
  Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());
  const int nelems = m_geometry.num_elems();
  const int nteams = get_num_concurrent_teams(m_policy_pre_exchange); 

  // Ref states (persistent views => nelems)
  m_buffers.dp_ref = decltype(m_buffers.dp_ref)(mem,nelems);
  mem += size_mid_scalar*nelems;

#ifdef HV_USE_THETA_REF
  m_buffers.theta_ref = decltype(m_buffers.theta_ref)(mem,nelems);
  mem += size_mid_scalar*nelems;
#endif

  m_buffers.phi_i_ref = decltype(m_buffers.phi_i_ref)(mem,nelems);
  mem += size_int_scalar*nelems;

  // Tens quantities (persistent views => nelems)
  m_buffers.dptens = decltype(m_buffers.dptens)(mem,nelems);
  mem += size_mid_scalar*nelems;

  m_buffers.ttens = decltype(m_buffers.ttens)(mem,nelems);
  mem += size_mid_scalar*nelems;

  m_buffers.wtens = decltype(m_buffers.wtens)(mem,nelems);
#ifdef XX_NONBFB_COMING
  mem += size_int_scalar*nelems;
#else
  mem += size_mid_scalar*nelems;
#endif

  m_buffers.phitens = decltype(m_buffers.phitens)(mem,nelems);
  mem += size_mid_scalar*nelems;

  m_buffers.vtens = decltype(m_buffers.vtens)(mem,nelems);
  mem += size_mid_vector*nelems;

  // Biharmonic views (non-persistent views => nteams)
  m_buffers.lapl_dp = decltype(m_buffers.lapl_dp)(mem,nteams);
  mem += size_bhm_scalar*nteams;

  m_buffers.lapl_w = decltype(m_buffers.lapl_w)(mem,nteams);
  mem += size_bhm_scalar*nteams;

  m_buffers.lapl_phi = decltype(m_buffers.lapl_phi)(mem,nteams);
  mem += size_bhm_scalar*nteams;

  m_buffers.lapl_theta = decltype(m_buffers.lapl_theta)(mem,nteams);
  mem += size_bhm_scalar*nteams;

  m_buffers.lapl_v = decltype(m_buffers.lapl_v)(mem,nteams);
  mem += size_bhm_vector*nteams;

  // Service views (non-persistent views => nteams)
  m_buffers.p = decltype(m_buffers.p)(mem,nteams);
  mem += size_mid_scalar*nteams;

  m_buffers.p_i = decltype(m_buffers.p_i)(mem,nteams);
  mem += size_int_scalar*nteams;

  // ps_ref can alias anything (except dp_ref), since it's used to compute dp_ref, then tossed
  m_buffers.ps_ref = decltype(m_buffers.ps_ref)(reinterpret_cast<Real*>(m_buffers.p.data()),nteams);

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
  m_be = std::make_shared<BoundaryExchange>();
  auto bm_exchange = Context::singleton().get<MpiBuffersManagerMap>()[MPI_EXCHANGE];

  m_be->set_buffers_manager(bm_exchange);
  if (m_theta_hydrostatic_mode) {
    m_be->set_num_fields(0, 0, 4);
  } else {
#ifdef XX_NONBFB_COMING
    m_be->set_num_fields(0, 0, 5, 1);
#else
    m_be->set_num_fields(0, 0, 6);
#endif
  }
  m_be->register_field(m_buffers.dptens);
  m_be->register_field(m_buffers.ttens);
  if (!m_theta_hydrostatic_mode) {
    m_be->register_field(m_buffers.wtens);
    m_be->register_field(m_buffers.phitens);
  }
  m_be->register_field(m_buffers.vtens, 2, 0);
  m_be->registration_completed();
}

void HyperviscosityFunctorImpl::run (const int np1, const Real dt, const Real eta_ave_w)
{
  m_data.np1 = np1;
  m_data.dt = dt/m_data.hypervis_subcycle;
  m_data.eta_ave_w = eta_ave_w;

  Kokkos::parallel_for(m_policy_ref_states,*this);

  for (int icycle = 0; icycle < m_data.hypervis_subcycle; ++icycle) {
    GPTLstart("hvf-bhwk");
    biharmonic_wk_theta ();
    GPTLstop("hvf-bhwk");

    // dispatch parallel_for for first kernel
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
  }

  // Finally, convert theta back to vtheta, and adjust w at surface
  auto state = m_state;
  auto geo = m_geometry;
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
      Kokkos::single(Kokkos::PerThread(team),[&](){
        using InfoI = ColInfo<NUM_INTERFACE_LEV>;
        using InfoM = ColInfo<NUM_PHYSICAL_LEV>;
        constexpr Real g = PhysicalConstants::g;

        const auto& grad_x = geo.m_gradphis(ie,0,igp,jgp);
        const auto& grad_y = geo.m_gradphis(ie,1,igp,jgp);
        const auto& u = state.m_v(ie,np1,0,igp,jgp,InfoM::LastPack)[InfoM::LastVecEnd];
        const auto& v = state.m_v(ie,np1,1,igp,jgp,InfoM::LastPack)[InfoM::LastVecEnd];

        auto& w = state.m_w_i(ie,np1,igp,jgp,InfoI::LastPack)[InfoI::LastVecEnd];

        w = (u*grad_x+v*grad_y) / g;
      });
    });
  });
}

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
}

} // namespace Homme
