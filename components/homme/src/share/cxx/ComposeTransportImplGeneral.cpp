/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#include "ComposeTransportImpl.hpp"
#include "compose_hommexx.hpp"

extern "C" void
sl_get_params(double* nu_q, double* hv_scaling, int* hv_q, int* hv_subcycle_q,
              int* limiter_option, int* cdr_check, int* geometry_type,
              int* trajectory_nsubstep);

namespace Homme {

static int calc_nslot (const int nelemd) {
  const auto tp = Homme::get_default_team_policy<ExecSpace>(nelemd);
  const auto tu = TeamUtils<ExecSpace>(tp);
  return std::min(nelemd, tu.get_num_ws_slots());
}

ComposeTransportImpl::ComposeTransportImpl ()
  : m_tp_ne(1,1,1), m_tp_ne_qsize(1,1,1), m_tp_ne_hv_q(1,1,1), // throwaway settings
    m_tu_ne(m_tp_ne), m_tu_ne_qsize(m_tp_ne_qsize), m_tu_ne_hv_q(m_tp_ne_hv_q)
{
  setup();
}

ComposeTransportImpl::ComposeTransportImpl (const int num_elems)
  : m_tp_ne(1,1,1), m_tp_ne_qsize(1,1,1), m_tp_ne_hv_q(1,1,1), // throwaway settings
    m_tu_ne(m_tp_ne), m_tu_ne_qsize(m_tp_ne_qsize), m_tu_ne_hv_q(m_tp_ne_hv_q)
{
  nslot = calc_nslot(m_geometry.num_elems());
}

void ComposeTransportImpl::setup () {
  m_hvcoord = Context::singleton().get<HybridVCoord>();
  m_elements = Context::singleton().get<Elements>();
  m_state = m_elements.m_state;
  m_derived = m_elements.m_derived;
  m_geometry = Context::singleton().get<ElementsGeometry>();
  m_tracers = Context::singleton().get<Tracers>();
  m_sphere_ops = Context::singleton().get<SphereOperators>();
  
  set_dp_tol();
  setup_enhanced_trajectory();
  
  nslot = calc_nslot(m_geometry.num_elems());
}

void ComposeTransportImpl::reset (const SimulationParams& params) {
  const auto num_elems = Context::singleton().get<Connectivity>().get_num_local_elements();

  const bool independent_time_steps = params.dt_tracer_factor > params.dt_remap_factor;

  sl_get_params(&m_data.nu_q, &m_data.hv_scaling, &m_data.hv_q, &m_data.hv_subcycle_q,
                &m_data.limiter_option, &m_data.cdr_check, &m_data.geometry_type,
                &m_data.trajectory_nsubstep);

  if (independent_time_steps != m_data.independent_time_steps ||
      m_data.nelemd != num_elems || m_data.qsize != params.qsize) {
    const auto& g = m_geometry;
    const auto& t = m_tracers;
    const auto& s = m_state;
    const auto& d = m_derived;
    const auto nel = num_elems;
    const auto nlev = NUM_LEV*packn;
    const int ndim = (m_data.trajectory_nsubstep == 0 ?
                      3 :
                      (independent_time_steps ? 4 : 3));
    m_data.dep_pts = DeparturePoints("dep_pts", nel, num_phys_lev, np, np, ndim);
    if (m_data.trajectory_nsubstep > 0)
      m_data.vnode = DeparturePoints("vnode", nel, num_phys_lev, np, np, ndim);
    if (m_data.trajectory_nsubstep > 1)
      m_data.vdep  = DeparturePoints("vdep" , nel, num_phys_lev, np, np, ndim);
    homme::compose::set_views(
      g.m_spheremp,
      homme::compose::SetView<Real****>  (reinterpret_cast<Real*>(d.m_dp.data()),
                                          nel, np, np, nlev),
      homme::compose::SetView<Real*****> (
        reinterpret_cast<Real*>(independent_time_steps ?
                                d.m_divdp.data() :
                                s.m_dp3d.data()),
        nel, (independent_time_steps ? 1 : NUM_TIME_LEVELS), np, np, nlev),
      homme::compose::SetView<Real******>(
        reinterpret_cast<Real*>(t.qdp.data()),
        nel, t.qdp.extent_int(1), t.qdp.extent_int(2), np, np, nlev),
      homme::compose::SetView<Real*****> (reinterpret_cast<Real*>(t.Q.data()),
                                          nel, t.Q.extent_int(1), np, np, nlev),
      m_data.dep_pts, m_data.vnode, m_data.vdep, ndim);
  }

  m_data.independent_time_steps = independent_time_steps;
  if (m_data.nelemd == num_elems && m_data.qsize == params.qsize) return;

  m_data.qsize = params.qsize;
  Errors::runtime_check(m_data.qsize > 0,
                        "SL transport requires qsize > 0; if qsize == 0, use Eulerian.");
  m_data.nelemd = num_elems;

  Errors::runtime_check(m_data.hv_q >= 0 && m_data.hv_q <= m_data.qsize,
                        "semi_lagrange_hv_q should be in [0, qsize].");
  Errors::runtime_check(m_data.hv_subcycle_q >= 0,
                        "hypervis_subcycle_q should be >= 0.");

  m_tp_ne = Homme::get_default_team_policy<ExecSpace>(m_data.nelemd);
  m_tp_ne_qsize = Homme::get_default_team_policy<ExecSpace>(m_data.nelemd * m_data.qsize);
  m_tu_ne = TeamUtils<ExecSpace>(m_tp_ne);
  m_tu_ne_qsize = TeamUtils<ExecSpace>(m_tp_ne_qsize);
  if (m_data.nu_q > 0 && m_data.hv_q > 0) {
    m_tp_ne_hv_q = Homme::get_default_team_policy<ExecSpace>(m_data.nelemd * m_data.hv_q);
    m_tu_ne_hv_q = TeamUtils<ExecSpace>(m_tp_ne_hv_q);
  }

  m_sphere_ops.allocate_buffers(m_tu_ne_qsize);

  if (Context::singleton().get<Connectivity>().get_comm().root())
    printf("compose> nelemd %d qsize %d hv_q %d hv_subcycle_q %d lim %d "
           "independent_time_steps %d\n",
           m_data.nelemd, m_data.qsize, m_data.hv_q, m_data.hv_subcycle_q,
           m_data.limiter_option, (int) m_data.independent_time_steps);
}

int ComposeTransportImpl::requested_buffer_size () const {
  // FunctorsBuffersManager wants the size in terms of sizeof(Real).
  return (m_data.n_buf1*Buf1Alloc::shmem_size(nslot) +
          m_data.n_buf2*Buf2::shmem_size(nslot))/sizeof(Real);
}

void ComposeTransportImpl::init_buffers (const FunctorsBuffersManager& fbm) {
  Scalar* mem = reinterpret_cast<Scalar*>(fbm.get_memory());
  for (int i = 0; i < m_data.n_buf1; ++i) {
    m_data.buf1o[i] = Buf1o(mem, nslot);
    m_data.buf1e[i] = Buf1e(mem, nslot); // use the same memory
    mem += Buf1Alloc::shmem_size(nslot)/sizeof(Scalar);
  }
  for (int i = 0; i < m_data.n_buf2; ++i) {
    m_data.buf2[i] = Buf2(mem, nslot);
    mem += Buf2::shmem_size(nslot)/sizeof(Scalar);
  }
}

void ComposeTransportImpl::init_boundary_exchanges () {
  assert(m_data.qsize > 0); // after reset() called

  auto bm_exchange = Context::singleton().get<MpiBuffersManagerMap>()[MPI_EXCHANGE];
  const auto& sp = Context::singleton().get<SimulationParams>();

  // For qdp DSS at end of transport step.
  for (int i = 0; i < Q_NUM_TIME_LEVELS; ++i) {
    m_qdp_dss_be[i] = std::make_shared<BoundaryExchange>();
    auto be = m_qdp_dss_be[i];
    be->set_label(std::string("ComposeTransport-qdp-DSS-" + std::to_string(i)));
    be->set_diagnostics_level(sp.internal_diagnostics_level);
    be->set_buffers_manager(bm_exchange);
    be->set_num_fields(0, 0, m_data.qsize + 1);
    be->register_field(m_tracers.qdp, i, m_data.qsize, 0);
    be->register_field(m_derived.m_omega_p);
    be->registration_completed();
  }

  if (m_data.trajectory_nsubstep == 0) {
    for (int i = 0; i < 2; ++i) {
      m_v_dss_be[i] = std::make_shared<BoundaryExchange>();
      auto be = m_v_dss_be[i];
      be->set_label(std::string("ComposeTransport-v-DSS-" + std::to_string(i)));
      be->set_diagnostics_level(sp.internal_diagnostics_level);
      be->set_buffers_manager(bm_exchange);
      be->set_num_fields(0, 0, 2+i);
      be->register_field(m_derived.m_vstar, 2, 0);
      if (i) be->register_field(m_derived.m_divdp);
      be->registration_completed();
    }
  } else {
    for (int i = 0; i < 2; ++i) {
      m_v_dss_be[i] = std::make_shared<BoundaryExchange>();
      auto be = m_v_dss_be[i];
      be->set_label(std::string("ComposeTransport-v-DSS-" + std::to_string(i)));
      be->set_diagnostics_level(sp.internal_diagnostics_level);
      be->set_buffers_manager(bm_exchange);
      be->set_num_fields(0, 0, 3+i);
      be->register_field(m_tracers.qtens_biharmonic, 3+i, 0);
      be->registration_completed();
    }
  }

  // For optional HV applied to q.
  if (m_data.hv_q > 0 && m_data.nu_q > 0) {
    for (int i = 0; i < 2; ++i) {
      m_hv_dss_be[i] = std::make_shared<BoundaryExchange>();
      auto be = m_hv_dss_be[i];
      be->set_label(std::string("ComposeTransport-q-HV-" + std::to_string(i)));
      be->set_diagnostics_level(sp.internal_diagnostics_level);
      be->set_buffers_manager(bm_exchange);
      be->set_num_fields(0, 0, m_data.hv_q);
      if (i == 0) 
        be->register_field(m_tracers.qtens_biharmonic, m_data.hv_q, 0);
      else
        be->register_field(m_tracers.Q, m_data.hv_q, 0);
      be->registration_completed();
    }
  }
}

void ComposeTransportImpl::run (const TimeLevel& tl, const Real dt) {
  GPTLstart("compose_transport");

  if (m_data.trajectory_nsubstep == 0)
    calc_trajectory(tl.np1, dt);
  else
    calc_enhanced_trajectory(tl.np1, dt);
  
  GPTLstart("compose_isl");
  homme::compose::advect(tl.np1, tl.n0_qdp, tl.np1_qdp);
  Kokkos::fence();
  GPTLstop("compose_isl");
  
  if (m_data.hv_q > 0 && m_data.nu_q > 0) {
    GPTLstart("compose_hypervis_scalar");
    advance_hypervis_scalar(dt);
    Kokkos::fence();
    GPTLstop("compose_hypervis_scalar");
  }
  
  GPTLstart("compose_cedr_global");
  homme::compose::set_dp3d_np1(m_data.independent_time_steps ?
                               0 : // dp3d is actually divdp
                               tl.np1);
  const auto run_cedr = homme::compose::property_preserve_global();
  if (run_cedr) Kokkos::fence();
  GPTLstop("compose_cedr_global");
  GPTLstart("compose_cedr_local");
  if (run_cedr) {
    homme::compose::property_preserve_local(m_data.limiter_option);
    Kokkos::fence();
  }
  GPTLstop("compose_cedr_local");    

  const auto np1 = tl.np1;
  const auto np1_qdp = tl.np1_qdp;
  const auto qsize = m_data.qsize;
  
  if ( ! run_cedr) {
    // For analysis purposes, property preservation was not run. Need to convert
    // Q to qdp.
    const auto qdp = m_tracers.qdp;
    const auto Q = m_tracers.Q;
    const auto dp3d = m_state.m_dp3d;
    const auto spheremp = m_geometry.m_spheremp;
    const auto f = KOKKOS_LAMBDA (const int idx) {
      int ie, q, i, j, lev;
      idx_ie_q_ij_nlev<num_lev_pack>(qsize, idx, ie, q, i, j, lev);
      qdp(ie,np1_qdp,q,i,j,lev) = Q(ie,q,i,j,lev)/dp3d(ie,np1,i,j,lev);
    };
    launch_ie_q_ij_nlev<num_lev_pack>(qsize, f);
  }
  
  { // DSS qdp and omega
    GPTLstart("compose_dss_q");
    const auto qdp = m_tracers.qdp;
    const auto spheremp = m_geometry.m_spheremp;
    const auto f1 = KOKKOS_LAMBDA (const int idx) {
      int ie, q, i, j, lev;
      idx_ie_q_ij_nlev<num_lev_pack>(qsize, idx, ie, q, i, j, lev);
      qdp(ie,np1_qdp,q,i,j,lev) *= spheremp(ie,i,j);
    };
    launch_ie_q_ij_nlev<num_lev_pack>(qsize, f1);
    const auto omega = m_derived.m_omega_p;
    const auto f2 = KOKKOS_LAMBDA (const int idx) {
      int ie, i, j, lev;
      idx_ie_ij_nlev<num_lev_pack>(idx, ie, i, j, lev);
      omega(ie,i,j,lev) *= spheremp(ie,i,j);
    };
    launch_ie_ij_nlev<num_lev_pack>(f2);
    m_qdp_dss_be[tl.np1_qdp]->exchange(m_geometry.m_rspheremp);
    Kokkos::fence();
    GPTLstop("compose_dss_q");
  }
  
  if (m_data.cdr_check) {
    GPTLstart("compose_cedr_check");
    homme::compose::property_preserve_check();
    Kokkos::fence();
    GPTLstop("compose_cedr_check");
  }
  
  GPTLstop("compose_transport");
}

} // namespace Homme

#endif // HOMME_ENABLE_COMPOSE
