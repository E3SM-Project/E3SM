/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#include "ComposeTransportImplEnhancedTrajectoryImpl.hpp"

namespace Homme {

// For limit_etam.
void ComposeTransportImpl::setup_enhanced_trajectory () {
  const auto etai = cmvdc(m_hvcoord.etai);
  const Real deta_ave = (etai(num_phys_lev) - etai(0)) / num_phys_lev;
  m_data.deta_tol = 10*std::numeric_limits<Real>::epsilon()*deta_ave;

  // diff(etai)
  m_data.hydetai = decltype(m_data.hydetai)("hydetai");
  const auto detai_pack = Kokkos::create_mirror_view(m_data.hydetai);
  HostViewUnmanaged<Real[NUM_PHYSICAL_LEV]> detai(pack2real(detai_pack));
  for (int k = 0; k < num_phys_lev; ++k)
    detai(k) = etai(k+1) - etai(k);
  Kokkos::deep_copy(m_data.hydetai, detai_pack);

  const auto etamp = cmvdc(m_hvcoord.etam);
  HostViewUnmanaged<Real[NUM_PHYSICAL_LEV]> etam(pack2real(etamp));
  
  // hydetam_ref.
  m_data.hydetam_ref = decltype(m_data.hydetam_ref)("hydetam_ref");
  const auto m = Kokkos::create_mirror_view(m_data.hydetam_ref);
  const int nlev = num_phys_lev;
  m(0) = etam(0) - etai(0);
  for (int k = 1; k < nlev; ++k) m(k) = etam(k) - etam(k-1);
  m(nlev) = etai(nlev) - etam(nlev-1);
  Kokkos::deep_copy(m_data.hydetam_ref, m);

  // etam
  homme::compose::set_hvcoord(etai(0), etai(num_phys_lev), etam.data());
}

void ComposeTransportImpl::calc_enhanced_trajectory (const int np1, const Real dt) {
  GPTLstart("compose_calc_enhanced_trajectory");

  const auto& dep_pts = m_data.dep_pts;
  const auto& vnode = m_data.vnode;
  const auto& vdep = m_data.vdep;

  init_dep_points(*this, dep_pts);

  const int nelemd = m_data.nelemd;
  const Real dtsub = dt / m_data.trajectory_nsubstep;
  const int nsubstep = m_data.trajectory_nsubstep;
  for (int step = 0; step < nsubstep; ++step) {
    {
      Kokkos::fence();
      GPTLstart("compose_vnode");
      const Real alpha[] = {Real(nsubstep-step-1)/nsubstep,
                            Real(nsubstep-step  )/nsubstep};
      const CVSlot v1(m_derived.m_vstar.data(), nelemd, 1);
      const CDpSlot dp1(m_derived.m_dp.data(), nelemd, 1);
      const auto& v2 = m_state.m_v;
      const auto& dp2 = m_state.m_dp3d;
      calc_nodal_velocities(*this, dtsub, alpha,
                            v1, dp1, 0, v2, dp2, np1,
                            vnode);
      Kokkos::fence();
      GPTLstop("compose_vnode");
    }

    GPTLstart("compose_v_bexchv");
    dss_vnode(*this, vnode);
    Kokkos::fence();
    GPTLstop("compose_v_bexchv");

    if (step == 0) {
      update_dep_points(*this, dtsub, vnode, dep_pts);
    } else {
      GPTLstart("compose_vdep");
      homme::compose::calc_v_departure(step, dtsub);
      Kokkos::fence();
      GPTLstop("compose_vdep");

      update_dep_points(*this, dtsub, vdep, dep_pts);
    }
  }
  Kokkos::fence();

  if (m_data.independent_time_steps) {
    GPTLstart("compose_floating_dep_pts");
    interp_departure_points_to_floating_level_midpoints(*this, np1);
    Kokkos::fence();
    GPTLstop("compose_floating_dep_pts");
  }

  GPTLstop("compose_calc_enhanced_trajectory");
}

} // namespace Homme

#endif // HOMME_ENABLE_COMPOSE
