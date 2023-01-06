/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#include "ComposeTransportImpl.hpp"

namespace Homme {

void ComposeTransportImpl::advance_hypervis_scalar (const Real dt_q) {
  const auto dt = dt_q / m_data.hv_subcycle_q;
  const auto hv_q = m_data.hv_q;
  const auto nu_q = m_data.nu_q;
  const auto Qtens = m_tracers.qtens_biharmonic;
  const auto Q = m_tracers.Q;
  const auto spheremp = m_geometry.m_spheremp;
  const auto tu_ne_hv_q = m_tu_ne_hv_q;
  const auto sphere_ops = m_sphere_ops;
  for (int it = 0; it < m_data.hv_subcycle_q; ++it) {
    { // Qtens = Q
      const auto f = KOKKOS_LAMBDA (const int idx) {
        int ie, q, i, j, lev;
        idx_ie_q_ij_nlev<num_lev_pack>(hv_q, idx, ie, q, i, j, lev);
        Qtens(ie,q,i,j,lev) = Q(ie,q,i,j,lev);
      };
      Kokkos::fence();
      launch_ie_q_ij_nlev<num_lev_pack>(hv_q, f);
    }
    // biharmonic_wk_scalar
    const auto laplace_simple_Qtens = [&] () {
      const auto f = KOKKOS_LAMBDA (const MT& team) {
        KernelVariables kv(team, hv_q, tu_ne_hv_q);
        const auto Qtens_ie = Homme::subview(Qtens, kv.ie, kv.iq);
        sphere_ops.laplace_simple(kv, Qtens_ie, Qtens_ie);
      };
      Kokkos::fence();
      Kokkos::parallel_for(m_tp_ne_hv_q, f);
    };
    laplace_simple_Qtens();
    m_hv_dss_be[0]->exchange(m_geometry.m_rspheremp);
    if (m_data.hv_scaling == 0) {
      Kokkos::fence();
      laplace_simple_Qtens();
    } else {
      const auto tensorvisc = m_geometry.m_tensorvisc;
      const auto f = KOKKOS_LAMBDA (const MT& team) {
        KernelVariables kv(team, hv_q, tu_ne_hv_q);
        const auto Qtens_ie = Homme::subview(Qtens, kv.ie, kv.iq);
        sphere_ops.laplace_tensor(kv, Homme::subview(tensorvisc, kv.ie),
                                  Qtens_ie, Qtens_ie);
      };
      Kokkos::fence();
      Kokkos::parallel_for(m_tp_ne_hv_q, f);
    }
    { // Compute Q = Q spheremp - dt nu_q Qtens. N.B. spheremp is already in
      // Qtens from divergence_sphere_wk.
      const auto f = KOKKOS_LAMBDA (const int idx) {
        int ie, q, i, j, lev;
        idx_ie_q_ij_nlev<num_lev_pack>(hv_q, idx, ie, q, i, j, lev);
        Q(ie,q,i,j,lev) = (Q(ie,q,i,j,lev) * spheremp(ie,i,j)
                           - dt * nu_q * Qtens(ie,q,i,j,lev));
      };
      Kokkos::fence();
      launch_ie_q_ij_nlev<num_lev_pack>(hv_q, f);
    }
    // Halo exchange Q and apply rspheremp.
    Kokkos::fence();
    m_hv_dss_be[1]->exchange(m_geometry.m_rspheremp);
  }
}

} // namespace Homme

#endif // HOMME_ENABLE_COMPOSE
