/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#include "Config.hpp"
#ifdef HOMME_ENABLE_COMPOSE

#include "ComposeTransportImpl.hpp"
#include "Context.hpp"
#include "VerticalRemapManager.hpp"
#include "RemapFunctor.hpp"

namespace Homme {
using cti = ComposeTransportImpl;

void ComposeTransportImpl
::remap_v (const ExecViewUnmanaged<const Scalar*[NUM_TIME_LEVELS][NP][NP][NUM_LEV]>& dp3d,
           const int np1, const ExecViewUnmanaged<const Scalar*[NP][NP][NUM_LEV]>& dp,
           const ExecViewUnmanaged<Scalar*[2][NP][NP][NUM_LEV]>& v) {
  using Kokkos::parallel_for;
  const auto& vrm = Context::singleton().get<VerticalRemapManager>();
  const auto r = vrm.get_remapper();
  const auto policy = Kokkos::RangePolicy<ExecSpace>(0, dp3d.extent_int(0)*NP*NP*NUM_LEV*2);
  const auto pre = KOKKOS_LAMBDA (const int idx) {
    int ie, q, i, j, k;
    cti::idx_ie_q_ij_nlev<NUM_LEV>(2, idx, ie, q, i, j, k);
    v(ie,q,i,j,k) *= dp3d(ie,np1,i,j,k);
  };
  parallel_for(policy, pre);
  Kokkos::fence();
  r->remap1(dp3d, np1, dp, v, 2);
  Kokkos::fence();
  const auto post = KOKKOS_LAMBDA (const int idx) {
    int ie, q, i, j, k;
    cti::idx_ie_q_ij_nlev<NUM_LEV>(2, idx, ie, q, i, j, k);
    v(ie,q,i,j,k) /= dp(ie,i,j,k);
  };
  parallel_for(policy, post);
}

void ComposeTransportImpl::remap_q (const TimeLevel& tl) {
  GPTLstart("compose_vertical_remap");
  const auto np1 = tl.np1;
  const auto np1_qdp = tl.np1_qdp;
  const auto dp = m_derived.m_divdp;
  const auto dp3d = m_state.m_dp3d;
  const auto qdp = m_tracers.qdp;
  const auto q = m_tracers.Q;
  const int nq = m_tracers.num_tracers();
  const auto& vrm = Context::singleton().get<VerticalRemapManager>();
  const auto r = vrm.get_remapper();
  r->remap1(dp, dp3d, np1, qdp, np1_qdp, nq);
  const auto post = KOKKOS_LAMBDA (const int idx) {
    int ie, iq, i, j, k;
    cti::idx_ie_q_ij_nlev<NUM_LEV>(nq, idx, ie, iq, i, j, k);
    q(ie,iq,i,j,k) = qdp(ie,np1_qdp,iq,i,j,k)/dp3d(ie,np1,i,j,k);
  };
  const auto policy = Kokkos::RangePolicy<ExecSpace>(0, dp3d.extent_int(0)*NP*NP*NUM_LEV*nq);
  Kokkos::fence();
  Kokkos::parallel_for(policy, post);
  Kokkos::fence();
  GPTLstop("compose_vertical_remap");
}

} // namespace Homme

#endif // HOMME_ENABLE_COMPOSE
