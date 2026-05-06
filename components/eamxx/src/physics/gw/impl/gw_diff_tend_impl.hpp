#ifndef GW_GW_DIFF_TEND_IMPL_HPP
#define GW_GW_DIFF_TEND_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gw_diff_tend. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gw_diff_tend(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const Int& pver,
  const Int& kbot,
  const Int& ktop,
  const uview_1d<const Real>& q,
  const Real& dt,
  const uview_1d<const Real>& decomp_ca,
  const uview_1d<const Real>& decomp_cc,
  const uview_1d<const Real>& decomp_dnom,
  const uview_1d<const Real>& decomp_ze,
  // Outputs
  const uview_1d<Real>& dq)
{
  auto qnew = workspace.take("qnew");

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, pver), [&] (const int k) {
      qnew(k) = q(k);
    });

  // Solve the diffusion matrix.
  vd_lu_solve(team, workspace, pver, decomp_ca, decomp_cc, decomp_dnom, decomp_ze, ktop+1, kbot+1, 0, qnew);

  // Evaluate tendency to be reported back.
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, pver), [&] (const int k) {
      dq(k) = (qnew(k)-q(k)) / dt;
    });

  workspace.release(qnew);
}

} // namespace gw
} // namespace scream

#endif
