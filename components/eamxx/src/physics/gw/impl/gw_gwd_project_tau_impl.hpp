#ifndef GW_GWD_PROJECT_TAU_IMPL_HPP
#define GW_GWD_PROJECT_TAU_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU

namespace scream {
namespace gw {

/*
 * Implementation of gw gwd_project_tau. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gwd_project_tau(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const GwCommonInit& init,
  const Int& pver,
  const Int& pgwv,
  const Int& tend_level,
  const uview_2d<const Real>& tau,
  const uview_1d<const Real>& ubi,
  const uview_1d<const Real>& c,
  const Real& xv,
  const Real& yv,
  // Outputs
  const uview_2d<Real>& taucd)
{
  uview_1d<Real> taub, tauf;
  workspace.template take_many_contiguous_unsafe<2>(
    {"taub", "tauf"},
    {&taub, &tauf});

  const Real ubi_tend = ubi(tend_level+1);

  const int num_pgwv = 2*pgwv + 1;

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, taub.size()), [&] (const int k) {
      taub(k) = 0;
      tauf(k) = 0;
    }
  );

  team.team_barrier();

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, init.ktop+1, tend_level + 2), [&] (const int k) {

      for (int l = 0; l < num_pgwv; ++l) {
        const Real tausg = Kokkos::copysign(tau(l,k), c(l)-ubi(k));
        if (c(l) < ubi_tend) {
          taub(k) += tausg;
        }
        else if (c(l) > ubi_tend) {
          tauf(k) += tausg;
        }
      }
    });

  team.team_barrier();

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, init.ktop+1, tend_level + 2), [&] (const int k) {
      if (xv > 0) {
        taucd(k, GWC::east) = tauf(k) * xv;
        taucd(k, GWC::west) = taub(k) * xv;
      }
      else if (xv < 0) {
        taucd(k, GWC::east) = taub(k) * xv;
        taucd(k, GWC::west) = tauf(k) * xv;
      }

      if (yv > 0) {
        taucd(k, GWC::north) = tauf(k) * yv;
        taucd(k, GWC::south) = taub(k) * yv;
      }
      else if (yv < 0) {
        taucd(k, GWC::north) = taub(k) * yv;
        taucd(k, GWC::south) = tauf(k) * yv;
      }
    });

  workspace.template release_many_contiguous<2>(
    {&taub, &tauf});
}

} // namespace gw
} // namespace scream

#endif
