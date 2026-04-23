#ifndef GW_GWD_COMPUTE_TENDENCIES_FROM_STRESS_DIVERGENCE_IMPL_HPP
#define GW_GWD_COMPUTE_TENDENCIES_FROM_STRESS_DIVERGENCE_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU
#include "share/util/eamxx_utils.hpp"

#include <ekat_math_utils.hpp>
#include <cmath>

namespace scream {
namespace gw {

/*
 * Implementation of gw gwd_compute_tendencies_from_stress_divergence. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gwd_compute_tendencies_from_stress_divergence(
  // Inputs
  const MemberType& team,
  const Workspace& workspace,
  const GwCommonInit& init,
  const Int& pver,
  const Int& pgwv,
  const bool& do_taper,
  const Real& dt,
  const Real& effgw,
  const Int& tend_level,
  const Int& max_level,
  const Real& lat,
  const uview_1d<const Real>& dpm,
  const uview_1d<const Real>& rdpm,
  const uview_1d<const Real>& c,
  const uview_1d<const Real>& ubm,
  const uview_1d<const Real>& t,
  const uview_1d<const Real>& nm,
  const Real& xv,
  const Real& yv,
  // Inputs/Outputs
  const uview_2d<Real>& tau,
  // Outputs
  const uview_2d<Real>& gwut,
  const uview_1d<Real>& utgw,
  const uview_1d<Real>& vtgw)
{
  const Real ptaper = do_taper ? std::cos(lat) : 1;

  // Get temporary workspaces and change them to desired dimensions
  auto work_1d = workspace.take("work_1d");
  uview_2d<Real> work(work_1d.data(), pver, 2*pgwv + 1);

  // Loop waves
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, -pgwv, pgwv+1), [&] (const int l) {
    //  Accumulate the mean wind tendency over wavenumber.
    const int pl_idx = l + pgwv; // 0-based idx for -pgwv:pgwv arrays

    // Force tau at the top of the model to zero, if requested.
    if (init.tau_0_ubc) {
      tau(pl_idx,0) = 0;
    }

    // Loop over levels from top to bottom. Each level reads and writes to
    // the next level, so this loop must be serialized.
    for (int k = init.ktop+1; k <= max_level; ++k) {

      // Determine the wind tendency, including excess stress carried down
      // from above.
      Real ubtl = C::gravit.value * (tau(pl_idx,k+1)-tau(pl_idx,k)) * rdpm(k);

      // if (Kokkos::isnan(tau(pl_idx,k  ))) { printf("[gwd_compute_tend 0] tau k+0 NaN @ k=%d  pl_idx=%d\n", k, pl_idx); }
      // if (Kokkos::isinf(tau(pl_idx,k  ))) { printf("[gwd_compute_tend 0] tau k+0 Inf @ k=%d  pl_idx=%d\n", k, pl_idx); }

      // if (Kokkos::isnan(tau(pl_idx,k+1))) { printf("[gwd_compute_tend 0] tau k+1 NaN @ k=%d  pl_idx=%d\n", k, pl_idx); }
      // if (Kokkos::isinf(tau(pl_idx,k+1))) { printf("[gwd_compute_tend 0] tau k+1 Inf @ k=%d  pl_idx=%d\n", k, pl_idx); }

      // if (Kokkos::isnan(rdpm(k))) { printf("[gwd_compute_tend 0] rdpm NaN @ k=%d\n", k); }
      // if (Kokkos::isinf(rdpm(k))) { printf("[gwd_compute_tend 0] rdpm Inf @ k=%d\n", k); }

      // if (Kokkos::isnan(ubtl)) { printf("[gwd_compute_tend 0] ubtl NaN @ k=%d\n", k); }
      // if (Kokkos::isinf(ubtl)) { printf("[gwd_compute_tend 0] ubtl Inf @ k=%d\n", k); }

      // Real tau_diff = tau(pl_idx,k+1)-tau(pl_idx,k);
      // if (Kokkos::isnan(ubtl)) { printf("[gct 0] ubtl NaN @ [%d,%d] => %f - %f = %f\n", k, pl_idx, tau(pl_idx,k+1), tau(pl_idx,k), tau_diff);}
      // if (Kokkos::isinf(ubtl)) { printf("[gct 0] ubtl Inf @ [%d,%d] => %f - %f = %f\n", k, pl_idx, tau(pl_idx,k+1), tau(pl_idx,k), tau_diff);}

      if (init.orographic_only) {
        // Require that the tendency be no larger than the analytic
        // solution for a saturated region [proportional to (u-c)^3].
        Real ubtlsat = init.effkwv * std::abs(bfb_cube(c(pl_idx)-ubm(k))) / (2*GWC::rog*t(k)*nm(k));
        ubtl = ekat::impl::min(ubtl, ubtlsat);
      }

      // if (Kokkos::isnan(t(k))) { printf("[gwd_compute_tend 1] t NaN @ k=%d\n", k); }
      // if (Kokkos::isinf(t(k))) { printf("[gwd_compute_tend 1] t Inf @ k=%d\n", k); }

      // if (Kokkos::isnan(nm(k))) { printf("[gwd_compute_tend 1] nm NaN @ k=%d\n", k); }
      // if (Kokkos::isinf(nm(k))) { printf("[gwd_compute_tend 1] nm Inf @ k=%d\n", k); }

      // if (Kokkos::isnan(ubtl)) { printf("[gwd_compute_tend 1] ubtl NaN @ k=%d\n", k); }
      // if (Kokkos::isinf(ubtl)) { printf("[gwd_compute_tend 1] ubtl Inf @ k=%d\n", k); }

      // Apply tendency limits to maintain numerical stability.
      // 1. du/dt < |c-u|/dt  so u-c cannot change sign
      //    (u^n+1 = u^n + du/dt * dt)
      // 2. du/dt < tndmax    so that ridicuously large tendencies are not
      //    permitted
      ubtl = ekat::impl::min(ubtl, GWC::umcfac * std::abs(c(pl_idx)-ubm(k)) / dt);
      ubtl = ekat::impl::min(ubtl, init.tndmax);

      // if (Kokkos::isnan(gwut(k,pl_idx))) { printf("[gwd_compute_tend] (outside) gwut(k,pl_idx) NaN @ k=%d  pl_idx=%d\n", k, pl_idx); }
      // if (Kokkos::isinf(gwut(k,pl_idx))) { printf("[gwd_compute_tend] (outside) gwut(k,pl_idx) Inf @ k=%d  pl_idx=%d\n", k, pl_idx); }

      if (k <= tend_level) {

        // if (Kokkos::isnan(effgw)) { printf("[gwd_compute_tend] effgw NaN @ k=%d\n", k); }
        // if (Kokkos::isinf(effgw)) { printf("[gwd_compute_tend] effgw Inf @ k=%d\n", k); }

        // if (Kokkos::isnan(ptaper)) { printf("[gwd_compute_tend] ptaper NaN @ k=%d\n", k); }
        // if (Kokkos::isinf(ptaper)) { printf("[gwd_compute_tend] ptaper Inf @ k=%d\n", k); }

        // if (Kokkos::isnan(ubtl)) { printf("[gwd_compute_tend 2] ubtl NaN @ k=%d\n", k); }
        // if (Kokkos::isinf(ubtl)) { printf("[gwd_compute_tend 2] ubtl Inf @ k=%d\n", k); }

        // if (Kokkos::isnan(effgw)) { printf("[gwd_compute_tend] effgw NaN @ k=%d\n", k); }
        // if (Kokkos::isinf(effgw)) { printf("[gwd_compute_tend] effgw Inf @ k=%d\n", k); }

        // if (Kokkos::isnan(c(pl_idx))) { printf("[gwd_compute_tend] c(pl_idx) NaN @ k=%d  pl_idx=%d\n", k, pl_idx); }
        // if (Kokkos::isinf(c(pl_idx))) { printf("[gwd_compute_tend] c(pl_idx) Inf @ k=%d  pl_idx=%d\n", k, pl_idx); }

        // if (Kokkos::isnan(ubm(k))) { printf("[gwd_compute_tend] ubm NaN @ k=%d\n", k); }
        // if (Kokkos::isinf(ubm(k))) { printf("[gwd_compute_tend] ubm Inf @ k=%d\n", k); }

        // if (Kokkos::isnan(gwut(k,pl_idx))) { printf("[gwd_compute_tend] (before) gwut(k,pl_idx) NaN @ k=%d  pl_idx=%d\n", k, pl_idx); }
        // if (Kokkos::isinf(gwut(k,pl_idx))) { printf("[gwd_compute_tend] (before) gwut(k,pl_idx) Inf @ k=%d  pl_idx=%d\n", k, pl_idx); }

        // Save tendency for each wave (for later computation of kzz),
        // applying efficiency and taper:
        gwut(k,pl_idx) = Kokkos::copysign(ubtl, c(pl_idx)-ubm(k)) * effgw * ptaper;

        // if (Kokkos::isnan(gwut(k,pl_idx))) { printf("[gwd_compute_tend] (after) gwut(k,pl_idx) NaN @ k=%d  pl_idx=%d\n", k, pl_idx); }
        // if (Kokkos::isinf(gwut(k,pl_idx))) { printf("[gwd_compute_tend] (after) gwut(k,pl_idx) Inf @ k=%d  pl_idx=%d\n", k, pl_idx); }

        // atomic_sum for a workspace item ubt(k) are another option here. It works
        // but, since the order of operations is non-deterministic, there are
        // non-deterministic round-off differences from run to run.
        if (!init.orographic_only) {
          work(k, pl_idx) = gwut(k,pl_idx);
        }
        else {
          work(k, pl_idx) = Kokkos::copysign(ubtl, c(pl_idx)-ubm(k));
        }

        // Redetermine the effective stress on the interface below from
        // the wind tendency. If the wind tendency was limited above,
        // then the new stress will be smaller than the old stress,
        // causing stress divergence in the next layer down. This
        // smoothes large stress divergences downward while conserving
        // total stress.
        tau(pl_idx,k+1) = tau(pl_idx,k) + ubtl * dpm(k) / C::gravit.value;
      }
    }
  });

  team.team_barrier();

  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, init.ktop+1, tend_level+1), [&] (const int k) {
    // Serialize the sum so it's repeatable
    Real ubt = 0;
    for (size_t i = 0; i < work.extent(1); ++i) {
      ubt += work(k, i);
    }

    // Project the mean wind tendency onto the components.
    if (!init.orographic_only) {
      utgw(k) = ubt * xv;
      vtgw(k) = ubt * yv;
    }
    else {
      utgw(k) = ubt * xv * effgw * ptaper;
      vtgw(k) = ubt * yv * effgw * ptaper;
    }
  });

  team.team_barrier();

  // Release temporary variables from the workspace
  workspace.release(work_1d);
}

} // namespace gw
} // namespace scream

#endif
