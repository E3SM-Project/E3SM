#ifndef GW_GWD_COMPUTE_TENDENCIES_FROM_STRESS_DIVERGENCE_IMPL_HPP
#define GW_GWD_COMPUTE_TENDENCIES_FROM_STRESS_DIVERGENCE_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU
#include "util/eamxx_utils.hpp"
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
  const Int& ngwv,
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
  // Temporary ubar tendencies (overall, at wave l, and at saturation).
  uview_1d<Real> ubt = workspace.take("ubt");

  const Real ptaper = do_taper ? std::cos(lat) : 1.;

  // Force tau at the top of the model to zero, if requested.
  if (s_common_init.tau_0_ubc) {
    for (size_t i = 0; i < tau.extent(0); ++i) {
      tau(i,0) = 0.;
    }
  }

  // Loop over levels from top to bottom
  Kokkos::parallel_for(
    Kokkos::TeamVectorRange(team, s_common_init.ktop, max_level+1), [&] (const int k) {
    //  Accumulate the mean wind tendency over wavenumber.
    ubt(k) = 0.;

    for (int l = -ngwv; l <= ngwv; ++l) { // loop over wave
      int nl_idx = l + ngwv; // 0-based idx for -ngwv:ngwv arrays
      int pl_idx = nl_idx + (pgwv - ngwv); // 0-based idx -pgwv:pgwv arrays

      // Determine the wind tendency, including excess stress carried down
      // from above.
      Real ubtl = C::gravit * (tau(pl_idx,k+1)-tau(pl_idx,k)) * rdpm(k);

      if (s_common_init.orographic_only) {
        // Require that the tendency be no larger than the analytic
        // solution for a saturated region [proportional to (u-c)^3].
        Real temp = c(pl_idx)-ubm(k);
        temp = temp * temp * temp; // BFB with fortran **3
        Real ubtlsat = s_common_init.effkwv * std::abs(temp) / (2*GWC::rog*t(k)*nm(k));
        ubtl = std::min(ubtl, ubtlsat);
      }

      // Apply tendency limits to maintain numerical stability.
      // 1. du/dt < |c-u|/dt  so u-c cannot change sign
      //    (u^n+1 = u^n + du/dt * dt)
      // 2. du/dt < tndmax    so that ridicuously large tendencies are not
      //    permitted
      ubtl = std::min(ubtl, GWC::umcfac * std::abs(c(pl_idx)-ubm(k)) / dt);
      ubtl = std::min(ubtl, s_common_init.tndmax);

      if (k <= tend_level) {
        // Save tendency for each wave (for later computation of kzz),
        // applying efficiency and taper:
        gwut(k,nl_idx) = sign(ubtl, c(pl_idx)-ubm(k)) * effgw * ptaper;

        if (!s_common_init.orographic_only) {
          ubt(k) = ubt(k) + gwut(k,nl_idx);
        }
        else {
          ubt(k) = ubt(k) + sign(ubtl, c(pl_idx)-ubm(k));
        }

        // Redetermine the effective stress on the interface below from
        // the wind tendency. If the wind tendency was limited above,
        // then the new stress will be smaller than the old stress,
        // causing stress divergence in the next layer down. This
        // smoothes large stress divergences downward while conserving
        // total stress.
        tau(pl_idx,k+1) = tau(pl_idx,k) + ubtl * dpm(k) / C::gravit;
      }
    }

    // Project the mean wind tendency onto the components.
    if (k <= tend_level) {
      if (!s_common_init.orographic_only) {
        utgw(k) = ubt(k) * xv;
        vtgw(k) = ubt(k) * yv;
      }
      else {
        utgw(k) = ubt(k) * xv * effgw * ptaper;
        vtgw(k) = ubt(k) * yv * effgw * ptaper;
      }
    }
  });

  workspace.release(ubt);
}

} // namespace gw
} // namespace scream

#endif
