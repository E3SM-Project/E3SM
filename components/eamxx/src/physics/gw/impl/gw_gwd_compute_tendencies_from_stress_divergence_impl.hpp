#ifndef GW_GWD_COMPUTE_TENDENCIES_FROM_STRESS_DIVERGENCE_IMPL_HPP
#define GW_GWD_COMPUTE_TENDENCIES_FROM_STRESS_DIVERGENCE_IMPL_HPP

#include "gw_functions.hpp" // for ETI only but harmless for GPU
#include <cmath>

namespace scream {
namespace gw {

/*
 * Implementation of gw gwd_compute_tendencies_from_stress_divergence. Clients should NOT
 * #include this file, but include gw_functions.hpp instead.
 */

// TODO: Move to common area
template <typename T>
T sign(T a, T b) {
    return (b >= 0) ? std::abs(a) : -std::abs(a);
}

template<typename S, typename D>
KOKKOS_FUNCTION
void Functions<S,D>::gwd_compute_tendencies_from_stress_divergence(
  // Inputs
  const Int& ncol,
  const Int& pver,
  const Int& pgwv,
  const Int& ngwv,
  const bool& do_taper,
  const Real& dt,
  const Real& effgw,
  const uview_1d<const Int>& tend_level,
  const uview_1d<const Real>& lat,
  const uview_2d<const Real>& dpm,
  const uview_2d<const Real>& rdpm,
  const uview_2d<const Real>& c,
  const uview_2d<const Real>& ubm,
  const uview_2d<const Real>& t,
  const uview_2d<const Real>& nm,
  const uview_1d<const Real>& xv,
  const uview_1d<const Real>& yv,
  // Inputs/Outputs
  const uview_3d<Real>& tau,
  // Outputs
  const uview_3d<Real>& gwut,
  const uview_2d<Real>& utgw,
  const uview_2d<Real>& vtgw)
{
  // Temporary ubar tendencies (overall, at wave l, and at saturation).
  // These should probably come from workspace manager
  view_2d<Real> ubt("ubt", ncol, s_common_init.pver);
  view_1d<Real> ubtl("ubtl", ncol), ubtlsat("ubtlsat", ncol);

  // Polar taper.
  view_1d<Real> ptaper("ptaper", ncol);

  if (do_taper) { // taper CM only
    for (int i = 0; i < ncol; ++i) {
      ptaper(i) = std::cos(lat(i));
    }
  }
  else { //  do not taper other cases
    Kokkos::deep_copy(ptaper, 1.);
  }

  // Force tau at the top of the model to zero, if requested.
  if (s_common_init.tau_0_ubc) {
    for (size_t i = 0; i < tau.extent(0); ++i) {
      for (size_t j = 0; j < tau.extent(1); ++j) {
        tau(i,j,0) = 0.;
      }
    }
  }

  // Find max
  int maxval = -999999;
  for (int i = 0; i < ncol; ++i) {
    if (tend_level(i) > maxval) maxval = tend_level(i);
  }

  // Loop over levels from top to bottom
  for (int k = s_common_init.ktop; k <= maxval; ++k) {
    //  Accumulate the mean wind tendency over wavenumber.
    for (int i = 0; i < ncol; ++i) {
      ubt(i,k) = 0.;
    }

    for (int l = -ngwv; l <= ngwv; ++l) { // loop over wave
      int nl_idx = l + ngwv; // 0-based idx for -ngwv:ngwv arrays
      int pl_idx = nl_idx + (pgwv - ngwv); // 0-based idx -pgwv:pgwv arrays

      // Determine the wind tendency, including excess stress carried down
      // from above.
      for (int i = 0; i < ncol; ++i) {
        ubtl(i) = C::gravit * (tau(i,pl_idx,k+1)-tau(i,pl_idx,k)) * rdpm(i,k);
      }

      if (s_common_init.orographic_only) {
        // Require that the tendency be no larger than the analytic
        // solution for a saturated region [proportional to (u-c)^3].
        for (int i = 0; i < ncol; ++i) {
          Real temp = c(i,pl_idx)-ubm(i,k);
          temp = temp * temp * temp; // BFB with fortran **3
          ubtlsat(i) = s_common_init.effkwv * std::abs(temp) / (2*GWC::rog*t(i,k)*nm(i,k));
          ubtl(i) = std::min(ubtl(i), ubtlsat(i));
        }
      }

      // Apply tendency limits to maintain numerical stability.
      // 1. du/dt < |c-u|/dt  so u-c cannot change sign
      //    (u^n+1 = u^n + du/dt * dt)
      // 2. du/dt < tndmax    so that ridicuously large tendencies are not
      //    permitted
      for (int i = 0; i < ncol; ++i) {
        ubtl(i) = std::min(ubtl(i), GWC::umcfac * std::abs(c(i,pl_idx)-ubm(i,k)) / dt);
        ubtl(i) = std::min(ubtl(i), s_common_init.tndmax);
      }

      for (int i = 0; i < ncol; ++i) {
        if (k <= tend_level(i)) {
          // Save tendency for each wave (for later computation of kzz),
          // applying efficiency and taper:
          gwut(i,k,nl_idx) = sign(ubtl(i), c(i,pl_idx)-ubm(i,k)) * effgw * ptaper(i);
        }
      }

      if (!s_common_init.orographic_only) {
        for (int i = 0; i < ncol; ++i) {
          if (k <= tend_level(i)) {
            ubt(i,k) = ubt(i,k) + gwut(i,k,nl_idx);
          }
        }
      }
      else {
        for (int i = 0; i < ncol; ++i) {
          if (k <= tend_level(i)) {
            ubt(i,k) = ubt(i,k) + sign(ubtl(i), c(i,pl_idx)-ubm(i,k));
          }
        }
      }

      for (int i = 0; i < ncol; ++i) {
        if (k <= tend_level(i)) {
          // Redetermine the effective stress on the interface below from
          // the wind tendency. If the wind tendency was limited above,
          // then the new stress will be smaller than the old stress,
          // causing stress divergence in the next layer down. This
          // smoothes large stress divergences downward while conserving
          // total stress.
          tau(i,pl_idx,k+1) = tau(i,pl_idx,k) + ubtl(i) * dpm(i,k) / C::gravit;
        }
      }
    }

    // Project the mean wind tendency onto the components.
    if (!s_common_init.orographic_only) {
      for (int i = 0; i < ncol; ++i) {
        if (k <= tend_level(i)) {
          utgw(i,k) = ubt(i,k) * xv(i);
          vtgw(i,k) = ubt(i,k) * yv(i);
        }
      }
    }
    else {
      for (int i = 0; i < ncol; ++i) {
        if (k <= tend_level(i)) {
          utgw(i,k) = ubt(i,k) * xv(i) * effgw * ptaper(i);
          vtgw(i,k) = ubt(i,k) * yv(i) * effgw * ptaper(i);
        }
      }
    }

  }
}

} // namespace gw
} // namespace scream

#endif
